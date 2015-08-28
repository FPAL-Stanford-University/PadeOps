module initialization
    use kind_parameters, only: rkind
    use exits, only: GracefulExit 
    use mpi
    use variables
    use decomp_2d
    use io, only: get_input, getHit3d_uvw, write_matlab_header,closeMatlabFile,dumpData4Matlab 
    use constants, only: zero, one, two, three, four,pi, imi 
    use spectralOps, only: divergence 
    
    implicit none

    private
    public :: initialize, finalize  

    character(len=64) :: inputFilename 
    integer :: ierr


contains

    subroutine initialize(tidx)
        integer, intent(out) :: tidx 
        character(len=64) :: tmp
        !integer :: dims(2),dummy_periods(2), dummy_coords(2)
        real(rkind) :: kmax_sq
        logical, parameter :: fixOddball = .false. 
        integer :: i, j, k 

        ! Begin MPI 
        call MPI_Init(ierr)
        
        ! Get the input file 
        if(iargc().eq.0) then
            call getarg(0,tmp)
            call GracefulExit("Usage format:"//trim(tmp)//' {input File}', 00)
        end if
        call getarg(1,inputFilename)

        ! Read input data 
        call get_input(trim(inputFilename))

        dx = two*pi/Nx
        dy = two*pi/Ny
        dz = two*pi/Nz

        kdealias_x = real(floor(Nx/three),rkind) 
        kdealias_y = real(floor(Ny/three),rkind)
        kdealias_z = real(floor(Nz/three),rkind)

        nu = one/REY

        kmax_sq = (min(kdealias_x,kdealias_y,kdealias_z))**2 
        
        ! Get 2d decomposition 
        call decomp_2d_init(Nx, Ny, Nz, 0, 0)
        allocate(gp)
        call get_decomp_info(gp)
       
        ! Allocate physical space fields  
        allocate(fieldsPhys(gp%xsz(1),gp%xsz(2),gp%xsz(3),3))
        allocate(Buff_real(gp%xsz(1),gp%xsz(2),gp%xsz(3)))

        ! Allocate 3d FFT 
        allocate(FT)
        ierr = FT%init(Nx,Ny,Nz,"x", dx, dy,dz,FFT_PLAN_EXHAUSTIVE,.true.,.true.)
        if (ierr .ne. 0) then
            call GracefulExit("FFT_3d derived could not be initialized", 01)
        end if

        ! Allocate Spectral space variables (See fft_3d class to see how this done)
        call FT%alloc_output(rhs,3)
        call FT%alloc_output(rhsold,3)
        call FT%alloc_output(fieldsSpec,3)
        call FT%alloc_output(DealiasMat)
        call FT%alloc_output(oneBykSq)
        call FT%alloc_output(Buff_Cmplx)
  
        ! Store one-by-k_squared in memory for fast poisson inversion
        do k = 1,size(oneBykSq,3)
            do j = 1,size(oneBykSq,2)
                do i = 1,size(oneBykSq,1)
                    if (FT%kabs_sq(i,j,k) .lt. 1.d-18) then
                        oneBykSq(i,j,k) = zero
                    else
                        oneBykSq(i,j,k) = one/FT%kabs_sq(i,j,k)
                    end if 
                end do 
            end do 
         end do 

        !call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, ierr) 

        !if ((dims(1) == 1) .or. (dims(2) == 1)) then
        !    if (nrank == 0) then
        !        print*, "============================================================================="
        !        print*, "WARNING: Effectively 1D decomposition is being used. &
        !        & The code hasn't been optimized for 1D decomposition."
        !        print*, "============================================================================="
        !    end if 
        !end if 
        
        if (.not. RESTART) then  
            if (useHit3dinit) then
                call getHit3d_uvw(HIT3dinitDir) 
            else
                call GracefulExit("Stochastic initialization is currenlty not supported", 02)
            end if 
            if (nrank == 0) then
                print*, "Initialization complete!"
                print*, "============================================================================="
                print*, "Problem Summary:"
                print*, "---------------"
                write(*,"(A40,3I4)"), "(Nx, Ny, Nz) GLOBAL:", Nx, Ny, Nz
                write(*,"(A40,I4)"), "Total Processors:",     nproc
                !write(*,"(A40,I4,A4,I4)"), "Processor Decomposition:", dims(1),"x", dims(2)
                if (dealias == 2) then
                    write(*,"(A40,A32)") "Dealiasing method:", "Boxed dealiasing (2/3rd rule)"
                    write(*,"(A40,3F10.4)"), "Dealiasing Wavenumbers (kx,ky,kz)_max:", kdealias_x,kdealias_y,kdealias_z
                elseif (dealias == 1) then
                    write(*,"(A40,A36)") "Dealiasing method:", "Spherical dealiasing (2/3rd rule)"
                    write(*,"(A40,3F10.4)"), "Dealiasing Wavenumber Magnitude, k_max:", sqrt(kmax_sq)
                end if 
            end if
            time = zero 
            tidx = 0
        else
            call gracefulExit("Restart capability not yet supported", 03)
            
        endif

        ! Generate dealiasing matrix 
        select case (dealias)
        case(0) ! No dealiasing
            if (nrank == 0) then
                print*, "=========================================================="
                print*, "WARNING: You have selected to not perform dealiasing.&
                        & Solution likely to blow up."
                print*, "=========================================================="
            end if 
            DealiasMat = one
        case(1) ! Spherical 
            do k = 1,size(DealiasMat,3)
                do j = 1,size(DealiasMat,2)
                    do i = 1,size(DealiasMat,1)
                        if (FT%kabs_sq(i,j,k) .gt. kmax_sq) then
                            DealiasMat(i,j,k) = zero
                        else
                            DealiasMat(i,j,k) = one
                        end if 
                    end do 
                end do 
            end do 
        case(2) ! Box 
            do k = 1,size(DealiasMat,3)
                do j = 1,size(DealiasMat,2)
                    do i = 1,size(DealiasMat,1)
                        if ((FT%k1(i,j,k) .gt. kdealias_x) .or. &
                            (FT%k2(i,j,k) .gt. kdealias_y) .or. &
                            (FT%k3(i,j,k) .gt. kdealias_z)) then

                            DealiasMat(i,j,k) = zero
                        else
                            DealiasMat(i,j,k) = one
                        end if 
                    end do 
                end do 
            end do
        case default
           call GracefulExit("Invalid choice for Dealiasing method", 04) 
        end select
      

       ! Compute the Forward transforms 
       call FT%fft3_x2z(fieldsPhys(:,:,:,1),fieldsSpec(:,:,:,1))
       call FT%fft3_x2z(fieldsPhys(:,:,:,2),fieldsSpec(:,:,:,2))
       call FT%fft3_x2z(fieldsPhys(:,:,:,3),fieldsSpec(:,:,:,3))

       ! Perform explit dealiasing on input
       fieldsSpec(:,:,:,1) = DealiasMat*fieldsSpec(:,:,:,1)
       fieldsSpec(:,:,:,2) = DealiasMat*fieldsSpec(:,:,:,2)
       fieldsSpec(:,:,:,3) = DealiasMat*fieldsSpec(:,:,:,3)

       ! Check divergence 
       call divergence

       ! Return the the dealised fields back in Physical space 
       call FT%ifft3_z2x(fieldsSpec(:,:,:,1),fieldsPhys(:,:,:,1))
       call FT%ifft3_z2x(fieldsSpec(:,:,:,2),fieldsPhys(:,:,:,2))
       call FT%ifft3_z2x(fieldsSpec(:,:,:,3),fieldsPhys(:,:,:,3))
       
       if (nrank == 0) then
            write(*,"(A40,E13.4)"), "Maximum divergence of the initial field:", MaxDivergence
            write(*,"(A40,E13.4)"), "Time step:", dt
       end if 

        call write_matlab_header
        ! Create first dump for restart/initialization
        call dumpData4Matlab(tidx)


    end subroutine 

    subroutine finalize
   
        call closeMatlabFile 
        if ( allocated (fieldsPhys)) deallocate (fieldsPhys)
        if ( allocated (fieldsSpec)) deallocate (fieldsSpec)
        if ( allocated (DealiasMat)) deallocate (DealiasMat)
        if ( allocated (rhs)) deallocate (rhs)
        if ( allocated (rhsold)) deallocate (rhsold)

        if ( allocated (gp)) then
            call decomp_info_finalize(gp)
            deallocate(gp)
        end if 

        if (allocated(FT)) then
            call FT%destroy
            deallocate(FT)
        end if 

        call decomp_2d_finalize
        call MPI_Finalize(ierr)

    end subroutine

end module 
