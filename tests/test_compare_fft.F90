module getHitinput
    implicit none
    contains
    subroutine getHit3d_uvw(Nx,Ny,Nz,fieldsPhys,gp,dir)
        use kind_parameters, only: rkind, mpirkind
        use mpi
        use decomp_2d,       only: decomp_info,nrank,nproc
        use exits, only: GracefulExit, message 
        class( decomp_info ), intent(in)                  :: gp
        integer, dimension(:,:), allocatable        :: xst,xen,xsz
        integer, intent(in)                               :: Nx,Ny,Nz
        real(rkind), dimension(:,:,:,:), intent(inout)    :: fieldsPhys
        character(len=*), intent(in)                      :: dir
        character(len=8),parameter                        :: uFile = "U.000000"
        character(len=8),parameter                        :: vFile = "V.000000"
        character(len=8),parameter                        :: wFile = "W.000000"
        character(:), allocatable                         :: fname
        real(rkind), dimension(:,:,:), allocatable        :: full_Field
        integer :: fid = 1234, i, j, k
        integer :: tag, idx, status(MPI_STATUS_SIZE), ierr
        integer :: sizes(3), chunksize

        call message("Reading in Hit3D U, V, W data from directory: "//trim(dir))

        ! Create data sharing info
        if (nrank == 0) then
            allocate(xst(0:nproc-1,3),xen(0:nproc-1,3),xsz(0:nproc-1,3))
            allocate(full_Field(Nx,Ny,Nz))
        end if


        ! communicate local processor grid info (Assume x-decomposition)
        if (nrank == 0) then
            xst(0,:) = gp%xst
            xen(0,:) = gp%xen
            
            tag = 0
            do idx = 1,nproc-1
                call MPI_RECV(xst(idx,:), 3, MPI_INTEGER, idx, tag,&
                    MPI_COMM_WORLD, status, ierr)
            end do 
            tag = 1
            do idx = 1,nproc-1
                call MPI_RECV(xen(idx,:), 3, MPI_INTEGER, idx, tag,&
                    MPI_COMM_WORLD, status, ierr)
            end do
           tag = 2
            do idx = 1,nproc-1
                call MPI_RECV(xsz(idx,:), 3, MPI_INTEGER, idx, tag,&
                    MPI_COMM_WORLD, status, ierr)
            end do

        else
            tag = 0
            call MPI_SEND(gp%xst, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)
            tag = 1
            call MPI_SEND(gp%xen, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)
            tag = 2
            call MPI_SEND(gp%xsz, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)

        end if 
        
        ! Read the u data  
        if (nrank == 0) then
            ! U data 
            fname = dir(:len_trim(dir))//'/'//uFile
            open(fid,file=fname,form='unformatted', access='stream')  
            read(fid) sizes(1:3)
            if ((sizes(1) == Nx) .and. (sizes(2) == Ny) .and. (sizes(3) == Nz)) then
                read(fid) (((full_field(i,j,k),i=1,Nx),j=1,Ny),k=1,Nz)
            else
                close(fid)
                print*, "Problem size      :", Nx, Ny, Nz
                print*, "Hit 3d input size :", sizes 
                call GracefulExit("Problem size mismatch with HIT_3d input file:"//uFile,1001) 
            end if
            tag = 10
            fieldsPhys(:,:,:,1) = full_field(xst(0,1):xen(0,1), xst(0,2):xen(0,2), xst(0,3):xen(0,3))
            do idx = 1,nproc-1
                chunksize = xsz(idx,1)*xsz(idx,2)*xsz(idx,3)
                call MPI_SEND(full_field(xst(idx,1):xen(idx,1), xst(idx,2):xen(idx,2),xst(idx,3):xen(idx,3)), chunksize,&
                mpirkind, idx, tag, MPI_COMM_WORLD, ierr)
            end do

            ! V data 
            fname = dir(:len_trim(dir))//'/'//vFile
            open(fid,file=fname,form='unformatted', access='stream')  
            read(fid) sizes(1:3)
            if ((sizes(1) == Nx) .and. (sizes(2) == Ny) .and. (sizes(3) == Nz)) then
                read(fid) (((full_field(i,j,k),i=1,Nx),j=1,Ny),k=1,Nz)
            else
                close(fid)
                print*, "Problem size      :", Nx, Ny, Nz
                print*, "Hit 3d input size :", sizes 
                call GracefulExit("Problem size mismatch with HIT_3d input file:"//vFile,1002) 
            end if
            tag = 20
            fieldsPhys(:,:,:,2) = full_field(xst(0,1):xen(0,1), xst(0,2):xen(0,2), xst(0,3):xen(0,3))
            do idx = 1,nproc-1
                chunksize = xsz(idx,1)*xsz(idx,2)*xsz(idx,3)
                call MPI_SEND(full_field(xst(idx,1):xen(idx,1), xst(idx,2):xen(idx,2),xst(idx,3):xen(idx,3)), chunksize,&
                mpirkind, idx, tag, MPI_COMM_WORLD, ierr)
            end do

            ! W data 
            fname = dir(:len_trim(dir))//'/'//wFile
            open(fid,file=fname,form='unformatted', access='stream')  
            read(fid) sizes(1:3)
            if ((sizes(1) == Nx) .and. (sizes(2) == Ny) .and. (sizes(3) == Nz)) then
                read(fid) (((full_field(i,j,k),i=1,Nx),j=1,Ny),k=1,Nz)
            else
                close(fid)
                print*, "Problem size      :", Nx, Ny, Nz
                print*, "Hit 3d input size :", sizes 
                call GracefulExit("Problem size mismatch with HIT_3d input file:"//wFile,1003) 
            end if
            tag = 30
            fieldsPhys(:,:,:,3) = full_field(xst(0,1):xen(0,1), xst(0,2):xen(0,2), xst(0,3):xen(0,3))
            do idx = 1,nproc-1
                chunksize = xsz(idx,1)*xsz(idx,2)*xsz(idx,3)
                call MPI_SEND(full_field(xst(idx,1):xen(idx,1), xst(idx,2):xen(idx,2),xst(idx,3):xen(idx,3)), chunksize,&
                mpirkind, idx, tag, MPI_COMM_WORLD, ierr)
            end do

        else
            tag = 10
            chunksize = gp%xsz(1)*gp%xsz(2)*gp%xsz(3)
            call MPI_RECV(fieldsPhys(:,:,:,1), chunksize, mpirkind, 0, tag,&
                    MPI_COMM_WORLD, status, ierr)
            tag = 20
            call MPI_RECV(fieldsPhys(:,:,:,2), chunksize, mpirkind, 0, tag,&
                    MPI_COMM_WORLD, status, ierr)
            tag = 30
            call MPI_RECV(fieldsPhys(:,:,:,3), chunksize, mpirkind, 0, tag,&
                    MPI_COMM_WORLD, status, ierr)
            
        end if 

        if (nrank == 0) then
            deallocate(full_Field)
        end if
    end subroutine 

end module 

module fft3ders
    use fft_3d_stuff, only: fft_3d 
    use kind_parameters, only: rkind, clen 
    use constants, only: imi
    
    implicit none
    character(len=1) :: fft_dir = "x"
    type(fft_3d) :: FT 
    complex(rkind), dimension(:,:,:), allocatable :: fhat

    contains 
    subroutine get_ddx(arr_in,arr_out)
        real(rkind), dimension(:,:,:), intent(in) :: arr_in
        real(rkind), dimension(:,:,:), intent(out) :: arr_out

        select case (fft_dir)
        case("x")
            call FT% fft3_x2z(arr_in,fhat)
            fhat = imi*FT%k1*fhat
            call FT%ifft3_z2x(fhat,arr_out) 
        case("y")
            call FT% fft3_y2y(arr_in,fhat)
            fhat = imi*FT%k1*fhat
            call FT%ifft3_y2y(fhat,arr_out) 
        case("z")
            call FT% fft3_z2x(arr_in,fhat)
            fhat = imi*FT%k1*fhat
            call FT%ifft3_x2z(fhat,arr_out) 
        end select
    end subroutine

    subroutine get_ddy(arr_in,arr_out)
        real(rkind), dimension(:,:,:), intent(in) :: arr_in
        real(rkind), dimension(:,:,:), intent(out) :: arr_out

        select case (fft_dir)
        case("x")
            call FT% fft3_x2z(arr_in,fhat)
            fhat = imi*FT%k2*fhat
            call FT%ifft3_z2x(fhat,arr_out) 
        case("y")
            call FT% fft3_y2y(arr_in,fhat)
            fhat = imi*FT%k2*fhat
            call FT%ifft3_y2y(fhat,arr_out) 
        case("z")
            call FT% fft3_z2x(arr_in,fhat)
            fhat = imi*FT%k2*fhat
            call FT%ifft3_x2z(fhat,arr_out) 
        end select
    end subroutine
    
    subroutine get_ddz(arr_in,arr_out)
        real(rkind), dimension(:,:,:), intent(in) :: arr_in
        real(rkind), dimension(:,:,:), intent(out) :: arr_out

        select case (fft_dir)
        case("x")
            call FT% fft3_x2z(arr_in,fhat)
            fhat = imi*FT%k3*fhat
            call FT%ifft3_z2x(fhat,arr_out) 
        case("y")
            call FT% fft3_y2y(arr_in,fhat)
            fhat = imi*FT%k3*fhat
            call FT%ifft3_y2y(fhat,arr_out) 
        case("z")
            call FT% fft3_z2x(arr_in,fhat)
            fhat = imi*FT%k3*fhat
            call FT%ifft3_x2z(fhat,arr_out) 
        end select
    end subroutine



end module 

module derWrapper
    use kind_parameters, only: rkind
    use decomp_2d
    use DerivativesMod,  only: derivatives
    implicit none 
    
    type(derivatives) :: der
    

contains
    subroutine der_ddx(arr_in,arr_out,gp)
        class(decomp_info), intent(in) :: gp
        real(rkind), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)), intent(in) :: arr_in 
        real(rkind), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)), intent(out) :: arr_out
        
        call der%ddx(arr_in,arr_out)

    end subroutine 
    
    subroutine der_ddy(arr_in,arr_out,gp)
        class(decomp_info), intent(in) :: gp
        real(rkind), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)), intent(in) :: arr_in 
        real(rkind), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)), intent(out) :: arr_out
        real(rkind), dimension(gp%ysz(1),gp%ysz(2),gp%ysz(3)) :: tmp1, tmp2 
        
        call transpose_x_to_y(arr_in,tmp1,gp)
        call der%ddy(tmp1,tmp2)
        call transpose_y_to_x(tmp2,arr_out)

    end subroutine 

    subroutine der_ddz(arr_in,arr_out,gp)
        class(decomp_info), intent(in) :: gp
        real(rkind), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)), intent(in) :: arr_in 
        real(rkind), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)), intent(out) :: arr_out
        real(rkind), dimension(gp%ysz(1),gp%ysz(2),gp%ysz(3)) :: tmp
        real(rkind), dimension(gp%zsz(1),gp%zsz(2),gp%zsz(3)) :: tmp1, tmp2 
        
        call transpose_x_to_y(arr_in,tmp,gp)
        call transpose_y_to_z(tmp,tmp1,gp)
        call der%ddz(tmp1,tmp2)
        call transpose_z_to_y(tmp2,tmp,gp)
        call transpose_y_to_x(tmp,arr_out)

    end subroutine 

end module

program compare_ffts
    use fft_3d_stuff, only: fft_3d 
    use kind_parameters, only: rkind, clen 
    use decomp_2d, only: decomp_info, nrank, decomp_2d_init, get_decomp_info 
    use mpi 
    use constants, only: pi
    use fft3ders
    use derWrapper
    use spectralMod, only: spectral  
    use exits, only: GracefulExit, message
    use getHITinput, only: getHit3d_uvw 
    use timer, only: tic, toc
    implicit none

    type(spectral) :: SPECT
    integer :: nx, ny, nz
    real(rkind), dimension(:,:,:), allocatable :: buff1_R,buff2_R
    real(rkind), target, dimension(:,:,:,:), allocatable :: fields 
    complex(rkind), dimension(:,:,:), allocatable :: buff1_C,buff2_C
    real(rkind), dimension(:,:,:), pointer :: u, v, w
    type(decomp_info) :: gp
    character(len=clen) :: directory
    character(len=8),parameter                        :: uFile = "U.000000"
    character(len=8),parameter                        :: vFile = "V.000000"
    character(len=8),parameter                        :: wFile = "W.000000"
    integer ::  ierr, idx   
    real(Rkind) :: dx, dy, dz
    integer :: nrow = 0, ncols = 0
    integer :: ntimer = 200
    real(rkind) :: ttime

    directory = "../../data/hit/input"
    nx = 128; ny = 128; nz = 128;
    
    ! Do the decomp shit
    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, nrow, ncols)
    call get_decomp_info(gp)
     

    ! Allocate stuff 
    allocate(buff1_R(gp%xsz(1),gp%xsz(2),gp%xsz(3)),buff2_R(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(fields(gp%xsz(1),gp%xsz(2),gp%xsz(3),3)) 

    u => fields(:,:,:,1)
    v => fields(:,:,:,2)
    w => fields(:,:,:,3)

    ! Read in the data 
    call getHit3d_uvw(Nx,Ny,Nz,fields,gp,directory) 

    dx = 2.*pi/nx 
    dy = 2.*pi/ny
    dz = 2.*pi/nz

    ! 1st kind - use the derivatives Module

    ! Initialize other shit
    call     der%init(          gp,& 
                                dx,     dy,    dz, &
                            .TRUE., .TRUE., .TRUE., &
                            "four", "four", "four" )
    
    ! 2nd kind - use the spectral module
    call SPECT%init(fft_dir, nx, ny, nz, dx, dy, dz, "four", "2/3rd",3, .true., .true.)
    call SPECT%alloc_r2c_out(buff1_c) 
    call SPECT%alloc_r2c_out(buff2_c) 
    
    ! 3rd kind - use the 3d FFT module 
    ierr = FT%init(nx,ny,nz,fft_dir,dx, dy,dz, .false., .true.)
    call FT%alloc_output(fhat) 

  
    ! Create input 
    buff1_R = u*v

    !!!! Test derivative !!!!
    ! 1st kind
    call tic
    do idx = 1,ntimer 
        call der_ddz(buff1_R,buff2_R,gp)
    end do 
    call toc("DERIVATIVES Mod:",ttime)
    !call sleep(nrank)
    !print*, buff2_R(1:2,1,1)
    call mpi_barrier(mpi_comm_world,ierr)

    ! 3rd kind
    call tic
    do idx = 1,ntimer 
        call get_ddz(buff1_R,buff2_R)
    end do 
    call toc("FFT_3d Mod:",ttime)
    !call sleep(nrank)
    !print*, buff2_R(1:2,1,1)
    call mpi_barrier(mpi_comm_world,ierr)
    
    ! 2nd kind 
    buff2_c = 0._rkind
    buff2_r = 0._rkind
    call tic 
    do idx = 1,ntimer
        call SPECT%fft(buff1_R,buff1_C)
        buff2_C = imi*SPECT%k3*buff1_C
        call SPECT%ifft(buff2_C,buff2_R)
    end do 
    call toc("SPECTRAL Mod:",ttime)
    !call sleep(nrank)
    !print*, buff2_R(1:2,1,1)
    call mpi_barrier(mpi_comm_world,ierr)
   
    
    call mpi_finalize(ierr)

end program 
