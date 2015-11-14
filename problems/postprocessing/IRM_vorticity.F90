module IRM_vorticity_mod
    use mpi
    use kind_parameters, only: rkind, clen, mpirkind, mpickind
    use constants,       only: zero, one
    use fftstuff,        only: ffts
    use miranda_tools,   only: miranda_reader
    use io_VTK_stuff,    only: io_VTK
    use DerivativesMod,  only: derivatives
    use decomp_2d,       only: nrank, nproc, transpose_y_to_z
    use exits,           only: message, GracefulExit

    implicit none

    type(miranda_reader) :: mir
    type(io_VTK)         :: viz
    type(derivatives)    :: der
    type(ffts)           :: fftz

    character(len=clen)  :: inputdir, outputdir                                                 ! Input and output dirs
    integer              :: prow, pcol                                                          ! Procs in x and z
    logical              :: periodicx = .FALSE., periodicy = .FALSE., periodicz = .FALSE.       ! Periodic?
    character(len=4)     :: derivative_x = 'cd10', derivative_y = 'cd10', derivative_z = 'cd10' ! Derivative method to use

    logical              :: writeviz                                                            ! Write out VTK viz files?

    ! This part is from the Miranda problem file to recalculate the physical species diffusivity
    ! Material properties list ------ Property --------  Air ------------ CO2 --------
    real(rkind), DIMENSION(2) :: mySig    = (/ 3.681D0            , 3.952D0   /)
    real(rkind), DIMENSION(2) :: myTeff   = (/ 91.46D0            , 200.0D0   /)
    real(rkind), DIMENSION(2) :: myMolwts = (/ 28.013D0           , 44.010D0  /)
    real(rkind), DIMENSION(2) :: myGammas = (/ 1.4D0              , 1.28D0    /)
    real(rkind), DIMENSION(2) :: myPr     = (/ 0.72D0             , 0.77D0    /)
    real(rkind), DIMENSION(2) :: myMix    = (/ 1.0D0              , 1.0D0     /)
    real(rkind), DIMENSION(2) :: myNuA    = (/ 1.45D-5            , 1.59D-5   /) ! in g/(cm-s-K^0.5)
    real(rkind), DIMENSION(2) :: myNuB    = (/ 1.105D2            , 2.438D2   /) ! in K
    real(rkind), DIMENSION(2) :: myNuC    = (/ 0.0D0              , 0.0D0     /) ! in g/(cm-s)
    real(rkind), DIMENSION(2) :: myNuN    = (/ 1.5D0              , 1.5D0     /)

    real(rkind), PARAMETER :: Runiv  = 8.314472D+7              ! Gas constant              [g*cm^2/(s^2*mol*K)]=[erg/K/mol]

contains

    subroutine read_inputs(inputfile)
        character(len=clen), intent(in) :: inputfile
        integer :: iounit = 67

        namelist /INPUT/ inputdir, outputdir, &
                         periodicx, periodicy, periodicz, &
                         derivative_x, derivative_y, derivative_z, &
                         prow, pcol, writeviz

        open(unit=iounit, file=trim(inputfile), form='FORMATTED')
        read(unit=iounit, NML=INPUT)
        close(iounit)
    end subroutine

    subroutine write_post2d(step, vort_avgz, TKE_avg, div_avg, rho_avg, chi_avg, MMF_avg, CO2_avg, density_self_correlation)
        integer, intent(in) :: step
        real(rkind), dimension(:,:,:), intent(in) :: vort_avgz
        real(rkind), dimension(:,:),   intent(in) :: TKE_avg, div_avg, rho_avg, chi_avg, MMF_avg, CO2_avg, density_self_correlation

        character(len=clen) :: post2d_file
        integer :: i, ierr, iounit2d=91

        write(post2d_file,'(2A,I4.4,A)') trim(outputdir), "/post2d_", step, ".dat"
        call message(1,"Writing 2D vorticity file to "//trim(post2d_file))
        if(nrank == 0) then
            open(unit=iounit2d, file=trim(post2d_file), form='UNFORMATTED', status='REPLACE')
            write(iounit2d) mir%nx, mir%ny
            write(iounit2d) mir%gp%ysz(1), mir%gp%ysz(2)
            close(iounit2d)
        end if
        do i = 0,nproc-1
            if (nrank == i) then
                if( mir%gp%yst(3) == 1 ) then
                    open(unit=iounit2d, file=trim(post2d_file), form='UNFORMATTED', status='OLD', position='APPEND')
                    write(iounit2d) mir%gp%yst(1), mir%gp%yen(1), mir%gp%yst(2), mir%gp%yen(2)
                    write(iounit2d) mir%x(:,:,1)
                    write(iounit2d) mir%y(:,:,1)
                    write(iounit2d) vort_avgz(:,:,3)
                    write(iounit2d) TKE_avg
                    write(iounit2d) div_avg
                    write(iounit2d) rho_avg
                    write(iounit2d) chi_avg
                    write(iounit2d) MMF_avg
                    write(iounit2d) CO2_avg
                    write(iounit2d) density_self_correlation
                    close(iounit2d)
                end if
            end if
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        end do
    end subroutine

    subroutine get_zfft(input,output,mask)
        real(rkind), dimension(:,:,:), intent(in) :: input
        real(rkind), dimension(:,:,:), intent(in) :: mask
        complex(rkind), dimension(:), intent(out) :: output

        real(rkind), dimension(mir%gp%zsz(1),mir%gp%zsz(2),mir%gp%zsz(3)) :: zbuffer
        complex(rkind), dimension(mir%gp%zsz(1),mir%gp%zsz(2),mir%gp%zsz(3)/2+1) :: zfft
        complex(rkind), dimension(mir%gp%zsz(3)/2+1) :: proc_output

        real(rkind) :: masksum, proc_masksum
        integer :: i,j,ierr

        ! Get FFT of input
        call transpose_y_to_z(input, zbuffer, mir%gp)
        call fftz%fftz(zbuffer, zfft)

        ! Transpose mask to Z to average the FFT
        call transpose_y_to_z( mask, zbuffer, mir%gp)
        
        proc_output = zero
        proc_masksum = zero
        do j = 1,mir%gp%zsz(2)
            do i = 1,mir%gp%zsz(1)
                proc_output  = proc_output + zfft(i,j,:) * zbuffer(i,j,1)
                proc_masksum = proc_masksum + zbuffer(i,j,1)
            end do
        end do

        call MPI_ALLREDUCE(proc_output,  output,  mir%gp%zsz(3)/2+1, mpickind, MPI_SUM, MPI_COMM_WORLD, ierr) ! Add proc_output  across processors
        call MPI_ALLREDUCE(proc_masksum, masksum,                 1, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr) ! Add proc_masksum across processors

        output = ( output / masksum ) / mir%nz ! Divide by nz to get physical FFT

    end subroutine

    subroutine write_spectrum(step, varname, fft)
        integer, intent(in) :: step
        character(len=*), intent(in) :: varname
        complex(rkind), dimension(:), intent(in) :: fft

        character(len=clen) :: filename
        integer :: i,n,iounit=21

        ! Assume all processors have the same fft array, so only master proc can write out the file

        n = size(fft)

        write(filename,'(4A,I4.4,A)') trim(outputdir), "/", trim(varname), "_", step, ".dat"
        call message(1,"Writing " // trim(varname) // " spectrum to " // trim(filename))
        if (nrank == 0) then
            open(unit=iounit, file=trim(filename), form='FORMATTED', status="REPLACE")

            write(iounit,'(2A)') "# Spectrum for ", trim(varname)
            write(iounit,'(A9,2A26)') "# wav no.", "Real part", "Imag part"
            do i=1,n
                write(iounit,'(I9,2ES26.16)') (i-1), real(fft(i)), aimag(fft(i))
            end do

            close(iounit)
        end if

    end subroutine

    subroutine get_volumefractions(Y,X)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), mir%ns), intent(in)  :: Y    ! Species mass fractions
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), mir%ns), intent(out) :: X    ! Species volume fractions
        
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) :: M
        
        integer :: n

        ! Effective Molecular weight for EOS
        M = zero
        DO n = 1,mir%ns
           M = M + Y(:,:,:,n)/myMolwts(n)
        END DO
        M = one/M
        
        ! Get volume fractons
        do n = 1,mir%ns
            X(:,:,:,n) = Y(:,:,:,n)*M/myMolwts(n)
        end do
   
    end subroutine

! ------------------------------------------------------------------------------
! Assign problem-specific viscosity, thermal conductivity,
! species diffusivities and magnetic diffusivity.
! ------------------------------------------------------------------------------ 
    SUBROUTINE prob_properties(mu,ktc,Diff)
     USE constants, ONLY : half,three,two

     real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)),                   intent(out) :: mu,ktc    ! Shear visc and conductivity
     real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), mir%ns),           intent(out) :: Diff      ! Species diffusivity

     real(rkind), dimension(:,:,:),   pointer :: T,rho,p
     real(rkind), dimension(:,:,:,:), pointer :: Y
     
     real(rkind), DIMENSION(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) :: M,CpMix,Rgas,Ts,Omega
     real(rkind), DIMENSION(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) :: mu_i,tmp,dum,Dij,Xii,Xjj
     real(rkind) :: Acc,Bcc,Ccc,Dcc,Ecc,Fcc,Gcc,Hcc
     real(rkind) :: Mij,sigij,Teij
     INTEGER :: n,ni,nj
  
     ! Point T, Y, rho and p to mir variables so we can use Miranda code here directly
     T   => mir%T
     Y   => mir%Ys
     rho => mir%rho
     p   => mir%p

     ! Effective Molecular weight for EOS
     M = zero
     DO n = 1,mir%ns
        M = M + Y(:,:,:,n)/myMolwts(n)
     END DO
     M = one/M
     Rgas = Runiv/M
   
   
     ! Chapman-Cowling Viscosity
     mu = zero
     ktc = zero
     tmp = zero
     Acc = 1.16145D0
     Bcc = -0.14874D0
     Ccc = 0.52487D0
     Dcc = -0.7732D0
     Ecc = 2.16178D0
     Fcc = -2.43787D0
     DO n=1,mir%ns
         Ts = T/myTeff(n)
         Omega = Acc*Ts**Bcc + Ccc*exp(Dcc*Ts) + Ecc*exp(Fcc*Ts)     ! Edit: Ecc** -> Ecc*
         mu_i = 2.6693D-6*sqrt( myMolwts(n)*T)/(Omega*mySig(n)**2)
         dum = Y(:,:,:,n)/(myMolwts(n)**half)
         mu = mu + mu_i * dum
         tmp = tmp + dum
   
        ! Thermal diff part
         CpMix = myGammas(n)*Y(:,:,:,n)/(myMolwts(n)*(myGammas(n)-one))*Runiv
         ktc = ktc + CpMix*mu_i*dum/myPr(n)
   
     END DO
     mu = mu/tmp * 10.0D0  ! convert kg/m s to g/cm s
     ktc = ktc/tmp
   
     ! Diffusivity
     Diff = zero
     tmp = zero
     Acc = 1.06036D0
     Bcc = -0.1561D0
     Ccc = 0.19300D0
     Dcc = -0.47635D0
     Ecc = 1.03587D0
     Fcc = -1.52996D0
     Gcc = 1.76474D0
     Hcc = -3.89411D0
     DO ni=1,mir%ns
         Xii = Y(:,:,:,ni)*M/myMolwts(ni)
         tmp = zero
         DO nj=1,mir%ns
            IF (ni .NE. nj) THEN
               Xjj = Y(:,:,:,nj)*M/myMolwts(nj)
               Mij = two/(one/myMolwts(ni)+one/myMolwts(nj))
               sigij = (mySig(ni)+mySig(nj))/two
               Teij = sqrt( myTeff(ni)*myTeff(nj) )
               Ts = T/Teij                                         ! Edit: Add omega calculation
               Omega = Acc*Ts**Bcc + Ccc*exp(Dcc*Ts) + Ecc*exp(Fcc*Ts) + Gcc*exp(Hcc*Ts)
               Dij = 0.0266/Omega * T**(three/two)/(p*sqrt(Mij)*sigij**2 )
               tmp = tmp + Xjj/Dij
            END IF
         END DO
         Diff(:,:,:,ni) = (one-Xii)/(tmp + 1.0D-16)
     END DO
   
     Diff = Diff * 10.0D4 / 10.0D0 ! Convert m^2/s to cm^2/s  (extra /10 is from pressure)
   
    END SUBROUTINE prob_properties


end module

program IRM_vorticity
    use mpi
    use kind_parameters, only: rkind, clen
    use constants,       only: zero, half, one, eps
    use miranda_tools,   only: miranda_reader
    use io_VTK_stuff,    only: io_VTK
    use DerivativesMod,  only: derivatives
    use gridtools,       only: alloc_buffs
    use operators,       only: curl, gradient, divergence
    use reductions,      only: P_AVGZ, P_MAXVAL, P_MINVAL, P_SUM
    use timer,           only: tic, toc
    use exits,           only: message, GracefulExit

    use IRM_vorticity_mod

    implicit none

    character(len=clen) :: inputfile

    character(len=clen), dimension(:), allocatable :: varnames

    real(rkind) :: Ymin = 0.01_rkind  ! Cutoff for interface region

    real(rkind), dimension(:,:,:,:), allocatable, target :: buffer
    real(rkind), dimension(:,:,:,:), pointer :: vort, Diff, Xs
    real(rkind), dimension(:,:,:), pointer :: upp, vpp, wpp, TKE, div, region, chi, mu, ktc

    real(rkind), dimension(:,:,:), allocatable, target :: buffer2d
    real(rkind), dimension(:,:,:), pointer :: vort_avgz
    real(rkind), dimension(:,:), pointer :: rho_avg, u_avg, v_avg, w_avg, CO2_avg, TKE_avg, div_avg, chi_avg, MMF_avg, density_self_correlation

    real(rkind), dimension(:), allocatable :: Y_air_x, Y_CO2_x

    complex(rkind), dimension(:), allocatable :: TKEfft

    real(rkind) :: MMF_int, chi_int, vortz_int, vortz_pos, vortz_neg, TKE_int, mwidth

    logical :: fftw_exhaustive = .FALSE.

    character(len=clen) :: time_message

    character(len=clen) :: outputfile
    integer :: iounit = 93

    integer :: nVTKvars
    integer :: ierr, step, i

    !======================================================================================================!
    ! Note about usage: All 3D arrays are in Y decomposition.
    !                   All 2D arrays are also in Y decomposition (effectively 1D decomp)
    !                   All spectra are in Z, averaged over a region and are stored on every processor
    !                   It is assumed in some places that the number of species is 2
    !======================================================================================================!

    call MPI_Init(ierr)

    call tic()

    if( command_argument_count() .LT. 1 ) then
        call GracefulExit("Usage: "//NEW_LINE('A')//"    mpiexec -n 8 ./test_Miranda_reader <input file>", 1729)
    end if

    call get_command_argument(1,inputfile)
    call read_inputs(inputfile)

    ! Initialize miranda_reader object
    call mir%init(inputdir, prow, pcol, periodicx, periodicy, periodicz)

    call message("Jobdir is " // adjustl(trim(inputdir)) )
    call message("Initialized the Miranda object" )
    
    call message("Reading in the grid" )
    call mir%read_grid()
    
    call message("Initializing the derivative object" )
    call der%init(                                      mir%gp, &
                          mir%dx,        mir%dy,        mir%dz, &
                   mir%periodicx, mir%periodicy, mir%periodicz, &
                    derivative_x,  derivative_y,  derivative_z  )

    
    ! Initialize the FFT object (in the Z decomposition since that is where FFT will be performed)
    call message("Initializing the FFT object" )
    ierr = fftz%init( mir%nz, 'z', mir%gp%zsz(1), mir%gp%zsz(2), mir%dz, fftw_exhaustive )

    ! Allocate 3D buffer and associate pointers for convenience
    call alloc_buffs(buffer, 13+2*mir%ns-1, 'y', mir%gp)
    vort => buffer(:,:,:,1:3)
    upp => buffer(:,:,:,4)
    vpp => buffer(:,:,:,5)
    wpp => buffer(:,:,:,6)
    TKE => buffer(:,:,:,7)
    div => buffer(:,:,:,8)
    region => buffer(:,:,:,9)
    chi => buffer(:,:,:,10)
    mu => buffer(:,:,:,11)
    ktc => buffer(:,:,:,12)
    Diff => buffer(:,:,:,13:13+mir%ns-1)
    Xs => buffer(:,:,:,13+mir%ns:13+2*mir%ns-1)

    ! Allocate 2D buffer and associate pointers for convenience
    allocate( buffer2d( mir%gp%ysz(1), mir%gp%ysz(2), 13 ) )
    vort_avgz => buffer2d(:,:,1:3)
    rho_avg => buffer2d(:,:,4)
    u_avg => buffer2d(:,:,5)
    v_avg => buffer2d(:,:,6)
    w_avg => buffer2d(:,:,7)
    CO2_avg => buffer2d(:,:,8)
    TKE_avg => buffer2d(:,:,9)
    div_avg => buffer2d(:,:,10)
    chi_avg => buffer2d(:,:,11)
    MMF_avg => buffer2d(:,:,12)
    density_self_correlation => buffer2d(:,:,13)

    allocate( TKEfft(mir%nz/2+1) )

    allocate( Y_air_x( mir%gp%ysz(1) ) )
    allocate( Y_CO2_x( mir%gp%ysz(1) ) )

    ! Initialize visualization stuff
    if ( writeviz ) then
        nVTKvars = 10
        allocate( varnames(nVTKvars) )
        varnames( 1) = 'X-vorticity'
        varnames( 2) = 'Y-vorticity'
        varnames( 3) = 'Z-vorticity'
        varnames( 4) = 'upp'
        varnames( 5) = 'vpp'
        varnames( 6) = 'wpp'
        varnames( 7) = 'TKE'
        varnames( 8) = 'divergence'
        varnames( 9) = 'region'
        varnames(10) = 'scalar_dissipation'
        call message("Initializing the VTK I/O object" )
        call viz%init(outputdir, 'vorticity', nVTKvars, varnames)
    end if
    
    ! Open file for scalar outputs
    write(outputfile,'(2A)') trim(outputdir), "/post_scalar.dat"
    if(nrank == 0) then
        open(unit=iounit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
        write(iounit,'(A,A25,7A26)') "#", "Time", "Mixed width", "Z vorticity", "+ve Z vorticity", "-ve Z vorticity", "TKE", "Scalar dissipation", "MMF"
    end if

    call toc("Time to initialize everything: ")

    ! Loop through viz dumps
    do step = 0,mir%nsteps-1
        call tic()  ! Start timer

        call message(0,"Reading in data for step",step)

        ! Read in data for this viz dump
        call mir%read_data(step)

        ! Get the physical mu, ktc and Diff
        call prob_properties(mu,ktc,Diff)
        
        ! Get volume fractions
        call get_volumefractions(mir%Ys,Xs)
        
        ! Get the interface region
        region = zero
        where( (mir%Ys(:,:,:,1) .GT. Ymin) .AND. (mir%Ys(:,:,:,1) .LT. (one-Ymin)) )
            region = one
        end where
        call P_AVGZ( mir%gp, region, div_avg ) ! Put avg of region in div avg for now
        where ( div_avg .GT. one/(real(2*mir%nz,rkind)) )
            div_avg = one
        elsewhere
            div_avg = zero
        end where
        do i = 1,mir%gp%ysz(3)
            region(:,:,i) = div_avg
        end do

        ! Clip volume fractions above and below bounds
        where ( Xs .GT. one )
            Xs = one - eps
        elsewhere ( Xs .LT. zero )
            Xs = zero + eps
        end where

        ! Now get the molecular mixing fraction
        call P_AVGZ( mir%gp, Xs(:,:,:,1)*Xs(:,:,:,2), MMF_avg)
        call P_AVGZ( mir%gp, Xs(:,:,:,1), div_avg)  ! Temporarily store in div_avg
        call P_AVGZ( mir%gp, Xs(:,:,:,2), TKE_avg)  ! Temporarily store in TKE_avg

        ! Get integrated MMF (integrated inside region)
        MMF_int = P_SUM( MMF_avg*region(:,:,1) ) / P_SUM( div_avg*TKE_avg*region(:,:,1) )
        call message(1,"MMF",MMF_int)
        
        where ( div_avg*TKE_avg == zero )
            MMF_avg = one
        elsewhere
            MMF_avg = MMF_avg / ( div_avg*TKE_avg ) ! MMF = <X_air X_CO2>/( <X_air> <X_CO2> )
        end where

        ! Get mixing width
        call P_AVGZ( mir%gp, mir%Ys(:,:,:,1), div_avg) ! Temporarily store in div_avg
        Y_air_x = SUM(div_avg, 2)/real(mir%ny,rkind)   ! Average along the Y direction

        call P_AVGZ( mir%gp, mir%Ys(:,:,:,2), TKE_avg) ! Temporarily store in TKE_avg
        Y_CO2_x = SUM(TKE_avg, 2)/real(mir%ny,rkind)   ! Average along the Y direction
        
        mwidth = P_SUM( Y_air_x*Y_CO2_x ) * mir%dx     ! Integrate in X (\int <Y_air><Y_CO2> dx)
        mwidth = mwidth / real(pcol,rkind)      ! Divide by the number of Z processors due to multiple counting

        ! Get vorticity
        call curl( mir%gp, der, mir%u, mir%v, mir%w, vort)
        do i = 1,3
            call P_AVGZ( mir%gp, vort(:,:,:,i), vort_avgz(:,:,i) )
        end do
        
        ! Get integrated Z vorticity (integrated inside region)
        vortz_int = P_SUM(vort(:,:,:,3))*mir%dx*mir%dy*mir%dz
        call message(1,"vortz_int",vortz_int)
        chi = vort(:,:,:,3)
        where ( chi .LT. zero )
            chi = zero
        end where
        vortz_pos = P_SUM(chi)*mir%dx*mir%dy*mir%dz
        chi = vort(:,:,:,3)
        where ( chi .GT. zero )
            chi = zero
        end where
        vortz_neg = P_SUM(chi)*mir%dx*mir%dy*mir%dz

        ! Get rho average
        call P_AVGZ( mir%gp, mir%rho, rho_avg )

        ! Get density self correlation
        call P_AVGZ( mir%gp, one/mir%rho, density_self_correlation)
        density_self_correlation = rho_avg*density_self_correlation - one

        ! Get scalar dissipation rate for CO2 (use upp, vpp, wpp to temporarily store grad{Y_CO2})
        associate( dY_x => upp, dY_y => vpp, dY_z => wpp )
            call gradient(mir%gp, der, mir%Ys(:,:,:,2), dY_x, dY_y, dY_z)

            chi = dY_x*dY_x + dY_y*dY_y + dY_z*dY_z ! grad(Y_CO2).grad(Y_CO2)
            chi = Diff(:,:,:,2)*chi ! Scalar dissipation rate = Diff_CO2 * grad(Y_CO2).grad(Y_CO2)

            ! Get dissipation rate spectrum (put in TKEfft)
            call get_zfft(chi,TKEfft,region)
            call write_spectrum(step, "scalar_dissipation", TKEfft)

            ! Get Scalar dissipation average
            call P_AVGZ( mir%gp, chi, chi_avg )
        end associate
        
        ! Get integrated chi (integrated inside region)
        chi_int = P_SUM(chi)*mir%dx*mir%dy*mir%dz
        call message(1,"chi_int",chi_int)

        ! Get Favre averaged velocities
        upp = mir%rho * mir%u
        vpp = mir%rho * mir%v
        wpp = mir%rho * mir%w
        call P_AVGZ( mir%gp, upp, u_avg )
        call P_AVGZ( mir%gp, vpp, v_avg )
        call P_AVGZ( mir%gp, wpp, w_avg )
        u_avg = u_avg / rho_avg         ! <rho*u> / <rho>
        v_avg = v_avg / rho_avg         ! <rho*v> / <rho>
        w_avg = w_avg / rho_avg         ! <rho*w> / <rho>

        ! Get velocity fluctuations
        do i = 1,mir%gp%ysz(3)
            upp(:,:,i) = mir%u(:,:,i) - u_avg
            vpp(:,:,i) = mir%v(:,:,i) - v_avg
            wpp(:,:,i) = mir%w(:,:,i) - w_avg
        end do

        ! Get TKE = rho*u_i*u_i/2
        TKE = half * mir%rho * ( upp*upp + vpp*vpp + wpp*wpp )
        
        ! Get integrated TKE (integrated inside region)
        TKE_int = P_SUM(TKE)*mir%dx*mir%dy*mir%dz
        call message(1,"TKE_int",TKE_int)

        ! Get TKE spectrum
        call get_zfft(TKE,TKEfft,region)
        call write_spectrum(step, "TKE", TKEfft)
        
        ! Get TKE average
        call P_AVGZ( mir%gp, TKE, TKE_avg )

        ! Get Y_CO2 average
        call P_AVGZ( mir%gp, mir%Ys(:,:,:,2), CO2_avg )

        ! Get Divergence
        call divergence( mir%gp, der, mir%u, mir%v, mir%w, div )

        ! Get Divergence average
        call P_AVGZ( mir%gp, div, div_avg )

        if(nrank == 0) then
            write(iounit,'(8ES26.16)') real(step*1.0d-4,rkind), mwidth, vortz_int, vortz_pos, vortz_neg, TKE_int, chi_int, MMF_int
        end if

        ! Write out visualization files
        if ( writeviz ) then
            ! Set vizcount to be same as Miranda step
            call viz%SetVizcount(step)
            call viz%WriteViz(mir%gp, mir%mesh, buffer(:,:,:,1:nVTKvars), Xs, ['volume_fraction_1','volume_fraction_2'])
        end if

        ! Write out 2D postprocessing file
        call write_post2d(step, vort_avgz, TKE_avg, div_avg, rho_avg, chi_avg, MMF_avg, CO2_avg, density_self_correlation)

        write(time_message,'(A,I,A)') "Time to postprocess step ", step, " :"
        call toc(trim(time_message))
    end do

    ! Close file for scalar inputs
    if(nrank == 0) close(iounit)

    ! Destroy all variables and exit cleanly
    if (allocated(buffer  )) deallocate( buffer   )
    if (allocated(buffer2d)) deallocate( buffer2d )
    if (allocated( TKEfft )) deallocate( TKEfft   )
    if (allocated( Y_air_x)) deallocate( Y_air_x  )
    if (allocated( Y_CO2_x)) deallocate( Y_CO2_x  )
    call der%destroy()
    call mir%destroy()
    if ( writeviz ) then
        deallocate( varnames )
        call viz%destroy()
    end if
    call MPI_Finalize(ierr)

end program 
