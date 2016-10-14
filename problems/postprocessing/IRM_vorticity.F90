module IRM_vorticity_mod
    use mpi
    use kind_parameters,  only: rkind, clen, mpirkind, mpickind
    use constants,        only: zero, eps, half, one, three
    use fftstuff,         only: ffts
    use miranda_tools,    only: miranda_reader
    use io_VTK_stuff,     only: io_VTK
    use DerivativesMod,   only: derivatives
    use FiltersMod,       only: filters
    use decomp_2d,        only: nrank, nproc, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y, &
                                decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize
    use mytranspose2DMod, only: mytranspose2D
    use exits,            only: message, warning, GracefulExit
    use reductions,       only: P_MAXVAL, P_AVGZ

    implicit none

    type(miranda_reader) :: mir
    type(io_VTK)         :: viz
    type(derivatives)    :: der, der2D
    type(filters)        :: gfil
    type(ffts)           :: fftz
    type(mytranspose2D)  :: gp2D

    character(len=clen)  :: inputdir, outputdir                                                 ! Input and output dirs
    integer              :: prow, pcol                                                          ! Procs in x and z
    logical              :: periodicx = .FALSE., periodicy = .FALSE., periodicz = .FALSE.       ! Periodic?
    character(len=4)     :: derivative_x = 'cd10', derivative_y = 'cd10', derivative_z = 'cd10' ! Derivative method to use

    logical              :: writeviz                                                            ! Write out VTK viz files?

    integer              :: XY_COMM
    integer              :: xyrank, xyproc

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

    real(rkind), PARAMETER :: Cd = real(3.0D-2,rkind)
    real(rkind), PARAMETER :: CY = real(1.0D2 ,rkind)

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

    subroutine write_post2d(step, vort_avgz, TKE_avg, div_avg, rho_avg, chi_avg, MMF_avg, CO2_avg, density_self_correlation, &
                            R11, R12, R13, R22, R23, R33, a_x, a_y, a_z, rhop_sq, eta, xi, &
                            pdil, meandiss, production, dissipation, pdil_fluct, tke_visc_transport, mix_pdil, mix_pdil_fluct, &
                            duidxj2D, gradp2D, PSij, diss_tensor)
        integer, intent(in) :: step
        real(rkind), dimension(:,:,:), intent(in) :: vort_avgz
        real(rkind), dimension(:,:),   intent(in) :: TKE_avg, div_avg, rho_avg, chi_avg, MMF_avg, CO2_avg, density_self_correlation
        real(rkind), dimension(:,:),   intent(in) :: R11, R12, R13, R22, R23, R33, a_x, a_y, a_z, rhop_sq, eta, xi
        real(rkind), dimension(:,:),   intent(in) :: pdil, meandiss, production, dissipation, pdil_fluct, tke_visc_transport, mix_pdil, mix_pdil_fluct
        real(rkind), dimension(:,:,:), intent(in) :: duidxj2D, gradp2D, PSij, diss_tensor

        character(len=clen) :: post2d_file
        integer :: i, idx, ierr, iounit2d=91

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
                    write(iounit2d) R11
                    write(iounit2d) R12
                    write(iounit2d) R13
                    write(iounit2d) R22
                    write(iounit2d) R23
                    write(iounit2d) R33
                    write(iounit2d) a_x
                    write(iounit2d) a_y
                    write(iounit2d) a_z
                    write(iounit2d) rhop_sq
                    write(iounit2d) eta
                    write(iounit2d) xi
                    write(iounit2d) pdil
                    write(iounit2d) meandiss
                    write(iounit2d) production
                    write(iounit2d) dissipation
                    write(iounit2d) pdil_fluct
                    write(iounit2d) tke_visc_transport
                    write(iounit2d) mix_pdil
                    write(iounit2d) mix_pdil_fluct

                    do idx = 1,9
                        write(iounit2d) duidxj2D(:,:,idx)
                    end do

                    do idx = 1,3
                        write(iounit2d) gradp2D(:,:,idx)
                    end do

                    do idx=1,6
                        write(iounit2d) PSij(:,:,idx)
                    end do
                    do idx=1,6
                        write(iounit2d) diss_tensor(:,:,idx)
                    end do
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

        write(filename,'(4A,I4.4,A)') trim(outputdir), "/spectrum_", trim(varname), "_", step, ".dat"
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

    subroutine get_conditional_zfft(input,Ys,nbins,step)
        real(rkind), dimension(:,:,:), intent(in)  :: input
        real(rkind), dimension(:,:,:), intent(in)  :: Ys
        integer,                       intent(in)  :: nbins, step

        complex(rkind), dimension(mir%gp%zsz(3)/2+1) :: output
        real(rkind), dimension(mir%gp%zsz(1),mir%gp%zsz(2),mir%gp%zsz(3)) :: zbuffer
        complex(rkind), dimension(mir%gp%zsz(1),mir%gp%zsz(2),mir%gp%zsz(3)/2+1) :: zfft
        complex(rkind), dimension(mir%gp%zsz(3)/2+1) :: proc_output

        real(rkind), dimension(mir%gp%zsz(1),mir%gp%zsz(2)) :: Ysavg
        real(rkind) :: Yssum, proc_Yssum, Ysmin, Ysmax, dYs, Ysbin
        integer :: i,j,ierr,bin
        character(len=clen) :: varname

        Ysmin = zero; Ysmax = one; dYs = (Ysmax - Ysmin) / (nbins - 1)

        ! Get FFT of input
        call transpose_y_to_z(input, zbuffer, mir%gp)
        call fftz%fftz(zbuffer, zfft)

        ! Transpose Ys to Z to average the FFT
        call transpose_y_to_z( Ys, zbuffer, mir%gp)
        Ysavg = sum(zbuffer,3) / mir%gp%zsz(3) ! Ys average
        
        Ysbin = zero
        
        do bin = 1,nbins
            proc_output = zero
            proc_Yssum = zero

            do j = 1,mir%gp%zsz(2)
                do i = 1,mir%gp%zsz(1)
                    if ( (Ysavg(i,j) .LE. Ysbin+half*dYs) .AND. (Ysavg(i,j) .GE. Ysbin-half*dYs) ) then
                        proc_output  = proc_output + zfft(i,j,:)
                        proc_Yssum = proc_Yssum + one
                    end if
                end do
            end do

            call MPI_ALLREDUCE(proc_output, output, mir%gp%zsz(3)/2+1, mpickind, MPI_SUM, MPI_COMM_WORLD, ierr) ! Add proc_output  across processors
            call MPI_ALLREDUCE(proc_Yssum,   Yssum,                 1, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr) ! Add proc_Yssum across processors

            output = ( output / Yssum ) / mir%nz ! Divide by nz to get physical FFT

            write(varname,'(A,F3.1)') "TKE_Ys", Ysbin
            call write_spectrum(step, varname, output)

            Ysbin = Ysbin + dYs
        end do

    end subroutine

    subroutine get_pdf(input,nbins,bins,min_input,max_input,pdf,mask)
        use reductions, only: P_MAXVAL, P_MINVAL, P_SUM

        real(rkind), dimension(:,:,:), intent(in) :: input
        integer, intent(in) :: nbins
        real(rkind), intent(in) :: max_input, min_input
        real(rkind), dimension(:,:,:), intent(in) :: mask  ! This should only have zeros and ones
        real(rkind), dimension(nbins+1), intent(out) :: bins
        real(rkind), dimension(nbins), intent(out) :: pdf

        real(rkind), dimension(nbins) :: proc_count
        real(rkind) ::binsize

        integer :: i, j, k, nb, ierr

        ! Create equispaced bins between maximum and minimum
        binsize = (max_input - min_input)/real(nbins,rkind)
        bins(1) = min_input
        do nb = 2,nbins+1
            bins(nb) = bins(nb-1) + binsize
        end do

        proc_count = zero
        ! Get the counts per bin for this proc
        do k=1,mir%gp%ysz(3)
            do j=1,mir%gp%ysz(2)
                do i=1,mir%gp%ysz(1)
                    if ( mask(i,j,k) .LT. half ) cycle ! No need to go to the bin loop if mask is zero
                    do nb=1,nbins
                        if ( ( input(i,j,k) .GE. bins(nb) ) .AND. ( input(i,j,k) .LE. bins(nb+1) ) ) then
                            proc_count(nb) = proc_count(nb) + one
                            exit ! No need to go over the other bins
                        end if
                    end do
                end do
            end do
        end do

        call MPI_ALLREDUCE(proc_count, pdf, nbins, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr) ! Now all procs will have the total counts per bin

        ! Now normalize the pdf to make sum of probabilities equal to unity
        pdf = pdf / (binsize*sum(pdf))

    end subroutine

    subroutine write_pdf(step, varname, bins, pdf)
        integer, intent(in) :: step
        character(len=*), intent(in) :: varname
        real(rkind), dimension(:), intent(in) :: bins, pdf

        character(len=clen) :: filename
        integer :: i,n,iounit=51

        ! Assume all processors have same pdf array, so only master proc can
        ! write out the file
        n = size(pdf)

        write(filename,'(4A,I4.4,A)') trim(outputdir), "/pdf_", trim(varname), "_", step, ".dat"
        call message(1,"Writing " // trim(varname) // " pdf to " // trim(filename))
        if (nrank == 0) then
            open(unit=iounit, file=trim(filename), form='FORMATTED', status="REPLACE")

            write(iounit,'(2A)') "# pdf for ", trim(varname)
            write(iounit,'(2A26)') "# bin", "PDF"
            do i=1,n
                write(iounit,'(2ES26.16)') half*(bins(i)+bins(i+1)), pdf(i)
            end do

            close(iounit)
        end if

    end subroutine

    subroutine write_tpc(step, varname, tpc)
        integer, intent(in) :: step
        character(len=*), intent(in) :: varname
        real(rkind), dimension(:), intent(in) :: tpc

        character(len=clen) :: filename
        integer :: i,n,iounit=51

        ! Assume all processors have same pdf array, so only master proc can
        ! write out the file
        n = size(tpc)

        write(filename,'(4A,I4.4,A)') trim(outputdir), "/tpc_", trim(varname), "_", step, ".dat"
        call message(1,"Writing " // trim(varname) // " two-point correlation to " // trim(filename))
        if (nrank == 0) then
            open(unit=iounit, file=trim(filename), form='FORMATTED', status="REPLACE")

            write(iounit,'(2A)') "# Two-point correlation for ", trim(varname)
            write(iounit,'(A,ES26.16)') "# dx ", mir%dx
            do i=1,n
                write(iounit,'(ES26.16)') tpc(i)
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

    SUBROUTINE sgs_diffusivity(Y,Diff)
     real(rkind), DIMENSION(:,:,:), INTENT(IN) :: Y
     real(rkind), DIMENSION(:,:,:), INTENT(OUT) :: Diff
     real(rkind), DIMENSION(SIZE(Y,1),SIZE(Y,2),SIZE(Y,3)) :: Lart,Lden,art,tmp,dum,Dsgs,Ysgs

     real(rkind), dimension(mir%gp%xsz(1),mir%gp%xsz(2),mir%gp%xsz(3)) :: xbuf1, xbuf2
     real(rkind), dimension(mir%gp%zsz(1),mir%gp%zsz(2),mir%gp%zsz(3)) :: zbuf1, zbuf2

     Diff = zero
   ! Oscillation control (Dsgs)========================================================================
        Lart = SQRT(mir%dx**2 + mir%dy**2 + mir%dz**2) ! grid scale
   ! Kawai's formula ----------------------------------
          dum = one
          call transpose_y_to_x( Y, xbuf1, mir%gp )
          CALL der%d2dx2(xbuf1,xbuf2)
          CALL der%d2dx2(xbuf2,xbuf1)
          call transpose_x_to_y( xbuf1, tmp, mir%gp )
          Dsgs = dum**4 * tmp * mir%dx**5
          Lart = ABS(tmp)*mir%dx     ! Numerator for Y weighted grid length
          Lden = tmp**2              ! Denominator for Y weighted grid length
   
          dum = one
          CALL der%d2dy2(Y,art)
          CALL der%d2dy2(art,tmp)
          Dsgs = Dsgs + dum**4 * tmp * mir%dy**5
          Lart = Lart + ABS(tmp)*mir%dy
          Lden = Lden + tmp**2
   
          dum = one
          call transpose_y_to_z( Y, zbuf1, mir%gp )
          CALL der%d2dz2(zbuf1,zbuf2)
          CALL der%d2dz2(zbuf2,zbuf1)
          call transpose_z_to_y( zbuf1, tmp, mir%gp )
          Dsgs = Dsgs + dum**4 * tmp * mir%dz**5
          Lart = Lart + ABS(tmp)*mir%dz
          Lden = Lden + tmp**2
   
          Dsgs = Cd*mir%c*Dsgs
          Lart = Lart/( SQRT(Lden) + eps)
   ! Kawai's formula ----------------------------------
   
   ! Bounds control (Ysgs)=============================================================================
          Ysgs = CY*mir%c*0.5*(ABS(Y)  -one+ABS(one-Y  ))*Lart   !    0<=Y<=1
   ! ceiling + filter==================================================================================
        art = MAX(Dsgs,Ysgs)
        tmp = art
        ! CALL filter(gfilter,tmp,art) 
        call transpose_y_to_x( tmp, xbuf1, mir%gp )
        call gfil%filterx(xbuf1,xbuf2)
        call transpose_x_to_y( xbuf2, art, mir%gp )
        call gfil%filtery( art, tmp )
        call transpose_y_to_z( tmp, zbuf1, mir%gp )
        call gfil%filterz(zbuf1,zbuf2)
        call transpose_z_to_y( zbuf2, art, mir%gp )
        Diff = MAX(Diff,art)                   ! cm^2/s
    END SUBROUTINE sgs_diffusivity

    subroutine favre(f,ff,rho_avg)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)), intent(in)  :: f
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)), intent(out) :: ff
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)), intent(in), optional  :: rho_avg
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)) :: tmp

        call P_AVGZ( mir%gp, mir%rho*f, ff )
        if ( present(rho_avg) ) then
            ff = ff / rho_avg
        else
            call P_AVGZ( mir%gp, mir%rho, tmp )
            ff = ff / tmp
        end if
    end subroutine

    subroutine invariants(R11,R12,R13,R22,R23,R33,eta,xi)
        ! Get the invariants of the Reynolds Stress anisotropy tensor
        use constants, only: six
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)), intent(in)  :: R11,R12,R13,R22,R23,R33
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)), intent(out) :: eta, xi

        real(rkind), dimension(3,3) :: b, b2, b3
        real(rkind) :: rtmp,rmax
        integer :: i,j

        real(rkind), parameter :: inv_cutoff = real(1.e-2,rkind)

        ! Find global max of R11, R22 and R33
        rmax = P_MAXVAL(R11)
        rtmp = P_MAXVAL(R22)
        if ( rtmp .GT. rmax ) rmax = rtmp
        rtmp = P_MAXVAL(R33)
        if ( rtmp .GT. rmax ) rmax = rtmp

        do j = 1,mir%gp%ysz(2)
            do i = 1,mir%gp%ysz(1)
                rtmp = R11(i,j) + R22(i,j) + R33(i,j)
                
                if ( rtmp/rmax .GE. inv_cutoff ) then
                    ! Anisotropy tensor
                    b(1,1) = R11(i,j)/rtmp - one/three
                    b(1,2) = R12(i,j)/rtmp
                    b(1,3) = R13(i,j)/rtmp
                    b(2,1) = b(1,2)
                    b(2,2) = R22(i,j)/rtmp - one/three
                    b(2,3) = R23(i,j)/rtmp
                    b(3,1) = b(1,3)
                    b(3,2) = b(3,2)
                    b(3,3) = R33(i,j)/rtmp - one/three
                else
                    ! If below the noise threshold, set equal to isotropic state
                    b = zero
                end if

                b2 = matmul(b,b)
                b3 = matmul(b,b2)

                eta(i,j) = sqrt( ( b2(1,1) + b2(2,2) + b2(3,3) )/six ) ! sqrt( tr(b^2)/6 )
                xi (i,j) = ( ( b3(1,1) + b3(2,2) + b3(3,3) )/six )**(one/three) ! cbrt( tr(b^3)/6 )

            end do
        end do

    end subroutine

    subroutine two_point_correlation(f,g,tpc,region)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)), intent(in)  :: f, g, region
        real(rkind), dimension(mir%nz/2 + 1), intent(out) :: tpc

        real(rkind), dimension(mir%gp%zsz(1), mir%gp%zsz(2), mir%gp%zsz(3)) :: fz, gz, regionz
        real(rkind), dimension(mir%nz/2 + 1) :: tpc_proc

        integer :: i, j, k1, k2, dk, nregion, ierr

        ! Transpose everything to Z before doing two-point correlation stuff
        call transpose_y_to_z(f, fz, mir%gp)
        call transpose_y_to_z(g, gz, mir%gp)
        call transpose_y_to_z(region, regionz, mir%gp)

        if (.NOT. mir%periodicz) then
            call warning("WARNING: two_point_correlation assumes periodicity in Z, but that's not the case here.")
        end if

        nregion = 0
        tpc_proc = zero

        do j = 1,mir%gp%zsz(2)
            do i = 1,mir%gp%zsz(1)

                ! Compute two-point correlation only for (x,y) in region
                if (regionz(i,j,1) .GT. zero) then
                    nregion = nregion + 1

                    do dk = 0,mir%nz/2
                        do k1 = 1,mir%nz
                            k2 = k1 + dk
                            if (k2 .GT. mir%nz) then
                                k2 = k2 - mir%nz ! Assuming periodicity in Z
                            end if
                            tpc_proc(dk+1) = tpc_proc(dk+1) + fz(i,j,k1)*gz(i,j,k2)/real(mir%nz,rkind)
                        end do
                    end do

                end if

            end do
        end do

        if (nregion .GT. 0) then
            tpc_proc = tpc_proc / real(nregion,rkind)
        end if
       
        ! Reduce tpc from all processors to one mean tpc
        call MPI_Allreduce(tpc_proc, tpc, mir%nz/2 + 1, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
        tpc = tpc / real(nproc,rkind)

    end subroutine

    subroutine get_duidxj(u,v,w,duidxj)
        use operators, only: gradient
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)),            intent(in)  :: u, v, w
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 9), target, intent(out) :: duidxj

        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

        call gradient(mir%gp,der,u,dudx,dudy,dudz)
        call gradient(mir%gp,der,v,dvdx,dvdy,dvdz)
        call gradient(mir%gp,der,w,dwdx,dwdy,dwdz)

    end subroutine

    subroutine get_duidxj2D(u,v,w,duidxj2D)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)),            intent(in)  :: u, v, w
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), 9), target, intent(out) :: duidxj2D

        real(rkind), dimension(mir%gp%xsz(1),mir%gp%xsz(2),1) :: xbuf1, xbuf2
        real(rkind), dimension(:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        dudx => duidxj2D(:,:,1); dudy => duidxj2D(:,:,2); dudz => duidxj2D(:,:,3);
        dvdx => duidxj2D(:,:,4); dvdy => duidxj2D(:,:,5); dvdz => duidxj2D(:,:,6);
        dwdx => duidxj2D(:,:,7); dwdy => duidxj2D(:,:,8); dwdz => duidxj2D(:,:,9);

        ! Set Z derivatives to zero
        dudz = zero; dvdz = zero; dwdz = zero

        ! dudx
        call gp2D%transpose_y_to_x(u,xbuf1)
        call der2D%ddx(xbuf1,xbuf2)
        call gp2D%transpose_x_to_y(xbuf2,dudx)
        
        ! dvdx
        call gp2D%transpose_y_to_x(v,xbuf1)
        call der2D%ddx(xbuf1,xbuf2)
        call gp2D%transpose_x_to_y(xbuf2,dvdx)
        
        ! dwdx
        call gp2D%transpose_y_to_x(w,xbuf1)
        call der2D%ddx(xbuf1,xbuf2)
        call gp2D%transpose_x_to_y(xbuf2,dwdx)

        call der2D%ddy(u,dudy)
        call der2D%ddy(v,dvdy)
        call der2D%ddy(w,dwdy)
        
    end subroutine

    subroutine get_gradp2D(p_avg,gradp2D)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)),            intent(in)  :: p_avg
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), 3), target, intent(out) :: gradp2D

        real(rkind), dimension(mir%gp%xsz(1),mir%gp%xsz(2),1) :: xbuf1, xbuf2
        real(rkind), dimension(:,:), pointer :: dpdx, dpdy, dpdz

        dpdx => gradp2D(:,:,1); dpdy => gradp2D(:,:,2); dpdz => gradp2D(:,:,3);

        ! Set Z derivatives to zero
        dpdz = zero

        ! dpdx
        call gp2D%transpose_y_to_x(p_avg,xbuf1)
        call der2D%ddx(xbuf1,xbuf2)
        call gp2D%transpose_x_to_y(xbuf2,dpdx)
        
        call der2D%ddy(p_avg,dpdy)
        
    end subroutine

    subroutine get_tauij(duidxj,mu,beta,tauij)
        use constants, only: two, three, four
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 9), target, intent(in)  :: duidxj
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)   ),         intent(in)  :: mu, beta
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 6), target, intent(out) :: tauij
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        tauxx => tauij(:,:,:,1); tauxy => tauij(:,:,:,2); tauxz => tauij(:,:,:,3);
                                 tauyy => tauij(:,:,:,4); tauyz => tauij(:,:,:,5);
                                                          tauzz => tauij(:,:,:,6);

        tauxx = ( (four/three)*mu + beta )*dudx + (beta - (two/three)*mu)*(dvdy+dwdz)
        tauxy = mu * (dudy + dvdx)
        tauxz = mu * (dudz + dwdx)
        tauyy = ( (four/three)*mu + beta )*dvdy + (beta - (two/three)*mu)*(dudx+dwdz)
        tauyz = mu * (dvdz + dwdy)
        tauzz = ( (four/three)*mu + beta )*dwdz + (beta - (two/three)*mu)*(dudx+dvdy)
    end subroutine

    subroutine get_tauij_avg(tauij, tauij_avg)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 6), intent(in)  :: tauij
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2),                6), intent(out) :: tauij_avg

        integer :: i

        do i=1,6
            call P_AVGZ( mir%gp, tauij(:,:,:,i), tauij_avg(:,:,i) )
        end do
    end subroutine

    subroutine get_dissipation(upp,vpp,wpp,tauij,dissipation,pdil_fluct, tke_visc_transport)
        use operators, only: gradient, divergence
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)   ),         intent(in)  :: upp,vpp,wpp
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 6), target, intent(in)  :: tauij
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)),                           intent(out) :: dissipation, pdil_fluct, tke_visc_transport
        
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) :: ytmp1,ytmp2,ytmp3,diss,pdil

        tauxx => tauij(:,:,:,1); tauxy => tauij(:,:,:,2); tauxz => tauij(:,:,:,3);
                                 tauyy => tauij(:,:,:,4); tauyz => tauij(:,:,:,5);
                                                          tauzz => tauij(:,:,:,6);

        ytmp1 = tauxx*upp + tauxy*vpp + tauxz*wpp
        ytmp2 = tauxy*upp + tauyy*vpp + tauyz*wpp
        ytmp3 = tauxz*upp + tauyz*vpp + tauzz*wpp
        call divergence( mir%gp, der, ytmp1, ytmp2, ytmp3, diss )
        call P_AVGZ( mir%gp, diss, tke_visc_transport )

        ! Dissipation = <tau_ij dupp_i/dx_j>
        call gradient(mir%gp,der,upp,ytmp1,ytmp2,ytmp3)
        diss = tauxx*ytmp1 + tauxy*ytmp2 + tauxz*ytmp3
        pdil = ytmp1

        call gradient(mir%gp,der,vpp,ytmp1,ytmp2,ytmp3)
        diss = diss + tauxy*ytmp1 + tauyy*ytmp2 + tauyz*ytmp3
        pdil = pdil + ytmp2

        call gradient(mir%gp,der,wpp,ytmp1,ytmp2,ytmp3)
        diss = diss + tauxz*ytmp1 + tauyz*ytmp2 + tauzz*ytmp3
        pdil = pdil + ytmp3

        call P_AVGZ( mir%gp, diss, dissipation )

        call P_AVGZ( mir%gp, mir%p*pdil, pdil_fluct )

    end subroutine

    subroutine get_production(R11,R12,R13,R22,R23,R33,duidxj2D,rho_avg,production)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)   ),         intent(in)  :: R11,R12,R13,R22,R23,R33,rho_avg
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), 9), target, intent(in)  :: duidxj2D
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)   ),         intent(out) :: production
        
        real(rkind), dimension(:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        dudx => duidxj2D(:,:,1); dudy => duidxj2D(:,:,2); dudz => duidxj2D(:,:,3);
        dvdx => duidxj2D(:,:,4); dvdy => duidxj2D(:,:,5); dvdz => duidxj2D(:,:,6);
        dwdx => duidxj2D(:,:,7); dwdy => duidxj2D(:,:,8); dwdz => duidxj2D(:,:,9);

        production = R11*dudx + R12*dudy + R13*dudz &
                   + R12*dvdx + R22*dvdy + R23*dvdz &
                   + R13*dwdx + R23*dwdy + R33*dwdz

        production = -rho_avg*production

    end subroutine

    subroutine get_meandiss(tauij_avg,duidxj2D,meandiss)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), 6), target, intent(in)  :: tauij_avg
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), 9), target, intent(in)  :: duidxj2D
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)   ),         intent(out) :: meandiss
        
        real(rkind), dimension(:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz

        dudx => duidxj2D(:,:,1); dudy => duidxj2D(:,:,2); dudz => duidxj2D(:,:,3);
        dvdx => duidxj2D(:,:,4); dvdy => duidxj2D(:,:,5); dvdz => duidxj2D(:,:,6);
        dwdx => duidxj2D(:,:,7); dwdy => duidxj2D(:,:,8); dwdz => duidxj2D(:,:,9);
        
        tauxx => tauij_avg(:,:,1); tauxy => tauij_avg(:,:,2); tauxz => tauij_avg(:,:,3);
                                   tauyy => tauij_avg(:,:,4); tauyz => tauij_avg(:,:,5);
                                                              tauzz => tauij_avg(:,:,6);

        meandiss = tauxx*dudx + tauxy*dudy + tauxz*dudz &
                 + tauxy*dvdx + tauyy*dvdy + tauyz*dvdz &
                 + tauxz*dwdx + tauyz*dwdy + tauzz*dwdz

    end subroutine

    subroutine get_mix_dilatation(Diff,mix_pdil,mix_pdil_fluct)
        use operators, only: gradient, divergence
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 2), intent(in)  :: Diff
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)                  ), intent(out) :: mix_pdil,mix_pdil_fluct

        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)   )  :: mix_dilatation
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 2)  :: art_Diff
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)   )  :: ytmp1, ytmp2, ytmp3
        
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)   ) :: p_avg, dil_avg
       
        integer :: i

        call sgs_diffusivity(mir%Ys(:,:,:,1),art_Diff(:,:,:,1))
        call sgs_diffusivity(mir%Ys(:,:,:,2),art_Diff(:,:,:,2))

        art_Diff = art_Diff + Diff
        art_Diff(:,:,:,1) = mir%Ys(:,:,:,2)*art_Diff(:,:,:,1) + mir%Ys(:,:,:,1)*art_Diff(:,:,:,2)  ! Put D_mix in art_Diff(:,:,:,1)

        call gradient(mir%gp,der,mir%rho,ytmp1,ytmp2,ytmp3)
        ytmp1 = -art_Diff(:,:,:,1) * ytmp1 / mir%rho ! - (D_mix / rho) * grad(rho)
        ytmp2 = -art_Diff(:,:,:,1) * ytmp2 / mir%rho ! - (D_mix / rho) * grad(rho)
        ytmp3 = -art_Diff(:,:,:,1) * ytmp3 / mir%rho ! - (D_mix / rho) * grad(rho)
        
        call divergence( mir%gp, der, ytmp1, ytmp2, ytmp3, mix_dilatation ) ! -div( (D_mix/rho) * grad(rho) )
        
        ! Get average pressure and dilatation
        call P_AVGZ( mir%gp, mir%p, p_avg )
        call P_AVGZ( mir%gp, mix_dilatation, dil_avg )

        mix_pdil = p_avg * dil_avg

        ! Get fluctuation dilatation
        do i = 1,mir%gp%ysz(3)
            mix_dilatation(:,:,i) = mix_dilatation(:,:,i) - dil_avg ! Get fluctuation dilatation
        end do

        call P_AVGZ( mir%gp, mir%p * mix_dilatation, mix_pdil_fluct )

    end subroutine

    subroutine get_PS_diss(upp,vpp,wpp,p,p_avg,tauij,rho,rho_avg,PSij,diss)
        use operators, only: gradient, divergence
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)   ),         intent(in)  :: upp,vpp,wpp,p,rho
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 6),         intent(in)  :: tauij
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2)),                           intent(in)  :: p_avg,rho_avg
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), 6), target,                intent(out) :: PSij, diss
        
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) :: p_prime, ytmp1
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) :: tauxx_pp, tauxy_pp, tauxz_pp, tauyy_pp, tauyz_pp, tauzz_pp
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 6), target :: duidxj
        real(rkind), dimension(:,:), pointer :: PSxx, PSxy, PSxz, PSyy, PSyz, PSzz
        real(rkind), dimension(:,:), pointer :: Dxx, Dxy, Dxz, Dyy, Dyz, Dzz
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        integer :: i

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

        PSxx => PSij(:,:,1); PSxy => PSij(:,:,2); PSxz => PSij(:,:,3);
                             PSyy => PSij(:,:,4); PSyz => PSij(:,:,5);
                                                  PSzz => PSij(:,:,6);

        Dxx => diss(:,:,1); Dxy => diss(:,:,2); Dxz => diss(:,:,3);
                            Dyy => diss(:,:,4); Dyz => diss(:,:,5);
                                                Dzz => diss(:,:,6);

        call get_duidxj(upp,vpp,wpp,duidxj)

        do i = 1,mir%gp%ysz(3)
            p_prime(:,:,i) = p(:,:,i) - p_avg
        end do

        ytmp1 = p_prime * ( dudx + dudx )
        call P_AVGZ( mir%gp, ytmp1, PSxx )

        ytmp1 = p_prime * ( dudy + dvdx )
        call P_AVGZ( mir%gp, ytmp1, PSxy )

        ytmp1 = p_prime * ( dudz + dwdx )
        call P_AVGZ( mir%gp, ytmp1, PSxz )

        ytmp1 = p_prime * ( dvdy + dvdy )
        call P_AVGZ( mir%gp, ytmp1, PSyy )

        ytmp1 = p_prime * ( dvdz + dwdy )
        call P_AVGZ( mir%gp, ytmp1, PSyz )

        ytmp1 = p_prime * ( dwdz + dwdz )
        call P_AVGZ( mir%gp, ytmp1, PSzz )

        ! Get favre( tauij )
        call P_AVGZ( mir%gp, rho*tauij(:,:,:,1), Dxx )
        call P_AVGZ( mir%gp, rho*tauij(:,:,:,2), Dxy )
        call P_AVGZ( mir%gp, rho*tauij(:,:,:,3), Dxz )
        call P_AVGZ( mir%gp, rho*tauij(:,:,:,4), Dyy )
        call P_AVGZ( mir%gp, rho*tauij(:,:,:,5), Dyz )
        call P_AVGZ( mir%gp, rho*tauij(:,:,:,6), Dzz )

        do i = 1,mir%gp%ysz(3)
            tauxx_pp(:,:,i) = tauij(:,:,i,1) - Dxx/rho_avg
            tauxy_pp(:,:,i) = tauij(:,:,i,2) - Dxy/rho_avg
            tauxz_pp(:,:,i) = tauij(:,:,i,3) - Dxz/rho_avg
            tauyy_pp(:,:,i) = tauij(:,:,i,4) - Dyy/rho_avg
            tauyz_pp(:,:,i) = tauij(:,:,i,5) - Dyz/rho_avg
            tauzz_pp(:,:,i) = tauij(:,:,i,6) - Dzz/rho_avg
        end do

        ! Dxx
        ytmp1 = -(tauxx_pp*dudx + tauxy_pp*dudy + tauxz_pp*dudz) - (tauxx_pp*dudx + tauxy_pp*dudy + tauxz_pp*dudz)
        call P_AVGZ( mir%gp, ytmp1, Dxx )

        ! Dxy
        ytmp1 = -(tauxy_pp*dudx + tauyy_pp*dudy + tauyz_pp*dudz) - (tauxx_pp*dvdx + tauxy_pp*dvdy + tauxz_pp*dvdz)
        call P_AVGZ( mir%gp, ytmp1, Dxy )

        ! Dxz
        ytmp1 = -(tauxz_pp*dudx + tauyz_pp*dudy + tauzz_pp*dudz) - (tauxx_pp*dwdx + tauxy_pp*dwdy + tauxz_pp*dwdz)
        call P_AVGZ( mir%gp, ytmp1, Dxz )

        ! Dyy
        ytmp1 = -(tauxy_pp*dvdx + tauyy_pp*dvdy + tauyz_pp*dvdz) - (tauxy_pp*dvdx + tauyy_pp*dvdy + tauyz_pp*dvdz)
        call P_AVGZ( mir%gp, ytmp1, Dyy )

        ! Dyz
        ytmp1 = -(tauxz_pp*dvdx + tauyz_pp*dvdy + tauzz_pp*dvdz) - (tauxy_pp*dwdx + tauyy_pp*dwdy + tauyz_pp*dwdz)
        call P_AVGZ( mir%gp, ytmp1, Dyz )

        ! Dzz
        ytmp1 = -(tauxz_pp*dwdx + tauyz_pp*dwdy + tauzz_pp*dwdz) - (tauxz_pp*dwdx + tauyz_pp*dwdy + tauzz_pp*dwdz)
        call P_AVGZ( mir%gp, ytmp1, Dzz )

    end subroutine

end module

program IRM_vorticity
    use mpi
    use kind_parameters, only: rkind, clen
    use constants,       only: zero, half, one, four, eps
    use miranda_tools,   only: miranda_reader
    use io_VTK_stuff,    only: io_VTK
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use gridtools,       only: alloc_buffs
    use operators,       only: curl, gradient, divergence
    use reductions,      only: P_AVGZ, P_MAXVAL, P_MINVAL, P_SUM
    use timer,           only: tic, toc
    use exits,           only: message, GracefulExit

    use IRM_vorticity_mod

    implicit none

    character(len=clen) :: inputfile

    character(len=clen), dimension(:), allocatable :: varnames

    real(rkind) :: cutoff = 0.5_rkind  ! Cutoff for interface region

    real(rkind), dimension(:,:,:,:), allocatable, target :: buffer
    real(rkind), dimension(:,:,:,:), pointer :: vort, Diff, Xs
    real(rkind), dimension(:,:,:), pointer :: upp, vpp, wpp, TKE, div, region, chi, mu, ktc, Rij

    real(rkind), dimension(:,:,:), allocatable, target :: buffer2d
    real(rkind), dimension(:,:,:), pointer :: vort_avgz
    real(rkind), dimension(:,:), pointer :: rho_avg, u_avg, v_avg, w_avg, CO2_avg, TKE_avg, div_avg, chi_avg, MMF_avg, density_self_correlation
    real(rkind), dimension(:,:), pointer :: R11, R12, R13, R22, R23, R33
    real(rkind), dimension(:,:), pointer :: a_x, a_y, a_z, rhop_sq, eta, xi
    real(rkind), dimension(:,:), pointer :: p_avg, pdil, meandiss, production, dissipation, pdil_fluct, tke_visc_transport, mix_pdil, mix_pdil_fluct

    real(rkind), dimension(:,:,:,:), allocatable :: duidxj, tauij
    real(rkind), dimension(:,:,:),   allocatable :: duidxj2D, tauij_avg, gradp2D, PSij, diss_tensor

    real(rkind), dimension(:,:,:), allocatable :: xtmp1, xtmp2, ytmp1

    real(rkind), dimension(:), allocatable :: Y_air_x, Y_CO2_x

    integer :: nbins = 64
    real(rkind), dimension(:), allocatable :: bins, pdf, tpc

    complex(rkind), dimension(:), allocatable :: TKEfft

    real(rkind) :: MMF_int, chi_int, chi_art, vortz_int, vortz_pos, vortz_neg, TKE_int, mwidth, rhop_sq_int

    logical :: fftw_exhaustive = .FALSE.

    character(len=clen) :: time_message

    character(len=clen) :: outputfile
    integer :: iounit = 93

    integer :: nVTKvars
    integer :: ierr, step, i, j

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
    call der2D%init(                                    mir%gp, &
                          mir%dx,        mir%dy,        mir%dz, &
                   mir%periodicx, mir%periodicy, mir%periodicz, &
                    derivative_x,  derivative_y,  derivative_z  )
    call der2D%set_xsz( [mir%gp%xsz(1), mir%gp%xsz(2), 1] )                 ! Set xsz to make arrays 2D
    call der2D%set_ysz( [mir%gp%ysz(1), mir%gp%ysz(2), 1] )                 ! Set ysz to make arrays 2D

    ! Split MPI_COMM_WORLD into individual XY plane communicators
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, mir%gp%yst(3), nrank, XY_COMM, ierr)
    call MPI_COMM_RANK(XY_COMM, xyrank, ierr)
    call MPI_COMM_SIZE(XY_COMM, xyproc, ierr)

    ! Initialize gp2D
    call gp2D%init(mir%nx,mir%ny,XY_COMM)

    call message("Initializing the gaussian filter object" )
    call gfil%init(                                     mir%gp, &
                   mir%periodicx, mir%periodicy, mir%periodicz, &
                      'gaussian',    'gaussian',    'gaussian'  )
    
    ! Initialize the FFT object (in the Z decomposition since that is where FFT will be performed)
    call message("Initializing the FFT object" )
    ierr = fftz%init( mir%nz, 'z', mir%gp%zsz(1), mir%gp%zsz(2), mir%dz, fftw_exhaustive )

    ! Allocate 3D buffer and associate pointers for convenience
    call alloc_buffs(buffer, 13+2*mir%ns, 'y', mir%gp)
    i = 1
    vort   => buffer(:,:,:,i:i+2);        i = i+3
    upp    => buffer(:,:,:,i);            i = i+1
    vpp    => buffer(:,:,:,i);            i = i+1
    wpp    => buffer(:,:,:,i);            i = i+1
    TKE    => buffer(:,:,:,i);            i = i+1
    div    => buffer(:,:,:,i);            i = i+1
    region => buffer(:,:,:,i);            i = i+1
    chi    => buffer(:,:,:,i);            i = i+1
    mu     => buffer(:,:,:,i);            i = i+1
    ktc    => buffer(:,:,:,i);            i = i+1
    Diff   => buffer(:,:,:,i:i+mir%ns-1); i = i+mir%ns
    Xs     => buffer(:,:,:,i:i+mir%ns-1); i = i+mir%ns
    Rij    => buffer(:,:,:,i)

    ! Allocate 2D buffer and associate pointers for convenience
    allocate( buffer2d( mir%gp%ysz(1), mir%gp%ysz(2), 34 ) )
    vort_avgz                => buffer2d(:,:,1:3)
    rho_avg                  => buffer2d(:,:,4 )
    u_avg                    => buffer2d(:,:,5 )
    v_avg                    => buffer2d(:,:,6 )
    w_avg                    => buffer2d(:,:,7 )
    CO2_avg                  => buffer2d(:,:,8 )
    TKE_avg                  => buffer2d(:,:,9 )
    div_avg                  => buffer2d(:,:,10)
    chi_avg                  => buffer2d(:,:,11)
    MMF_avg                  => buffer2d(:,:,12)
    density_self_correlation => buffer2d(:,:,13)
    R11                      => buffer2d(:,:,14)
    R12                      => buffer2d(:,:,15)
    R13                      => buffer2d(:,:,16)
    R22                      => buffer2d(:,:,17)
    R23                      => buffer2d(:,:,18)
    R33                      => buffer2d(:,:,19)
    a_x                      => buffer2d(:,:,20)
    a_y                      => buffer2d(:,:,21)
    a_z                      => buffer2d(:,:,22)
    rhop_sq                  => buffer2d(:,:,23)
    eta                      => buffer2d(:,:,24)
    xi                       => buffer2d(:,:,25)
    p_avg                    => buffer2d(:,:,26)
    pdil                     => buffer2d(:,:,27)
    meandiss                 => buffer2d(:,:,28)
    production               => buffer2d(:,:,29)
    dissipation              => buffer2d(:,:,30)
    pdil_fluct               => buffer2d(:,:,31)
    tke_visc_transport       => buffer2d(:,:,32)
    mix_pdil                 => buffer2d(:,:,33)
    mix_pdil_fluct           => buffer2d(:,:,34)
    
    allocate( xtmp1( mir%gp%xsz(1), mir%gp%xsz(2), 1 ) )
    allocate( xtmp2( mir%gp%xsz(1), mir%gp%xsz(2), 1 ) )
    
    allocate( ytmp1( mir%gp%ysz(1), mir%gp%ysz(2), 1 ) )
    
    allocate( TKEfft(mir%nz/2+1) )
    allocate( tpc   (mir%nz/2+1) )

    allocate( Y_air_x( mir%gp%ysz(1) ) )
    allocate( Y_CO2_x( mir%gp%ysz(1) ) )

    allocate( bins(nbins+1) )
    allocate( pdf (nbins  ) )
    
    allocate( duidxj(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 9) )
    allocate( tauij (mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 6) )
    allocate( duidxj2D (mir%gp%ysz(1), mir%gp%ysz(2), 9) )
    allocate( tauij_avg(mir%gp%ysz(1), mir%gp%ysz(2), 6) )
    allocate( gradp2D (mir%gp%ysz(1), mir%gp%ysz(2), 3) )
    allocate( PSij(mir%gp%ysz(1), mir%gp%ysz(2), 6) )
    allocate( diss_tensor(mir%gp%ysz(1), mir%gp%ysz(2), 6) )


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
        write(iounit,'(A,A25,9A26)') "#", "Time", "Mixed width", "Z vorticity", "+ve Z vorticity", "-ve Z vorticity", "TKE", "Scalar dissipation", "Art sca dissipation", "MMF", "< rho' rho' >"
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
        
        ! Get Y_CO2 average
        call P_AVGZ( mir%gp, mir%Ys(:,:,:,2), CO2_avg )

        ! Get the interface region
        region = zero
        do j=1,mir%gp%ysz(2)
            do i=1,mir%gp%ysz(1)
                if ( four * (one-CO2_avg(i,j)) * CO2_avg(i,j) .GT. cutoff ) then
                    region(i,j,:) = one
                end if
            end do
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
        
        ! Get integrated Z vorticity
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

        ! Get < rho' rho' >
        do i = 1,mir%gp%ysz(3)
            chi(:,:,i) = ( mir%rho(:,:,i) - rho_avg )**2 ! Put (rho' rho') in chi temporarily
        end do
        call P_AVGZ( mir%gp, chi, rhop_sq )
        rhop_sq_int = P_SUM(chi*region)/P_SUM(region)

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

            ! Get integrated chi
            chi_int = P_SUM(chi)*mir%dx*mir%dy*mir%dz
            call message(1,"chi_int",chi_int)
            
            call sgs_diffusivity( mir%Ys(:,:,:,2), Diff(:,:,:,2) )
            chi = Diff(:,:,:,2)*( dY_x*dY_x + dY_y*dY_y + dY_z*dY_z ) ! Artificial scalar dissipation
            
            ! Get integrated artificial chi
            chi_art = P_SUM(chi)*mir%dx*mir%dy*mir%dz
            
            ! Get artificial scalar dissipation average
            call P_AVGZ( mir%gp, chi, chi_avg )
        end associate
        
        ! Recompute physical Diff
        call prob_properties(mu,ktc,Diff)
        ! Get dilatation due to mixing
        call get_mix_dilatation(Diff,mix_pdil,mix_pdil_fluct)

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

        ! Get Reynolds stress R11
        Rij = mir%rho*upp*upp
        call P_AVGZ( mir%gp, Rij, R11 )
        R11 = R11 / rho_avg

        ! Get Reynolds stress R12
        Rij = mir%rho*upp*vpp
        call P_AVGZ( mir%gp, Rij, R12 )
        R12 = R12 / rho_avg

        ! Get Reynolds stress R13
        Rij = mir%rho*upp*wpp
        call P_AVGZ( mir%gp, Rij, R13 )
        R13 = R13 / rho_avg

        ! Get Reynolds stress R22
        Rij = mir%rho*vpp*vpp
        call P_AVGZ( mir%gp, Rij, R22 )
        R22 = R22 / rho_avg

        ! Get Reynolds stress R23
        Rij = mir%rho*vpp*wpp
        call P_AVGZ( mir%gp, Rij, R23 )
        R23 = R23 / rho_avg

        ! Get Reynolds stress R33
        Rij = mir%rho*wpp*wpp
        call P_AVGZ( mir%gp, Rij, R33 )
        R33 = R33 / rho_avg

        ! Get Reynolds stress anisotropy invariants
        call  invariants(R11,R12,R13,R22,R23,R33,eta,xi)

        ! Get turbulent mass flux
        call P_AVGZ( mir%gp, upp, a_x )
        a_x = -a_x ! a_x = -<upp>
        call P_AVGZ( mir%gp, vpp, a_y )
        a_y = -a_y ! a_y = -<vpp>
        call P_AVGZ( mir%gp, wpp, a_z )
        a_z = -a_z ! a_z = -<wpp>

        ! Get TKE = (1/2) * (rho*u_i*u_i)/<rho>
        TKE = half * mir%rho * ( upp*upp + vpp*vpp + wpp*wpp )
        ! do i = 1,mir%gp%ysz(3)
        !     TKE(:,:,i) = TKE(:,:,i) / rho_avg
        ! end do
        
        ! Get integrated TKE
        TKE_int = P_SUM(TKE)*mir%dx*mir%dy*mir%dz
        call message(1,"TKE_int",TKE_int)

        ! Get TKE spectrum
        call get_zfft(TKE,TKEfft,region)
        call write_spectrum(step, "TKE", TKEfft)
        call get_conditional_zfft(TKE,mir%Ys(:,:,:,2),11,step)
        
        ! Get TKE average
        call P_AVGZ( mir%gp, TKE, TKE_avg )

        ! Get Divergence
        call divergence( mir%gp, der, mir%u, mir%v, mir%w, div )

        ! Get Divergence average
        call P_AVGZ( mir%gp, div, div_avg )
    
        Rij = one
        ! Get Y_CO2 PDF
        call get_pdf(mir%Ys(:,:,:,2),nbins,bins,0.1_rkind,0.9_rkind,pdf,Rij)
        call write_pdf(step, "CO2", bins, pdf)
    
        call two_point_correlation(wpp,wpp,tpc,region)
        call write_tpc(step, "wpp", tpc)

        !!!! =============================================
        !!!! Energetics ----------------------------------
        
        call get_duidxj(mir%u,mir%v,mir%w,duidxj)
        call get_tauij(duidxj,mir%mu,mir%bulk,tauij)
        call get_tauij_avg(tauij, tauij_avg)

        call get_duidxj2D(u_avg,v_avg,w_avg,duidxj2D)

        ! Pressure-dilatation correlation
        call P_AVGZ( mir%gp, mir%p, p_avg )
        pdil = p_avg*(duidxj2D(:,:,1) + duidxj2D(:,:,5))
        
        call get_gradp2D(p_avg,gradp2D)

        ! Mean-dissipation
        call get_meandiss(tauij_avg,duidxj2D,meandiss)

        ! Shear-production
        call get_production(R11,R12,R13,R22,R23,R33,duidxj2D,rho_avg,production)

        ! Turbulent dissipation and fluctuation pressure-dilatation
        call get_dissipation(upp,vpp,wpp,tauij,dissipation,pdil_fluct, tke_visc_transport)

        ! Get model form of pressure-strain tensor and dissipation tensor
        call get_PS_diss(upp,vpp,wpp,mir%p,p_avg,tauij,mir%rho,rho_avg,PSij,diss_tensor)

        !!!! Energetics ----------------------------------
        !!!! =============================================

        if(nrank == 0) then
            write(iounit,'(10ES26.16)') real(step*1.0d-4,rkind), mwidth, vortz_int, vortz_pos, vortz_neg, TKE_int, chi_int, chi_art, MMF_int, rhop_sq_int
        end if

        ! Write out visualization files
        if ( writeviz ) then
            ! Set vizcount to be same as Miranda step
            call viz%SetVizcount(step)
            call viz%WriteViz(mir%gp, mir%mesh, buffer(:,:,:,1:nVTKvars), secondary=Xs, secondary_names=['volume_fraction_1','volume_fraction_2'])
        end if

        ! Write out 2D postprocessing file
        call write_post2d(step, vort_avgz, TKE_avg, div_avg, rho_avg, chi_avg, MMF_avg, CO2_avg, density_self_correlation, &
                          R11, R12, R13, R22, R23, R33, a_x, a_y, a_z, rhop_sq, eta, xi, &
                          pdil, meandiss, production, dissipation, pdil_fluct, tke_visc_transport, mix_pdil, mix_pdil_fluct, &
                          duidxj2D, gradp2D, PSij, diss_tensor)

        write(time_message,'(A,I4,A)') "Time to postprocess step ", step, " :"
        call toc(trim(time_message))
    end do

    ! Close file for scalar inputs
    if(nrank == 0) close(iounit)

    ! Destroy all variables and exit cleanly
    if (allocated(buffer  )) deallocate( buffer   )
    if (allocated(buffer2d)) deallocate( buffer2d )
    if (allocated( xtmp1  )) deallocate( xtmp1    )
    if (allocated( xtmp2  )) deallocate( xtmp2    )
    if (allocated( ytmp1  )) deallocate( ytmp1    )
    if (allocated( TKEfft )) deallocate( TKEfft   )
    if (allocated( tpc    )) deallocate( tpc      )
    if (allocated( Y_air_x)) deallocate( Y_air_x  )
    if (allocated( Y_CO2_x)) deallocate( Y_CO2_x  )
    if (allocated( bins   )) deallocate( bins     )
    if (allocated( pdf    )) deallocate( pdf      )
    if (allocated( duidxj )) deallocate( duidxj )
    if (allocated( tauij )) deallocate( tauij )
    if (allocated( duidxj2D )) deallocate( duidxj2D )
    if (allocated( tauij_avg )) deallocate( tauij_avg )
    if (allocated( gradp2D )) deallocate( gradp2D )
    if (allocated( PSij )) deallocate( PSij )
    if (allocated( diss_tensor )) deallocate( diss_tensor )

    call der%destroy()
    call gfil%destroy()
    call mir%destroy()
    if ( writeviz ) then
        deallocate( varnames )
        call viz%destroy()
    end if
    call MPI_Finalize(ierr)

end program 
