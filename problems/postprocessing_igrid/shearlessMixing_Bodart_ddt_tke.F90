#include "../incompressible/shearlessMixing_Bodart_files/initialize.F90"
program shearlessMixing_Bodart_ddt_tke
    use kind_parameters,     only: rkind, clen, mpirkind
    use mpi
    use incompressibleGrid, only: igrid
    use basic_io,           only: write_1d_ascii
    use decomp_2d,          only: nrank, transpose_x_to_y, transpose_y_to_z
    use reductions,         only: p_sum, p_maxval, p_minval
    use exits,              only: message_min_max, message
    use timer,              only: tic, toc
    implicit none

    integer :: ierr, tid, tidst, tiden, tstep, ioUnit
    type(igrid), allocatable :: SM
    character(len=clen) :: inputfile, outputdir, fname, tempname
    real(rkind), dimension(:,:,:), allocatable :: u, v, w, T
    real(rkind), dimension(:), allocatable :: dkdt, ddt_MKE, uavg_n, vavg_n, wavg_n, Tavg_n
    real(rkind), dimension(:), allocatable :: dTTdt, ddt_MPE, kn, knp2
    real(rkind), dimension(:), allocatable :: uavg_np2, vavg_np2, wavg_np2, Tavg_np2
    real(rkind) :: dtForce = 1.d-3

    ! Initialize MPI 
    call MPI_Init(ierr)                                               
    
    ! Read input file to get parameters for post processing
    call GETARG(1,inputfile)                                            

    namelist/postProcess/ tidst, tiden, tstep, outputdir, dtForce 
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=postProcess)
    close(ioUnit)    

    allocate(SM)                                                     
    
    call SM%init(inputfile, .true.)                              
    call SM%start_io(.true.)                                          
    call SM%printDivergence()

    ! Allocate memory
    allocate(u(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(v(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(w(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(T(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))

    allocate(dkdt(SM%gpC%zsz(3)))
    allocate(kn(SM%gpC%zsz(3)))
    allocate(knp2(SM%gpC%zsz(3)))
    allocate(ddt_MKE(SM%gpC%zsz(3)))
    allocate(dTTdt(SM%gpC%zsz(3)))
    allocate(ddt_MPE(SM%gpC%zsz(3)))
    allocate(uavg_n(SM%gpC%zsz(3)))
    allocate(vavg_n(SM%gpC%zsz(3)))
    allocate(wavg_n(SM%gpC%zsz(3)))
    allocate(Tavg_n(SM%gpC%zsz(3)))
    allocate(uavg_np2(SM%gpC%zsz(3)))
    allocate(vavg_np2(SM%gpC%zsz(3)))
    allocate(wavg_np2(SM%gpC%zsz(3)))
    allocate(Tavg_np2(SM%gpC%zsz(3)))

    do tid = tidst, tiden, tstep
        call tic()
        ! Read data
        call SM%reinit(tid,restartFromViz=.true.)

        ! Store u^n
        u = SM%u
        v = SM%v
        w = SM%wC
        if (SM%isStratified) T = SM%T

        call message(1,"TID",tid)

        ! Time advance to u^{n+1}
        call SM%timeAdvance(dtForce)

        ! Compute bugdet RHS

        ! Compute ddt(Total KE) with central difference
        SM%rbuffxC(:,:,:,1) = 0.5*(u*u + v*v + w*w)
        SM%rbuffxC(:,:,:,2) = 0.5*(SM%u*SM%u + SM%v*SM%v + SM%wC*SM%wC)
        
        ! Take x-y average
        call get_xy_average(SM%rbuffxC(:,:,:,1), SM, kn)
        call get_xy_average(SM%rbuffxC(:,:,:,2), SM, knp2)
        dkdt = (knp2 - kn)/(dtForce)
        call message_min_max(2,"dkdt bounds",p_maxval(maxval(abs(dkdt))),p_minval(minval(abs(dkdt))))

        ! Compute ddt(MKE)
        call get_xy_average(u,SM,uavg_n)
        call get_xy_average(v,SM,vavg_n)
        call get_xy_average(w,SM,wavg_n)
        call get_xy_average(SM%u,SM,uavg_np2)
        call get_xy_average(SM%v,SM,vavg_np2)
        call get_xy_average(SM%wC,SM,wavg_np2)

        if (SM%isStratified) then
          SM%rbuffxC(:,:,:,1) = (0.5*(SM%T*SM%T) - 0.5*(T*T))/(dtForce)
          call get_xy_average(SM%rbuffxC(:,:,:,1), SM, dTTdt)
          call message_min_max(2,"dTTdt bounds",p_maxval(maxval(abs(dTTdt))),p_minval(minval(abs(dTTdt))))
          call get_xy_average(T,SM,Tavg_n)
          call get_xy_average(SM%T,SM,Tavg_np2)
        end if

        if (nrank == 0) then
            ddt_MKE = (0.5*(uavg_np2*uavg_np2 + vavg_np2*vavg_np2 + wavg_np2*wavg_np2) - &
              0.5*(uavg_n*uavg_n + vavg_n*vavg_n + wavg_n*wavg_n))/(dtForce)

            ! Subtract ddt(MKE) to get ddt(TKE)
            dkdt = dkdt - ddt_MKE

            ! Dump dkdt to disk
            write(tempname,"(A3,I2.2,A7,I6.6,A4)")"Run",SM%runID,"_dkdt_t",tid,".stt"
            fname = trim(outputdir)//'/'//trim(tempname)
            call write_1d_ascii(dkdt,trim(fname))

            if (SM%isStratified) then
              ddt_MPE = (0.5*Tavg_np2*Tavg_np2 - 0.5*Tavg_n*Tavg_n)/(dtForce)
              dTTdt = dTTdt - ddt_MPE
              write(tempname,"(A3,I2.2,A8,I6.6,A4)")"Run",SM%runID,"_dTTdt_t",tid,".stt"
              fname = trim(outputdir)//'/'//trim(tempname)
              call write_1d_ascii(dTTdt,trim(fname))

            end if
        end if

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if (tiden - tidst > 0) then
          call toc(tid,tiden-tidst,whichMssg=2)
        else
          call toc()
        end if
    end do
    contains

      subroutine get_xy_average(f, sim, fmean)
        real(rkind), dimension(:,:,:), intent(in) :: f
        type(igrid), intent(inout) :: sim
        real(rkind), dimension(:), intent(out) :: fmean
        integer :: nx, ny, nz, ierr, k
        real(rkind) :: avgFact

        nx = sim%gpC%xsz(1)
        ny = sim%gpC%ysz(2)
        nz = sim%gpC%zsz(3)
        avgFact = 1.d0/(real(nx,rkind)*real(ny,rkind))

        !call transpose_x_to_y(f,sim%rbuffyC(:,:,:,1),sim%gpC)
        !call transpose_y_to_z(sim%rbuffyC(:,:,:,1),sim%rbuffzC(:,:,:,1),sim%gpC)
        !do k = 1,nz
        !  fmean(k) = p_sum(sum(sim%rbuffzC(:,:,k,1)))*avgFact
        !end do
        
        call sim%spectC%fft(f, sim%cbuffyC(:,:,:,1))
        call transpose_y_to_z(sim%cbuffyC(:,:,:,1), sim%cbuffzC(:,:,:,1), sim%sp_gpC)
        if (nrank == 0) then
          fmean = real(sim%cbuffzC(1,1,:,1),rkind)*avgFact
        else
          fmean = 0.d0
        end if
        call MPI_Bcast(fmean,nz,mpirkind,0,MPI_COMM_WORLD, ierr)

      end subroutine
end program shearlessMixing_Bodart_ddt_tke
