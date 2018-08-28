program EkmanLayerDNSMomentumBalance
    use kind_parameters, only: rkind, clen
    use igrid_Operators, only: igrid_ops
    use constants, only: pi, two
    use mpi
    use decomp_2d, only: nrank 
    use timer, only: tic, toc
    use exits, only: message
    implicit none

    real(rkind), dimension(:,:,:), allocatable :: &
                &buff1, buff2, buff3, buff4, &
                &ufluct, vfluct, u, v, w
    real(rkind), dimension(:,:), allocatable :: &
                &data2write,&
                &umean_t, vmean_t, &
                &uu_t, uv_t, uw_t, vv_t, vw_t, ww_t 
    real(rkind), dimension(:), allocatable :: &
                &umean, vmean, &
                &uu, uv, uw, vv, vw, ww, &
                &d2Udz2, d2Vdz2,&
                &dUdz, dVdz, duwdz, dvwdz
    real(rkind) :: &
                &time, dx, dy, dz, &
                &Re=400.d0, Lx=26.d0, Ly=26.d0, Lz=24.d0
    integer ::  nx, ny, nz, nt, RunID, TIDX, botBC=0, topBC=1, &
                &idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
    type(igrid_ops) :: ops
    character(len=clen) ::  inputdir, outputdir, inputfile
    logical :: isZPeriodic = .false. 

    namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, Re, NumericalSchemeVert 
    
    call MPI_Init(ierr)               
    call GETARG(1,inputfile)          
    open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=99, NML=INPUT)
    close(unit=99)

    dx = Lx/real(nx,rkind) 
    dy = Ly/real(ny,rkind) 
    dz = Lz/real(nz,rkind)
    
    ! Initialize the operator class
    call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NumericalSchemeVert)
     
    ! Allocate all the needed memory 
    call ops%allocate3DField(buff1)
    call ops%allocate3DField(buff2)
    call ops%allocate3DField(buff3)
    call ops%allocate3DField(buff4)
    call ops%allocate3DField(u)
    call ops%allocate3DField(v)
    call ops%allocate3DField(w)
    call ops%allocate3DField(ufluct)
    call ops%allocate3DField(vfluct)
    
    ! Get number of timesteps and allocate arrays
    nt = (tstop - tstart)/tstep
    allocate(umean_t(nz,nt))
    allocate(vmean_t(nz,nt))
    allocate(umean(nz))
    allocate(vmean(nz))
    allocate(uu_t(nz,nt))
    allocate(uv_t(nz,nt))
    allocate(uw_t(nz,nt))
    allocate(vv_t(nz,nt))
    allocate(vw_t(nz,nt))
    allocate(ww_t(nz,nt))
    allocate(uu(nz))
    allocate(uv(nz))
    allocate(uw(nz))
    allocate(vv(nz))
    allocate(vw(nz))
    allocate(ww(nz))
    allocate(dUdz(nz))
    allocate(dVdz(nz))
    allocate(duwdz(nz))
    allocate(dvwdz(nz))
    allocate(d2Udz2(nz))
    allocate(d2Vdz2(nz))
 
    idx = 1
    do while(idx <= nt)
        TIDX = tstart + tstep * (idx - 1)
        call message(0, "Reading fields for tid:", TIDX)
        call tic()
        call ops%ReadField3D(u,"uVel",TIDX)
        call ops%ReadField3D(v,"vVel",TIDX)
        call ops%ReadField3D(w,"wVel",TIDX)
        time = ops%getSimTime(tidx)
        call message(0, "Read simulation data at time:", time)

        ! Fluctuations and means
        call ops%getFluct_from_MeanZ(u,ufluct)
        call ops%getFluct_from_MeanZ(v,vfluct)
        call ops%TakeMean_xy(u-ufluct,umean_t(:,idx))
        call ops%TakeMean_xy(v-vfluct,vmean_t(:,idx))
 
        ! Reynolds stresses
        call ops%TakeMean_xy(ufluct*ufluct, uu_t(:,idx))        
        call ops%TakeMean_xy(ufluct*vfluct, uv_t(:,idx))        
        call ops%TakeMean_xy(ufluct*w,      uw_t(:,idx))        
        call ops%TakeMean_xy(vfluct*vfluct, vv_t(:,idx))        
        call ops%TakeMean_xy(vfluct*w,      vw_t(:,idx))       
        call ops%TakeMean_xy(w*w,           ww_t(:,idx))       
        
        call toc()
        idx = idx + 1
    end do 
    
    ! Average in time
    umean = sum(umean_t,2)/nt 
    vmean = sum(vmean_t,2)/nt 
    uu = sum(uu_t,2)/nt 
    uv = sum(uv_t,2)/nt 
    uw = sum(uw_t,2)/nt 
    vv = sum(vv_t,2)/nt 
    vw = sum(vw_t,2)/nt 
    ww = sum(ww_t,2)/nt 
       
    ! Derivatives on the mean: dU/dz, d2U/dz2
    call ops%ddz_1d(umean, dUdz, -1, 1 )
    call ops%ddz_1d(vmean, dVdz, -1, 1 )
    call ops%d2dz2_1d(umean, d2Udz2, 0, 0 )
    call ops%d2dz2_1d(vmean, d2Vdz2, 0, 0 )
 
    ! Derivatives on the stresses: 
    call ops%ddz_1d(uw, duwdz, 0, 0)
    call ops%ddz_1d(vw, dvwdz, 0, 0)
    
    ! Viscous stresses
    call ops%d2dz2_1d(umean, d2Udz2,0,0)
    call ops%d2dz2_1d(vmean, d2Vdz2,0,0)
    
    ! Write out time averages
    if (nrank == 0) then
        call message(0,"Writing files")
        allocate(data2write(Nz,14))
        data2write(:,1) = umean
        data2write(:,2) = vmean 
        data2write(:,3) = uu
        data2write(:,4) = uv
        data2write(:,5) = uw
        data2write(:,6) = vv
        data2write(:,7) = vw
        data2write(:,8) = ww
        data2write(:,9) = dUdz
        data2write(:,10) = dVdz
        data2write(:,11) = d2Udz2
        data2write(:,12) = d2Vdz2
        data2write(:,13) = duwdz
        data2write(:,14) = dvwdz
        call ops%WriteASCII_2D(data2write, "stat")
    end if 
    
    call ops%destroy()
    call MPI_Finalize(ierr)           
    call message(0,"done") 
end program 

