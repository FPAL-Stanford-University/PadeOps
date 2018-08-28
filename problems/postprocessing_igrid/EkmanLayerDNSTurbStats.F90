program EkmanLayerDNSTurbStats
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
                &data2write, &
                &umean_t, vmean_t, &
                &uu_t, uv_t, uw_t, vv_t, vw_t, ww_t, &
                &d2Udz2_t, d2Vdz2_t,&
                &dUdz_t, dVdz_t, duwdz_t, dvwdz_t
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
    allocate(uu_t(nz,nt))
    allocate(uv_t(nz,nt))
    allocate(uw_t(nz,nt))
    allocate(vv_t(nz,nt))
    allocate(vw_t(nz,nt))
    allocate(ww_t(nz,nt))
    allocate(dUdz_t(nz,nt))
    allocate(dVdz_t(nz,nt))
    allocate(duwdz_t(nz,nt))
    allocate(dvwdz_t(nz,nt))
    allocate(d2Udz2_t(nz,nt))
    allocate(d2Vdz2_t(nz,nt))
 
    ! Compute for each timestep the z profiles 
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
       
        ! Viscous stresses
        call ops%d2dz2(u-ufluct, buff1,-1,-1)
        call ops%TakeMean_xy(buff1, d2Udz2_t(:,idx))
        call ops%d2dz2(v-vfluct, buff2,-1,-1)
        call ops%TakeMean_xy(buff2, d2Vdz2_t(:,idx))
        call toc()
        idx = idx + 1
    end do 

    ! d/dz terms: time avgd, fn of height
    call ops%ddz_1d(sum(umean_t,2)/nt, dUdz,0,0)
    call ops%ddz_1d(sum(vmean_t,2)/nt, dVdz,0,0)
    call ops%ddz_1d(sum(uw_t,2)/nt, duwdz,0,0)
    call ops%ddz_1d(sum(vw_t,2)/nt, dvwdz,0,0)
    
    
    ! Write out time averages
    if (nrank == 0) then
       allocate(data2write(nz,10))
       !data2write(:,1) = sum(umean_t,2)/nt
       !data2write(:,2) = sum(vmean_t,2)/nt
       !data2write(:,3) = sum(uu_t,2)/nt 
       !data2write(:,4) = sum(uv_t,2)/nt
       !data2write(:,5) = sum(uw_t,2)/nt
       !data2write(:,6) = sum(vv_t,2)/nt
       !data2write(:,7) = sum(vw_t,2)/nt 
       !data2write(:,8) = sum(ww_t,2)/nt 
       !data2write(:,9) = dUdz
       !data2write(:,10) = dVdz
       !call ops%WriteASCII_2D(data2write, "tavg")
       !data2write(:,1) = dUdz 
       !data2write(:,2) = dVdz
       !data2write(:,3) = sum(d2Udz2_t,2)/nt
       !data2write(:,4) = sum(d2Vdz2_t,2)/nt
       !data2write(:,5) = duwdz
       !data2write(:,6) = dvwdz
       !data2write(:,7) = 0
       !data2write(:,8) = 0
       !data2write(:,9) = 0
       !data2write(:,10) = 0 
       !call ops%WriteASCII_2D(data2write, "d_dz")
    end if 
        
    call ops%destroy()
    call MPI_Finalize(ierr)           

end program 

