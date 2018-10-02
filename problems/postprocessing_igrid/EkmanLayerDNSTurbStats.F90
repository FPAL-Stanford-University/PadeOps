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
                &ufluct, vfluct, u, v, w, Rij_t
    real(rkind), dimension(:,:), allocatable :: &
                &data2write, Rij, &
                &umean_t, vmean_t
    real(rkind), dimension(:), allocatable :: &
                &umean, vmean,&
                &d2Udz2, d2Vdz2,&
                &dUdz, dVdz, duwdz, dvwdz
    real(rkind) :: &
                &time, dx, dy, dz, &
                &Re=400.d0, Lx=26.d0, Ly=26.d0, Lz=24.d0
    integer ::  nx, ny, nz, nt, RunID, TIDX,&
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
    allocate(Rij_t(nz,6,nt))
    allocate(umean(nz))
    allocate(vmean(nz))
    allocate(Rij(nz,6))
    allocate(dUdz(nz))
    allocate(dVdz(nz))
    allocate(duwdz(nz))
    allocate(dvwdz(nz))
    allocate(d2Udz2(nz))
    allocate(d2Vdz2(nz))
 
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
        call ops%TakeMean_xy(ufluct*ufluct, Rij_t(:,1,idx))        
        call ops%TakeMean_xy(ufluct*vfluct, Rij_t(:,2,idx))        
        call ops%TakeMean_xy(ufluct*w,      Rij_t(:,3,idx))        
        call ops%TakeMean_xy(vfluct*vfluct, Rij_t(:,4,idx))        
        call ops%TakeMean_xy(vfluct*w,      Rij_t(:,5,idx))       
        call ops%TakeMean_xy(w*w,           Rij_t(:,6,idx))       
       
        idx = idx + 1
    end do 

    ! Averages in time:
    umean = sum(umean_t,2)/nt
    vmean = sum(vmean_t,2)/nt
    Rij = sum(Rij_t,3)/nt
        
    ! d/dz terms:
    call ops%ddz_1d(umean, dUdz,0,1)
    call ops%ddz_1d(vmean, dVdz,0,1)
    call ops%ddz_1d(Rij(:,3), duwdz,-1,-1)
    call ops%ddz_1d(Rij(:,5), dvwdz,-1,-1)
    
    ! Viscous stresses
    call ops%d2dz2_1d(umean, d2Udz2,0,-1)
    call ops%d2dz2_1d(vmean, d2Vdz2,0,-1)

    ! Write out time averages
    if (nrank == 0) then
        allocate(data2write(nz,8))
        data2write(:,1) = umean
        data2write(:,2) = vmean
        data2write(:,3) = dUdz 
        data2write(:,4) = dVdz
        data2write(:,5) = d2Udz2
        data2write(:,6) = d2Vdz2
        data2write(:,7) = duwdz
        data2write(:,8) = dvwdz
        call ops%WriteASCII_2D(data2write, "tavg")
        call ops%WriteASCII_2D(Rij, "uiuj")
    end if 
        
    call ops%destroy()
    call MPI_Finalize(ierr)           

end program 

