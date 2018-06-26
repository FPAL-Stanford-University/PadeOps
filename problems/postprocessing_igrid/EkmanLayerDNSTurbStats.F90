program EkmanLayerDNSTurbStats
    use kind_parameters, only: rkind, clen
    use igrid_Operators, only: igrid_ops
    use constants, only: pi, two
    use mpi
    use decomp_2d, only: nrank 
    use timer, only: tic, toc
    use exits, only: message
    implicit none

    real(rkind), dimension(:,:,:), allocatable :: buff1,buff2,buff3,buff4, u,v,w, ufluct,vfluct
    real(rkind), dimension(:,:), allocatable :: umean_t,vmean_t,data2write, data2write2, uu_t,uv_t,uw_t, vv_t,vw_t,ww_t, uuw_t,vvw_t,www_t, Diss_t
    real(rkind), dimension(:), allocatable :: buff1d, umean, vmean, uu,uv,uw,vv,vw,ww, Tran,Prod,Diff,Diss,Cori
    real(rkind) :: time, dx, dy, dz, Re=400.d0, Lx=26.d0, Ly=26.d0, Lz=24.d0!, nu=1.145d-5, G=0.0115, D=3.9626d-1
    integer :: nx, ny, nz, nt, RunID, TIDX, botBC=0, topBC=1, idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
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
    allocate(uuw_t(nz,nt))
    allocate(vvw_t(nz,nt))
    allocate(www_t(nz,nt))
    allocate(Diss_t(nz,nt))
    allocate(uu(nz))
    allocate(uv(nz))
    allocate(uw(nz))
    allocate(vv(nz))
    allocate(vw(nz))
    allocate(ww(nz))
    allocate(Tran(nz))
    allocate(Prod(nz))
    allocate(Diff(nz))
    allocate(Diss(nz))
    allocate(Cori(nz))
    allocate(buff1D(nz))
 
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
       
        ! TKE Budget terms:
        ! Transport terms:      Tran = d/dz <uuw + vvw + www>
        call ops%TakeMean_xy(ufluct*ufluct*w, uuw_t(:,idx))       
        call ops%TakeMean_xy(vfluct*vfluct*w, vvw_t(:,idx))       
        call ops%TakeMean_xy(w**3,            www_t(:,idx))       
        ! Viscous dissip:       Diss = - <dui/dxk*dui/dxk>
        call ops%GetGradient(ufluct,buff1,buff2,buff3,botBC,topBC)
        buff4 = buff1**2 + buff2**2 + buff3**2
        call ops%GetGradient(vfluct,buff1,buff2,buff3,botBC,topBC)
        buff4 = buff4 + buff1**2 + buff2**2 + buff3**2
        call ops%GetGradient(w, buff1,buff2,buff3,botBC,topBC)
        buff4 = buff4 + buff1**2 + buff2**2 + buff3**2
        call ops%TakeMean_xy(buff4, Diss_t(:,idx))

        call toc()
        idx = idx + 1
    end do 

    ! Mean velocities, time avgd
    umean = sum(umean_t,2)/nt
    vmean = sum(vmean_t,2)/nt 

    ! Reynold stresses, time avgd
    uu = sum(uu_t,2)/nt
    uv = sum(uv_t,2)/nt
    uw = sum(uw_t,2)/nt
    vv = sum(vv_t,2)/nt
    vw = sum(vw_t,2)/nt
    ww = sum(ww_t,2)/nt

        ! Production
        buff1(1,1,:) = umean 
        call ops%ddz(buff1, buff2, botBC, topBC)
        Prod = uw*buff2(1,1,:)
        buff1(1,1,:) = vmean
        call ops%ddz(buff1, buff2, botBC, topBC)
        Prod = Prod + vw*buff2(1,1,:)
       
         ! Transport terms
        buff1D = sum(uuw_t,2)/nt
        buff1D = sum(vvw_t,2)/nt + buff1D
        buff1D = sum(www_t,2)/nt + buff1D
        call ops%ddz_1d( buff1D, Tran)
        
        ! Viscous diffusion 
        !call ops%d2dz2( uu+vv+ww, Diff )
        Diff = 0

        ! Viscous dissipation
        Diss = sum(Diss_t,2)/nt

    ! Write out time averages
    if (nrank == 0) then
        allocate(data2write(nz,8))
        data2write(:,1) = sum(umean_t,2)/nt
        data2write(:,2) = sum(vmean_t,2)/nt 
        data2write(:,3) = sum(uu_t,2)/nt 
        data2write(:,4) = sum(uv_t,2)/nt 
        data2write(:,5) = sum(uw_t,2)/nt 
        data2write(:,6) = sum(vv_t,2)/nt 
        data2write(:,7) = sum(vw_t,2)/nt 
        data2write(:,8) = sum(ww_t,2)/nt
        call ops%WriteASCII_2D(data2write, "tavg")
        
        allocate(data2write2(nz,4))
        data2write(:,1) = Tran
        data2write(:,2) = Prod
        data2write(:,3) = Diff
        data2write(:,4) = Diss
        call ops%WriteASCII_2D(data2write2, "tkeb")
    end if 
        
    call ops%destroy()
    call MPI_Finalize(ierr)           

end program 

