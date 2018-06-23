program EkmanLayerDNSTurbStats
    use kind_parameters, only: rkind, clen
    use igrid_Operators, only: igrid_ops
    use constants, only: pi, two
    use mpi
    use decomp_2d, only: nrank 
    use timer, only: tic, toc
    use exits, only: message
    implicit none

    real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3
    real(rkind), dimension(:,:,:), allocatable :: u,v,w, ufluct,vfluct, umean,vmean, dUdx,dUdy,dUdz,dVdx,dVdy,dVdz
    real(rkind), dimension(:,:), allocatable :: ustar, vstar, time, uu,vv,uv,ww, S11,S12,S22
    real(rkind) :: dx, dy, dz, Re = 400.d0, nu = 1.145d-5, Lx = 26.d0, Ly = 26.d0, Lz = 24.d0
    integer :: nx, ny, nz, nt, RunID, TIDX, botBC=0, topBC=1
    type(igrid_ops) :: ops
    character(len=clen) ::  inputdir, outputdir, inputfile
    integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
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
    call ops%allocate3DField(u)
    call ops%allocate3DField(v)
    call ops%allocate3DField(w)
    call ops%allocate3DField(ufluct)
    call ops%allocate3DField(vfluct)
    call ops%allocate3DField(umean)
    call ops%allocate3DField(vmean)
    call ops%allocate3DField(dUdx)
    call ops%allocate3DField(dUdy)
    call ops%allocate3DField(dUdz)
    call ops%allocate3DField(dVdx)
    call ops%allocate3DField(dVdy)
    call ops%allocate3DField(dVdz)
    call ops%allocate3DField(buff1)
    call ops%allocate3DField(buff2)
    call ops%allocate3DField(buff3)

    ! Get number of timesteps and allocate arrays
    nt = (tstop - tstart)/tstep
    allocate(ustar(nz,nt))
    allocate(vstar(nz,nt))
    allocate(time(nt,1))
    allocate(S11(nz,nt))
    allocate(S22(nz,nt))
    allocate(S12(nz,nt))
    !allocate(uu(nz,nt))
    !allocate(vv(nz,nt))
    !allocate(uv(nz,nt))
    !allocate(ww(nz,nt))
    
    ! Compute for each timestep: 
    idx = 1
    do while(idx <= nt)
        TIDX = tstart + tstep * (idx - 1)
        call message(0, "Reading fields for tid:", TIDX)
        call tic()
        call ops%ReadField3D(u,"uVel",TIDX)
        call ops%ReadField3D(v,"vVel",TIDX)
        call ops%ReadField3D(w,"wVel",TIDX)
        time(idx,1) = ops%getSimTime(tidx)
        call message(0, "Read simulation data at time:", time(idx,1))

        ! Fluctuations and means
        call ops%getFluct_from_MeanZ(u,ufluct)
        call ops%getFluct_from_MeanZ(v,vfluct)
        buff1 = u - ufluct
        buff2 = v - vfluct
        
        ! Get gradients of the mean: dUdx, dUdy, dUdz and Sij
        call ops%GetGradient(buff1,dUdx,dUdy,dUdz,botBC,topBC) 
        call ops%GetGradient(buff2,dVdx,dVdy,dVdz,botBC,topBC) 
        call ops%TakeMean_xy(dUdx, S11(:,idx))
        call ops%TakeMean_xy(dVdx, S22(:,idx))
        call ops%TakeMean_xy((dUdy+dVdx)/2, S12(:,idx))
        !
        !! Reynolds stresses ui*uj
        !call ops%TakeMean_xy(ufluct*ufluct, uu(:,idx))
        !call ops%TakeMean_xy(vfluct*vfluct, vv(:,idx))
        !call ops%TakeMean_xy(ufluct*vfluct, uv(:,idx))
        !call ops%TakeMean_xy(w*w, ww(:,idx))

        ! Friction velocity sqrt(nu*tau_wall)
        call ops%TakeMean_xy(dUdz,ustar(:,idx)) 
        call ops%TakeMean_xy(dVdz,vstar(:,idx)) 
        ustar(:,idx) = sqrt(nu*ustar(:,idx))
        vstar(:,idx) = sqrt(nu*ustar(:,idx))

        call toc()
        idx = idx + 1
    end do 
   
    if (nrank == 0) then
       call ops%WriteASCII_2D(ustar,"ustar")
       call ops%WriteASCII_2D(vstar,"vstar")
       call ops%WriteASCII_2D(time, "times")
    end if 

    call ops%destroy()
    call MPI_Finalize(ierr)           

end program 

