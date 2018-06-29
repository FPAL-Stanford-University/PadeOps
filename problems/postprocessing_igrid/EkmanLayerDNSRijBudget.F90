program EkmanLayerDNSTurbStats
    use kind_parameters, only: rkind, clen
    use igrid_Operators, only: igrid_ops
    use PadeDerOps, only: Pade6stagg
    use constants, only: pi, two
    use mpi
    use decomp_2d, only: nrank 
    use timer, only: tic, toc
    use exits, only: message
    implicit none

    real(rkind), dimension(:,:,:), allocatable :: buff1,buff2,buff3,buff4, u,v,w,p, ufluct,vfluct
    real(rkind), dimension(:,:), allocatable :: umean_t, vmean_t, data2write, &
                                                &Tran_t, Prod_t, Diff_t, Diss_t, PD_t
    real(rkind) :: time, dx, dy, dz, Re=400.d0, Lx=26.d0, Ly=26.d0, Lz=24.d0
    integer :: nx, ny, nz, nt, RunID, TIDX, botBC=0, topBC=1, &
                &idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
    type(igrid_ops) :: ops
    type(Pade6stagg) :: Pade6opZ    
    character(len=clen) ::  inputdir, outputdir, inputfile
    logical :: isZPeriodic = .false. 

    namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep,&
             nx, ny, nz, Re, NumericalSchemeVert 
    
    call MPI_Init(ierr)               
    call GETARG(1,inputfile)          
    open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=99, NML=INPUT)
    close(unit=99)

    dx = Lx/real(nx,rkind) 
    dy = Ly/real(ny,rkind) 
    dz = Lz/real(nz,rkind)
    
    ! Initialize the operator class and derivative class
    call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NumericalSchemeVert)
    !call decomp_2d_init(nx, ny, nz, 0, 0)
    !call get_decomp_info(gpC)
    !call decomp_info_init(nx, ny, nz+1, gpE)
    !call spectC%init("x", nx, ny, nz, dx, dy, dz, "four", '2/3rd', 2 , .false.)
    !call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", '2/3rd', 2 ,.false.)
    !sp_gpC => spectC%spectdecomp
    !sp_gpE => spectE%spectdecomp
    !call Pade6opz%init(gpC, sp_gpC, gpE, sp_gpE, dz, NumericalSchemeVert, isZperiodic)

    ! Allocate all the needed memory 
    call ops%allocate3DField(buff1)
    call ops%allocate3DField(buff2)
    call ops%allocate3DField(buff3)
    call ops%allocate3DField(buff4)
    call ops%allocate3DField(u)
    call ops%allocate3DField(v)
    call ops%allocate3DField(w)
    call ops%allocate3DField(p)
    call ops%allocate3DField(ufluct)
    call ops%allocate3DField(vfluct)

    ! Get number of timesteps and allocate arrays
    nt = (tstop - tstart)/tstep
    allocate(umean_t(nz,nt))
    allocate(vmean_t(nz,nt))
    allocate(Tran_t(nz,nt))
    allocate(Prod_t(nz,nt))
    allocate(Diff_t(nz,nt))
    allocate(Diss_t(nz,nt))
    allocate(PD_t(nz,nt))
 
    ! Compute for each timestep the z profiles 
    idx = 1
    do while(idx <= nt)
        TIDX = tstart + tstep * (idx - 1)
        call message(0, "Reading fields for tid:", TIDX)
        call tic()
        call ops%ReadField3D(u,"uVel",TIDX)
        call ops%ReadField3D(v,"vVel",TIDX)
        call ops%ReadField3D(w,"wVel",TIDX)
        call ops%ReadField3D(p,"prss",TIDX)
        time = ops%getSimTime(tidx)
        call message(0, "Read simulation data at time:", time)

        ! Fluctuations and means
        call ops%getFluct_from_MeanZ(u,ufluct)
        call ops%getFluct_from_MeanZ(v,vfluct)
        call ops%TakeMean_xy(u-ufluct,umean_t(:,idx))
        call ops%TakeMean_xy(v-vfluct,vmean_t(:,idx))
        
        ! Transport terms: d/dz <uuw + vvw + www>
        call ops%ddz(ufluct*ufluct*w + vfluct*vfluct*w + w**3, buff1, botBC, topBC)
        call ops%TakeMean_xy(-buff1, Tran_t(:,idx))            
        
        ! Production: -(<uw>*dUdz + <vw>*dVdz)
        call ops%ddz(u-ufluct, buff1, botBC, topBC)
        buff2 = ufluct*w*buff1
        call ops%ddz(v-vfluct, buff1, botBC, topBC)
        buff2 = buff2 + ufluct*w*buff1
        call ops%TakeMean_xy(-buff2,Prod_t(:,idx))

        ! Diffusion: d2dz2 <uiui>
        ! call Pade6opZ%d2dz2_C2C(ufluct**2+vfluct**2+w**2,buff1,botBC,topBC)
        Diff_t(:,idx) = 0

        ! Viscous dissip: -2<dui/dxk*dui/dxk>
        call ops%GetGradient(ufluct,buff1,buff2,buff3,botBC,topBC)
        buff4 = buff1**2 + buff2**2 + buff3**2
        call ops%GetGradient(vfluct,buff1,buff2,buff3,botBC,topBC)
        buff4 = buff4 + buff1**2 + buff2**2 + buff3**2
        call ops%GetGradient(w, buff1,buff2,buff3,botBC,topBC)
        buff4 = buff4 + buff1**2 + buff2**2 + buff3**2
        call ops%TakeMean_xy(-2*buff4, Diss_t(:,idx))

        ! Pressure diffusion: -d/dz<wp'>
        call ops%getFluct_from_MeanZ(p,buff1)
        call ops%ddz( w*buff1, buff2, botBC, topBC )
        call ops%TakeMean_xy(-buff2, PD_t(:,idx))

        call toc()
        idx = idx + 1
    end do 

    ! Write out time averages
    if (nrank == 0) then
        allocate(data2write(nz,5))
        data2write(:,1) = sum(Tran_t,2)/nt
        data2write(:,2) = sum(Prod_t,2)/nt
        data2write(:,3) = sum(Diff_t,2)/nt
        data2write(:,4) = sum(Diss_t,2)/nt
        data2write(:,5) = sum(PD_t,2)/nt
        call ops%WriteASCII_2D(data2write, "tkeb")
    end if 
        
    call ops%destroy()
    call MPI_Finalize(ierr)           

end program 

