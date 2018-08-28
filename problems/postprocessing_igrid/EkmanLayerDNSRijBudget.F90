program EkmanLayerDNSRijBudget
    use kind_parameters, only: rkind, clen
    use igrid_Operators, only: igrid_ops
    use constants, only: pi, two
    use mpi
    use decomp_2d
    use timer, only: tic, toc
    use exits, only: message
    implicit none

    real(rkind), dimension(:,:,:), allocatable :: &
                &buff1, buff2, buff3, buff4, buff5, buff6, &
                &u,v,w,p, ufluct, vfluct, pfluct
    real(rkind), dimension(:,:), allocatable :: &
                &umean_t, vmean_t, uw_t, vw_t, ww_t, &
                &Tran, Prod, Diff, Diss, PDif, PStr
    real(rkind) :: &
                &time, dx, dy, dz, &
                &Re=400.d0, Lx=26.d0, Ly=26.d0, Lz=24.d0
    integer ::  nx, ny, nz, nt, RunID, TIDX, botBC=0, topBC=1, &
                &idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
    type(igrid_ops) :: ops
    character(len=clen) ::  inputdir, outputdir, inputfile
    logical :: isZPeriodic = .false. 

    namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep,&
             nx, ny, nz, Re, NumericalSchemeVert 
    
    call GETARG(1,inputfile)          
    open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=99, NML=INPUT)
    close(unit=99)

    dx = Lx/real(nx,rkind) 
    dy = Ly/real(ny,rkind) 
    dz = Lz/real(nz,rkind)
    
    ! Initialize the operator class and derivative class
    call MPI_Init(ierr)   
    call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, &
                   & RunID, isZPeriodic, NumericalSchemeVert)
    call get_decomp_info(ops%gp)
    call mpi_barrier(mpi_comm_world, ierr)

    ! Allocate all the needed memory 
    call ops%allocate3DField(buff1)
    call ops%allocate3DField(buff2)
    call ops%allocate3DField(buff3)
    call ops%allocate3DField(buff4)
    call ops%allocate3DField(buff5)
    call ops%allocate3DField(buff6)
    call ops%allocate3DField(u)
    call ops%allocate3DField(v)
    call ops%allocate3DField(w)
    call ops%allocate3DField(p)
    call ops%allocate3DField(ufluct)
    call ops%allocate3DField(vfluct)
    call ops%allocate3DField(pfluct)

    ! Get number of timesteps and allocate arrays
    nt = (tstop - tstart)/tstep
    allocate(umean_t(nz,nt))
    allocate(vmean_t(nz,nt))
    allocate(Tran(nz,nt))
    allocate(Prod(nz,nt))
    allocate(Diff(nz,nt))
    allocate(Diss(nz,nt))
    allocate(PDif(nz,nt))
    allocate(PStr(nz,nt))
    allocate(uw_t(nz,nt))
    allocate(vw_t(nz,nt))
    allocate(ww_t(nz,nt))
 
    ! Compute for each timestep: 
    idx = 1
    do while(idx <= nt)
        call tic()
        TIDX = tstart + tstep * (idx - 1)
        call message(0, "Reading fields for tid:", TIDX)
        call ops%ReadField3D(u,"uVel",TIDX)
        call ops%ReadField3D(v,"vVel",TIDX)
        call ops%ReadField3D(w,"wVel",TIDX)
        call ops%ReadField3D(p,"prss",TIDX)
        time = ops%getSimTime(TIDX)
        call message(0, "Read simulation data at time:", time)

        ! Fluctuations and means
        call ops%getFluct_from_MeanZ(u,ufluct)
        call ops%getFluct_from_MeanZ(v,vfluct)
        call ops%getFluct_from_MeanZ(p,pfluct)
        call ops%TakeMean_xy(u-ufluct,umean_t(:,idx))
        call ops%TakeMean_xy(v-vfluct,vmean_t(:,idx))
       
        call ops%TakeMean_xy( ufluct*w, uw_t(:,idx))
        call ops%TakeMean_xy( vfluct*w, vw_t(:,idx))
        call ops%TakeMean_xy( w*w, ww_t(:,idx))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! uw budget
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Transport: -d/dz<uww>
        call ops%ddz( ufluct*w*w, buff1, botBC, topBC)
        call ops%TakeMean_xy( -buff1, Tran(:,idx))            
        
        ! Production: -<ww>dUdz
        call ops%ddz(u-ufluct, buff1, botBC, topBC)
        call ops%TakeMean_xy( -w*w*buff1,Prod(:,idx))

        ! Diffusion: d2dz2<uw>
        call ops%d2dz2( ufluct*w, buff1, -1, 1 )
        call ops%TakeMean_xy( buff1, Diff(:,idx) )

        ! Viscous Dissip: -2<du/dxk*dw/dxk>
        call ops%GetGradient( ufluct, buff1, buff2, buff3, 0,0)
        call ops%GetGradient( w,      buff4, buff5, buff6, 0,0)
        call ops%TakeMean_xy(-2*(buff1*buff4+buff2*buff5+buff3*buff6), Diss(:,idx))

        ! Pressure diffusion: -d/dz<up'>
        call ops%ddz( ufluct*pfluct, buff1, 0, 0 )
        call ops%TakeMean_xy(-buff1, PDif(:,idx))

        ! Pressure strain: p'( du'/dz + dw'/dx )
        call ops%ddz( ufluct, buff1,0 , topBC )
        call ops%ddz( w,      buff2,0 , topBC )
        call ops%TakeMean_xy(pfluct*(buff1+buff2), PStr(:,idx))
       
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ww budget
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Transport: -d/dz<www>
        !call ops%ddz( w**3, buff1, botBC, topBC)
        !call ops%TakeMean_xy( -buff1, Tran(:,idx))            
        !
        !! Production: -<ww>dWdz
        !Prod(:,idx) = 0

        !! Diffusion: d2dz2<ww>
        !call ops%d2dz2( w*w, buff1, -1, 1 )
        !call ops%TakeMean_xy( buff1, Diff(:,idx) )

        !! Viscous Dissip: -2<dw/dxk*dw/dxk>
        !call ops%GetGradient( w, buff1, buff2, buff3, 0,0)
        !call ops%TakeMean_xy( -2*(buff1**2+buff2**2+buff3**2), Diss(:,idx))

        !! Pressure diffusion: -2d/dz<wp'>
        !call ops%ddz( w*pfluct, buff1, 0, 0 )
        !call ops%TakeMean_xy( -2*buff1, PDif(:,idx))

        !! Pressure strain: 2p'dw/dz
        !call ops%ddz( w, buff1, 0 , topBC )
        !call ops%TakeMean_xy( pfluct*buff1, PStr(:,idx))

        !call toc()
        idx = idx + 1
    end do 

    ! Write out time averages
    if (nrank == 0) then
        call message(0,"Writing files...")
        call ops%WriteASCII_2D(Tran, "Tran")
        call ops%WriteASCII_2D(Prod, "Prod")
        call ops%WriteASCII_2D(Diff, "Diff")
        call ops%WriteASCII_2D(Diss, "Diss")
        call ops%WriteASCII_2D(PDif, "PDif")
        call ops%WriteASCII_2D(PStr, "PStr")
        call ops%WriteASCII_2D(uw_t, "uw_t")
        call ops%WriteASCII_2D(vw_t, "vw_t")
        call ops%WriteASCII_2D(vw_t, "ww_t")
    end if 
        
    call ops%destroy()
    call MPI_Finalize(ierr)           
    call message(0,"Done")

end program 

