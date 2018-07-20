program EkmanLayerDNSTurbStats
    use kind_parameters, only: rkind, clen
    use igrid_Operators, only: igrid_ops
    use DerivativesMod, only: derivatives
    use constants, only: pi, two
    use mpi
    use decomp_2d
    use timer, only: tic, toc
    use exits, only: message
    implicit none

    real(rkind), dimension(:,:,:), allocatable :: &
                &buff1, buff2, buff3, buff4, buffz1, buffz2, &
                &u,v,w,p, ufluct, vfluct
    real(rkind), dimension(:,:), allocatable :: &
                &umean_t, vmean_t, data2write, &
                &Tran_t, Prod_t, Diff_t, Diss_t, PDif_t
    real(rkind) :: &
                &time, dx, dy, dz, &
                &Re=400.d0, Lx=26.d0, Ly=26.d0, Lz=24.d0
    integer ::  nx, ny, nz, nt, RunID, TIDX, botBC=0, topBC=1, &
                &idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
    type(igrid_ops) :: ops
    type(derivatives) :: der  
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
    call der%init( ops%gp,   dx,     dy,    dz, &
                        .TRUE., .TRUE., .FALSE., &
                        "four", "four", "cd10" )
    call mpi_barrier(mpi_comm_world, ierr)

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
    allocate(buffz1(ops%gp%zsz(1),ops%gp%zsz(2),ops%gp%zsz(3)))
    allocate(buffz2(ops%gp%zsz(1),ops%gp%zsz(2),ops%gp%zsz(3)))

    ! Get number of timesteps and allocate arrays
    nt = (tstop - tstart)/tstep
    allocate(umean_t(nz,nt))
    allocate(vmean_t(nz,nt))
    allocate(Tran_t(nz,nt))
    allocate(Prod_t(nz,nt))
    allocate(Diff_t(nz,nt))
    allocate(Diss_t(nz,nt))
    allocate(PDif_t(nz,nt))
 
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
        time = ops%getSimTime(tidx)
        call mes:tabsage(0, "Read simulation data at time:", time)

        ! Fluctuations and means
        call ops%getFluct_from_MeanZ(u,ufluct)
        call ops%getFluct_from_MeanZ(v,vfluct)
        call ops%TakeMean_xy(u-ufluct,umean_t(:,idx))
        call ops%TakeMean_xy(v-vfluct,vmean_t(:,idx))
        
        ! Transport terms: T13 = d/dz <u'ww> 
        call ops%ddz(ufluct*w*w, buff1, botBC, topBC)
        call ops%TakeMean_xy(-buff1, Tran_13(:,idx))            
        
        ! Production: -(<uw>*dUdz + <vw>*dVdz)
        call ops%ddz(u-ufluct, buff1, botBC, topBC)
        buff2 = ufluct*w*buff1
        call ops%ddz(v-vfluct, buff1, botBC, topBC)
        buff2 = buff2 + vfluct*w*buff1
        call ops%TakeMean_xy(-buff2,Prod_t(:,idx))

        ! Diffusion: d2dz2 <uiui>
        buff1 = ufluct**2+vfluct**2+w**2
        call transpose_y_to_z(buff1,buffz1,ops%gp)
        call der % d2dz2(buffz1,buffz2)
        call transpose_z_to_y(buffz2,buff2,ops%gp)
        call ops%TakeMean_xy(buff2,Diff_t(:,idx))

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
        call ops%TakeMean_xy(-buff2, PDif_t(:,idx))

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
        data2write(:,5) = sum(PDif_t,2)/nt
        call ops%WriteASCII_2D(data2write, "test")
    end if 
        
    call ops%destroy()
    call MPI_Finalize(ierr)           

end program 

