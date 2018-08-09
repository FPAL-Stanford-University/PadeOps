program EkmanLayerDNSTKEBudget
    use kind_parameters, only: rkind, clen
    use igrid_Operators, only: igrid_ops
    use constants, only: pi, two
    use mpi
    use decomp_2d
    use timer, only: tic, toc
    use exits, only: message
    implicit none

    real(rkind), dimension(:,:,:), allocatable :: &
                &buff1, buff2, buff3, buff4, buff5, &
                &u,v,w,p, ufluct, vfluct
    real(rkind), dimension(:,:), allocatable :: &
                &umean_t, vmean_t, &
                &Tran_t, Prod_t, Diff_t, Diss_t, PDif_t, T_da_t, uiui_t, times
    real(rkind) :: &
                &time, dx, dy, dz, &
                &Re=400.d0, Lx=26.d0, Ly=26.d0, Lz=24.d0
    integer ::  nx, ny, nz, nt, RunID, TIDX, botBC=0, topBC=1, &
                &idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
    type(igrid_ops) :: ops
    character(len=clen) ::  inputdir, outputdir, inputfile
    logical :: isZPeriodic = .false. 

    namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, Re, NumericalSchemeVert 
    
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
    allocate(PDif_t(nz,nt))
    allocate(uiui_t(nz,nt))
    allocate(T_da_t(nz,nt))
    allocate(times(1,nt))
 
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
        times(1,idx) = TIDX

        ! Fluctuations and means
        call ops%getFluct_from_MeanZ(u,ufluct)
        call ops%getFluct_from_MeanZ(v,vfluct)
        call ops%TakeMean_xy(u-ufluct,umean_t(:,idx))
        call ops%TakeMean_xy(v-vfluct,vmean_t(:,idx))
       
        ! TKE: uiui = uu+vv+ww
        buff1 = ufluct**2 + vfluct**2 + w**2        
        call ops%TakeMean_xy(buff1, uiui_t(:,idx))
 
        ! Transport terms: d/dz <uuw + vvw + www>
        call ops%ddz(ufluct*ufluct*w + vfluct*vfluct*w + w**3, buff1, botBC, topBC)
        call ops%TakeMean_xy(-buff1, Tran_t(:,idx))            
       
        ! Corrections terms for transport term dealiasing
        call ops%GetGradient(ufluct,buff2,buff3,buff4,botBC,topBC)
        buff5 = ufluct*ufluct
        call ops%dealias(buff5)
        call ops%dealias(buff2)
        buff1 = buff5*buff2         !uu*dudx
        buff5 = ufluct*vfluct
        call ops%dealias(buff5)
        call ops%dealias(buff3)
        buff1 = buff5*buff3 + buff1 !uv*dudy
        buff5 = ufluct*w
        call ops%dealias(buff5)
        call ops%dealias(buff4)
        buff1 = buff5*buff4 + buff1 !uw*dudz
        call ops%GetGradient(vfluct,buff2,buff3,buff4,botBC,topBC)
        buff5 = vfluct*ufluct
        call ops%dealias(buff5)
        call ops%dealias(buff2)
        buff1 = buff5*buff2 + buff1 !vu*dvdz
        buff5 = vfluct*vfluct
        call ops%dealias(buff5)
        call ops%dealias(buff3)
        buff1 = buff5*buff3 + buff1 !vv*dvdz
        buff5 = vfluct*w
        call ops%dealias(buff5)
        call ops%dealias(buff4)
        buff1 = buff5*buff4 + buff1 !vw*dvdz
        call ops%GetGradient(w,buff2,buff3,buff4,botBC,topBC)
        buff5 = w*ufluct
        call ops%dealias(buff5)
        call ops%dealias(buff2)
        buff1 = buff5*buff2 + buff1 !wu*dwdz
        buff5 = w*vfluct
        call ops%dealias(buff5)
        call ops%dealias(buff3)
        buff1 = buff5*buff3 + buff1 !ww*dwdz
        buff5 = w*w
        call ops%dealias(buff5)
        call ops%dealias(buff4)
        buff1 = buff5*buff4 + buff1 !ww*dwdz
        call ops%TakeMean_xy(buff1,T_da_t(:,idx))
 
        
        ! Production: -(<uw>*dUdz + <vw>*dVdz)
        call ops%ddz(u-ufluct, buff1, botBC, topBC)
        buff2 = ufluct*w*buff1
        call ops%ddz(v-vfluct, buff1, botBC, topBC)
        buff2 = buff2 + vfluct*w*buff1
        call ops%TakeMean_xy(-buff2,Prod_t(:,idx))

        ! Diffusion: d2dz2 <uiui>
        call ops%d2dz2(ufluct**2+vfluct**2+w**2,buff1,-1,-1)
        call ops%TakeMean_xy(buff1,Diff_t(:,idx))

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
        call ops%ddz( w*buff1, buff2, 0, 0)
        call ops%TakeMean_xy(-buff2, PDif_t(:,idx))

        call toc()
        idx = idx + 1
    end do 

    ! Write out time averages
    if (nrank == 0) then
        call message(0,"Writing files...")
        call ops%WriteASCII_2D(Tran_t, "Tran")
        call ops%WriteASCII_2D(Prod_t, "Prod")
        call ops%WriteASCII_2D(Diff_t, "Diff")
        call ops%WriteASCII_2D(Diss_t, "Diss")
        call ops%WriteASCII_2D(PDif_t, "PDif")
        call ops%WriteASCII_2D(T_da_t, "T_da")
        call ops%WriteASCII_2D(uiui_t, "uiui")
        call ops%WriteASCII_2D(times, "time")
    end if 
        
    call ops%destroy()
    call MPI_Finalize(ierr)           
    call message(0,"Done")

end program 

