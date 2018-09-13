program EkmanLayerDNSTKEBudget
    use kind_parameters, only: rkind, clen
    use igrid_Operators, only: igrid_ops
    use DerivativesMod, only: derivatives
    use constants, only: pi, two
    use mpi
    use decomp_2d
    use timer, only: tic, toc
    use exits, only: message
    implicit none

    real(rkind), dimension(:,:,:,:), allocatable :: Rij
    real(rkind), dimension(:,:,:), allocatable :: &
         &buff1, buff2, buff3,& 
         &u,v,w,p, umean,vmean,pmean, ufluct,vfluct,pfluct, MKE
    real(rkind), dimension(:,:), allocatable :: &
         & Prod, Diss, Diff, PDif, RDif,&
         & umean_1d, vmean_1d, MKE_1d
    real(rkind) :: &
            & dx, dy, dz , &
            &Re=400.d0, Lx=26.d0, Ly=26.d0, Lz=24.d0
    integer :: & 
        &nx, ny, nz, nt, RunID, TIDX, & 
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
    call ops%allocate3DField(u)
    call ops%allocate3DField(v)
    call ops%allocate3DField(w)
    call ops%allocate3DField(p)
    call ops%allocate3DField(ufluct)
    call ops%allocate3DField(vfluct)
    call ops%allocate3DField(pfluct)

    ! Get number of timesteps and allocate arrays
    nt = (tstop - tstart)/tstep
    allocate(umean(nx,nx,nz))
    allocate(vmean(nx,nx,nz))
    allocate(pmean(nx,nx,nz))
    allocate(MKE(nx,nx,nz))
    allocate(Rij(nx,nx,nz,6))
    allocate(Prod(nz,1))
    allocate(Diss(nz,1))
    allocate(Diff(nz,1))
    allocate(PDif(nz,1))
    allocate(RDif(nz,1))
    allocate(umean_1d(nz,1))
    allocate(vmean_1d(nz,1))
    allocate(MKE_1d(nz,1))

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

        ! Fluctuations and means as 3D fields
        call ops%getFluct_from_MeanZ(u,ufluct)
        call ops%getFluct_from_MeanZ(v,vfluct)
        call ops%getFluct_from_MeanZ(p,pfluct)
        umean = umean + u-ufluct
        vmean = vmean + v-vfluct
        pmean = pmean + p-pfluct
        
        ! Reynolds stresses
        Rij(:,:,:,1) = Rij(:,:,:,1) + ufluct*ufluct
        Rij(:,:,:,2) = Rij(:,:,:,2) + ufluct*vfluct
        Rij(:,:,:,3) = Rij(:,:,:,3) + ufluct*w
        Rij(:,:,:,4) = Rij(:,:,:,4) + vfluct*vfluct
        Rij(:,:,:,5) = Rij(:,:,:,5) + vfluct*w
        Rij(:,:,:,6) = Rij(:,:,:,6) + w*w 
      
        call toc()
        idx = idx + 1
    end do 

    ! Averages in time, 3D fields
    umean = umean/nt
    vmean = vmean/nt
    pmean = pmean/nt
    Rij = Rij/nt
    MKE = 0.5*(umean**2 + vmean**2)

    ! Advection Uj*DE/dxj = 0
    ! call ops%getGradient(MKE,buff1,buff2,buff3,0,0)
    ! call ops%TakeMean_xy(umean*buff1+vmean*buff2+wmean*buff3, Advc)

    ! Pressure diffusion d/dxj(PUj) = Uj*dP/dxj= 0
    call ops%getGradient(pmean,buff1,buff2,buff3,0,0)
    call ops%TakeMean_xy(umean*buff1+vmean*buff2, PDif)

    ! Viscous diffusion d2/dxj2 MKE
    call ops%d2dz2(MKE,buff1,0,0)
    call ops%TakeMean_xy(buff1,Diff)

    ! Dissipation dUi/dxj*dUi/dxj = dUdz^2 + dVdz^2
    call ops%ddz(umean,buff1,0,1)
    call ops%ddz(vmean,buff2,0,1)
    call ops%TakeMean_xy( buff1**2 + buff2**2, Diss)

    ! Production uiuj*dUi/dxj = uwdUdz + vwdVdz
    call ops%TakeMean_xy(Rij(:,:,:,3)*buff1+Rij(:,:,:,5)*buff2, Prod)

    ! Reynolds diffusion terms d/dxj(Ui uiuj) = ddz(Uuw+Vvw)
    call ops%ddz( umean*Rij(:,:,:,3)+vmean*Rij(:,:,:,5),buff1,0,0)
    call ops%TakeMean_xy(buff1,RDif)

    ! Coriolis terms: omega3*(umean*G2-vmean*G1)
    ! 1D Averages
    call ops%TakeMean_xy(umean,umean_1d)
    call ops%TakeMean_xy(vmean,vmean_1d)
    call ops%TakeMean_xy(MKE,MKE_1d)
    
    ! Write out time averages
    if (nrank == 0) then
        call message(0,"Writing  files...")
        call ops%WriteASCII_2D(Prod, "Prod")
        call ops%WriteASCII_2D(Diff, "Diff")
        call ops%WriteASCII_2D(Diff, "PDif")
        call ops%WriteASCII_2D(Diff, "RDif")
        call ops%WriteASCII_2D(Diss, "Diss")
        call ops%WriteASCII_2D(umean_1d, "u_mn")
        call ops%WriteASCII_2D(vmean_1d, "v_mn")
        call ops%WriteASCII_2D(MKE_1d, "MKE_")
    end if 
        
    call ops%destroy()
    call MPI_Finalize(ierr)           
    call message(0,"Done")

end program 
