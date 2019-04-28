program EkmanLayerDNSTurbStats
    use kind_parameters, only: rkind, clen
    use igrid_Operators, only: igrid_ops
    use constants, only: pi, two
    use mpi
    use decomp_2d, only: nrank 
    use timer, only: tic, toc
    use exits, only: message
    implicit none

    real(rkind), dimension(:,:,:), allocatable ::buff1,buff2,buff3, u,v,w,ufluct,vfluct,omegax,omegay,omegaz
    real(rkind) :: time, dx, dy, dz, Re=400.d0, Lx=26.d0, Ly=26.d0, Lz=24.d0
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
    call ops%allocate3DField(u)
    call ops%allocate3DField(v)
    call ops%allocate3DField(w)
    call ops%allocate3DField(ufluct)
    call ops%allocate3DField(vfluct)
    call ops%allocate3DField(omegax)
    call ops%allocate3DField(omegay)
    call ops%allocate3DField(omegaz)

    ! Get number of timesteps and allocate arrays
    nt = (tstop - tstart)/tstep
 
    ! Compute for each timestep
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

        ! Fluctuations and curl
        call ops%getCurl(u, v, w, buff1, buff2, buff3,&
                & botBC, topBC, botBC, topBC)
        call ops%getFluct_from_MeanZ(buff1,omegax)
        call ops%getFluct_from_MeanZ(buff2,omegay)
        call ops%getFluct_from_MeanZ(buff3,omegaz)
        call message(0,"max omegax:",maxval(omegax))

        ! Write out files
        call ops%WriteField3D(omegax, "om_x", TIDX)
        call ops%WriteField3D(omegay, "om_y", TIDX)
        call ops%WriteField3D(omegaz, "om_z", TIDX)
        call toc()
        idx = idx + 1
    end do 
        
    call ops%destroy()
    call MPI_Finalize(ierr)           

end program 

