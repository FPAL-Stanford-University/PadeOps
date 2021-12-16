program test_gaborMode
    use domainSetup, only: xDom, yDom, zDom
    use gaborModeMod, only : gaborMode, useStrain 
    use kind_parameters, only: rkind, clen 
    use EulerG_mod, only: EulerG 
    use constants 
    use decomp_2d 
    use timer, only: tic, toc 
    use exits, only: message

    integer :: nx = 32, ny = 32, nz = 32
    real(rkind), dimension(2), parameter :: DomX = [0.d0, two*pi], DomY = [0.d0, two*pi], DomZ = [0.d0, two*pi]
    type(decomp_info) :: gp
    integer, parameter :: levels = 4
    real(rkind), dimension(:,:,:), pointer :: uF, vF, wF 
    type(gaborMode) :: gm
    type(EulerG) :: Egrid
    real(rkind) :: xRange(2), yRange(2), zRange(2), delta(3)
    complex(kind=8), parameter :: uhat = (0.25d0,2.d0), vhat = (-2.42d0,1.5d0), what = (0.35d0,-2.65d0)
    real(kind=4), dimension(:,:,:), allocatable :: u_sp, v_sp, w_sp 
    character(len=clen) :: datadir = '/work2/06632/ryanhass/stampede2/Enrichment/Gabor_class_tests/Vadvect'
    character(len=18) :: fname
    integer :: n, ierr
    
    ! Test with these inputs  
    real(rkind), parameter :: wSupport = 1.5d0  
    real(rkind), parameter :: kx = 26.d0, ky = 26.d0, kz = 40.d0 
    real(rkind), parameter :: xMode = 1.d0, yMode = 1.d0, zMode = 1.d0   
    integer, parameter :: modeAtLevel = 4
    real(rkind) :: ules = 0.d0, vles = two*pi, wles = 0.d0
    real(rkind) :: dudx = 0.d0, dudy = 2.d0, dudz = 0.d0 ! Shear
    real(rkind) :: dvdx = 0.d0, dvdy = 0.d0, dvdz = 0.d0
    real(rkind) :: dwdx = 0.d0, dwdy = 0.d0, dwdz = 0.d0
    real(rkind) :: dt = 0.01, epsKE = 0.d0, kmin = 11.d0
    integer :: nsteps = 100
    character(len=clen) :: mssg

    call MPI_Init(ierr)

    call decomp_2d_init(nx+1, ny+1, nz+1, 0, 0)
    call get_decomp_info(gp)
   
    xDom = DomX; yDom = DomY; zDom = DomZ

    call Egrid%init( gp, levels)
    call Egrid%getPointersToFields(uF, vF, wF, modeAtLevel)
    call Egrid%getDomainRangeForLevel(xRange, yRange, zRange, modeAtLevel)
    call Egrid%getDeltaForLevel(delta, modeAtLevel)
    allocate(u_sp(size(uF,1), size(uF,2), size(uF,3)))
    allocate(v_sp(size(uF,1), size(uF,2), size(uF,3)))
    allocate(w_sp(size(uF,1), size(uF,2), size(uF,3)))

    ! Do gabor mode stuff     
    call gm%init(uhat, vhat, what, xMode, yMode, zMode, kx, ky, kz, wSupport)

    call tic
    ! Render 
    call gm%render(u_sp, v_sp, w_sp, xRange, yRange, zRange, delta)

    ! Before agglomerating, cast single precision fields to double precision 
    uF = real(u_sp, rkind)
    vF = real(v_sp, rkind)
    wF = real(w_sp, rkind)

    ! See what it looks like on the finest level and mode's level   
    call Egrid%agglomerate()
    call toc

    call Egrid%writeFields(trim(datadir)//"/Run01_init", modeAtLevel)

    if (modeAtLevel .ne. levels) then 
        call Egrid%writeFields(trim(datadir)//"/Run01_init", levels)
    end if
    
    ! Zero velocity field
    u_sp = 0.e0
    v_sp = 0.e0
    w_sp = 0.e0 
  
    ! --- Evolve the mode (Simple advection) --- !

    call Egrid%ResetVelocityToZero()

    useStrain = .false.
    do n = 1,nsteps
        call gm%evolve(ules,vles,wles,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,epsKE,dt,kmin)
        if (n == nsteps/2) then
          call gm%render(u_sp, v_sp, w_sp, xRange, yRange, zRange, delta)
          uF = real(u_sp, rkind)
          vF = real(v_sp, rkind)
          wF = real(w_sp, rkind)
          call Egrid%agglomerate()
          write(fname,'(A,I3.3)')'Run01_adv_step',n
          call Egrid%writeFields(trim(datadir)//'/'//trim(fname), modeAtLevel)
          call Egrid%ResetVelocityToZero()
    
          ! Zero velocity field
          u_sp = 0.e0
          v_sp = 0.e0
          w_sp = 0.e0 
        end if
        if (mod(n,10) == 0) then
          write(mssg,'(A,I4,A,F10.5)')"Time step ", n, '; gm%uhatR: ', gm%uhatR
          call message(trim(mssg))
        end if
    end do
    
    call tic
    
    ! Render 
    write(mssg,'(A,F5.2,A,F5.2,A,F5.2,A)')'mode location: (', gm%x,',',gm%y,',',gm%z,')'
    call message(trim(mssg))
    call gm%render(u_sp, v_sp, w_sp, xRange, yRange, zRange, delta)

    ! Before agglomerating, cast single precision fields to double precision 
    uF = real(u_sp, rkind)
    vF = real(v_sp, rkind)
    wF = real(w_sp, rkind)

    ! See what it looks like on the finest level and mode's level   
    call Egrid%agglomerate()
    call toc

    call Egrid%writeFields(trim(datadir)//"/Run01_adv", modeAtLevel)

    if (modeAtLevel .ne. levels) then 
        call Egrid%writeFields(trim(datadir)//"/Run01_adv", levels)
    end if 
    
    ! Zero velocity field
    u_sp = 0.e0
    v_sp = 0.e0
    w_sp = 0.e0 

    ! --- Evolve the mode (w/ straining) --- !

    call Egrid%ResetVelocityToZero()
    
    useStrain = .true.
    do n = 1,nsteps
        call gm%evolve(ules,vles,wles,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,epsKE,dt,kmin)
        if (n == nsteps/2) then
          call gm%render(u_sp, v_sp, w_sp, xRange, yRange, zRange, delta)
          uF = real(u_sp, rkind)
          vF = real(v_sp, rkind)
          wF = real(w_sp, rkind)
          call Egrid%agglomerate()
          write(fname,'(A,I3.3)')'Run01_strn_step',n
          call Egrid%writeFields(trim(datadir)//'/'//trim(fname), modeAtLevel)
          call Egrid%ResetVelocityToZero()
    
          ! Zero velocity field
          u_sp = 0.e0
          v_sp = 0.e0
          w_sp = 0.e0 
        end if
        if (mod(n,10) == 0) then
          write(mssg,'(A,I4,A,F10.5)')"Time step ", n, '; gm%uhatR: ', gm%uhatR
          call message(trim(mssg))
        end if
    end do
    
    call tic
    ! Render 
    write(mssg,'(A,F5.2,A,F5.2,A,F5.2,A)')'mode location: (', gm%x,',',gm%y,',',gm%z,')'
    call message(trim(mssg))
    call gm%render(u_sp, v_sp, w_sp, xRange, yRange, zRange, delta)

    ! Before agglomerating, cast single precision fields to double precision 
    uF = real(u_sp, rkind)
    vF = real(v_sp, rkind)
    wF = real(w_sp, rkind)

    ! See what it looks like on the finest level and mode's level   
    call Egrid%agglomerate()
    call toc

    call Egrid%writeFields(trim(datadir)//"/Run01_strn", modeAtLevel)

    if (modeAtLevel .ne. levels) then 
        call Egrid%writeFields(trim(datadir)//"/Run01_strn", levels)
    end if

!BEGIN_DEBUG****************************************************************
write(mssg,'(A,I3,A,I3,A)')'(yst,yen) = (',Egrid%fields(4)%gp%xst(2),',',Egrid%fields(4)%gp%xen(2),')'
print*, trim(mssg)
write(mssg,'(A,F20.16,A,F20.16,A)')'Y domain = [',yRange(1),',',yRange(2),']'
print*, trim(mssg)
!END_DEBUG******************************************************************
    call Egrid%destroy()
    call decomp_2d_finalize
    call MPI_Finalize(ierr)
end program 
