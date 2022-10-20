program test_gaborMode
    use domainSetup, only: get_window_function_bounds, xDom, yDom, zDom
    use gaborModeMod, only : gaborMode, useStrain 
    use kind_parameters, only: rkind, clen 
    use EulerG_mod, only: EulerG 
    use constants 
    use decomp_2d
    use mpi
    use timer, only: tic, toc 
    use exits, only: message, gracefulExit

    integer :: nx = 32, ny = 32, nz = 32
    integer, dimension(3) :: ng
    real(rkind), dimension(2), parameter :: DomX = [0.d0, two*pi], DomY = [0.d0, two*pi], DomZ = [0.d0, two*pi]
    type(decomp_info) :: gp, gpGM
    integer :: nxsupp, nysupp, nzsupp
    logical, dimension(3), parameter :: periodic_bcs = [.true.,.true.,.true.]
    integer, parameter :: levels = 4
    real(rkind), dimension(:,:,:), pointer :: uF, vF, wF 
    type(gaborMode) :: gm
    type(EulerG) :: Egrid
    real(rkind) :: xRange(2), yRange(2), zRange(2), delta(3)
    complex(kind=8), parameter :: uhat = (0.25d0,2.d0), vhat = (-2.42d0,1.5d0), what = (0.35d0,-2.65d0)
    real(kind=4), dimension(:,:,:), allocatable :: u_sp, v_sp, w_sp
    real(rkind), dimension(:,:,:), allocatable :: haloVel
    character(len=clen) :: datadir = '/work2/06632/ryanhass/stampede2/Enrichment/Gabor_class_tests/Vadvect'
    character(len=18) :: fname
    integer :: n, ierr
    
    ! Test with these inputs  
    real(rkind) :: wSupport 
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

    integer :: istg, ieng, jstg, jeng, kstg, keng

    call MPI_Init(ierr)

    call get_window_function_bounds(DomX,DomY,DomZ,nx,ny,nz,modeAtLevel,2,2,2,nxsupp,nysupp,nzsupp,wSupport)
    xDom = DomX; yDom = DomY; zDom = DomZ
   
    call decomp_2d_init(nx+1, ny+1, nz+1, 0, 0, periodic_bcs)
    call get_decomp_info(gp)
print*, "A"
    call Egrid%init( gp, levels)
    call Egrid%getPointersToFields(uF, vF, wF, modeAtLevel)
    call Egrid%getDomainRangeForLevel(xRange, yRange, zRange, modeAtLevel)
    call Egrid%getDeltaForLevel(delta, modeAtLevel)

print*, "B"
    ! Define Gabor mode grid partition and velocity arrays. Eventually this will
    ! be absorbed in the constructor for the Gabor mode population class
    call Egrid%fields(modeAtLevel)%getMeshSize(ng)
    call decomp_info_init(ng(1),ng(2),ng(3),gpGM)
print*, "C"

    ! Define global index range for velocity fields
    istg = gpGM%xst(1)-nxsupp/2
    ieng = gpGM%xen(1)+nxsupp/2
    
    jstg = gpGM%xst(2)-nysupp/2
    jeng = gpGM%xen(2)+nysupp/2

    kstg = gpGM%xst(3)-nzsupp/2
    keng = gpGM%xen(3)+nzsupp/2

    ! Allocate memory for Gabor-induced velocity
    allocate(u_sp(istg:ieng,jstg:jeng,kstg:keng))
    allocate(v_sp(istg:ieng,jstg:jeng,kstg:keng))
    allocate(w_sp(istg:ieng,jstg:jeng,kstg:keng))
    !allocate(u_sp(size(uF,1)+nxsupp, size(uF,2)+nysupp, size(uF,3)+nzsupp))
    !allocate(v_sp(size(uF,1)+nxsupp, size(uF,2)+nysupp, size(uF,3)+nzsupp))
    !allocate(w_sp(size(uF,1)+nxsupp, size(uF,2)+nysupp, size(uF,3)+nzsupp))

    ! Do gabor mode stuff     
    call gm%init(uhat, vhat, what, xMode, yMode, zMode, kx, ky, kz, wSupport)
print*, "D"

    call tic
    ! Render
    ! modify Ranges to account for halo regions
    xRange(1) = xRange(1)-wSupport/2
    xRange(2) = xRange(2)+wSupport/2
     
    yRange(1) = yRange(1)-wSupport/2
    yRange(2) = yRange(2)+wSupport/2
     
    zRange(1) = zRange(1)-wSupport/2
    zRange(2) = zRange(2)+wSupport/2
     
print*, "D2"
    call gm%render(u_sp, v_sp, w_sp, xRange, yRange, zRange, delta)

print*, "E"
print*, "shape(u_sp) = ", shape(u_sp)
    ! Now exchange halo regions
    call update_halo(real(u_sp,rkind),haloVel,nxsupp/2)
print*, "F"
    ! y-halo contribution
    u_sp(:,gpGM%xst(2):gpGM%xst(2)+nysupp/2,:) = &
      u_sp(:,gpGM%xst(2):gpGM%xst(2)+nysupp/2,:) + real(haloVel(:,jstg-nysupp/2:jstg-1,kstg:keng),kind=4)
    u_sp(:,gpGM%xen(2)-nysupp/2+1:gpGM%xen(2)+nysupp/2,:) = &
      u_sp(:,gpGM%xen(2)-nysupp/2+1:gpGM%xen(2),:) + real(haloVel(:,jeng+1:jeng+nysupp/2,kstg:keng),kind=4)

    ! z-halo contribution
    u_sp(:,:,gpGM%xst(3):gpGM%xst(3)+nzsupp/2) = &
      u_sp(:,:,gpGM%xst(3):gpGM%xst(3)+nzsupp/2) + real(haloVel(:,jstg:jeng,kstg-nzsupp/2:kstg-1),kind=4)
    u_sp(:,:,gpGM%xen(3)-nzsupp/2+1:gpGM%xen(3)+nzsupp/2) = &
      u_sp(:,:,gpGM%xen(3)-nzsupp/2+1:gpGM%xen(3)) + real(haloVel(:,jstg:jeng,keng+1:keng+nzsupp/2),kind=4)

    ! Add x-periodic contribution
    u_sp(1:nxsupp/2,:,:) = u_sp(ng(1)+1:ng(1)+nxsupp/2,:,:) + u_sp(1:nxsupp/2,:,:)
    u_sp(ng(1)-nxsupp/2+1:ng(1),:,:) = u_sp(-nxsupp/2+1:0,:,:) + u_sp(ng(1)-nxsupp/2+1:ng(1),:,:)

    call update_halo(real(v_sp,rkind),haloVel,nxsupp/2)
    ! y-halo contribution
    v_sp(:,gpGM%xst(2):gpGM%xst(2)+nysupp/2,:) = &
      v_sp(:,gpGM%xst(2):gpGM%xst(2)+nysupp/2,:) + real(haloVel(:,jstg-nysupp/2:jstg-1,kstg:keng),kind=4)
    v_sp(:,gpGM%xen(2)-nysupp/2+1:gpGM%xen(2)+nysupp/2,:) = &
      v_sp(:,gpGM%xen(2)-nysupp/2+1:gpGM%xen(2),:) + real(haloVel(:,jeng+1:jeng+nysupp/2,kstg:keng),kind=4)

    ! z-halo contribution
    v_sp(:,:,gpGM%xst(3):gpGM%xst(3)+nzsupp/2) = &
      v_sp(:,:,gpGM%xst(3):gpGM%xst(3)+nzsupp/2) + real(haloVel(:,jstg:jeng,kstg-nzsupp/2:kstg-1),kind=4)
    v_sp(:,:,gpGM%xen(3)-nzsupp/2+1:gpGM%xen(3)+nzsupp/2) = &
      v_sp(:,:,gpGM%xen(3)-nzsupp/2+1:gpGM%xen(3)) + real(haloVel(:,jstg:jeng,keng+1:keng+nzsupp/2),kind=4)

    ! Add x-periodic contribution
    v_sp(1:nxsupp/2,:,:) = v_sp(ng(1)+1:ng(1)+nxsupp/2,:,:) + v_sp(1:nxsupp/2,:,:)
    v_sp(ng(1)-nxsupp/2+1:ng(1),:,:) = v_sp(-nxsupp/2+1:0,:,:) + v_sp(ng(1)-nxsupp/2+1:ng(1),:,:)
    
    call update_halo(real(w_sp,rkind),haloVel,nxsupp/2)
print*, "G"
    ! y-halo contribution
    w_sp(:,gpGM%xst(2):gpGM%xst(2)+nysupp/2,:) = &
      w_sp(:,gpGM%xst(2):gpGM%xst(2)+nysupp/2,:) + real(haloVel(:,jstg-nysupp/2:jstg-1,kstg:keng),kind=4)
    w_sp(:,gpGM%xen(2)-nysupp/2+1:gpGM%xen(2)+nysupp/2,:) = &
      w_sp(:,gpGM%xen(2)-nysupp/2+1:gpGM%xen(2),:) + real(haloVel(:,jeng+1:jeng+nysupp/2,kstg:keng),kind=4)

    ! z-halo contribution
    w_sp(:,:,gpGM%xst(3):gpGM%xst(3)+nzsupp/2) = &
      w_sp(:,:,gpGM%xst(3):gpGM%xst(3)+nzsupp/2) + real(haloVel(:,jstg:jeng,kstg-nzsupp/2:kstg-1),kind=4)
    w_sp(:,:,gpGM%xen(3)-nzsupp/2+1:gpGM%xen(3)+nzsupp/2) = &
      w_sp(:,:,gpGM%xen(3)-nzsupp/2+1:gpGM%xen(3)) + real(haloVel(:,jstg:jeng,keng+1:keng+nzsupp/2),kind=4)

    ! Add x-periodic contribution
    w_sp(1:nxsupp/2,:,:) = w_sp(ng(1)+1:ng(1)+nxsupp/2,:,:) + w_sp(1:nxsupp/2,:,:)
    w_sp(ng(1)-nxsupp/2+1:ng(1),:,:) = w_sp(-nxsupp/2+1:0,:,:) + w_sp(ng(1)-nxsupp/2+1:ng(1),:,:)
    
    ! Before agglomerating, cast single precision fields to double precision
    uF = real(u_sp(gpGM%xst(1):gpGM%xen(1),gpGM%xst(2):gpGM%xen(2),gpGM%xst(3):gpGM%xen(3)), rkind)
    vF = real(v_sp(gpGM%xst(1):gpGM%xen(1),gpGM%xst(2):gpGM%xen(2),gpGM%xst(3):gpGM%xen(3)), rkind)
    wF = real(w_sp(gpGM%xst(1):gpGM%xen(1),gpGM%xst(2):gpGM%xen(2),gpGM%xst(3):gpGM%xen(3)), rkind)

    ! See what it looks like on the finest level and mode's level   
    call Egrid%agglomerate()
print*, "H"
    call toc

    call Egrid%writeFields(trim(datadir)//"/Run01_init", modeAtLevel)

print*, "I"
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
