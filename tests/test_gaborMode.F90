program test_gaborMode
    use gaborModeMod, only : gaborMode 
    use kind_parameters, only: rkind, clen 
    use EulerG_mod, only: EulerG 
    use constants 
    use decomp_2d 
    use timer, only: tic, toc  

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
    character(len=clen) :: datadir = '/work2/06632/ryanhass/stampede2/Enrichment/Gabor_class_tests'
    
    ! Test with these inputs  
    real(rkind), parameter :: wSupport = 0.5d0  
    real(rkind), parameter :: kx = 26.d0, ky = 18.d0, kz = 40.d0 
    real(rkind), parameter :: xMode = 1.d0, yMode = 1.d0, zMode = 1.d0   
    integer, parameter :: modeAtLevel = 4

    call MPI_Init(ierr)

    call decomp_2d_init(nx+1, ny+1, nz+1, 0, 0)
    call get_decomp_info(gp)
   
    call Egrid%init( gp, DomX, DomY, DomZ, levels)
    call Egrid%getPointersToFields(uF, vF, wF, modeAtLevel)
    call Egrid%getDomainRangeForLevel(xRange, yRange, zRange, modeAtLevel)
    call Egrid%getDeltaForLevel(delta, modeAtLevel)
    allocate(u_sp(size(uF,1), size(uF,2), size(uF,3)))
    allocate(v_sp(size(uF,1), size(uF,2), size(uF,3)))
    allocate(w_sp(size(uF,1), size(uF,2), size(uF,3)))

    ! Do gabor mode stuff     
    call gm%init(uhat, vhat, what, xMode, yMode, zMode, kx, ky, kz, wSupport)
  
    ! Evolve if needed ...

    ! call Egrid%ResetVelocityToZero() << depends on whether velocity has been rendered before 

    call tic()
    ! Render 
    call gm%render(u_sp, v_sp, w_sp, xRange, yRange, zRange, delta)

    ! Before agglomerating, cast single precision fields to double precision 
    uF = real(u_sp, rkind)
    vF = real(v_sp, rkind)
    wF = real(w_sp, rkind)

    ! See what it looks like on the finest level and mode's level   
    call Egrid%agglomerate()
    call toc

    call Egrid%writeFields(trim(datadir)//"/Run01", modeAtLevel)

    if (modeAtLevel .ne. levels) then 
        call Egrid%writeFields(trim(datadir)//"/Run01", levels)
    end if 

end program 
