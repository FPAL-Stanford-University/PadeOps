module SolidGrid
    use kind_parameters, only: rkind, clen
    use constants,       only: zero,eps,third,half,one,two,three,four,pi
    use FiltersMod,      only: filters
    use GridMod,         only: grid
    use gridtools,       only: alloc_buffs, destroy_buffs
    use sgrid_hooks,     only: meshgen, initfields, hook_output, hook_bc, hook_timestep, hook_mixture_source
    use decomp_2d,       only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                               transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,  only: derivatives
    use LADMod,          only: ladobject
    use StiffGasEOS,     only: stiffgas
    use Sep1SolidEOS,    only: sep1solid
    use SolidMixtureMod, only: solid_mixture
   
    implicit none

    integer, parameter :: rho_index    = 1 
    integer, parameter :: u_index      = 2
    integer, parameter :: v_index      = 3
    integer, parameter :: w_index      = 4
    integer, parameter :: p_index      = 5
    integer, parameter :: T_index      = 6
    integer, parameter :: e_index      = 7
    integer, parameter :: sos_index    = 8
    integer, parameter :: mu_index     = 9
    integer, parameter :: bulk_index   = 10
    integer, parameter :: kap_index    = 11
    integer, parameter :: sxx_index    = 12
    integer, parameter :: sxy_index    = 13
    integer, parameter :: sxz_index    = 14
    integer, parameter :: syy_index    = 15
    integer, parameter :: syz_index    = 16
    integer, parameter :: szz_index    = 17

    integer, parameter :: nfields = 17

    integer, parameter :: mom_index = 1
    integer, parameter :: TE_index = mom_index+3
    integer, parameter :: ncnsrv  = TE_index

    ! These indices are for data management, do not change if you're not sure of what you're doing
    integer, parameter :: tauxyidx = 2
    integer, parameter :: tauxzidx = 3
    integer, parameter :: tauyzidx = 6
    integer, parameter :: tauxxidx = 4
    integer, parameter :: tauyyidx = 7
    integer, parameter :: tauzzidx = 8

    integer, parameter :: qxidx = 1
    integer, parameter :: qyidx = 5
    integer, parameter :: qzidx = 9


    ! Number of buffers to create
    integer, parameter :: nbufsx = 2
    integer, parameter :: nbufsy = 6
    integer, parameter :: nbufsz = 2

    type, extends(grid) :: sgrid
       
        type(filters),       allocatable :: gfil
        type(solid_mixture), allocatable :: mix
        type(ladobject),     allocatable :: LAD

        logical     :: plastic
        logical     :: explPlast
        real(rkind) :: tau0
        real(rkind) :: invtau0

        real(rkind), dimension(:,:,:,:), allocatable :: Wcnsrv                               ! Conserved variables
        real(rkind), dimension(:,:,:,:), allocatable :: xbuf, ybuf, zbuf   ! Buffers
       
        ! real(rkind) :: Cmu, Cbeta, Ckap, Cdiff, CY

        real(rkind) :: rho0

        real(rkind), dimension(:,:,:), pointer :: x 
        real(rkind), dimension(:,:,:), pointer :: y 
        real(rkind), dimension(:,:,:), pointer :: z 
        
        real(rkind), dimension(:,:,:), pointer :: rho 
        real(rkind), dimension(:,:,:), pointer :: u 
        real(rkind), dimension(:,:,:), pointer :: v 
        real(rkind), dimension(:,:,:), pointer :: w 
        real(rkind), dimension(:,:,:), pointer :: p 
        real(rkind), dimension(:,:,:), pointer :: T 
        real(rkind), dimension(:,:,:), pointer :: e 
        real(rkind), dimension(:,:,:), pointer :: sos 
        real(rkind), dimension(:,:,:), pointer :: mu 
        real(rkind), dimension(:,:,:), pointer :: bulk 
        real(rkind), dimension(:,:,:), pointer :: kap
        real(rkind), dimension(:,:,:), pointer :: tauaiidivu

        real(rkind), dimension(:,:,:,:), pointer :: devstress
        real(rkind), dimension(:,:,:), pointer :: sxx
        real(rkind), dimension(:,:,:), pointer :: sxy
        real(rkind), dimension(:,:,:), pointer :: sxz
        real(rkind), dimension(:,:,:), pointer :: syy
        real(rkind), dimension(:,:,:), pointer :: syz
        real(rkind), dimension(:,:,:), pointer :: szz
        
        contains
            procedure          :: init
            procedure          :: destroy
            procedure          :: laplacian
            procedure          :: gradient 
            procedure          :: advance_RK45
            procedure          :: simulate
            procedure, private :: get_dt
            procedure, private :: get_primitive
            procedure, private :: get_conserved
            procedure, private :: post_bc
            procedure, private :: getRHS
            procedure, private :: getRHS_x
            procedure, private :: getRHS_y
            procedure, private :: getRHS_z
            ! procedure          :: getLAD
            procedure          :: filter
            procedure          :: getPhysicalProperties
            procedure, private :: get_tau
            procedure, private :: get_q
            ! procedure, private :: getPlasticSources
    end type

contains
    subroutine init(this, inputfile )
        use reductions, only: P_MAXVAL
        use exits,      only: message, warning, nancheck, GracefulExit
        class(sgrid),target, intent(inout) :: this
        character(len=clen), intent(in) :: inputfile  

        integer :: nx, ny, nz
        integer :: ns = 1
        character(len=clen) :: outputdir
        character(len=clen) :: inputdir
        character(len=clen) :: vizprefix = "sgrid"
        real(rkind) :: tviz = zero
        character(len=clen), dimension(nfields) :: varnames
        logical :: periodicx = .true. 
        logical :: periodicy = .true. 
        logical :: periodicz = .true.
        character(len=clen) :: derivative_x = "cd10"  
        character(len=clen) :: derivative_y = "cd10" 
        character(len=clen) :: derivative_z = "cd10"
        character(len=clen) :: filter_x = "cf90"  
        character(len=clen) :: filter_y = "cf90" 
        character(len=clen) :: filter_z = "cf90"
        integer :: prow = 0, pcol = 0 
        integer :: i, j, k 
        integer :: ioUnit
        real(rkind) :: gam = 1.4_rkind
        real(rkind) :: Rgas = one
        real(rkind) :: PInf = zero
        real(rkind) :: shmod = zero
        real(rkind) :: rho0 = one
        integer :: nsteps = -1
        real(rkind) :: dt = -one
        real(rkind) :: tstop = one
        real(rkind) :: CFL = -one
        logical :: SkewSymm = .FALSE.
        real(rkind) :: Cmu = 0.002_rkind
        real(rkind) :: Cbeta = 0.5_rkind
        real(rkind) :: Ckap = 0.01_rkind
        real(rkind) :: Cdiff = 0.003_rkind
        real(rkind) :: CY = 100._rkind
        logical     :: plastic = .FALSE.
        real(rkind) :: yield = real(1.D30,rkind)
        logical     :: explPlast = .FALSE.
        real(rkind) :: tau0 = one

        !real(rkind), dimension(:,:,:,:), allocatable :: finger, fingersq
        !real(rkind), dimension(:,:,:),   allocatable :: trG, trG2, detG

        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                             inputdir, outputdir, vizprefix, tviz, &
                                  periodicx, periodicy, periodicz, &
                         derivative_x, derivative_y, derivative_z, &
                                     filter_x, filter_y, filter_z, &
                                                       prow, pcol, &
                                                         SkewSymm  
        namelist /SINPUT/  gam, Rgas, PInf, shmod, rho0, plastic, yield, &
                           explPlast, tau0, ns, Cmu, Cbeta, Ckap, Cdiff, CY

        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=INPUT)
        read(unit=ioUnit, NML=SINPUT)
        close(ioUnit)

        this%nx = nx
        this%ny = ny
        this%nz = nz

        this%tsim = zero
        this%tstop = tstop
        this%dtfixed = dt
        this%dt = dt
        this%CFL = CFL

        this%step = 0
        this%nsteps = nsteps

        this%rho0 = rho0

        ! this%Cmu = Cmu
        ! this%Cbeta = Cbeta
        ! this%Ckap = Ckap
        ! this%Cdiff = Cdiff
        ! this%CY = CY

        this%plastic = plastic
        
        this%explPlast = explPlast
        this%tau0 = tau0

        ! Allocate decomp
        if ( allocated(this%decomp) ) deallocate(this%decomp)
        allocate(this%decomp)
        
        ! Initialize decomp
        call decomp_2d_init(nx, ny, nz, prow, pcol)
        call get_decomp_info(this%decomp)
        
        ! Set all the attributes of the abstract grid type         
        this%outputdir = outputdir 
        
        this%periodicx = periodicx
        this%periodicy = periodicy
        this%periodicz = periodicz

        this%derivative_x = derivative_x    
        this%derivative_y = derivative_y    
        this%derivative_z = derivative_z  

        this%filter_x = filter_x    
        this%filter_y = filter_y    
        this%filter_z = filter_z  
        
        ! Finally, set the local array dimensions
        this%nxp = this%decomp%ysz(1)
        this%nyp = this%decomp%ysz(2)
        this%nzp = this%decomp%ysz(3)

        ! Allocate mesh
        if ( allocated(this%mesh) ) deallocate(this%mesh) 
        call alloc_buffs(this%mesh,3,'y',this%decomp)
        
        ! Associate pointers for ease of use
        this%x    => this%mesh  (:,:,:, 1) 
        this%y    => this%mesh  (:,:,:, 2) 
        this%z    => this%mesh  (:,:,:, 3)

        ! Generate default mesh: X \in [-1, 1), Y \in [-1, 1), Z \in [-1, 1)
        this%dx = two/nx
        this%dy = two/ny
        this%dz = two/nz

        ! Generate default mesh 
        do k = 1,size(this%mesh,3)
            do j = 1,size(this%mesh,2)
                do i = 1,size(this%mesh,1)
                    this%mesh(i,j,k,1) = -one + (this%decomp%yst(1) - 1 + i - 1)*this%dx           
                    this%mesh(i,j,k,2) = -one + (this%decomp%yst(2) - 1 + j - 1)*this%dy           
                    this%mesh(i,j,k,3) = -one + (this%decomp%yst(3) - 1 + k - 1)*this%dz           
                end do 
            end do 
        end do  

        ! Go to hooks if a different mesh is desired 
        call meshgen(this%decomp, this%dx, this%dy, this%dz, this%mesh) 

        ! Allocate der
        if ( allocated(this%der) ) deallocate(this%der)
        allocate(this%der)

        ! Initialize derivatives 
        call this%der%init(                           this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
                      derivative_x,  derivative_y,   derivative_z, &
                           .false.,       .false.,        .false., &
                           .false.)      

        ! Allocate fil and gfil
        if ( allocated(this%fil) ) deallocate(this%fil)
        allocate(this%fil)
        if ( allocated(this%gfil) ) deallocate(this%gfil)
        allocate(this%gfil)
        
        ! Initialize filters
        call this%fil%init(                           this%decomp, &
                         periodicx,     periodicy,      periodicz, &
                          filter_x,      filter_y,       filter_z  )      
        call this%gfil%init(                          this%decomp, &
                         periodicx,     periodicy,      periodicz, &
                        "gaussian",    "gaussian",     "gaussian"  )      

        ! Allocate LAD object
        if ( allocated(this%LAD) ) deallocate(this%LAD)
        allocate(this%LAD)
        call this%LAD%init(this%decomp,this%der,this%gfil,2,this%dx,this%dy,this%dz,Cbeta,Cmu,Ckap,Cdiff,CY)

        ! Allocate mixture
        if ( allocated(this%mix) ) deallocate(this%mix)
        allocate(this%mix)
        call this%mix%init(this%decomp,this%der,this%fil,this%LAD,ns)
        !allocate(this%mix, source=solid_mixture(this%decomp,this%der,this%fil,this%LAD,ns))

        ! Allocate fields
        if ( allocated(this%fields) ) deallocate(this%fields) 
        call alloc_buffs(this%fields,nfields,'y',this%decomp)
        
        if ( allocated(this%Wcnsrv) ) deallocate(this%Wcnsrv) 
        call alloc_buffs(this%Wcnsrv,ncnsrv,'y',this%decomp)
        
        ! Associate pointers for ease of use
        this%rho  => this%fields(:,:,:, rho_index) 
        this%u    => this%fields(:,:,:,   u_index) 
        this%v    => this%fields(:,:,:,   v_index) 
        this%w    => this%fields(:,:,:,   w_index)  
        this%p    => this%fields(:,:,:,   p_index)  
        this%T    => this%fields(:,:,:,   T_index)  
        this%e    => this%fields(:,:,:,   e_index)  
        this%sos  => this%fields(:,:,:, sos_index)  
        this%mu   => this%fields(:,:,:,  mu_index)  
        this%bulk => this%fields(:,:,:,bulk_index)  
        this%kap  => this%fields(:,:,:, kap_index)   
       
        this%devstress => this%fields(:,:,:,sxx_index:szz_index)
        this%sxx  => this%fields(:,:,:, sxx_index)   
        this%sxy  => this%fields(:,:,:, sxy_index)   
        this%sxz  => this%fields(:,:,:, sxz_index)   
        this%syy  => this%fields(:,:,:, syy_index)   
        this%syz  => this%fields(:,:,:, syz_index)   
        this%szz  => this%fields(:,:,:, szz_index)   
        
        ! Initialize everything to a constant Zero
        this%fields = zero  

        ! Go to hooks if a different initialization is desired (Set mixture p, Ys, VF, u, v, w, rho)
        call initfields(this%decomp, this%dx, this%dy, this%dz, inputfile, this%mesh, this%fields, &
                        this%mix, this%tstop, this%dtfixed, tviz)
        ! Get hydrodynamic and elastic energies, stresses
        call this%mix%get_rhoYs_from_gVF(this%rho)  ! Get mixture rho and species Ys from species deformations and volume fractions
        call this%post_bc()
        
        !if (P_MAXVAL(abs( this%rho/this%rho0/(detG)**half - one )) > 10._rkind*eps) then
        !    call warning("Inconsistent initialization: rho/rho0 and g are not compatible")
        !end if

        ! Allocate 2 buffers for each of the three decompositions
        call alloc_buffs(this%xbuf,nbufsx,"x",this%decomp)
        call alloc_buffs(this%ybuf,nbufsy,"y",this%decomp)
        call alloc_buffs(this%zbuf,nbufsz,"z",this%decomp)

        this%SkewSymm = SkewSymm

        varnames( 1) = 'density'
        varnames( 2) = 'u'
        varnames( 3) = 'v'
        varnames( 4) = 'w'
        varnames( 5) = 'p'
        varnames( 6) = 'T'
        varnames( 7) = 'e'
        varnames( 8) = 'sos'
        varnames( 9) = 'mu'
        varnames(10) = 'bulk'
        varnames(11) = 'kap'
        varnames(12) = 'Sxx'
        varnames(13) = 'Sxy'
        varnames(14) = 'Sxz'
        varnames(15) = 'Syy'
        varnames(16) = 'Syz'
        varnames(17) = 'Szz'

        allocate(this%viz)
        call this%viz%init(this%outputdir, vizprefix, nfields, varnames)
        this%tviz = tviz

        ! Do this here to keep tau0 and invtau0 compatible at the end of this subroutine
        this%invtau0 = one/this%tau0

    end subroutine


    subroutine destroy(this)
        class(sgrid), intent(inout) :: this

        ! Nullify pointers
        nullify(this%x); nullify(this%y); nullify(this%z)

        nullify(this%rho      )
        nullify(this%u        )
        nullify(this%v        )
        nullify(this%w        )
        nullify(this%p        )
        nullify(this%T        )
        nullify(this%e        )
        nullify(this%sos      )
        nullify(this%mu       )
        nullify(this%bulk     )
        nullify(this%kap      )
        nullify(this%devstress)
        nullify(this%sxx      )
        nullify(this%sxy      )
        nullify(this%sxz      )
        nullify(this%syy      )
        nullify(this%syz      )
        nullify(this%szz      )
        if (allocated(this%mesh)) deallocate(this%mesh) 
        if (allocated(this%fields)) deallocate(this%fields) 
        
        call this%der%destroy()
        if (allocated(this%der)) deallocate(this%der) 
        
        call this%fil%destroy()
        if (allocated(this%fil)) deallocate(this%fil) 
        
        call this%gfil%destroy()
        if (allocated(this%gfil)) deallocate(this%gfil) 
        
        call destroy_buffs(this%xbuf)
        call destroy_buffs(this%ybuf)
        call destroy_buffs(this%zbuf)

        if ( allocated(this%mix) ) deallocate(this%mix)

        if ( allocated(this%LAD) ) deallocate(this%LAD)
        
        if (allocated(this%Wcnsrv)) deallocate(this%Wcnsrv) 
        
        call this%viz%destroy()
        if (allocated(this%viz)) deallocate(this%viz)

        call decomp_2d_finalize
        if (allocated(this%decomp)) deallocate(this%decomp) 

    end subroutine

    subroutine gradient(this, f, dfdx, dfdy, dfdz)
        class(sgrid),target, intent(inout) :: this
        real(rkind), intent(in), dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdx
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdy
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdz

        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        real(rkind), dimension(:,:,:), pointer :: xtmp,xdum,ztmp,zdum
        
        der => this%der
        decomp => this%decomp
        xtmp => this%xbuf(:,:,:,1)
        xdum => this%xbuf(:,:,:,2)
        ztmp => this%zbuf(:,:,:,1)
        zdum => this%zbuf(:,:,:,2)

        ! Get Y derivatives
        call der%ddy(f,dfdy)

        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,decomp)
        call der%ddx(xtmp,xdum)
        call transpose_x_to_y(xdum,dfdx)

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,decomp)
        call der%ddz(ztmp,zdum)
        call transpose_z_to_y(zdum,dfdz)

    end subroutine 

    subroutine laplacian(this, f, lapf)
        use timer
        class(sgrid),target, intent(inout) :: this
        real(rkind), intent(in), dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: lapf
        
        real(rkind), dimension(:,:,:), pointer :: xtmp,xdum,ztmp,zdum, ytmp
        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        

        der => this%der
        decomp => this%decomp
        xtmp => this%xbuf(:,:,:,1)
        xdum => this%xbuf(:,:,:,2)
        ztmp => this%zbuf(:,:,:,1)
        zdum => this%zbuf(:,:,:,2)
        ytmp => this%ybuf(:,:,:,1)

        ! Get Y derivatives
        call der%d2dy2(f,lapf)
        
        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,this%decomp) 
        call this%der%d2dx2(xtmp,xdum)
        call transpose_x_to_y(xdum,ytmp,this%decomp)

        lapf = lapf + ytmp

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,this%decomp)
        call this%der%d2dz2(ztmp,zdum)
        call transpose_z_to_y(zdum,ytmp,this%decomp)
        
        lapf = lapf + ytmp

    end subroutine

    subroutine simulate(this)
        use reductions, only: P_MEAN
        use timer,      only: tic, toc
        use exits,      only: GracefulExit, message
        class(sgrid), target, intent(inout) :: this

        logical :: tcond, vizcond, stepcond
        character(len=clen) :: stability
        real(rkind) :: cputime
        real(rkind), dimension(:,:,:,:), allocatable, target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: ehmix
        integer :: i, imat

        allocate( duidxj(this%nxp, this%nyp, this%nzp, 9) )
        ! Get artificial properties for initial conditions
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
            !write(*,'(a,2(e21.14,1x))') 'inf', maxval(this%p)-0.1d0,  minval(this%p)-0.1d0
            !write(*,'(a,2(e21.14,1x))') 'enf', maxval(this%e)-3.75d0,  minval(this%e)-3.75d0
        call this%gradient(this%u,dudx,dudy,dudz)
        call this%gradient(this%v,dvdx,dvdy,dvdz)
        call this%gradient(this%w,dwdx,dwdy,dwdz)
            !write(*,'(a,2(e21.14,1x))') 'inf', maxval(this%p)-0.1d0,  minval(this%p)-0.1d0
            !write(*,'(a,2(e21.14,1x))') 'enf', maxval(this%e)-3.75d0,  minval(this%e)-3.75d0

        do i=1,this%mix%ns
            call this%gradient(this%mix%material(i)%Ys,this%mix%material(i)%Ji(:,:,:,1),&
                               this%mix%material(i)%Ji(:,:,:,2),this%mix%material(i)%Ji(:,:,:,3))
        end do

        ! compute artificial shear and bulk viscosities
        call this%getPhysicalProperties()
            !write(*,'(a,2(e21.14,1x))') 'inf', maxval(this%p)-0.1d0,  minval(this%p)-0.1d0
            !write(*,'(a,2(e21.14,1x))') 'enf', maxval(this%e)-3.75d0,  minval(this%e)-3.75d0
        call this%LAD%get_viscosities(this%rho,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc)
            !write(*,'(a,2(e21.14,1x))') 'inf', maxval(this%p)-0.1d0,  minval(this%p)-0.1d0
            !write(*,'(a,2(e21.14,1x))') 'enf', maxval(this%e)-3.75d0,  minval(this%e)-3.75d0
        ehmix => duidxj(:,:,:,4) ! use some storage space
        ehmix = this%e
        do imat = 1, this%mix%ns
            ehmix = ehmix - this%mix%material(imat)%Ys * this%mix%material(imat)%eel
        enddo
        call this%LAD%get_conductivity(this%rho,ehmix,this%T,this%sos,this%kap,this%x_bc,this%y_bc,this%z_bc) ! Change to eh instead of e
            !write(*,'(a,2(e21.14,1x))') 'inf', maxval(this%p)-0.1d0,  minval(this%p)-0.1d0
            !write(*,'(a,2(e21.14,1x))') 'enf', maxval(this%e)-3.75d0,  minval(this%e)-3.75d0

        nullify(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,ehmix)
        deallocate( duidxj )

        ! compute species artificial conductivities and diffusivities
        ! call this%mix%getSOS(this%rho,this%p,this%sos)
        call this%mix%getLAD(this%rho,this%sos,this%x_bc,this%y_bc,this%z_bc)  ! Compute species LAD (kap, diff)
        ! ------------------------------------------------

        call this%get_dt(stability)

        ! Write out initial conditions
        call hook_output(this%decomp, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
        call this%viz%WriteViz(this%decomp, this%mesh, this%fields, this%tsim)
        vizcond = .FALSE.
        
        ! Check for visualization condition and adjust time step
        if ( (this%tviz > zero) .AND. (this%tsim + this%dt > this%tviz * this%viz%vizcount) ) then
            this%dt = this%tviz * this%viz%vizcount - this%tsim
            vizcond = .TRUE.
            stability = 'vizdump'
        end if

        tcond = .TRUE.
        ! Check tstop condition
        if ( (this%tstop > zero) .AND. (this%tsim >= this%tstop) ) then
            tcond = .FALSE.
        else if ( (this%tstop > zero) .AND. (this%tsim + this%dt >= this%tstop) ) then
            this%dt = this%tstop - this%tsim
        end if

        ! Check nsteps condition
        if ( (this%nsteps <= 0) .OR. (this%step < this%nsteps) ) then
            stepcond = .TRUE.
        else
            stepcond = .FALSE.
        end if

        if ( (this%tstop <= zero) .AND. (this%nsteps <= 0) ) then
            call GracefulExit('No stopping criterion set. Set either tstop or nsteps to be positive.', 345)
        end if

        ! Start the simulation while loop
        do while ( tcond .AND. stepcond )
            ! Advance time
            call tic()
            call this%advance_RK45()
            call toc(cputime)
            
            call message(1,"Time",this%tsim)
            call message(2,"Time step",this%dt)
            call message(2,"Stability limit: "//trim(stability))
            call message(2,"CPU time (in seconds)",cputime)
            call hook_timestep(this%decomp, this%mesh, this%fields, this%mix, this%step, this%tsim)
          
            ! Write out vizualization dump if vizcond is met 
            if (vizcond) then
                call hook_output(this%decomp, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
                call this%viz%WriteViz(this%decomp, this%mesh, this%fields, this%tsim)
                vizcond = .FALSE.
            end if
            
            ! Get the new time step
            call this%get_dt(stability)
            
            ! Check for visualization condition and adjust time step
            if ( (this%tviz > zero) .AND. (this%tsim + this%dt >= this%tviz * this%viz%vizcount) ) then
                this%dt = this%tviz * this%viz%vizcount - this%tsim
                vizcond = .TRUE.
            end if

            ! Check tstop condition
            if ( (this%tstop > zero) .AND. (this%tsim >= this%tstop - eps) ) then
                tcond = .FALSE.
            else if ( (this%tstop > zero) .AND. (this%tsim + this%dt + eps >= this%tstop) ) then
                this%dt = this%tstop - this%tsim
                stability = 'stop'
                vizcond = .TRUE.
            end if

            ! Check nsteps condition
            if ( (this%nsteps <= 0) .OR. (this%step < this%nsteps) ) then
                stepcond = .TRUE.
            else
                stepcond = .FALSE.
            end if

        end do

    end subroutine

    subroutine advance_RK45(this)
        use RKCoeffs,   only: RK45_steps,RK45_A,RK45_B
        use exits,      only: message,nancheck,GracefulExit
        use reductions, only: P_MAXVAL, P_MINVAL
        use mpi
        class(sgrid), target, intent(inout) :: this

        real(rkind)                                               :: Qtmpt     ! Temporary variable for RK45
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv) :: rhs       ! RHS for conserved variables
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv) :: Qtmp      ! Temporary variable for RK45
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: divu      ! Velocity divergence for species energy eq
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: viscwork  ! Viscous work term for species energy eq
        integer :: isub,i,j,k,l,ierr

        character(len=clen) :: charout

        Qtmp  = zero
        Qtmpt = zero

        do isub = 1,RK45_steps
            !write(*,*) '------RK45 substep', isub, '--------'
            ! hard-coded velocity for material interface straining problem
            !this%u = -sin(pi*this%x)**2*sin(two*pi*this%y)*cos(pi*this%tsim/4.0_rkind) ! u
            !this%v =  sin(pi*this%y)**2*sin(two*pi*this%x)*cos(pi*this%tsim/4.0_rkind)  ! v

            call this%get_conserved()

            if ( nancheck(this%Wcnsrv,i,j,k,l) ) then
                call message("Wcnsrv: ",this%Wcnsrv(i,j,k,l))
                write(charout,'(A,I1,A,I5,A,4(I5,A))') "NaN encountered in solution (Wcnsrv) at &
                    &substep ", isub, " of step ", this%step+1, " at (",i,", ",j,", ",k,", ",l,") of Wcnsrv"
                call GracefulExit(trim(charout), 999)
            end if
            call this%mix%checkNaN()

            ! Pre-compute stress, LAD, J, etc.
            ! call this%mix%getSOS(this%rho,this%p,this%sos)
            call this%mix%getLAD(this%rho,this%sos,this%x_bc,this%y_bc,this%z_bc)  ! Compute species LAD (kap, diff)
            call this%mix%get_J(this%rho)                                          ! Compute diffusive mass fluxes
            call this%mix%get_q(this%x_bc,this%y_bc,this%z_bc)                     ! Compute diffusive thermal fluxes (including enthalpy diffusion)

            ! Update total mixture conserved variables
            call this%getRHS(rhs,divu,viscwork)
              !rhs = zero
            Qtmp  = this%dt*rhs  + RK45_A(isub)*Qtmp
              !write(*,'(a,3(e19.12,1x))') 'wcne', this%Wcnsrv(1,1,1,4), this%Wcnsrv(1,1,1,4)
              !write(*,'(a,3(e19.12,1x))') 'Ycne', this%rho(51,1,1), this%mix%material(1)%Ys(51,1,1), this%mix%material(2)%Ys(51,1,1)
              !write(*,'(a,3(e19.12,1x))') 'Ycne', this%mix%material(2)%consrv(51,1,1,1), this%mix%material(1)%consrv(51,1,1,1)!, this%mix%material(2)%consrv(51,1,1,1) + this%mix%material(1)%consrv(51,1,1,1)
              !write(*,'(a,3(e19.12,1x))') 'Ycne', this%mix%material(2)%consrv(51,1,1,1) + this%mix%material(1)%consrv(51,1,1,1)!, this%mix%material(2)%consrv(51,1,1,1) + this%mix%material(1)%consrv(51,1,1,1)
              !write(*,'(a,3(e19.12,1x))')
            this%Wcnsrv = this%Wcnsrv + RK45_B(isub)*Qtmp
            !write(*,'(a,2(e21.14,1x))') '-rhsu-', maxval(rhs(:,:,:,1)),  minval(rhs(:,:,:,1))
            !write(*,'(a,2(e21.14,1x))') '-rhsv-', maxval(rhs(:,:,:,2)),  minval(rhs(:,:,:,2))
            !write(*,'(a,2(e21.14,1x))') '-rhsw-', maxval(rhs(:,:,:,3)),  minval(rhs(:,:,:,3))
            !write(*,'(a,2(e21.14,1x))') '-rhse-', maxval(rhs(:,:,:,4)),  minval(rhs(:,:,:,4))

            ! Now update all the individual species variables
            call this%mix%update_g (isub,this%dt,this%rho,this%u,this%v,this%w,this%x,this%y,this%z,this%tsim)               ! g tensor
            call this%mix%update_Ys(isub,this%dt,this%rho,this%u,this%v,this%w,this%x,this%y,this%z,this%tsim)               ! Volume Fraction
            !call this%mix%update_eh(isub,this%dt,this%rho,this%u,this%v,this%w,this%x,this%y,this%z,this%tsim,divu,viscwork) ! Hydrodynamic energy
            !call this%mix%update_VF(isub,this%dt,this%u,this%v,this%w,this%x,this%y,this%z,this%tsim)                        ! Volume Fraction

            ! Integrate simulation time to keep it in sync with RK substep
            Qtmpt = this%dt + RK45_A(isub)*Qtmpt
            this%tsim = this%tsim + RK45_B(isub)*Qtmpt
              !write(*,'(a,3(e19.12,1x))') 'wcne', this%Wcnsrv(40,1,1,4), this%Wcnsrv(51,1,1,4)
              !write(*,'(a,3(e19.12,1x))') 'Ycn1', this%mix%material(1)%consrv(1,1,1,1), this%mix%material(1)%consrv(1,1,1,1)
              !write(*,'(a,3(e19.12,1x))') 'Ycn2', this%mix%material(2)%consrv(40,1,1,1), this%mix%material(2)%consrv(51,1,1,1)
              !write(*,'(a,3(e19.12,1x))') 'Ycns', this%mix%material(2)%consrv(40,1,1,1) + this%mix%material(1)%consrv(40,1,1,1), this%mix%material(2)%consrv(51,1,1,1) + this%mix%material(1)%consrv(51,1,1,1)
              !write(*,'(a,3(e19.12,1x))')

            ! Filter the conserved variables
            do i = 1,ncnsrv
                call this%filter(this%Wcnsrv(:,:,:,i), this%fil, 1)
            end do
             call mpi_barrier(mpi_comm_world,ierr) 
              !write(*,'(a,3(e19.12,1x))') 'wcnf', this%Wcnsrv(1,1,1,4), this%Wcnsrv(1,1,1,4)
             call mpi_barrier(mpi_comm_world,ierr) 
           

            ! Filter the individual species variables
            call this%mix%filter(1)
              !write(*,'(a,3(e19.12,1x))') 'Ys1', this%mix%material(1)%consrv(40,1,1,1), this%mix%material(1)%consrv(51,1,1,1)
              !write(*,'(a,3(e19.12,1x))') 'Ys2', this%mix%material(2)%consrv(1,1,1,1), this%mix%material(2)%consrv(1,1,1,1)
            
            call this%get_primitive()

            ! if (.NOT. this%explPlast) then
            !     if (this%plastic) then
            !         ! Effect plastic deformations
            !         ! call this%elastic%plastic_deformation(this%g)
            !         call this%mix%plastic_deformation()
            !         call this%get_primitive()

            !         ! Filter the conserved variables
            !         do i = 1,5
            !             call this%filter(this%Wcnsrv(:,:,:,i), this%fil, 1)
            !         end do
            !         ! Filter the g tensor
            !         do i = 1,9
            !             call this%filter(this%g(:,:,:,i), this%fil, 1)
            !         end do
            !     end if
            ! end if

             !call this%mix%relaxPressure(this%rho, this%e, this%p)
            !write(*,'(a,2(e21.14,1x))') 'pbe', maxval(this%p)-0.1d0,  minval(this%p)-0.1d0
            !write(*,'(a,2(e21.14,1x))') 'ebe', maxval(this%e)-3.75d0, minval(this%e)-3.75d0
            !write(*,'(a,3(e19.12,1x))') '-e-', this%e(40,1,1), this%e(51,1,1)

             call this%mix%equilibratePressureTemperature(this%rho, this%e, this%p, this%T)
            
            call hook_bc(this%decomp, this%mesh, this%fields, this%mix, this%tsim)
            !write(*,'(a,2(e21.14,1x))') '-p-', maxval(this%p)-0.1d0,  minval(this%p)-0.1d0
            !write(*,'(a,2(e21.14,1x))') '-e-', maxval(this%e)-3.75d0, minval(this%e)-3.75d0
            call this%post_bc()
        end do

        ! this%tsim = this%tsim + this%dt
        this%step = this%step + 1
            
    end subroutine

    subroutine get_dt(this,stability)
        use reductions, only : P_MAXVAL
        class(sgrid), target, intent(inout) :: this
        character(len=*), intent(out) :: stability
        real(rkind) :: dtCFL, dtmu, dtbulk, dtkap, dtdiff, dtplast, delta

        delta = min(this%dx, this%dy, this%dz)

        ! call this%sgas%get_sos(this%rho,this%p,cs)  ! Speed of sound - hydrodynamic part
        ! call this%elastic%get_sos(this%rho0,cs)     ! Speed of sound - elastic part

        ! continuum
        dtCFL  = this%CFL / P_MAXVAL( ABS(this%u)/this%dx + ABS(this%v)/this%dy + ABS(this%w)/this%dz &
               + this%sos*sqrt( one/(this%dx**2) + one/(this%dy**2) + one/(this%dz**2) ))
        dtmu   = 0.2_rkind * delta**2 / (P_MAXVAL( this%mu  / this%rho ) + eps)
        dtbulk = 0.2_rkind * delta**2 / (P_MAXVAL( this%bulk/ this%rho ) + eps)
        !dtkap  = 0.2_rkind * delta**2 / (P_MAXVAL( this%kap*this%T/(this%rho*this%sos**2)) + eps)     ! Cook (2007) formulation
        dtkap  = 0.2_rkind * one / ( (P_MAXVAL(this%kap*this%T/(this%rho*delta**4)))**(third) + eps)    ! Cook (2009) formulation

        ! species specific
        call this%mix%get_dt(this%rho, delta, dtkap, dtdiff, dtplast)
        !dtkap = 0.2_rkind * dtkap
        dtdiff = 0.2_rkind * dtdiff

        ! Use fixed time step if CFL <= 0
        if ( this%CFL .LE. zero ) then
            this%dt = this%dtfixed
            stability = 'fixed'
        else
            stability = 'convective'
            this%dt = dtCFL
            if ( this%dt > dtmu ) then
                this%dt = dtmu
                stability = 'shear'
            else if ( this%dt > dtbulk ) then
                this%dt = dtbulk
                stability = 'bulk'
            else if ( this%dt > dtkap ) then
                this%dt = dtkap
                stability = 'conductive'
            else if ( this%dt > dtdiff ) then
                this%dt = dtdiff
                stability = 'diffusive'
            else if ( this%dt > dtplast ) then
                this%dt = dtplast
                stability = 'plastic'
            end if

            if (this%step .LE. 10) then
                this%dt = this%dt / 10._rkind
                stability = 'startup'
            end if
        end if

    end subroutine

    subroutine get_primitive(this)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(:,:,:), pointer :: onebyrho
        real(rkind), dimension(:,:,:), pointer :: rhou,rhov,rhow,TE

        onebyrho => this%ybuf(:,:,:,1)

        ! this%rho  =  this%Wcnsrv(:,:,:,1)
        
        call this%mix%get_rho(this%rho)

        rhou => this%Wcnsrv(:,:,:,mom_index  )
        rhov => this%Wcnsrv(:,:,:,mom_index+1)
        rhow => this%Wcnsrv(:,:,:,mom_index+2)
        TE   => this%Wcnsrv(:,:,:, TE_index  )

        onebyrho = one/this%rho
        this%u = rhou * onebyrho
        this%v = rhov * onebyrho
        this%w = rhow * onebyrho
        this%e = (TE*onebyrho) - half*( this%u*this%u + this%v*this%v + this%w*this%w )
       
        call this%mix%get_primitive(this%rho, this%devstress, this%p, this%sos, this%e)                  ! Get primitive variables for individual species

    end subroutine

    !pure subroutine get_conserved(this)
    subroutine get_conserved(this)
        class(sgrid), intent(inout) :: this

        ! Assume rho is already available
        this%Wcnsrv(:,:,:,mom_index  ) = this%rho * this%u
        this%Wcnsrv(:,:,:,mom_index+1) = this%rho * this%v
        this%Wcnsrv(:,:,:,mom_index+2) = this%rho * this%w
            !write(*,'(a,2(e21.14,1x))') 'ecn', maxval(this%e)-3.75d0,  minval(this%e)-3.75d0
        this%Wcnsrv(:,:,:, TE_index  ) = this%rho * ( this%e + half*( this%u*this%u + this%v*this%v + this%w*this%w ) )
            !this%e = this%Wcnsrv(:,:,:, TE_index  ) / this%rho - ( half*( this%u*this%u + this%v*this%v + this%w*this%w ) )
            !write(*,'(a,2(e21.14,1x))') 'ecn', maxval(this%e)-3.75d0,  minval(this%e)-3.75d0

        ! add 2M (mass fraction and hydrodynamic energy) variables here
        call this%mix%get_conserved(this%rho)

    end subroutine

    subroutine post_bc(this)
        class(sgrid), intent(inout) :: this

        call this%mix%get_eelastic_devstress(this%devstress)   ! Get species elastic energies, and mixture and species devstress
        call this%mix%get_ehydro_from_p(this%rho,this%p,this%T)            ! Get species hydrodynamic energy, temperature; and mixture pressure, temperature
        !write(*,'(a,4(e19.12,1x))') 'rho,p,e,t  ', this%rho(71,1,1), this%p(71,1,1), this%e(71,1,1), this%T(71,1,1)
        !write(*,'(a,4(e19.12,1x))') 'rho,p,e,t-1', this%mix%material(1)%Ys(71,1,1)*this%rho(71,1,1)/this%mix%material(1)%VF(71,1,1), this%mix%material(1)%p(71,1,1), this%mix%material(1)%eh(71,1,1), this%mix%material(1)%T(71,1,1)
        !write(*,'(a,4(e19.12,1x))') 'rho,p,e,t-2', this%mix%material(2)%Ys(71,1,1)*this%rho(71,1,1)/this%mix%material(2)%VF(71,1,1), this%mix%material(2)%p(71,1,1), this%mix%material(2)%eh(71,1,1), this%mix%material(2)%T(71,1,1)
        !write(*,'(a,4(e19.12,1x))') 'Cv1', (this%mix%material(1)%hydro%Cv), (this%mix%material(1)%hydro%Cv), (this%mix%material(2)%hydro%Cv), (this%mix%material(2)%hydro%Cv)
        !write(*,'(a,4(e19.12,1x))') 'Pin', (this%mix%material(1)%hydro%PInf), (this%mix%material(1)%hydro%Pinf), (this%mix%material(2)%hydro%Pinf), (this%mix%material(2)%hydro%PInf)
        !write(*,'(a,4(e19.12,1x))') 'e-1', maxval(this%mix%material(1)%eh), minval(this%mix%material(1)%eh), maxval(this%mix%material(2)%eh), minval(this%mix%material(2)%eh)
        call this%mix%getSOS(this%rho,this%p,this%sos)

        ! assuming pressures have relaxed and sum( (Ys*(ehydro + eelastic) ) over all
        ! materials equals e
        call this%mix%get_emix(this%e)
       
            !write(*,'(a,2(e21.14,1x))') 'pbc', maxval(this%p)-0.1d0,  minval(this%p)-0.1d0
            !write(*,'(a,2(e21.14,1x))') 'ebc', maxval(this%e)-3.75d0,  minval(this%e)-3.75d0
        !call this%mix%get_T(this%T) ! Get mixture temperature (?)
        !call this%sgas%get_T(this%e,this%T)  ! Get updated temperature

    end subroutine

    subroutine getRHS(this, rhs, divu, viscwork)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,ncnsrv), intent(out) :: rhs
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),     intent(out) :: divu
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),     intent(out) :: viscwork
        !real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), intent(out) :: rhsg
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        real(rkind), dimension(:,:,:), pointer :: qx,qy,qz
        real(rkind), dimension(:,:,:), pointer :: ehmix
        integer :: imat

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        call this%gradient(this%u,dudx,dudy,dudz)
        call this%gradient(this%v,dvdx,dvdy,dvdz)
        call this%gradient(this%w,dwdx,dwdy,dwdz)

        divu = dudx + dvdy + dwdz

        call this%getPhysicalProperties()
        call this%LAD%get_viscosities(this%rho,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc)

        ! subtract elastic energies to determine mixture hydrostatic energy. conductivity 
        ! is assumed a function of only hydrostatic energy
        ehmix => duidxj(:,:,:,tauxxidx) ! use some storage space
        ehmix = this%e
        do imat = 1, this%mix%ns
            ehmix = ehmix - this%mix%material(imat)%Ys * this%mix%material(imat)%eel
        enddo
        call this%LAD%get_conductivity(this%rho,ehmix,this%T,this%sos,this%kap,this%x_bc,this%y_bc,this%z_bc) ! Change to eh instead of e

        ! Get tau tensor and q (heat conduction) vector. Put in components of duidxj
        call this%get_tau( duidxj, viscwork )
        ! Now, associate the pointers to understand what's going on better
        tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx); tauxz => duidxj(:,:,:,tauxzidx);
                                         tauyy => duidxj(:,:,:,tauyyidx); tauyz => duidxj(:,:,:,tauyzidx);
                                                                          tauzz => duidxj(:,:,:,tauzzidx);

        ! for use in species hydrodynamic equations --- check carefully that
        ! dudx, dvdy, dwdz have not been destroyed in get_tau
        ! viscwork = tauxx*dudx + tauyy*dvdy + tauzz*dwdz
       
        ! Add the deviatoric stress to the tau for use in fluxes 
        tauxx = tauxx + this%sxx; tauxy = tauxy + this%sxy; tauxz = tauxz + this%sxz
                                  tauyy = tauyy + this%syy; tauyz = tauyz + this%syz
                                                            tauzz = tauzz + this%szz
        
        qx => duidxj(:,:,:,qxidx); qy => duidxj(:,:,:,qyidx); qz => duidxj(:,:,:,qzidx);
        call this%mix%get_qmix(qx, qy, qz)     ! get only inter-species enthalpy fluxes
            !write(*,'(a,2(e21.14,1x))') 'qxb', maxval(qx),  minval(qx)
            !write(*,'(a,2(e21.14,1x))') 'qyb', maxval(qy),  minval(qy)
            !write(*,'(a,2(e21.14,1x))') 'qzb', maxval(qz),  minval(qz)
            !write(*,'(a,2(e21.14,1x))') '---'
        call this%get_q(qx, qy, qz)            ! add artificial thermal conduction fluxes
            !write(*,'(a,2(e21.14,1x))') 'qxa', maxval(qx),  minval(qx)
            !write(*,'(a,2(e21.14,1x))') 'qya', maxval(qy),  minval(qy)
            !write(*,'(a,2(e21.14,1x))') 'qza', maxval(qz),  minval(qz)
            !write(*,'(a,2(e21.14,1x))') '---'

        rhs = zero
        call this%getRHS_x(              rhs,&
                           tauxx,tauxy,tauxz,&
                               qx )

        call this%getRHS_y(              rhs,&
                           tauxy,tauyy,tauyz,&
                               qy )

        call this%getRHS_z(              rhs,&
                           tauxz,tauyz,tauzz,&
                               qz )

        ! Call problem source hook
        call hook_mixture_source(this%decomp, this%mesh, this%fields, this%mix, this%tsim, rhs)
 
    end subroutine

    subroutine getRHS_x(       this,  rhs,&
                        tauxx,tauxy,tauxz,&
                            qx )
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxx,tauxy,tauxz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qx

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: flux
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2

        ! flux => this%ybuf(:,:,:,1:4)
        xtmp1 => this%xbuf(:,:,:,1); xtmp2 => this%xbuf(:,:,:,2)

        ! flux(:,:,:,1) = this%Wcnsrv(:,:,:,2)   ! rho*u
        ! flux(:,:,:,1) = this%Wcnsrv(:,:,:,1)*this%u + this%p - tauxx
        ! flux(:,:,:,2) = this%Wcnsrv(:,:,:,1)*this%v          - tauxy
        ! flux(:,:,:,3) = this%Wcnsrv(:,:,:,1)*this%w          - tauxz
        ! flux(:,:,:,4) = (this%Wcnsrv(:,:,:,4) + this%p - tauxx)*this%u - this%v*tauxy - this%w*tauxz - qx

        ! ! Now, get the x-derivative of the fluxes
        ! do i=1,4
        !     call transpose_y_to_x(flux(:,:,:,i),xtmp1,this%decomp)
        !     call this%der%ddx(xtmp1,xtmp2)
        !     call transpose_x_to_y(xtmp2,flux(:,:,:,i),this%decomp)
        ! end do

        ! ! Add to rhs
        ! rhs = rhs - flux

        flux = this%Wcnsrv(:,:,:,mom_index  )*this%u + this%p - tauxx ! x-momentum
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2, this%x_bc(1), this%x_bc(2)) ! Symmetric for x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux

        flux = this%Wcnsrv(:,:,:,mom_index  )*this%v          - tauxy ! y-momentum
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux

        flux = this%Wcnsrv(:,:,:,mom_index  )*this%w          - tauxz ! z-momentum
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux

        flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauxx)*this%u - this%v*tauxy - this%w*tauxz + qx ! Total Energy
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux

    end subroutine

    subroutine getRHS_y(       this,  rhs,&
                        tauxy,tauyy,tauyz,&
                            qy )
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxy,tauyy,tauyz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qy

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: flux
        real(rkind), dimension(:,:,:), pointer :: ytmp1

        ! flux => this%ybuf(:,:,:,1:5)
        ytmp1 => this%ybuf(:,:,:,6)

        ! flux(:,:,:,1) = this%Wcnsrv(:,:,:,3)   ! rho*v
        ! flux(:,:,:,2) = this%Wcnsrv(:,:,:,3)*this%u          - tauxy
        ! flux(:,:,:,3) = this%Wcnsrv(:,:,:,3)*this%v + this%p - tauyy
        ! flux(:,:,:,4) = this%Wcnsrv(:,:,:,3)*this%w          - tauyz
        ! flux(:,:,:,5) = (this%Wcnsrv(:,:,:,5) + this%p - tauyy)*this%v - this%u*tauxy - this%w*tauyz - qy

        ! ! Now, get the x-derivative of the fluxes
        ! do i=1,5
        !     call this%der%ddy(flux(:,:,:,i),ytmp1)

        !     ! Add to rhs
        !     rhs(:,:,:,i) = rhs(:,:,:,i) - ytmp1
        ! end do

        flux = this%Wcnsrv(:,:,:,mom_index+1)*this%u          - tauxy ! x-momentum
        call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - ytmp1

        flux = this%Wcnsrv(:,:,:,mom_index+1)*this%v + this%p - tauyy ! y-momentum
        call this%der%ddy(flux,ytmp1, this%y_bc(1), this%y_bc(2)) ! Symmetric for y-momentum
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - ytmp1

        flux = this%Wcnsrv(:,:,:,mom_index+1)*this%w          - tauyz ! z-momentum
        call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - ytmp1

        flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauyy)*this%v - this%u*tauxy - this%w*tauyz + qy ! Total Energy
        call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - ytmp1


    end subroutine

    subroutine getRHS_z(       this,  rhs,&
                        tauxz,tauyz,tauzz,&
                            qz )
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxz,tauyz,tauzz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qz

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: flux
        real(rkind), dimension(:,:,:), pointer :: ztmp1,ztmp2

        ! flux => this%ybuf(:,:,:,1:5)
        ztmp1 => this%zbuf(:,:,:,1); ztmp2 => this%zbuf(:,:,:,2)

        ! flux(:,:,:,1) = this%Wcnsrv(:,:,:,4)   ! rho*w
        ! flux(:,:,:,2) = this%Wcnsrv(:,:,:,4)*this%u          - tauxz
        ! flux(:,:,:,3) = this%Wcnsrv(:,:,:,4)*this%v          - tauyz
        ! flux(:,:,:,4) = this%Wcnsrv(:,:,:,4)*this%w + this%p - tauzz
        ! flux(:,:,:,5) = (this%Wcnsrv(:,:,:,5) + this%p - tauzz)*this%w - this%u*tauxz - this%v*tauyz - qz

        ! ! Now, get the x-derivative of the fluxes
        ! do i=1,5
        !     call transpose_y_to_z(flux(:,:,:,i),ztmp1,this%decomp)
        !     call this%der%ddz(ztmp1,ztmp2)
        !     call transpose_z_to_y(ztmp2,flux(:,:,:,i),this%decomp)
        ! end do

        ! ! Add to rhs
        ! rhs = rhs - flux

        flux = this%Wcnsrv(:,:,:,mom_index+2)*this%u          - tauxz ! x-momentum
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux

        flux = this%Wcnsrv(:,:,:,mom_index+2)*this%v          - tauyz ! y-momentum
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux

        flux = this%Wcnsrv(:,:,:,mom_index+2)*this%w + this%p - tauzz ! z-momentum
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2, this%z_bc(1), this%z_bc(2)) ! Symmetric for z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux

        flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauzz)*this%w - this%u*tauxz - this%v*tauyz + qz ! Total Energy
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux

    end subroutine

    subroutine filter(this,arr,myfil,numtimes)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(inout) :: arr
        type(filters), target, optional, intent(in) :: myfil
        integer, optional, intent(in) :: numtimes
        
        type(filters), pointer :: fil2use
        integer :: times2fil
        real(rkind), dimension(:,:,:), pointer :: tmp_in_y, tmp1_in_x, tmp1_in_z, tmp2_in_x, tmp2_in_z
        integer :: lastx, lasty, lastz, idx


        if (present(myfil)) then
            fil2use => myfil
        else
            fil2use => this%fil
        end if 

        if (present(numtimes)) then
            times2fil = numtimes
        else
            times2fil = 1
        end if

        ! Allocate pointers for the needed buffers 
        ! Atleast 2 buffers in x and z are assumed
        ! Last two buffers are occupied

        lastx = size(this%xbuf,4)
        lasty = size(this%ybuf,4)
        lastz = size(this%zbuf,4)

        tmp1_in_x => this%xbuf(:,:,:,lastx)
        tmp2_in_x => this%xbuf(:,:,:,lastx-1)
        tmp_in_y => this%ybuf(:,:,:,lasty)
        tmp1_in_z => this%zbuf(:,:,:,lastz)
        tmp2_in_z => this%zbuf(:,:,:,lastz-1)
       
        
        ! First filter in y
        call fil2use%filtery(arr,tmp_in_y)
        ! Subsequent refilters 
        do idx = 1,times2fil-1
            arr = tmp_in_y
            call fil2use%filtery(arr,tmp_in_y)
        end do
        
        ! Then transpose to x
        call transpose_y_to_x(tmp_in_y,tmp1_in_x,this%decomp)

        ! First filter in x
        call fil2use%filterx(tmp1_in_x,tmp2_in_x)
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_x = tmp2_in_x
            call fil2use%filterx(tmp1_in_x,tmp2_in_x)
        end do 

        ! Now transpose back to y
        call transpose_x_to_y(tmp2_in_x,tmp_in_y,this%decomp)

        ! Now transpose to z
        call transpose_y_to_z(tmp_in_y,tmp1_in_z,this%decomp)

        !First filter in z
        call fil2use%filterz(tmp1_in_z,tmp2_in_z)
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_z = tmp2_in_z
            call fil2use%filterz(tmp1_in_z,tmp2_in_z)
        end do 

        ! Now transpose back to y
        call transpose_z_to_y(tmp2_in_z,arr,this%decomp)

        ! Finished
    end subroutine
   
    subroutine getPhysicalProperties(this)
        class(sgrid), intent(inout) :: this

        ! If inviscid set everything to zero (otherwise use a model)
        this%mu = zero
        this%bulk = zero
        this%kap = zero

    end subroutine  

    ! Get tau_ij for momentum equation and simulataneously calculate the viscous work
    ! for the material hydrodynamic equations
    subroutine get_tau(this,duidxj,viscwork)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), target, intent(inout) :: duidxj
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),           intent(out)   :: viscwork

        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: lambda, bambda

        lambda => this%ybuf(:,:,:,1)
        bambda => this%ybuf(:,:,:,2)

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
       
        ! Compute the multiplying factors (thermo-shit)
        bambda = (four/three)*this%mu + this%bulk
        lambda = this%bulk - (two/three)*this%mu

        ! Step 1: Get tau_12  (dudy is destroyed)
        dudy =  dudy + dvdx
        viscwork  = (this%mu*dudy) * (dudy)  ! tau_12 * (S_12+S_21) (Since symmetric)
        dudy = this%mu*dudy
        !tauxyidz = 2
    
        ! Step 2: Get tau_13 (dudz is destroyed)
        dudz = dudz + dwdx
        viscwork  = viscwork + (this%mu*dudz) * (dudz)  ! tau_13 *(S_13 + S_31)
        dudz = this%mu*dudz
        !tauxzidx = 3

        ! Step 3: Get tau_23 (dvdz is destroyed)
        dvdz = dvdz + dwdy
        viscwork  = viscwork + (this%mu*dvdz) * (half*dvdz)  ! tau_23 * (S_23 + S_32)
        dvdz = this%mu*dvdz
        !tauyzidx = 6

        ! Step 4: Get tau_11 (dvdx is destroyed)
        dvdx = bambda*dudx + lambda*(dvdy + dwdz)
        viscwork  = viscwork + dvdx * dudx  ! tau_11 * S_11
        !tauxxidx = 4

        ! Step 5: Get tau_22 (dwdx is destroyed)
        dwdx = bambda*dvdy + lambda*(dudx + dwdz)
        viscwork  = viscwork + dwdx * dvdy  ! tau_22 * S_22
        !tauyyidx = 7

        ! Step 6: Get tau_33 (dwdy is destroyed)
        dwdy = bambda*dwdz + lambda*(dudx + dvdy)
        viscwork  = viscwork + dwdy * dwdz  ! tau_33 * S_33
        !tauzzidx = 8

        ! Done 
    end subroutine 

    subroutine get_q(this,qx,qy,qz)
        use exits, only: nancheck
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(inout) :: qx,qy,qz

        integer :: i
        real(rkind), dimension(:,:,:), pointer :: tmp1_in_x, tmp2_in_x, tmp1_in_y, tmp1_in_z, tmp2_in_z
        type(derivatives), pointer :: der

        der => this%der

        tmp1_in_x => this%xbuf(:,:,:,1)
        tmp2_in_x => this%xbuf(:,:,:,2)

        tmp1_in_z => this%zbuf(:,:,:,1)
        tmp2_in_z => this%zbuf(:,:,:,2)

        tmp1_in_y => this%ybuf(:,:,:,1)

        ! Step 1: Get qy (dvdy is destroyed)
        call der%ddy(this%T,tmp1_in_y,this%y_bc(1),this%y_bc(2))
        qy = qy - this%kap*tmp1_in_y

        ! Step 2: Get qx (dudx is destroyed)
        call transpose_y_to_x(this%T,tmp1_in_x,this%decomp)
        call der%ddx(tmp1_in_x,tmp2_in_x,this%x_bc(1),this%x_bc(2))
        call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
        qx = qx - this%kap*tmp1_in_y

        ! Step 3: Get qz (dwdz is destroyed)
        call transpose_y_to_z(this%T,tmp1_in_z,this%decomp)
        call der%ddz(tmp1_in_z,tmp2_in_z,this%z_bc(1),this%z_bc(2))
        call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
        qz = qz - this%kap*tmp1_in_y

        ! done inside SolidMixture
        !! If multispecies, add the inter-species enthalpy flux
        !if (this%mix%ns .GT. 1) then
        !    do i = 1,this%mix%ns
        !        call this%mix%material(i)%mat%get_enthalpy(this%T,tmp1_in_y)
        !        duidxj(:,:,:,qxidx) = duidxj(:,:,:,qxidx) + ( tmp1_in_y * Jx(:,:,:,i) )
        !        duidxj(:,:,:,qyidx) = duidxj(:,:,:,qyidx) + ( tmp1_in_y * Jy(:,:,:,i) )
        !        duidxj(:,:,:,qzidx) = duidxj(:,:,:,qzidx) + ( tmp1_in_y * Jz(:,:,:,i) )
        !    end do
        !end if

        ! Done
    end subroutine


end module 
