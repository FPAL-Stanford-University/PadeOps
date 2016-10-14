module SolidGrid
    use kind_parameters, only: rkind, clen
    use constants, only: zero,eps,third,half,one,two,three,four
    use FiltersMod, only: filters
    use GridMod, only: grid
    use gridtools, only: alloc_buffs, destroy_buffs
    use sgrid_hooks, only: meshgen, initfields, hook_output, hook_bc, hook_timestep, hook_source
    use decomp_2d, only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,  only: derivatives
    use StiffGasEOS,     only: stiffgas
    use Sep1SolidEOS,    only: sep1solid
    use GeneralMatEOS,   only: generaleos
   
    implicit none

    integer, parameter :: rho_index    = 1 
    integer, parameter :: u_index      = 2
    integer, parameter :: v_index      = 3
    integer, parameter :: w_index      = 4
    integer, parameter :: p_index      = 5
    integer, parameter :: T_index      = 6
    integer, parameter :: e_index      = 7
    integer, parameter :: mu_index     = 8
    integer, parameter :: bulk_index   = 9
    integer, parameter :: kap_index    = 10
    integer, parameter :: g11_index    = 11
    integer, parameter :: g12_index    = 12
    integer, parameter :: g13_index    = 13
    integer, parameter :: g21_index    = 14
    integer, parameter :: g22_index    = 15
    integer, parameter :: g23_index    = 16
    integer, parameter :: g31_index    = 17
    integer, parameter :: g32_index    = 18
    integer, parameter :: g33_index    = 19
    integer, parameter :: eel_index    = 20
    integer, parameter :: sxx_index    = 21
    integer, parameter :: sxy_index    = 22
    integer, parameter :: sxz_index    = 23
    integer, parameter :: syy_index    = 24
    integer, parameter :: syz_index    = 25
    integer, parameter :: szz_index    = 26
    integer, parameter :: Ent_index    = 27
    integer, parameter :: sos_index    = 28

    integer, parameter :: nfields = 28

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
       
        type(filters),   allocatable :: gfil

        type(stiffgas),  allocatable :: sgas
        type(sep1solid), allocatable :: elastic

        type(generaleos),  allocatable :: geneos
        integer                        :: eostype

        logical     :: plastic
        logical     :: explPlast
        real(rkind) :: tau0
        real(rkind) :: invtau0

        real(rkind), dimension(:,:,:,:), allocatable :: Wcnsrv                               ! Conserved variables
        real(rkind), dimension(:,:,:,:), allocatable :: xbuf, ybuf, zbuf   ! Buffers
       
        real(rkind) :: Cmu, Cbeta, Ckap

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
        real(rkind), dimension(:,:,:), pointer :: mu 
        real(rkind), dimension(:,:,:), pointer :: bulk 
        real(rkind), dimension(:,:,:), pointer :: kap

        real(rkind), dimension(:,:,:,:), pointer :: g
        real(rkind), dimension(:,:,:), pointer :: g11
        real(rkind), dimension(:,:,:), pointer :: g12
        real(rkind), dimension(:,:,:), pointer :: g13
        real(rkind), dimension(:,:,:), pointer :: g21
        real(rkind), dimension(:,:,:), pointer :: g22
        real(rkind), dimension(:,:,:), pointer :: g23
        real(rkind), dimension(:,:,:), pointer :: g31
        real(rkind), dimension(:,:,:), pointer :: g32
        real(rkind), dimension(:,:,:), pointer :: g33
       
        real(rkind), dimension(:,:,:,:), pointer :: devstress
        real(rkind), dimension(:,:,:), pointer :: sxx
        real(rkind), dimension(:,:,:), pointer :: sxy
        real(rkind), dimension(:,:,:), pointer :: sxz
        real(rkind), dimension(:,:,:), pointer :: syy
        real(rkind), dimension(:,:,:), pointer :: syz
        real(rkind), dimension(:,:,:), pointer :: szz
        
        real(rkind), dimension(:,:,:), pointer :: eel
        real(rkind), dimension(:,:,:), pointer :: Ent
        real(rkind), dimension(:,:,:), pointer :: sos
         
        logical     :: gfilttimes
        real(rkind) :: etafac
         
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
            procedure          :: getLAD
            procedure          :: filter
            procedure          :: getPhysicalProperties
            procedure, private :: get_tau
            procedure, private :: get_q
            procedure, private :: getPlasticSources
    end type

contains
    subroutine init(this, inputfile )
        use reductions, only: P_MAXVAL
        use exits,      only: message, warning, nancheck, GracefulExit
        class(sgrid),target, intent(inout) :: this
        character(len=clen), intent(in) :: inputfile  

        integer :: nx, ny, nz
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
        integer :: i, j, k, l 
        integer :: ioUnit
        real(rkind) :: gam = 1.4_rkind  !--- replace this block by eosparams
        real(rkind) :: Rgas = one
        real(rkind) :: PInf = zero
        real(rkind) :: shmod = zero
        real(rkind) :: yield = real(1.D30,rkind)
        real(rkind) :: tau0 = one       !--- since these are specific to separable EOS
        integer     :: eostype = 1
        real(rkind) :: eosparams(20)
        logical     :: explPlast = .FALSE.
        real(rkind) :: rho0 = one
        integer :: nsteps = -1
        real(rkind) :: dt = -one
        real(rkind) :: tstop = one
        real(rkind) :: CFL = -one
        logical :: SkewSymm = .FALSE.
        real(rkind) :: Cmu = 0.002_rkind
        real(rkind) :: Cbeta = 1.75_rkind
        real(rkind) :: Ckap = 0.01_rkind
        logical     :: plastic = .FALSE.
        integer     :: x_bc1 = 0, x_bcn = 0, y_bc1 = 0, y_bcn = 0, z_bc1 = 0, z_bcn = 0    ! 0: general, 1: symmetric/anti-symmetric
        real(rkind) :: etafac = zero
        logical     :: gfilttimes = .TRUE.

        character(len=clen) :: charout
        real(rkind), dimension(:,:,:,:), allocatable :: finger, fingersq
        real(rkind), dimension(:,:,:),   allocatable :: trG, trG2, detG

        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                             inputdir, outputdir, vizprefix, tviz, &
                                  periodicx, periodicy, periodicz, &
                         derivative_x, derivative_y, derivative_z, &
                                     filter_x, filter_y, filter_z, &
                                                       prow, pcol, &
                                                         SkewSymm  
        !namelist /SINPUT/  gam, Rgas, PInf, shmod, rho0, plastic, yield, &
        !                   explPlast, tau0, Cmu, Cbeta, Ckap,            &
        !                   x_bc1, x_bcn, y_bc1, y_bcn, z_bc1, z_bcn,     &
        !                   gfilttimes, etafac
        namelist /SINPUT/  rho0, eostype, eosparams, plastic, &
                           explPlast, Cmu, Cbeta, Ckap,            &
                           x_bc1, x_bcn, y_bc1, y_bcn, z_bc1, z_bcn,     &
                           gfilttimes, etafac

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

        this%eostype = eostype

        this%Cmu = Cmu
        this%Cbeta = Cbeta
        this%Ckap = Ckap

        this%plastic = plastic
        
        this%explPlast = explPlast

        this%gfilttimes = gfilttimes
        this%etafac = etafac

        ! Allocate decomp
        if ( allocated(this%decomp) ) deallocate(this%decomp)
        allocate(this%decomp)
        
        ! Initialize decomp
        call decomp_2d_init(nx, ny, nz, prow, pcol)
        call get_decomp_info(this%decomp)
        
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

        ! Allocate sgas
        if(this%eostype == 1) then
            ! separable eos
            gam   = eosparams(1);   Rgas  = eosparams(2);   PInf = eosparams(3);
            shmod = eosparams(4);    yield = eosparams(5);   tau0 = eosparams(6);
            this%tau0 = tau0    ! maybe move this to elastic or geneos
            if ( allocated(this%sgas) ) deallocate(this%sgas)
            allocate(this%sgas)
            call this%sgas%init(gam,Rgas,PInf)

            ! Allocate elastic
            if ( allocated(this%elastic) ) deallocate(this%elastic)
            allocate(this%elastic)
            call this%elastic%init(shmod,yield)
        else
            ! general eos
            this%tau0 = eosparams(9)    ! maybe move this to elastic or geneos
            if ( allocated(this%geneos) ) deallocate(this%geneos)
            allocate(this%geneos)
            call this%geneos%init(this%decomp,this%eostype,eosparams)
        endif

        ! Go to hooks if a different mesh is desired 
        call meshgen(this%decomp, this%dx, this%dy, this%dz, this%mesh) 

        ! Allocate fields
        if ( allocated(this%fields) ) deallocate(this%fields) 
        call alloc_buffs(this%fields,nfields,'y',this%decomp)
        call alloc_buffs(this%Wcnsrv,5,'y',this%decomp)
        
        ! Associate pointers for ease of use
        this%rho  => this%fields(:,:,:, rho_index) 
        this%u    => this%fields(:,:,:,   u_index) 
        this%v    => this%fields(:,:,:,   v_index) 
        this%w    => this%fields(:,:,:,   w_index)  
        this%p    => this%fields(:,:,:,   p_index)  
        this%T    => this%fields(:,:,:,   T_index)  
        this%e    => this%fields(:,:,:,   e_index)  
        this%mu   => this%fields(:,:,:,  mu_index)  
        this%bulk => this%fields(:,:,:,bulk_index)  
        this%kap  => this%fields(:,:,:, kap_index)   
       
        this%g    => this%fields(:,:,:,g11_index:g33_index)
        this%g11  => this%fields(:,:,:, g11_index)   
        this%g12  => this%fields(:,:,:, g12_index)   
        this%g13  => this%fields(:,:,:, g13_index)   
        this%g21  => this%fields(:,:,:, g21_index)   
        this%g22  => this%fields(:,:,:, g22_index)   
        this%g23  => this%fields(:,:,:, g23_index)   
        this%g31  => this%fields(:,:,:, g31_index)   
        this%g32  => this%fields(:,:,:, g32_index)   
        this%g33  => this%fields(:,:,:, g33_index)   
        
        this%devstress => this%fields(:,:,:,sxx_index:szz_index)
        this%sxx  => this%fields(:,:,:, sxx_index)   
        this%sxy  => this%fields(:,:,:, sxy_index)   
        this%sxz  => this%fields(:,:,:, sxz_index)   
        this%syy  => this%fields(:,:,:, syy_index)   
        this%syz  => this%fields(:,:,:, syz_index)   
        this%szz  => this%fields(:,:,:, szz_index)   
        
        this%eel  => this%fields(:,:,:, eel_index)   
        this%Ent  => this%fields(:,:,:, Ent_index)   
        this%sos  => this%fields(:,:,:, sos_index)   
       
        ! Initialize everything to a constant Zero
        this%fields = zero  

        ! Go to hooks if a different initialization is derired 
        !call initfields(this%decomp, this%dx, this%dy, this%dz, inputfile, this%mesh, this%fields, &
        !                rho0=this%rho0, mu=this%elastic%mu, gam=this%sgas%gam, PInf=this%sgas%PInf, &
        !                tstop=this%tstop, dt=this%dtfixed, tviz=tviz, yield=this%elastic%yield, &
        !                tau0=this%tau0)
        call initfields(this%decomp, this%dx, this%dy, this%dz, inputfile, this%mesh, this%fields, &
                        this%eostype, eosparams, rho0=this%rho0,&
                        tstop=this%tstop, dt=this%dtfixed, tviz=tviz)
       
        ! Get hydrodynamic and elastic energies 
        if(this%eostype == 1) then
            call this%sgas%get_e_from_p(this%rho,this%p,this%e)
        else
            ! energy already specified - not specifying primitive variables, so
            ! energy does not have to be computed from primitive variables
            call this%geneos%get_e_from_rhoT(this%rho0, this%g, this%rho, this%T, this%e)
        endif
        
        ! Check if the initialization was okay
        if ( nancheck(this%e) ) then
            call GracefulExit("NaN encountered at initialization in the hydrodynamic energy", 999)
        end if

        if(this%eostype == 1) then
            call alloc_buffs(finger,  6,"y",this%decomp)
            call alloc_buffs(fingersq,6,"y",this%decomp)
            allocate( trG (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )
            allocate( trG2(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )
            allocate( detG(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )

            call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG)
            call this%elastic%get_eelastic(this%rho0,trG,trG2,detG,this%eel) 
            call this%elastic%get_devstress(finger, fingersq, trG, trG2, detG, this%devstress)
        else
            allocate( detG(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )

            detG = this%g11*(this%g22*this%g33-this%g23*this%g32) &
                 - this%g12*(this%g21*this%g33-this%g31*this%g23) &
                 + this%g13*(this%g21*this%g32-this%g31*this%g22)

            call this%geneos%get_p_devstress_T_sos(this%rho0, this%g, this%rho, this%e, this%Ent, this%p, this%T, this%devstress, this%sos)
        endif
        
        ! Check if the initialization was okay
        if ( nancheck(this%eel) ) then
            call GracefulExit("NaN encountered at initialization in the elastic energy", 999)
        end if
        if ( nancheck(this%p) ) then
            call GracefulExit("NaN encountered at initialization in the pressure", 999)
        end if
        if ( nancheck(this%devstress) ) then
            call GracefulExit("NaN encountered at initialization in the deviatoric stress", 999)
        end if

        if(this%eostype == 1) then
            this%e = this%e + this%eel
            call this%sgas%get_T(this%e,this%T)
        else
            !call this%geneos%get_T(this%e, this%T) -- all clubbed together in get_p_devstress_T_sos
        endif

        if (P_MAXVAL(abs( this%rho/this%rho0/(detG) - one )) > 10._rkind*eps) then
            call warning("Inconsistent initialization: rho/rho0 and g are not compatible")
        end if

        if(this%eostype == 1) then
            deallocate( finger   )
            deallocate( fingersq )
            deallocate( trG      )
            deallocate( trG2     )
            deallocate( detG     )
        elseif(this%eostype == 3) then
            deallocate( detG     )
        endif

        ! Set all the attributes of the abstract grid type         
        this%outputdir = outputdir 
        
        this%periodicx = periodicx
        this%periodicy = periodicy
        this%periodicz = periodicz

        if ( ((x_bc1 /= 0) .AND. (x_bc1 /= 1)) ) call GracefulExit("x_bc1 can be only 0 (general) or 1 (symmetric)",4634)
        if ( ((x_bcn /= 0) .AND. (x_bcn /= 1)) ) call GracefulExit("x_bcn can be only 0 (general) or 1 (symmetric)",4634)
        if ( ((y_bc1 /= 0) .AND. (y_bc1 /= 1)) ) call GracefulExit("y_bc1 can be only 0 (general) or 1 (symmetric)",4634)
        if ( ((y_bcn /= 0) .AND. (y_bcn /= 1)) ) call GracefulExit("y_bcn can be only 0 (general) or 1 (symmetric)",4634)
        if ( ((z_bc1 /= 0) .AND. (z_bc1 /= 1)) ) call GracefulExit("z_bc1 can be only 0 (general) or 1 (symmetric)",4634)
        if ( ((z_bcn /= 0) .AND. (z_bcn /= 1)) ) call GracefulExit("z_bcn can be only 0 (general) or 1 (symmetric)",4634)
        this%x_bc = [x_bc1, x_bcn]
        this%y_bc = [y_bc1, y_bcn]
        this%z_bc = [z_bc1, z_bcn]

        this%derivative_x = derivative_x    
        this%derivative_y = derivative_y    
        this%derivative_z = derivative_z  

        this%filter_x = filter_x    
        this%filter_y = filter_y    
        this%filter_z = filter_z  
        
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

        ! Finally, set the local array dimensions
        this%nxp = this%decomp%ysz(1)
        this%nyp = this%decomp%ysz(2)
        this%nzp = this%decomp%ysz(3)


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
        varnames( 8) = 'mu'
        varnames( 9) = 'bulk'
        varnames(10) = 'kap'
        varnames(11) = 'g11'
        varnames(12) = 'g12'
        varnames(13) = 'g13'
        varnames(14) = 'g21'
        varnames(15) = 'g22'
        varnames(16) = 'g23'
        varnames(17) = 'g31'
        varnames(18) = 'g32'
        varnames(19) = 'g33'
        varnames(20) = 'e_elastic'
        varnames(21) = 'Sxx'
        varnames(22) = 'Sxy'
        varnames(23) = 'Sxz'
        varnames(24) = 'Syy'
        varnames(25) = 'Syz'
        varnames(26) = 'Szz'
        varnames(27) = 'Entr'

        allocate(this%viz)
        call this%viz%init(this%outputdir, vizprefix, nfields, varnames)
        this%tviz = tviz

        ! Do this here to keep tau0 and invtau0 compatible at the end of this subroutine
        this%invtau0 = one/this%tau0

        ! Check if the initialization was okay
        if ( nancheck(this%fields(:,:,:,8:27),i,j,k,l) ) then
            call message("fields: ",this%fields(i,j,k,l))
            write(charout,'(A,4(I5,A))') "NaN encountered in initialization ("//trim(varnames(l+7))//")  at (",i,", ",j,", ",k,", ",l,") of fields"
            call GracefulExit(trim(charout), 999)
        else
            call message("Initialization successful")
        end if
        
    end subroutine


    subroutine destroy(this)
        class(sgrid), intent(inout) :: this

        if(this%eostype == 1) then
        else
            call this%geneos%destroy()
        endif

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

        if (allocated(this%sgas)) deallocate(this%sgas) 
        if (allocated(this%elastic)) deallocate(this%elastic) 
        if (allocated(this%geneos))  deallocate(this%geneos) 
        
        if (allocated(this%Wcnsrv)) deallocate(this%Wcnsrv) 
        
        call this%viz%destroy()
        if (allocated(this%viz)) deallocate(this%viz)

        call decomp_2d_finalize
        if (allocated(this%decomp)) deallocate(this%decomp) 

        call hook_finalize()

    end subroutine

    subroutine gradient(this, f, dfdx, dfdy, dfdz, x_bc, y_bc, z_bc)
        class(sgrid),target, intent(inout) :: this
        real(rkind), intent(in), dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdx
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdy
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdz
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc

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
        call der%ddy(f,dfdy,y_bc(1),y_bc(2))

        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,decomp)
        call der%ddx(xtmp,xdum,x_bc(1),x_bc(2))
        call transpose_x_to_y(xdum,dfdx)

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,decomp)
        call der%ddz(ztmp,zdum,z_bc(1),z_bc(2))
        call transpose_z_to_y(zdum,dfdz)

    end subroutine 

    subroutine laplacian(this, f, lapf, x_bc, y_bc, z_bc)
        use timer
        class(sgrid),target, intent(inout) :: this
        real(rkind), intent(in), dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: lapf
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc
        
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
        call der%d2dy2(f,lapf,y_bc(1),y_bc(2))
        
        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,this%decomp) 
        call this%der%d2dx2(xtmp,xdum,x_bc(1),x_bc(2))
        call transpose_x_to_y(xdum,ytmp,this%decomp)

        lapf = lapf + ytmp

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,this%decomp)
        call this%der%d2dz2(ztmp,zdum,z_bc(1),z_bc(2))
        call transpose_z_to_y(zdum,ytmp,this%decomp)
        
        lapf = lapf + ytmp

    end subroutine

    subroutine simulate(this)
        use reductions, only: P_MEAN
        use timer,      only: tic, toc
        use exits,      only: GracefulExit, message
        class(sgrid), target, intent(inout) :: this

        logical :: tcond, vizcond, stepcond, hookcond = .FALSE.
        character(len=clen) :: stability
        real(rkind) :: cputime
        real(rkind), dimension(:,:,:,:), allocatable, target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        allocate( duidxj(this%nxp, this%nyp, this%nzp, 9) )
        ! Get artificial properties for initial conditions
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        call this%gradient(this%u,dudx,dudy,dudz,-this%x_bc, this%y_bc, this%z_bc)
        call this%gradient(this%v,dvdx,dvdy,dvdz, this%x_bc,-this%y_bc, this%z_bc)
        call this%gradient(this%w,dwdx,dwdy,dwdz, this%x_bc, this%y_bc,-this%z_bc)

        call this%getPhysicalProperties()

        call this%getLAD(dudx,dudy,dudz,&
                         dvdx,dvdy,dvdz,&
                         dwdx,dwdy,dwdz )
        deallocate( duidxj )
        ! ------------------------------------------------

        call this%get_dt(stability)

        ! Write out initial conditions
        call hook_output(this%decomp, this%der, this%fil, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%tsim, this%viz%vizcount, this%x_bc, this%y_bc, this%z_bc)
        call this%viz%WriteViz(this%decomp, this%mesh, this%fields, this%tsim)
        vizcond = .FALSE.
        
        ! Check for visualization condition and adjust time step
        if ( (this%tviz > zero) .AND. (this%tsim + this%dt >= this%tviz * this%viz%vizcount) ) then
            this%dt = this%tviz * this%viz%vizcount - this%tsim
            vizcond = .TRUE.
            stability = 'vizdump'
        end if

        tcond = .TRUE.
        ! Check tstop condition
        if ( (this%tstop > zero) .AND. (this%tsim >= this%tstop*(one-eps)) ) then
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
        do while ( tcond .AND. stepcond .AND. (.NOT. hookcond) )
            ! Advance time
            call tic()
            call this%advance_RK45()
            call toc(cputime)
            if (hookcond) stability = "hook"
            call message(1,"Time",this%tsim)
            call message(2,"Time step",this%dt)
            call message(2,"Stability limit: "//trim(stability))
            call message(2,"CPU time (in seconds)",cputime)
            call hook_timestep(this%decomp, this%der, this%mesh, this%fields, this%step, this%tsim, this%dt, this%x_bc, this%y_bc, this%z_bc, hookcond)
          
            ! Write out vizualization dump if vizcond is met 
            if (vizcond) then
                call hook_output(this%decomp, this%der, this%fil, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%tsim, this%viz%vizcount, this%x_bc, this%y_bc, this%z_bc)
                call this%viz%WriteViz(this%decomp, this%mesh, this%fields, this%tsim)
                vizcond = .FALSE.
            end if
            
            ! Get the new time step
            call this%get_dt(stability)
            
            ! Check for visualization condition and adjust time step
            if ( (this%tviz > zero) .AND. ((this%tsim + this%dt)*(one + eps) >= this%tviz * this%viz%vizcount) ) then
                this%dt = this%tviz * this%viz%vizcount - this%tsim
                vizcond = .TRUE.
            end if

            ! Check tstop condition
            if ( (this%tstop > zero) .AND. (this%tsim >= this%tstop*(one - eps) ) ) then
                tcond = .FALSE.
            else if ( (this%tstop > zero) .AND. (this%tsim + this%dt >= this%tstop) ) then
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
        class(sgrid), target, intent(inout) :: this

        real(rkind)                                          :: Qtmpt ! Temporary variable for RK45
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,5) :: rhs   ! RHS for conserved variables
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,5) :: Qtmp  ! Temporary variable for RK45
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: rhsg  ! RHS for g tensor equation
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: Qtmpg ! Temporary variable for RK45
        integer :: isub,i,j,k,l

        character(len=clen) :: charout

        Qtmp  = zero
        Qtmpg = zero
        Qtmpt = zero

        do isub = 1,RK45_steps
            call this%get_conserved()

            if ( nancheck(this%Wcnsrv,i,j,k,l) ) then
                call message("Wcnsrv: ",this%Wcnsrv(i,j,k,l))
                write(charout,'(A,I1,A,I5,A,4(I5,A))') "NaN encountered in solution (Wcnsrv) at substep ", isub, " of step ", this%step+1, " at (",i,", ",j,", ",k,", ",l,") of Wcnsrv"
                call GracefulExit(trim(charout), 999)
            end if
            if ( nancheck(this%g,i,j,k,l) ) then
                call message("g: ",this%g(i,j,k,l))
                write(charout,'(A,I1,A,I5,A,4(I5,A))') "NaN encountered in solution (g) at substep ", isub, " of step ", this%step+1, " at (",i,", ",j,", ",k,", ",l,") of Wcnsrv"
                call GracefulExit(trim(charout), 999)
            end if

            call this%getRHS(rhs, rhsg)
            Qtmp  = this%dt*rhs  + RK45_A(isub)*Qtmp
            Qtmpg = this%dt*rhsg + RK45_A(isub)*Qtmpg
            Qtmpt = this%dt + RK45_A(isub)*Qtmpt
            this%Wcnsrv = this%Wcnsrv + RK45_B(isub)*Qtmp
            this%g      = this%g      + RK45_B(isub)*Qtmpg
            this%tsim = this%tsim + RK45_B(isub)*Qtmpt

            ! Filter the conserved variables
            call this%filter(this%Wcnsrv(:,:,:,1), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)   ! continuity
            call this%filter(this%Wcnsrv(:,:,:,2), this%fil, 1,-this%x_bc, this%y_bc, this%z_bc)   ! x_mom
            call this%filter(this%Wcnsrv(:,:,:,3), this%fil, 1, this%x_bc,-this%y_bc, this%z_bc)   ! y_mom
            call this%filter(this%Wcnsrv(:,:,:,4), this%fil, 1, this%x_bc, this%y_bc,-this%z_bc)   ! z_mom
            call this%filter(this%Wcnsrv(:,:,:,5), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)   ! tot_energy

            ! Filter the g tensor
            if(this%gfilttimes) then
                call this%filter(this%g11(:,:,:), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
                call this%filter(this%g12(:,:,:), this%fil, 1,-this%x_bc,-this%y_bc, this%z_bc)
                call this%filter(this%g13(:,:,:), this%fil, 1,-this%x_bc, this%y_bc,-this%z_bc)
                call this%filter(this%g21(:,:,:), this%fil, 1,-this%x_bc,-this%y_bc, this%z_bc)
                call this%filter(this%g22(:,:,:), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
                call this%filter(this%g23(:,:,:), this%fil, 1, this%x_bc,-this%y_bc,-this%z_bc)
                call this%filter(this%g31(:,:,:), this%fil, 1,-this%x_bc, this%y_bc,-this%z_bc)
                call this%filter(this%g32(:,:,:), this%fil, 1, this%x_bc,-this%y_bc,-this%z_bc)
                call this%filter(this%g33(:,:,:), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
            end if
            
            call this%get_primitive()

            if (.NOT. this%explPlast .and. isub==RK45_steps) then
                if (this%plastic) then
                    ! Effect plastic deformations
                    if(this%eostype == 1) call this%elastic%plastic_deformation(this%devstress, this%dt, this%invtau0, this%g)
                    call this%get_primitive()      ! --- shouldn't this be after filtering?

                    !! Filter the conserved variables
                    !do i = 1,5
                    !    call this%filter(this%Wcnsrv(:,:,:,i), this%fil, 1)
                    !end do
                    !! Filter the g tensor
                    !do i = 1,9
                    !    call this%filter(this%g(:,:,:,i), this%fil, 1)
                    !end do

                    ! Filter the conserved variables
                    call this%filter(this%Wcnsrv(:,:,:,1), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)   ! continuity
                    call this%filter(this%Wcnsrv(:,:,:,2), this%fil, 1,-this%x_bc, this%y_bc, this%z_bc)   ! x_mom
                    call this%filter(this%Wcnsrv(:,:,:,3), this%fil, 1, this%x_bc,-this%y_bc, this%z_bc)   ! y_mom
                    call this%filter(this%Wcnsrv(:,:,:,4), this%fil, 1, this%x_bc, this%y_bc,-this%z_bc)   ! z_mom
                    call this%filter(this%Wcnsrv(:,:,:,5), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)   ! tot_energy

                    ! Filter the g tensor
                    if(this%gfilttimes) then
                        call this%filter(this%g11(:,:,:), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
                        call this%filter(this%g12(:,:,:), this%fil, 1,-this%x_bc,-this%y_bc, this%z_bc)
                        call this%filter(this%g13(:,:,:), this%fil, 1,-this%x_bc, this%y_bc,-this%z_bc)
                        call this%filter(this%g21(:,:,:), this%fil, 1,-this%x_bc,-this%y_bc, this%z_bc)
                        call this%filter(this%g22(:,:,:), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
                        call this%filter(this%g23(:,:,:), this%fil, 1, this%x_bc,-this%y_bc,-this%z_bc)
                        call this%filter(this%g31(:,:,:), this%fil, 1,-this%x_bc, this%y_bc,-this%z_bc)
                        call this%filter(this%g32(:,:,:), this%fil, 1, this%x_bc,-this%y_bc,-this%z_bc)
                        call this%filter(this%g33(:,:,:), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
                    end if
                end if
            end if
            
            call hook_bc(this%decomp, this%mesh, this%fields, this%tsim, this%x_bc, this%y_bc, this%z_bc)
            call this%post_bc()
        end do

        ! this%tsim = this%tsim + this%dt
        this%step = this%step + 1
            
    end subroutine

    subroutine get_dt(this,stability)
        use reductions, only : P_MAXVAL
        class(sgrid), target, intent(inout) :: this
        character(len=clen), intent(out) :: stability
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: cs
        real(rkind) :: dtCFL, dtmu, dtbulk, dtkap, dtplast

        if(this%eostype == 1) then
            call this%sgas%get_sos(this%rho,this%p,cs)  ! Speed of sound - hydrodynamic part
            call this%elastic%get_sos(this%rho0,cs)     ! Speed of sound - elastic part
        else
            !call this%geneos%get_sos(this%rho0,this%rho,this%devstress,cs)     ! Speed of sound -- all clubbed together in get_p_devstress_T_sos
            cs = this%sos
        endif

        dtCFL  = this%CFL / P_MAXVAL( ABS(this%u)/this%dx + ABS(this%v)/this%dy + ABS(this%w)/this%dz &
               + cs*sqrt( one/(this%dx**two) + one/(this%dy**two) + one/(this%dz**two) ))
        dtmu   = 0.2_rkind * min(this%dx,this%dy,this%dz)**2 / (P_MAXVAL( this%mu  / this%rho ) + eps)
        dtbulk = 0.2_rkind * min(this%dx,this%dy,this%dz)**2 / (P_MAXVAL( this%bulk/ this%rho ) + eps)
        dtkap  = 0.2_rkind * one / ( (P_MAXVAL( this%kap*this%T/(this%rho* (min(this%dx,this%dy,this%dz)**4))))**(third) + eps)
        dtplast = this%tau0

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
            else if ( this%dt > dtplast .and. this%explPlast .and. this%plastic) then
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

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6) :: finger, fingersq
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: trG, trG2, detG

        onebyrho => this%ybuf(:,:,:,1)

        this%rho  =  this%Wcnsrv(:,:,:,1)
        
        rhou => this%Wcnsrv(:,:,:,2)
        rhov => this%Wcnsrv(:,:,:,3)
        rhow => this%Wcnsrv(:,:,:,4)
        TE   => this%Wcnsrv(:,:,:,5)

        onebyrho = one/this%rho
        this%u = rhou * onebyrho
        this%v = rhov * onebyrho
        this%w = rhow * onebyrho
        this%e = (TE*onebyrho) - half*( this%u*this%u + this%v*this%v + this%w*this%w )
       
        if(this%eostype == 1) then 
            call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG)
            call this%elastic%get_eelastic(this%rho0,trG,trG2,detG,this%eel)
        
            call this%sgas%get_T(this%e,this%T)
            call this%sgas%get_p(this%rho,(this%e-this%eel),this%p)

            call this%elastic%get_devstress(finger, fingersq, trG, trG2, detG, this%devstress)
        else
            call this%geneos%get_p_devstress_T_sos(this%rho0, this%g, this%rho, this%e, this%Ent, this%p, this%T, this%devstress, this%sos)
            !call this%geneos%get_T(this%e, this%T) -- all clubbed together in get_p_devstress_T_sos
        endif

    end subroutine

    pure subroutine get_conserved(this)
        class(sgrid), intent(inout) :: this

        this%Wcnsrv(:,:,:,1) = this%rho
        this%Wcnsrv(:,:,:,2) = this%rho * this%u
        this%Wcnsrv(:,:,:,3) = this%rho * this%v
        this%Wcnsrv(:,:,:,4) = this%rho * this%w
        this%Wcnsrv(:,:,:,5) = this%rho * ( this%e + half*( this%u*this%u + this%v*this%v + this%w*this%w ) )

    end subroutine

    subroutine post_bc(this)
        class(sgrid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6) :: finger, fingersq
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: trG, trG2, detG

        if(this%eostype == 1) then 
            this%e = this%e - this%eel ! Get only hydrodynamic part

            call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG)
            call this%elastic%get_eelastic(this%rho0,trG,trG2,detG,this%eel)  ! Update elastic energy
            
            call this%sgas%get_e_from_p(this%rho,this%p,this%e)  ! Update hydrodynamic energy
            this%e = this%e + this%eel ! Combine both energies
            call this%sgas%get_T(this%e,this%T)  ! Get updated temperature

            call this%elastic%get_devstress(finger, fingersq, trG, trG2, detG, this%devstress)  ! Get updated stress
        else
            ! passive boundaries for now
        endif

    end subroutine

    subroutine getRHS(this, rhs, rhsg)
        use operators, only: curl
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,5), intent(out) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), intent(out) :: rhsg
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        real(rkind), dimension(:,:,:), pointer :: qx,qy,qz
        real(rkind), dimension(:,:,:), pointer :: penalty, tmp, detg
        real(rkind), dimension(:,:,:,:), pointer :: curlg
        !real(rkind) :: etafac = one/6._rkind

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        call this%gradient(this%u,dudx,dudy,dudz,-this%x_bc, this%y_bc, this%z_bc)
        call this%gradient(this%v,dvdx,dvdy,dvdz, this%x_bc,-this%y_bc, this%z_bc)
        call this%gradient(this%w,dwdx,dwdy,dwdz, this%x_bc, this%y_bc,-this%z_bc)

        call this%getPhysicalProperties()

        call this%getLAD(dudx,dudy,dudz,&
                         dvdx,dvdy,dvdz,&
                         dwdx,dwdy,dwdz )

        ! Get tau tensor and q (heat conduction) vector. Put in components of duidxj
        call this%get_tau( duidxj )
        call this%get_q  ( duidxj )

        ! Now, associate the pointers to understand what's going on better
        tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx); tauxz => duidxj(:,:,:,tauxzidx);
                                         tauyy => duidxj(:,:,:,tauyyidx); tauyz => duidxj(:,:,:,tauyzidx);
                                                                          tauzz => duidxj(:,:,:,tauzzidx);
       
        ! Add the deviatoric stress to the tau for use in fluxes 
        tauxx = tauxx + this%sxx; tauxy = tauxy + this%sxy; tauxz = tauxz + this%sxz
                                  tauyy = tauyy + this%syy; tauyz = tauyz + this%syz
                                                            tauzz = tauzz + this%szz
        
        qx => duidxj(:,:,:,qxidx); qy => duidxj(:,:,:,qyidx); qz => duidxj(:,:,:,qzidx);

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

        ! inverse deformation gradient tensor
        penalty => this%ybuf(:,:,:,1)
        tmp => this%ybuf(:,:,:,2)
        curlg => this%ybuf(:,:,:,3:5)
        detg => this%ybuf(:,:,:,6)

        detg = this%g11*(this%g22*this%g33-this%g23*this%g32) &
             - this%g12*(this%g21*this%g33-this%g31*this%g23) &
             + this%g13*(this%g21*this%g32-this%g31*this%g22)

        penalty = (this%etafac/this%dt)*(this%rho/detg/this%rho0-one)

        tmp = -this%u*this%g11-this%v*this%g12-this%w*this%g13
        call this%gradient(tmp,rhsg(:,:,:,1),rhsg(:,:,:,2),rhsg(:,:,:,3),-this%x_bc,this%y_bc,this%z_bc)
        
        call curl(this%decomp, this%der, this%g11, this%g12, this%g13, curlg,-this%x_bc,this%y_bc,this%z_bc)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + this%v*curlg(:,:,:,3) - this%w*curlg(:,:,:,2) + penalty*this%g11
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + this%w*curlg(:,:,:,1) - this%u*curlg(:,:,:,3) + penalty*this%g12
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + this%u*curlg(:,:,:,2) - this%v*curlg(:,:,:,1) + penalty*this%g13
 
        tmp = -this%u*this%g21-this%v*this%g22-this%w*this%g23
        call this%gradient(tmp,rhsg(:,:,:,4),rhsg(:,:,:,5),rhsg(:,:,:,6),this%x_bc,-this%y_bc,this%z_bc)
        
        call curl(this%decomp, this%der, this%g21, this%g22, this%g23, curlg,this%x_bc,-this%y_bc,this%z_bc)
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + this%v*curlg(:,:,:,3) - this%w*curlg(:,:,:,2) + penalty*this%g21
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + this%w*curlg(:,:,:,1) - this%u*curlg(:,:,:,3) + penalty*this%g22
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + this%u*curlg(:,:,:,2) - this%v*curlg(:,:,:,1) + penalty*this%g23
 
        tmp = -this%u*this%g31-this%v*this%g32-this%w*this%g33
        call this%gradient(tmp,rhsg(:,:,:,7),rhsg(:,:,:,8),rhsg(:,:,:,9),this%x_bc,this%y_bc,-this%z_bc)

        call curl(this%decomp, this%der, this%g31, this%g32, this%g33, curlg,this%x_bc,this%y_bc,-this%z_bc)
        rhsg(:,:,:,7) = rhsg(:,:,:,7) + this%v*curlg(:,:,:,3) - this%w*curlg(:,:,:,2) + penalty*this%g31
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + this%w*curlg(:,:,:,1) - this%u*curlg(:,:,:,3) + penalty*this%g32
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + this%u*curlg(:,:,:,2) - this%v*curlg(:,:,:,1) + penalty*this%g33

        if (this%explPlast) then
            call this%getPlasticSources(detg,rhsg)
        end if

        ! Call problem source hook
        call hook_source(this%decomp, this%mesh, this%fields, this%tsim, rhs, rhsg)
 
    end subroutine

    subroutine getRHS_x(       this,  rhs,&
                        tauxx,tauxy,tauxz,&
                            qx )
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, 5), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxx,tauxy,tauxz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qx

        real(rkind), dimension(:,:,:,:), pointer :: flux
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2

        flux => this%ybuf(:,:,:,1:5)
        xtmp1 => this%xbuf(:,:,:,1); xtmp2 => this%xbuf(:,:,:,2)

        flux(:,:,:,1) = this%Wcnsrv(:,:,:,2)   ! rho*u
        call transpose_y_to_x(flux(:,:,:,1),xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2))
        call transpose_x_to_y(xtmp2,flux(:,:,:,1),this%decomp)

        flux(:,:,:,2) = this%Wcnsrv(:,:,:,2)*this%u + this%p - tauxx
        call transpose_y_to_x(flux(:,:,:,2),xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
        call transpose_x_to_y(xtmp2,flux(:,:,:,2),this%decomp)

        flux(:,:,:,3) = this%Wcnsrv(:,:,:,2)*this%v          - tauxy
        call transpose_y_to_x(flux(:,:,:,3),xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2))
        call transpose_x_to_y(xtmp2,flux(:,:,:,3),this%decomp)

        flux(:,:,:,4) = this%Wcnsrv(:,:,:,2)*this%w          - tauxz
        call transpose_y_to_x(flux(:,:,:,4),xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2))
        call transpose_x_to_y(xtmp2,flux(:,:,:,4),this%decomp)

        flux(:,:,:,5) = (this%Wcnsrv(:,:,:,5) + this%p - tauxx)*this%u - this%v*tauxy - this%w*tauxz + qx
        call transpose_y_to_x(flux(:,:,:,5),xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2))
        call transpose_x_to_y(xtmp2,flux(:,:,:,5),this%decomp)

        ! Add to rhs
        rhs = rhs - flux

    end subroutine

    subroutine getRHS_y(       this,  rhs,&
                        tauxy,tauyy,tauyz,&
                            qy )
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, 5), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxy,tauyy,tauyz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qy

        real(rkind), dimension(:,:,:,:), pointer :: flux
        real(rkind), dimension(:,:,:), pointer :: ytmp1

        flux => this%ybuf(:,:,:,1:5)
        ytmp1 => this%ybuf(:,:,:,6)

        flux(:,:,:,1) = this%Wcnsrv(:,:,:,3)   ! rho*v
        call this%der%ddy(flux(:,:,:,1),ytmp1,-this%y_bc(1),-this%y_bc(2))
        rhs(:,:,:,1) = rhs(:,:,:,1) - ytmp1

        flux(:,:,:,2) = this%Wcnsrv(:,:,:,3)*this%u          - tauxy
        call this%der%ddy(flux(:,:,:,2),ytmp1,-this%y_bc(1),-this%y_bc(2))
        rhs(:,:,:,2) = rhs(:,:,:,2) - ytmp1

        flux(:,:,:,3) = this%Wcnsrv(:,:,:,3)*this%v + this%p - tauyy
        call this%der%ddy(flux(:,:,:,3),ytmp1,this%y_bc(1),this%y_bc(2))
        rhs(:,:,:,3) = rhs(:,:,:,3) - ytmp1

        flux(:,:,:,4) = this%Wcnsrv(:,:,:,3)*this%w          - tauyz
        call this%der%ddy(flux(:,:,:,4),ytmp1,-this%y_bc(1),-this%y_bc(2))
        rhs(:,:,:,4) = rhs(:,:,:,4) - ytmp1

        flux(:,:,:,5) = (this%Wcnsrv(:,:,:,5) + this%p - tauyy)*this%v - this%u*tauxy - this%w*tauyz + qy
        call this%der%ddy(flux(:,:,:,5),ytmp1,-this%y_bc(1),-this%y_bc(2))
        rhs(:,:,:,5) = rhs(:,:,:,5) - ytmp1

    end subroutine

    subroutine getRHS_z(       this,  rhs,&
                        tauxz,tauyz,tauzz,&
                            qz )
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, 5), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxz,tauyz,tauzz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qz

        real(rkind), dimension(:,:,:,:), pointer :: flux
        real(rkind), dimension(:,:,:), pointer :: ztmp1,ztmp2

        flux => this%ybuf(:,:,:,1:5)
        ztmp1 => this%zbuf(:,:,:,1); ztmp2 => this%zbuf(:,:,:,2)

        flux(:,:,:,1) = this%Wcnsrv(:,:,:,4)   ! rho*w
        call transpose_y_to_z(flux(:,:,:,1),ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2))
        call transpose_z_to_y(ztmp2,flux(:,:,:,1),this%decomp)

        flux(:,:,:,2) = this%Wcnsrv(:,:,:,4)*this%u          - tauxz
        call transpose_y_to_z(flux(:,:,:,2),ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2))
        call transpose_z_to_y(ztmp2,flux(:,:,:,2),this%decomp)

        flux(:,:,:,3) = this%Wcnsrv(:,:,:,4)*this%v          - tauyz
        call transpose_y_to_z(flux(:,:,:,3),ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2))
        call transpose_z_to_y(ztmp2,flux(:,:,:,3),this%decomp)

        flux(:,:,:,4) = this%Wcnsrv(:,:,:,4)*this%w + this%p - tauzz
        call transpose_y_to_z(flux(:,:,:,4),ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
        call transpose_z_to_y(ztmp2,flux(:,:,:,4),this%decomp)

        flux(:,:,:,5) = (this%Wcnsrv(:,:,:,5) + this%p - tauzz)*this%w - this%u*tauxz - this%v*tauyz + qz
        call transpose_y_to_z(flux(:,:,:,5),ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2))
        call transpose_z_to_y(ztmp2,flux(:,:,:,5),this%decomp)

        ! Add to rhs
        rhs = rhs - flux

    end subroutine

    subroutine getLAD(this,dudx,dudy,dudz,&
                           dvdx,dvdy,dvdz,&
                           dwdx,dwdy,dwdz )
        use reductions, only: P_MAXVAL
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: mustar,bulkstar,kapstar,func
        
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        real(rkind), dimension(:,:,:), pointer :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5
        real(rkind), dimension(:,:,:), pointer :: ztmp1,ztmp2

        xtmp1 => this%xbuf(:,:,:,1); xtmp2 => this%xbuf(:,:,:,2)
        
        ytmp1 => this%ybuf(:,:,:,1); ytmp2 => this%ybuf(:,:,:,2)
        ytmp3 => this%ybuf(:,:,:,3); ytmp4 => this%ybuf(:,:,:,4)
        ytmp5 => this%ybuf(:,:,:,5)
        
        ztmp1 => this%zbuf(:,:,:,1); ztmp2 => this%zbuf(:,:,:,2)

        ! -------- Artificial Shear Viscosity --------

        ! Magnitude of strain rate
        func = sqrt(dudx**2 + half*(dvdx+dudy)**2 + half*(dwdx+dudz)**2 &
                            +             dvdy**2 + half*(dwdy+dvdz)**2 &
                                                  +             dwdz**2 )
        
        ! Get 4th derivative in X
        call transpose_y_to_x(func,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,this%x_bc(1),this%x_bc(2))
        xtmp2 = xtmp1*this%dx**6
        call transpose_x_to_y(xtmp2,mustar,this%decomp)
        
        ! Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
        ztmp2 = ztmp1*this%dz**6
        call transpose_z_to_y(ztmp2,ytmp1,this%decomp)
        mustar = mustar + ytmp1
        
        ! Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp1,this%y_bc(1),this%y_bc(2))
        call this%der%d2dy2(ytmp1,ytmp2,this%y_bc(1),this%y_bc(2))
        ytmp1 = ytmp2*this%dy**6
        mustar = mustar + ytmp1

        mustar = this%Cmu*this%rho*abs(mustar)
        
        ! Filter mustar
        call this%filter(mustar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)
        
        ! -------- Artificial Bulk Viscosity --------
        
        func = dudx + dvdy + dwdz      ! dilatation
        
        ! Step 1: Get components of grad(rho) squared individually
        call this%gradient(this%rho,ytmp1,ytmp2,ytmp3,this%x_bc,this%y_bc,this%z_bc) ! Does not use any Y buffers
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(func,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,this%x_bc(1),this%x_bc(2))
        xtmp2 = xtmp1*this%dx**4
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        bulkstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
        ztmp2 = ztmp1*this%dz**4
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        bulkstar = bulkstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp4,this%y_bc(1),this%y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,this%y_bc(1),this%y_bc(2))
        ytmp4 = ytmp5*this%dy**4
        bulkstar = bulkstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2

        ! Now, all ytmps are free to use
        ytmp1 = dwdy-dvdz; ytmp2 = dudz-dwdx; ytmp3 = dvdx-dudy
        ytmp4 = ytmp1*ytmp1 + ytmp2*ytmp2 + ytmp3*ytmp3 ! |curl(u)|^2
        ytmp2 = func*func ! dilatation^2

        ! Calculate the switching function
        ytmp1 = ytmp2 / (ytmp2 + ytmp4 + real(1.0D-32,rkind)) ! Switching function f_sw
        where (func .GE. zero)
            ytmp1 = zero
        end where

        bulkstar = this%Cbeta*this%rho*ytmp1*abs(bulkstar)

        ! Filter bulkstar
        call this%filter(bulkstar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)

        ! -------- Artificial Conductivity --------

        ! Step 1: Get components of grad(e) squared individually
        call this%gradient(this%e,ytmp1,ytmp2,ytmp3,this%x_bc,this%y_bc,this%z_bc) ! Does not use any Y buffers
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(this%e,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,this%x_bc(1),this%x_bc(2))
        xtmp2 = xtmp1*this%dx**4
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        kapstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(this%e,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
        ztmp2 = ztmp1*this%dz**4
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        kapstar = kapstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(this%e,ytmp4,this%y_bc(1),this%y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,this%y_bc(1),this%y_bc(2))
        ytmp4 = ytmp5*this%dy**4
        kapstar = kapstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Now, all ytmps are free to us
        if(this%eostype == 1) then
            call this%sgas%get_sos(this%rho,this%p,ytmp1)  ! Speed of sound - hydrodynamic part
            call this%elastic%get_sos(this%rho0,ytmp1)     ! Speed of sound - elastic part
        else
            !call this%geneos%get_sos(this%rho0,this%rho,this%devstress,ytmp1)     ! Speed of sound -- all clubbed together in get_p_devstress_T_sos
        endif

        kapstar = this%Ckap*this%rho*ytmp1*abs(kapstar)/this%T

        ! Filter kapstar
        call this%filter(kapstar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)

        ! Now, add to physical fluid properties
        this%mu   = this%mu   + mustar
        this%bulk = this%bulk + bulkstar
        this%kap  = this%kap  + kapstar

    end subroutine

    subroutine filter(this,arr,myfil,numtimes,x_bc_,y_bc_,z_bc_)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(inout) :: arr
        type(filters), target, optional, intent(in) :: myfil
        integer, optional, intent(in) :: numtimes
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc
        
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
       
        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_
        
        ! First filter in y
        call fil2use%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
        ! Subsequent refilters 
        do idx = 1,times2fil-1
            arr = tmp_in_y
            call fil2use%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
        end do
        
        ! Then transpose to x
        call transpose_y_to_x(tmp_in_y,tmp1_in_x,this%decomp)

        ! First filter in x
        call fil2use%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_x = tmp2_in_x
            call fil2use%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
        end do 

        ! Now transpose back to y
        call transpose_x_to_y(tmp2_in_x,tmp_in_y,this%decomp)

        ! Now transpose to z
        call transpose_y_to_z(tmp_in_y,tmp1_in_z,this%decomp)

        !First filter in z
        call fil2use%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_z = tmp2_in_z
            call fil2use%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
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

    subroutine get_tau(this,duidxj)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), target, intent(inout) :: duidxj

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
        dudy = this%mu*dudy
        !tauxyidz = 2
    
        ! Step 2: Get tau_13 (dudz is destroyed)
        dudz = dudz + dwdx
        dudz = this%mu*dudz
        !tauxzidx = 3

        ! Step 3: Get tau_23 (dvdz is destroyed)
        dvdz = dvdz + dwdy
        dvdz = this%mu*dvdz
        !tauyzidx = 6

        ! Step 4: Get tau_11 (dvdx is destroyed)
        dvdx = bambda*dudx + lambda*(dvdy + dwdz)
        !tauxxidx = 4

        ! Step 5: Get tau_22 (dwdx is destroyed)
        dwdx = bambda*dvdy + lambda*(dudx + dwdz)
        !tauyyidx = 7

        ! Step 6: Get tau_33 (dwdy is destroyed)
        dwdy = bambda*dwdz + lambda*(dudx + dvdy)
        !tauzzidx = 8

        ! Done 
    end subroutine 

    subroutine get_q(this,duidxj)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(inout) :: duidxj

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
        duidxj(:,:,:,qyidx) = -this%kap*tmp1_in_y

        ! Step 2: Get qx (dudx is destroyed)
        call transpose_y_to_x(this%T,tmp1_in_x,this%decomp)
        call der%ddx(tmp1_in_x,tmp2_in_x,this%x_bc(1),this%x_bc(2))
        call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
        duidxj(:,:,:,qxidx) = -this%kap*tmp1_in_y

        ! Step 3: Get qz (dwdz is destroyed)
        call transpose_y_to_z(this%T,tmp1_in_z,this%decomp)
        call der%ddz(tmp1_in_z,tmp2_in_z,this%z_bc(1),this%z_bc(2))
        call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
        duidxj(:,:,:,qzidx) = -this%kap*tmp1_in_y

        ! Done
    end subroutine

    subroutine getPlasticSources(this, detg, rhsg)
        use constants, only: twothird
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)    :: detg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(inout) :: rhsg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: invtaurel

        ! Get S'S'
        invtaurel = this%sxx*this%sxx + two*this%sxy*this%sxy + two*this%sxz*this%sxz &
                                      +     this%syy*this%syy + two*this%syz*this%syz &
                                                              +     this%szz*this%szz

        ! 1/tau_rel
        invtaurel = this%invtau0 * ( invtaurel - (twothird)*this%elastic%yield**2 ) / this%elastic%mu**2
        where (invtaurel .LE. zero)
            invtaurel = zero
        end where
        invtaurel = invtaurel / (two * this%elastic%mu * detg)

        ! Add (1/tau_rel)*g*S to the rhsg (explicit plastic source terms)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + invtaurel * ( this%g11*this%sxx + this%g12*this%sxy + this%g13*this%sxz ) ! g11 
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + invtaurel * ( this%g11*this%sxy + this%g12*this%syy + this%g13*this%syz ) ! g12 
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + invtaurel * ( this%g11*this%sxz + this%g12*this%syz + this%g13*this%szz ) ! g13 
 
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + invtaurel * ( this%g21*this%sxx + this%g22*this%sxy + this%g23*this%sxz ) ! g21 
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + invtaurel * ( this%g21*this%sxy + this%g22*this%syy + this%g23*this%syz ) ! g22 
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + invtaurel * ( this%g21*this%sxz + this%g22*this%syz + this%g23*this%szz ) ! g23 

        rhsg(:,:,:,7) = rhsg(:,:,:,7) + invtaurel * ( this%g31*this%sxx + this%g32*this%sxy + this%g33*this%sxz ) ! g31 
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + invtaurel * ( this%g31*this%sxy + this%g32*this%syy + this%g33*this%syz ) ! g32 
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + invtaurel * ( this%g31*this%sxz + this%g32*this%syz + this%g33*this%szz ) ! g33 

    end subroutine
end module 
