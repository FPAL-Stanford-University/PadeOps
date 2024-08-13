module CurvilCompressibleGrid
    use kind_parameters,       only: rkind, clen
    use constants,             only: zero,eps,third,half,one,two,three,four
    use FiltersMod,            only: filters
    use GridMod,               only: grid
    use gridtools,             only: alloc_buffs, destroy_buffs
    use cvlgrid_hooks,         only: meshgen, initfields, hook_output, hook_bc, hook_timestep, hook_source
    use decomp_2d,             only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                                     transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,        only: derivatives
    use io_hdf5_stuff,         only: io_hdf5
    use IdealGasEOS,           only: idealgas
    use MixtureEOSMod,         only: mixture
    use ShearViscosityMod,     only: shearViscosity
    use TKEBudgetMod,          only: tkeBudget
    use ScaleDecompositionMod, only: scaleDecomposition
    use sgsmod_cgrid,          only: sgs_cgrid
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
    integer, parameter :: Ys_index     = 11

    integer            :: nfields      = 12
    integer            :: ncnsrv       = 5

    integer            :: mass_index = 1
    integer            ::  mom_index = 2
    integer            ::   TE_index = 5

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

    type, extends(grid) :: cvlgrid
       
        type(filters), allocatable :: gfil
        type(mixture), allocatable :: mix

        type(io_hdf5), allocatable :: viz
        type(io_hdf5), allocatable :: restart

        logical         :: compute_tke_budget
        type(tkeBudget) :: budget

        logical                  :: compute_scale_decomposition
        type(scaleDecomposition) :: scaledecomp
        
        logical                  :: forcing_mat 
        logical                  :: forcing_cha 
        real(rkind)              :: tsim_0, dtheta_0
        real(rkind), dimension(:,:,:,:), allocatable :: Wcnsrv             ! Conserved variables
        real(rkind), dimension(:,:,:,:), allocatable :: xbuf, ybuf, zbuf   ! Buffers
        real(rkind), dimension(:,:,:,:), allocatable :: meshcvl            ! Curvilinear mesh variables
        real(rkind), dimension(:,:,:,:), allocatable :: metric_multipliers ! Curvilinear mesh metrics
       
        real(rkind) :: Cmu, Cbeta, Ckap, Cdiff, CY

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
        real(rkind), dimension(:,:,:,:), pointer :: Ys
        real(rkind), dimension(:,:,:,:), pointer :: diff

        integer :: nrestart = 0
        integer :: vizramp  = 5

        ! SGS model
        logical :: useSGS
        real(rkind), allocatable, dimension(:,:,:,:) :: tausgs, Qjsgs
        type(sgs_cgrid), allocatable :: sgsmodel

        ! stretched and curvilinear meshes
        logical     :: xmetric=.false., ymetric=.false., zmetric=.false.
        real(rkind), dimension(:,:,:), pointer     :: xi, eta, zeta
        real(rkind), dimension(:,:,:), pointer     :: dxidx, detadx, dzetadx
        real(rkind), dimension(:,:,:), pointer     :: dxidy, detady, dzetady
        real(rkind), dimension(:,:,:), pointer     :: dxidz, detadz, dzetadz
        real(rkind), dimension(:,:,:), allocatable :: dxs, dys, dzs, InvJ
        real(rkind)                                :: dxi, deta, dzeta

        contains
            procedure          :: init
            procedure          :: init_curvilinear
            procedure          :: init_metric
            procedure          :: destroy_grid
            procedure          :: laplacian
            procedure          :: gradient 
            procedure          :: gradient_cvl 
            procedure          :: advance_RK45
            procedure          :: simulate
            procedure          :: get_dt
            ! procedure, private :: get_dt
            procedure, private :: get_primitive
            procedure, private :: get_conserved
            procedure, private :: post_bc
            procedure, private :: getRHS
            procedure, private :: getRHS_xi
            procedure, private :: getRHS_eta
            procedure, private :: getRHS_zeta
            procedure          :: getLAD
            procedure          :: filter
            procedure          :: getPhysicalProperties
            procedure, private :: get_tau
            procedure, private :: get_q
            procedure, private :: get_J
            procedure, private :: write_viz
            procedure, private :: write_restart
            procedure          :: read_restart
    end type

contains
    subroutine init(this, inputfile )
        use mpi
        use reductions, only: P_MAXVAL, P_MINVAL
        use exits, only: message, nancheck, GracefulExit
        class(cvlgrid),target, intent(inout) :: this
        character(len=clen), intent(in) :: inputfile  

        integer :: nx, ny, nz
        integer :: ns = 1
        character(len=clen) :: outputdir
        character(len=clen) :: inputdir
        character(len=clen) :: vizprefix = "cvlgrid"
        logical :: reduce_precision = .true.
        real(rkind) :: tviz = zero, Pr, Cp !! -- Pr, Cp need to be fixed
        character(len=clen), dimension(:), allocatable :: varnames
        logical :: periodicx = .true. 
        logical :: periodicy = .true. 
        logical :: periodicz = .true.
        integer :: x_bc1 = 0, x_bcn = 0
        integer :: y_bc1 = 0, y_bcn = 0
        integer :: z_bc1 = 0, z_bcn = 0
        character(len=clen) :: derivative_x = "cd10"  
        character(len=clen) :: derivative_y = "cd10" 
        character(len=clen) :: derivative_z = "cd10"
        character(len=clen) :: filter_x = "cf90"  
        character(len=clen) :: filter_y = "cf90" 
        character(len=clen) :: filter_z = "cf90"
        integer :: prow = 0, pcol = 0 
        integer :: i, j, k, l
        integer :: ioUnit
        real(rkind) :: gam = 1.4_rkind
        real(rkind) :: Rgas = one
        integer :: nsteps = -1
        real(rkind) :: dt = -one
        real(rkind) :: tstop = one
        real(rkind) :: CFL = -one
        logical :: SkewSymm = .FALSE.
        real(rkind) :: Cmu = 0.002_rkind
        real(rkind) :: Cbeta = 1.75_rkind
        real(rkind) :: Ckap = 0.01_rkind
        real(rkind) :: Cdiff = 0.003_rkind
        real(rkind) :: CY = 100.0_rkind
        character(len=clen) :: charout
        real(rkind) :: Ys_error
        logical     :: inviscid = .true.
        integer     :: nrestart = 0
        logical     :: rewrite_viz = .true.
        integer     :: vizramp = 5
        logical     :: compute_tke_budget = .false.
        logical     :: compute_scale_decomposition = .false.
        logical     :: forcing_mat = .false. ! KVM 2021
        logical     :: forcing_cha = .false. ! Vishwaja 2024
        logical     :: useSGS = .false.
        logical     :: xmetric=.false., ymetric=.false., zmetric=.false.

        namelist /INPUT/ nx, ny, nz, tstop, dt, CFL, nsteps, inputdir, &
                         outputdir, vizprefix, tviz, reduce_precision, &
                                      periodicx, periodicy, periodicz, &
                             derivative_x, derivative_y, derivative_z, &
                                         filter_x, filter_y, filter_z, &
                                          xmetric,  ymetric,  zmetric, &
                                                           prow, pcol, &
                                                             SkewSymm  
        namelist /CINPUT/  ns, gam, Rgas, Cmu, Cbeta, Ckap, Cdiff, CY, &
                             inviscid, nrestart, rewrite_viz, vizramp, &
                      compute_tke_budget, compute_scale_decomposition, &
                             x_bc1, x_bcn, y_bc1, y_bcn, z_bc1, z_bcn, &
                                                      forcing_mat, forcing_cha, useSGS


        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=INPUT)
        read(unit=ioUnit, NML=CINPUT)
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

        this%Cmu = Cmu
        this%Cbeta = Cbeta
        this%Ckap = Ckap
        this%Cdiff = Cdiff
        this%CY = CY

        this%compute_tke_budget = compute_tke_budget
        this%compute_scale_decomposition = compute_scale_decomposition
        this%forcing_mat = forcing_mat 
        this%forcing_cha = forcing_cha

        ! Allocate decomp
        if ( allocated(this%decomp) ) deallocate(this%decomp)
        allocate(this%decomp)
        
        ! Initialize decomp
        call decomp_2d_init(nx, ny, nz, prow, pcol)
        call get_decomp_info(this%decomp)
        
        ! Allocate mesh
        if ( allocated(this%mesh) ) deallocate(this%mesh) 
        call alloc_buffs(this%mesh,3,'y',this%decomp)
        call alloc_buffs(this%meshcvl,3,'y',this%decomp)
        
        ! Associate pointers for ease of use
        this%x    => this%mesh  (:,:,:, 1) 
        this%y    => this%mesh  (:,:,:, 2) 
        this%z    => this%mesh  (:,:,:, 3)

        this%xi   => this%meshcvl(:,:,:, 1) 
        this%eta  => this%meshcvl(:,:,:, 2) 
        this%zeta => this%meshcvl(:,:,:, 3)

        ! Generate default mesh: X \in [-1, 1), Y \in [-1, 1), Z \in [-1, 1)
        this%dxi   = two/nx
        this%deta  = two/ny
        this%dzeta = two/nz

        ! Generate default mesh
        do k = 1,size(this%mesh,3)
            do j = 1,size(this%mesh,2)
                do i = 1,size(this%mesh,1)
                    this%mesh(i,j,k,1) = -one + (this%decomp%yst(1) - 1 + i - 1)*this%dxi 
                    this%mesh(i,j,k,2) = -one + (this%decomp%yst(2) - 1 + j - 1)*this%deta
                    this%mesh(i,j,k,3) = -one + (this%decomp%yst(3) - 1 + k - 1)*this%dzeta
                end do
            end do
        end do

        ! setup metrics for curvilinear grids
        allocate(this%dxs(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)))
        allocate(this%dys(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)))
        allocate(this%dzs(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)))
        allocate(this%InvJ(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)))

        call alloc_buffs(this%metric_multipliers, 9, "y", this%decomp)
        this%dxidx   => this%metric_multipliers(:,:,:,1)
        this%dxidy   => this%metric_multipliers(:,:,:,2)
        this%dxidz   => this%metric_multipliers(:,:,:,3)
        this%detadx  => this%metric_multipliers(:,:,:,4)
        this%detady  => this%metric_multipliers(:,:,:,5)
        this%detadz  => this%metric_multipliers(:,:,:,6)
        this%dzetadx => this%metric_multipliers(:,:,:,7)
        this%dzetady => this%metric_multipliers(:,:,:,8)
        this%dzetadz => this%metric_multipliers(:,:,:,9)

        ! Allocate 2 buffers for each of the three decompositions
        call alloc_buffs(this%xbuf,nbufsx,"x",this%decomp)
        call alloc_buffs(this%ybuf,nbufsy,"y",this%decomp)
        call alloc_buffs(this%zbuf,nbufsz,"z",this%decomp)

        ! Go to hooks if a different mesh is desired
        call meshgen(this%decomp,this%dxi, this%deta, this%dzeta, this%mesh, inputfile, &
                     this%meshcvl, this%dxs, this%dys, this%dzs, this%xbuf, this%zbuf)

        !call this%init_curvilinear(inputfile)

        if (ns .LT. 1) call GracefulExit("Cannot have less than 1 species. Must have ns >= 1.",4568)

        ! Allocate mixture
        if (allocated(this%mix)) deallocate(this%mix)
        allocate(this%mix)
        call this%mix%init(this%decomp,ns,inviscid)
        ! allocate(this%mix , source=mixture(this%decomp,ns))

        ! Set default materials with the same gam and Rgas
        if (this%mix%inviscid) then
            do i = 1,ns
                call this%mix%set_material(i,idealgas(gam,Rgas))
            end do
        end if

        nfields = kap_index + 2*ns   ! Add ns massfractions to fields
        ncnsrv  = ncnsrv   + ns - 1  ! Add ns-1 conserved variables for the massfraction equations
        

        ! Set mass, momentum and energy indices in Wcnsrv
        mass_index = 1
        mom_index  = mass_index + ns
        TE_index   = mom_index + 3

        ! Allocate fields
        if ( allocated(this%fields) ) deallocate(this%fields) 
        call alloc_buffs(this%fields,nfields,'y',this%decomp)
        call alloc_buffs(this%Wcnsrv,ncnsrv ,'y',this%decomp)
        
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
        this%Ys   => this%fields(:,:,:,Ys_index:Ys_index+ns-1)
        this%diff => this%fields(:,:,:,Ys_index+ns:Ys_index+2*ns-1)

        !!! ----- block moved from here --------

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
        if ((this%filter_x=='none') .AND. (this%filter_y=='none') .AND. (this%filter_z=='none')) then
            call message("Note: Filtering turned off")
        end if

        ! Allocate der
        if ( allocated(this%der) ) deallocate(this%der)
        allocate(this%der)

        if ( allocated(this%der_metric) ) deallocate(this%der_metric)
        allocate(this%der_metric)
        ! Initialize derivatives 
        call this%der%init(                           this%decomp, &
                           this%dxi,    this%deta,     this%dzeta, &
                         periodicx,     periodicy,      periodicz, &
                      derivative_x,  derivative_y,   derivative_z)
                           ! this%x,        this%y,         this%z) !, &
                           !.false.,       .false.,        .false., &
                           !.false.,     inputfile,                 &
                           !this%xi,      this%eta,      this%zeta, this%xbuf, this%zbuf)

        call this%der_metric%init(                           this%decomp, &
                           this%dxi,    this%deta,     this%dzeta, &
                         .false.,     .false.,      .false., &
                      derivative_x,  derivative_y,   derivative_z)

        ! Allocate fil and gfil
        if ( allocated(this%fil) ) deallocate(this%fil)
        allocate(this%fil)
        if ( allocated(this%gfil) ) deallocate(this%gfil)
        allocate(this%gfil)
        
        ! Initialize filters
        if ((this%filter_x=='none') .AND. (this%filter_y=='none') .AND. (this%filter_z=='none')) then
            call this%fil%init( this%decomp, &
                         periodicx,   periodicy,    periodicz, &
                         'cf90',      'cf90',       'cf90'  ) 
        else
            call this%fil%init( this%decomp, &
                         periodicx,     periodicy,      periodicz, &
                          filter_x,      filter_y,       filter_z  ) 
        end if
        call this%gfil%init(                          this%decomp, &
                         periodicx,     periodicy,      periodicz, &
                        "gaussian",    "gaussian",     "gaussian"  )      

        this%SkewSymm = SkewSymm

        ! Initialize everything to a constant Zero
        this%fields = zero
        this%Ys(:,:,:,1) = one   ! So that all massfractions add up to unity

        !Cp = 3.5_rkind
        !Pr = 0.7_rkind  !!!needed to be fixed
        
        this%useSGS = useSGS

        if(this%useSGS) then
            call alloc_buffs(this%tausgs,6,'y',this%decomp)
            call alloc_buffs(this%Qjsgs, 3,'y',this%decomp)
           
            allocate(this%sgsmodel)
            
            call this%sgsmodel%init(this%der, this%decomp, Cp, Pr, this%dxi, this%deta, this%dzeta, inputfile, this%xbuf, this%ybuf, this%zbuf, periodicx, periodicy, periodicz, x_bc1, x_bcn, y_bc1, y_bcn, z_bc1, z_bcn, this%xmetric, this%ymetric, this%zmetric, this%dxs, this%dys,this%dzs)
         
        endif

        ! Finally, set the local array dimensions
        this%nxp = this%decomp%ysz(1)
        this%nyp = this%decomp%ysz(2)
        this%nzp = this%decomp%ysz(3)

        ! Find metrics for curvilinear
        call this%init_curvilinear(inputfile)
        ! Go to hooks if a different initialization is derired 
        call initfields(this%decomp, this%dxi, this%deta, this%dzeta, inputfile, this%mesh, this%fields, &
                        this%mix, this%tsim, this%tstop, this%dtfixed, tviz)
        
        ! Check for correct initialization of the mixture object
        call this%mix%check_initialization()

        ! Update mix
        call this%mix%update(this%Ys)
        call this%mix%get_e_from_p(this%rho,this%p,this%e)
        call this%mix%get_T(this%e,this%T)

        allocate(varnames(nfields))
        varnames(rho_index ) = 'density'
        varnames(u_index   ) = 'u'
        varnames(v_index   ) = 'v'
        varnames(w_index   ) = 'w'
        varnames(p_index   ) = 'p'
        varnames(T_index   ) = 'T'
        varnames(e_index   ) = 'e'
        varnames(mu_index  ) = 'mu'
        varnames(bulk_index) = 'bulk'
        varnames(kap_index ) = 'kap'
        do i = 1,ns
            write(charout,'(I2.2)') i
            varnames(Ys_index + i - 1     ) = 'Massfraction_'//trim(charout)
            varnames(Ys_index + i - 1 + ns) = 'Diffusivity_'//trim(charout)
        end do

        allocate(this%viz)
        call this%viz%init( mpi_comm_world, this%decomp, 'y', this%outputdir, vizprefix, &
                            reduce_precision=reduce_precision, read_only=.false., jump_to_last=.true.)
        this%tviz = tviz
        this%vizramp = vizramp
        ! Write mesh coordinates to file
        call this%viz%write_coords(this%mesh)

        ! Always write data in full precision for restart
        allocate(this%restart)
        call this%restart%init( mpi_comm_world, this%decomp, 'y', this%outputdir, 'restart', &
                                reduce_precision=.false., write_xdmf=.false., read_only=.true., jump_to_last=.true.)
        if (this%restart%vizcount >= 0) call this%read_restart(this%restart%vizcount)
        call this%restart%destroy()
        call this%restart%init( mpi_comm_world, this%decomp, 'y', this%outputdir, 'restart', &
                                reduce_precision=.false., write_xdmf=.false., read_only=.false., jump_to_last=.true.)
        this%nrestart = nrestart

        if ( ((this%step > 0) .or. (this%tsim > zero)) .and. (rewrite_viz)) then
            this%viz%vizcount = int(this%tsim / this%tviz + 100._rkind*eps) + 1
        end if

        ! Check for consistency of massfractions
        Ys_error = P_MAXVAL( abs(sum(this%Ys,4) - one) )
        if ( Ys_error .GT. 10._rkind*eps ) then
            write(charout,'(A,ES8.1E2)') "Inconsistency in massfractions (do not sum to unity). Maximum deviation from unity = ", Ys_error
            call GracefulExit(trim(charout),4387)
        end if

        ! Check if the initialization was okay
        if ( nancheck(this%fields,i,j,k,l) ) then
            call message("fields: ",this%fields(i,j,k,l))
            write(charout,'(A,4(I5,A))') "NaN encountered in initialization ("//trim(varnames(l))//") at (",i,", ",j,", ",k,", ",l,") of fields"
            call GracefulExit(trim(charout),999)
        end if

        deallocate(varnames)

        ! Instantiate TKE budget object
        if (this%compute_tke_budget) then
            ! if ( this%periodicx .or. this%periodicy .or. (.not. this%periodicz) ) then
            !     call GracefulExit("TKE budgets only supported for Z direction averaging. Sorry :(", 2083)
            ! end if
            
            ! if (allocated(this%budget)) deallocate(this%budget)
            ! allocate(this%budget , source=tkeBudget(this%decomp, this%der, this%mesh, this%dxi, this%deta, this%dzeta, &
            !                                         [this%periodicx, this%periodicy, this%periodicz], this%outputdir, &
            !                                         this%x_bc, this%y_bc, this%z_bc, reduce_precision))
            call this%budget%init(this%decomp, this%der, this%mesh, this%dxi, this%deta, this%dzeta, this%mix%ns, &
                                  [this%periodicx, this%periodicy, this%periodicz], this%outputdir, &
                                  this%x_bc, this%y_bc, this%z_bc, this%x, this%y, this%z, inputfile, &
                                  this%xi, this%eta, this%zeta, reduce_precision=.false.)
        end if

        ! Instantiate scale decomposition object
        if (this%compute_scale_decomposition) then
            call this%scaledecomp%init(this%decomp, this%der, this%gfil, this%mesh, this%dxi, this%deta, this%dzeta, this%mix%ns, &
                                       this%x_bc, this%y_bc, this%z_bc, inputfile)
        end if
        
        ! KVM 2021 Initialize for ddt(dtheta)
        if (this%forcing_mat) then
            !call this%budget%get_dtheta(this%decomp,this%mesh(:,:,:,2), &
            !    this%rho,this%u,this%dtheta_0)
            !this%tsim_0 = 0.d0
            call message(0,"Forcing terms ON")
        endif

    end subroutine


    subroutine destroy_grid(this)
        class(cvlgrid), intent(inout) :: this

        if(associated(this%dzetadz)) nullify(this%dzetadz)
        if(associated(this%dzetady)) nullify(this%dzetady)
        if(associated(this%dzetadx)) nullify(this%dzetadx)
        if(associated(this%detadz))  nullify(this%detadz)
        if(associated(this%detady))  nullify(this%detady)
        if(associated(this%detadx))  nullify(this%detadx)
        if(associated(this%dxidz))   nullify(this%dxidz)
        if(associated(this%dxidy))   nullify(this%dxidy)
        if(associated(this%dxidx))   nullify(this%dxidx)
        if (allocated(this%metric_multipliers)) deallocate(this%metric_multipliers)

        if (allocated(this%meshcvl)) deallocate(this%meshcvl) 
        if (allocated(this%mesh)) deallocate(this%mesh) 
        if (allocated(this%fields)) deallocate(this%fields) 
        
        if (this%useSGS) then
          call this%sgsmodel%destroy()
          deallocate(this%sgsmodel)

          if(allocated(this%tausgs)) deallocate(this%tausgs)
          if(allocated(this%Qjsgs )) deallocate(this%Qjsgs)
        end if

        if(allocated(this%dxs) ) deallocate(this%dxs)
        if(allocated(this%dys) ) deallocate(this%dys)
        if(allocated(this%dzs) ) deallocate(this%dzs)
        if(allocated(this%InvJ) ) deallocate(this%InvJ)

        call this%der%destroy()
        if (allocated(this%der)) deallocate(this%der) 
        
        call this%der_metric%destroy()
        if (allocated(this%der_metric)) deallocate(this%der_metric)
 
        call this%fil%destroy()
        if (allocated(this%fil)) deallocate(this%fil) 
        
        call this%gfil%destroy()
        if (allocated(this%gfil)) deallocate(this%gfil) 
        
        call destroy_buffs(this%xbuf)
        call destroy_buffs(this%ybuf)
        call destroy_buffs(this%zbuf)

        call this%mix%destroy()
        if (allocated(this%mix)) deallocate(this%mix)
        
        if (allocated(this%Wcnsrv)) deallocate(this%Wcnsrv) 
        
        ! if (allocated(this%budget)) deallocate(this%budget)

        call this%viz%destroy()
        if (allocated(this%viz)) deallocate(this%viz)

        call this%restart%destroy()
        if (allocated(this%restart)) deallocate(this%restart)

        call decomp_2d_finalize
        if (allocated(this%decomp)) deallocate(this%decomp) 

    end subroutine

    subroutine init_curvilinear(this, inputfile)
        use reductions,       only : P_MAXVAL,P_MINVAL
        use decomp_2d,        only: decomp_info, nrank
        use constants,        only: zero,eps,third,half,one,two,three,four,pi,eight
        class(cvlgrid), intent(inout) :: this
        character(len=* ) ,intent(in) :: inputfile

        real(rkind), allocatable, dimension(:,:) :: metric_params
        integer                                  ::  xmetric_flag=0,  ymetric_flag=0, zmetric_flag=0, curvil_flag=0, ierr
        real(rkind), dimension(:,:,:,:), allocatable, target :: dxidxij
        real(rkind), dimension(:,:,:), allocatable:: Jacob_inv
        real(rkind), dimension(:,:,:),   pointer  :: dxdxi,dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta
        integer :: i,j,k
        character(len=clen) :: outputfile,str
        namelist /METRICS/ xmetric_flag, ymetric_flag, zmetric_flag, metric_params, curvil_flag

        allocate(metric_params(3,5))    ! 3 :: (x,y,z); 5 :: max no of parameters
        metric_params = zero

        open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=123, NML=METRICS)
        close(123)


        if(curvil_flag==0) then !! uniform mesh
            this%metric_multipliers(:,:,:,1) = one
            this%metric_multipliers(:,:,:,2) = zero
            this%metric_multipliers(:,:,:,3) = zero
            this%metric_multipliers(:,:,:,4) = zero
            this%metric_multipliers(:,:,:,5) = one
            this%metric_multipliers(:,:,:,6) = zero
            this%metric_multipliers(:,:,:,7) = zero
            this%metric_multipliers(:,:,:,8) = zero
            this%metric_multipliers(:,:,:,9) = one
        elseif(curvil_flag==1)then  !! only stretching
            call this%init_metric(this%x, this%xi,  xmetric_flag, metric_params(1,:), this%dxidx)
            call this%init_metric(this%y, this%eta, ymetric_flag, metric_params(2,:), this%detady)
            call this%init_metric(this%z, this%zeta,zmetric_flag, metric_params(3,:), this%dzetadz)
            this%metric_multipliers(:,:,:,2) = zero
            this%metric_multipliers(:,:,:,3) = zero
            this%metric_multipliers(:,:,:,4) = zero
            this%metric_multipliers(:,:,:,6) = zero
            this%metric_multipliers(:,:,:,7) = zero
            this%metric_multipliers(:,:,:,8) = zero
        elseif(curvil_flag==2) then !! fully curvilinear
            allocate( dxidxij(this%decomp%ysz(1), this%decomp%ysz(2), this%decomp%ysz(3), 9) )
            allocate( Jacob_inv(this%decomp%ysz(1), this%decomp%ysz(2), this%decomp%ysz(3)) )
            dxdxi => dxidxij(:,:,:,1); dxdeta => dxidxij(:,:,:,2); dxdzeta => dxidxij(:,:,:,3);
            dydxi => dxidxij(:,:,:,4); dydeta => dxidxij(:,:,:,5); dydzeta => dxidxij(:,:,:,6);
            dzdxi => dxidxij(:,:,:,7); dzdeta => dxidxij(:,:,:,8); dzdzeta => dxidxij(:,:,:,9);
            !!!! one-sided stencil for metric terms !!!!        
            call this%gradient_cvl(this%x,dxdxi,dxdeta,dxdzeta,-[0,0], [0,0], [0,0])
            call this%gradient_cvl(this%y,dydxi,dydeta,dydzeta, [0,0],-[0,0], [0,0])
            call this%gradient_cvl(this%z,dzdxi,dzdeta,dzdzeta, [0,0], [0,0],-[0,0])
            
            !Jacob_inv  = dxdxi  *(dydeta *dzdzeta- dzdeta *dydzeta) + &
            !             dxdeta *(dydzeta*dzdxi  - dzdzeta*dydxi)   + &
            !             dxdzeta*(dydxi  *dzdeta - dzdxi  *dydeta)    !!Inverse of Jacobian

            !this%metric_multipliers(:,:,:,1) =  (dydeta*dzdzeta-dydzeta*dzdeta)/ (Jacob_inv + real(1.0D-32,rkind) )
            !this%metric_multipliers(:,:,:,2) = -(dxdeta*dzdzeta-dxdzeta*dzdeta)/ (Jacob_inv + real(1.0D-32,rkind) )
            !this%metric_multipliers(:,:,:,3) =  (dxdeta*dydzeta-dxdzeta*dydeta)/ (Jacob_inv + real(1.0D-32,rkind) )
            !this%metric_multipliers(:,:,:,4) = -(dydxi*dzdzeta-dydzeta*dzdxi)  / (Jacob_inv + real(1.0D-32,rkind) )
            !this%metric_multipliers(:,:,:,5) =  (dxdxi*dzdzeta-dxdzeta*dzdxi)  / (Jacob_inv + real(1.0D-32,rkind) )
            !this%metric_multipliers(:,:,:,6) = -(dxdxi*dydzeta-dxdzeta*dydxi)  / (Jacob_inv + real(1.0D-32,rkind) )
            !this%metric_multipliers(:,:,:,7) =  (dydxi*dzdeta-dydeta*dzdxi)    / (Jacob_inv + real(1.0D-32,rkind) )
            !this%metric_multipliers(:,:,:,8) = -(dxdxi*dzdeta-dxdeta*dzdxi)    / (Jacob_inv + real(1.0D-32,rkind) )
            !this%metric_multipliers(:,:,:,9) =  (dxdxi*dydeta-dxdeta*dydxi)    / (Jacob_inv + real(1.0D-32,rkind) )

            !!!!=============== Metrics for 2D ==============!!!!
            Jacob_inv  = dxdxi  * dydeta - dxdeta * dydxi
            
            this%metric_multipliers(:,:,:,1) =  dydeta/ (Jacob_inv + real(1.0D-32,rkind) )
            this%metric_multipliers(:,:,:,2) = -dxdeta/ (Jacob_inv + real(1.0D-32,rkind) )
            this%metric_multipliers(:,:,:,3) =  zero
            this%metric_multipliers(:,:,:,4) = -dydxi/ (Jacob_inv + real(1.0D-32,rkind) )
            this%metric_multipliers(:,:,:,5) =  dxdxi / (Jacob_inv + real(1.0D-32,rkind) )
            this%metric_multipliers(:,:,:,6) =  zero
            this%metric_multipliers(:,:,:,7) =  zero      
            this%metric_multipliers(:,:,:,8) =  zero      
            this%metric_multipliers(:,:,:,9) =  one     
  
            !!!!=============== Metrics for 2D ==============!!!! 

            !write(outputfile, '(a)') 'der_num_x.dat'
            !open(10,file=outputfile,status='unknown')
            !do i=1,this%decomp%xsz(1)
            !   write(10,'(e19.12,2x,e19.12,2x,e19.12)') this%x(i,1,1), dxidxij(i,this%decomp%ysz(2)/2,1,1), dxidxij(i,this%decomp%ysz(2)/2,1,2)
            !enddo
            !close(10)

            !write(outputfile, '(a)') 'der_num_y.dat'
            !open(10,file=outputfile,status='unknown')
            !do i=1,this%decomp%ysz(2)
            !   write(10,'(e19.12,2x,e19.12,2x,e19.12)') this%y(1,i,1), dxidxij(this%decomp%xsz(1)/2,i,1,4), dxidxij(this%decomp%xsz(1)/2,i,1,5)
            !enddo
            !close(10)
 
            deallocate(dxidxij)
            deallocate(Jacob_inv)
        endif
        
        this%InvJ = this%dxidx*(this%detady*this%dzetadz-this%dzetady*this%detadz) + &
                    this%dxidy*(this%detadz*this%dzetadx-this%dzetadz*this%detadx) + &
                    this%dxidz*(this%detadx*this%dzetady-this%dzetadx*this%detady)
        this%InvJ = 1 / (this%InvJ + real(1.0D-32,rkind) )
        
            !!!!========= Check whether Identities are satisfied ==========!!!!
            allocate( dxidxij(this%decomp%ysz(1), this%decomp%ysz(2), this%decomp%ysz(3), 9) )
            dxdxi => dxidxij(:,:,:,1); dxdeta => dxidxij(:,:,:,2); dxdzeta => dxidxij(:,:,:,3);
            dydxi => dxidxij(:,:,:,4); dydeta => dxidxij(:,:,:,5); dydzeta => dxidxij(:,:,:,6);
            dzdxi => dxidxij(:,:,:,7); dzdeta => dxidxij(:,:,:,8); dzdzeta => dxidxij(:,:,:,9);
            !! Check whether Identities are satisfied
            call this%gradient_cvl(this%dxidx*this%InvJ,dxdxi,dxdeta,dxdzeta,-[0,0], [0,0], [0,0])
            call this%gradient_cvl(this%detadx*this%InvJ,dydxi,dydeta,dydzeta, [0,0],-[0,0], [0,0])
            call this%gradient_cvl(this%dzetadx*this%InvJ,dzdxi,dzdeta,dzdzeta, [0,0], [0,0],-[0,0])
            !! Identity - 1
            dxdeta = dxdxi + dydeta !+ dzdzeta 
            print*, '>>> Max value of I1 <<<', P_MAXVAL(abs(dxidxij(:,:,:,2)))

            !!!!========= Check whether Identities are satisfied ==========!!!!
            call this%gradient_cvl(this%dxidy*this%InvJ,dxdxi,dxdeta,dxdzeta,-[0,0], [0,0], [0,0])
            call this%gradient_cvl(this%detady*this%InvJ,dydxi,dydeta,dydzeta, [0,0],-[0,0], [0,0])
            call this%gradient_cvl(this%dzetady*this%InvJ,dzdxi,dzdeta,dzdzeta, [0,0], [0,0],-[0,0])
            !! Identity - 2
            dxdeta = dxdxi + dydeta !+ dzdzeta
            print*, '>>> Max value of I2 <<<', P_MAXVAL(abs(dxidxij(:,:,:,2)))

            !!!!========= Check whether Identities are satisfied ==========!!!!
            call this%gradient_cvl(this%dxidz*this%InvJ,dxdxi,dxdeta,dxdzeta,-[0,0], [0,0], [0,0])
            call this%gradient_cvl(this%detadz*this%InvJ,dydxi,dydeta,dydzeta, [0,0],-[0,0], [0,0])
            call this%gradient_cvl(this%dzetadz*this%InvJ,dzdxi,dzdeta,dzdzeta, [0,0], [0,0],-[0,0])
            !! Identity - 3
            dxdeta = dxdxi + dydeta !+ dzdzeta
            print*, '>>> Max value of I3 <<<', P_MAXVAL(abs(dxidxij(:,:,:,2)))
            deallocate(dxidxij)

        deallocate(metric_params)   
    end subroutine


    subroutine init_metric(this, xstretch, xuniform, flag, params, dxudxs)
        use exits, only: GracefulExit
        class(cvlgrid),             intent(in) :: this
        integer,                    intent(in) :: flag
        real(rkind), dimension(5), intent(in)  :: params
        real(rkind), dimension(:,:,:), intent(in)  :: xstretch, xuniform
        real(rkind), dimension(:,:,:), intent(out) :: dxudxs
        real(rkind) :: xfocus, xfocus_adj, xtau, hh, xstart, num, num1, den, BB, BB2, xstretch_loc
        real(rkind) :: xbyxfocm1, xuniform_adj, alpha, beta
        integer     :: i, j, k

        if(flag==0) then
           do k = 1, this%decomp%ysz(3)
              do j = 1, this%decomp%ysz(2)
                 do i = 1, this%decomp%ysz(1)
                    xstretch_loc = xuniform(i,j,k)
                    dxudxs(i,j,k) = one
                    ! compare with xstretch specified in meshgen
                    if(abs(xstretch(i,j,k)-xstretch_loc) > 1.0d-12) then
                        !print '(3i5,1x,2(e19.12,1x))', i,j,k, xstretch(i,j,k), xstretch_loc
                        call GracefulExit("flag = 0; metric is not consistent with meshgen. Check details.", 21)
                    endif
                 enddo
              enddo
           enddo
        elseif(flag==1) then
           ! concentrate towards the center -- Pletcher, Tannehill, Anderson
           ! (Section 5.6, Transformation 3, pg. 332) 
           xfocus = params(1);  xtau   = params(2);  xstart = params(3); hh = params(4)
           !print '(5(e19.12,1x))', params(:)
           !print '(5(e19.12,1x))', xfocus, xtau, xstart, hh
           xfocus_adj = xfocus - xstart
           num = one + (xfocus_adj/hh) * (exp( xtau) - one)
           den = one + (xfocus_adj/hh) * (exp(-xtau) - one)
           BB = half/xtau*log(num/den)
           do k = 1, this%decomp%ysz(3)
              do j = 1, this%decomp%ysz(2)
                 do i = 1, this%decomp%ysz(1)
                    ! adjust for starting point
                    xuniform_adj = (xuniform(i,j,k) - xstart) !/ hh

                    ! stretched location
                    num = sinh(xtau*BB)
                    xstretch_loc = xfocus_adj * (one + sinh(xtau * (xuniform_adj/hh-BB))/num) + xstart
                    xbyxfocm1 = (xstretch_loc-xstart)/xfocus_adj - one

                    ! metric for first derivative
                    dxudxs(i,j,k) = num * hh / (xtau * xfocus_adj * sqrt(one + (xbyxfocm1 * num)**2))

                    !print '(i5,1x,4(e19.12,1x))', j, xstretch(1,j,1), xuniform(1,j,1), dxudxs(1,j,1)
                    ! compare with xstretch specified in meshgen
                    if(abs(xstretch(i,j,k)-xstretch_loc) > 1.0d-12) then
                        !print '(3i5,1x,2(e19.12,1x))', i,j,k, xstretch(i,j,k), xstretch_loc
                        call GracefulExit("flag = 1; metric is not consistent with meshgen. Check details.", 21)
                    endif
                 enddo
              enddo
           enddo
        elseif(flag==2) then
           ! concentrate towards the two ends(alpha = 0.5) or one side at xend(alpha = 0)
           ! concentrate towards the center -- Pletcher, Tannehill, Anderson
           ! (Section 5.6, Transformation 2, pg. 335) 
           alpha  = params(1);  beta   = params(2);  xstart = params(3); hh = params(4)
           !print '(5(e19.12,1x))', params(:)
           !print '(5(e19.12,1x))', alpha, beta, xstart, hh
           BB = (beta + 1) / (beta - 1)
           num1 = 2*beta * (2*alpha+1) * (1-alpha) / log(BB)
           do k = 1, this%decomp%ysz(3)
              do j = 1, this%decomp%ysz(2)
                 do i = 1, this%decomp%ysz(1)
                    ! adjust for starting point
                    xuniform_adj = (xuniform(i,j,k) - xstart) / hh
                    ! stretched location
                    BB2 = BB ** ( (xuniform_adj-alpha) / (1-alpha) )
                    num = ((beta+2*alpha)*BB2 - beta + 2*alpha ) * hh
                    xstretch_loc = num/( (2*alpha+1)*(1+BB2) )   + xstart

                    ! metric for first derivative
                    BB2 = 2*alpha - ((xstretch_loc-xstart)/hh) * (2*alpha+1)
                    dxudxs(i,j,k) = num1 / (beta**2 - BB2**2)

                    !print '(i5,1x,4(e19.12,1x))', j, xstretch(1,j,1), xuniform(1,j,1), dxudxs(1,j,1)
                    ! compare with xstretch specified in meshgen
                    if(abs(xstretch(i,j,k)-xstretch_loc) > 1.0d-12) then
                        print '(3i5,1x,2(e19.12,1x))', i,j,k, xstretch(i,j,k), xstretch_loc
                        call GracefulExit("flag = 2; metric is not consistent with meshgen. Check details.", 21)
                    endif
                 enddo
              enddo
           enddo
        elseif(flag==3) then
           ! concentrate at one side at xstart
           ! (Section 5.6, Transformation 1, pg. 334) 
           alpha  = params(1);  beta   = params(2);  xstart = params(3); hh = params(4)
           !print '(5(e19.12,1x))', params(:)
           !print '(5(e19.12,1x))', alpha, beta, xstart, hh
           BB = (beta + 1) / (beta - 1)
           num1 = (2*beta)/log(BB)
           do k = 1, this%decomp%ysz(3)
              do j = 1, this%decomp%ysz(2)
                 do i = 1, this%decomp%ysz(1)
                    ! adjust for starting point
                    xuniform_adj = xuniform(i,j,k) - xstart
                    ! stretched location
                    num = (beta+1) - (beta-1)*(BB**(1-xuniform_adj/hh))
                    den = (BB**(1-xuniform_adj/hh)) + 1
                    xstretch_loc = hh*(num/den) + xstart

                    ! metric for first derivative
                    BB2 = 1 - ((xstretch_loc-xstart)/hh)
                    dxudxs(i,j,k) = num1/(beta**2 - BB2**2)

                    if(abs(xstretch(i,j,k)-xstretch_loc) > 1.0d-12) then
                        print '(3i5,1x,2(e19.12,1x))', i,j,k, xstretch(i,j,k), xstretch_loc
                        call GracefulExit("flag = 3; metric is not consistent with meshgen. Check details.", 21)
                    endif
                 enddo
              enddo
           enddo

        elseif(flag==4) then
           ! finite-difference evaluation of metrics (reduces order of accuracy)
            call GracefulExit("flag = 4 (finite-difference evaluation of metrics) is incomplete right now",21)
        else
            call GracefulExit("metric flag must be 0-4. Check details.",21)
        endif
    end subroutine 

    subroutine gradient_cvl(this, f, dfdx, dfdy, dfdz, x_bc, y_bc, z_bc)
        class(cvlgrid),target, intent(inout) :: this
        real(rkind), intent(in), dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdx
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdy
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdz
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc

        type(derivatives), pointer :: der_metric
        type(decomp_info), pointer :: decomp
        real(rkind), dimension(:,:,:), pointer :: xtmp,xdum,ztmp,zdum
        real(rkind), dimension(:,:,:), pointer :: ytmp1, ytmp2, ytmp3
        
        der_metric => this%der_metric
        decomp => this%decomp
        xtmp   => this%xbuf(:,:,:,1)
        xdum   => this%xbuf(:,:,:,2)
        ztmp   => this%zbuf(:,:,:,1)
        zdum   => this%zbuf(:,:,:,2)

        ! Get Y derivatives
        call der_metric%ddeta(f,dfdy,y_bc(1),y_bc(2))

        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,decomp)
        call der_metric%ddxi(xtmp,xdum,x_bc(1),x_bc(2))
        call transpose_x_to_y(xdum,dfdx)

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,decomp)
        call der_metric%ddzeta(ztmp,zdum,z_bc(1),z_bc(2))
        call transpose_z_to_y(zdum,dfdz)

    end subroutine 

    subroutine gradient(this, f, dfdx, dfdy, dfdz, x_bc, y_bc, z_bc)
        class(cvlgrid),target, intent(inout) :: this
        real(rkind), intent(in), dimension(this%nxp, this%nyp, this%nzp) :: f
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdx
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdy
        real(rkind), intent(out), dimension(this%nxp, this%nyp, this%nzp) :: dfdz
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc

        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        real(rkind), dimension(:,:,:), pointer :: xtmp,xdum,ztmp,zdum
        real(rkind), dimension(:,:,:), pointer :: ytmp1, ytmp2, ytmp3
        
        der    => this%der
        decomp => this%decomp
        xtmp   => this%xbuf(:,:,:,1)
        xdum   => this%xbuf(:,:,:,2)
        ztmp   => this%zbuf(:,:,:,1)
        zdum   => this%zbuf(:,:,:,2)

        ytmp1  => this%ybuf(:,:,:,1)
        ytmp2  => this%ybuf(:,:,:,2)
        ytmp3  => this%ybuf(:,:,:,3)

        ! Get Y derivatives
        !call der%ddeta(f,dfdy,y_bc(1),y_bc(2))
        call der%ddeta(f,ytmp2,y_bc(1),y_bc(2))

        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,decomp)
        call der%ddxi(xtmp,xdum,x_bc(1),x_bc(2))
        !call transpose_x_to_y(xdum,dfdx)
        call transpose_x_to_y(xdum,ytmp1)

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,decomp)
        call der%ddzeta(ztmp,zdum,z_bc(1),z_bc(2))
        !call transpose_z_to_y(zdum,dfdz)
        call transpose_z_to_y(zdum,ytmp3)

        dfdx = this%dxidx*ytmp1 + this%detadx*ytmp2 + this%dzetadx*ytmp3 
        dfdy = this%dxidy*ytmp1 + this%detady*ytmp2 + this%dzetady*ytmp3 
        dfdz = this%dxidz*ytmp1 + this%detadz*ytmp2 + this%dzetadz*ytmp3 
    end subroutine 

    subroutine laplacian(this, f, lapf, x_bc, y_bc, z_bc)
        use timer
        class(cvlgrid),target, intent(inout) :: this
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
        call der%d2deta2(f,lapf,y_bc(1),y_bc(2))
        
        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,this%decomp) 
        call this%der%d2dxi2(xtmp,xdum,x_bc(1),x_bc(2))
        call transpose_x_to_y(xdum,ytmp,this%decomp)

        lapf = lapf + ytmp

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,this%decomp)
        call this%der%d2dzeta2(ztmp,zdum,z_bc(1),z_bc(2))
        call transpose_z_to_y(zdum,ytmp,this%decomp)
        
        lapf = lapf + ytmp

    end subroutine

    subroutine simulate(this)
        use reductions, only: P_MEAN
        use timer,      only: tic, toc
        use exits,      only: GracefulExit, message
        use reductions, only: P_MAXVAL, P_MINVAL
        use RKCoeffs,   only: RK45_steps
        class(cvlgrid), target, intent(inout) :: this

        logical :: tcond, vizcond, stepcond
        character(len=clen) :: stability
        real(rkind) :: cputime
        real(rkind), dimension(:,:,:,:), allocatable, target :: duidxj
        real(rkind), dimension(:,:,:,:), allocatable, target :: gradYs
        real(rkind), dimension(:,:,:,:), allocatable, target :: tauij
        real(rkind), dimension(:,:,:),   allocatable :: tke_old ! Temporary variable for tke rate
        real(rkind), dimension(:,:,:,:), allocatable :: tke_prefilter, tke_postfilter ! Temporary variable for tke dissipation
        real(rkind), dimension(:,:,:,:),   allocatable :: rhoPsi_old ! Temporary variable for Psi rate
        real(rkind), dimension(:,:,:,:,:), allocatable :: rhoPsi_prefilter, rhoPsi_postfilter ! Temporary variable for Psi dissipation
        real(rkind), dimension(:,:,:),     allocatable :: KE_L_old ! Temporary variable for large scale KE rate
        real(rkind), dimension(:,:,:,:),   allocatable :: KE_L_prefilter, KE_L_postfilter ! Temporary variable for large scale KE dissipation
        real(rkind), dimension(:,:,:),     allocatable :: KE_S_old ! Temporary variable for small scale KE rate
        real(rkind), dimension(:,:,:,:),   allocatable :: KE_S_prefilter, KE_S_postfilter ! Temporary variable for small scale KE dissipation
        real(rkind), dimension(:,:,:),     allocatable :: rhoPsi_SD_old        ! Temporary variable for large scale KE rate
        real(rkind), dimension(:,:,:,:),   allocatable :: rhoPsi_SD_prefilter  ! Temporary variable for large scale KE dissipation
        real(rkind), dimension(:,:,:,:),   allocatable :: rhoPsi_SD_postfilter ! Temporary variable for large scale KE dissipation

        real(rkind), dimension(:,:,:),   pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:,:), pointer :: dYsdx, dYsdy, dYsdz
        ! real(rkind) :: dt_tke
        integer :: i
        integer :: stepramp = 0

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv) :: rhs  ! RHS for conserved variables

        ! call this%get_dt(stability)

        allocate( duidxj(this%nxp, this%nyp, this%nzp, 9) )
        allocate( gradYs(this%nxp, this%nyp, this%nzp,3*this%mix%ns) )
        
        ! Get artificial properties for initial condition output file
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        call this%gradient(this%u,dudx,dudy,dudz,-this%x_bc, this%y_bc, this%z_bc)
        call this%gradient(this%v,dvdx,dvdy,dvdz, this%x_bc,-this%y_bc, this%z_bc)
        call this%gradient(this%w,dwdx,dwdy,dwdz, this%x_bc, this%y_bc,-this%z_bc)

        ! call this%getPhysicalProperties()
        call this%mix%get_transport_properties(this%p, this%T, this%Ys, this%mu, this%bulk, this%kap, this%diff)

        if (this%mix%ns .GT. 1) then
            dYsdx => gradYs(:,:,:,              1:  this%mix%ns)
            dYsdy => gradYs(:,:,:,  this%mix%ns+1:2*this%mix%ns)
            dYsdz => gradYs(:,:,:,2*this%mix%ns+1:3*this%mix%ns)

            do i = 1,this%mix%ns
                call this%gradient(this%Ys(:,:,:,i),dYsdx(:,:,:,i),dYsdy(:,:,:,i),dYsdz(:,:,:,i), this%x_bc, this%y_bc, this%z_bc)
            end do
        
            call this%getLAD(dudx, dudy, dudz,&
                             dvdx, dvdy, dvdz,&
                             dwdx, dwdy, dwdz,&
                            dYsdx,dYsdy,dYsdz )
        else
            call this%getLAD(dudx, dudy, dudz,&
                             dvdx, dvdy, dvdz,&
                             dwdx, dwdy, dwdz )
        end if

        deallocate( gradYs )

        ! Write out initial conditions
        if (this%viz%vizcount == 0) then 
            call hook_output(this%decomp, this%der, this%dxi, this%deta, this%dzeta, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
            call this%write_viz()

            if ((this%compute_tke_budget) .or. (this%compute_scale_decomposition)) then
                if (this%compute_tke_budget) then
                    allocate( tke_old       (this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3)) )
                    allocate( tke_prefilter (this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3),RK45_steps) )
                    allocate( tke_postfilter(this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3),RK45_steps) )
                    tke_old = zero
                    tke_prefilter = zero
                    tke_postfilter = zero

                    if (this%mix%ns > 1) then
                        allocate( rhoPsi_old       (this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3),this%mix%ns) )
                        allocate( rhoPsi_prefilter (this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3),this%mix%ns,RK45_steps) )
                        allocate( rhoPsi_postfilter(this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3),this%mix%ns,RK45_steps) )

                        rhoPsi_old = zero
                        rhoPsi_prefilter = zero
                        rhoPsi_postfilter = zero
                    end if
                end if

                if (this%compute_scale_decomposition) then
                    allocate( KE_L_old       (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )
                    allocate( KE_L_prefilter (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )
                    allocate( KE_L_postfilter(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )
                    allocate( KE_S_old       (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )
                    allocate( KE_S_prefilter (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )
                    allocate( KE_S_postfilter(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )

                    KE_L_old = zero
                    KE_L_prefilter = zero
                    KE_L_postfilter = zero
                    KE_S_old = zero
                    KE_S_prefilter = zero
                    KE_S_postfilter = zero

                    if (this%mix%ns > 1) then
                        allocate( rhoPsi_SD_old       (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )
                        allocate( rhoPsi_SD_prefilter (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )
                        allocate( rhoPsi_SD_postfilter(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )

                        rhoPsi_SD_old = zero
                        rhoPsi_SD_prefilter = zero
                        rhoPsi_SD_postfilter = zero
                    end if
                end if

                ! Get tau tensor and q (heat conduction) vector. Put in components of duidxj
                call this%get_tau( duidxj )      
                allocate( tauij (this%nxp,this%nyp,this%nzp,6) )
                ! Now, associate the pointers to understand what's going on better
                tauij(:,:,:,1) = duidxj(:,:,:,tauxxidx)
                tauij(:,:,:,2) = duidxj(:,:,:,tauxyidx)
                tauij(:,:,:,3) = duidxj(:,:,:,tauxzidx)
                tauij(:,:,:,4) = duidxj(:,:,:,tauyyidx)
                tauij(:,:,:,5) = duidxj(:,:,:,tauyzidx)
                tauij(:,:,:,6) = duidxj(:,:,:,tauzzidx)

                ! dt_tke = one
                if (this%compute_tke_budget) then
                    call this%budget%tke_budget(this%rho, this%u, this%v, this%w, this%p, tauij, &
                                                tke_old, tke_prefilter, tke_postfilter, this%tsim, this%dt)
                    if (this%mix%ns > 1) then
                        call this%budget%mixing_budget(this%rho, this%u, this%v, this%w, this%Ys, this%diff, &
                                                       rhoPsi_old, rhoPsi_prefilter, rhoPsi_postfilter, this%tsim, this%dt)
                    end if
                end if

                if (this%compute_scale_decomposition) then
                    call this%scaledecomp%tke_budget(this%rho, this%u, this%v, this%w, this%p, tauij, &
                         KE_L_old, KE_S_old, KE_L_prefilter, KE_L_postfilter, KE_S_prefilter, KE_S_postfilter, this%tsim, this%dt)
                    if (this%mix%ns > 1) then
                        call this%scaledecomp%mix_budget(this%rho, this%u, this%v, this%w, this%Ys, this%diff, &
                             rhoPsi_SD_old, rhoPsi_SD_prefilter, rhoPsi_SD_postfilter, this%tsim, this%dt)
                    end if
                end if

                deallocate( tauij  )

                if ( allocated(tke_old) )        deallocate( tke_old )
                if ( allocated(tke_prefilter) )  deallocate( tke_prefilter )
                if ( allocated(tke_postfilter) ) deallocate( tke_postfilter )
                if ( allocated(rhoPsi_old) )        deallocate( rhoPsi_old )
                if ( allocated(rhoPsi_prefilter) )  deallocate( rhoPsi_prefilter )
                if ( allocated(rhoPsi_postfilter) ) deallocate( rhoPsi_postfilter )

                if ( allocated( KE_L_old        ) ) deallocate( KE_L_old        )
                if ( allocated( KE_L_prefilter  ) ) deallocate( KE_L_prefilter  )
                if ( allocated( KE_L_postfilter ) ) deallocate( KE_L_postfilter )
                if ( allocated( KE_S_old        ) ) deallocate( KE_S_old        )
                if ( allocated( KE_S_prefilter  ) ) deallocate( KE_S_prefilter  )
                if ( allocated( KE_S_postfilter ) ) deallocate( KE_S_postfilter )

                if ( allocated( rhoPsi_SD_old        ) ) deallocate( rhoPsi_SD_old        )
                if ( allocated( rhoPsi_SD_prefilter  ) ) deallocate( rhoPsi_SD_prefilter  )
                if ( allocated( rhoPsi_SD_postfilter ) ) deallocate( rhoPsi_SD_postfilter )
            end if

        end if
        vizcond = .FALSE.
       
        deallocate( duidxj )
        ! ------------------------------------------------

        ! Write initial restart file
        if (this%restart%vizcount == 0) call this%write_restart()

        ! Get time step (after computing all artificial properties)
        call this%get_dt(stability)

        ! Check for visualization condition and adjust time step
        if ( (this%tviz > zero) .AND. (this%tsim + this%dt > this%tviz * this%viz%vizcount) ) then
            this%dt = this%tviz * this%viz%vizcount - this%tsim
            vizcond = .TRUE.
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
            call this%advance_RK45(rhs, vizcond)
            call toc(cputime)
            
            call message(1,"Time",this%tsim)
            call message(2,"Step",this%step)
            call message(2,"Time step",this%dt)
            call message(2,"Stability limit: "//trim(stability))
            call message(2,"CPU time (in seconds)",cputime)
            if(this%useSGS) then
              call hook_timestep(this%decomp, this%mesh, this%fields, this%mix, this%step, this%tsim, this%sgsmodel)
            else
              call hook_timestep(this%decomp, this%mesh, this%fields, this%mix, this%step, this%tsim)
            endif
          
            ! Write out vizualization dump if vizcond is met 
            if (vizcond) then
                call hook_output(this%decomp, this%der, this%dxi, this%deta, this%dzeta, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
                ! call this%viz%WriteViz(this%decomp, this%mesh, this%fields, this%tsim)
                call this%getRHS(rhs, .false.)
                call this%write_viz()
                vizcond = .FALSE.
            end if
            
            ! Write restart file
            if ( (this%nrestart > 0) .and. (mod(this%step,this%nrestart) == 0) ) call this%write_restart()

            ! Get the new time step
            call this%get_dt(stability)
                
            ! Check for visualization condition and adjust time step
            if (stepramp == 0) then
                if ( (this%tviz > zero) .AND. (this%tsim + this%vizramp*this%dt >= this%tviz * this%viz%vizcount) ) then
                    this%dt = (this%tviz * this%viz%vizcount - this%tsim)/real(this%vizramp,rkind)
                    stability = 'viz dump'
                    stepramp = stepramp + 1
                end if
            else if (stepramp == this%vizramp-1) then
                this%dt = (this%tviz * this%viz%vizcount - this%tsim)/real(this%vizramp-stepramp,rkind)
                stepramp = 0
                vizcond = .TRUE.
                stability = 'viz dump'
            else
                this%dt = (this%tviz * this%viz%vizcount - this%tsim)/real(this%vizramp-stepramp,rkind)
                stepramp = stepramp + 1
                stability = 'viz dump'
            end if

            ! Check tstop condition
            if ( (this%tstop > zero) .AND. (this%tsim >= this%tstop*(one - eps) ) ) then
                tcond = .FALSE.
            else if ( (this%tstop > zero) .AND. (this%tsim + this%dt >= this%tstop) ) then
                this%dt = this%tstop - this%tsim
                vizcond = .TRUE.
                stability = 'final'
            end if

            ! Check nsteps condition
            if ( (this%nsteps <= 0) .OR. (this%step < this%nsteps) ) then
                stepcond = .TRUE.
            else
                stepcond = .FALSE.
            end if

        end do

        ! Write final restart file
        call this%write_restart()

    end subroutine

    subroutine advance_RK45(this, rhs, vizcond)
        use RKCoeffs,   only: RK45_steps,RK45_A,RK45_B
        use exits,      only: message,nancheck,GracefulExit
        use reductions, only: P_MAXVAL, P_MINVAL
        use decomp_2d,        only: decomp_info, nrank
        class(cvlgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv), intent(inout) :: rhs  ! RHS for conserved variables
        logical,              intent(in)    :: vizcond

        real(rkind)                                               :: Qtmpt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv) :: Qtmp ! Temporary variable for RK45

        real(rkind), dimension(:,:,:,:), target, allocatable :: duidxj, tauij, gradYs
        real(rkind), dimension(:,:,:),     allocatable :: tke_old ! Temporary variable for tke rate
        real(rkind), dimension(:,:,:,:),   allocatable :: tke_prefilter, tke_postfilter ! Temporary variable for tke dissipation
        real(rkind), dimension(:,:,:,:),   allocatable :: rhoPsi_old ! Temporary variable for Psi rate
        real(rkind), dimension(:,:,:,:,:), allocatable :: rhoPsi_prefilter, rhoPsi_postfilter ! Temporary variable for Psi dissipation
        real(rkind), dimension(:,:,:),     allocatable :: KE_L_old ! Temporary variable for large scale KE rate
        real(rkind), dimension(:,:,:,:),   allocatable :: KE_L_prefilter, KE_L_postfilter ! Temporary variable for large scale KE dissipation
        real(rkind), dimension(:,:,:),     allocatable :: KE_S_old ! Temporary variable for small scale KE rate
        real(rkind), dimension(:,:,:,:),   allocatable :: KE_S_prefilter, KE_S_postfilter ! Temporary variable for small scale KE dissipation
        real(rkind), dimension(:,:,:),     allocatable :: rhoPsi_SD_old        ! Temporary variable for large scale KE rate
        real(rkind), dimension(:,:,:,:),   allocatable :: rhoPsi_SD_prefilter  ! Temporary variable for large scale KE dissipation
        real(rkind), dimension(:,:,:,:),   allocatable :: rhoPsi_SD_postfilter ! Temporary variable for large scale KE dissipation
        ! real(rkind)                                    :: dt_tke

        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:,:), pointer :: dYsdx, dYsdy, dYsdz

        integer :: isub,i,j,k,l
        logical :: newTimeStep

        character(len=clen) :: charout
        real(rkind), dimension(:,:,:), pointer :: ytmp1,ytmp2

        ytmp1 => this%ybuf(:,:,:,1); ytmp2 => this%ybuf(:,:,:,2)

        if (this%step == 0) then
            call this%mix%update(this%Ys)
            call this%mix%get_e_from_p(this%rho,this%p,this%e)
            call this%mix%get_T(this%e,this%T)
        end if

        if ( (vizcond) .and. (this%compute_tke_budget) ) then
            allocate( tke_old       (this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3)) )
            allocate( tke_prefilter (this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3),RK45_steps) )
            allocate( tke_postfilter(this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3),RK45_steps) )
            if (this%mix%ns > 1) then
                allocate( rhoPsi_old       (this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3),this%mix%ns) )
                allocate( rhoPsi_prefilter (this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3),this%mix%ns,RK45_steps) )
                allocate( rhoPsi_postfilter(this%budget%avg%avg_size(1),this%budget%avg%avg_size(2),this%budget%avg%avg_size(3),this%mix%ns,RK45_steps) )
            end if
        end if

        if ( (vizcond) .and. (this%compute_scale_decomposition) ) then
            allocate( KE_L_old       (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )
            allocate( KE_L_prefilter (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )
            allocate( KE_L_postfilter(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )
            allocate( KE_S_old       (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )
            allocate( KE_S_prefilter (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )
            allocate( KE_S_postfilter(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )

            if (this%mix%ns > 1) then
                allocate( rhoPsi_SD_old       (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )
                allocate( rhoPsi_SD_prefilter (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )
                allocate( rhoPsi_SD_postfilter(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),RK45_steps) )
            end if
        end if

        Qtmp = zero
        Qtmpt = zero

        if ((vizcond) .and. (this%compute_tke_budget)) then
            call this%budget%get_tke(this%rho, this%u, this%v, this%w, tke_old)
            if (this%mix%ns > 1) then
                call this%budget%get_rhoPsi_bar(this%rho, this%Ys, rhoPsi_old)
            end if
        end if

        if ((vizcond) .and. (this%compute_scale_decomposition)) then
            call this%scaledecomp%get_kinetic_energies(this%rho, this%u, this%v, this%w, KE_L_old, KE_S_old)
            if (this%mix%ns > 1) then
                call this%scaledecomp%get_rhoPsi(this%rho, this%Ys, rhoPsi_SD_old)
            end if
        end if


        do isub = 1,RK45_steps

            if(isub==1) then
                newTimeStep = .true.
            else
                newTimeStep = .false.
            endif

            call this%get_conserved()        

            if ( nancheck(this%Wcnsrv,i,j,k,l) ) then
                call message("Wcnsrv: ",this%Wcnsrv(i,j,k,l))
                write(charout,'(A,I0,A,I0,A,4(I0,A))') "NaN encountered in solution at substep ", isub, " of step ", this%step+1, " at Wcnsrv(",i,",",j,",",k,",",l,")"
                print *, charout
                call GracefulExit(trim(charout), 999)
            end if

            call this%getRHS(rhs, newTimeStep)       
            do i = 1,ncnsrv
               rhs(:,:,:,i) = rhs(:,:,:,i)/( this%InvJ(:,:,:) + real(1.0D-32,rkind) )
            end do
            Qtmp = this%dt*rhs + RK45_A(isub)*Qtmp
            Qtmpt = this%dt + RK45_A(isub)*Qtmpt
            this%Wcnsrv = this%Wcnsrv + RK45_B(isub)*Qtmp    
            this%tsim = this%tsim + RK45_B(isub)*Qtmpt

            ! if ( (vizcond) .and. (this%compute_tke_budget) .and. (isub == RK45_steps) ) then
            if ( (vizcond) .and. ((this%compute_tke_budget) .or. (this%compute_scale_decomposition)) ) then
                call this%get_primitive()
                if (this%compute_tke_budget) then
                    call this%budget%get_tke(this%rho, this%u, this%v, this%w, tke_prefilter(:,:,:,isub))
                    call this%budget%get_rhoPsi_bar(this%rho, this%Ys, rhoPsi_prefilter(:,:,:,:,isub))
                end if
                if (this%compute_scale_decomposition) then
                    call this%scaledecomp%get_kinetic_energies(this%rho, this%u, this%v, this%w, &
                                                  KE_L_prefilter(:,:,:,isub), KE_S_prefilter(:,:,:,isub))
                    if (this%mix%ns > 1) then
                        call this%scaledecomp%get_rhoPsi(this%rho, this%Ys, rhoPsi_SD_prefilter(:,:,:,isub))
                    end if
                end if
                ! dt_tke = RK45_B(isub)*Qtmpt
            end if

            ! Filter the conserved variables
            if (.not. ((this%filter_x=='none') .AND. (this%filter_y=='none') .AND.(this%filter_z=='none') )) then
                do i = 1,this%mix%ns
                    call this%filter(this%Wcnsrv(:,:,:,i), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
                end do
                call this%filter(this%Wcnsrv(:,:,:,mom_index  ), this%fil, 1,-this%x_bc, this%y_bc, this%z_bc)
                call this%filter(this%Wcnsrv(:,:,:,mom_index+1), this%fil, 1, this%x_bc,-this%y_bc, this%z_bc)
                call this%filter(this%Wcnsrv(:,:,:,mom_index+2), this%fil, 1, this%x_bc, this%y_bc,-this%z_bc)
                call this%filter(this%Wcnsrv(:,:,:, TE_index  ), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
            end if

            call this%get_primitive()      
            call hook_bc(this%decomp, this%mesh, this%fields, this%mix, this%tsim, &
                          this%x_bc, this%y_bc, this%z_bc, newTimeStep, this%step)
            call this%post_bc()
            !!!============Vishwaja channel flow isothermal walls=========!!!
            if (this%forcing_cha) then
            call hook_bc(this%decomp, this%mesh, this%fields, this%mix, this%tsim, &
                          this%x_bc, this%y_bc, this%z_bc, newTimeStep, this%step)
            end if
            ! Compute TKE budgets
            if ((vizcond) .and. ((this%compute_tke_budget) .or. (this%compute_scale_decomposition))) then
                if (this%compute_tke_budget) then
                    ! if (isub == RK45_steps-1) then
                    !     call this%budget%get_tke(this%rho, this%u, this%v, this%w, tke_old)
                    ! end if
                    call this%budget%get_tke(this%rho, this%u, this%v, this%w, tke_postfilter(:,:,:,isub))
                    call this%budget%get_rhoPsi_bar(this%rho, this%Ys, rhoPsi_postfilter(:,:,:,:,isub))
                end if
                if (this%compute_scale_decomposition) then
                    call this%scaledecomp%get_kinetic_energies(this%rho, this%u, this%v, this%w, &
                                                  KE_L_postfilter(:,:,:,isub), KE_S_postfilter(:,:,:,isub))
                    if (this%mix%ns > 1) then
                        call this%scaledecomp%get_rhoPsi(this%rho, this%Ys, rhoPsi_SD_postfilter(:,:,:,isub))
                    end if
                end if

                if (isub == RK45_steps)then
                    allocate( duidxj(this%nxp,this%nyp,this%nzp,9) )
                    allocate( gradYs(this%nxp, this%nyp, this%nzp,3*this%mix%ns) )

                    dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
                    dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
                    dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
                    
                    call this%gradient(this%u,dudx,dudy,dudz,-this%x_bc, this%y_bc, this%z_bc)
                    call this%gradient(this%v,dvdx,dvdy,dvdz, this%x_bc,-this%y_bc, this%z_bc)
                    call this%gradient(this%w,dwdx,dwdy,dwdz, this%x_bc, this%y_bc,-this%z_bc)

                    ! call this%getPhysicalProperties()
                    call this%mix%get_transport_properties(this%p, this%T, this%Ys, this%mu, this%bulk, this%kap, this%diff)

                    if (this%mix%ns .GT. 1) then
                        dYsdx => gradYs(:,:,:,1:this%mix%ns); dYsdy => gradYs(:,:,:,this%mix%ns+1:2*this%mix%ns); dYsdz => gradYs(:,:,:,2*this%mix%ns+1:3*this%mix%ns);
                        do i = 1,this%mix%ns
                            call this%gradient(this%Ys(:,:,:,i),dYsdx(:,:,:,i),dYsdy(:,:,:,i),dYsdz(:,:,:,i), this%x_bc, this%y_bc, this%z_bc)
                        end do
                        call this%getLAD(dudx, dudy, dudz,&
                                         dvdx, dvdy, dvdz,&
                                         dwdx, dwdy, dwdz,&
                                        dYsdx,dYsdy,dYsdz )
                    else
                        call this%getLAD(dudx, dudy, dudz,&
                                         dvdx, dvdy, dvdz,&
                                         dwdx, dwdy, dwdz )
                    end if
                    deallocate( gradYs )

                    ! Get tau tensor and q (heat conduction) vector. Put in components of duidxj
                    call this%get_tau( duidxj )
                    allocate( tauij (this%nxp,this%nyp,this%nzp,6) )
                    ! Now, associate the pointers to understand what's going on better
                    tauij(:,:,:,1) = duidxj(:,:,:,tauxxidx)
                    tauij(:,:,:,2) = duidxj(:,:,:,tauxyidx)
                    tauij(:,:,:,3) = duidxj(:,:,:,tauxzidx)
                    tauij(:,:,:,4) = duidxj(:,:,:,tauyyidx)
                    tauij(:,:,:,5) = duidxj(:,:,:,tauyzidx)
                    tauij(:,:,:,6) = duidxj(:,:,:,tauzzidx)
                    deallocate( duidxj )

                    if (this%compute_tke_budget) then
                        ! call this%budget%tke_budget(this%rho, this%u, this%v, this%w, this%p, tauij, &
                        !                             tke_old, tke_prefilter, this%tsim, dt_tke)
                        call this%budget%tke_budget(this%rho, this%u, this%v, this%w, this%p, tauij, &
                                                    tke_old, tke_prefilter, tke_postfilter, this%tsim, this%dt)
                        if (this%mix%ns > 1) then
                            call this%budget%mixing_budget(this%rho, this%u, this%v, this%w, this%Ys, this%diff, &
                                                           rhoPsi_old, rhoPsi_prefilter, rhoPsi_postfilter, this%tsim, this%dt)
                        end if
                    end if
                    if (this%compute_scale_decomposition) then
                        call this%scaledecomp%tke_budget(this%rho, this%u, this%v, this%w, this%p, tauij, &
                             KE_L_old, KE_S_old, KE_L_prefilter, KE_L_postfilter, KE_S_prefilter, KE_S_postfilter, this%tsim, this%dt)
                        if (this%mix%ns > 1) then
                            call this%scaledecomp%mix_budget(this%rho, this%u, this%v, this%w, this%Ys, this%diff, &
                                 rhoPsi_SD_old, rhoPsi_SD_prefilter, rhoPsi_SD_postfilter, this%tsim, this%dt)
                        end if
                    end if

                    deallocate( tauij  )

                end if
            end if
                
        end do
        
        ! KVM 2021 Update tsim_0,dthet_0
        if (this%forcing_mat) then
        !if (.true.) then
            call this%budget%get_dtheta(this%decomp,this%mesh(:,:,:,2), &
                this%rho,this%u,this%dtheta_0)
            this%tsim_0 = this%tsim 
            call message(1,"dtheta = ",this%dtheta_0)
        endif
            
        !this%tsim = this%tsim + this%dt
        this%step = this%step + 1
            
        if ( allocated(tke_old) )        deallocate( tke_old )
        if ( allocated(tke_prefilter) )  deallocate( tke_prefilter )
        if ( allocated(tke_postfilter) ) deallocate( tke_postfilter )
        if ( allocated(rhoPsi_old) )        deallocate( rhoPsi_old )
        if ( allocated(rhoPsi_prefilter) )  deallocate( rhoPsi_prefilter )
        if ( allocated(rhoPsi_postfilter) ) deallocate( rhoPsi_postfilter )

        if ( allocated( KE_L_old        ) ) deallocate( KE_L_old        )
        if ( allocated( KE_L_prefilter  ) ) deallocate( KE_L_prefilter  )
        if ( allocated( KE_L_postfilter ) ) deallocate( KE_L_postfilter )
        if ( allocated( KE_S_old        ) ) deallocate( KE_S_old        )
        if ( allocated( KE_S_prefilter  ) ) deallocate( KE_S_prefilter  )
        if ( allocated( KE_S_postfilter ) ) deallocate( KE_S_postfilter )

        if ( allocated( rhoPsi_SD_old        ) ) deallocate( rhoPsi_SD_old        )
        if ( allocated( rhoPsi_SD_prefilter  ) ) deallocate( rhoPsi_SD_prefilter  )
        if ( allocated( rhoPsi_SD_postfilter ) ) deallocate( rhoPsi_SD_postfilter )
    end subroutine


    subroutine get_dt(this,stability)
        use reductions, only : P_MAXVAL,P_MINVAL
        use decomp_2d,        only: decomp_info, nrank
        class(cvlgrid), target, intent(inout) :: this
        character(len=*), intent(out) :: stability
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: cs
        real(rkind) :: dtCFL, dtmu, dtbulk, dtkap, dtdiff

        call this%mix%get_sos(this%rho,this%p,cs)  ! Speed of sound - hydrodynamic part

        dtCFL  = this%CFL / P_MAXVAL( ABS(this%u)/this%dxs + ABS(this%v)/this%dys + ABS(this%w)/this%dzs &
               + cs*sqrt( one/(this%dxs**two) + one/(this%dys**two) + one/(this%dzs**two) ))
        dtmu   = 0.2_rkind * min(P_MINVAL(this%dxs),P_MINVAL(this%dys),P_MINVAL(this%dzs))**2 / (P_MAXVAL( this%mu  / this%rho ) + eps)
        dtbulk = 0.2_rkind * min(P_MINVAL(this%dxs),P_MINVAL(this%dys),P_MINVAL(this%dzs))**2 / (P_MAXVAL( this%bulk/ this%rho ) + eps)
        dtkap  = 0.2_rkind * one / ( (P_MAXVAL( this%kap*this%T/(this%rho* (min(P_MINVAL(this%dxs),P_MINVAL(this%dys),P_MINVAL(this%dzs))**4))))**(third) + eps)
        dtdiff = 0.2_rkind * min(P_MINVAL(this%dxs),P_MINVAL(this%dys),P_MINVAL(this%dzs))**2 / (P_MAXVAL( this%diff ) + eps)

        !print*, nrank, P_MINVAL(this%dxs), P_MINVAL(this%dys), P_MINVAL(this%dzs)

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
            end if

            if (this%step .LE. 10) then
                this%dt = this%dt / 10._rkind
                stability = 'startup'
            end if
        end if

    end subroutine

    pure subroutine get_primitive(this)
        class(cvlgrid), target, intent(inout) :: this
        real(rkind), dimension(:,:,:), pointer :: onebyrho
        real(rkind), dimension(:,:,:), pointer :: rhou,rhov,rhow,TE
        integer :: i

        onebyrho => this%ybuf(:,:,:,1)

        select case(this%mix%ns)
        case(1)
            this%rho = this%Wcnsrv(:,:,:,1)
        case default
            this%rho = sum( this%Wcnsrv(:,:,:,1:this%mix%ns),4 )
        end select
        
        rhou => this%Wcnsrv(:,:,:,mom_index  )
        rhov => this%Wcnsrv(:,:,:,mom_index+1)
        rhow => this%Wcnsrv(:,:,:,mom_index+2)
        TE   => this%Wcnsrv(:,:,:, TE_index  )

        onebyrho = one/this%rho
        do i = 1,this%mix%ns
            this%Ys(:,:,:,i) = this%Wcnsrv(:,:,:,i) * onebyrho
        end do
        this%u = rhou * onebyrho
        this%v = rhov * onebyrho
        this%w = rhow * onebyrho
        this%e = (TE*onebyrho) - half*( this%u*this%u + this%v*this%v + this%w*this%w )
        
        call this%mix%update(this%Ys)
        call this%mix%get_p(this%rho,this%e,this%p)
        call this%mix%get_T(this%e,this%T)

    end subroutine

    pure subroutine get_conserved(this)
        class(cvlgrid), intent(inout) :: this
        integer :: i

        do i = 1,this%mix%ns
            this%Wcnsrv(:,:,:,i) = this%rho*this%Ys(:,:,:,i)
        end do
        this%Wcnsrv(:,:,:,mom_index  ) = this%rho * this%u
        this%Wcnsrv(:,:,:,mom_index+1) = this%rho * this%v
        this%Wcnsrv(:,:,:,mom_index+2) = this%rho * this%w
        this%Wcnsrv(:,:,:, TE_index  ) = this%rho*( this%e + half*( this%u*this%u + this%v*this%v + this%w*this%w ) )

    end subroutine

    subroutine post_bc(this)
        class(cvlgrid), intent(inout) :: this

        call this%mix%update(this%Ys)
        call this%mix%get_e_from_p(this%rho,this%p,this%e)
        call this%mix%get_T(this%e,this%T)

    end subroutine

    subroutine getRHS(this, rhs, newTimeStep)
        class(cvlgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,ncnsrv), intent(out) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,3), target :: gradT
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,3*this%mix%ns), target :: gradYs
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:,:), pointer :: dYsdx, dYsdy, dYsdz
        real(rkind), dimension(:,:,:), pointer :: dTdx, dTdy, dTdz
        real(rkind), dimension(:,:,:), pointer :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        real(rkind), dimension(:,:,:), pointer :: qx,qy,qz
        real(rkind), dimension(:,:,:,:), pointer :: Jx,Jy,Jz
        integer :: i
        logical :: newTimeStep

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        dTdx =>  gradT(:,:,:,1); dTdy =>  gradT(:,:,:,2); dTdz =>  gradT(:,:,:,3);
        
        call this%gradient(this%u,dudx,dudy,dudz,-this%x_bc, this%y_bc, this%z_bc)
        call this%gradient(this%v,dvdx,dvdy,dvdz, this%x_bc,-this%y_bc, this%z_bc)
        call this%gradient(this%w,dwdx,dwdy,dwdz, this%x_bc, this%y_bc,-this%z_bc)
        call this%gradient(this%T,dTdx,dTdy,dTdz, this%x_bc, this%y_bc, this%z_bc)

        ! call this%getPhysicalProperties()
        call this%mix%get_transport_properties(this%p, this%T, this%Ys, this%mu, this%bulk, this%kap, this%diff)

        if(this%useSGS) then
            !call this%sgsmodel%getTauSGS(newTimeStep, duidxj, this%rho, this%u, this%v, this%w, this%tausgs)
            !call this%sgsmodel%getQjSGS (newTimeStep, duidxj, this%rho, this%u, this%v, this%w, this%T, gradT, this%Qjsgs )
            call this%sgsmodel%getTauijQjSGS (newTimeStep, duidxj, gradT, this%rho, &
                            this%u, this%v, this%w, this%T, this%tausgs, this%Qjsgs )
        endif

        if (this%mix%ns .GT. 1) then
            dYsdx => gradYs(:,:,:,1:this%mix%ns); dYsdy => gradYs(:,:,:,this%mix%ns+1:2*this%mix%ns); dYsdz => gradYs(:,:,:,2*this%mix%ns+1:3*this%mix%ns);
            do i = 1,this%mix%ns
                call this%gradient(this%Ys(:,:,:,i),dYsdx(:,:,:,i),dYsdy(:,:,:,i),dYsdz(:,:,:,i), this%x_bc, this%y_bc, this%z_bc)
            end do
            call this%getLAD(dudx, dudy, dudz,&
                             dvdx, dvdy, dvdz,&
                             dwdx, dwdy, dwdz,&
                            dYsdx,dYsdy,dYsdz )
        else
            call this%getLAD(dudx, dudy, dudz,&
                             dvdx, dvdy, dvdz,&
                             dwdx, dwdy, dwdz )
        end if
        
        ! Get tau tensor and q (heat conduction) vector. Put in components of duidxj
        call this%get_tau( duidxj )
        ! Now, associate the pointers to understand what's going on better
        tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx); tauxz => duidxj(:,:,:,tauxzidx);
                                         tauyy => duidxj(:,:,:,tauyyidx); tauyz => duidxj(:,:,:,tauyzidx);
                                                                          tauzz => duidxj(:,:,:,tauzzidx);

        ! Get species mass fluxes. Put in components on gradYs
        if (this%mix%ns .GT. 1) then
            call this%get_J( gradYs )
        end if
        Jx => gradYs(:,:,:,1:this%mix%ns); Jy => gradYs(:,:,:,this%mix%ns+1:2*this%mix%ns); Jz => gradYs(:,:,:,2*this%mix%ns+1:3*this%mix%ns);

        call this%get_q  ( duidxj, Jx, Jy, Jz )
        qx => duidxj(:,:,:,qxidx); qy => duidxj(:,:,:,qyidx); qz => duidxj(:,:,:,qzidx);

        if(this%useSGS) then
            tauxx = tauxx - this%tausgs(:,:,:,1)
            tauxy = tauxy - this%tausgs(:,:,:,2)
            tauxz = tauxz - this%tausgs(:,:,:,3)
            tauyy = tauyy - this%tausgs(:,:,:,4)
            tauyz = tauyz - this%tausgs(:,:,:,5)
            tauzz = tauzz - this%tausgs(:,:,:,6)

            qx = qx + this%Qjsgs(:,:,:,1)
            qy = qy + this%Qjsgs(:,:,:,2)
            qz = qz + this%Qjsgs(:,:,:,3)
        endif

        rhs = zero
        call this%getRHS_xi  (              rhs,&
                           tauxx,tauxy,tauxz,tauyy,tauyz,tauzz,&
                               qx,qy,qz,Jx,Jy,Jz )

        call this%getRHS_eta (              rhs,&
                           tauxx,tauxy,tauxz,tauyy,tauyz,tauzz,&
                               qx,qy,qz,Jx,Jy,Jz )

        call this%getRHS_zeta(              rhs,&
                           tauxx,tauxy,tauxz,tauyy,tauyz,tauzz,&
                               qx,qy,qz,Jx,Jy,Jz )

        ! KVM 2021 Call problem source hook
        if (this%forcing_mat) then
            !call hook_source(this%decomp, this%mesh, this%fields, &
            !    this%mix, this%tsim, rhs, &
            !    this%Wcnsrv,this%budget,this%tsim_0,this%dtheta_0)
        elseif (this%forcing_cha) then
            call hook_source(this%decomp, this%mesh, this%fields, this%mix, this%tsim, rhs, this%der, this%dt, this%step, this%dys, this%detady)
        else
            call hook_source(this%decomp, this%mesh, this%fields, this%mix, this%tsim, rhs)
        endif
    end subroutine

    subroutine getRHS_xi(       this,  rhs,&
                           tauxx,tauxy,tauxz,tauyy,tauyz,tauzz,&
                               qx,qy,qz,Jx,Jy,Jz )
        class(cvlgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qx,qy,qz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, this%mix%ns), intent(in) :: Jx,Jy,Jz

        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: flux
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        real(rkind), dimension(:,:,:), pointer :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5,ytmp6
        integer :: i

        xtmp1 => this%xbuf(:,:,:,1); xtmp2 => this%xbuf(:,:,:,2)
        
        ytmp1 => this%ybuf(:,:,:,1); ytmp2 => this%ybuf(:,:,:,2)
        ytmp3 => this%ybuf(:,:,:,3); ytmp4 => this%ybuf(:,:,:,4)
        ytmp5 => this%ybuf(:,:,:,5); ytmp6 => this%ybuf(:,:,:,6)


        ytmp1 = this%dxidx*this%u+ this%dxidy*this%v+ this%dxidz*this%w   !   U^hat
        !ytmp2 is included in species equation
        ytmp3 = this%dxidx*tauxx + this%dxidy*tauxy + this%dxidz*tauxz    ! \tau_xx^hat
        ytmp4 = this%dxidx*tauxy + this%dxidy*tauyy + this%dxidz*tauyz    ! \tau_xy^hat
        ytmp5 = this%dxidx*tauxz + this%dxidy*tauyz + this%dxidz*tauzz    ! \tau_xz^hat
        ytmp6 = this%dxidx*qx    + this%dxidy*qy    + this%dxidz*qz       ! \q_x^hat

        ytmp1 = ytmp1*this%InvJ; ytmp3 = ytmp3*this%InvJ
        ytmp4 = ytmp4*this%InvJ; ytmp5 = ytmp5*this%InvJ; ytmp6 = ytmp6*this%InvJ

        select case(this%mix%ns)
        case(1)
            !flux  = this%Wcnsrv(:,:,:,mom_index  )   ! mass
            flux  =  this%rho * ytmp1                 ! mass
            call transpose_y_to_x(flux,xtmp1,this%decomp)
            call this%der%ddxi(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
            call transpose_x_to_y(xtmp2,flux,this%decomp)
            rhs(:,:,:,1) = rhs(:,:,:,1) - flux
        case default
            do i = 1,this%mix%ns
                !flux = this%Wcnsrv(:,:,:,i)*this%u + Jx(:,:,:,i)   ! mass
                ytmp2 = this%dxidx*Jx(:,:,:,i) + this%dxidy*Jy(:,:,:,i) + this%dxidz*Jz(:,:,:,i)
                ytmp2 = ytmp2*this%InvJ
                flux  = this%Wcnsrv(:,:,:,i)*ytmp1 + ytmp2
                call transpose_y_to_x(flux,xtmp1,this%decomp)
                call this%der%ddxi(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
                call transpose_x_to_y(xtmp2,flux,this%decomp)
                rhs(:,:,:,i) = rhs(:,:,:,i) - flux
            end do
        end select


        !flux = this%Wcnsrv(:,:,:,mom_index  )*this%u + this%p - tauxx ! x-momentum
        flux = this%Wcnsrv(:,:,:,mom_index  )*ytmp1 + this%dxidx*this%InvJ*this%p - ytmp3 ! x-momentum
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddxi(xtmp1,xtmp2, this%x_bc(1), this%x_bc(2)) ! Symmetric for x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux

        !flux = this%Wcnsrv(:,:,:,mom_index  )*this%v          - tauxy ! y-momentum
        flux = this%Wcnsrv(:,:,:,mom_index+1)*ytmp1 + this%dxidy*this%InvJ*this%p - ytmp4 ! y-momentum
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddxi(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux

        !flux = this%Wcnsrv(:,:,:,mom_index  )*this%w          - tauxz ! z-momentum
        flux = this%Wcnsrv(:,:,:,mom_index+2 )*ytmp1 + this%dxidz*this%InvJ*this%p - ytmp5 ! z-momentum
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddxi(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux

        !flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauxx)*this%u - this%v*tauxy - this%w*tauxz + qx ! Total Energy
        flux  = (this%Wcnsrv(:,:,:, TE_index  ) + this%p)*ytmp1 - this%u*ytmp3 - this%v*ytmp4 - this%w*ytmp5 + ytmp6 
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddxi(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux

    end subroutine

    subroutine getRHS_eta(       this,  rhs,&
                           tauxx,tauxy,tauxz,tauyy,tauyz,tauzz,&
                               qx,qy,qz,Jx,Jy,Jz )
        class(cvlgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qx,qy,qz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, this%mix%ns), intent(in) :: Jx,Jy,Jz

        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: flux
        real(rkind), dimension(:,:,:), pointer :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5,ytmp6
        integer :: i

        ytmp1 => this%ybuf(:,:,:,1); ytmp2 => this%ybuf(:,:,:,2)
        ytmp3 => this%ybuf(:,:,:,3); ytmp4 => this%ybuf(:,:,:,4)
        ytmp5 => this%ybuf(:,:,:,5); ytmp6 => this%ybuf(:,:,:,6)


        ytmp1 = this%detadx*this%u+ this%detady*this%v+ this%detadz*this%w
        !ytmp2 is included in species equation and for transpose
        ytmp3 = this%detadx*tauxx + this%detady*tauxy + this%detadz*tauxz
        ytmp4 = this%detadx*tauxy + this%detady*tauyy + this%detadz*tauyz
        ytmp5 = this%detadx*tauxz + this%detady*tauyz + this%detadz*tauzz
        ytmp6 = this%detadx*qx    + this%detady*qy    + this%detadz*qz

        ytmp1 = ytmp1*this%InvJ; ytmp3 = ytmp3*this%InvJ
        ytmp4 = ytmp4*this%InvJ; ytmp5 = ytmp5*this%InvJ; ytmp6 = ytmp6*this%InvJ

        select case(this%mix%ns)
        case(1)
            !flux = this%Wcnsrv(:,:,:,mom_index+1)   ! mass
            flux  =  this%rho * ytmp1                 ! mass
            call this%der%ddeta(flux,ytmp2,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
            rhs(:,:,:,1) = rhs(:,:,:,1) - ytmp2
        case default
            do i = 1,this%mix%ns
                !flux = this%Wcnsrv(:,:,:,i)*this%v + Jy(:,:,:,i) ! mass
                ytmp2 = this%detadx*Jx(:,:,:,i) + this%detady*Jy(:,:,:,i) + this%detadz*Jz(:,:,:,i)
                ytmp2 = ytmp2*this%InvJ
                flux  = this%Wcnsrv(:,:,:,i)*ytmp1 + ytmp2
                call this%der%ddeta(flux,ytmp2,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
                rhs(:,:,:,i) = rhs(:,:,:,i) - ytmp2
            end do
        end select

        !flux = this%Wcnsrv(:,:,:,mom_index+1)*this%u          - tauxy ! x-momentum
        flux = this%Wcnsrv(:,:,:,mom_index  )*ytmp1 + this%detadx*this%InvJ*this%p - ytmp3 ! x-momentum
        call this%der%ddeta(flux,ytmp2,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - ytmp2

        !flux = this%Wcnsrv(:,:,:,mom_index+1)*this%v + this%p - tauyy ! y-momentum
        flux = this%Wcnsrv(:,:,:,mom_index+1)*ytmp1 + this%detady*this%InvJ*this%p - ytmp4 ! y-momentum
        call this%der%ddeta(flux,ytmp2, this%y_bc(1), this%y_bc(2)) ! Symmetric for y-momentum
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - ytmp2

        !flux = this%Wcnsrv(:,:,:,mom_index+1)*this%w          - tauyz ! z-momentum
        flux = this%Wcnsrv(:,:,:,mom_index+2 )*ytmp1 + this%detadz*this%InvJ*this%p - ytmp5 ! z-momentum
        call this%der%ddeta(flux,ytmp2,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - ytmp2

        !flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauyy)*this%v - this%u*tauxy - this%w*tauyz + qy ! Total Energy
        flux  = (this%Wcnsrv(:,:,:, TE_index  ) + this%p)*ytmp1 - this%u*ytmp3 - this%v*ytmp4 - this%w*ytmp5 + ytmp6 
        call this%der%ddeta(flux,ytmp2,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - ytmp2

    end subroutine

    subroutine getRHS_zeta(       this,  rhs,&
                           tauxx,tauxy,tauxz,tauyy,tauyz,tauzz,&
                               qx,qy,qz,Jx,Jy,Jz )
        class(cvlgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qx,qy,qz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, this%mix%ns), intent(in) :: Jx,Jy,Jz

        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: flux
        real(rkind), dimension(:,:,:), pointer :: ztmp1,ztmp2
        real(rkind), dimension(:,:,:), pointer :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5,ytmp6
        integer :: i

        ztmp1 => this%zbuf(:,:,:,1); ztmp2 => this%zbuf(:,:,:,2)
        
        ytmp1 => this%ybuf(:,:,:,1); ytmp2 => this%ybuf(:,:,:,2)
        ytmp3 => this%ybuf(:,:,:,3); ytmp4 => this%ybuf(:,:,:,4)
        ytmp5 => this%ybuf(:,:,:,5); ytmp6 => this%ybuf(:,:,:,6)


        ytmp1 = this%dzetadx*this%u+ this%dzetady*this%v+ this%dzetadz*this%w
        !ytmp2 is included in species equation
        ytmp3 = this%dzetadx*tauxx + this%dzetady*tauxy + this%dzetadz*tauxz
        ytmp4 = this%dzetadx*tauxy + this%dzetady*tauyy + this%dzetadz*tauyz
        ytmp5 = this%dzetadx*tauxz + this%dzetady*tauyz + this%dzetadz*tauzz
        ytmp6 = this%dzetadx*qx    + this%dzetady*qy    + this%dzetadz*qz

        ytmp1 = ytmp1*this%InvJ; ytmp3 = ytmp3*this%InvJ
        ytmp4 = ytmp4*this%InvJ; ytmp5 = ytmp5*this%InvJ; ytmp6 = ytmp6*this%InvJ
        
        select case(this%mix%ns)
        case(1)
            !flux = this%Wcnsrv(:,:,:,mom_index+2)   ! mass
            flux  =  this%rho * ytmp1                 ! mass
            call transpose_y_to_z(flux,ztmp1,this%decomp)
            call this%der%ddzeta(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
            call transpose_z_to_y(ztmp2,flux,this%decomp)
            rhs(:,:,:,1) = rhs(:,:,:,1) - flux
        case default
            do i = 1,this%mix%ns
                !flux = this%Wcnsrv(:,:,:,i)*this%w + Jz(:,:,:,i)   ! mass
                ytmp2 = this%dzetadx*Jx(:,:,:,i) + this%dzetady*Jy(:,:,:,i) + this%dzetadz*Jz(:,:,:,i)
                ytmp2 = ytmp2*this%InvJ
                flux  = this%Wcnsrv(:,:,:,i)*ytmp1 + ytmp2
                call transpose_y_to_z(flux,ztmp1,this%decomp)
                call this%der%ddzeta(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
                call transpose_z_to_y(ztmp2,flux,this%decomp)
                rhs(:,:,:,i) = rhs(:,:,:,i) - flux
            end do
        end select

        !flux = this%Wcnsrv(:,:,:,mom_index+2)*this%u          - tauxz ! x-momentum
        flux = this%Wcnsrv(:,:,:,mom_index   )*ytmp1 + this%dzetadx*this%InvJ*this%p - ytmp3 ! x-momentum
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddzeta(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux

        !flux = this%Wcnsrv(:,:,:,mom_index+2)*this%v          - tauyz ! y-momentum
        flux = this%Wcnsrv(:,:,:,mom_index+1 )*ytmp1 + this%dzetady*this%InvJ*this%p - ytmp4 ! y-momentum
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddzeta(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux

        !flux = this%Wcnsrv(:,:,:,mom_index+2)*this%w + this%p - tauzz ! z-momentum
        flux = this%Wcnsrv(:,:,:,mom_index+2 )*ytmp1 + this%dzetadz*this%InvJ*this%p - ytmp5 ! z-momentum
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddzeta(ztmp1,ztmp2, this%z_bc(1), this%z_bc(2)) ! Symmetric for z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux

        !flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauzz)*this%w - this%u*tauxz - this%v*tauyz + qz ! Total Energy
        flux  = (this%Wcnsrv(:,:,:, TE_index  ) + this%p)*ytmp1 - this%u*ytmp3 - this%v*ytmp4 - this%w*ytmp5 + ytmp6 
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddzeta(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux

    end subroutine

    subroutine getLAD(this, dudx, dudy, dudz,&
                            dvdx, dvdy, dvdz,&
                            dwdx, dwdy, dwdz,&
                           dYsdx,dYsdy,dYsdz )
        use reductions, only: P_MAXVAL
        class(cvlgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, this%mix%ns), optional, intent(in) :: dYsdx,dYsdy,dYsdz
        
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: mustar,bulkstar,kapstar,diffstar,func
        integer :: i, j, k, p
        
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        real(rkind), dimension(:,:,:), pointer :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5,ytmp6
        real(rkind), dimension(:,:,:), pointer :: ztmp1,ztmp2

        xtmp1 => this%xbuf(:,:,:,1); xtmp2 => this%xbuf(:,:,:,2)
        
        ytmp1 => this%ybuf(:,:,:,1); ytmp2 => this%ybuf(:,:,:,2)
        ytmp3 => this%ybuf(:,:,:,3); ytmp4 => this%ybuf(:,:,:,4)
        ytmp5 => this%ybuf(:,:,:,5); ytmp6 => this%ybuf(:,:,:,6)
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
        !xtmp2(i,j,k) = xtmp1(i,j,k)*this%dxs(i)**6
        !call transpose_x_to_y(xtmp2,mustar,this%decomp)
        !!! To include dxs
        call transpose_x_to_y(xtmp1,ytmp1,this%decomp)
        do k = 1, this%nzp
           do j = 1, this%nyp
              do i = 1, this%nxp
                 mustar(i,j,k) = ytmp1(i,j,k)*this%dxs(i,j,k)**6
              end do
           end do
        end do
        
        ! Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
        !ztmp2(:,:,k) = ztmp1(:,:,k)*this%dzs(k)**6
        !call transpose_z_to_y(ztmp2,ytmp1,this%decomp)
        !!! To include dzs
        call transpose_z_to_y(ztmp1,ytmp1,this%decomp)
        do k = 1, this%nzp
           ytmp2(:,:,k) = ytmp1(:,:,k)*this%dzs(:,:,k)**6
        end do
        mustar = mustar + ytmp2
        
        ! Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp1,this%y_bc(1),this%y_bc(2))
        call this%der%d2dy2(ytmp1,ytmp2,this%y_bc(1),this%y_bc(2))
        !ytmp1 = ytmp2*this%dy**6
        !!! To include dys
        do k = 1, this%nzp
           do j = 1, this%nyp
              ytmp1(:,j,k) = ytmp2(:,j,k)*this%dys(:,j,k)**6
           end do
        end do
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
        !xtmp2 = xtmp1*this%dx**4
        !call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        !bulkstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2
        !!! To include dxs
        call transpose_x_to_y(xtmp1,ytmp4,this%decomp)
        do k = 1, this%nzp
           do j = 1, this%nyp
              do i = 1, this%nxp
                 ytmp5(i,j,k)    = ytmp4(i,j,k) * this%dxs(i,j,k)**4
                 bulkstar(i,j,k) = ytmp5(i,j,k) * ( this%dxs(i,j,k) * ytmp1(i,j,k) / (ytmp1(i,j,k) + &
                                      ytmp2(i,j,k) + ytmp3(i,j,k) + real(1.0D-32,rkind)) )**2
              end do
           end do
        end do

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
        !ztmp2 = ztmp1*this%dz**4
        !call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        !bulkstar = bulkstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2
        !!! To include dzs
        call transpose_z_to_y(ztmp1,ytmp4,this%decomp)
        do k = 1, this%nzp
           ytmp5(:,:,k) = ytmp4(:,:,k)*this%dzs(:,:,k)**4
           bulkstar(:,:,k) = bulkstar(:,:,k) + ytmp5(:,:,k) * ( this%dzs(:,:,k) * ytmp3(:,:,k) / (ytmp1(:,:,k) + &
                               ytmp2(:,:,k) + ytmp3(:,:,k) + real(1.0D-32,rkind)) )**2
        end do

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp4,this%y_bc(1),this%y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,this%y_bc(1),this%y_bc(2))
        !ytmp4 = ytmp5*this%dy**4
        !bulkstar = bulkstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2
        !!! To include dys
        do k = 1, this%nzp
           do j = 1, this%nyp
              ytmp4(:,j,k)    = ytmp5(:,j,k)*this%dys(:,j,k)**4
              bulkstar(:,j,k) = bulkstar(:,j,k) + ytmp4(:,j,k) * ( this%dys(:,j,k) * ytmp2(:,j,k) / (ytmp1(:,j,k) &
                                  + ytmp2(:,j,k) + ytmp3(:,j,k) + real(1.0D-32,rkind)) )**2
           end do
        end do

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
        !xtmp2 = xtmp1*this%dx**4
        !call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        !kapstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
        !!! To include dxs
        call transpose_x_to_y(xtmp1,ytmp4,this%decomp)
        do k = 1, this%nzp
           do j = 1, this%nyp
              do i = 1, this%nxp
                 ytmp5(i,j,k) = ytmp4(i,j,k)*this%dxs(i,j,k)**4
                 kapstar(i,j,k) = ytmp5(i,j,k) * ( this%dxs(i,j,k) * ytmp1(i,j,k) / (ytmp1(i,j,k) &
                                    + ytmp2(i,j,k) + ytmp3(i,j,k) + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
              end do
           end do
        end do

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(this%e,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
        !ztmp2 = ztmp1*this%dz**4
        !call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        !kapstar = kapstar + ytmp4 * ( this%dz * ytmp3/(ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)))! Add eps in case denominator is zero
        !!! To include dzs 
        call transpose_z_to_y(ztmp1,ytmp4,this%decomp)
        do k = 1, this%nzp
           ytmp5(:,:,k) = ytmp4(:,:,k)*this%dzs(:,:,k)**4
           kapstar(:,:,k) = kapstar(:,:,k) + ytmp5(:,:,k) * ( this%dzs(:,:,k) * ytmp3(:,:,k) / (ytmp1(:,:,k) &
                              + ytmp2(:,:,k) + ytmp3(:,:,k) + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
        end do

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(this%e,ytmp4,this%y_bc(1),this%y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,this%y_bc(1),this%y_bc(2))
        !ytmp4 = ytmp5*this%dy**4
        !kapstar = kapstar + ytmp4 * (this%dy * ytmp2/(ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind))) ! Add eps in case denominator is zero
        !!! To include dys
        do k = 1, this%nzp
           do j = 1, this%nyp
              ytmp4(:,j,k) = ytmp5(:,j,k)*this%dys(:,j,k)**4
              kapstar(:,j,k) = kapstar(:,j,k) + ytmp4(:,j,k) * ( this%dys(:,j,k) * ytmp2(:,j,k) / &
                          (ytmp1(:,j,k) + ytmp2(:,j,k) + ytmp3(:,j,k) + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
           end do
        end do

        ! Now, all ytmps are free to use
        call this%mix%get_sos(this%rho,this%p,ytmp1)  ! Speed of sound

        kapstar = this%Ckap*this%rho*ytmp1*abs(kapstar)/this%T

        ! Filter kapstar
        call this%filter(kapstar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)

        ! -------- Artificial Diffusivity ---------
        if (present(dYsdx) .AND. present(dYsdy) .AND. present(dYsdz)) then
            if (this%mix%ns .GT. 1) then
                do p = 1,this%mix%ns
                    ! Step 1: Get components of grad(Ys) squared individually
                    ytmp1 = dYsdx(:,:,:,p)*dYsdx(:,:,:,p)
                    ytmp2 = dYsdy(:,:,:,p)*dYsdy(:,:,:,p)
                    ytmp3 = dYsdz(:,:,:,p)*dYsdz(:,:,:,p)

                    call this%mix%get_sos(this%rho,this%p,ytmp6)  ! Speed of sound  

                    ! Step 2: Get 4th derivative in X
                    call transpose_y_to_x(this%Ys(:,:,:,p),xtmp1,this%decomp)
                    call this%der%d2dx2(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
                    call this%der%d2dx2(xtmp2,xtmp1,this%x_bc(1),this%x_bc(2))
                    !xtmp2 = xtmp1*this%dx**4
                    !call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
                    !diffstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
                    !!! To include dxs
                    call transpose_x_to_y(xtmp1,ytmp4,this%decomp)
                    do k = 1, this%nzp
                       do j = 1, this%nyp
                          do i = 1, this%nxp
                             ytmp5(i,j,k) = ytmp4(i,j,k)*this%dxs(i,j,k)**4
                             diffstar(i,j,k) = ytmp5(i,j,k) * ( this%dxs(i,j,k) * ytmp1(i,j,k) / (ytmp1(i,j,k) + &
                                           ytmp2(i,j,k) + ytmp3(i,j,k) + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
                          end do
                       end do
                    end do

                    ! Step 3: Get 4th derivative in Z
                    call transpose_y_to_z(this%Ys(:,:,:,p),ztmp1,this%decomp)
                    call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
                    call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
                    !ztmp2 = ztmp1*this%dz**4
                    !call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
                    !diffstar = diffstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
                    !!! To include dzs
                    call transpose_z_to_y(ztmp1,ytmp4,this%decomp)
                    do k = 1, this%nzp
                       ytmp5(:,:,k) = ytmp4(:,:,k)*this%dzs(:,:,k)**4
                       diffstar(:,:,k) = diffstar(:,:,k) + ytmp5(:,:,k) * ( this%dzs(:,:,k) * ytmp3(:,:,k) / (ytmp1(:,:,k) &
                                  + ytmp2(:,:,k) + ytmp3(:,:,k) + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
                    end do

                    ! Step 4: Get 4th derivative in Y
                    call this%der%d2dy2(this%Ys(:,:,:,p),ytmp4,this%y_bc(1),this%y_bc(2))
                    call this%der%d2dy2(ytmp4,ytmp5,this%y_bc(1),this%y_bc(2))
                    !ytmp4 = ytmp5*this%dy**4
                    !diffstar = diffstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
                    !!! To include dys
                    do k = 1, this%nzp
                       do j = 1, this%nyp
                          ytmp4(:,j,k)    = ytmp5(:,j,k)*this%dys(:,j,k)**4
                          diffstar(:,j,k) = diffstar(:,j,k) + ytmp4(:,j,k) * ( this%dys(:,j,k) * ytmp2(:,j,k) / (ytmp1(:,j,k) + &
                                           ytmp2(:,j,k) + ytmp3(:,j,k) + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
                       end do
                    end do

                    diffstar = this%Cdiff*ytmp6*abs(diffstar)     ! CD part of diff
                    !ytmp4 = (sqrt(ytmp1)*this%dx + sqrt(ytmp2)*this%dy + sqrt(ytmp3)*this%dz) &
                    !               / ( sqrt(ytmp1+ytmp2+ytmp3) + real(1.0D-32,rkind) ) ! grid scale 
                    !!! To include dxs, dys, dzs
                    do k = 1, this%nzp
                       do j = 1, this%nyp
                          do i = 1, this%nxp
                              ytmp4(i,j,k) = (sqrt(ytmp1(i,j,k))*this%dxs(i,j,k) + sqrt(ytmp2(i,j,k))*this%dys(i,j,k) + &
                                  sqrt(ytmp3(i,j,k))*this%dzs(i,j,k)) / ( sqrt(ytmp1(i,j,k)+ytmp2(i,j,k)+ytmp3(i,j,k)) + real(1.0D-32,rkind) )
                          end do
                       end do
                    end do

                    ytmp5 = this%CY*ytmp6*( half*(abs(this%Ys(:,:,:,p))-one + abs(this%Ys(:,:,:,p)-one)) )*ytmp4 ! CY part of diff

                    diffstar = max(diffstar, ytmp5) ! Take max of both terms instead of add to minimize the dissipation

                    ! Filter diffstar
                    call this%filter(diffstar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)

                    ! Add to physical diffusivity
                    this%diff(:,:,:,p) = this%diff(:,:,:,p) + diffstar
                end do
            end if
        end if

        ! Now, add to physical fluid properties
        this%mu   = this%mu   + mustar
        this%bulk = this%bulk + bulkstar
        this%kap  = this%kap  + kapstar

    end subroutine

    subroutine filter(this,arr,myfil,numtimes, x_bc, y_bc, z_bc)
        class(cvlgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(inout) :: arr
        type(filters), target, optional, intent(in) :: myfil
        integer, optional, intent(in) :: numtimes
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc
        
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
        class(cvlgrid), intent(inout) :: this

        ! TODO
        ! Hard code these values for now. Need to make a better interface for this later
        real(rkind) :: mu_ref = 0.6_rkind * half / (100._rkind * sqrt(3._rkind))
        real(rkind) :: T_ref = 1._rkind / 1.4_rkind
        real(rkind) :: Pr = 0.70_rkind

        this%mu   = mu_ref * (this%T / T_ref)**(three/four)
        this%bulk = zero
        this%kap  = this%mix%material(1)%mat%gam / (this%mix%material(1)%mat%gam - one) * this%mix%material(1)%mat%Rgas * this%mu / Pr
        this%diff = zero

        ! If inviscid set everything to zero (otherwise use a model)
        ! this%mu = zero
        ! this%bulk = zero
        ! this%kap = zero
        ! this%diff = zero

    end subroutine  

    subroutine get_tau(this,duidxj)
        class(cvlgrid), target, intent(inout) :: this
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

    subroutine get_q(this,duidxj,Jx,Jy,Jz)
        use exits, only: nancheck
        class(cvlgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(inout) :: duidxj
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,this%mix%ns), intent(in) :: Jx,Jy,Jz

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

        ! If multispecies, add the inter-species enthalpy flux
        if (this%mix%ns .GT. 1) then
            do i = 1,this%mix%ns
                call this%mix%material(i)%mat%get_enthalpy(this%T,tmp1_in_y)
                duidxj(:,:,:,qxidx) = duidxj(:,:,:,qxidx) + ( tmp1_in_y * Jx(:,:,:,i) )
                duidxj(:,:,:,qyidx) = duidxj(:,:,:,qyidx) + ( tmp1_in_y * Jy(:,:,:,i) )
                duidxj(:,:,:,qzidx) = duidxj(:,:,:,qzidx) + ( tmp1_in_y * Jz(:,:,:,i) )
            end do
        end if

        ! Done
    end subroutine 

    subroutine get_J(this,gradYs)
        class(cvlgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,3*this%mix%ns), target, intent(in) :: gradYs
        real(rkind), dimension(:,:,:,:), pointer :: dYsdx, dYsdy, dYsdz
        real(rkind), dimension(:,:,:), pointer :: sumJx, sumJy, sumJz
        integer :: i

        sumJx => this%ybuf(:,:,:,1)
        sumJy => this%ybuf(:,:,:,2)
        sumJz => this%ybuf(:,:,:,3)

        dYsdx => gradYs(:,:,:,1:this%mix%ns); dYsdy => gradYs(:,:,:,this%mix%ns+1:2*this%mix%ns); dYsdz => gradYs(:,:,:,2*this%mix%ns+1:3*this%mix%ns);

        ! Get diff*gradYs
        do i=1,this%mix%ns
            dYsdx(:,:,:,i) = dYsdx(:,:,:,i)*this%diff(:,:,:,i)
            dYsdy(:,:,:,i) = dYsdy(:,:,:,i)*this%diff(:,:,:,i)
            dYsdz(:,:,:,i) = dYsdz(:,:,:,i)*this%diff(:,:,:,i)
        end do

        ! Get sum of diff*gradYs and correct so this sum becomes zero (No net diffusive flux)
        sumJx = sum(dYsdx,4); sumJy = sum(dYsdy,4); sumJz = sum(dYsdz,4);

        ! Put the fluxes in dYsdx itself
        do i=1,this%mix%ns
            dYsdx(:,:,:,i) = -this%rho*( dYsdx(:,:,:,i) - this%Ys(:,:,:,i)*sumJx )
            dYsdy(:,:,:,i) = -this%rho*( dYsdy(:,:,:,i) - this%Ys(:,:,:,i)*sumJy )
            dYsdz(:,:,:,i) = -this%rho*( dYsdz(:,:,:,i) - this%Ys(:,:,:,i)*sumJz )
        end do

    end subroutine

    subroutine write_viz(this)
        use exits, only: message
        use timer, only: tic, toc
        class(cvlgrid), intent(inout) :: this
        character(len=clen) :: charout
        real(rkind) :: cputime
        integer :: i

        call tic()

        ! Start visualization dump
        call this%viz%start_viz(this%tsim)

        write(charout,'(A,I0,A,A)') "Writing visualization dump ", this%viz%vizcount, " to ", adjustl(trim(this%viz%filename))
        call message(charout)

        ! Write variables
        call this%viz%write_variable(this%rho , 'rho' )
        call this%viz%write_variable(this%u   , 'u'   )
        call this%viz%write_variable(this%v   , 'v'   )
        call this%viz%write_variable(this%w   , 'w'   )
        call this%viz%write_variable(this%p   , 'p'   )
        call this%viz%write_variable(this%T   , 'T'   )
        call this%viz%write_variable(this%e   , 'e'   )
        call this%viz%write_variable(this%mu  , 'mu'  )
        call this%viz%write_variable(this%bulk, 'bulk')
        call this%viz%write_variable(this%kap , 'kap' )

        if (this%mix%ns > 1) then
            do i = 1,this%mix%ns
                write(charout,'(I2.2)') i
                call this%viz%write_variable(this%Ys(:,:,:,i)  , 'Massfraction_'//adjustl(trim(charout)) )
                call this%viz%write_variable(this%diff(:,:,:,i), 'Diffusivity_'//adjustl(trim(charout)) )
            end do
        end if
        
        if (this%useSGS ) then
           call this%viz%write_variable(this%tausgs(:,:,:,1) ,'t11')
           call this%viz%write_variable(this%tausgs(:,:,:,2) ,'t12')
           call this%viz%write_variable(this%tausgs(:,:,:,3) ,'t13')
           call this%viz%write_variable(this%tausgs(:,:,:,4) ,'t22')
           call this%viz%write_variable(this%tausgs(:,:,:,5) ,'t23')
           call this%viz%write_variable(this%tausgs(:,:,:,6) ,'t33')
           call this%viz%write_variable(this%Qjsgs(:,:,:,1) ,'q1')
           call this%viz%write_variable(this%Qjsgs(:,:,:,2) ,'q2')
           call this%viz%write_variable(this%Qjsgs(:,:,:,3) ,'q3')
          ! call this%viz%write_variable(this%sgsmodel%nusgs,'nusgs')
        end if
 
        ! TODO: Add hook_viz here

        ! End visualization dump
        call this%viz%end_viz()

        call toc(cputime)
        write(charout,'(A,ES11.3,A)') "Finished writing visualization dump in ", cputime, " seconds"
        call message(charout)

    end subroutine

    subroutine write_restart(this)
        use exits, only: message
        use timer, only: tic, toc
        class(cvlgrid), intent(inout) :: this
        character(len=clen) :: charout
        real(rkind) :: cputime
        integer :: i

        ! Get updated conserved variables
        call this%get_conserved()

        call tic()

        ! Start visualization dump
        call this%restart%start_viz(this%tsim)

        write(charout,'(A,I0,A,A)') "Writing restart dump ", this%restart%vizcount, " to ", adjustl(trim(this%restart%filename))
        call message(charout)

        ! Write conserved variables
        if (this%mix%ns > 1) then
            do i = 1,this%mix%ns
                write(charout,'(I4.4)') i
                call this%restart%write_variable(this%Wcnsrv(:,:,:,i), 'rhoY_'//adjustl(trim(charout)) )
            end do
        else
            call this%restart%write_variable(this%Wcnsrv(:,:,:,mass_index), 'rho' )
        end if
        call this%restart%write_variable(this%Wcnsrv(:,:,:,mom_index  ), 'rhou')
        call this%restart%write_variable(this%Wcnsrv(:,:,:,mom_index+1), 'rhov')
        call this%restart%write_variable(this%Wcnsrv(:,:,:,mom_index+2), 'rhow')
        call this%restart%write_variable(this%Wcnsrv(:,:,:, TE_index  ), 'TE')

        call this%restart%write_attribute(1, [this%step], 'step', '/')

        ! End visualization dump
        call this%restart%end_viz()

        call toc(cputime)
        write(charout,'(A,ES11.3,A)') "Finished writing restart dump in ", cputime, " seconds"
        call message(charout)

    end subroutine

    subroutine read_restart(this, vizcount)
        use exits, only: message
        use timer, only: tic, toc
        class(cvlgrid), intent(inout) :: this
        integer,      intent(in)    :: vizcount
        character(len=clen) :: charout
        real(rkind) :: cputime
        integer :: i
        integer, dimension(1) :: tmp_int
        real(rkind), dimension(1) :: tmp_real

        ! Start reading restart file
        call this%restart%start_reading(vizcount)

        call tic()

        write(charout,'(A,A)') "Reading restart dump from ", adjustl(trim(this%restart%filename))
        call message(charout)

        ! Read conserved variables
        if (this%mix%ns > 1) then
            do i = 1,this%mix%ns
                write(charout,'(I4.4)') i
                call this%restart%read_dataset(this%Wcnsrv(:,:,:,i), 'rhoY_'//adjustl(trim(charout)) )
            end do
        else
            call this%restart%read_dataset(this%Wcnsrv(:,:,:,mass_index), 'rho' )
        end if
        call this%restart%read_dataset(this%Wcnsrv(:,:,:,mom_index  ), 'rhou')
        call this%restart%read_dataset(this%Wcnsrv(:,:,:,mom_index+1), 'rhov')
        call this%restart%read_dataset(this%Wcnsrv(:,:,:,mom_index+2), 'rhow')
        call this%restart%read_dataset(this%Wcnsrv(:,:,:, TE_index  ), 'TE')

        call this%restart%read_attribute(1, tmp_int, 'step')
        this%step = tmp_int(1)

        ! Reset nsteps so that simulation runs for no. of steps in inputfile
        if (this%nsteps > 0) this%nsteps = this%step + this%nsteps

        call this%restart%read_attribute(1, tmp_real, 'Time')
        this%tsim = tmp_real(1)

        ! End visualization dump
        call this%restart%end_reading()

        call toc(cputime)
        write(charout,'(A,ES11.3,A)') "Finished reading restart dump in ", cputime, " seconds"
        call message(charout)

        ! Get primitive variables from conserved
        call this%get_primitive()

    end subroutine
    
    subroutine setup_postprocessing(this, nrestarts)
        use mpi
        class(cvlgrid), intent(inout) :: this
        integer,      intent(out)   :: nrestarts

        ! Destroy old restart object
        call this%restart%destroy()

        ! Initialize new restart object in read only mode
        call this%restart%init( mpi_comm_world, this%decomp, 'y', this%outputdir, 'restart', &
                                reduce_precision=.false., write_xdmf=.false., read_only=.true., jump_to_last=.true.)
        nrestarts = this%restart%vizcount

    end subroutine

end module 
