module CompressibleGrid
    use kind_parameters,      only: rkind, clen
    use constants,            only: zero,eps,third,half,one,two,three,four
    use FiltersMod,           only: filters
    use GridMod,              only: grid
    use gridtools,            only: alloc_buffs, destroy_buffs
    use cgrid_hooks,          only: meshgen, initfields, hook_output, hook_bc, hook_timestep, hook_source
    use decomp_2d,            only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,       only: derivatives
    use io_hdf5_stuff,        only: io_hdf5
    use IdealGasEOS,          only: idealgas
    use MixtureEOSMod,        only: mixture
    use PowerLawViscosityMod, only: powerLawViscosity
   
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

    type, extends(grid) :: cgrid
       
        type(filters), allocatable :: gfil
        type(mixture), allocatable :: mix

        type(io_hdf5), allocatable :: viz
        type(io_hdf5), allocatable :: restart

        real(rkind), dimension(:,:,:,:), allocatable :: Wcnsrv                               ! Conserved variables
        real(rkind), dimension(:,:,:,:), allocatable :: xbuf, ybuf, zbuf   ! Buffers
       
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
            procedure, private :: get_J
            procedure, private :: write_viz
            procedure, private :: write_restart
            procedure, private :: read_restart
    end type

contains
    subroutine init(this, inputfile )
        use mpi
        use reductions, only: P_MAXVAL, P_MINVAL
        use exits, only: message, nancheck, GracefulExit
        class(cgrid),target, intent(inout) :: this
        character(len=clen), intent(in) :: inputfile  

        integer :: nx, ny, nz
        integer :: ns = 1
        character(len=clen) :: outputdir
        character(len=clen) :: inputdir
        character(len=clen) :: vizprefix = "cgrid"
        logical :: reduce_precision = .true.
        real(rkind) :: tviz = zero
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

        namelist /INPUT/ nx, ny, nz, tstop, dt, CFL, nsteps, inputdir, &
                         outputdir, vizprefix, tviz, reduce_precision, &
                                      periodicx, periodicy, periodicz, &
                             derivative_x, derivative_y, derivative_z, &
                                         filter_x, filter_y, filter_z, &
                                                           prow, pcol, &
                                                             SkewSymm  
        namelist /CINPUT/  ns, gam, Rgas, Cmu, Cbeta, Ckap, Cdiff, CY, &
                           inviscid, nrestart, rewrite_viz, vizramp,   &
                           x_bc1, x_bcn, y_bc1, y_bcn, z_bc1, z_bcn


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

        ! Go to hooks if a different mesh is desired 
        call meshgen(this%decomp, this%dx, this%dy, this%dz, this%mesh) 

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

        ! Initialize everything to a constant Zero
        this%fields = zero
        this%Ys(:,:,:,1) = one   ! So that all massfractions add up to unity

        ! Finally, set the local array dimensions
        this%nxp = this%decomp%ysz(1)
        this%nyp = this%decomp%ysz(2)
        this%nzp = this%decomp%ysz(3)

        ! Go to hooks if a different initialization is derired 
        call initfields(this%decomp, this%dx, this%dy, this%dz, inputfile, this%mesh, this%fields, &
                        this%mix, this%tsim, this%tstop, this%dtfixed, tviz)
        
        ! Check for correct initialization of the mixture object
        call this%mix%check_initialization()

        ! Update mix
        call this%mix%update(this%Ys)
        call this%mix%get_e_from_p(this%rho,this%p,this%e)
        call this%mix%get_T(this%e,this%T)

        ! print *, "Cp(1,1,1) = ", this%mix%Cp(1,1,1)
        ! print *, "Cp(-1,1,1) = ", this%mix%Cp(this%nxp,1,1)

        ! print *, "T(1,1,1) = ", this%T(1,1,1)
        ! print *, "T(-1,1,1) = ", this%T(this%nxp,1,1)

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


        ! Allocate 2 buffers for each of the three decompositions
        call alloc_buffs(this%xbuf,nbufsx,"x",this%decomp)
        call alloc_buffs(this%ybuf,nbufsy,"y",this%decomp)
        call alloc_buffs(this%zbuf,nbufsz,"z",this%decomp)

        this%SkewSymm = SkewSymm

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
        if (this%restart%vizcount > 0) call this%read_restart(this%restart%vizcount)
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

    end subroutine


    subroutine destroy(this)
        class(cgrid), intent(inout) :: this

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

        if (allocated(this%mix)) deallocate(this%mix)
        
        if (allocated(this%Wcnsrv)) deallocate(this%Wcnsrv) 
        
        call this%viz%destroy()
        if (allocated(this%viz)) deallocate(this%viz)

        call this%restart%destroy()
        if (allocated(this%restart)) deallocate(this%restart)

        call decomp_2d_finalize
        if (allocated(this%decomp)) deallocate(this%decomp) 

    end subroutine

    subroutine gradient(this, f, dfdx, dfdy, dfdz, x_bc, y_bc, z_bc)
        class(cgrid),target, intent(inout) :: this
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
        class(cgrid),target, intent(inout) :: this
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
        use reductions, only: P_MAXVAL, P_MINVAL
        class(cgrid), target, intent(inout) :: this

        logical :: tcond, vizcond, stepcond
        character(len=clen) :: stability
        real(rkind) :: cputime
        real(rkind), dimension(:,:,:,:), allocatable, target :: duidxj
        real(rkind), dimension(:,:,:,:), allocatable, target :: gradYs
        real(rkind), dimension(:,:,:),   pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:,:), pointer :: dYsdx, dYsdy, dYsdz
        integer :: i
        integer :: stepramp = 0

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
        deallocate( duidxj )
        ! ------------------------------------------------

        ! Write out initial conditions
        if (this%viz%vizcount == 0) then 
            call hook_output(this%decomp, this%der, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
            call this%write_viz()
        end if
        vizcond = .FALSE.
       
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
            call this%advance_RK45()
            call toc(cputime)
            
            call message(1,"Time",this%tsim)
            call message(2,"Time step",this%dt)
            call message(2,"Stability limit: "//trim(stability))
            call message(2,"CPU time (in seconds)",cputime)
            call hook_timestep(this%decomp, this%mesh, this%fields, this%mix, this%step, this%tsim)
          
            ! Write out vizualization dump if vizcond is met 
            if (vizcond) then
                call hook_output(this%decomp, this%der, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
                ! call this%viz%WriteViz(this%decomp, this%mesh, this%fields, this%tsim)
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

    end subroutine

    subroutine advance_RK45(this)
        use RKCoeffs,   only: RK45_steps,RK45_A,RK45_B
        use exits,      only: message,nancheck,GracefulExit
        class(cgrid), target, intent(inout) :: this

        real(rkind)                                          :: Qtmpt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv) :: rhs  ! RHS for conserved variables
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv) :: Qtmp ! Temporary variable for RK45
        integer :: isub,i,j,k,l

        character(len=clen) :: charout

        if (this%step == 0) then
            call this%mix%update(this%Ys)
            call this%mix%get_e_from_p(this%rho,this%p,this%e)
            call this%mix%get_T(this%e,this%T)
        end if

        Qtmp = zero
        Qtmpt = zero

        do isub = 1,RK45_steps
            call this%get_conserved()

            if ( nancheck(this%Wcnsrv,i,j,k,l) ) then
                call message("Wcnsrv: ",this%Wcnsrv(i,j,k,l))
                write(charout,'(A,I0,A,I0,A,4(I0,A))') "NaN encountered in solution at substep ", isub, " of step ", this%step+1, " at Wcnsrv(",i,",",j,",",k,",",l,")"
                print *, charout
                call GracefulExit(trim(charout), 999)
            end if

            call this%getRHS(rhs)
            Qtmp = this%dt*rhs + RK45_A(isub)*Qtmp
            Qtmpt = this%dt + RK45_A(isub)*Qtmpt
            this%Wcnsrv = this%Wcnsrv + RK45_B(isub)*Qtmp
            this%tsim = this%tsim + RK45_B(isub)*Qtmpt

            ! Filter the conserved variables
            do i = 1,this%mix%ns
                call this%filter(this%Wcnsrv(:,:,:,i), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
            end do
            call this%filter(this%Wcnsrv(:,:,:,mom_index  ), this%fil, 1,-this%x_bc, this%y_bc, this%z_bc)
            call this%filter(this%Wcnsrv(:,:,:,mom_index+1), this%fil, 1, this%x_bc,-this%y_bc, this%z_bc)
            call this%filter(this%Wcnsrv(:,:,:,mom_index+2), this%fil, 1, this%x_bc, this%y_bc,-this%z_bc)
            call this%filter(this%Wcnsrv(:,:,:, TE_index  ), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
            
            call this%get_primitive()
            call hook_bc(this%decomp, this%mesh, this%fields, this%mix, this%tsim, this%x_bc, this%y_bc, this%z_bc)
            call this%post_bc()
        end do

        !this%tsim = this%tsim + this%dt
        this%step = this%step + 1
            
    end subroutine

    subroutine get_dt(this,stability)
        use reductions, only : P_MAXVAL
        class(cgrid), target, intent(inout) :: this
        character(len=*), intent(out) :: stability
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: cs
        real(rkind) :: dtCFL, dtmu, dtbulk, dtkap, dtdiff

        call this%mix%get_sos(this%rho,this%p,cs)  ! Speed of sound - hydrodynamic part

        dtCFL  = this%CFL / P_MAXVAL( ABS(this%u)/this%dx + ABS(this%v)/this%dy + ABS(this%w)/this%dz &
               + cs*sqrt( one/(this%dx**two) + one/(this%dy**two) + one/(this%dz**two) ))
        dtmu   = 0.2_rkind * min(this%dx,this%dy,this%dz)**2 / (P_MAXVAL( this%mu  / this%rho ) + eps)
        dtbulk = 0.2_rkind * min(this%dx,this%dy,this%dz)**2 / (P_MAXVAL( this%bulk/ this%rho ) + eps)
        dtkap  = 0.2_rkind * one / ( (P_MAXVAL( this%kap*this%T/(this%rho* (min(this%dx,this%dy,this%dz)**4))))**(third) + eps)
        dtdiff = 0.2_rkind * min(this%dx,this%dy,this%dz)**2 / (P_MAXVAL( this%diff ) + eps)

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
        class(cgrid), target, intent(inout) :: this
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
        class(cgrid), intent(inout) :: this
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
        class(cgrid), intent(inout) :: this

        call this%mix%update(this%Ys)
        call this%mix%get_e_from_p(this%rho,this%p,this%e)
        call this%mix%get_T(this%e,this%T)

    end subroutine

    subroutine getRHS(this, rhs)
        class(cgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,ncnsrv), intent(out) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,3*this%mix%ns), target :: gradYs
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:,:), pointer :: dYsdx, dYsdy, dYsdz
        real(rkind), dimension(:,:,:), pointer :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        real(rkind), dimension(:,:,:), pointer :: qx,qy,qz
        real(rkind), dimension(:,:,:,:), pointer :: Jx,Jy,Jz
        integer :: i

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


        rhs = zero
        call this%getRHS_x(              rhs,&
                           tauxx,tauxy,tauxz,&
                               qx, Jx )

        call this%getRHS_y(              rhs,&
                           tauxy,tauyy,tauyz,&
                               qy, Jy )

        call this%getRHS_z(              rhs,&
                           tauxz,tauyz,tauzz,&
                               qz, Jz )

        ! Call problem source hook
        call hook_source(this%decomp, this%mesh, this%fields, this%mix, this%tsim, rhs)

    end subroutine

    subroutine getRHS_x(       this,  rhs,&
                        tauxx,tauxy,tauxz,&
                            qx, Jx )
        class(cgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxx,tauxy,tauxz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qx
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, this%mix%ns), intent(in) :: Jx

        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: flux
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        integer :: i

        xtmp1 => this%xbuf(:,:,:,1); xtmp2 => this%xbuf(:,:,:,2)

        select case(this%mix%ns)
        case(1)
            flux = this%Wcnsrv(:,:,:,mom_index  )   ! mass
            call transpose_y_to_x(flux,xtmp1,this%decomp)
            call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
            call transpose_x_to_y(xtmp2,flux,this%decomp)
            rhs(:,:,:,1) = rhs(:,:,:,1) - flux
        case default
            do i = 1,this%mix%ns
                flux = this%Wcnsrv(:,:,:,i)*this%u + Jx(:,:,:,i)   ! mass
                call transpose_y_to_x(flux,xtmp1,this%decomp)
                call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
                call transpose_x_to_y(xtmp2,flux,this%decomp)
                rhs(:,:,:,i) = rhs(:,:,:,i) - flux
            end do
        end select

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
                            qy, Jy )
        class(cgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxy,tauyy,tauyz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qy
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, this%mix%ns), intent(in) :: Jy

        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: flux
        real(rkind), dimension(:,:,:), pointer :: ytmp1
        integer :: i

        ytmp1 => this%ybuf(:,:,:,6)

        select case(this%mix%ns)
        case(1)
            flux = this%Wcnsrv(:,:,:,mom_index+1)   ! mass
            call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
            rhs(:,:,:,1) = rhs(:,:,:,1) - ytmp1
        case default
            do i = 1,this%mix%ns
                flux = this%Wcnsrv(:,:,:,i)*this%v + Jy(:,:,:,i) ! mass
                call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
                rhs(:,:,:,i) = rhs(:,:,:,i) - ytmp1
            end do
        end select

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
                            qz, Jz )
        class(cgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxz,tauyz,tauzz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, this%mix%ns), intent(in) :: Jz

        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: flux
        real(rkind), dimension(:,:,:), pointer :: ztmp1,ztmp2
        integer :: i

        ztmp1 => this%zbuf(:,:,:,1); ztmp2 => this%zbuf(:,:,:,2)

        select case(this%mix%ns)
        case(1)
            flux = this%Wcnsrv(:,:,:,mom_index+2)   ! mass
            call transpose_y_to_z(flux,ztmp1,this%decomp)
            call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
            call transpose_z_to_y(ztmp2,flux,this%decomp)
            rhs(:,:,:,1) = rhs(:,:,:,1) - flux
        case default
            do i = 1,this%mix%ns
                flux = this%Wcnsrv(:,:,:,i)*this%w + Jz(:,:,:,i)   ! mass
                call transpose_y_to_z(flux,ztmp1,this%decomp)
                call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
                call transpose_z_to_y(ztmp2,flux,this%decomp)
                rhs(:,:,:,i) = rhs(:,:,:,i) - flux
            end do
        end select

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

    subroutine getLAD(this, dudx, dudy, dudz,&
                            dvdx, dvdy, dvdz,&
                            dwdx, dwdy, dwdz,&
                           dYsdx,dYsdy,dYsdz )
        use reductions, only: P_MAXVAL
        class(cgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, this%mix%ns), optional, intent(in) :: dYsdx,dYsdy,dYsdz
        
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: mustar,bulkstar,kapstar,diffstar,func
        integer :: i
        
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

        ! Now, all ytmps are free to use
        call this%mix%get_sos(this%rho,this%p,ytmp1)  ! Speed of sound

        kapstar = this%Ckap*this%rho*ytmp1*abs(kapstar)/this%T

        ! Filter kapstar
        call this%filter(kapstar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)

        ! -------- Artificial Diffusivity ---------
        if (present(dYsdx) .AND. present(dYsdy) .AND. present(dYsdz)) then
            if (this%mix%ns .GT. 1) then
                do i = 1,this%mix%ns
                    ! Step 1: Get components of grad(Ys) squared individually
                    ytmp1 = dYsdx(:,:,:,i)*dYsdx(:,:,:,i)
                    ytmp2 = dYsdy(:,:,:,i)*dYsdy(:,:,:,i)
                    ytmp3 = dYsdz(:,:,:,i)*dYsdz(:,:,:,i)

                    call this%mix%get_sos(this%rho,this%p,ytmp6)  ! Speed of sound

                    ! Step 2: Get 4th derivative in X
                    call transpose_y_to_x(this%Ys(:,:,:,i),xtmp1,this%decomp)
                    call this%der%d2dx2(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
                    call this%der%d2dx2(xtmp2,xtmp1,this%x_bc(1),this%x_bc(2))
                    xtmp2 = xtmp1*this%dx**4
                    call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
                    diffstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

                    ! Step 3: Get 4th derivative in Z
                    call transpose_y_to_z(this%Ys(:,:,:,i),ztmp1,this%decomp)
                    call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
                    call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
                    ztmp2 = ztmp1*this%dz**4
                    call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
                    diffstar = diffstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

                    ! Step 4: Get 4th derivative in Y
                    call this%der%d2dy2(this%Ys(:,:,:,i),ytmp4,this%y_bc(1),this%y_bc(2))
                    call this%der%d2dy2(ytmp4,ytmp5,this%y_bc(1),this%y_bc(2))
                    ytmp4 = ytmp5*this%dy**4
                    diffstar = diffstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero


                    diffstar = this%Cdiff*ytmp6*abs(diffstar) ! CD part of diff

                    ytmp4 = (sqrt(ytmp1)*this%dx + sqrt(ytmp2)*this%dy + sqrt(ytmp3)*this%dz) &
                                    / ( sqrt(ytmp1+ytmp2+ytmp3) + real(1.0D-32,rkind) ) ! grid scale
                    ytmp5 = this%CY*ytmp6*( half*(abs(this%Ys(:,:,:,i))-one + abs(this%Ys(:,:,:,i)-one)) )*ytmp4 ! CY part of diff

                    diffstar = max(diffstar, ytmp5) ! Take max of both terms instead of add to minimize the dissipation

                    ! Filter diffstar
                    call this%filter(diffstar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)

                    ! Add to physical diffusivity
                    this%diff(:,:,:,i) = this%diff(:,:,:,i) + diffstar
                end do
            end if
        end if

        ! Now, add to physical fluid properties
        this%mu   = this%mu   + mustar
        this%bulk = this%bulk + bulkstar
        this%kap  = this%kap  + kapstar

    end subroutine

    subroutine filter(this,arr,myfil,numtimes, x_bc, y_bc, z_bc)
        class(cgrid), target, intent(inout) :: this
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
        class(cgrid), intent(inout) :: this

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
        class(cgrid), target, intent(inout) :: this
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
        class(cgrid), target, intent(inout) :: this
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
        class(cgrid), target, intent(inout) :: this
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
        class(cgrid), intent(inout) :: this
        character(len=clen) :: charout
        integer :: i

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

        ! TODO: Add hook_viz here

        ! End visualization dump
        call this%viz%end_viz()
    end subroutine

    subroutine write_restart(this)
        use exits, only: message
        class(cgrid), intent(inout) :: this
        character(len=clen) :: charout
        integer :: i

        ! Get updated conserved variables
        call this%get_conserved()

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
    end subroutine

    subroutine read_restart(this, vizcount)
        use exits, only: message
        class(cgrid), intent(inout) :: this
        integer,      intent(in)    :: vizcount
        character(len=clen) :: charout
        integer :: i
        integer, dimension(1) :: tmp_int
        real(rkind), dimension(1) :: tmp_real

        ! Start reading restart file
        call this%restart%start_reading(vizcount)

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

        call this%restart%read_attribute(1, tmp_real, 'Time')
        this%tsim = tmp_real(1)

        ! End visualization dump
        call this%restart%end_reading()

        ! Get primitive variables from conserved
        call this%get_primitive()
    end subroutine

end module 
