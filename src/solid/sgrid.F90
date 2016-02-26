module SolidGrid
    use kind_parameters, only: rkind, clen
    use constants, only: zero,half,one,two,three,four
    use FiltersMod, only: filters
    use GridMod, only: grid
    use gridtools, only: alloc_buffs, destroy_buffs
    use hooks, only: meshgen, initfields, hook_output, hook_bc, hook_timestep
    use decomp_2d, only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,  only: derivatives
    use StiffGasEOS,     only: stiffgas
    use Sep1SolidEOS,    only: sep1solid
   
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

    integer, parameter :: nfields = 20

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
        
        real(rkind), dimension(:,:,:), pointer :: eel
         
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
            procedure, private :: getRHS
            procedure, private :: getRHS_x
            procedure, private :: getRHS_y
            procedure, private :: getRHS_z
            procedure          :: getSGS
            procedure          :: filter
            procedure          :: getPhysicalProperties
            procedure, private :: get_tau
            procedure, private :: get_q
    end type

contains
    subroutine init(this, inputfile )
        use reductions, only: P_MAXVAL
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
        real(rkind) :: Cbeta = 1.75_rkind
        real(rkind) :: Ckap = 0.01_rkind

        real(rkind), dimension(:,:,:,:), allocatable :: finger, fingersq
        real(rkind), dimension(:,:,:),   allocatable :: trG, trG2, detG

        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                             inputdir, outputdir, vizprefix, tviz, &
                                  periodicx, periodicy, periodicz, &
                         derivative_x, derivative_y, derivative_z, &
                                     filter_x, filter_y, filter_z, &
                                                       prow, pcol, &
                                                         SkewSymm  
        namelist /SINPUT/  gam, Rgas, PInf, shmod, rho0, Cmu, Cbeta, Ckap

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

        this%step = 0
        this%nsteps = nsteps

        this%rho0 = rho0

        this%Cmu = Cmu
        this%Cbeta = Cbeta
        this%Ckap = Ckap

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
        if ( allocated(this%sgas) ) deallocate(this%sgas)
        allocate(this%sgas)
        call this%sgas%init(gam,Rgas,PInf)

        ! Allocate elastic
        if ( allocated(this%elastic) ) deallocate(this%elastic)
        allocate(this%elastic)
        call this%elastic%init(shmod)

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
        
        this%eel  => this%fields(:,:,:, eel_index)   
       
        ! Initialize everything to a constant Zero
        this%fields = zero  

        ! Go to hooks if a different initialization is derired 
        call initfields(this%decomp, this%dx, this%dy, this%dz, inputdir, this%mesh, this%fields, this%rho0)
       
        ! Get hydrodynamic and elastic energies 
        call this%sgas%get_e_from_p(this%rho,this%p,this%e)

        call alloc_buffs(finger,  6,"y",this%decomp)
        call alloc_buffs(fingersq,6,"y",this%decomp)
        allocate( trG (this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )
        allocate( trG2(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )
        allocate( detG(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) )

        call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG)
        call this%elastic%get_eelastic(this%rho0,trG,trG2,detG,this%eel)
       
        this%e = this%e + this%eel
        call this%sgas%get_T(this%e,this%T)

        if (P_MAXVAL(abs( this%rho/this%rho0/(detG)**half - one )) > 10._rkind*eps) then
            call warning("Inconsistent initialization: rho/rho0 and g are not compatible")
        end if

        deallocate( finger   )
        deallocate( fingersq )
        deallocate( trG      )
        deallocate( trG2     )
        deallocate( detG     )

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

        allocate(this%viz)
        call this%viz%init(this%outputdir, vizprefix, nfields, varnames)
        this%tviz = tviz

    end subroutine


    subroutine destroy(this)
        class(sgrid), intent(inout) :: this

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
        real(rkind) :: cputime
        real(rkind), dimension(:,:,:,:), allocatable, target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        call this%get_dt()

        allocate( duidxj(this%nxp, this%nyp, this%nzp, 9) )
        ! Get artificial properties for initial conditions
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        call this%gradient(this%u,dudx,dudy,dudz)
        call this%gradient(this%v,dvdx,dvdy,dvdz)
        call this%gradient(this%w,dwdx,dwdy,dwdz)

        call this%getPhysicalProperties()

        call this%getSGS(dudx,dudy,dudz,&
                         dvdx,dvdy,dvdz,&
                         dwdx,dwdy,dwdz )
        deallocate( duidxj )
        ! ------------------------------------------------

        ! Write out initial conditions
        call hook_output(this%decomp, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%tsim, this%viz%vizcount)
        call this%viz%WriteViz(this%decomp, this%mesh, this%fields, this%tsim)
        vizcond = .FALSE.
        
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
            call message(2,"CPU time (in seconds)",cputime)
          
            ! Write out vizualization dump if vizcond is met 
            if (vizcond) then
                call hook_output(this%decomp, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%tsim, this%viz%vizcount)
                call this%viz%WriteViz(this%decomp, this%mesh, this%fields, this%tsim)
                vizcond = .FALSE.
            end if
            
            ! Get the new time step
            call this%get_dt()
            
            ! Check for visualization condition and adjust time step
            if ( (this%tviz > zero) .AND. (this%tsim + this%dt >= this%tviz * this%viz%vizcount) ) then
                this%dt = this%tviz * this%viz%vizcount - this%tsim
                vizcond = .TRUE.
            end if

            ! Check tstop condition
            if ( (this%tstop > zero) .AND. (this%tsim >= this%tstop) ) then
                tcond = .FALSE.
            else if ( (this%tstop > zero) .AND. (this%tsim + this%dt >= this%tstop) ) then
                this%dt = this%tstop - this%tsim
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

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,5) :: rhs   ! RHS for conserved variables
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,5) :: Qtmp  ! Temporary variable for RK45
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: rhsg  ! RHS for g tensor equation
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: Qtmpg ! Temporary variable for RK45
        integer :: isub,i,j,k,l

        character(len=clen) :: charout

        ! if (this%step == 0) then
        !     call this%gas%get_e_from_p(this%rho,this%p,this%e)
        !     call this%gas%get_T(this%e,this%T)
        ! end if

        Qtmp  = zero
        Qtmpg = zero

        do isub = 1,RK45_steps
            call this%get_conserved()

            if ( nancheck(this%Wcnsrv,i,j,k,l) ) then
                call message("Wcnsrv: ",this%Wcnsrv(i,j,k,l))
                write(charout,'(A,I1,A,I5,A,4(I5,A))') "NaN encountered in solution at substep ", isub, " of step ", this%step+1, " at (",i,", ",j,", ",k,", ",l,") of Wcnsrv"
                call GracefulExit(trim(charout), 999)
            end if

            call this%getRHS(rhs, rhsg)
            Qtmp  = this%dt*rhs  + RK45_A(isub)*Qtmp
            Qtmpg = this%dt*rhsg + RK45_A(isub)*Qtmpg
            this%Wcnsrv = this%Wcnsrv + RK45_B(isub)*Qtmp
            this%g      = this%g      + RK45_B(isub)*Qtmpg

            ! Filter the conserved variables
            do i = 1,5
                call this%filter(this%Wcnsrv(:,:,:,i), this%fil, 1)
            end do
            
            ! Filter the g tensor
            do i = 1,9
                call this%filter(this%g(:,:,:,i), this%fil, 1)
            end do
            
            call this%get_primitive()
            call hook_bc(this%decomp, this%mesh, this%fields, this%tsim)
        end do

        this%tsim = this%tsim + this%dt
        this%step = this%step + 1
            
        call message(1,"Time",this%tsim)
        call message(2,"Minimum density",P_MINVAL(this%rho))
        call message(2,"Maximum velocity",sqrt(P_MAXVAL(this%u*this%u+this%v*this%v+this%w*this%w)))
        call message(2,"Minimum pressure",P_MINVAL(this%p))
        call hook_timestep(this%decomp, this%mesh, this%fields, this%tsim)

    end subroutine

    subroutine get_dt(this)
        class(sgrid), target, intent(inout) :: this

        if (this%CFL > zero) then
            return
        end if

        ! Use fixed time step if CFL <= 0
        this%dt = this%dtfixed

    end subroutine

    pure subroutine get_primitive(this)
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
        
        call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG)
        call this%elastic%get_eelastic(this%rho0,trG,trG2,detG,this%eel)
        
        call this%sgas%get_T(this%e,this%T)
        call this%sgas%get_p(this%rho,(this%e-this%eel),this%p)

    end subroutine

    pure subroutine get_conserved(this)
        class(sgrid), intent(inout) :: this

        this%Wcnsrv(:,:,:,1) = this%rho
        this%Wcnsrv(:,:,:,2) = this%rho * this%u
        this%Wcnsrv(:,:,:,3) = this%rho * this%v
        this%Wcnsrv(:,:,:,4) = this%rho * this%w
        this%Wcnsrv(:,:,:,5) = this%rho * ( this%e + half*( this%u*this%u + this%v*this%v + this%w*this%w ) )

    end subroutine

    subroutine getRHS(this, rhs)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,5), intent(out) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        real(rkind), dimension(:,:,:), pointer :: qx,qy,qz

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        call this%gradient(this%u,dudx,dudy,dudz)
        call this%gradient(this%v,dvdx,dvdy,dvdz)
        call this%gradient(this%w,dwdx,dwdy,dwdz)

        call this%getPhysicalProperties()

        call this%getSGS(dudx,dudy,dudz,&
                         dvdx,dvdy,dvdz,&
                         dwdx,dwdy,dwdz )

        ! Get tau tensor and q (heat conduction) vector. Put in components of duidxj
        call this%get_tau( duidxj )
        call this%get_q  ( duidxj )

        ! Now, associate the pointers to understand what's going on better
        tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx); tauxz => duidxj(:,:,:,tauxzidx);
                                         tauyy => duidxj(:,:,:,tauyyidx); tauyz => duidxj(:,:,:,tauyzidx);
                                                                          tauzz => duidxj(:,:,:,tauzzidx);
        
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
        integer :: i

        flux => this%ybuf(:,:,:,1:5)
        xtmp1 => this%xbuf(:,:,:,1); xtmp2 => this%xbuf(:,:,:,2)

        flux(:,:,:,1) = this%Wcnsrv(:,:,:,2)   ! rho*u
        flux(:,:,:,2) = this%Wcnsrv(:,:,:,2)*this%u + this%p - tauxx
        flux(:,:,:,3) = this%Wcnsrv(:,:,:,2)*this%v          - tauxy
        flux(:,:,:,4) = this%Wcnsrv(:,:,:,2)*this%w          - tauxz
        flux(:,:,:,5) = (this%Wcnsrv(:,:,:,5) + this%p - tauxx)*this%u - this%v*tauxy - this%w*tauxz - qx

        ! Now, get the x-derivative of the fluxes
        do i=1,5
            call transpose_y_to_x(flux(:,:,:,i),xtmp1,this%decomp)
            call this%der%ddx(xtmp1,xtmp2)
            call transpose_x_to_y(xtmp2,flux(:,:,:,i),this%decomp)
        end do

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
        integer :: i

        flux => this%ybuf(:,:,:,1:5)
        ytmp1 => this%ybuf(:,:,:,6)

        flux(:,:,:,1) = this%Wcnsrv(:,:,:,3)   ! rho*v
        flux(:,:,:,2) = this%Wcnsrv(:,:,:,3)*this%u          - tauxy
        flux(:,:,:,3) = this%Wcnsrv(:,:,:,3)*this%v + this%p - tauyy
        flux(:,:,:,4) = this%Wcnsrv(:,:,:,3)*this%w          - tauyz
        flux(:,:,:,5) = (this%Wcnsrv(:,:,:,5) + this%p - tauyy)*this%v - this%u*tauxy - this%w*tauyz - qy

        ! Now, get the x-derivative of the fluxes
        do i=1,5
            call this%der%ddy(flux(:,:,:,i),ytmp1)

            ! Add to rhs
            rhs(:,:,:,i) = rhs(:,:,:,i) - ytmp1
        end do

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
        integer :: i

        flux => this%ybuf(:,:,:,1:5)
        ztmp1 => this%zbuf(:,:,:,1); ztmp2 => this%zbuf(:,:,:,2)

        flux(:,:,:,1) = this%Wcnsrv(:,:,:,4)   ! rho*w
        flux(:,:,:,2) = this%Wcnsrv(:,:,:,4)*this%u          - tauxz
        flux(:,:,:,3) = this%Wcnsrv(:,:,:,4)*this%v          - tauyz
        flux(:,:,:,4) = this%Wcnsrv(:,:,:,4)*this%w + this%p - tauzz
        flux(:,:,:,5) = (this%Wcnsrv(:,:,:,5) + this%p - tauzz)*this%w - this%u*tauxz - this%v*tauyz - qz

        ! Now, get the x-derivative of the fluxes
        do i=1,5
            call transpose_y_to_z(flux(:,:,:,i),ztmp1,this%decomp)
            call this%der%ddz(ztmp1,ztmp2)
            call transpose_z_to_y(ztmp2,flux(:,:,:,i),this%decomp)
        end do

        ! Add to rhs
        rhs = rhs - flux

    end subroutine

    subroutine getSGS(this,dudx,dudy,dudz,&
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
        call this%der%d2dx2(xtmp1,xtmp2)
        call this%der%d2dx2(xtmp2,xtmp1)
        xtmp2 = xtmp1*this%dx**6
        call transpose_x_to_y(xtmp2,mustar,this%decomp)
        
        ! Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2)
        call this%der%d2dz2(ztmp2,ztmp1)
        ztmp2 = ztmp1*this%dz**6
        call transpose_z_to_y(ztmp2,ytmp1,this%decomp)
        mustar = mustar + ytmp1
        
        ! Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp1)
        call this%der%d2dy2(ytmp1,ytmp2)
        ytmp1 = ytmp2*this%dy**6
        mustar = mustar + ytmp1

        mustar = this%Cmu*this%rho*abs(mustar)
        
        ! Filter mustar
        call this%filter(mustar, this%gfil, 2)
        
        ! -------- Artificial Bulk Viscosity --------
        
        func = dudx + dvdy + dwdz      ! dilatation
        
        ! Step 1: Get components of grad(rho) squared individually
        call this%gradient(this%rho,ytmp1,ytmp2,ytmp3) ! Does not use any Y buffers
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(func,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2)
        call this%der%d2dx2(xtmp2,xtmp1)
        xtmp2 = xtmp1*this%dx**6
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        bulkstar = ytmp4 * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind))

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2)
        call this%der%d2dz2(ztmp2,ztmp1)
        ztmp2 = ztmp1*this%dz**6
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        bulkstar = bulkstar + ytmp4 * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind))

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp4)
        call this%der%d2dy2(ytmp4,ytmp5)
        ytmp4 = ytmp5*this%dy**6
        bulkstar = bulkstar + ytmp4 * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind))

        ! Now, all ytmps are free to use
        ytmp1 = dwdy-dvdz; ytmp2 = dudz-dwdx; ytmp3 = dvdx-dudy
        ytmp4 = ytmp1*ytmp1 + ytmp2+ytmp2 + ytmp3*ytmp3 ! |curl(u)|^2
        ytmp2 = func*func

        ! Calculate the switching function
        ytmp1 = ytmp2 / (ytmp2 + ytmp4 + real(1.0D-32,rkind)) ! Switching function f_sw
        where (func .GE. zero)
            ytmp1 = zero
        end where

        bulkstar = this%Cbeta*this%rho*ytmp1*abs(bulkstar)

        ! Filter bulkstar
        call this%filter(bulkstar, this%gfil, 2)

        ! -------- Artificial Conductivity --------

        ! Step 1: Get components of grad(e) squared individually
        call this%gradient(this%e,ytmp1,ytmp2,ytmp3) ! Does not use any Y buffers
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(this%e,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2)
        call this%der%d2dx2(xtmp2,xtmp1)
        xtmp2 = xtmp1*this%dx**6
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        kapstar = ytmp4 * ytmp1 / (ytmp1 + ytmp2 + ytmp3)

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(this%e,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2)
        call this%der%d2dz2(ztmp2,ztmp1)
        ztmp2 = ztmp1*this%dz**6
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        kapstar = kapstar + ytmp4 * ytmp3 / (ytmp1 + ytmp2 + ytmp3)

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(this%e,ytmp4)
        call this%der%d2dy2(ytmp4,ytmp5)
        ytmp4 = ytmp5*this%dy**6
        kapstar = kapstar + ytmp4 * ytmp2 / (ytmp1 + ytmp2 + ytmp3)

        ! Now, all ytmps are free to use
        call this%gas%get_sos(this%rho,this%p,ytmp1)  ! Speed of sound

        kapstar = this%Ckap*this%rho*ytmp1*abs(kapstar)/this%T

        ! Filter kapstar
        call this%filter(kapstar, this%gfil, 2)

        ! Now, add to physical fluid properties
        this%mu   = this%mu   + mustar
        this%bulk = this%bulk + bulkstar
        this%kap  = this%kap  + kapstar

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
        call der%ddy(this%T,tmp1_in_y)
        duidxj(:,:,:,qyidx) = -this%kap*tmp1_in_y

        ! Step 2: Get qx (dudx is destroyed)
        call transpose_y_to_x(this%T,tmp1_in_x,this%decomp)
        call der%ddx(tmp1_in_x,tmp2_in_x)
        call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
        duidxj(:,:,:,qxidx) = -this%kap*tmp1_in_y

        ! Step 3: Get qz (dwdz is destroyed)
        call transpose_y_to_z(this%T,tmp1_in_z,this%decomp)
        call der%ddz(tmp1_in_z,tmp2_in_z)
        call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
        duidxj(:,:,:,qzidx) = -this%kap*tmp1_in_y

        ! Done
    end subroutine 
end module 
