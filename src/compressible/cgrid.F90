module CompressibleGrid
    use kind_parameters, only: rkind, clen
    use constants, only: zero,one,two
    use FiltersMod, only: filters
    use GridMod, only: grid, alloc_buffs, destroy_buffs
    use hooks, only: meshgen, initfields
    use decomp_2d, only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,  only: derivatives
    use IdealGasEOS,     only: idealgas
   
    implicit none

    integer :: rho_index    = 1 
    integer :: u_index      = 2
    integer :: v_index      = 3
    integer :: w_index      = 4
    integer :: p_index      = 5
    integer :: T_index      = 6
    integer :: e_index      = 7
    integer :: mu_index     = 8
    integer :: bulk_index   = 9
    integer :: kap_index    = 10

    type, extends(grid) :: cgrid
       
        type(filters) :: gfil
        type(idealgas), allocatable :: gas

        real(rkind), dimension(:,:,:,:) :: W                               ! Conserved variables
        real(rkind), dimension(:,:,:,:), allocatable :: xbuf, ybuf, zbuf   ! Buffers
        
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
         
        contains
            procedure :: init
            procedure :: destroy
            procedure :: laplacian
            procedure :: gradient 
    end type

contains
    subroutine init(this, inputfile )
        class(cgrid),target, intent(inout) :: this
        character(len=clen), intent(in) :: inputfile  

        integer :: nx, ny, nz
        character(len=clen) :: outputdir
        character(len=clen) :: inputdir
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
        integer :: nsteps = -1
        real(rkind) :: dt = -one
        real(rkind) :: tstop = one
        real(rkind) :: CFL = -one
        logical :: SkewSymm = .FALSE.

        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                                              inputdir, outputdir, &
                                  periodicx, periodicy, periodicz, &
                         derivative_x, derivative_y, derivative_z, &
                                     filter_x, filter_y, filter_z, &
                                                       prow, pcol, &
                                                         SkewSymm  &
        namelist /CINPUT/  gam, Rgas


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

        this%step = 0
        this%nsteps = nsteps

        ! Initialize decomp
        call decomp_2d_init(nx, ny, nz, prow, pcol)
        call get_decomp_info(this%decomp)
        
        ! Allocate mesh
        if ( allocated(this%mesh) ) deallocate(this%mesh) 
        call alloc_buffs(this%mesh,3,'y',this%decomp)

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

        allocate(this%gas)
        call this%gas%init(gam,Rgas)

        ! Go to hooks if a different mesh is desired 
        call meshgen(this%decomp, this%dx, this%dy, this%dz, this%mesh) 

   
        ! Allocate fields
        if ( allocated(this%fields) ) deallocate(this%fields) 
        call alloc_buffs(this%fields,10,'y',this%decomp)
        call alloc_buffs(this%W,5,'y',this%decomp)
       
        ! Initialize everything to a constant Zero
        this%fields = zero  

        ! Go to hooks if a different initialization is derired 
        call initfields(this%decomp, this%dx, this%dy, this%dz, inputdir, this%mesh, this%fields)

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
        
        ! Initialize derivatives 
        call this%der%init(                           this%decomp, &
                           this%dx,       this%dy,        this%dz, &
                         periodicx,     periodicy,      periodicz, &
                      derivative_x,  derivative_y,   derivative_z, &
                           .false.,       .false.,        .false., &
                           .false.)      

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
        call alloc_buffs(this%xbuf,2,"x",this%decomp)
        call alloc_buffs(this%ybuf,1,"y",this%decomp)
        call alloc_buffs(this%zbuf,2,"z",this%decomp)

        ! Associate pointers for ease of use
        x    => this%mesh  (:,:,:, 1) 
        y    => this%mesh  (:,:,:, 2) 
        z    => this%mesh  (:,:,:, 3)

        rho  => this%fields(:,:,:, 1) 
        u    => this%fields(:,:,:, 2) 
        v    => this%fields(:,:,:, 3) 
        w    => this%fields(:,:,:, 4)  
        p    => this%fields(:,:,:, 5)  
        T    => this%fields(:,:,:, 6)  
        e    => this%fields(:,:,:, 7)  
        mu   => this%fields(:,:,:, 8)  
        bulk => this%fields(:,:,:, 9)  
        kap  => this%fields(:,:,:,10)   

        this%SkewSymm = SkewSymm

    end subroutine


    subroutine destroy(this)
        class(cgrid), intent(inout) :: this

        if (allocated(this%mesh)) deallocate(this%mesh) 
        if (allocated(this%fields)) deallocate(this%fields) 
        call this%der%destroy()
        call this%fil%destroy()
        call this%gfil%destroy()
        call destroy_buffs(this%xbuf)
        call destroy_buffs(this%ybuf)
        call destroy_buffs(this%zbuf)
        if (allocated(this%gas)) deallocate(this%gas) 
        if (allocated(this%W)) deallocate(this%W) 
        call decomp_2d_finalize

    end subroutine

    subroutine gradient(this, f, dfdx, dfdy, dfdz)
        class(cgrid),target, intent(inout) :: this
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
        class(cgrid),target, intent(inout) :: this
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

    subroutine advance_RK45(this)
        use RKCoeffs, only: RK45_steps,RK45_A,RK45_B
        class(cgrid), target, intent(inout) :: this

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,5) :: rhs  ! RHS for conserved variables
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,5) :: Qtmp ! Temporary variable for RK45
        integer :: isub

        Qtmp = zero

        do isub = 1,RK45_steps
            call this%get_conservative()

            call this%getRHS(rhs)
            Qtmp = this%dt*rhs + A(isub)*Qtmp
            W = W + B(isub)*Qtmp

            call this%filter()
            call this%get_primitive()
        end do

    end subroutine

    subroutine getRHS(this, rhs)
        class(cgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,5), intent(out) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9) :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        call this%gradient(u,dudx,dudy,dudz)
        call this%gradient(v,dvdx,dvdy,dvdz)
        call this%gradient(w,dwdx,dwdy,dwdz)

        call this%getSGS(dudx,dudy,dudz,&
                         dvdx,dvdy,dvdz,&
                         dwdx,dwdy,dwdz )

        call GetInviscidRHS(rhs)
        call GetViscousFluxes(rhs)

    end subroutine

    subroutine getSGS(this,dudx,dudy,dudz,&
                           dvdx,dvdy,dvdz,&
                           dwdx,dwdy,dwdz )
        class(cgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: mustar,bulkstar,kapstar,func
        
        real(rkind), dimension(:,:,:) :: xtmp1,xtmp2,xtmp3
        real(rkind), dimension(:,:,:) :: ytmp1,ytmp2,ytmp3,ytmp4,grho1,grho2,grho3
        real(rkind), dimension(:,:,:) :: ztmp1,ztmp2,ztmp3

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

        mustar = Cmu*rho*abs(mustar)

        ! Filter mustar
        call this%gfilter(mustar,ytmp1)
        call this%gfilter(ytmp1,mustar)

        ! -------- Artificial Bulk Viscosity --------
        
        func = dudx + dvdy + dwdz      ! dilatation
        
        ! Step 1: Get components of grad(rho) squared individually
        call this%gradient(rho,ytmp1,ytmp2,ytmp3)
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(func,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2)
        call this%der%d2dx2(xtmp2,xtmp1)
        xtmp2 = xtmp1*this%dx**6
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        bulkstar = ytmp4 * ytmp1 / (ytmp1 + ytmp2 + ytmp3)

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2)
        call this%der%d2dz2(ztmp2,ztmp1)
        ztmp2 = ztmp1*this%dz**6
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        bulkstar = bulkstar + ytmp4 * ytmp3 / (ytmp1 + ytmp2 + ytmp3)

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp4)
        call this%der%d2dy2(ytmp4,ytmp5)
        ytmp4 = ytmp5*this%dy**6
        bulkstar = bulkstar + ytmp4 * ytmp2 / (ytmp1 + ytmp2 + ytmp3)

        ! Now, all ytmps are free to use
        ytmp1 = dwdy-dvdz; ytmp2 = dudz-dwdx; ytmp3 = dvdx-dudy
        ytmp4 = ytmp1*ytmp1 + ytmp2+ytmp2 + ytmp3*ytmp3 ! |curl(u)|^2
        ytmp2 = func*func

        ! Calculate the switching function
        ytmp1 = ytmp2 / (ytmp2 + ytmp4 + 1.0D-32_rkind) ! Switching function f_sw
        where (func .GE. zero)
            ytmp1 = zero
        end where

        bulkstar = Cbeta*rho*ytmp1*abs(bulkstar)

        ! Filter bulkstar
        call this%gfilter(bulkstar,ytmp2)
        call this%gfilter(ytmp2,bulkstar)

        ! -------- Artificial Conductivity --------

    end subroutine

    subroutine GetInviscidRHS(this,rhs)
        class(cgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,5), intent(inout) :: rhs


    end subroutine

    subroutine filter(this)
        class(cgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp
        integer :: i

        do i=1,5
            call this%fil%filterx(W(:,:,:,i),tmp)
            W(:,:,:,i) = tmp
            call this%fil%filtery(W(:,:,:,i),tmp)
            W(:,:,:,i) = tmp
            call this%fil%filterz(W(:,:,:,i),tmp)
            W(:,:,:,i) = tmp
        end do

    end subroutine
    
    subroutine gfilter(this)
        class(cgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp
        integer :: i

        do i=1,5
            call this%fil%filterx(W(:,:,:,i),tmp)
            W(:,:,:,i) = tmp
            call this%fil%filtery(W(:,:,:,i),tmp)
            W(:,:,:,i) = tmp
            call this%fil%filterz(W(:,:,:,i),tmp)
            W(:,:,:,i) = tmp
        end do

    end subroutine

end module 
