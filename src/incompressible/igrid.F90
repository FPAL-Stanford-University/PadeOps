module IncompressibleGrid
    use kind_parameters, only: rkind, clen
    use constants, only: zero,one,two,three,half 
    use GridMod, only: grid, alloc_buffs, destroy_buffs
    use hooks, only: meshgen, initfields
    use decomp_2d, only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,  only: derivatives
  
    use exits, only: GracefulExit, message

    implicit none

    integer :: u_index      = 1
    integer :: v_index      = 2
    integer :: w_index      = 3
    integer :: nu_index     = 4

    type, extends(grid) :: igrid
       
        real(rkind), dimension(:,:,:,:), allocatable :: xbuf, ybuf, zbuf 
        real(rkind), dimension(:,:,:), pointer :: x
        real(rkind), dimension(:,:,:), pointer :: y
        real(rkind), dimension(:,:,:), pointer :: z
        real(rkind), dimension(:,:,:), pointer :: u
        real(rkind), dimension(:,:,:), pointer :: v
        real(rkind), dimension(:,:,:), pointer :: w
        real(rkind), dimension(:,:,:), pointer :: nu
        real(rkind), dimension(:,:,:), pointer :: u_in_x 
        real(rkind), dimension(:,:,:), pointer :: v_in_x
        real(rkind), dimension(:,:,:), pointer :: w_in_x
        real(rkind), dimension(:,:,:), pointer :: u_in_z 
        real(rkind), dimension(:,:,:), pointer :: v_in_z
        real(rkind), dimension(:,:,:), pointer :: w_in_z
        real(rkind), dimension(:,:,:), pointer :: tmp1_in_x
        real(rkind), dimension(:,:,:), pointer :: tmp1_in_y
        real(rkind), dimension(:,:,:), pointer :: tmp1_in_z
        real(rkind), dimension(:,:,:), pointer :: tmp2_in_x
        real(rkind), dimension(:,:,:), pointer :: tmp2_in_y
        real(rkind), dimension(:,:,:), pointer :: tmp2_in_z
        
        real(rkind), dimension(:,:,:,:,:), allocatable :: rhs 
        logical :: useSGS

        contains
            procedure :: init
            procedure :: destroy
            procedure :: laplacian
            procedure :: gradient
            procedure,private :: CnsrvFrm 
            procedure,private :: NonCnsrvFrm
            procedure,private :: getRHS
            procedure,private :: getDT
            procedure,private :: addVisc
            procedure :: AdamsBashforth
            procedure :: printDivergence

    end type

contains
    subroutine init(this, inputfile )
        class(igrid),target, intent(inout) :: this
        character(len=clen), intent(in) :: inputfile  

        integer :: nx, ny, nz
        character(len=clen) :: outputdir
        character(len=clen) :: inputdir
        logical :: SkewSymm = .false. 
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
        logical :: useSGS = .false. 
        integer :: i, j, k 
        integer :: ioUnit
        integer :: nsteps = -1
        real(rkind) :: dt = -one
        real(rkind) :: tstop = one
        real(rkind) :: CFL = -one
        real(rkind) :: nu = 0.02_rkind

        namelist /INPUT/       nx, ny, nz, tstop, dt, CFL, nsteps, &
                                              inputdir, outputdir, &
                                  periodicx, periodicy, periodicz, &
                         derivative_x, derivative_y, derivative_z, &
                                     filter_x, filter_y, filter_z, &
                                                       prow, pcol, &
                                                         SkewSymm
        namelist /IINPUT/  nu, useSGS


        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=INPUT)
        read(unit=ioUnit, NML=IINPUT)
        close(ioUnit)

        this%nx = nx
        this%ny = ny
        this%nz = nz

        this%SkewSymm = SkewSymm 
        this%tsim = zero
        this%tstop = tstop
        this%CFL = CFL

        this%useSGS = useSGS
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

        ! Go to hooks if a different mesh is desired 
        call meshgen(this%decomp, this%dx, this%dy, this%dz, this%mesh) 
   
        ! Allocate fields
        if ( allocated(this%fields) ) deallocate(this%fields) 
        call alloc_buffs(this%fields,4,'y',this%decomp)
       
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

        ! Finally, set the local array dimensions
        this%nxp = this%decomp%ysz(1)
        this%nyp = this%decomp%ysz(2)
        this%nzp = this%decomp%ysz(3)


        ! Allocate 2 buffers for each of the three decompositions
        call alloc_buffs(this%xbuf,5,"x",this%decomp)
        call alloc_buffs(this%ybuf,2,"y",this%decomp)
        call alloc_buffs(this%zbuf,5,"z",this%decomp)

        ! Link Pointers
        this%u => this%fields(:,:,:,u_index)
        this%v => this%fields(:,:,:,v_index)
        this%w => this%fields(:,:,:,w_index)
        this%nu => this%fields(:,:,:,nu_index)
        this%x => this%mesh(:,:,:,1)
        this%y => this%mesh(:,:,:,2)
        this%z => this%mesh(:,:,:,3)
        
        this%u_in_x => this%xbuf(:,:,:,1)
        this%v_in_x => this%xbuf(:,:,:,2)
        this%w_in_x => this%xbuf(:,:,:,3)
        this%tmp1_in_x => this%xbuf(:,:,:,4)
        this%tmp2_in_x => this%xbuf(:,:,:,5)

        this%u_in_z => this%zbuf(:,:,:,1)
        this%v_in_z => this%zbuf(:,:,:,2)
        this%w_in_z => this%zbuf(:,:,:,3)
        this%tmp1_in_z => this%zbuf(:,:,:,4)
        this%tmp2_in_z => this%zbuf(:,:,:,5)

        this%tmp1_in_y => this%ybuf(:,:,:,1)
        this%tmp2_in_y => this%ybuf(:,:,:,2)


        ! Initialize time
        this%tsim = zero
        this%step = 0

        ! Create the storage RHS storage
        ! Assuming Adams-Bashforth
        allocate(this%rhs(this%nxp,this%nyp,this%nzp,3,2))

        if ((this%CFL < 0) .and. (this%dt < 0)) then
            call GracefulExit("Neither CFL not dt were specified in the input file", 123)
        end if 
    end subroutine


    subroutine destroy(this)
        class(igrid), intent(inout) :: this

        if (allocated(this%mesh)) deallocate(this%mesh) 
        if (allocated(this%fields)) deallocate(this%fields) 
        call this%der%destroy()
        call this%fil%destroy()
        call destroy_buffs(this%xbuf)
        call destroy_buffs(this%ybuf)
        call destroy_buffs(this%zbuf)
        call decomp_2d_finalize

    end subroutine

    subroutine gradient(this, f, dfdx, dfdy, dfdz)
        class(igrid),target, intent(inout) :: this
        real(rkind), intent(in),  dimension(this%nxp, this%nyp, this%nzp) :: f
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
        class(igrid),target, intent(inout) :: this
        real(rkind), intent(in),  dimension(this%nxp, this%nyp, this%nzp) :: f
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

        call der%d2dy2(f,lapf)
        
        call transpose_y_to_x(f,xtmp,this%decomp) 
        call this%der%d2dx2(xtmp,xdum)
        call transpose_x_to_y(xdum,ytmp,this%decomp)

        lapf = lapf + ytmp

        call transpose_y_to_z(f,ztmp,this%decomp)
        call this%der%d2dz2(ztmp,zdum)
        call transpose_z_to_y(zdum,ytmp,this%decomp)
        
        lapf = lapf + ytmp

    end subroutine


    subroutine AdamsBashforth(this)
        class(igrid), target, intent(inout) :: this
        real(rkind), dimension(:,:,:,:), pointer :: rhsold, rhs
        real(rkind) :: dtby2

        call this%getDT()
        dtby2 = this%dt/two

        rhs    => this%rhs(:,:,:,:,1)
        rhsold => this%rhs(:,:,:,:,2)

        call this%getRHS(rhs)
        this%u = this%u + dtby2*(three*rhs(:,:,:,1) - rhsold(:,:,:,1))   
        this%v = this%v + dtby2*(three*rhs(:,:,:,2) - rhsold(:,:,:,2))   
        this%w = this%w + dtby2*(three*rhs(:,:,:,3) - rhsold(:,:,:,3))   
        
        rhsold = rhs 

    end subroutine 

    subroutine getDT(this)
        class(igrid), intent(inout) :: this 

        if (this%CFL < 0) then
            return
        end if 
       
        ! Otherwise Compute dt
        this%dt = 0.00001  
    end subroutine 

    subroutine getRHS(this,rhs)
        class(igrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3), intent(out) :: rhs
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: rhs_tmp
        

        ! Now transpose fields from y->x
        call transpose_y_to_x(this%u,this%u_in_x,this%decomp)
        call transpose_y_to_x(this%v,this%v_in_x,this%decomp)
        call transpose_y_to_x(this%w,this%w_in_x,this%decomp)

        ! Now transpose fields from y->z
        call transpose_y_to_z(this%u,this%u_in_z,this%decomp)
        call transpose_y_to_z(this%v,this%v_in_z,this%decomp)
        call transpose_y_to_z(this%w,this%w_in_z,this%decomp)

        call this%CnsrvFrm(rhs)       
        if (this%SkewSymm) then
            call this%NonCnsrvFrm(rhs_tmp)
            rhs = half*rhs + half*rhs_tmp
        end if 
        
        call this%addVisc(rhs)

    end subroutine


    subroutine NonCnsrvFrm(this,rhs_tmp)
        class(igrid),target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3), intent(out) :: rhs_tmp
        type(derivatives), pointer :: der
        
        der => this%der
        
        ! Step 1: Get the x derivatives and add to rhs
        call der%ddx(this%u_in_x,this%tmp1_in_x)
        call transpose_x_to_y(this%tmp1_in_x,this%tmp1_in_y,this%decomp) 
        rhs_tmp(:,:,:,1) = -this%u*this%tmp1_in_y

        call der%ddx(this%v_in_x,this%tmp1_in_x)
        call transpose_x_to_y(this%tmp1_in_x,this%tmp1_in_y,this%decomp) 
        rhs_tmp(:,:,:,2) = -this%u*this%tmp1_in_y

        call der%ddx(this%w_in_x,this%tmp1_in_x)
        call transpose_x_to_y(this%tmp1_in_x,this%tmp1_in_y,this%decomp) 
        rhs_tmp(:,:,:,3) = -this%u*this%tmp1_in_y

        ! Step 2: Get the y derivatives and add to rhs
        call der%ddx(this%u,this%tmp1_in_y)
        rhs_tmp(:,:,:,1) = rhs_tmp(:,:,:,1) - this%v*this%tmp1_in_y

        call der%ddx(this%v,this%tmp1_in_y)
        rhs_tmp(:,:,:,2) = rhs_tmp(:,:,:,3) - this%v*this%tmp1_in_y

        call der%ddx(this%w,this%tmp1_in_y)
        rhs_tmp(:,:,:,3) = rhs_tmp(:,:,:,3) - this%v*this%tmp1_in_y
        
        ! Step 3: Get the z derivatives and add to rhs
        call der%ddz(this%u_in_z,this%tmp1_in_z)
        call transpose_z_to_y(this%tmp1_in_z,this%tmp1_in_y,this%decomp)
        rhs_tmp(:,:,:,1) = rhs_tmp(:,:,:,1) - this%w*this%tmp1_in_y

        call der%ddz(this%v_in_z,this%tmp1_in_z)
        call transpose_z_to_y(this%tmp1_in_z,this%tmp1_in_y,this%decomp)
        rhs_tmp(:,:,:,2) = rhs_tmp(:,:,:,2) - this%w*this%tmp1_in_y

        call der%ddz(this%w_in_z,this%tmp1_in_z)
        call transpose_z_to_y(this%tmp1_in_z,this%tmp1_in_y,this%decomp)
        rhs_tmp(:,:,:,3) = rhs_tmp(:,:,:,3) - this%w*this%tmp1_in_y

    end subroutine

    subroutine CnsrvFrm(this,rhs)
        class(igrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3), intent(out) :: rhs
        type(derivatives), pointer :: der

        der => this%der 
        
        ! Step 1: Get x derivatives and add to the rhs
        this%tmp1_in_x = this%u_in_x*this%u_in_x
        call der%ddx(this%tmp1_in_x,this%tmp2_in_x)
        call transpose_x_to_y(this%tmp2_in_x,this%tmp1_in_y,this%decomp)
        rhs(:,:,:,1) = - this%tmp1_in_y

        this%tmp1_in_x = this%u_in_x*this%v_in_x
        call der%ddx(this%tmp1_in_x,this%tmp2_in_x)
        call transpose_x_to_y(this%tmp2_in_x,this%tmp1_in_y,this%decomp)
        rhs(:,:,:,2) = - this%tmp1_in_y

        this%tmp1_in_x = this%u_in_x*this%w_in_x
        call der%ddx(this%tmp1_in_x,this%tmp2_in_x)
        call transpose_x_to_y(this%tmp2_in_x,this%tmp1_in_y,this%decomp)
        rhs(:,:,:,3) = - this%tmp1_in_y

        ! Step 2: Get y derivatives and add to the rhs
        this%tmp1_in_y = this%u*this%v
        call der%ddy(this%tmp1_in_y,this%tmp2_in_y)
        rhs(:,:,:,1) = rhs(:,:,:,1) - this%tmp2_in_y

        this%tmp1_in_y = this%v*this%v
        call der%ddy(this%tmp1_in_y,this%tmp2_in_y)
        rhs(:,:,:,2) = rhs(:,:,:,2) - this%tmp2_in_y

        this%tmp1_in_y = this%v*this%w
        call der%ddy(this%tmp1_in_y,this%tmp2_in_y)
        rhs(:,:,:,3) = rhs(:,:,:,3) - this%tmp2_in_y

        ! Step 3: Get z derivatives and add to the rhs
        this%tmp1_in_z = this%u_in_z*this%w_in_z
        call der%ddz(this%tmp1_in_z,this%tmp2_in_z)
        call transpose_z_to_y(this%tmp2_in_z,this%tmp1_in_y,this%decomp)
        rhs(:,:,:,1) = rhs(:,:,:,1) - this%tmp1_in_y

        this%tmp1_in_z = this%v_in_z*this%w_in_z
        call der%ddz(this%tmp1_in_z,this%tmp2_in_z)
        call transpose_z_to_y(this%tmp2_in_z,this%tmp1_in_y,this%decomp)
        rhs(:,:,:,2) = rhs(:,:,:,2) - this%tmp1_in_y

        this%tmp1_in_z = this%w_in_z*this%w_in_z
        call der%ddz(this%tmp1_in_z,this%tmp2_in_z)
        call transpose_z_to_y(this%tmp2_in_z,this%tmp1_in_y,this%decomp)
        rhs(:,:,:,3) = rhs(:,:,:,3) - this%tmp1_in_y


    end subroutine

    subroutine addVisc(this,rhs)
        class(igrid),target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3), intent(inout) :: rhs
        type(derivatives), pointer :: der

        der => this%der 

        if (this%useSGS) then
            ! Call the SGS model of choice to update "nu"
        end if
         
        ! Step 1: get 1st equation derivatives and add to its rhs
        call der%d2dx2(this%u_in_x,this%tmp1_in_x)
        call transpose_x_to_y(this%tmp1_in_x,this%tmp1_in_y,this%decomp)
        call der%d2dy2(this%u,this%tmp2_in_y)
        this%tmp2_in_y = this%tmp2_in_y  + this%tmp1_in_y
        call der%d2dz2(this%u_in_z,this%tmp1_in_z)
        call transpose_z_to_y(this%tmp1_in_z,this%tmp1_in_y,this%decomp)
        this%tmp2_in_y = this%tmp2_in_y  + this%tmp1_in_y
        rhs(:,:,:,1) = rhs(:,:,:,1) + this%nu*this%tmp2_in_y

        ! Step 2: get 2nd equation derivatives and add to its rhs
        call der%d2dx2(this%v_in_x,this%tmp1_in_x)
        call transpose_x_to_y(this%tmp1_in_x,this%tmp1_in_y,this%decomp)
        call der%d2dy2(this%v,this%tmp2_in_y)
        this%tmp2_in_y = this%tmp2_in_y  + this%tmp1_in_y
        call der%d2dz2(this%v_in_z,this%tmp1_in_z)
        call transpose_z_to_y(this%tmp1_in_z,this%tmp1_in_y,this%decomp)
        this%tmp2_in_y = this%tmp2_in_y  + this%tmp1_in_y
        rhs(:,:,:,2) = rhs(:,:,:,2) + this%nu*this%tmp2_in_y

        ! Step 3: get 2nd equation derivatives and add to its rhs
        call der%d2dx2(this%w_in_x,this%tmp1_in_x)
        call transpose_x_to_y(this%tmp1_in_x,this%tmp1_in_y,this%decomp)
        call der%d2dy2(this%w,this%tmp2_in_y)
        this%tmp2_in_y = this%tmp2_in_y  + this%tmp1_in_y
        call der%d2dz2(this%w_in_z,this%tmp1_in_z)
        call transpose_z_to_y(this%tmp1_in_z,this%tmp1_in_y,this%decomp)
        this%tmp2_in_y = this%tmp2_in_y  + this%tmp1_in_y
        rhs(:,:,:,3) = rhs(:,:,:,3) + this%nu*this%tmp2_in_y


    end subroutine 

    subroutine printDivergence(this)
        use reductions, only: p_maxval
        class(igrid),target, intent(inout) :: this
        type(derivatives), pointer :: der

        der => this%der 
        

        ! Get du_dx 
        call transpose_y_to_x(this%u,this%tmp1_in_x,this%decomp)
        call der%ddx(this%tmp1_in_x,this%tmp2_in_x) 
        call transpose_x_to_y(this%tmp2_in_x,this%tmp1_in_y,this%decomp)

        ! Get dv_dy
        call der%ddy(this%v,this%tmp2_in_y)
        this%tmp1_in_y = this%tmp1_in_y + this%tmp2_in_y

        ! Get dw_dz
        call transpose_y_to_z(this%w,this%tmp1_in_z,this%decomp)
        call der%ddz(this%tmp1_in_z,this%tmp2_in_z)
        call transpose_z_to_y(this%tmp2_in_z,this%tmp2_in_y,this%decomp)

        this%tmp1_in_y = this%tmp1_in_y + this%tmp2_in_y

        call message("Maximum divergence",p_maxval(this%tmp1_in_y))

    end subroutine

end module 
