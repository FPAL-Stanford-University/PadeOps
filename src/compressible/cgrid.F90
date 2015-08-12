module CompressibleGrid
    use kind_parameters, only: rkind, clen
    use constants, only: zero,one,two
    use GridMod, only: grid
    use hooks, only: meshgen, initfields
    use decomp_2d, only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize
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
        contains
        procedure :: init
        procedure :: destroy
        procedure :: laplacian
        procedure :: gradient 
    end type

contains
    subroutine init(this, inputfile )
        class(cgrid), intent(inout) :: this
        character(len=clen), intent(in) :: inputfile  

        integer :: nx, ny, nz
        character(len=clen) :: outputdir
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

        namelist /CINPUT/ nx, ny, nz, outputdir,periodicx, periodicy, periodicz, &
                                     derivative_x, derivative_y, derivative_z,   &
                                     filter_x, filter_y, filter_z, prow, pcol    


        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=CINPUT)
        close(ioUnit)

        this%nx = nx
        this%ny = ny
        this%nz = nz

        ! Initialize decomp
        call decomp_2d_init(nx, ny, nz, prow, pcol)
        call get_decomp_info(this%decomp)
        
        ! Allocate mesh
        if ( allocated(this%mesh) ) deallocate(this%mesh) 
        allocate(this%mesh(this%decomp%ysz(1),this%decomp%ysz(1),this%decomp%ysz(1),3))

        ! Generate default mesh: X \in [-1, 1), Y \in [-1, 1), Z \in [-1, 1)
        this%dx = two/nx
        this%dy = two/ny
        this%dz = two/nz

        ! Generate default mesh 
        do k = 1,this%decomp%ysz(3)
            do j = 1,this%decomp%ysz(2)
                do i = 1,this%decomp%ysz(1)
                    this%mesh(i,j,k,1) = -one + (i - 1)*this%dx           
                    this%mesh(i,j,k,2) = -one + (j - 1)*this%dy           
                    this%mesh(i,j,k,3) = -one + (k - 1)*this%dz           
                end do 
            end do 
        end do  

        ! Go to hooks if a different mesh is desired 
        call meshgen(nx, ny, nz, this%decomp%yst, this%decomp%yen, this%decomp%ysz, &
                    this%dx, this%dy, this%dz, this%mesh) 
        

   
        ! Allocate fields
        if ( allocated(this%fields) ) deallocate(this%fields) 
        allocate(this%fields(this%decomp%ysz(1),this%decomp%ysz(1),this%decomp%ysz(1),10))
       
        ! Initialize everything to a constant Zero
        this%fields = zero  

        ! Go to hooks if a different initialization is derired 
        call initfields(nx, ny, nz, this%decomp%yst, this%decomp%yen, this%decomp%ysz, &
                    this%dx, this%dy, this%dz, size(this%fields,4), this%mesh, this%fields) 

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
        this%nx_proc = this%decomp%ysz(1)
        this%ny_proc = this%decomp%ysz(2)
        this%nz_proc = this%decomp%ysz(3)

    end subroutine


    subroutine destroy(this)
        class(cgrid), intent(inout) :: this

        if (allocated(this%mesh)) deallocate(this%mesh) 
        if (allocated(this%fields)) deallocate(this%fields) 
        call this%der%destroy()
        call this%fil%destroy()
        call decomp_2d_finalize
    end subroutine

    subroutine gradient(this, f, dfdx, dfdy, dfdz)
        class(cgrid), intent(in) :: this
        real(rkind), intent(in), dimension(this%nx_proc, this%ny_proc, this%nz_proc) :: f
        real(rkind), intent(out), dimension(this%nx_proc, this%ny_proc, this%nz_proc) :: dfdx
        real(rkind), intent(out), dimension(this%nx_proc, this%ny_proc, this%nz_proc) :: dfdy
        real(rkind), intent(out), dimension(this%nx_proc, this%ny_proc, this%nz_proc) :: dfdz

        ! Set some shit to avoid warnings
        dfdx = f - f
        dfdy = f - f
        dfdz = f - f
    end subroutine 

    subroutine laplacian(this, f, lapf)
        class(cgrid), intent(in) :: this
        real(rkind), intent(in), dimension(this%nx_proc, this%ny_proc, this%nz_proc) :: f
        real(rkind), intent(out), dimension(this%nx_proc, this%ny_proc, this%nz_proc) :: lapf

        ! Set some shit to avoid warnings
        lapf = f - f

    end subroutine 
end module 
