module CompressibleGrid
    use GridMod, only: grid
    use hooks, only: meshgen, initfields
    use decomp_2d, only: decomp_info 
    use constants, only: one

    type, extends(grid) :: cgrid 
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

        contains
        procedure :: init
        procedure :: destroy
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
        character(len=clen) :: filter_x = "cd10"  
        character(len=clen) :: filter_y = "cd10" 
        character(len=clen) :: filter_z = "cd10"
        integer :: prow = 0, pcol = 0 
        integer :: i, j, k 
        integer :: ioUnit

        namelist /CINPUT/ nx, ny, nz, outputdir,periodicx, periodicy, periodicz,
                                     derivative_x, derivative_y, derivative_z, 
                                     filter_x, filter_y, filter_z, prow, pcol


        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=CINPUT)
        close(ioUnit)

        this%nx = nx
        this%ny = ny
        this%nz = nz

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

        ! Initialize decomp
        call decomp_2d_init(nx, ny, nz, prow, pcol)
        call get_decomp_info(this%decomp)
   
        ! Initialize derivatives 
        call this%der%init(                                    gp, &
                                dx,            dy,             dz, &
                         periodicx,     periodicy,      periodicz, &
                      derivative_x,  derivative_y,   derivative_z, &
                           .false.,       .false.,        .false., &
                           .false.)      

        ! Initialize filters
        call this%fil%init(                                    gp, &
                         periodicx,     periodicy,      periodicz, &
                          filter_x,      filter_y,       filter_z, &
                          )      

        ! Allocate mesh
        if ( allocated(this%mesh) ) deallocate(this%mesh) 
        allocate(this%mesh(this%decomp%ysz(1),this%decomp%ysz(1),this%decomp%ysz(1),3))

        ! Allocate fields
        if ( allocated(this%fields) ) deallocate(this%fields) 
        allocate(this%fields(this%decomp%ysz(1),this%decomp%ysz(1),this%decomp%ysz(1),10))


        ! Generate default mesh: X \in [-1, 1), Y \in [-1, 1), Z \in [-1, 1)

        this%dx = two/nx
        this%dy = two/ny
        this%dz = two/nz

        ! Generate default mesh 
        do k = 1,this%ysz(3)
            do j = 1,this%ysz(2)
                do i = 1,this%ysz(1)
                    this%mesh(i,j,k,1) = -one + (i - 1)*this%dx           
                    this%mesh(i,j,k,2) = -one + (j - 1)*this%dy           
                    this%mesh(i,j,k,3) = -one + (k - 1)*this%dz           
                end do 
            end do 
        end do  

        ! Go to hooks if a different mesh is desired 
        call meshgen(this) 
       
        ! Initialize everything to a constant Zero
        this%fields = zero  

        ! Go to hooks if a different initialization is derired 
        call initfields(this)

    end subroutine


    subroutine destroy(this)
        class(cgrid), intent(inout) :: this

        if (allocated(this%mesh)) deallocate(this%mesh) 
        if (allocated(this%fields)) deallocate(this%fields) 
        call this%der%destroy()
        call this%fil%destroy()
        call decomp_2d_finalize
    end subroutine 

end module 
