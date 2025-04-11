module MultiBlockTopology

    use kind_parameters, only: rkind
    use constants,       only: half,one

    implicit none

    private
    public :: multiblocktopol

    type :: multiblocktopol

        private 
        integer :: x_num_blocks, y_num_blocks, z_num_blocks
        integer, allocatable, dimension(:,:) :: xst, yst, zst, xen, yen, zen
        
        contains

            procedure :: init
            procedure :: destroy

    end type

contains

    subroutine init(this)
        class(multiblocktopol), intent(inout) :: this
        integer :: ib

        this%x_num_blocks = 1
        this%y_num_blocks = 1
        this%z_num_blocks = 1

        allocate(this%xst(this%x_num_blocks, 3))
        allocate(this%xen(this%x_num_blocks, 3))
        allocate(this%yst(this%y_num_blocks, 3))
        allocate(this%yen(this%y_num_blocks, 3))
        allocate(this%zst(this%z_num_blocks, 3))
        allocate(this%zen(this%z_num_blocks, 3))

        ! only 1 block for now
        ! x-decomp
        do ib = 1, this%x_num_blocks
          this%xst(ib, 1) = 1;          this%xen(ib, 1) = 10
          this%xst(ib, 2) = 1;          this%xen(ib, 2) = 10
          this%xst(ib, 3) = 1;          this%xen(ib, 3) = 10
        end do

        ! y-decomp
        do ib = 1, this%y_num_blocks
          this%yst(ib, 1) = 1;          this%yen(ib, 1) = 10
          this%yst(ib, 2) = 1;          this%yen(ib, 2) = 10
          this%yst(ib, 3) = 1;          this%yen(ib, 3) = 10
        end do

        ! z-decomp
        do ib = 1, this%z_num_blocks
          this%zst(ib, 1) = 1;          this%zen(ib, 1) = 10
          this%zst(ib, 2) = 1;          this%zen(ib, 2) = 10
          this%zst(ib, 3) = 1;          this%zen(ib, 3) = 10
        end do
    end subroutine

    subroutine destroy(this)
        class(multiblocktopol), intent(inout) :: this

        deallocate(this%zen)
        deallocate(this%zst)
        deallocate(this%yen)
        deallocate(this%yst)
        deallocate(this%xen)
        deallocate(this%xst)

    end subroutine

end module
