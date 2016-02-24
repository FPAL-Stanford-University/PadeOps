module staggOpsMod
    use kind_parameters, only: rkind
    use exits, only: GracefulExit
    use decomp_2d
    use constants, only: half

    implicit none

    private

    public :: staggOps

    type :: staggOps
        private
        integer :: nxcell, nycell, nzcell
        integer :: nxedge, nyedge, nzedge
        type(decomp_info), pointer :: cellDecomp
        type(decomp_info), allocatable :: edgeDecomp
        integer :: stagg_scheme
        real(rkind) :: dx, dy, dz
        contains
            procedure :: init
            procedure :: destroy
            procedure :: InterpZ_Edge2Cell 
            
            !procedure :: InterpZ_Cell2Edge
            !procedure :: Stagg_ddz
    end type

contains

    subroutine init(this, gp, stagg_scheme , dx, dy, dz)
        class(staggOps), intent(inout) :: this
        class(decomp_info), intent(in), target:: gp
        integer, intent(in) :: stagg_scheme
        real(rkind), intent(in) :: dx, dy, dz

        this%nxcell = gp%zsz(1)
        this%nycell = gp%zsz(2)
        this%nzcell = gp%zsz(3)

        this%nzedge = this%nzcell + 1
        this%nxedge = this%nxcell 
        this%nyedge = this%nycell 

        this%stagg_scheme = stagg_scheme
        this%cellDecomp => gp
        call decomp_info_init(this%nxedge, this%nyedge, this%nzedge, this%edgedecomp)

        this%dx = dx
        this%dy = dy
        this%dz = dz

    end subroutine


    subroutine destroy(this)
        class(staggOps), intent(inout) :: this

        if (allocated(this%edgeDecomp)) deallocate(this%edgeDecomp)

    end subroutine

    pure subroutine InterpZ_Edge2Cell(this, edgeArr, cellArr)
        class(staggOps), intent(in) :: this
        complex(rkind), intent(in), dimension(this%nxedge, this%nyedge, this%nzedge) :: edgeArr
        complex(rkind), intent(out), dimension(this%nxcell, this%nycell, this%nzcell) :: cellArr
        integer :: k, j, i 

        cellArr = edgeArr(1:this%nxcell,1:this%nycell,1:this%nzcell)
        do k = 1,this%nzedge-1
            do j = 1,this%nycell
                do i = 1,this%nxcell
                    cellArr(i,j,k) = cellArr(i,j,k) + edgeArr(i,j,k+1)
                 end do 
            end do 
        end do
        cellArr = half*cellArr   

    end subroutine

end module 
