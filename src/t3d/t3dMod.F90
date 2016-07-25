module t3dMod
    use kind_parameters, only: rkind 
    
    use mpi

    implicit none
    private
    public :: t3d

    type :: t3d
        private
        integer, public :: px, py, pz
        integer, public :: comm3d, commX, commY, commZ
        integer, public :: commXY, commYZ, commXZ
        integer, public :: rank3d, rankX, rankY, rankZ
        integer, public :: nprocs3d, nprocsX, nprocsY, nprocsZ
        integer, dimension(3) :: coords3d
        integer, dimension(2) :: coordsX, coordsY, coordsZ
        integer :: px_y, px_z
        integer :: py_x, py_z
        integer :: pz_y, pz_x
        integer, dimension(3), public :: sz3d, st3d, en3d 
        integer, dimension(3), public :: szX , stX , enX  
        integer, dimension(3), public :: szY , stY , enY  
        integer, dimension(3), public :: szZ , stZ , enZ  
    
        
    contains
        procedure :: transpose_3d_to_x
        procedure :: transpose_x_to_3d
        final :: destroy
    end type 
    
    interface t3d
        module procedure init
    end interface

contains 

    function init(comm3d, nx, ny, nz, px, py, pz) result(this)
        type(t3d) :: this
        integer, intent(in) :: comm3d, nx, ny, nz, px, py, pz
        integer :: ierr

        ! SAFEGUARD 

        if ((px == 0) .or. (py == 0) .or. (pz == 0)) then
            call GracefulExit("Incomplete. Optimization is not done.", 0)
        end if 
        call mpi_comm_size(comm3d, this%nprocs3d, ierr)
        call mpi_comm_rank(comm3d, this%rank3d, ierr)
         
        if (px*py*pz .ne. this%nprocs3d) then
            call GracefulExit("(PX x PY x PZ) must equal the number of processes in &
            & comm3d sent in.", 1)
        end if

    end function
   

    impure elemental subroutine destroy(this)
        type(t3d), intent(inout) :: this

    end subroutine

    subroutine transpose_3d_to_x(this, input, output)
        class(t3d), intent(in) :: this
        real(rkind), dimension(this%sz3d(1),this%sz3d(2),this%sz3d(3)), intent(in)  :: input
        real(rkind), dimension(this%szX (1),this%szX (2),this%szX (3)), intent(out) :: output

    end subroutine 

    subroutine transpose_x_to_3d(this, input, output)
        class(t3d), intent(in) :: this
        real(rkind), dimension(this%sz3d(1),this%sz3d(2),this%sz3d(3)), intent(out)  :: output
        real(rkind), dimension(this%szX (1),this%szX (2),this%szX (3)), intent(in)   :: input

    end subroutine 

end module 
