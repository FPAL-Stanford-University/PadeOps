module turbineMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa 
    use decomp_2d
    use StaggOpsMod, only: staggOps  
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use mpi 
    use reductions, only: p_maxval, p_sum
    use timer, only: tic, toc

    implicit none

    private
    public :: TurbineArray


    type :: TurbineArray
        integer :: myProc
        integer, dimension(:), allocatable  :: xst, xen, yst, yen
        integer :: nTurbines
        real(rkind), dimension(:), allocatable :: xLoc, yLoc, zLoc
        type(decomp_info), pointer :: gpC, sp_gpC, gpE, sp_gpE
        type(spectral), pointer :: spectC, spectE
        integer :: myLeftNeigh, myRightNeigh, myTopNeigh, myBotNeigh 
   

    contains
        procedure :: init
        procedure :: destroy
        procedure :: getForceRHS 
    end type
contains
subroutine init(this, inputFile, gpC, gpE, spectC, spectE)
    class(TurbineArray), intent(inout) :: this
    character(len=*), intent(in) :: inputFile
    type(spectral), target :: spectC, spectE
    type(decomp_info), target :: gpC, gpE

    print*, inputFile
    this%gpC => gpC
    this%spectC => this%spectC
    this%sp_gpC => this%spectC%spectdecomp

    this%gpE => gpE
    this%spectE => this%spectE
    this%sp_gpE => this%spectE%spectdecomp

    call GracefulExit("Wind Turbine stuff is incomplete", 423)

end subroutine


subroutine destroy(this)
    class(TurbineArray), intent(inout) :: this
    nullify(this%gpC, this%gpE, this%spectC, this%sp_gpC)
end subroutine


subroutine getForceRHS(this, u, v, wC, urhs, vrhs, wrhs)
    class(TurbineArray), intent(inout) :: this
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: u, v, wC
    complex(rkind), dimension(this%sp_gpC%xsz(1),this%sp_gpC%xsz(2),this%sp_gpC%xsz(3)), intent(inout) :: urhs, vrhs
    complex(rkind), dimension(this%sp_gpE%xsz(1),this%sp_gpE%xsz(2),this%sp_gpE%xsz(3)), intent(inout) :: wrhs 

    urhs = urhs + zero*u
    vrhs = vrhs + zero*v
    wrhs = wrhs + zero*wC


end subroutine 

end module
