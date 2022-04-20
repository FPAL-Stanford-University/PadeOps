module GaborModPopMod
  use kind_parameters, only: rkind
  use GaborModeMod, only: GaborMode
  implicit none
  private
  public :: GaborModePop

  type :: GaborModePop
    class(GaborMode), dimension(:), allocatable :: modes

    contains
      procedure :: init
      procedure :: dest

  end type

contains 
  
  subroutine init(this)
    class(GaborModePop), intent(inout) :: this
  end subroutine

  subroutine dest(this)
    class(GaborModePop), intent(inout) :: this
  end subroutine

end module
