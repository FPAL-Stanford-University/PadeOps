module gaborHooks
  use kind_parameters, only: rkind
  use IncompressibleGrid, only: igrid
  implicit none

  interface getLargeScaleParams
    subroutine getLargeScaleParams(KE, L, LES, datadir)
      import :: rkind
      import :: igrid
      real(rkind), dimension(:,:,:), intent(out) :: KE, L
      type(igrid), intent(inout) :: LES
      ! Make the option to simply read the large scale parameters from disk
      character(len=*), intent(in) :: datadir

    end subroutine
  end interface

  interface getDomainBoundaries
    subroutine getDomainBoundaries(xDom,yDom,zDom,mesh)
      import ::  rkind
      real(rkind), dimension(2), intent(out) :: xDom, yDom, zDom
      real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    end subroutine
  end interface
end module gaborHooks
