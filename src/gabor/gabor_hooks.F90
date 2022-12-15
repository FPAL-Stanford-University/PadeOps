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
      character(len=*), intent(in), optional :: datadir

    end subroutine
  end interface

  interface getDomainBounds
    subroutine getDomainBounds()
  end interface
end module gaborHooks
