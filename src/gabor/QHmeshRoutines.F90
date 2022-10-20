module QHmeshRoutines
    use kind_parameters, only: rkind
    implicit none 

contains
    
    pure function getMeshEdgeValues(st,en,h) result(x)
      ! Define a 1D vector based on global indices and uniform grid spacing
      ! The global mesh begins at the domain boundary
      ! Inputs:
      !   st, en --> start and ending indices of vector
      !   h --> grid spacing
      integer, intent(in) :: st, en
      real(rkind), intent(in) :: h
      real(rkind), dimension(st:en) :: x
      integer :: i
         
      do i = st,en
        x(i) = (i - 1)*h
      end do
    end function

end module QHmeshRoutines