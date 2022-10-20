module QHmeshRoutines
    use kind_parameters, only: rkind
    implicit none 

contains
    
    pure function getMeshEdgeValues(st,sz,h,xst) result(x)
      ! Define a 1D vector based on global indices and uniform grid spacing
      ! The global mesh begins at the domain boundary
      ! Inputs:
      !   st --> global start indicex of vector
      !   sz --> size of the vector
      !   h --> grid spacing
      !   xst --> domain boundary, i.e. x(1) = xst
      integer, intent(in) :: st, sz
      real(rkind), intent(in) :: h
      real(rkind), intent(in) :: xst
      real(rkind), dimension(sz) :: x
      integer :: i
         
      do i = 1,sz
        x(i) = xst + real(st + i - 2,rkind)*h
      end do
    end function

end module QHmeshRoutines
