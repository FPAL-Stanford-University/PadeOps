subroutine doDebugChecks(this)
  class(enrichmentOperator), intent(inout) :: this
  integer :: n
  real(rkind) :: Lx, Ly, Lz

  Lx = this%largeScales%dx*this%largeScales%nx
  Ly = this%largeScales%dy*this%largeScales%ny
  Lz = this%largeScales%dz*this%largeScales%nz

  do n = 1,this%nmodes
    ! Confirm the mode locations are within the domain bounds
    call assert(this%x(n) >= 0.d0, 'Mode x-location less than 0')
    call assert(this%x(n) <= Lx,   'Mode x-location greater than Lx')
    call assert(this%y(n) >= 0.d0, 'Mode y-location less than 0')
    call assert(this%y(n) <= Ly,   'Mode y-location greater than Ly')
    call assert(this%z(n) >= 0.d0, 'Mode z-location less than 0')
    call assert(this%z(n) <= Lz,   'Mode z-location greater than Lz')
  end do
end subroutine
