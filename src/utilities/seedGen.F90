module seedGen
    use kind_parameters, only: rkind
    use constants, only: rmaxInt
    use exits, only: warning
    use fortran_assert, only: assert
   
    implicit none 
    
    real(rkind), parameter :: initFact = 8.145d0
    real(rkind), parameter :: logMapFact = 4.d0
contains 

function get_seed_from_location(xloc, yloc, zloc, offset, globalIdx, seedNumber) result(seed)  
    real(rkind), intent(in) :: xloc, yloc, zloc 
    integer, intent(in), optional :: globalIdx, seedNumber
    real(rkind), intent(in), optional :: offset 
    real(kind=4) :: location, addThis
    real(rkind) :: newSeed
    integer :: seed
    
    addThis = 8.4241E32
    if (present(offset)) addThis = real(offset, kind=4)

    location = real(xloc,kind=4)*1.314D32 + real(yloc,kind=4)*3.14341D32 + real(zloc,kind=4)*7.542351D32 + real(3.4234D32,kind=4) + addThis 
    seed = abs(transfer(location,1))

    if (present(globalIdx)) then
      !num = 1
      !if (present(seedNumber)) num = seedNumber
      !if (real(num,rkind) < initFact) then
      !  oldSeed = real(num,rkind)/initFact
      !  do n = 1,globalIdx
      !    call incrementLogisticMap(oldSeed,newSeed)
      !    oldSeed = newSeed
      !  end do
      !else
      !  call warning('Seed issue!')
      !end if

      if (present(seedNumber)) then
        newSeed = initializeLogMap(globalIdx,seedNumber=seedNumber)
      else
        newSeed = initializeLogMap(globalIdx)
      end if

      seed = nint((0.5d0*newSeed + 0.25d0)*rmaxInt)
    end if
end function

function initializeLogMap(globalIdx,seedNumber,x0) result(logMapOut)
  integer, intent(in) :: globalIdx
  integer, intent(in), optional :: seedNumber
  real(rkind), intent(in), optional :: x0
  integer :: num, n
  real(rkind) :: logMapOut, old, new, x0use
  
  num = 1
  if (present(seedNumber)) num = seedNumber
  x0use = real(num,rkind)/initFact ! Initial condition for logistic map
  
  if (present(x0)) x0use = x0      ! User-prescribed initial condition

  if (x0use < 1.d0) then
    old = x0use
    do n = 1,globalIdx
      call incrementLogisticMap(old,new)
      old = new
    end do
  else
    call assert(.false.,'Seed issue!')
  end if
  logMapOut = new
end function

subroutine incrementLogisticMap(old,new)
    real(rkind), intent(in) :: old
    real(rkind), intent(out) :: new
            
    new = logMapFact*old*(1.d0 - old)

end subroutine

end module 
