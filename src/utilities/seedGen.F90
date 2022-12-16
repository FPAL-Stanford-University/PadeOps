module seedGen
    use kind_parameters, only: rkind
    use constants, only: rmaxInt
    use exits, only: warning
    implicit none 
contains 

function get_seed_from_location(xloc, yloc, zloc, offset, globalIdx, seedNumber) result(seed)  
    real(rkind), intent(in) :: xloc, yloc, zloc 
    integer, intent(in), optional :: globalIdx, seedNumber
    real(rkind), intent(in), optional :: offset 
    real(kind=4) :: location, addThis
    real(rkind) :: newSeed, oldSeed
    integer :: seed, n, num
    real(rkind), parameter :: initFact = 8.145d0
    real(rkind), parameter :: logMapFact = 4.d0
    
    addThis = 8.4241E32
    if (present(offset)) addThis = real(offset, kind=4)

    location = real(xloc,kind=4)*1.314D32 + real(yloc,kind=4)*3.14341D32 + real(zloc,kind=4)*7.542351D32 + real(3.4234D32,kind=4) + addThis 
    seed = abs(transfer(location,1))

    if (present(globalIdx)) then
      num = 1
      if (present(seedNumber)) num = seedNumber
      if (real(num,rkind) < initFact) then
        oldSeed = real(seedNumber,rkind)/initFact
        do n = 1,globalIdx
          newSeed = logMapFact*oldSeed*(1.d0 - oldSeed)
          oldSeed = newSeed
        end do
      else
        call warning('Seed issue!')
      end if

      seed = nint((0.5d0*newSeed + 0.25d0)*rmaxInt)
    end if
end function 

end module 
