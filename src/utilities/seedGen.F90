module seedGen
    use kind_parameters, only: rkind

    implicit none 
contains 

function get_seed_from_location(xloc, yloc, zloc, offset) result(seed)  
    real(rkind), intent(in) :: xloc, yloc, zloc 
    real(rkind), intent(in), optional :: offset 
    real(rkind) :: location, addThis 
    integer :: seed
    
    addThis = 8.4241D14
    if (present(offset)) addThis = offset 

    location = xloc*1.314D18 + yloc*3.14341D18 + zloc*7.542351D18 + 3.4234D15 + addThis 
    seed = abs(transfer(location,1))

end function 

end module 
