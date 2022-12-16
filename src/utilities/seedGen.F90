module seedGen
    use kind_parameters, only: rkind

    implicit none 
contains 

function get_seed_from_location(xloc, yloc, zloc, offset) result(seed)  
    real(rkind), intent(in) :: xloc, yloc, zloc 
    real(rkind), intent(in), optional :: offset 
    real(kind=4) :: location, addThis 
    integer :: seed
    
    addThis = 8.4241E32
    if (present(offset)) addThis = real(offset, kind=4)

    location = real(xloc,kind=4)*1.314D32 + real(yloc,kind=4)*3.14341D32 + real(zloc,kind=4)*7.542351D32 + real(3.4234D32,kind=4) + addThis 
    seed = abs(transfer(location,1))

end function 

end module 
