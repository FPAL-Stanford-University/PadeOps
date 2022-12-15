module temporalHook
    use IncompressibleGrid, only: igrid
    implicit none 
contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout), target :: gp 

    end subroutine

end module 
