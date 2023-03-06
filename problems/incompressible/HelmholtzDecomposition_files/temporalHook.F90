module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval, p_sum
    use exits,              only: message, message_min_max
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    use decomp_2d
    use HelmholtzDecomposition_parameters, only: A 

    implicit none 

    integer :: nt_print2screen = 1
    real(rkind) :: DomMaxDiv
    integer :: ierr 

contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout), target :: gp 
      
        call gp%dumpFullField(A(:,:,:,1),'vecx',gp%gpC)
        call gp%dumpFullField(A(:,:,:,2),'vecy',gp%gpC)
        call gp%dumpFullField(A(:,:,:,3),'vecz',gp%gpC)
      
        call message(0, "Just dumped potential fields")
        
    end subroutine

end module 
