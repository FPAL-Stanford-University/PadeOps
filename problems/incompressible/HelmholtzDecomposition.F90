! Template for PadeOps

#include "HelmholtzDecomposition_files/initialize.F90"       
#include "HelmholtzDecomposition_files/temporalHook.F90"  

module doHelmholtzDecomposition
  use incompressibleGrid, only: igrid, wBC_bottom, wBC_top
  use kind_parameters,    only: rkind
  use reductions,         only: p_maxval
  use exits,              only: message
  use fortran_assert,     only: assert
  implicit none
  
contains
 
  subroutine getVectorPotential(igp,A)
    use HelmholtzDecomposition_parameters, only: divergence
    ! Solve for the scalar (phi) and vector potentials (A) which yield u =
    ! uphi + upsi = -grad(phi) + curl(A), div(A)=0. See Appendix 1 in Davidson'a "turbulence
    ! - an introduction for scientists and engineers (2d edition)"
  
    ! This routine assumes we have the all three velocities, u, v, and (w and wC)
    type(igrid), intent(inout), target :: igp
    real(rkind), dimension(:,:,:,:), intent(inout), target :: A
    real(rkind), parameter :: tol = 1.d-12
    real(rkind), dimension(:,:,:), pointer :: Ax, Ay, AzC
    real(rkind), dimension(:,:,:), pointer :: Az
    complex(rkind), dimension(:,:,:), pointer :: Axhat, Ayhat, Azhat
  
    ! Compute the divergence
    call igp%padepoiss%DivergenceCheck(igp%uhat, igp%vhat, igp%what, igp%divergence)
    
    ! Make sure the velocity field is divergence free
    if (p_maxval(maxval(abs(igp%divergence))) > tol) then
      call message('WARNING. Max divergence:',p_maxval(maxval(abs(igp%divergence))))
      call assert(p_maxval(maxval(abs(igp%divergence))) < tol,'Non-zero divergence')
    end if
  
    call igp%padepoiss%PoissonSolver_HomogeneousNeumannBCz(-igp%ox, A(:,:,:,1))
    call igp%padepoiss%PoissonSolver_HomogeneousNeumannBCz(-igp%oy, A(:,:,:,2))
    call igp%padepoiss%PoissonSolver_HomogeneousNeumannBCz(-igp%oz, A(:,:,:,3))

    ! Confirm A is solenoidal
    Axhat => igp%cbuffyC(:,:,:,1)
    Ayhat => igp%cbuffyC(:,:,:,2)
    Azhat => igp%cbuffyE(:,:,:,1)
    Az    => igp%rbuffxE(:,:,:,1)

    Ax    => A(:,:,:,1)
    Ay    => A(:,:,:,2)
    AzC   => A(:,:,:,3)

    call igp%spectC%fft(Ax, Axhat)
    call igp%spectC%fft(Ay, Ayhat)
    call igp%interpolate_cellField_to_edgeField(AzC, Az, -1, -1)
    call igp%spectE%fft(Az, Azhat)
    call igp%padepoiss%divergenceCheck(Axhat, Ayhat, Azhat, divergence)
    nullify(Axhat, Ayhat, Azhat, Az) 
    nullify(Ax, Ay, AzC) 
  end subroutine 
end module

program HelmholtzDecomposition
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use exits, only: message
    use HelmholtzDecomposition_parameters, only: tidst, tiden, tidStride, &
      A, finalizeProblem, allocateMemory
    use doHelmholtzDecomposition, only: getVectorPotential

    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr, tid

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)

    call allocateMemory(igp%gpC,igp%sp_gpC)

    ! Compute the vector potential
    call getVectorPotential(igp,A)

    ! Write the vector potential to disk
    call doTemporalStuff(igp)
 
    call igp%destroy()                !<-- Destroy the IGRID derived type

    deallocate(igp)
    
    call finalizeProblem()

    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
