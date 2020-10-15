!#include "test_budget2_files/io.F90"
#include "test_budget2_files/initialize.F90"       

program test_budget2
    use mpi
    use kind_parameters, only: rkind, clen
    use IncompressibleGrid, only: igrid
    use budgets_time_avg_mod, only: budgets_time_avg  

    implicit none 

    
    type(igrid), allocatable, target :: adsim
    type(budgets_time_avg) :: budg_tavg
    character(len=clen) :: inputfile, AD_InputFile
    integer :: ierr

    call MPI_Init(ierr)                                                

    call GETARG(1,AD_Inputfile)

    allocate(adsim)                                                     
    call adsim%init(AD_InputFile, .true.)                               
    call budg_tavg%init(AD_Inputfile, adsim)   !<-- Budget class initialization 

    ! read budget0, budget1
    ! construct budget2_3 from budget0, budget1
    ! construct budget2_3 from its definition

    call budg_tavg%destroy()
    call adsim%finalize_io()
    call adsim%destroy()
    deallocate(adsim)

    call MPI_Finalize(ierr)

end program test_budget2
