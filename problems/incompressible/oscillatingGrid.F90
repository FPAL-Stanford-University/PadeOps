! Template for PadeOps

#include "oscillatingGrid_files/initialize.F90"       
#include "oscillatingGrid_files/temporalHook.F90"  

module gridFuncMod
    use kind_parameters, only: rkind 
    use IncompressibleGrid, only: igrid
    use constants 
    implicit none 

contains 

    subroutine setGridMask(xG, yG, zG, time, mask) 
        real(rkind), dimension(:,:,:), intent(in) :: xG, yG, zG
        real(rkind), dimension(:,:,:), intent(out) :: mask
        real(rkind), intent(in) :: time 
        integer :: i, j, k
        real(rkind) :: x, y, z

        ! TODO: Ryan - figure out what the grid function actually is. 
        do k = 1,size(mask,3) 
            z = zG(1,1,k)
            do j = 1,size(mask,2) 
                y = yG(1,j,1) 
                do i = 1,size(mask,1) 
                    x = xG(i, 1, 1)
                end do 
            end do 
        end do 

        mask = max(x*y*z*time, one)

    end subroutine 

    subroutine set_grid_position(igp, rkStage) 
        class(igrid), intent(inout) :: igp 
        integer :: rkStage
        real(rkind) :: tnow 

        select case (rkStage) 
        ! TO: Figure out substage times 
        case (1) 
            tnow = igp%tsim
        
        case default 
            tnow = igp%tsim 
        end select 

        call setGridMask(igp%mesh(:,:,:,1), igp%mesh(:,:,:,2), igp%mesh(:,:,:,3), tnow, & 
                igp%immersedBodies(1)%RbodyC) 
        call setGridMask(igp%mesh(:,:,:,1), igp%mesh(:,:,:,2), igp%zE           , tnow, & 
                igp%immersedBodies(1)%RbodyE) 

        igp%immersedBodies(1)%uTarget = zero 
        igp%immersedBodies(1)%vTarget = zero 
        igp%immersedBodies(1)%wTarget = zero ! TODO: Ryan Finish this  

    end subroutine 

end module 

program oscillatingGrid
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use gridFuncMod
    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    call igp%start_io(.false.)                !<-- Start I/O by creating a header file (see io.F90)
    
    call igp%printDivergence()

    if (igp%gpC%xsz(2) .ne. igp%gpE%xsz(2)) then 
        print*, "GP issue:" 
        print*, igp%gpC%xsz(2), igp%gpE%xsz(2)
        stop 
    end if 

    call tic() 
    do while (igp%tsim < igp%tstop) 

       call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igridWallM.F90)
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
       call set_grid_position(igp, 1)
       call igp%advance_SSP_RK45_Stage_1()
       
       call set_grid_position(igp, 2)
       call igp%advance_SSP_RK45_Stage_2()

       call set_grid_position(igp, 3)
       call igp%advance_SSP_RK45_Stage_3()

       call set_grid_position(igp, 4)
       call igp%advance_SSP_RK45_Stage_4()

       call set_grid_position(igp, 5)
       call igp%advance_SSP_RK45_Stage_5()
       
       call igp%wrapup_timestep()
       
       call doTemporalStuff(igp)                                        

    end do 
 
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
