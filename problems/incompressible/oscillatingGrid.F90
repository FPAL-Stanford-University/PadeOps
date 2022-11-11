! Template for PadeOps

#include "oscillatingGrid_files/initialize.F90"       
#include "oscillatingGrid_files/temporalHook.F90"  

module gridFuncMod
    use kind_parameters,    only: rkind 
    use IncompressibleGrid, only: igrid
    use igrid_Operators,    only: igrid_ops
    use constants
    use decomp_2d,          only: decomp_info
    use exits,              only: gracefulExit, message
    use reductions,         only: p_maxval
    use oscillating_grid_parameters, only: omega, stroke, lxbar, lybar, lzbar, z0, &
      Nbx, Nby, Lx, Ly, Lz, filterMask
    implicit none 

contains 

    subroutine setGridMask(x, y, z, time, mask, gp) 
        real(rkind), dimension(:,:,:), intent(in) :: x, y, z
        real(rkind), dimension(:,:,:), intent(out) :: mask
        real(rkind), intent(in) :: time
        class(decomp_info), intent(in) :: gp
        integer :: n, i, j, k
        integer :: ist, ien, jst, jen, kst, ken 
        real(rkind) :: xmin, xmax, ymin, ymax, zmin, zmax, dxgrid, dygrid
        real(rkind) :: dx, dy, dz

        dx = x(2,1,1) - x(1,1,1)
        dy = y(1,2,1) - y(1,1,1)
        dz = z(1,1,2) - z(1,1,1)

        ! dxgrid, dygrid: Size of grid gaps + bar width
        dxgrid = Lx/real(Nbx,rkind) - lxbar 
        dygrid = Ly/real(Nby,rkind) - lybar
        zmin = z0 - 0.5d0*lzbar + stroke*sin(omega*time)
        zmax = z0 + 0.5d0*lzbar + stroke*sin(omega*time)
        
        mask = zero

        kst = max(ceiling(zmin/dz),gp%xst(3))
        ken = min(floor(  zmax/dz),gp%xen(3))
          
        kst = kst - gp%xst(3) + 1
        ken = ken - gp%xst(3) + 1

        do n = 1,Nbx
          xmin = (n-1)*(dxgrid + lxbar) 
          xmax = xmin + lxbar - 1.d-14
          ist = max(ceiling(xmin/dx),gp%xst(1))
          ien = min(floor(  xmax/dx),gp%xen(1))

          ist = ist - gp%xst(1) + 1
          ien = ien - gp%xst(1) + 1

          do k = kst, ken
            do j = 1, gp%xsz(2)
              do i = ist, ien
                mask(i,j,k) = one 
              end do
            end do
          end do
        end do
        
        do n = 1,Nby
          ymin = (n-1)*(dygrid + lybar) 
          ymax = ymin + lybar - 1.d-14
          jst = max(ceiling(ymin/dy),gp%xst(2))
          jen = min(floor(  ymax/dy),gp%xen(2))
          
          jst = jst - gp%xst(2) + 1
          jen = jen - gp%xst(2) + 1
          
          do k = kst, ken
            do j = jst, jen
              mask(:,j,k) = one
            end do
          end do
        end do

    end subroutine 

    subroutine set_grid_position(igp, igpOp, rkStage)
        class(igrid), intent(inout) :: igp 
        class(igrid_ops), intent(inout) :: igpOp 
        integer, intent(in) :: rkStage
        real(rkind) :: tnow 
        integer :: ierr

        select case (rkStage) 
        ! TO: Figure out substage times 
        case (1) 
            tnow = igp%tsim
        case (2) 
            tnow = igp%tsim + 0.5d0*igp%dt
        case (3) 
            tnow = igp%tsim + 0.5d0*igp%dt
        case (4) 
            tnow = igp%tsim +       igp%dt
        case default 
            call gracefulExit('Invalid RK substep',ierr) 
        end select 

        call setGridMask(igp%mesh(:,:,:,1), igp%mesh(:,:,:,2), igp%mesh(:,:,:,3), &
          tnow, igp%immersedBodies(1)%RbodyC, igp%gpC) 

        if (filterMask) then
          igp%rbuffxC(:,:,:,1) = igp%immersedBodies(1)%RbodyC
          igp%rbuffxE(:,:,:,1) = igp%immersedBodies(1)%RbodyE
          
          call igpOp%filterField(igp%rbuffxC(:,:,:,1),&
            igp%immersedBodies(1)%RbodyC)
        end if
!DEBUG
if (rkStage == 1) call igp%dumpFullField(igp%immersedBodies(1)%RbodyC,'mask',igp%gpC)
!END DEBUG
        
        call igp%interpolate_cellField_to_edgeField(&
          igp%immersedBodies(1)%RbodyC,igp%immersedBodies(1)%RbodyE,0,0)

        igp%immersedBodies(1)%uTarget = zero 
        igp%immersedBodies(1)%vTarget = zero 
        igp%immersedBodies(1)%wTarget = stroke*omega*cos(omega*tnow) ! TODO: Ryan Finish this  

    end subroutine 

end module 

program oscillatingGrid
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
    use igrid_Operators, only: igrid_ops
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use gridFuncMod, only: set_grid_position
    use oscillating_grid_parameters, only: filterMask
    use fortran_assert, only: assert

    implicit none

    type(igrid), allocatable, target :: igp
    type(igrid_ops) :: igpOp
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

    call igpOp%init(igp%nx, igp%ny, igp%nz, igp%dx, igp%dy, igp%dz, &
      igp%inputdir, igp%outputdir, igp%runID, igp%periodicInZ, &
      igp%numericalSchemeVert)
    if (filterMask) then
      call assert(mod(igp%nx,2) == 0,'choose different domain or filter size')
      call assert(mod(igp%ny,2) == 0,'choose different domain or filter size')
      call igpOp%initFilter(igp%nx/2,igp%ny/2,1)
    endif

    call tic() 
    do while (igp%tsim < igp%tstop .and. igp%step < igp%nsteps) 

       call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igridWallM.F90)
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
       call set_grid_position(igp, igpOp, 1)
       call igp%advance_RK4_Stage(stage=1)
       
       call set_grid_position(igp, igpOp, 2)
       call igp%advance_RK4_Stage(stage=2)

       call set_grid_position(igp, igpOp, 3)
       call igp%advance_RK4_Stage(stage=3)

       call set_grid_position(igp, igpOp, 4)
       call igp%advance_RK4_Stage(stage=4)

       call igp%wrapup_timestep()
       
       call doTemporalStuff(igp)                                        

    end do 
 
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
