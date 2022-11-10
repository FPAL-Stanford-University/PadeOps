! Template for PadeOps

#include "oscillatingGrid_files/initialize.F90"       
#include "oscillatingGrid_files/temporalHook.F90"  

module gridFuncMod
    use kind_parameters,    only: rkind 
    use IncompressibleGrid, only: igrid
    use constants
    use decomp_2d,          only: decomp_info
    use exits,              only: gracefulExit, message
    use reductions,         only: p_maxval
    implicit none 

contains 

    subroutine setGridMask(x, y, z, Lx, Ly, Lz, time, omega, stroke, mask, &
        lxbar, lybar, lzbar, Nbx, Nby, z0, gp) 
        real(rkind), dimension(:,:,:), intent(in) :: x, y, z
        real(rkind), intent(in) :: Lx, Ly, Lz
        real(rkind), dimension(:,:,:), intent(out) :: mask
        real(rkind), intent(in) :: time
        real(rkind), intent(in) :: omega      ! Oscillating frequency of the grid
        real(rkind), intent(in) :: stroke     ! Amplitude of the grid (i.e. "stroke length")
        real(rkind), intent(in) :: lxbar, lybar, lzbar ! The width of each "bar"
        integer, intent(in) :: Nbx, Nby       ! The number of bars in x and y
        real(rkind), intent(in) :: z0         ! z-location of the center of the grid at t=0
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
          xmin = x(1,1,1) + (n-1)*(dxgrid + lxbar) 
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
          ymin = y(1,1,1) + (n-1)*(dygrid + lybar) 
          ymax = ymin + lybar - 1.d-14
          jst = max(ceiling(ymin/dy),gp%xst(2))
          jen = min(floor(  ymax/dy),gp%xen(2))
          
          jst = jst - gp%xst(2) + 1
          jen = jen - gp%xst(2) + 1
          
          do k = kst, ken
            do j = jst, jen
              do i = 1, gp%xsz(1)
                mask(i,j,k) = one
              end do
            end do
          end do
        end do

    end subroutine 

    subroutine set_grid_position(igp, rkStage, Lx, Ly, Lz, omega, stroke, &
        lxbar, lybar, lzbar, Nbx, Nby, z0) 
        class(igrid), intent(inout) :: igp 
        integer, intent(in) :: rkStage
        real(rkind), intent(in) :: Lx, Ly, Lz
        real(rkind), intent(in) :: omega      ! Oscillating frequency of the grid
        real(rkind), intent(in) :: stroke     ! Amplitude of the grid (i.e. "stroke length")
        real(rkind), intent(in) :: lxbar, lybar, lzbar ! The width of each "bar"
        integer, intent(in) :: Nbx, Nby       ! The number of bars in x and y
        real(rkind), intent(in) :: z0         ! z-location of the center of the grid at t=0
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
          Lx, Ly, Lz, tnow, omega, stroke, igp%immersedBodies(1)%RbodyC, lxbar, lybar, &
          lzbar, Nbx, Nby, z0, igp%gpC) 
        call setGridMask(igp%mesh(:,:,:,1), igp%mesh(:,:,:,2), igp%zE           , &
          Lx, Ly, Lz, tnow, omega, stroke, igp%immersedBodies(1)%RbodyE, lxbar, lybar, &
          lzbar, Nbx, Nby, z0, igp%gpE) 

        igp%immersedBodies(1)%uTarget = zero 
        igp%immersedBodies(1)%vTarget = zero 
        igp%immersedBodies(1)%wTarget = stroke*omega*cos(omega*tnow) ! TODO: Ryan Finish this  

    end subroutine 

end module 

program oscillatingGrid
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use gridFuncMod, only: set_grid_position
    use oscillating_grid_parameters, only: omega, stroke, lxbar, lybar, lzbar, z0, &
      Nbx, Nby, Lx, Ly, Lz

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
       
       call set_grid_position(igp, 1, Lx, Ly, Lz, omega, stroke, lxbar, lybar, lzbar, &
         Nbx, Nby, z0)
       call igp%advance_RK4_Stage(stage=1)
       
       call set_grid_position(igp, 2, Lx, Ly, Lz, omega, stroke, lxbar, lybar, lzbar, &
         Nbx, Nby, z0)
       call igp%advance_RK4_Stage(stage=2)

       call set_grid_position(igp, 3, Lx, Ly, Lz, omega, stroke, lxbar, lybar, lzbar, &
         Nbx, Nby, z0)
       call igp%advance_RK4_Stage(stage=3)

       call set_grid_position(igp, 4, Lx, Ly, Lz, omega, stroke, lxbar, lybar, lzbar, &
         Nbx, Nby, z0)
       call igp%advance_RK4_Stage(stage=4)

       call igp%wrapup_timestep()
       
       call doTemporalStuff(igp)                                        

    end do 
 
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
