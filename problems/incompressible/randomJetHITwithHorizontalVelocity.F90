! Template for PadeOps

#include "randomJetHIT_files/initialize.F90"       
#include "randomJetHIT_files/temporalHook.F90"  
#include "randomJetHIT_files/jetMod.F90"  

module gridFuncMod
    use kind_parameters,    only: rkind 
    use IncompressibleGrid, only: igrid
    use constants
    use decomp_2d,          only: decomp_info, transpose_x_to_y, &
      transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use exits,              only: gracefulExit, message
    use reductions,         only: p_maxval, p_minval
    use gaussianstuff,      only: gaussian
    use randomJetHIT_parameters, only: filterMaskHowManyTimes, Lx, Ly
    use fortran_assert,     only: assert
    use jetMod,             only: jet 
    implicit none 
    
    type(gaussian) :: gauss_x, gauss_y, gauss_zC, gauss_zE

contains 

    subroutine initializeJetArray(jetArray, NjetsX, NjetsY, jetXsz, jetYsz, &
        jetZsz, jetZst, inputfile, tInit)
      type(jet), dimension(:,:), allocatable, intent(inout) :: jetArray
      integer, intent(in) :: NjetsX, NjetsY
      real(rkind), intent(in) :: jetXsz, jetYsz, jetZsz, jetZst, tInit
      character(len=*), intent(in) :: inputfile
      real(rkind) :: dxJet, dyJet, yloc, xloc
      integer :: i, j, gID, NjetsTotal

      call assert(NjetsX > 0 .and. NjetsY > 0,'Must use at least one jet')
      allocate(jetArray(NjetsX,NjetsY))
    
      dxJet = Lx/real(NjetsX,rkind)
      dyJet = Ly/real(NjetsY,rkind)
    
      NjetsTotal = NjetsY*NjetsX

      do j = 1,NjetsY
        yloc = 0.5d0*dyJet + real(j - 1,rkind)*dyJet
        do i = 1,NjetsX
          xloc = 0.5d0*dxJet + real(i - 1,rkind)*dxJet
          gID = (j-1)*NjetsX + i
          call jetArray(i,j)%init(xloc,yloc,jetZst,jetXsz,jetYsz,jetZsz,gID,&
            inputfile,NjetsTotal,tInit)
        end do
      end do
    end subroutine

    subroutine getStEnIndices(xmin,xmax,dx,gpxst,gpxen,ist,ien)
      real(rkind), intent(in) :: xmin, xmax, dx
      integer, intent(in) :: gpxst, gpxen
      integer, intent(out) :: ist, ien
      ist = max(ceiling(xmin/dx),gpxst)
      ien = min(floor(  xmax/dx),gpxen)

      ist = ist - gpxst + 1
      ien = ien - gpxst + 1
    end subroutine

    subroutine setMask(ist,ien,jst,jen,kst,ken,gp,mask,maskVal)
      integer, intent(in) :: ist, ien, jst, jen, kst, ken
      class(decomp_info), intent(in) :: gp
      real(rkind), dimension(:,:,:), intent(inout) :: mask
      real(rkind), intent(in) :: maskVal
      integer :: i, j, k

      if (ist == 1 .and. ien == gp%xsz(1)) then
        do k = kst, ken
          do j = jst, jen
            mask(:,j,k) = mask(i,j,k) + maskVal
          end do
        end do
      else
        do k = kst, ken
          do j = jst, jen
            do i = ist, ien
              mask(i,j,k) = mask(i,j,k) + maskVal
            end do
          end do
        end do
      end if
    end subroutine

    subroutine setJetMask(x, y, z, time, arr, gp, thisJet, arrID)
        real(rkind), dimension(:,:,:), intent(in) :: x, y, z
        real(rkind), dimension(:,:,:), intent(out) :: arr
        real(rkind), intent(in) :: time
        integer, intent(in) :: arrID
        integer, parameter :: u = 1, v = 2, w = 3, mask = 4
        class(decomp_info), intent(in) :: gp
        class(jet), intent(inout) :: thisJet
        integer :: ist, ien, jst, jen, kst, ken 
        real(rkind) :: xmin, xmax, ymin, ymax, zmin, zmax
        real(rkind) :: dx, dy, dz
        
        if (thisJet%isOn(time)) then
          dx = x(2,1,1) - x(1,1,1)
          dy = y(1,2,1) - y(1,1,1)
          dz = z(1,1,2) - z(1,1,1)

          zmin = thisJet%zst 
          zmax = zmin + thisJet%zsz 
          
          kst = max(ceiling(zmin/dz),gp%xst(3))
          ken = min(floor(  zmax/dz),gp%xen(3))
            
          kst = kst - gp%xst(3) + 1
          ken = ken - gp%xst(3) + 1
          
          xmin = thisJet%xloc - 0.5d0*thisJet%xsz 
          xmax = thisJet%xloc + 0.5d0*thisJet%xsz
          ymin = thisJet%yloc - 0.5d0*thisJet%ysz 
          ymax = thisJet%yloc + 0.5d0*thisJet%ysz

          call getStEnIndices(xmin,xmax,dx,gp%xst(1),gp%xen(1),ist,ien)
          call getStEnIndices(ymin,ymax,dy,gp%xst(2),gp%xen(2),jst,jen)
          select case (arrID)
          case (u)
            call setMask(ist,ien,jst,jen,kst,ken,gp,arr,thisJet%weight*thisJet%u)
          case (v)
            call setMask(ist,ien,jst,jen,kst,ken,gp,arr,thisJet%weight*thisJet%v)
          case (w)
            call setMask(ist,ien,jst,jen,kst,ken,gp,arr,thisJet%weight*thisJet%w)
          case (mask)
            call setMask(ist,ien,jst,jen,kst,ken,gp,arr,1.d0)
          end select
        end if

    end subroutine

    subroutine filterField(field, gauss_x, gauss_y, gauss_z, gp, buffx, &
        buffy, buffz)
      real(rkind), dimension(:,:,:), intent(inout) :: field
      real(rkind), dimension(:,:,:,:), intent(inout) :: buffx, buffy, buffz
      class(gaussian), intent(inout) :: gauss_x, gauss_y, gauss_z
      class(decomp_info), intent(in) :: gp
      
      buffx(:,:,:,1) = field
      call gauss_x%filter1(buffx(:,:,:,1),buffx(:,:,:,2), gp%xsz(2), gp%xsz(3))
      call transpose_x_to_y(buffx(:,:,:,2),buffy(:,:,:,1),gp)
      call gauss_y%filter2(buffy(:,:,:,1),buffy(:,:,:,2), gp%ysz(1), gp%ysz(3))
      call transpose_y_to_z(buffy(:,:,:,2),buffz(:,:,:,1),gp)
      call gauss_z%filter3(buffz(:,:,:,1),buffz(:,:,:,2), gp%zsz(1), gp%zsz(2))
      call transpose_z_to_y(buffz(:,:,:,2),buffy(:,:,:,1),gp)
      call transpose_y_to_x(buffy(:,:,:,1),field,gp)
    end subroutine

    subroutine set_grid_position(igp, rkStage, jetArray)
        use randomJetHIT_parameters, only: dumpMaskFreq, NjetsX, NjetsY
        use jetMod, only: wJet
        class(igrid), intent(inout) :: igp 
        integer, intent(in) :: rkStage
        class(jet), dimension(:,:), intent(inout) :: jetArray
        real(rkind) :: tnow 
        integer :: i, j, ierr

        select case (rkStage) 
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

        igp%immersedBodies(1)%RbodyC  = zero 
        igp%immersedBodies(1)%RbodyE  = zero 
        igp%immersedBodies(1)%uTarget = zero 
        igp%immersedBodies(1)%vTarget = zero 
        igp%immersedBodies(1)%wTarget = zero

        do j = 1,NjetsY
          do i = 1,NjetsX
            call setJetMask(igp%mesh(:,:,:,1), igp%mesh(:,:,:,2), igp%mesh(:,:,:,3), &
              tnow, igp%immersedBodies(1)%RbodyC, igp%gpC, jetArray(i,j), 4)
            call setJetMask(igp%mesh(:,:,:,1), igp%mesh(:,:,:,2), igp%mesh(:,:,:,3), &
              tnow, igp%immersedBodies(1)%RbodyE, igp%gpE, jetArray(i,j), 4)
            
            call setJetMask(igp%mesh(:,:,:,1), igp%mesh(:,:,:,2), igp%mesh(:,:,:,3), &
              tnow, igp%immersedBodies(1)%uTarget, igp%gpC, jetArray(i,j), 1)
            call setJetMask(igp%mesh(:,:,:,1), igp%mesh(:,:,:,2), igp%mesh(:,:,:,3), &
              tnow, igp%immersedBodies(1)%vTarget, igp%gpC, jetArray(i,j), 2)
            call setJetMask(igp%mesh(:,:,:,1), igp%mesh(:,:,:,2), igp%mesh(:,:,:,3), &
              tnow, igp%immersedBodies(1)%wTarget, igp%gpE, jetArray(i,j), 3)
          end do
        end do
       
        do i = 1,filterMaskHowManyTimes
          call filterField(igp%immersedBodies(1)%RbodyC,gauss_x,gauss_y,gauss_zC,&
            igp%gpC,igp%rbuffxC(:,:,:,1:2), igp%rbuffyC(:,:,:,1:2), &
            igp%rbuffzC(:,:,:,1:2))
          call filterField(igp%immersedBodies(1)%RbodyE,gauss_x,gauss_y,gauss_zE,&
            igp%gpE,igp%rbuffxE(:,:,:,1:2), igp%rbuffyE(:,:,:,1:2), &
            igp%rbuffzE(:,:,:,1:2))
          
          call filterField(igp%immersedBodies(1)%uTarget,gauss_x,gauss_y,gauss_zC,&
            igp%gpC,igp%rbuffxC(:,:,:,1:2), igp%rbuffyC(:,:,:,1:2), &
            igp%rbuffzC(:,:,:,1:2))
          call filterField(igp%immersedBodies(1)%vTarget,gauss_x,gauss_y,gauss_zC,&
            igp%gpC,igp%rbuffxC(:,:,:,1:2), igp%rbuffyC(:,:,:,1:2), &
            igp%rbuffzC(:,:,:,1:2))
          call filterField(igp%immersedBodies(1)%wTarget,gauss_x,gauss_y,gauss_zE,&
            igp%gpE,igp%rbuffxE(:,:,:,1:2), igp%rbuffyE(:,:,:,1:2), &
            igp%rbuffzE(:,:,:,1:2))
        end do

        where (igp%immersedBodies(1)%RbodyC > 0.d0) igp%immersedBodies(1)%RbodyC = 1.d0
        where (igp%immersedBodies(1)%RbodyE > 0.d0) igp%immersedBodies(1)%RbodyE = 1.d0
        
        if (rkStage == 1 .and. mod(igp%step,dumpMaskFreq) == 0) then
          call igp%dumpFullField(igp%immersedBodies(1)%uTarget,'utgt',igp%gpC)
          call igp%dumpFullField(igp%immersedBodies(1)%vTarget,'vtgt',igp%gpC)
          call igp%dumpFullField(igp%immersedBodies(1)%wTarget,'wtgt',igp%gpE)
          call igp%dumpFullField(igp%immersedBodies(1)%RbodyC,'mask',igp%gpC)
        end if


    end subroutine 

end module 

program randomJetHITwithHoriontalVelocity
    use mpi
    use kind_parameters,  only: clen, rkind
    use IncompressibleGrid, only: igrid
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use gridFuncMod, only: set_grid_position, gauss_x, gauss_y, gauss_zC, &
      gauss_zE, initializeJetArray
    use randomJetHIT_parameters, only: filterMaskHowManyTimes, &
      Lx, Ly, Lz, NjetsX, NjetsY, jetXsz, jetYsz, jetZsz, jetZst
    use fortran_assert, only: assert
    use jetMod, only: jet

    implicit none

    type(igrid), allocatable, target :: igp
    type(jet), dimension(:,:), allocatable :: jetArray
    character(len=clen) :: inputfile
    integer :: ierr
    !integer :: i, j, gID
    !real(rkind) :: xloc, yloc, dxJet, dyJet

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

    if (filterMaskHowManyTimes > 0) then
      ierr = gauss_x%init( igp%nx,  .true.)
      call assert(ierr == 0,'gauss_x initialization')
      ierr = gauss_y%init( igp%ny,  .true.)
      call assert(ierr == 0,'gauss_y initialization')
      ierr = gauss_zC%init(igp%nz,  .false.)
      call assert(ierr == 0,'gauss_zC initialization')
      ierr = gauss_zE%init(igp%nz+1,.false.)
      call assert(ierr == 0,'gauss_zE initialization')
    endif

    call initializeJetArray(jetArray, NjetsX, NjetsY, jetXsz, jetYsz, jetZsz, jetZst, &
      inputfile, igp%tsim)

    call tic() 
    do while (igp%tsim < igp%tstop .and. igp%step < igp%nsteps) 

       call set_grid_position(igp, 1, jetArray)
       call igp%advance_RK4_Stage(stage=1)
       
       call set_grid_position(igp, 2, jetArray)
       call igp%advance_RK4_Stage(stage=2)

       call set_grid_position(igp, 3, jetArray)
       call igp%advance_RK4_Stage(stage=3)

       call set_grid_position(igp, 4, jetArray)
       call igp%advance_RK4_Stage(stage=4)
       
       call igp%wrapup_timestep()
       
       call doTemporalStuff(igp)                                        

    end do 
 
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults

    !call finalizeProblem(jetArray)

    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
