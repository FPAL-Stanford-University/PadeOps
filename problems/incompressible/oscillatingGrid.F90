! Template for PadeOps

#include "oscillatingGrid_files/initialize.F90"       
#include "oscillatingGrid_files/temporalHook.F90"  

module gridFuncMod
<<<<<<< HEAD
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
=======
    use kind_parameters,    only: rkind 
    use IncompressibleGrid, only: igrid
    use igrid_Operators,    only: igrid_ops
    use constants
    use decomp_2d,          only: decomp_info, transpose_x_to_y, &
      transpose_y_to_z, transpose_z_to_y, transpose_y_to_x
    use exits,              only: gracefulExit, message
    use reductions,         only: p_maxval, p_minval
    use gaussianstuff, only: gaussian
    use oscillating_grid_parameters, only: widthFact, lengthFact, filterMask, &
      omega, stroke, NbarsTotal
    use fortran_assert,     only: assert 
    implicit none 
    
    type(gaussian) :: gauss_x, gauss_y, gauss_zC, gauss_zE

    interface addBars
      module procedure addBarsOnly, addBarsAndGetInfo 
    end interface 

contains 

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
            mask(:,j,k) = maskVal
          end do
        end do
      else
        do k = kst, ken
          do j = jst, jen
            do i = ist, ien
              mask(i,j,k) = maskVal
            end do
          end do
        end do
      end if
    end subroutine

    subroutine setGridMask(x, y, z, time, mask, gp) 
      use oscillating_grid_parameters, only: lxbar, lybar, lzbar, z0, &
        Nbx, Nby, Lx, Ly, Lz, gridType, width1, length1, &
        levels, useGridBreaks, phi1, phi2, omega1, omega2, &
        breakSzY, breakSzX, side, fractal, homogeneous
        real(rkind), dimension(:,:,:), intent(in) :: x, y, z
        real(rkind), dimension(:,:,:), intent(out) :: mask
        real(rkind), intent(in) :: time
        class(decomp_info), intent(in) :: gp
        integer :: n, ierr 
        integer :: ist, ien, jst, jen, kst, ken 
        real(rkind) :: xmin, xmax, ymin, ymax, zmin, zmax, dxgrid, dygrid
        real(rkind) :: dx, dy, dz
        real(rkind), dimension(2) :: center
        real(rkind) :: y1
        !integer, dimension(NbarsTotal,4) :: barStEn
        !real(rkind), dimension(NbarsTotal) :: barWidths
        real(rkind),dimension(NbarsTotal,5) :: barInfo
        integer, dimension(NbarsTotal) :: myCorner
        
        y1 = p_minval(minval(y))
        dx = x(2,1,1) - x(1,1,1)
        dy = y(1,2,1) - y(1,1,1)
        dz = z(1,1,2) - z(1,1,1)

        zmin = z0 - 0.5d0*lzbar + 0.5d0*stroke*sin(omega*time)
        zmax = z0 + 0.5d0*lzbar + 0.5d0*stroke*sin(omega*time)
        
        mask = zero

        kst = max(ceiling(zmin/dz),gp%xst(3))
        ken = min(floor(  zmax/dz),gp%xen(3))
          
        kst = kst - gp%xst(3) + 1
        ken = ken - gp%xst(3) + 1
        
        select case (gridType)
        case (homogeneous)
          ! dxgrid, dygrid: Size of grid gaps + bar width
          dxgrid = Lx/real(Nbx,rkind) - lxbar 
          dygrid = Ly/real(Nby,rkind) - lybar

          do n = 1,Nbx
            xmin = x(1,1,1) + (n-1)*(dxgrid + lxbar) 
            xmax = xmin + lxbar
            call getStEnIndices(xmin,xmax,dx,gp%xst(1),gp%xen(1),ist,ien)
            call setMask(ist,ien,1,gp%xsz(2),kst,ken,gp,mask,one)
          end do
          
          do n = 1,Nby
            ymin = y1 + (n-1)*(dygrid + lybar) 
            ymax = ymin + lybar
            call getStEnIndices(ymin,ymax,dy,gp%xst(2),gp%xen(2),jst,jen)
            call setMask(1,gp%xsz(1),jst,jen,kst,ken,gp,mask,one)
          end do
          if (useGridBreaks) then
            do n = 1,Nbx
              ymin = (0.5d0*Ly-breakSzY)*(1.d0 + sin(omega1(n)*time + phi1(n)))
              ymax = ymin + breakSzY
              xmin = x(1,1,1) + (n-1)*(dxgrid + lxbar)
              xmax = xmin + lxbar
              call getStEnIndices(xmin,xmax,dx,gp%xst(1),gp%xen(1),ist,ien)
              call getStEnIndices(ymin,ymax,dy,gp%xst(2),gp%xen(2),jst,jen)
              call setMask(ist,ien,jst,jen,kst,ken,gp,mask,zero)
            end do
            do n = 1,Nby
              xmin = (0.5d0*Lx-breakSzX)*(1.d0 + sin(omega2(n)*time + phi2(n)))
              xmax = xmin + breakSzX
              ymin = y1 + (n-1)*(dygrid + lybar)
              ymax = ymin + lybar
              call getStEnIndices(xmin,xmax,dx,gp%xst(1),gp%xen(1),ist,ien)
              call getStEnIndices(ymin,ymax,dy,gp%xst(2),gp%xen(2),jst,jen)
              call setMask(ist,ien,jst,jen,kst,ken,gp,mask,zero)
            end do
          end if
        case (fractal)
          center = [0.5d0*Lx, 0.5d0*Ly]
          if (useGridBreaks) then
            call addBars(1,levels,width1,length1,center,dx,dy,kst,ken,gp,mask,barInfo,myCorner)
            call addBreaks(mask,barInfo,dx,dy,kst,ken,time,omega1,phi1,side,breakSzX,breakSzY,gp)
          else
            call addBars(1,levels,width1,length1,center,dx,dy,kst,ken,gp,mask)
          end if
        case default
        end select

    end subroutine

    subroutine addBreaks(mask,barInfo,dx,dy,kst,ken,time,omega,phi,side,breakSzX,breakSzY,gp)
      real(rkind), dimension(:,:,:), intent(inout) :: mask
      real(rkind), dimension(:,:), intent(in) :: barInfo
      real(rkind), dimension(:), intent(in) :: omega, phi 
      real(rkind), intent(in) :: dx, dy, time, breakSzX, breakSzY
      integer, intent(in) :: kst, ken
      integer, dimension(:), intent(in) :: side
      class(decomp_info), intent(in) :: gp
      integer :: n, ist, ien, jst, jen
      real(rkind) :: xmin, xmax, ymin, ymax, xbreakMin, xbreakMax, ybreakMin, &
        ybreakMax
      real(rkind) :: width, A, start

      do n = 1,size(barInfo,1)
        xmin = barInfo(n,1)
        xmax = barInfo(n,2)
        ymin = barInfo(n,3)
        ymax = barInfo(n,4)
        width = barInfo(n,5)
        select case (side(n))
          case (1)
            A = 0.5d0*(ymax - ymin)
            start = ymin + A
            ybreakMin = start + A*sin(omega(n)*time + phi(n))
            ybreakMax = ybreakMin + breakSzY
            xbreakMin = xmin
            xbreakMax = xmin + width
          case (2) 
            A = 0.5d0*(xmax - xmin)
            start = xmin + A
            xbreakMin = start + A*sin(omega(n)*time + phi(n))
            xbreakMax = xbreakMin + breakSzX
            ybreakMax = ymax
            ybreakMin = ymax - width
          case (3)
            A = 0.5d0*(ymax - ymin)
            start = ymin + A
            ybreakMin = start + A*sin(omega(n)*time + phi(n))
            ybreakMax = ybreakMin + breakSzY
            xbreakMax = xmax
            xbreakMin = xmax - width
          case (4)
            A = 0.5d0*(xmax - xmin)
            start = xmin + A
            xbreakMin = start + A*sin(omega(n)*time + phi(n))
            xbreakMax = xbreakMin + breakSzX
            ybreakMin = ymin
            ybreakMax = ymin + width
          case default
            call assert(.false.,'Invalide side')
        end select
        call getStEnIndices(xbreakMin,xbreakMax,dx,gp%xst(1),gp%xen(1),ist,ien)
        call getStEnIndices(ybreakMin,ybreakMax,dy,gp%xst(2),gp%xen(2),jst,jen)
        call setMask(ist,ien,jst,jen,kst,ken,gp,mask,zero)
      end do
    end subroutine
    
    recursive subroutine addBarsAndGetInfo(level,nLevels,D,L,center,dx,dy,kst,ken,gp,mask,barInfo,myCorner)
      integer, intent(in) :: level, nLevels, kst, ken
      real(rkind), intent(in) :: D, L, dx, dy
      real(rkind), dimension(2), intent(in) :: center
      type(decomp_info), intent(in) :: gp
      real(rkind), dimension(:,:,:), intent(inout) :: mask
      real(rkind), dimension(:,:), intent(inout) :: barInfo
      integer, dimension(:), intent(inout) :: myCorner
      integer :: i, j, k, ist, ien, jst, jen, gID
      real(rkind) :: xmin, xmax, ymin, ymax, xminStore, xmaxStore, yminStore, ymaxStore
      real(rkind), dimension(2) :: corner

      xminStore = 1.d9
      xmaxStore = 1.d-9
      yminStore = 1.d9
      ymaxStore = 1.d-9

      if (level == 1) then
        myCorner = 0
        myCorner(1) = 1
      end if
      
      gID = getGlobalID(level,myCorner)
      call assert(gID > 0,'gID > 0')
      call assert(gId < NbarsTotal + 1,'gID < NbarsTotal + 1')

      ! Left side of square
      xmin = center(1) - 0.5d0*L
      xmax = xmin + D
      ymin = center(2) - 0.5d0*L
      ymax = center(2) + 0.5d0*L
#include "oscillatingGrid_files/gridFunc_files/setMaskCommon.F90"
#include "oscillatingGrid_files/gridFunc_files/updateMinMax.F90"
      
      ! Right side of square
      xmax = center(1) + 0.5d0*L
      xmin = xmax - D
      ymin = center(2) - L/2
      ymax = center(2) + L/2
#include "oscillatingGrid_files/gridFunc_files/setMaskCommon.F90"
#include "oscillatingGrid_files/gridFunc_files/updateMinMax.F90"
      
      ! Bottom of square
      xmin = center(1) - 0.5d0*L + D
      xmax = center(1) + 0.5d0*L - D
      ymin = center(2) - 0.5d0*L
      ymax = ymin + D
#include "oscillatingGrid_files/gridFunc_files/setMaskCommon.F90"
#include "oscillatingGrid_files/gridFunc_files/updateMinMax.F90"
      
      ! Top of square
      xmin = center(1) - 0.5d0*L + D
      xmax = center(1) + 0.5d0*L - D
      ymax = center(2) + 0.5d0*L
      ymin = ymax - D
#include "oscillatingGrid_files/gridFunc_files/setMaskCommon.F90"
#include "oscillatingGrid_files/gridFunc_files/updateMinMax.F90"
      
      call updateBarInfo(barInfo,gID,D,xminStore,yminStore,xmaxStore,ymaxStore,dx,dy,gp)

      if (level < nLevels) then
        call getBounds(center,D,L,'NW',corner)
        myCorner(level+1) = 1
        call addBars(level+1,nLevels,widthFact*D,lengthFact*L,corner,dx,dy,kst,ken,gp,mask,barInfo,myCorner)
        
        call getBounds(center,D,L,'NE',corner)
        myCorner(level+1) = 2
        call addBars(level+1,nLevels,widthFact*D,lengthFact*L,corner,dx,dy,kst,ken,gp,mask,barInfo,myCorner)
        
        call getBounds(center,D,L,'SE',corner)
        myCorner(level+1) = 3
        call addBars(level+1,nLevels,widthFact*D,lengthFact*L,corner,dx,dy,kst,ken,gp,mask,barInfo,myCorner)
        
        call getBounds(center,D,L,'SW',corner)
        myCorner(level+1) = 4
        call addBars(level+1,nLevels,widthFact*D,lengthFact*L,corner,dx,dy,kst,ken,gp,mask,barInfo,myCorner)
      end if 
    end subroutine
    
    recursive subroutine addBarsOnly(level,nLevels,D,L,center,dx,dy,kst,ken,gp,mask)
      integer, intent(in) :: level, nLevels, kst, ken
      real(rkind), intent(in) :: D, L, dx, dy
      real(rkind), dimension(2), intent(in) :: center
      type(decomp_info), intent(in) :: gp
      real(rkind), dimension(:,:,:), intent(inout) :: mask
      integer :: i, j, k
      integer :: ist, ien, jst, jen
      real(rkind) :: xmin, xmax, ymin, ymax
      real(rkind), dimension(2) :: corner

      ! Left side of square
      xmin = center(1) - 0.5d0*L
      xmax = xmin + D
      ymin = center(2) - 0.5d0*L
      ymax = center(2) + 0.5d0*L
#include "oscillatingGrid_files/gridFunc_files/setMaskCommon.F90"
      
      ! Right side of square
      xmax = center(1) + 0.5d0*L
      xmin = xmax - D
      ymin = center(2) - L/2
      ymax = center(2) + L/2
#include "oscillatingGrid_files/gridFunc_files/setMaskCommon.F90" 
      
      ! Bottom of square
      xmin = center(1) - 0.5d0*L + D
      xmax = center(1) + 0.5d0*L - D
      ymin = center(2) - 0.5d0*L
      ymax = ymin + D
#include "oscillatingGrid_files/gridFunc_files/setMaskCommon.F90" 
      
      ! Top of square
      xmin = center(1) - 0.5d0*L + D
      xmax = center(1) + 0.5d0*L - D
      ymax = center(2) + 0.5d0*L
      ymin = ymax - D
#include "oscillatingGrid_files/gridFunc_files/setMaskCommon.F90"

      if (level < nLevels) then
        call getBounds(center,D,L,'NW',corner)
        call addBars(level+1,nLevels,widthFact*D,lengthFact*L,corner,dx,dy,kst,ken,gp,mask)
        call getBounds(center,D,L,'NE',corner)
        call addBars(level+1,nLevels,widthFact*D,lengthFact*L,corner,dx,dy,kst,ken,gp,mask)
        call getBounds(center,D,L,'SE',corner)
        call addBars(level+1,nLevels,widthFact*D,lengthFact*L,corner,dx,dy,kst,ken,gp,mask)
        call getBounds(center,D,L,'SW',corner)
        call addBars(level+1,nLevels,widthFact*D,lengthFact*L,corner,dx,dy,kst,ken,gp,mask)
      end if 
    end subroutine
 
    subroutine getBounds(center,D,L,cornerDesc,corner)
      real(rkind), dimension(2), intent(in) :: center
      real(rkind), intent(in) :: D, L
      character(len=2), intent(in) :: cornerDesc
      real(rkind), dimension(2), intent(out) :: corner

      select case (cornerDesc)
      case ('NW')
        corner = [center(1) - 0.5d0*L + 0.5d0*D, center(2) + 0.5d0*L - 0.5d0*D] 
      case ('NE')
        corner = [center(1) + 0.5d0*L - 0.5d0*D, center(2) + 0.5d0*L - 0.5d0*D] 
      case ('SE')
        corner = [center(1) + 0.5d0*L - 0.5d0*D, center(2) - 0.5d0*L + 0.5d0*D] 
      case ('SW')
        corner = [center(1) - 0.5d0*L + 0.5d0*D, center(2) - 0.5d0*L + 0.5d0*D] 
      end select
    
    end subroutine
    
    subroutine updateBarinfo(barInfo,gID,D,xmin,ymin,xmax,ymax,dx,dy,gp)
      real(rkind), dimension(:,:), intent(inout) :: barInfo
      integer, intent(in) :: gID
      real(rkind), intent(in) :: D, xmin, ymin, xmax, ymax, dx, dy
      class(decomp_info), intent(in) :: gp
      
      !call getStEnIndices(xmin,xmax,dx,gp%xst(1),gp%xen(1),barStEn(gID,1),barStEn(gID,2))
      !call getStEnIndices(ymin,ymax,dy,gp%xst(2),gp%xen(2),barStEn(gID,3),barStEn(gID,4))
      barInfo(gID,1) = xmin
      barInfo(gID,2) = xmax
      barInfo(gID,3) = ymin
      barInfo(gID,4) = ymax
      barInfo(gID,5) = D
    end subroutine

    function getGlobalID(level,myCorner) result(gID)
      integer, intent(in) :: level
      integer, dimension(:), intent(in) :: myCorner
      integer :: gID, gIDfinalLastLevel, n
      
      if (level == 1) then
        gID = 1
      else
        gIDfinalLastLevel = 0
        do n = 1,level-1
          gIDfinalLastLevel = gIDfinalLastLevel + 4**(n-1)
        end do
        gID = gIDfinalLastLevel + 4*(myCorner(level-1) - 1) + myCorner(level)
      end if
    end function


    subroutine set_grid_position(igp, rkStage)
        use oscillating_grid_parameters, only: dumpMaskFreq
        class(igrid), intent(inout) :: igp 
        integer, intent(in) :: rkStage
        real(rkind) :: tnow 
        integer :: ierr
>>>>>>> bca0e56daec689d11e59fe4531fb699163924c45

        select case (rkStage) 
        ! TO: Figure out substage times 
        case (1) 
            tnow = igp%tsim
<<<<<<< HEAD
        
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
=======
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
        call setGridMask(igp%mesh(:,:,:,1), igp%mesh(:,:,:,2), igp%zE           , &
          tnow, igp%immersedBodies(1)%RbodyE, igp%gpE) 

        if (filterMask) then
          igp%rbuffxC(:,:,:,1) = igp%immersedBodies(1)%RbodyC
          call gauss_x%filter1(igp%rbuffxC(:,:,:,1),igp%rbuffxC(:,:,:,2), &
            igp%gpC%xsz(2), igp%gpC%xsz(3))
          call transpose_x_to_y(igp%rbuffxC(:,:,:,2),igp%rbuffyC(:,:,:,1),igp%gpC)
          call gauss_y%filter2(igp%rbuffyC(:,:,:,1),igp%rbuffyC(:,:,:,2), &
            igp%gpC%ysz(1), igp%gpC%ysz(3))
          call transpose_y_to_z(igp%rbuffyC(:,:,:,2),igp%rbuffzC(:,:,:,1),igp%gpC)
          call gauss_zC%filter3(igp%rbuffzC(:,:,:,1),igp%rbuffzC(:,:,:,2), &
            igp%gpC%zsz(1), igp%gpC%zsz(2))
          call transpose_z_to_y(igp%rbuffzC(:,:,:,2),igp%rbuffyC(:,:,:,1),igp%gpC)
          call transpose_y_to_x(igp%rbuffyC(:,:,:,1),igp%immersedBodies(1)%RbodyC,igp%gpC)
          
          igp%rbuffxE(:,:,:,1) = igp%immersedBodies(1)%RbodyE
          call gauss_x%filter1(igp%rbuffxE(:,:,:,1),igp%rbuffxE(:,:,:,2), &
            igp%gpE%xsz(2), igp%gpE%xsz(3))
          call transpose_x_to_y(igp%rbuffxE(:,:,:,2),igp%rbuffyE(:,:,:,1),igp%gpE)
          call gauss_y%filter2(igp%rbuffyE(:,:,:,1),igp%rbuffyE(:,:,:,2), &
            igp%gpE%ysz(1), igp%gpE%ysz(3))
          call transpose_y_to_z(igp%rbuffyE(:,:,:,2),igp%rbuffzE(:,:,:,1),igp%gpE)
          call gauss_zE%filter3(igp%rbuffzE(:,:,:,1),igp%rbuffzE(:,:,:,2), &
            igp%gpE%zsz(1), igp%gpE%zsz(2))
          call transpose_z_to_y(igp%rbuffzE(:,:,:,2),igp%rbuffyE(:,:,:,1),igp%gpE)
          call transpose_y_to_x(igp%rbuffyE(:,:,:,1),igp%immersedBodies(1)%RbodyE,igp%gpE)
        end if
        
        if (rkStage == 1 .and. mod(igp%step,dumpMaskFreq) == 0) then
          call igp%dumpFullField(igp%immersedBodies(1)%RbodyC,'mask',igp%gpC)
        end if

        igp%immersedBodies(1)%uTarget = zero 
        igp%immersedBodies(1)%vTarget = zero 
        igp%immersedBodies(1)%wTarget = 0.5d0*stroke*omega*cos(omega*tnow) ! TODO: Ryan Finish this  
>>>>>>> bca0e56daec689d11e59fe4531fb699163924c45

    end subroutine 

end module 

program oscillatingGrid
    use mpi
    use kind_parameters,  only: clen
    use IncompressibleGrid, only: igrid
<<<<<<< HEAD
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use gridFuncMod
=======
    use igrid_Operators, only: igrid_ops
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use gridFuncMod, only: set_grid_position, gauss_x, gauss_y, gauss_zC, &
      gauss_zE
    use oscillating_grid_parameters, only: filterMask, useGridBreaks, &
      getGridBreakParameters, finalizeProblem
    use fortran_assert, only: assert

>>>>>>> bca0e56daec689d11e59fe4531fb699163924c45
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

<<<<<<< HEAD
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
=======
    if (filterMask) then
      ierr = gauss_x%init( igp%nx,  .true.)
      call assert(ierr == 0,'gauss_x initialization')
      ierr = gauss_y%init( igp%ny,  .true.)
      call assert(ierr == 0,'gauss_y initialization')
      ierr = gauss_zC%init(igp%nz,  .false.)
      call assert(ierr == 0,'gauss_zC initialization')
      ierr = gauss_zE%init(igp%nz+1,.false.)
      call assert(ierr == 0,'gauss_zE initialization')
    endif

    if (useGridBreaks) then
      call getGridBreakParameters()
    end if

    call tic() 
    do while (igp%tsim < igp%tstop .and. igp%step < igp%nsteps) 

       !call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igridWallM.F90)
       !call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
       call set_grid_position(igp, 1)
       call igp%advance_RK4_Stage(stage=1)
       
       call set_grid_position(igp, 2)
       call igp%advance_RK4_Stage(stage=2)

       call set_grid_position(igp, 3)
       call igp%advance_RK4_Stage(stage=3)

       call set_grid_position(igp, 4)
       call igp%advance_RK4_Stage(stage=4)
>>>>>>> bca0e56daec689d11e59fe4531fb699163924c45
       
       call igp%wrapup_timestep()
       
       call doTemporalStuff(igp)                                        

    end do 
 
    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
<<<<<<< HEAD
=======

    call finalizeProblem()
>>>>>>> bca0e56daec689d11e59fe4531fb699163924c45
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
