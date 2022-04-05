module domainSetup
    use kind_parameters, only: rkind
    use decomp_2d
    use EulerG_mod, only: EulerG
    implicit none
    real(rkind), dimension(2) :: xDom, yDom, zDom
    real(rkind) :: Lx, Ly, Lz
    type(decomp_info), allocatable :: gpLES, gpQHcent
    type(EulerG) :: fineMesh
 
    ! LES discretization
      real(rkind), dimension(:), allocatable :: xLES, yLES, zLES
      real(rkind) :: dxLES, dyLES, dzLES
      integer :: nxLES, nyLES, nzLES
    
    ! QH mesh discretization
      real(rkind), dimension(:), allocatable :: xQHedge, yQHedge, zQHedge
      real(rkind), dimension(:), allocatable :: xQHcent, yQHcent, zQHcent
      real(rkind) :: dxQH, dyQH, dzQH
      integer :: nxQH, nyQH, nzQH
   
    ! High resolution discretization
      integer :: nxF, nyF, nzF
      ! Number of "levels" to fine mesh. One level constitutes a factor of 2
      ! refinment in each coordinate direction
      integer :: nlevels
    
    contains
      
      subroutine setupDomainXYperiodic(inputfile)
        ! This subroutine defines all relevant numerical meshes which include:
        !   - LES mesh (including boundary points)
        !   - QH edges
        !   - QH centers
        !   - High resolution mesh to render Gabor-induced fields
        character(len=*), intent(in) :: inputfile
        real(rkind) :: Lx, Ly, Lz
        integer :: ist, ien, jst, jen, kst, ken
        integer :: ierr, pcol, prow

        ! QH mesh
        integer :: nxLESperQH, nyLESperQH, nzLESperQH

        namelist /DOMAIN/ Lx, Ly, Lz, nxLES, nyLES, nzLES, &
          nxLESperQH, nyLESperQH, nzLESperQH, &
          nxF, nyF, nzF, pcol, prow

        ! Read inputfile
        open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=123, NML=DOMAIN)
        close(123)

        ! Define domain bounds
        xDom = [0.d0, Lx]
        yDom = [0.d0, Ly]
        zDom = [0.d0, Lz]
        
        ! Allocate memory for grid partitions
        allocate(gpLES,gpQHcent)

        ! Initialize grid partition for LES mesh
        call decomp_2d_init(nxLES+1,nyLES+1,nzLES+2,prow,pcol,[.true.,.true.,.false.])
        call get_decomp_info(gpLES)

        ! Initializ grid partition for QH edges and centers
        nxQH = nxLES/nxLESperQH
        nyQH = nyLES/nyLESperQH
        nzQH = nzLES/nzLESperQH
        call decomp_info_init(nxQH,nyQH,nzQH,gpQHcent)

        ! LES field
          call getMeshSpacing(Lx,Ly,Lz,nxLES,nyLES,nzLES,dxLES,dyLES,dzLES)

          ! Define start and end global indidces
            call getStartAndEndIndices(gpLES,ist,ien,jst,jen,kst,ken)
          
          ! Allocate memory for the 1D vectors defining coordinate axes
            allocate(xLES(ist:ien),yLES(jst:jen),zLES(kst:ken))

          ! Define the physical space grid
            xLES = getPeriodicNodeValues(ist,ien,dxLES)
            yLES = getPeriodicNodeValues(jst,jen,dyLES)
            if (kst == 1 .and. ken == nzLES+2) then
              zLES(1) = 0.d0
              zLES(nzLES+2) = Lz
              zLES(2:nzLES+1) = getWallNormalNodeValues(1,nzLES+1,dzLES)
            elseif (kst == 1) then
              zLES(1) = 0.d0
              zLES(2:ken) = getWallNormalNodeValues(1,ken,dzLES)
            elseif (ken == nzLES+2) then
              zLES(nzLES+2) = Lz
              zLES(kst:nzLES+1) = getWallNormalNodeValues(kst-1,nzLES+1,dzLES)
            else
              zLES = getWallNormalNodeValues(kst-1,ken,dzLES)
            end if
  
        ! QH mesh
          call getMeshSpacing(Lx,Ly,Lz,nxQH,nyQH,nzQH,dxQH,dyQH,dzQH)
          call getStartAndEndIndices(gpQHcent,ist,ien,jst,jen,kst,ken)
          allocate(xQHcent(ist:ien),yQHcent(jst:jen),zQHcent(kst:ken))
          allocate(xQHedge(ist:ien+1),yQHedge(jst:jen+1),zQHedge(kst:ken+1))
          
          xQHedge = getPeriodicNodeValues(ist,ien+1,dxQH)
          yQHedge = getPeriodicNodeValues(jst,jen+1,dyQH)
          ! While the z-direction is not periodic, the QH edges are equivalently defined by this routine:
          zQHedge = getPeriodicNodeValues(kst,ken+1,dzQH)
          
          xQHcent = 0.5d0*(xQHedge(ist+1:ien+1)+xQHedge(ist:ien))
          yQHcent = 0.5d0*(yQHedge(jst+1:jen+1)+yQHedge(jst:jen))
          zQHcent = 0.5d0*(zQHedge(kst+1:ken+1)+zQHedge(kst:ken))

        ! High resolution mesh
          nlevels = nint(log2(real(nxF/nxLES,rkind)))
          call fineMesh%init(gpLES,nlevels,xDom,yDom,zDom)
      end subroutine

      subroutine finalizeDomainSetup()
        if (allocated(gpLES)) then
          call decomp_info_finalize(gpLES)
          deallocate(gpLES)
        end if
        if (allocated(gpQHcent)) then
          call decomp_info_finalize(gpQHcent)
          deallocate(gpQHcent)
        end if

        if (allocated(xLES)) deallocate(xLES)
        if (allocated(yLES)) deallocate(yLES)
        if (allocated(zLES)) deallocate(zLES)
        if (allocated(xQHedge)) deallocate(xQHedge)
        if (allocated(yQHedge)) deallocate(yQHedge)
        if (allocated(zQHedge)) deallocate(zQHedge)
        if (allocated(xQHcent)) deallocate(xQHcent)
        if (allocated(yQHcent)) deallocate(yQHcent)
        if (allocated(zQHcent)) deallocate(zQHcent)
  
      end subroutine


      pure function getPeriodicNodeValues(st,en,h) result(x)
        ! Define a 1D vector based on global indices and uniform grid spacing
        ! The global mesh begins at the domain boundary
        ! Inputs:
        !   st, en --> start and ending indices of vector
        !   h --> grid spacing
        integer, intent(in) :: st, en
        real(rkind), intent(in) :: h
        real(rkind), dimension(st:en) :: x
        integer :: i
           
        do i = st,en
          x(i) = (i - 1)*h
        end do
      end function
      pure function getWallNormalNodeValues(st,en,h) result(x)
        ! Define a 1D vector based on global indices and uniform grid spacing
        ! The global mesh is defined as h/2:h:L-h/2
        ! Inputs:
        !   st, en --> start and ending indices of vector
        !   h --> grid spacing
        integer, intent(in) :: st, en
        real(rkind), intent(in) :: h
        real(rkind), dimension(st:en) :: x
        integer :: i
           
        do i = st,en
          x(i) = 0.5d0*h + (i - 1)*h
        end do

      end function

      pure subroutine getMeshSpacing(Lx,Ly,Lz,nx,ny,nz,dx,dy,dz)
        ! Routing to compute the uniform grid spacing
        ! Inputs:
        !     L --> Domain length in x, y, and z
        !     N --> Number of grid points in x, y, and z
        !     h --> Grid spacing in x, y, and z
        real(rkind), intent(in) :: Lx, Ly, Lz
        integer, intent(in) :: nx, ny, nz
        real(rkind), intent(out) :: dx, dy, dz
        
        dx = Lx/real(nx,rkind)
        dy = Ly/real(ny,rkind)
        dz = Lz/real(nz,rkind)
      end subroutine

      pure subroutine getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken)
        ! Grab the gp%xst(i) and gp%xen(i) values for more readable code
        type(decomp_info), intent(in) :: gp
        integer, intent(out) :: ist, ien, jst, jen, kst, ken
        
        ist = gp%xst(1); ien = gp%xen(1)
        jst = gp%xst(2); jen = gp%xen(2)
        kst = gp%xst(3); ken = gp%xen(3)
      end subroutine

      pure function log2(x) result(y)
        ! logarithm base 2
        real(rkind), intent(in) :: x
        real(rkind) :: y
        y = log(x)/log(2.d0)
      end function
      
      subroutine getWindowFunctionBounds(nxLES,nyLES,nzLES,level,&
          & nxLESperQH,nyLESperQH,nzLESperQH,nxsupp,nysupp,nzsupp,wSupport)
        integer, intent(in) :: nxLES, nyLES, nzLES, level
        integer, intent(in) :: nxLESperQH, nyLESperQH, nzLESperQH
        integer, intent(out) :: nxsupp, nysupp, nzsupp
        real(rkind), intent(out) :: wSupport
        integer :: nxFine, nyFine, nzFine
        real(rkind) :: dx, dy, dz
        integer :: nx_per_QH, ny_per_QH, nz_per_QH

        nxFine = nxLES*2**(level)
        nyFine = nyLES*2**(level)
        nzFine = nzLES*2**(level)

        dx = Lx/real(nxFine,rkind)
        dy = Ly/real(nyFine,rkind)
        dz = Lz/real(nzFine,rkind)

        nx_per_QH = nxFine/nxLES*nxLESperQH
        ny_per_QH = nyFine/nyLES*nyLESperQH
        nz_per_QH = nzFine/nzLES*nzLESperQH

        nxsupp = 2*nx_per_QH
        nysupp = 2*ny_per_QH
        nzsupp = 2*nz_per_QH

        wSupport = nxsupp*dx
      end subroutine
end module
