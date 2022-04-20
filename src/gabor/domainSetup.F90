module domainSetup
    ! This module stores key domain information and routines to defined the
    ! computational mesh. Here are descriptions of some of the variables:
    !      gpLESb --> grid partition for the LES field including the following boundaries:
    !                 x=Lx, y=Ly, z=0, z=Lz
    !      gpLES  --> grid partition for LES field used by LES solver (i.e.
    !                 excluding the boundaries listed above)
    !      gpQHcent --> grid partition for QH regions (the center of the QH regions)
    !      gpF --> grid partition for the high resolution mesh used to render
    !              the velocity
    !      xLESb --> LES mesh including the boundary points x=Lx, y=Ly, z=0, z=Lz
    !      xLES  --> LES mesh excluding the aforementioned boundary points

    use kind_parameters, only: rkind
    use decomp_2d
    use EulerG_mod, only: EulerG
    use constants, only: pi
    use mpi
    implicit none
    real(rkind), dimension(2) :: xDom, yDom, zDom
    real(rkind) :: Lx, Ly, Lz
    type(decomp_info), allocatable :: gpLESb
    type(decomp_info), allocatable :: gpLES
    type(decomp_info), allocatable :: gpQHcent
    type(decomp_info), allocatable :: gpF
    character(len=1) :: decomp2Dpencil = 'x'
    type(EulerG) :: fineMesh
    logical, dimension(3) :: periodic ! <-- Boundary conditions
    logical :: finishDomainSetup = .false.
 
    ! LES discretization
      real(rkind), dimension(:), allocatable :: xLESb, yLESb, zLESb
      real(rkind), dimension(:), allocatable :: xLES, yLES, zLES
      real(rkind) :: dxLES, dyLES, dzLES
      integer :: nxLES, nyLES, nzLES
    
    ! QH mesh discretization
      real(rkind), dimension(:), allocatable :: xQHedge, yQHedge, zQHedge
      real(rkind), dimension(:), allocatable :: xQHcent, yQHcent, zQHcent
      real(rkind) :: dxQH, dyQH, dzQH
      integer :: nxQH, nyQH, nzQH
      integer :: nxLESperQH, nyLESperQH, nzLESperQH
      real(rkind) :: Wx, Wy, Wz ! Size of the window function support for each Gabor mode
   
    ! High resolution discretization
      integer :: nxF, nyF, nzF
      integer :: nxsupp, nysupp, nzsupp
      real(rkind) :: dxF, dyF, dzF
      ! Number of "levels" to fine mesh. One level constitutes a factor of 2
      ! refinment in each coordinate direction
      integer :: nlevels
      ! Numerical mesh corresponding to high resolution fields
      real(rkind), dimension(:), allocatable :: xF, yF, zF
      ! Numerical mesh for velocity rendering including halo regions for each MPI rank
      real(rkind), dimension(:), allocatable :: xFh, yFh, zFh

    ! Max and min wavenumber based on LES and high resolution meshes
      real(rkind) :: kmin, kmax
    
    contains
      
      subroutine setupDomainXYperiodic(inputfile)
        ! This subroutine defines all relevant numerical meshes which include:
        !   - LES mesh (including boundary points)
        !   - QH edges
        !   - QH centers
        !   - High resolution mesh to render Gabor-induced fields
        character(len=*), intent(in) :: inputfile
        integer :: ist, ien, jst, jen, kst, ken
        integer :: isz, jsz, ksz
        integer :: ierr, ioUnit, pcol, prow
        real(rkind) :: wSupport

        namelist /DOMAIN/ Lx, Ly, Lz, nxLES, nyLES, nzLES, &
          nxLESperQH, nyLESperQH, nzLESperQH, &
          nxF, nyF, nzF, pcol, prow, decomp2Dpencil

        ! Read inputfile
        ioUnit = 1
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=DOMAIN)
        close(ioUnit)

        ! Define domain bounds
        xDom = [0.d0, Lx]
        yDom = [0.d0, Ly]
        zDom = [0.d0, Lz]
        
        ! Allocate memory for grid partitions
        allocate(gpLES,gpLESb,gpQHcent,gpF)

        ! Initialize grid partition for LES mesh
        periodic = [.true.,.true.,.false.]
        call decomp_2d_init(nxLES+1,nyLES+1,nzLES+2,prow,pcol,periodic)
        call get_decomp_info(gpLESb) ! Includes domain boundaries
        call decomp_info_init(nxLES,nyLES,nzLES,gpLES) ! Excludes wall-normal boundaries 
                                                       ! and redundant periodic boundary

        ! Initializ grid partition for QH edges and centers
        nxQH = nxLES/nxLESperQH
        nyQH = nyLES/nyLESperQH
        nzQH = nzLES/nzLESperQH
        call decomp_info_init(nxQH,nyQH,nzQH,gpQHcent)
        
        ! Initialize grid partition for high resolution mesh  
        call decomp_info_init(nxF,nyF,nzF,gpF)

        ! LES field
          call getMeshSpacing(Lx,Ly,Lz,nxLES,nyLES,nzLES,dxLES,dyLES,dzLES)

          ! Define start and end global indidces
            call getStartAndEndIndices(gpLESb,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
          
          ! Allocate memory for the 1D vectors defining coordinate axes
            allocate(xLESb(isz),yLESb(jsz),zLESb(ksz))

          ! Define the physical space grid
            xLESb = getPeriodicNodeValues(ist,ien,dxLES)
            yLESb = getPeriodicNodeValues(jst,jen,dyLES)
            if (kst == 1 .and. ken == nzLES+2) then
              zLESb(1) = 0.d0
              zLESb(nzLES+2) = Lz
              zLESb(2:nzLES+1) = getWallNormalNodeValues(1,nzLES+1,dzLES)
            elseif (kst == 1) then
              zLESb(1) = 0.d0
              zLESb(2:ken) = getWallNormalNodeValues(1,ken,dzLES)
            elseif (ken == nzLES+2) then
              zLESb(ksz) = Lz
              zLESb(1:ksz-1) = getWallNormalNodeValues(kst-1,nzLES+1,dzLES)
            else
              zLESb = getWallNormalNodeValues(kst-1,ken,dzLES)
            end if

          ! Repeat for LES domain without domain boundaries
          call getStartAndEndIndices(gpLES,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
          allocate(xLES(isz),yLES(jsz),zLES(ksz))
          xLES = getPeriodicNodeValues(ist,ien,dxLES)
          yLES = getPeriodicNodeValues(jst,jen,dyLES)
          zLES = getWallNormalNodeValues(kst,ken,dzLES)
  
        ! QH mesh
          call getMeshSpacing(Lx,Ly,Lz,nxQH,nyQH,nzQH,dxQH,dyQH,dzQH)
          call getStartAndEndIndices(gpQHcent,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
          allocate(xQHcent(gpQHcent%xsz(1))  ,yQHcent(gpQHcent%xsz(2))  ,zQHcent(gpQHcent%xsz(3)))
          allocate(xQHedge(gpQHcent%xsz(1)+1),yQHedge(gpQHcent%xsz(2)+1),zQHedge(gpQHcent%xsz(3)+1))
          
          xQHedge = getPeriodicNodeValues(ist,ien+1,dxQH)
          yQHedge = getPeriodicNodeValues(jst,jen+1,dyQH)
          ! While the z-direction is not periodic, the QH edges are equivalently defined by this routine:
          zQHedge = getPeriodicNodeValues(kst,ken+1,dzQH)
          
          xQHcent = 0.5d0*(xQHedge(2:isz+1)+xQHedge(1:isz))
          yQHcent = 0.5d0*(yQHedge(2:jsz+1)+yQHedge(1:jsz))
          zQHcent = 0.5d0*(zQHedge(2:ksz+1)+zQHedge(1:ksz))

          ! Define support window 
          Wx = 2.d0*dxQH
          Wy = 2.d0*dyQH
          Wz = 2.d0*dzQH

        ! High resolution mesh
          nlevels = nint(log2(real(nxF/nxLES,rkind)))
          ! TODO: implement hierarchy of fine meshes
          !call fineMesh%init(gpLESb,nlevels,xDom,yDom,zDom)
          call getMeshSpacing(Lx,Ly,Lz,nxF,nyF,nzF,dxF,dyF,dzF)
          call getStartAndEndIndices(gpF,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
          allocate(xF(isz),yF(jsz),zF(ksz))
          xF = getPeriodicNodeValues(ist,ien,dxF)
          yF = getPeriodicNodeValues(jst,jen,dyF)
          zF = getWallNormalNodeValues(kst,ken,dzF)

          ! High res mesh including halo regions for each MPI process 
          call getWindowFunctionBounds(nlevels,nxsupp,nysupp,nzsupp)
          ist = ist - nxsupp/2; ien = ien + nxsupp/2
          jst = jst - nysupp/2; jen = jen + nysupp/2
          kst = max(1,kst-nzsupp/2)
          ken = min(nzF+1,ken+nzsupp/2) 
        
          !allocate(xFh(isz+nxsupp),yFh(jsz+nysupp),zFh(ksz+nzsupp))  
          allocate(xFh(ist:ien),yFh(jst:jen),zFh(kst:ken))  
          xFh = getPeriodicNodeValues(ist,ien,dxF)
          yFh = getPeriodicNodeValues(jst,jen,dyF)
          zFh = getPeriodicNodeValues(kst,ken,dzF) ! <-- While not periodic in
                                                   ! z, the mesh is defined from 
                                                   ! the wall so uses the same routine
                                                   ! as the periodic directions

        ! kmin and kmax for enrichment
          call computeKminKmax()

        finishDomainSetup = .true.
      end subroutine

      subroutine finalizeDomainSetup()
        if (allocated(gpLESb)) then
          call decomp_info_finalize(gpLESb)
          deallocate(gpLESb)
        end if
        if (allocated(gpLES)) then
          call decomp_info_finalize(gpLES)
          deallocate(gpLES)
        end if
        if (allocated(gpQHcent)) then
          call decomp_info_finalize(gpQHcent)
          deallocate(gpQHcent)
        end if
        if (allocated(gpF)) then
          call decomp_info_finalize(gpF)
          deallocate(gpF)
        end if
        call decomp_2d_finalize

        if (allocated(xLES)) deallocate(xLES)
        if (allocated(yLES)) deallocate(yLES)
        if (allocated(zLES)) deallocate(zLES)
        if (allocated(xLESb)) deallocate(xLESb)
        if (allocated(yLESb)) deallocate(yLESb)
        if (allocated(zLESb)) deallocate(zLESb)
        if (allocated(xQHedge)) deallocate(xQHedge)
        if (allocated(yQHedge)) deallocate(yQHedge)
        if (allocated(zQHedge)) deallocate(zQHedge)
        if (allocated(xQHcent)) deallocate(xQHcent)
        if (allocated(yQHcent)) deallocate(yQHcent)
        if (allocated(zQHcent)) deallocate(zQHcent)
        if (allocated(xF)) deallocate(xF)
        if (allocated(yF)) deallocate(yF)
        if (allocated(zF)) deallocate(zF)
        if (allocated(xFh)) deallocate(xFh)
        if (allocated(yFh)) deallocate(yFh)
        if (allocated(zFh)) deallocate(zFh)
  
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

      subroutine getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
        ! Grab the gp%xst(i), gp%xen(i), and gp%xsz(i) values for more readable code

        use exits, only: gracefulExit
        type(decomp_info), intent(in) :: gp
        integer, intent(out) :: ist, ien, jst, jen, kst, ken
        integer, intent(out) :: isz, jsz, ksz
        integer :: ierr
        
        if (decomp2Dpencil == 'x') then
          ist = gp%xst(1); ien = gp%xen(1)
          jst = gp%xst(2); jen = gp%xen(2)
          kst = gp%xst(3); ken = gp%xen(3)
          
          isz = gp%xsz(1)
          jsz = gp%xsz(2)
          ksz = gp%xsz(3)
        elseif (decomp2Dpencil == 'y') then
          ist = gp%yst(1); ien = gp%yen(1)
          jst = gp%yst(2); jen = gp%yen(2)
          kst = gp%yst(3); ken = gp%yen(3)
          
          isz = gp%ysz(1)
          jsz = gp%ysz(2)
          ksz = gp%ysz(3)
        elseif (decomp2Dpencil == 'z') then
          ist = gp%zst(1); ien = gp%zen(1)
          jst = gp%zst(2); jen = gp%zen(2)
          kst = gp%zst(3); ken = gp%zen(3)
          
          isz = gp%zsz(1)
          jsz = gp%zsz(2)
          ksz = gp%zsz(3)
        else
          call gracefulExit("Must specify decomp2Dpencil as either 'x', 'y', or 'z'",ierr)
        end if
      end subroutine

      pure function log2(x) result(y)
        ! logarithm base 2
        real(rkind), intent(in) :: x
        real(rkind) :: y
        y = log(x)/log(2.d0)
      end function
      
      subroutine getWindowFunctionBounds(level,nxsupp,nysupp,nzsupp)
        integer, intent(in) :: level
        integer, intent(out) :: nxsupp, nysupp, nzsupp
        integer :: nxFine, nyFine, nzFine
        real(rkind) :: dx, dy, dz
        integer :: nxPerQH, nyPerQH, nzPerQH

        nxFine = nxLES*2**(level)
        nyFine = nyLES*2**(level)
        nzFine = nzLES*2**(level)

        dx = Lx/real(nxFine,rkind)
        dy = Ly/real(nyFine,rkind)
        dz = Lz/real(nzFine,rkind)

        nxPerQH = nxFine/nxLES*nxLESperQH
        nyPerQH = nyFine/nyLES*nyLESperQH
        nzPerQH = nzFine/nzLES*nzLESperQH

        nxsupp = 2*nxPerQH
        nysupp = 2*nyPerQH
        nzsupp = 2*nzPerQH
      end subroutine
      
      subroutine computeKminKmax()
        real(rkind) :: kxNyqLES, kyNyqLES, kzNyqLES
        real(rkind) :: kxNyqF, kyNyqF, kzNyqF

        kxNyqLES = getNyquist(Lx,nxLES)
        kyNyqLES = getNyquist(Ly,nyLES)
        kzNyqLES = getNyquist(Lz,nzLES)

        kxNyqF = getNyquist(Lx,nxF)
        kyNyqF = getNyquist(Ly,nyF)
        kzNyqF = getNyquist(Lz,nzF)

        kmin = minval([kxNyqLES, kyNyqLES, kzNyqLES])
        kmax = minval([kxNyqF,   kyNyqF,   kzNyqF])
      end subroutine

      pure function getNyquist(L,n) result(kNyq)
        real(rkind), intent(in) :: L
        integer, intent(in) :: n
        real(rkind) :: kNyq

        kNyq = real(n,rkind)*0.5d0*(2.d0*pi/L)
      end function

end module
