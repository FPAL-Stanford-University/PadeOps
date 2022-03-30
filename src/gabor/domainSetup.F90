module domainSetup
    use kind_parameters, only: rkind
    implicit none
    real(rkind), dimension(2) :: xDom, yDom, zDom
    
    contains
      subroutine get_window_function_bounds(DomX,DomY,DomZ,nxLES,nyLES,nzLES,level,&
          & nxLES_per_QH,nyLES_per_QH,nzLES_per_QH,nxsupp,nysupp,nzsupp,wSupport)
        real(rkind), dimension(2), intent(in) :: DomX, DomY, DomZ
        integer, intent(in) :: nxLES, nyLES, nzLES, level
        integer, intent(in) :: nxLES_per_QH, nyLES_per_QH, nzLES_per_QH
        integer, intent(out) :: nxsupp, nysupp, nzsupp
        real(rkind), intent(out) :: wSupport
        real(rkind) :: Lx, Ly, Lz
        integer :: nxFine, nyFine, nzFine
        real(rkind) :: dx, dy, dz
        integer :: nx_per_QH, ny_per_QH, nz_per_QH

        Lx = DomX(2) - DomX(1)
        Ly = DomY(2) - DomY(1)
        Lz = DomZ(2) - DomZ(1)
  
        nxFine = nxLES*2**(level)
        nyFine = nyLES*2**(level)
        nzFine = nzLES*2**(level)

        dx = Lx/real(nxFine,rkind)
        dy = Ly/real(nyFine,rkind)
        dz = Lz/real(nzFine,rkind)

        nx_per_QH = nxFine/nxLES*nxLES_per_QH
        ny_per_QH = nyFine/nyLES*nyLES_per_QH
        nz_per_QH = nzFine/nzLES*nzLES_per_QH

        nxsupp = 2*nx_per_QH
        nysupp = 2*ny_per_QH
        nzsupp = 2*nz_per_QH

        wSupport = nxsupp*dx
      end subroutine
end module
