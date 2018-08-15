program test_PadePoisson_NeumannBCs
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two 
    use spectralMod, only: spectral
    use cd06staggstuff, only: cd06stagg 
    use PadePoissonMod, only: padepoisson
    use PadeDerOps, only: Pade6stagg
    use decomp_2d
    use reductions, only: p_sum, p_maxval
    implicit none

    type(decomp_info) :: gpC, gpE
    integer :: nx = 128, ny = 128, nz = 128
    integer :: ierr, prow = 0, pcol = 0
    
    type(spectral), allocatable :: spect, spectE
    type(padepoisson), allocatable :: poiss 
    type(Pade6stagg) :: Pade6opZ
    integer :: scheme = 1
    integer :: ix1, ixn, iy1, iyn, iz1, izn, i, j, k
    logical :: computeStokesPressure = .false.
    real(rkind), dimension(:,:,:), allocatable :: x, y, z, g, f, ftrue
    real(rkind) :: dx, dy, dz, MaxError
    real(rkind), parameter :: Lx = 2.d0*pi, Ly = 2.d0*pi, Lz = 1.d0 

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, prow, pcol)

    call get_decomp_info(gpC)
    call decomp_info_init(nx, ny, nz+1, gpE)

    dx = Lx/real(nx,rkind)
    dy = Ly/real(ny,rkind)
    dz = Lz/real(nz,rkind)
    ix1 = gpC%xst(1); iy1 = gpC%xst(2); iz1 = gpC%xst(3)
    ixn = gpC%xen(1); iyn = gpC%xen(2); izn = gpC%xen(3)

    allocate(x(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(y(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(z(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(g(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(f(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(ftrue(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))

    do k=1,size(x,3)
        do j=1,size(x,2)
            do i=1,size(x,1)
                x(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                y(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                z(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
            end do
        end do
    end do

    ! Shift everything to the origin 
    x = x - dx
    y = y - dy
    z = z - dz 
    
    g = 34*cos(5*x)*sin(3*y)*(cos(16*pi*z)/(80*pi) - z**2/9.d0 - z**3/27.d0 + z**4/3.d0 &
     - z**5/5.d0 + 7.d0/540.d0) + cos(5*x)*sin(3*y)*((2*z)/9.d0 + (16*pi*cos(16*pi*z))/5.d0 - 4*z**2 + 4*z**3 + 2.d0/9.d0)
    ftrue = -cos(5*x)*sin(3*y)*(cos(16*pi*z)/(80*pi) - z**2/9.d0 - z**3/27.d0 + z**4/3.d0 - z**5/5.d0 + 7.d0/540.d0)

    allocate(spect)
    allocate(spectE)
    allocate(poiss)
    call spect%init("x", nx, ny, nz, dx, dy, dz, "four", "2/3rd", 2,.false.)
    call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", "2/3rd", 2,.false.)
    call Pade6opz%init(gpC, spect%spectdecomp, gpE, spectE%spectdecomp, dz, scheme,.false.)
    call poiss%init(dx, dy, dz, spect, spectE, computeStokesPressure, Lz, .false., gpC, Pade6opz, .false. )
   
    call poiss%PoissonSolver_HomogeneousNeumannBCz(g, f)

    maxError = p_maxval(abs(f - ftrue))/p_maxval(abs(ftrue))
    if (nrank == 0) print*, "Max error:", maxError

    call MPI_Finalize(ierr)

end program 
