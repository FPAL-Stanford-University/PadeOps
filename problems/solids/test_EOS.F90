module test_eos
    use kind_parameters, only: rkind
    use GeneralMatEOS,    only: generaleos
    use decomp_2d,       only: decomp_info, decomp_2d_init, get_decomp_info

    implicit none

    integer, parameter :: nx = 1, ny = 1, nz = 1
    type( decomp_info ):: decomp
    type(generaleos) :: geneos
    real(rkind) :: mu = real(77.D9,rkind)
    real(rkind) :: yield = real(2.49D9,rkind)

contains

end module

program test_eos_prog
    use kind_parameters, only: rkind
    use constants,       only: one, two, three
    use test_eos
    use mpi
    implicit none

    real(rkind), dimension(nx,ny,nz) :: rho, e, Entr, p, T, sos
    real(rkind), dimension(nx,ny,nz,6), target :: devstress
    real(rkind), dimension(nx,ny,nz,9), target :: g
    real(rkind), dimension(:,:,:), pointer :: g11,g12,g13,g21,g22,g23,g31,g32,g33
    real(rkind), dimension(:,:,:), pointer :: s11,s12,s13,s22,s23,s33
    real(rkind) :: eosparams(7), rho0
    integer     :: eostype, ierr

    call MPI_INIT(ierr)

    ! Initialize decomp
    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(decomp)

    eostype = 3
    eosparams = (/39.38d0, 3.0_rkind, 15.28d6, 1.0_rkind, 390.0_rkind, 300.0_rkind, 2.0_rkind/)

    call geneos%init(decomp,eostype,eosparams)

    g11 => g(:,:,:,1); g12 => g(:,:,:,2); g13 => g(:,:,:,3)
    g21 => g(:,:,:,4); g22 => g(:,:,:,5); g23 => g(:,:,:,6)
    g31 => g(:,:,:,7); g32 => g(:,:,:,8); g33 => g(:,:,:,9)

    g11 = real(1.0D0 ,rkind); g12 = real(1.0D-3,rkind); g13 = real(0.0D0 ,rkind)
    g21 = real(2.0D-3,rkind); g22 = real(0.0D0 ,rkind); g23 = real(0.0D0 ,rkind)
    g31 = real(0.0D0 ,rkind); g32 = real(0.0D0 ,rkind); g33 = real(0.0D0 ,rkind)

    g11 = g11 + one; g22 = g22 + one; g33 = g33 + one

    rho0 = 8930.0d0
    rho = 9000.0d0
    e = 100.0d9

    call geneos%get_p_devstress_T_sos(rho0, g, rho, e, Entr, p, T, devstress, sos)

    write(*,*) "Initial g:"
    write(*,*) g11(1,1,1), g12(1,1,1), g13(1,1,1)
    write(*,*) g21(1,1,1), g22(1,1,1), g23(1,1,1)
    write(*,*) g31(1,1,1), g32(1,1,1), g33(1,1,1)
    write(*,*) 

    write(*,*) 'rho0, rho, e'
    write(*,*) rho0, rho(1,1,1), e(1,1,1)
    write(*,*) 
    write(*,*) 'Entr, p, T, sos'
    write(*,*) Entr(1,1,1), p(1,1,1), T(1,1,1), sos(1,1,1)
    write(*,*) 
    write(*,*) 'devstress'
    write(*,*) s11(1,1,1), s12(1,1,1), s13(1,1,1)
    write(*,*) s22(1,1,1), s23(1,1,1)
    write(*,*) s33(1,1,1)

    call MPI_Finalize(ierr)


end program
