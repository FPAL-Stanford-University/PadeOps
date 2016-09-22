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
    use kind_parameters, only: rkind, clen
    use constants,       only: one, two, three
    use test_eos
    use mpi
    implicit none

    real(rkind), dimension(nx,ny,nz) :: rho, e, Entr, p, T, sos, detg
    real(rkind), dimension(nx,ny,nz,6), target :: devstress
    real(rkind), dimension(nx,ny,nz,9), target :: g
    real(rkind), dimension(:,:,:), pointer :: g11,g12,g13,g21,g22,g23,g31,g32,g33
    real(rkind), dimension(:,:,:), pointer :: s11,s12,s13,s22,s23,s33
    real(rkind) :: eosparams(7), rho0
    integer     :: eostype, ierr, IVPflag

    character(len=clen) :: inputfile

    ! Start MPI
    call MPI_Init(ierr)

    ! Get file location 
    call GETARG(1,inputfile)

    open(unit=11, file=trim(inputfile), form='FORMATTED')
    read(11,*) IVPflag

    ! Initialize decomp
    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(decomp)

    eostype = 3
   !eosparams = (/mu0,     beta,      K,       alpha,     C_V,         T_0,         gamma_s  /)
    eosparams = (/39.38d9, 3.0_rkind, 15.28d6, 1.0_rkind, 390.0_rkind, 300.0_rkind, 2.0_rkind/)
    rho0 = 8930.0d0

    call geneos%init(decomp,eostype,eosparams)

    g11 => g(:,:,:,1); g12 => g(:,:,:,2); g13 => g(:,:,:,3)
    g21 => g(:,:,:,4); g22 => g(:,:,:,5); g23 => g(:,:,:,6)
    g31 => g(:,:,:,7); g32 => g(:,:,:,8); g33 => g(:,:,:,9)

    if(IVPflag==1) then
      ! IVP 1 in Lopez-Ortega, JCP 257, 2014.3
      g11 = real(1.25D0 ,rkind); g12 = real(0.0D-3,rkind);  g13 = real(0.0D0 ,rkind)
      g21 = real(0.0D-3,rkind);  g22 = real(1.25D0 ,rkind); g23 = real(0.0D0 ,rkind)
      g31 = real(0.0D0 ,rkind);  g32 = real(0.0D0 ,rkind);  g33 = real(1.25D0 ,rkind)

    elseif(IVPflag==2) then

      ! IVP 2 in Lopez-Ortega, JCP 257, 2014.3
      ! inverse of [0.8 -0.1 0; -0.1 0.9 0; 0 0 0.9]
      g11 = 1.267605633802817_rkind; g12 = 0.140845070422535_rkind; g13 = 0.0_rkind
      g21 = 0.140845070422535_rkind; g22 = 1.126760563380282_rkind; g23 = 0.0_rkind
      g31 = 0.0_rkind;               g32 = 0.0_rkind;               g33 = 1.111111111111111_rkind

    endif

    detg = g11*(g22*g33-g23*g32) + g12*(g23*g31-g21*g33) + g13*(g21*g32-g22*g31)
    T = 300.0_rkind*detg**2*exp(0.5_rkind)

    s11 => devstress(:,:,:,1); s12 => devstress(:,:,:,2); s13 => devstress(:,:,:,3)
                               s22 => devstress(:,:,:,4); s23 => devstress(:,:,:,5)
                                                          s33 => devstress(:,:,:,6)

    rho = rho0*detg
    write(*,*) "Initial g:"
    write(*,'(3(e19.12,1x))') g11(1,1,1), g12(1,1,1), g13(1,1,1)
    write(*,'(3(e19.12,1x))') g21(1,1,1), g22(1,1,1), g23(1,1,1)
    write(*,'(3(e19.12,1x))') g31(1,1,1), g32(1,1,1), g33(1,1,1)
    write(*,*) 
    write(*,*) "Initial T:"
    write(*,'(3(e19.12,1x))') T(1,1,1)


    call geneos%get_e_from_rhoT(rho0, g, rho, T, e)
    call geneos%get_p_devstress_T_sos(rho0, g, rho, e, Entr, p, T, devstress, sos)

    write(*,*) 'rho0, rho, e'
    write(*,'(3(e19.12,1x))') rho0, rho(1,1,1), e(1,1,1)
    write(*,*) 
    write(*,*) 'Entr, p, T, sos'
    write(*,'(3(e19.12,1x))') Entr(1,1,1), p(1,1,1), T(1,1,1), sos(1,1,1)
    write(*,*) 
    write(*,*) 'devstress'
    write(*,'(3(e19.12,1x))') s11(1,1,1), s12(1,1,1), s13(1,1,1)
    write(*,'(3(e19.12,1x))') s12(1,1,1), s22(1,1,1), s23(1,1,1)
    write(*,'(3(e19.12,1x))') s13(1,1,1), s23(1,1,1), s33(1,1,1)
    write(*,*) 'totstress'
    write(*,'(3(e19.12,1x))') -p(1,1,1)+s11(1,1,1), s12(1,1,1),           s13(1,1,1)
    write(*,'(3(e19.12,1x))') s12(1,1,1),           -p(1,1,1)+s22(1,1,1), s23(1,1,1)
    write(*,'(3(e19.12,1x))') s13(1,1,1),           s23(1,1,1),           -p(1,1,1)+s33(1,1,1)

    call MPI_Finalize(ierr)


end program
