module test_godromeos
    use kind_parameters,        only: rkind
    use GodunovRomenskiiEOSMod, only: godromeos
    use decomp_2d,              only: decomp_info, decomp_2d_init, get_decomp_info

    implicit none

    integer, parameter :: nx = 1, ny = 1, nz = 1
    type( decomp_info ):: decomp
    type(godromeos) :: godrom
    real(rkind) :: mu = real(77.D9,rkind)
    real(rkind) :: yield = real(2.49D9,rkind)

contains

end module

program test_godromeos_prog
    use kind_parameters, only: rkind, clen
    use constants,       only: zero, one, two, three
    use AbstractEOSMod,  only: get_invariants
    use test_godromeos
    use mpi
    implicit none

    real(rkind), dimension(nx,ny,nz) :: rho, e, p, T, sos, detg
    real(rkind), dimension(nx,ny,nz,6), target :: devstress
    real(rkind), dimension(nx,ny,nz,9), target :: g
    real(rkind), dimension(:,:,:), pointer :: g11,g12,g13,g21,g22,g23,g31,g32,g33
    real(rkind), dimension(:,:,:), pointer :: s11,s12,s13,s22,s23,s33
    integer     :: ierr, IVPflag

    character(len=clen) :: input

    real(rkind) :: rho0 = 8900, K0 = real(1.528D7,rkind), Cv = 390, T0 = 300, B0 = real(4.41D6,rkind), alpha = one, beta = three, gam = two
    logical :: usegTg = .FALSE.

    ! Start MPI
    call MPI_Init(ierr)

    ! Get file location 
    call GETARG(1,input)

    read(input,*) IVPflag

    ! Initialize decomp
    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(decomp)

    ! Initialize godrom using the constructor
    godrom = godromeos(rho0,K0,Cv,T0,B0,alpha,beta,gam,usegTg)

    g11 => g(:,:,:,1); g12 => g(:,:,:,2); g13 => g(:,:,:,3)
    g21 => g(:,:,:,4); g22 => g(:,:,:,5); g23 => g(:,:,:,6)
    g31 => g(:,:,:,7); g32 => g(:,:,:,8); g33 => g(:,:,:,9)

    if(IVPflag==1) then
        ! IVP 1 (Left) in Barton 2009, Section 2.6
        g11 = one/real(0.98D0 ,rkind);    g12 = zero;  g13 = zero
        g21 = real(-0.02D0/0.98D0,rkind); g22 = one;   g23 = real(-0.1,rkind)
        g31 = zero;                       g32 = zero;  g33 = one

        detg = g11*(g22*g33-g23*g32) + g12*(g23*g31-g21*g33) + g13*(g21*g32-g22*g31)
        T = T0 * (detg**gam) * exp(1000._rkind/Cv)
    endif


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


    ! call godrom%get_finger(g,finger,fingersq)
    ! call get_invariants(finger,I1,I2,I3)

    call godrom%get_e_from_rho_g_T(rho, g, T, e)
    call godrom%get_p_devstress_T_sos2(g, rho, e, p, T, devstress, sos)

    write(*,*) 'rho0, rho, e'
    write(*,'(3(e19.12,1x))') rho0, rho(1,1,1), e(1,1,1)
    write(*,*) 
    write(*,*) 'p, T, sos'
    write(*,'(3(e19.12,1x))') p(1,1,1), T(1,1,1), sos(1,1,1)
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
