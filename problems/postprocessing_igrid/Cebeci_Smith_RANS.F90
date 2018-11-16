module NS_1d_rhs
    use kind_parameters, only: rkind
    use cd06staggstuff, only:cd06stagg
    implicit none

    real(rkind), parameter :: dz_match = 5.d0/1000.d0 !50.d0/1000.d0
    real(rkind) :: umatch = 0.010d0 !14.8439612264d0
    integer, parameter :: nz = 400
    real(rkind), dimension(:), allocatable :: nu_t, z, zE

    real(rkind) :: ustar = 1.d0 
    integer :: matching_mode = 1
    integer :: scheme = 1 ! 1: 2nd order, 2: CD06
    real(rkind), parameter :: Retau = 1000
    
    type(cd06stagg) :: der
    real(rkind) :: dz


contains

    subroutine get_RHS(nz, u, rhs, tmpC, dudz, uE, F, nu_t)
        integer, intent(in) :: nz
        real(rkind), intent(in) :: F
        real(rkind), intent(in) , dimension(1,1,nz) :: u
        real(rkind), intent(out), dimension(1,1,nz) :: rhs, dudz, tmpC
        real(rkind), intent(out), dimension(1,1,nz+1) :: uE
        real(rkind), intent(in), dimension(nz) :: nu_t
        integer :: idx

        rhs = 0.d0 
        select case(scheme)
        case(1)
            uE(1,1,2:nz) = 0.5d0*(u(1,1,2:nz) + u(1,1,1:nz-1))
        case(2)
            call der%InterpZ_C2E(u,uE,1,1)
        end select
        
        if (matching_mode == 1) then
            umatch = get_utarget(uE(1,1,2),(1+zE(2)),(1+zE(1)))
        end if 
        ! Dirichlet BCs
        uE(1,1,1) = umatch 
        uE(1,1,nz+1) = umatch

        select case(scheme)
        case(1)
            dudz = (uE(:,:,2:nz+1) - uE(:,:,1:nz))/dz
        case(2)
            call der%ddz_E2C(uE,dudz,1,1)
        end select 
        
        do idx = 1,nz
            tmpC(1,1,idx) = ((1./Retau) + nu_t(idx))*dudz(1,1,idx)
        end do 

        select case(scheme)
        case(1)
            rhs(1,1,2:nz-1) = (tmpC(1,1,3:nz) - tmpC(1,1,1:nz-2))/(2.d0*dz)
            rhs(1,1,1) = ((-3.d0/2.d0)*tmpC(1,1,1) + 2.d0*tmpC(1,1,2) - (1.d0/2.d0)*tmpC(1,1,3))/dz
            rhs(1,1,nz) = rhs(1,1,1) ! Symmetry about centerline 
        case(2)
            call der%ddz_C2C(tmpC,rhs,1,1)
        end select 
        rhs = rhs + F

    end subroutine 

    pure elemental function get_musker_profile(yp) result(up)
        real(rkind), intent(in) :: yp
        real(rkind) :: up

        up = 5.424d0*atan(0.11976047904*yp - 0.48802395209580833) + 0.434d0*log(((yp + 10.6)**9.6)/((yp**2 - 8.15*yp + 86)**2)) - 3.50727901936264842

    end function

    function find_utau(utau_guess, umatch, zmatch) result(utau)
        real(rkind), intent(in) :: utau_guess, umatch, zmatch
        real(rkind) :: zplus, utau_diff, up1, utaunew
        real(rkind) :: utau

        utau_diff = utau_guess
        utau = utau_guess
        do while (utau_diff > 1.d-8)
            zplus = zmatch*Retau*utau 
            up1 = get_musker_profile(zplus)                
            utaunew = umatch/up1
            utau_diff = abs(utaunew - utau)
            utau = utaunew
        end do 

        ! print*, tau
    end function 

        
    function get_utarget(u1,z1,z0) result(u0)
        real(rkind), intent(in) :: u1, z1, z0
        real(rkind) ::  up0
        real(rkind) :: u0, ustarnew

        ustarnew = find_utau(ustar, u1, z1)
        ustar = ustarnew

        up0 = get_musker_profile(z0*Retau*ustar)
        u0 = up0*ustar

    end function 
end module 


program Cebeci_smith_RANS
    use kind_parameters, only: rkind, clen 
!    use cd06staggstuff, only:cd06stagg
 !   use staggOps, 
    use basic_io
    use NS_1d_rhs
    use interpolation

    implicit none
    real(rkind), dimension(:,:), allocatable :: data2read

    real(rkind), parameter :: tstop = 5000.0d0 
    real(rkind), parameter :: cvisc =  0.4d0
    real(rkind), dimension(:,:,:), allocatable :: dt 
    real(rkind), dimension(:,:,:), allocatable :: u, uE, tmpC, rhs, dudz, uLam
    real(rkind), dimension(:,:,:), allocatable :: uin, k1, k2, k3, k4
    real(rkind), parameter :: F = 1.d0 
    character(len=clen) :: fname, input_fname = "/scratch/04076/tg833754/Cebecci_smith/Cebeci_smith_ReTau1000.txt"
    character(len=clen) :: restart_fname = "/scratch/04076/tg833754/Cebecci_smith/u_CS_WM0002400000.dat"
    logical :: restartSim = .false. 
    integer :: idx, tprint = 100000, tdatadump = 5000000 
    integer, parameter :: tscheme = 3
    real(rkind), dimension(:), allocatable :: bsp, csp, dsp 
    real(rkind), dimension(:,:), allocatable :: data2read2

    ! get the nu_t from disk
    call read_2d_ascii(data2read,input_fname)

    allocate(bsp(size(data2read,1)))
    allocate(csp(size(data2read,1)))
    allocate(dsp(size(data2read,1)))

    call spline (data2read(:,1), data2read(:,2), bsp, csp, dsp, size(bsp))


    !nz = size(data2read,1)
    print*, "Nz:", nz

    allocate(z(nz),nu_t(nz))
    allocate(u(1,1,nz))
    allocate(uE(1,1,nz+1))
    allocate(dudz(1,1,nz))
    allocate(tmpC(1,1,nz))
    allocate(uLam(1,1,nz))
    allocate(rhs(1,1,nz))

    allocate(zE(nz+1))

    allocate(dt(1,1,nz))
    dz = (2.d0 - 2.d0*dz_match)/nz
    z(1) = (-1.d0 + dz_match) + dz/2.d0;
    do idx = 2,nz
        z(idx) = z(idx-1)+dz
    end do 
   
    zE(1:nz) = z - dz/2.d0 
    zE(nz+1) = zE(nz) + dz

    ! Interpolate nu_t
    do idx = 1,nz
       nu_t(idx) = ispline(z(idx), data2read(:,1), data2read(:,2), bsp, csp, dsp, size(bsp)) 
    end do
    !nu_t = 0.d0 

    uLam(1,1,:) = 0.5*Retau*F*(1 - z*z)

    !deallocate(data2read)
    ! print*, z(45), nu_t(45)

    call der%init(nz, dz, isTopEven=.false., isBotEven=.false., isTopSided=.true., isBotSided=.true.)
    ! Initialize u
    if (restartSim) then
        deallocate(data2read)
        allocate(data2read(nz,2))
        call read_ascii_2d_2rows(data2read,restart_fname,nz)
        u(1,1,:) = data2read(:,2)
    else
        u = 1.d0 !am 
    end if 

    !nu_t = 0.d0 ! temporary


    do idx = 1,nz
        dt(1,1,idx) = cvisc*dz*dz/(1.d0/Retau + nu_t(idx))
    end do 
    dt = 1.d-4

    allocate(data2read2(nz+1,2))

    allocate(k1(1,1,nz))
    allocate(k2(1,1,nz))
    allocate(k3(1,1,nz))
    allocate(k4(1,1,nz))
    allocate(uin(1,1,nz))

    deallocate(data2read)
    allocate(data2read(nz,2))
    idx = 1
    do while(idx<100000000)
        ! Interpolate cell to edge

        select case (tscheme)
        case (0)
            call get_RHS(nz, u, k1, tmpC, dudz, uE, F, nu_t)
            u = u + k1*dt
        case (1)
            uin = u
            call get_RHS(nz, uin, k1, tmpC, dudz, uE, F, nu_t)
            k1 = dt*k1

            uin = u + 0.5d0*k1
            call get_RHS(nz, uin, k2, tmpC, dudz, uE, F, nu_t)
            k2 = dt*k2

            uin = u + 0.5d0*k2
            call get_RHS(nz, uin, k3, tmpC, dudz, uE, F, nu_t)
            k3 = dt*k3

            uin = u + k3
            call get_RHS(nz, uin, k4, tmpC, dudz, uE, F, nu_t)
            k4 = dt*k4

            u = u + (1.d0/6.d0)*(k1 + k2 + k3 + k4)
        case (2)
            call get_RHS(nz, u, k1, tmpC, dudz, uE, F, nu_t)
            k1 = u + dt*k1

            call get_RHS(nz, k1, k2, tmpC, dudz, uE, F, nu_t)
            k2 = (3.d0/4.d0)*u + (1.d0/4.d0)*k1 + (1.d0/4.d0)*dt*k2

            call get_RHS(nz, k2, k3, tmpC, dudz, uE, F, nu_t)
            
            u = (1.d0/3.d0)*u + (2.d0/3.d0)*k2 + (2.d0/3.d0)*dt*k3
        case (3)
            call get_RHS(nz, u, k1, tmpC, dudz, uE, F, nu_t)
            k1 = u + dt*k1

            call get_RHS(nz, k1, k2, tmpC, dudz, uE, F, nu_t)
            
            u = 0.5d0*u + 0.5d0*k1 + 0.5d0*dt*k2
        end select

        !t = t + dt
        idx = idx + 1
        if (mod(idx,tprint) == 0) then
            print*, "ubottom:", u(1,1,1)
            print*, "uquarte:", u(1,1,nz/4)
            print*, "ucenter:", u(1,1,nz/2)
            print*, "-----------" 
            print*, "utau:", ustar
            print*, "maxvalu:", maxval(u)
            print*, "minvalu:", minval(u)
            print*, "-----------" 
            print*, "-----------" 
            !print*,"Laminar error:", maxval(abs(u - ulam))
        end if

        if (mod(idx,tdatadump) == 0) then
            write(fname,"(A7,I10.10,A4)") "u_CS_2O",idx, ".dat"
            data2read(:,1) = z
            data2read(:,2) = u(1,1,:)
            call write_2d_ascii(data2read,trim(fname))
            
            write(fname,"(A7,I10.10,A4)") "u_CS_ED",idx, ".dat"
            data2read2(:,1) = (zE+1)*ustar*Retau
            data2read2(:,2) = uE(1,1,:)/ustar
            call write_2d_ascii(data2read2,trim(fname))
        end if 
    end do 

end program 
