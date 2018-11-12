module NS_1d_rhs
    use kind_parameters, only: rkind
    use cd06staggstuff, only:cd06stagg

    implicit none

    real(rkind), parameter :: dz_match = 50.47d0/1000.d0
    real(rkind), parameter :: umatch = 14.87d0
    integer, parameter :: nz = 32
contains

    subroutine get_RHS(nz, u, rhs, tmpC, dudz, uE, F, der, Retau, nu_t)
        integer, intent(in) :: nz
        real(rkind), intent(in) :: F, Retau
        real(rkind), intent(in) , dimension(1,1,nz) :: u
        real(rkind), intent(out), dimension(1,1,nz) :: rhs, dudz, tmpC
        real(rkind), intent(out), dimension(1,1,nz+1) :: uE
        real(rkind), intent(in), dimension(nz) :: nu_t
        type(cd06stagg), intent(inout) :: der
        integer :: idx

        rhs = 0.d0 
        call der%InterpZ_C2E(u,uE,1,1)
        ! Dirichlet BCs
        uE(1,1,1) = umatch 
        uE(1,1,nz+1) = umatch
        
        call der%ddz_E2C(uE,dudz,1,1)
        do idx = 1,nz
            tmpC(1,1,idx) = ((1./Retau) + nu_t(idx))*dudz(1,1,idx)
        end do 
        call der%ddz_C2C(tmpC,rhs,1,1)
        rhs = rhs + F

    end subroutine 
end module 


program Cebeci_smith_RANS
    use kind_parameters, only: rkind, clen 
    use cd06staggstuff, only:cd06stagg
    use basic_io, only: read_2d_ascii, write_2d_ascii
    use NS_1d_rhs
    use interpolation

    implicit none
    real(rkind), dimension(:,:), allocatable :: data2read
    real(rkind), dimension(:), allocatable :: nu_t, z

    real(rkind), parameter :: Retau = 1000
    real(rkind), parameter :: tstop = 5000.0d0 
    real(rkind), parameter :: cvisc =  0.4d0
    real(rkind), dimension(:,:,:), allocatable :: dt 
    real(rkind), dimension(:,:,:), allocatable :: u, uE, tmpC, rhs, dudz, uLam
    real(rkind), dimension(:,:,:), allocatable :: uin, k1, k2, k3, k4
    real(rkind), parameter :: F = 1.d0 
    character(len=clen) :: fname, input_fname = "/scratch/04076/tg833754/Cebecci_smith/Cebeci_smith_ReTau1000.txt"
    character(len=clen) :: restart_fname = "/home1/04076/tg833754/Codes/PadeOps/build/build_opt/problems/postprocessing_igrid/u_CS_NS0005460000.dat"
    logical :: restartSim = .false. 
    real(rkind) :: dz
    type(cd06stagg) :: der
    integer :: idx, tprint = 1000, tdatadump = 100000 
    integer, parameter :: tscheme = 3
    real(rkind), dimension(:), allocatable :: bsp, csp, dsp 


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

    allocate(dt(1,1,nz))
    dz = (2.d0 - 2.d0*dz_match)/nz
    z(1) = (-1.d0 + dz_match) + dz/2.d0;
    do idx = 2,nz
        z(idx) = z(idx-1)+dz
    end do 
   
    ! Interpolate nu_t
    do idx = 1,nz
       nu_t(idx) = ispline(z(idx), data2read(:,1), data2read(:,2), bsp, csp, dsp, size(bsp)) 
    end do
    !nu_t = data2read(:,2)
   
    uLam(1,1,:) = 0.5*Retau*F*(1 - z*z)

    !deallocate(data2read)
    ! print*, z(45), nu_t(45)

    call der%init(nz, dz, isTopEven=.false., isBotEven=.false., isTopSided=.true., isBotSided=.true.)
    ! Initialize u
    if (restartSim) then
        deallocate(data2read)
        call read_2d_ascii(data2read,restart_fname)
        print*, allocated(data2read)
        !print*, size(data2read)
        u(1,1,:) = data2read(2,:)
    else
        u = 0.d0 !am 
    end if 

    !nu_t = 0.d0 ! temporary

    do idx = 1,nz
        dt(1,1,idx) = cvisc*dz*dz/(1.d0/Retau + nu_t(idx))
    end do 
    dt = 1.d-4

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
            call get_RHS(nz, u, k1, tmpC, dudz, uE, F, der, Retau, nu_t)
            u = u + k1*dt
        case (1)
            uin = u
            call get_RHS(nz, uin, k1, tmpC, dudz, uE, F, der, Retau, nu_t)
            k1 = dt*k1

            uin = u + 0.5d0*k1
            call get_RHS(nz, uin, k2, tmpC, dudz, uE, F, der, Retau, nu_t)
            k2 = dt*k2

            uin = u + 0.5d0*k2
            call get_RHS(nz, uin, k3, tmpC, dudz, uE, F, der, Retau, nu_t)
            k3 = dt*k3

            uin = u + k3
            call get_RHS(nz, uin, k4, tmpC, dudz, uE, F, der, Retau, nu_t)
            k4 = dt*k4

            u = u + (1.d0/6.d0)*(k1 + k2 + k3 + k4)
        case (2)
            call get_RHS(nz, u, k1, tmpC, dudz, uE, F, der, Retau, nu_t)
            k1 = u + dt*k1

            call get_RHS(nz, k1, k2, tmpC, dudz, uE, F, der, Retau, nu_t)
            k2 = (3.d0/4.d0)*u + (1.d0/4.d0)*k1 + (1.d0/4.d0)*dt*k2

            call get_RHS(nz, k2, k3, tmpC, dudz, uE, F, der, Retau, nu_t)
            
            u = (1.d0/3.d0)*u + (2.d0/3.d0)*k2 + (2.d0/3.d0)*dt*k3
        case (3)
            call get_RHS(nz, u, k1, tmpC, dudz, uE, F, der, Retau, nu_t)
            k1 = u + dt*k1

            call get_RHS(nz, k1, k2, tmpC, dudz, uE, F, der, Retau, nu_t)
            
            u = 0.5d0*u + 0.5d0*k1 + 0.5d0*dt*k2
        end select

        !t = t + dt
        idx = idx + 1
        if (mod(idx,tprint) == 0) then
            print*, "ubottom:", u(1,1,1)
            print*, "uquarte:", u(1,1,nz/4)
            print*, "ucenter:", u(1,1,nz/2)
            print*, "-----------" 
            print*, "maxvalu:", maxval(u)
            print*, "minvalu:", minval(u)
            print*, "-----------" 
            print*, "-----------" 
        end if

        if (mod(idx,tdatadump) == 0) then
            write(fname,"(A7,I10.10,A4)") "u_CS_WM",idx, ".dat"
            data2read(:,1) = z
            data2read(:,2) = u(1,1,:)
            call write_2d_ascii(data2read,trim(fname))
        end if 
    end do 

    data2read(:,1) = z
    data2read(:,2) = u(1,1,:)
    call write_2d_ascii(data2read,"u_NS_CebiciSmith.dat")

end program 
