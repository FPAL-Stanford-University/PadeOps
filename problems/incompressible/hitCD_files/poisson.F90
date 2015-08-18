module poissonMod
    use kind_parameters, only: rkind
    use decomp_2d, only: decomp_info
    use fft_3d_stuff, only: fft_3d
    use constants, only: one, zero 
    use exits, only: GracefulExit 
    implicit none 

    type :: poisson
        private
        type(decomp_info), allocatable :: myGP
        real(rkind), dimension(:,:,:,:), allocatable :: k_tensor 
        type(fft_3d), allocatable :: FT
        real(rkind), dimension(:,:,:), allocatable :: Gfilt
        complex(rkind), dimension(:,:,:,:), allocatable :: ustar_hat
        complex(rkind), dimension(:,:,:,:), allocatable :: u_hat
        real(rkind), dimension(:,:,:,:), allocatable :: u_in_x

        integer :: nx_p_out, ny_p_out, nz_p_out
        integer :: nx_p_in, ny_p_in, nz_p_in

        contains
            procedure :: init
            procedure :: destroy

            procedure :: pressureProjection
    end type

contains
    
    pure elemental function GetFilterTransferFunction(kin,dx) result(T)
        use constants, only: two
        use cf90stuff, only: alpha90, beta90, a90, b90, c90, d90, e90
        real(rkind), intent(in) :: kin,dx
        real(rkind) :: k
        real(rkind) :: T
    
        k = kin*dx
        T = (a90 + two* b90*COS(k) + two*c90*COS(two*k) +two*d90*COS(3._rkind*k) + two*e90*COS(4._rkind*k) ) &
          / (1._rkind + two*alpha90*COS(k) + two*beta90*COS(two*k) )
    
    end function
    
    pure elemental function GetCD10ModWaveNum(kin,dx) result(kp)
        use constants, only: two
        use cd10stuff, only: alpha10d1, beta10d1, a10d1, b10d1, c10d1
        real(rkind), intent(in) :: kin,dx
        real(rkind) :: k
        real(rkind) :: kp
    
        k = kin*dx
        kp = ( two*a10d1*sin(k) + two*b10d1*sin(two*k) + two*c10d1*sin(3._rkind*k) ) / (1._rkind + two*alpha10d1*cos(k) + two*beta10d1*cos(two*k))
    
    end function
    
    pure elemental function GetCD06ModWaveNum(kin,dx) result(kp)
        use constants, only: two
        use cd06stuff, only: alpha06d1, a06d1, b06d1
        real(rkind), intent(in) :: kin,dx
        real(rkind) :: k
        real(rkind) :: kp
    
        k = kin*dx
        kp = ( two*a06d1*sin(k) + two*b06d1*sin(two*k) ) / (1._rkind + two*alpha06d1*cos(k))
    
    end function


    subroutine init(this, gp, dx, dy, dz, method,isBoxFilt)
        class(poisson), target, intent(inout) :: this
        class(decomp_info), intent(in) :: gp
        real(rkind), intent(in) :: dx, dy, dz
        integer :: ierr
        logical, intent(in) :: isBoxFilt
        character(len=4), intent(in) :: method
        real(rkind), dimension(:,:,:), allocatable :: G1,G2,G3,mk1, mk2, mk3, mkabs_sq
        real(rkind), dimension(:,:,:), pointer :: k1k1, k1k2, k1k3, k2k2, k2k3, k3k3
        integer :: nx, ny, nz
        real(rkind), dimension(:,:,:), allocatable :: dummy

        nx = gp%xsz(1)
        ny = gp%ysz(2)
        nz = gp%zsz(3)

        allocate(this%FT)
        ierr = this%FT%init(nx,ny,nz,"x",dx, dy,dz, .false.)
       

        call this%FT%alloc_output(mk1)
        call this%FT%alloc_output(mk2)
        call this%FT%alloc_output(mk3)
        call this%FT%alloc_output(mkabs_sq)
        call this%FT%alloc_output(this%ustar_hat,3)
        call this%FT%alloc_output(this%u_hat,3)

        allocate(this%k_tensor(size(mk1,1),size(mk1,2),size(mk1,3),6))
       
        select case (method)
        case ("cd10")
            mk1 = GetCD10ModWaveNum(this%FT%k1,dx)
            mk2 = GetCD10ModWaveNum(this%FT%k2,dy)
            mk3 = GetCD10ModWaveNum(this%FT%k3,dz)
        case ("cd06")
            mk1 = GetCD06ModWaveNum(this%FT%k1,dx)
            mk2 = GetCD06ModWaveNum(this%FT%k2,dy)
            mk3 = GetCD06ModWaveNum(this%FT%k3,dz)
        case default 
            call GracefulExit("Incorrect method selected in poisson init",1031)
        end select 
        mkabs_sq = mk1*mk1 + mk2*mk2 + mk3*mk3
        where (mkabs_sq .lt. 1.d-18) 
            mkabs_sq = 1.d-10
        end where

        k1k1 => this%k_tensor(:,:,:,1)
        k1k2 => this%k_tensor(:,:,:,2)
        k1k3 => this%k_tensor(:,:,:,3)
        k2k2 => this%k_tensor(:,:,:,4)
        k2k3 => this%k_tensor(:,:,:,5)
        k3k3 => this%k_tensor(:,:,:,6)


        call this%FT%alloc_output(this%Gfilt)
        if (isBoxFilt) then
            call this%FT%alloc_output(G1)
            call this%FT%alloc_output(G2)
            call this%FT%alloc_output(G3)
            G1 = GetFilterTransferFunction(this%FT%k1, dx)
            G2 = GetFilterTransferFunction(this%FT%k2, dy)
            G3 = GetFilterTransferFunction(this%FT%k3, dz)
            this%Gfilt = G1*G2*G3
        else
            this%Gfilt = GetFilterTransferFunction(sqrt(this%FT%kabs_sq), dx)
        end if 
        this%nx_p_out = size(mk1,1) 
        this%ny_p_out = size(mk1,2) 
        this%nz_p_out = size(mk1,3)

        k1k1 = this%Gfilt*(one  - mk1*mk1/mkabs_sq) 
        k1k2 = this%Gfilt*(zero - mk1*mk2/mkabs_sq)
        k1k3 = this%Gfilt*(zero - mk1*mk3/mkabs_sq)
        k2k2 = this%Gfilt*(one  - mk2*mk2/mkabs_sq)
        k2k3 = this%Gfilt*(zero - mk2*mk3/mkabs_sq)
        k3k3 = this%Gfilt*(one  - mk3*mk3/mkabs_sq)
        
        
        deallocate (mk1, mk2, mk3, mkabs_sq, G1, G2, G3)

        allocate(this%myGP,source=gp)

        allocate(this%u_in_x(this%myGP%xsz(1),this%myGP%xsz(2),this%myGP%xsz(3),3))
    end subroutine 

    subroutine destroy(this)
        class(poisson), intent(inout) :: this

        if (allocated(this%myGP)) deallocate(this%myGP)
        if (allocated(this%k_tensor)) deallocate(this%k_tensor)
        if (allocated(this%Gfilt)) deallocate(this%Gfilt)
        if (allocated(this%ustar_hat)) deallocate(this%ustar_hat)
        if (allocated(this%u_in_x)) deallocate(this%u_in_x)
        if (allocated(this%FT)) then
            call this%FT%destroy
            deallocate(this%k_tensor)
        end if 

    end subroutine 

    subroutine pressureProjection(this,ustar,u)
        class(poisson), target, intent(inout) :: this
        real(rkind), dimension(this%myGP%ysz(1),this%myGP%ysz(2),this%myGP%ysz(3),3), intent(in) :: ustar
        real(rkind), dimension(this%myGP%ysz(1),this%myGP%ysz(2),this%myGP%ysz(3),3), intent(out) :: u

        real(rkind), dimension(:,:,:), pointer :: k1k1, k1k2, k1k3, k2k2, k2k3, k3k3



        k1k1 => this%k_tensor(:,:,:,1)
        k1k2 => this%k_tensor(:,:,:,2)
        k1k3 => this%k_tensor(:,:,:,3)
        k2k2 => this%k_tensor(:,:,:,4)
        k2k3 => this%k_tensor(:,:,:,5)
        k3k3 => this%k_tensor(:,:,:,6)
       


        ! First move us, vs, ws from y -> x
        call transpose_y_to_x(ustar(:,:,:,1),this%u_in_y(:,:,:,1),this%myGP)
        call transpose_y_to_x(ustar(:,:,:,2),this%u_in_y(:,:,:,2),this%myGP)
        call transpose_y_to_x(ustar(:,:,:,3),this%u_in_y(:,:,:,3),this%myGP)



        call this%FT%fft3_x2z(u_in_x(:,:,:,1),this%ustar_hat(:,:,:,1)) 
        call this%FT%fft3_x2z(u_in_x(:,:,:,2),this%ustar_hat(:,:,:,2)) 
        call this%FT%fft3_x2z(u_in_x(:,:,:,3),this%ustar_hat(:,:,:,3)) 
        
        this%u_hat(:,:,:,1) = k1k1*this%ustar_hat(:,:,:,1) + k1k2*this%ustar_hat(:,:,:,2) &
                                    +  k1k3*this%ustar_hat(:,:,:,3)
        
        this%u_hat(:,:,:,2) = k1k2*this%ustar_hat(:,:,:,1) + k2k2*this%ustar_hat(:,:,:,2) &
                                    +  k2k3*this%ustar_hat(:,:,:,3)

        this%u_hat(:,:,:,3) = k1k3*this%ustar_hat(:,:,:,1) + k2k3*this%ustar_hat(:,:,:,2) &
                                    +  k3k3*this%ustar_hat(:,:,:,3)

        call this%FT%ifft3_z2x(this%u_hat(:,:,:,1),u_in_x(:,:,:,1))
        call this%FT%ifft3_z2x(this%u_hat(:,:,:,2),u_in_x(:,:,:,2))
        call this%FT%ifft3_z2x(this%u_hat(:,:,:,3),u_in_x(:,:,:,3))

        ! Now transpose back from x -> y
        call transpose_x_to_y(u_in_x(:,:,:,1),u(:,:,:,1),this%myGP)
        call transpose_x_to_y(u_in_x(:,:,:,2),u(:,:,:,2),this%myGP)
        call transpose_x_to_y(u_in_x(:,:,:,3),u(:,:,:,3),this%myGP)

    end subroutine

end module 
