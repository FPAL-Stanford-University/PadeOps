module poissonMod
    use kind_parameters, only: rkind
    use decomp_2d, only: decomp_info
    type :: poisson
        private
        real(rkind), dimension(:,:,:,:), allocatable :: k_tensor 
        type(fft_3d), allocatable :: FT
        real(rkind), dimension(:,:,:,:), allocatable :: Gfilt
        complex(rkind), dimension(:,:,:,:), allocatable :: ustar_hat
        complex(rkind), dimension(:,:,:,:), allocatable :: u_hat

        integer :: nx_p_out, ny_p_out, nz_p_out
        integer :: nx_p_in, ny_p_in, nz_p_in

        logical :: isProjfromZ
        contains
            procedure :: init
            procedure :: destroy

            procedure :: pressureProjectionX
            procedure :: pressureProjectionZ
    end type

contains
#include "hitCD_files/meshgen.F90"
#include "hitCD_files/initfields.F90"


    subroutine init(this, gp, base_dec, dx, dy, dz, method,isBoxFilt)
        class(poisson), intent(inout) :: this
        class(decomp_info), intent(in) :: gp
        real(rkind), intent(in) :: dx, dy, dz
        character(len=1), intent(in) :: base_dec
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
        ierr = this%FT%init(nx,ny,nz,base_dec,dx, dy,dz, .false.)
       

        call this%FT%alloc_output(mk1)
        call this%FT%alloc_output(mk2)
        call this%FT%alloc_output(mk3)
        call this%FT%alloc_output(mkabs_sq)
        call this%FT%alloc_output(this%ustar_hat,3)
        call this%FT%alloc_output(this%u_hat,3)

        allocate(k_tensor(size(mk1,1),size(mk1,2),size(mk1,3),6))
       
        select case (method)
        case ("cd10")
            mk1 = GetCD10ModWaveNum(this%k1,dx)
            mk2 = GetCD10ModWaveNum(this%k2,dy)
            mk3 = GetCD10ModWaveNum(this%k3,dz)
        case ("cd06")
            mk1 = GetCD06ModWaveNum(this%k1,dx)
            mk2 = GetCD06ModWaveNum(this%k2,dy)
            mk3 = GetCD06ModWaveNum(this%k3,dz)
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
            G1 = GetFilterTransferFunction(this%k1, dx)
            G2 = GetFilterTransferFunction(this%k2, dy)
            G3 = GetFilterTransferFunction(this%k3, dz)
            this%Gfilt = G1*G2*G3
        else
            this%Gfilt = GetFilterTransferFunction(sqrt(this%kabs_sq), dx)
        end if 
        this%nx_p_out = size(this%mk1,1) 
        this%ny_p_out = size(this%mk1,2) 
        this%nz_p_out = size(this%mk1,3)

        k1k1 = this%Gfilt*(one  - mk1*mk1/mkabs_sq) 
        k1k2 = this%Gfilt*(zero - mk1*mk2/mkabs_sq)
        k1k3 = this%Gfilt*(zero - mk1*mk3/mkabs_sq)
        k2k2 = this%Gfilt*(one  - mk2*mk2/mkabs_sq)
        k2k3 = this%Gfilt*(zero - mk2*mk3/mkabs_sq)
        k3k3 = this%Gfilt*(one  - mk3*mk3/mkabs_sq)
        
        
        deallocate (mk1, mk2, mk3, mkabs_sq, G1, G2, G3)

        call this%FT%alloc_input(dummy)
        this%nx_p_in = size(dummy,1) 
        this%ny_p_in = size(dummy,2) 
        this%nz_p_in = size(dummy,3) 
        deallocate(dummy)


        if (base_dec == "z") then
            this%isProjfromZ = .true. 
        else
            this%isProjfromZ = .false. 
        end if 
    end subroutine 

    subroutine destroy(this)
        class(poisson), intent(inout) :: this

        if (allocated(this%k_tensor)) deallocate(this%k_tensor)
        if (allocated(this%Gfilt)) deallocate(this%Gfilt)
        if (allocated(this%ustar_hat)) deallocate(this%ustar_hat)
        if (allocated(this%FT)) then
            call this%FT%destroy
            deallocate(this%k_tensor)
        end if 

    end subroutine 

    subroutine pressureProjection(this,ustar,u)
        class(poisson), intent(inout) :: this
        real(rkind), dimension(this%nx_p_in,this%nx_p_in,this%nx_p_in), intent(in) :: ustar
        real(rkind), dimension(this%nx_p_in,this%nx_p_in,this%nx_p_in), intent(out) :: u

        real(rkind), dimension(:,:,:), pointer :: k1k1, k1k2, k1k3, k2k2, k2k3, k3k3



        k1k1 => this%k_tensor(:,:,:,1)
        k1k2 => this%k_tensor(:,:,:,2)
        k1k3 => this%k_tensor(:,:,:,3)
        k2k2 => this%k_tensor(:,:,:,4)
        k2k3 => this%k_tensor(:,:,:,5)
        k3k3 => this%k_tensor(:,:,:,6)
        
        if (this%isProjfromZ) then
            call this%FT%fft3_z2x(ustar(:,:,:,1),this%ustar_hat(:,:,:,1)) 
            call this%FT%fft3_z2x(ustar(:,:,:,2),this%ustar_hat(:,:,:,2)) 
            call this%FT%fft3_z2x(ustar(:,:,:,3),this%ustar_hat(:,:,:,3)) 
            
            this%u_hat(:,:,:,1) = k1k1*ustar_hat(:,:,:,1) + k1k2*ustar_hat(:,:,:,2) &
                                        +  k1k3*ustar_hat(:,:,:,3)
            
            this%u_hat(:,:,:,2) = k1k2*ustar_hat(:,:,:,1) + k2k2*ustar_hat(:,:,:,2) &
                                        +  k2k3*ustar_hat(:,:,:,3)

            this%u_hat(:,:,:,3) = k1k3*ustar_hat(:,:,:,1) + k2k3*ustar_hat(:,:,:,2) &
                                        +  k3k3*ustar_hat(:,:,:,3)

            call this%FT%ifft3_x2z(this%u_hat(:,:,:,1),u(:,:,:,1))
            call this%FT%ifft3_x2z(this%u_hat(:,:,:,2),u(:,:,:,2))
            call this%FT%ifft3_x2z(this%u_hat(:,:,:,3),u(:,:,:,3))

        else
            call this%FT%fft3_x2z(ustar(:,:,:,1),this%ustar_hat(:,:,:,1)) 
            call this%FT%fft3_x2z(ustar(:,:,:,2),this%ustar_hat(:,:,:,2)) 
            call this%FT%fft3_x2z(ustar(:,:,:,3),this%ustar_hat(:,:,:,3)) 
            
            this%u_hat(:,:,:,1) = k1k1*ustar_hat(:,:,:,1) + k1k2*ustar_hat(:,:,:,2) &
                                        +  k1k3*ustar_hat(:,:,:,3)
            
            this%u_hat(:,:,:,2) = k1k2*ustar_hat(:,:,:,1) + k2k2*ustar_hat(:,:,:,2) &
                                        +  k2k3*ustar_hat(:,:,:,3)

            this%u_hat(:,:,:,3) = k1k3*ustar_hat(:,:,:,1) + k2k3*ustar_hat(:,:,:,2) &
                                        +  k3k3*ustar_hat(:,:,:,3)

            call this%FT%ifft3_z2x(this%u_hat(:,:,:,1),u(:,:,:,1))
            call this%FT%ifft3_z2x(this%u_hat(:,:,:,2),u(:,:,:,2))
            call this%FT%ifft3_z2x(this%u_hat(:,:,:,3),u(:,:,:,3))

        end if 

    end subroutine


end module 
