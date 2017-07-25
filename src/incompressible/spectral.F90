module spectralMod
    use kind_parameters, only: rkind
    use decomp_2d, only: decomp_info, decomp_info_init, &
                    transpose_x_to_y, transpose_y_to_x, &
                    transpose_y_to_z, transpose_z_to_y, nrank 
    use decomp_2d_fft, only: decomp_2d_fft_init, decomp_2d_fft_finalize, decomp_2d_fft_get_size
    use exits, only: GracefulExit, message 
    use constants, only: pi, one, zero, two, three, four, eight 
    use fft_3d_stuff, only: fft_3d
    use mpi
    use reductions, only: p_sum 
    use numerics, only: use3by2rule
 
    implicit none
    private
    public :: spectral, GetWaveNums, useExhaustiveFFT 
    
    include "fftw3.f"

    logical :: useExhaustiveFFT = .true. 

    type :: spectral
        private
        real(rkind), dimension(:,:,:), allocatable, public :: k1, k2, k3, kabs_sq, k1_der2, k2_der2, k3_der2, one_by_kabs_sq
        integer, dimension(3) :: fft_start, fft_end, fft_size
        real(rkind), dimension(:,:,:), allocatable, public :: Gdealias, GtestFilt, arr1Up, arr2Up, GksPrep1, GksPrep2
        !real(rkind), dimension(:,:,:), allocatable :: k3inZ
        complex(rkind), dimension(:,:,:), allocatable :: fhatz, ctmpz
        integer :: rPencil ! Pencil dimension for the real input
        logical :: is3dFFT = .true. ! use 3d FFTs
        logical :: isInitialized = .false.
        logical :: fixOddball = .true. 
        integer, public :: nx_g, ny_g, nz_g
        integer :: nx_r,ny_r,nz_r
        integer :: nx_c,ny_c,nz_c
        type(fft_3d), allocatable :: FT
        logical :: use2decompFFT = .false. 
        logical :: useConsrvD2 = .true. 
        real(rkind) :: normfact, normfactZ 
        type(decomp_info), allocatable, public :: spectdecomp, physdecomp, dealiasdecomp
        logical                                :: StoreK = .true., init_periodicInZ = .false. 
        logical, public :: carryingZeroK = .false.
        integer, public :: zeroK_i = 123456, zeroK_j = 123456
        real(rkind), dimension(:,:), allocatable :: GsurfaceFilter 
        real(rkind) :: dealiasFact = 2.d0/3.d0

        logical, dimension(:,:,:), allocatable :: G_bandpass
        integer, dimension(:,:,:), allocatable :: G_PostProcess
        complex(rkind), dimension(:,:,:), pointer :: cbuffz_bp, cbuffy_bp
        
        logical :: BandPassFilterInitialized = .false. 
        logical :: initPostProcessor = .false.
        integer(kind=8) :: plan_c2c_fwd_z_oop
        integer(kind=8) :: plan_c2c_fwd_z_ip
        integer(kind=8) :: plan_c2c_bwd_z_oop
        integer(kind=8) :: plan_c2c_bwd_z_ip
        integer(kind=8) :: plan_r2c_z, plan_c2r_z 
        complex(rkind), dimension(:), allocatable :: k3_C2Eshift, k3_E2Cshift, E2Cshift, C2Eshift
        real(rkind), dimension(:), allocatable :: mk3sq

        contains
            procedure           :: init
            procedure           :: init_bandpass_filter
            procedure           :: destroy
            procedure, private  :: alloc_r2c_out_Rank3
            procedure, private  :: alloc_r2c_out_Rank4
            generic             :: alloc_r2c_out => alloc_r2c_out_Rank4, alloc_r2c_out_Rank3
            procedure           :: fft_y2z
            procedure           :: ifft_z2y
            procedure, private  :: take_fftz
            procedure, private  :: take_ifftz
            procedure, private  :: init_periodic_inZ_procedures 
            procedure           :: fft
            procedure           :: ifft
            procedure           :: fft1_x2y
            procedure, private  :: initializeEverything
            procedure           :: dealias
            procedure           :: dealias_edgeField
            procedure           :: mTimes_ik1_oop
            procedure           :: mTimes_ik1_ip
            procedure           :: mTimes_ik2_oop
            procedure           :: mTimes_ik2_ip
            procedure           :: TestFilter_ip
            procedure           :: TestFilter_oop
            procedure           :: SurfaceFilter_ip
            procedure           :: SurfaceFilter_oop
            procedure           :: dealiasedMult_ip
            procedure           :: dealiasedMult_oop
            procedure           :: dealiasedDiv_ip
            procedure           :: dealiasedDiv_oop
            procedure           :: dealiasedSquare_ip
            procedure           :: dealiasedSquare_oop
            procedure           :: KSprepFilter2
            procedure           :: KSprepFilter1
            procedure           :: bandpassFilter

            procedure           :: take_fft1d_z2z_ip
            procedure           :: take_ifft1d_z2z_ip
            procedure           :: shiftz_E2C
            procedure           :: shiftz_C2E
            procedure, private  :: ddz_E2C_spect_cmplx
            procedure, private  :: ddz_C2E_spect_cmplx
            procedure, private  :: d2dz2_C2C_spect_cmplx
            procedure, private  :: d2dz2_E2E_spect_cmplx
            procedure, private  :: interp_E2C_spect_cmplx
            procedure, private  :: interp_C2E_spect_cmplx

            procedure, private  :: ddz_E2C_spect_real
            procedure, private  :: ddz_C2E_spect_real
            procedure, private  :: d2dz2_C2C_spect_real
            procedure, private  :: d2dz2_E2E_spect_real
            procedure, private  :: interp_E2C_spect_real
            procedure, private  :: interp_C2E_spect_real
            
            generic             :: ddz_E2C_spect => ddz_E2C_spect_cmplx, ddz_E2C_spect_real
            generic             :: ddz_C2E_spect => ddz_C2E_spect_cmplx, ddz_C2E_spect_real
            generic             :: interp_E2C_spect => interp_E2C_spect_cmplx, interp_E2C_spect_real
            generic             :: interp_C2E_spect => interp_C2E_spect_cmplx, interp_C2E_spect_real
            generic             :: d2dz2_C2C_spect => d2dz2_C2C_spect_cmplx, d2dz2_C2C_spect_real
            generic             :: d2dz2_E2E_spect => d2dz2_E2E_spect_cmplx, d2dz2_E2E_spect_real

            procedure           :: initPP
            procedure           :: destroyPP
            procedure           :: destroy_bandpass_filter
            procedure           :: spectralFilter_ip
            !procedure, private  :: upsample_Fhat
            !procedure, private  :: downsample_Fhat

    end type

contains

      subroutine bandpassFilter(this, uhat, uFilt) 
         class(spectral),  intent(inout)         :: this
         complex(rkind), dimension(this%spectdecomp%ysz(1),this%spectdecomp%ysz(2), this%spectdecomp%ysz(3)), intent(in) :: uhat 
         real(rkind)   , dimension(this%physdecomp%xsz(1),this%physdecomp%xsz(2), this%physdecomp%xsz(3)), intent(out) :: uFilt

         call transpose_y_to_z(uhat, this%cbuffz_bp, this%spectdecomp)
         call this%take_fft1d_z2z_ip(this%cbuffz_bp)

         where (this%G_bandpass == .false. ) 
            this%cbuffz_bp = zero
         end where

         call this%take_ifft1d_z2z_ip(this%cbuffz_bp)
         call transpose_z_to_y(this%cbuffz_bp, this%cbuffy_bp, this%spectdecomp)
         call this%ifft(this%cbuffy_bp, uFilt)

      end subroutine 

      subroutine init_bandpass_filter(this, kleft, kright, cbuffz, cbuffy)
         class(spectral),  intent(inout)         :: this
         real(rkind), intent(in) :: kleft, kright
         complex(rkind), dimension(this%spectdecomp%ysz(1),this%spectdecomp%ysz(2), this%spectdecomp%ysz(3)), intent(in), target :: cbuffy
         complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2), this%spectdecomp%zsz(3)), intent(in), target :: cbuffz
         real(rkind), dimension(:,:,:), allocatable :: rbuffz1, rbuffz2
        
         this%cbuffz_bp => cbuffz
         this%cbuffy_bp => cbuffy
         allocate(rbuffz1(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)))
         allocate(rbuffz2(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)))
         allocate(this%G_bandpass(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)))

         call transpose_y_to_z(this%k1, rbuffz1, this%spectdecomp)
         rbuffz2 = rbuffz1**2 

         call transpose_y_to_z(this%k2, rbuffz1, this%spectdecomp)
         rbuffz2 = rbuffz2 + rbuffz1**2

         call transpose_y_to_z(this%k3, rbuffz1, this%spectdecomp)
         rbuffz2 = rbuffz2 + rbuffz1**2

         rbuffz2 = sqrt(rbuffz2)

         this%G_bandpass = .true. 
         where (rbuffz2 < kleft)
            this%G_bandpass = .false. 
         end where

         where (rbuffz2 > kright) 
            this%G_bandpass = .false. 
         end where

         this%BandPassFilterInitialized = .true. 
         deallocate(rbuffz1, rbuffz2)
         call message(0, "Band pass filter initialized")
      end subroutine

    subroutine destroy_bandpass_filter(this)
        class(spectral),  intent(inout)         :: this

        if (this%BandPassFilterInitialized) then
            deallocate(this%G_bandpass)
        end if

    end subroutine

    subroutine KSprepFilter1(this, finout)
        class(spectral),  intent(in)         :: this
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(inout) :: finout

        finout = this%GksPrep1*finout
    end subroutine

    subroutine KSprepFilter2(this, finout)
        class(spectral),  intent(in)         :: this
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(inout) :: finout

        finout = this%GksPrep2*finout
    end subroutine

    subroutine mTimes_ik1_oop(this, fin, fout)
        class(spectral),  intent(in)         :: this
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(in) :: fin
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(out) :: fout
        integer :: i, j, k
        real(rkind) :: rpart, ipart
       
        do k = 1,this%fft_size(3)
            do j = 1,this%fft_size(2)
                !$omp simd 
                do i = 1,this%fft_size(1)
                    rpart = -this%k1(i,j,k)*dimag(fin(i,j,k))
                    ipart = this%k1(i,j,k)*dreal(fin(i,j,k))
                    fout(i,j,k) = dcmplx(rpart,ipart)
                end do 
            end do 
        end do 

    end subroutine
    
    subroutine mTimes_ik2_oop(this, fin, fout)
        class(spectral),  intent(in)         :: this
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(in) :: fin
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(out) :: fout
        integer :: i, j, k
        real(rkind) :: rpart, ipart
        
        do k = 1,this%fft_size(3)
            do j = 1,this%fft_size(2)
                !$omp simd 
                do i = 1,this%fft_size(1)
                    rpart = -this%k2(i,j,k)*dimag(fin(i,j,k))
                    ipart = this%k2(i,j,k)*dreal(fin(i,j,k))
                    fout(i,j,k) = dcmplx(rpart,ipart)
                end do 
            end do 
        end do 

    end subroutine


    subroutine mTimes_ik1_ip(this, fin)
        class(spectral),  intent(in)         :: this
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(inout) :: fin
        integer :: i, j, k
        real(rkind) :: rpart, ipart
        
        do k = 1,this%fft_size(3)
            do j = 1,this%fft_size(2)
                !$omp simd 
                do i = 1,this%fft_size(1)
                    rpart = -this%k1(i,j,k)*dimag(fin(i,j,k))
                    ipart = this%k1(i,j,k)*dreal(fin(i,j,k))
                    fin(i,j,k) = dcmplx(rpart,ipart)
                end do 
            end do 
        end do 

    end subroutine

    subroutine mTimes_ik2_ip(this, fin)
        class(spectral),  intent(in)         :: this
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(inout) :: fin
        integer :: i, j, k
        real(rkind) :: rpart, ipart
        
        do k = 1,this%fft_size(3)
            do j = 1,this%fft_size(2)
                !$omp simd 
                do i = 1,this%fft_size(1)
                    rpart = -this%k2(i,j,k)*dimag(fin(i,j,k))
                    ipart = this%k2(i,j,k)*dreal(fin(i,j,k))
                    fin(i,j,k) = dcmplx(rpart,ipart)
                end do 
            end do 
        end do 

    end subroutine

    subroutine dealias(this, fhat)
        class(spectral),  intent(inout)         :: this
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(inout) :: fhat
        integer :: i, j, k

        if (this%init_periodicinZ) then
            call this%take_fftz(fhat)
            do k = 1,this%spectdecomp%zsz(3)
               do j = 1,this%spectdecomp%zsz(2)
                  !$omp simd
                  do i = 1,this%spectdecomp%zsz(1)
                     this%ctmpz(i,j,k) = this%ctmpz(i,j,k)*this%Gdealias(i,j,k)
                  end do 
               end do 
            end do 
            call this%take_ifftz(fhat)
        else
         do k = 1,this%fft_size(3)
             do j = 1,this%fft_size(2)
                 !$omp simd 
                 do i = 1,this%fft_size(1)
                     fhat(i,j,k) = fhat(i,j,k)*this%Gdealias(i,j,k)
                 end do 
             end do 
         end do 
        end if 
         
    end subroutine
  
    subroutine dealias_edgeField(this, fhat)
        class(spectral),  intent(inout)         :: this
        complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)+1), intent(inout) :: fhat
        integer :: i, j, k

        if (this%init_periodicInZ) then
            call dfftw_execute_dft(this%plan_c2c_fwd_z_ip, fhat(:,:,1:this%spectdecomp%zsz(3)) , fhat(:,:,1:this%spectdecomp%zsz(3)))  
            do k = 1,this%spectdecomp%zsz(3)
               do j = 1,this%spectdecomp%zsz(2)
                  !$omp simd
                  do i = 1,this%spectdecomp%zsz(1)
                     fhat(i,j,k) = fhat(i,j,k)*this%Gdealias(i,j,k)
                  end do 
               end do 
            end do 
            call dfftw_execute_dft(this%plan_c2c_bwd_z_ip, fhat(:,:,1:this%spectdecomp%zsz(3)) , fhat(:,:,1:this%spectdecomp%zsz(3)))  
            fhat(:,:,1:this%spectdecomp%zsz(3)) = this%normfactz*fhat(:,:,1:this%spectdecomp%zsz(3))
            fhat(:,:,this%spectdecomp%zsz(3)+1) = fhat(:,:,1)
        end if 

    end subroutine 
   
    subroutine ddz_E2C_spect_real(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)+1), intent(in)  :: arr_in
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)),   intent(out) :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft_r2c(this%plan_r2c_z, arr_in(:,:,1:this%physdecomp%zsz(3)), this%fhatz)
         do k = 1,this%physdecomp%zsz(3)/2 ! Note that the oddball is ignored 
            do j = 1,this%physdecomp%zsz(2)
               !$omp simd
               do i = 1,this%physdecomp%zsz(1)
                  this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%k3_E2Cshift(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft_c2r(this%plan_c2r_z, this%fhatz, arr_out)
         arr_out = arr_out*this%normfactz
      end if

    end subroutine 
    
    subroutine ddz_E2C_spect_cmplx(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)+1), intent(in)  :: arr_in
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)),   intent(out) :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft(this%plan_c2c_fwd_z_oop, arr_in(:,:,1:this%spectdecomp%zsz(3)), arr_out(:,:,1:this%spectdecomp%zsz(3)))
         do k = 1,this%spectdecomp%zsz(3)
            do j = 1,this%spectdecomp%zsz(2)
               !$omp simd
               do i = 1,this%spectdecomp%zsz(1)
                  arr_out(i,j,k) = arr_out(i,j,k)*this%k3_E2Cshift(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft(this%plan_c2c_bwd_z_ip, arr_out(:,:,1:this%spectdecomp%zsz(3)), arr_out(:,:,1:this%spectdecomp%zsz(3)))
         arr_out = arr_out*this%normfactz
      end if

    end subroutine 

    subroutine shiftz_E2C(this, arr_inout)
      class(spectral), intent(in) :: this
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)), intent(inout)  :: arr_inout
      integer :: i, j, k

      do k = 1,this%spectdecomp%zsz(3) 
         do j = 1,this%spectdecomp%zsz(2)
            !$omp simd
            do i = 1,this%spectdecomp%zsz(1)
               arr_inout(i,j,k) = arr_inout(i,j,k)*this%E2Cshift(k) 
            end do 
         end do 
      end do
    end subroutine 

    subroutine shiftz_C2E(this, arr_inout)
      class(spectral), intent(in) :: this
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)), intent(inout)  :: arr_inout
      integer :: i, j, k

      do k = 1,this%spectdecomp%zsz(3) 
         do j = 1,this%spectdecomp%zsz(2)
            !$omp simd
            do i = 1,this%spectdecomp%zsz(1)
               arr_inout(i,j,k) = arr_inout(i,j,k)*this%C2Eshift(k) 
            end do 
         end do 
      end do
    end subroutine 

    subroutine interp_E2C_spect_real(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)+1), intent(in)  :: arr_in
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)),   intent(out) :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft_r2c(this%plan_r2c_z, arr_in(:,:,1:this%physdecomp%zsz(3)), this%fhatz)
         do k = 1,this%physdecomp%zsz(3)/2 ! Note that the oddball is ignored 
            do j = 1,this%physdecomp%zsz(2)
               !$omp simd
               do i = 1,this%physdecomp%zsz(1)
                  this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%E2Cshift(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft_c2r(this%plan_c2r_z, this%fhatz, arr_out)
         arr_out = arr_out*this%normfactz
      end if

    end subroutine 


    subroutine interp_E2C_spect_cmplx(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)+1), intent(in)  :: arr_in
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)),   intent(out) :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft(this%plan_c2c_fwd_z_oop, arr_in(:,:,1:this%spectdecomp%zsz(3)), arr_out(:,:,1:this%spectdecomp%zsz(3)))
         do k = 1,this%spectdecomp%zsz(3)
            do j = 1,this%spectdecomp%zsz(2)
               !$omp simd
               do i = 1,this%spectdecomp%zsz(1)
                  arr_out(i,j,k) = arr_out(i,j,k)*this%E2Cshift(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft(this%plan_c2c_bwd_z_ip, arr_out(:,:,1:this%spectdecomp%zsz(3)), arr_out(:,:,1:this%spectdecomp%zsz(3)))
         arr_out = arr_out*this%normfactz
      end if

    end subroutine 

    subroutine ddz_C2E_spect_real(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)  ), intent(in)  :: arr_in
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)+1), intent(out) :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft_r2c(this%plan_r2c_z, arr_in, this%fhatz)
         do k = 1,this%physdecomp%zsz(3)/2 ! Note that the oddball is ignored 
            do j = 1,this%physdecomp%zsz(2)
               !$omp simd
               do i = 1,this%physdecomp%zsz(1)
                  this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%k3_C2Eshift(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft_c2r(this%plan_c2r_z, this%fhatz, arr_out(:,:,1:this%physdecomp%zsz(3)))
         arr_out(:,:,1:this%physdecomp%zsz(3)) = arr_out(:,:,1:this%physdecomp%zsz(3))*this%normfactz
         arr_out(:,:,this%physdecomp%zsz(3)+1) = arr_out(:,:,1)
      end if

    end subroutine 

    subroutine ddz_C2E_spect_cmplx(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)  ), intent(in )  :: arr_in
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)+1), intent(out) :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft(this%plan_c2c_fwd_z_oop, arr_in, this%ctmpz)
         do k = 1,this%spectdecomp%zsz(3)
            do j = 1,this%spectdecomp%zsz(2)
               !$omp simd
               do i = 1,this%spectdecomp%zsz(1)
                  this%ctmpz(i,j,k) = this%ctmpz(i,j,k)*this%k3_C2Eshift(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft(this%plan_c2c_bwd_z_oop, this%ctmpz, arr_out(:,:,1:this%spectdecomp%zsz(3)))
         arr_out(:,:,1:this%spectdecomp%zsz(3)) =  this%normfactz*arr_out(:,:,1:this%spectdecomp%zsz(3))
         arr_out(:,:,this%spectdecomp%zsz(3)+1) =  arr_out(:,:,1)
      end if

    end subroutine 

    subroutine interp_C2E_spect_real(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)  ), intent(in)  :: arr_in
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)+1), intent(out) :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft_r2c(this%plan_r2c_z, arr_in, this%fhatz)
         do k = 1,this%physdecomp%zsz(3)/2 ! Note that the oddball is ignored 
            do j = 1,this%physdecomp%zsz(2)
               !$omp simd
               do i = 1,this%physdecomp%zsz(1)
                  this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%C2Eshift(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft_c2r(this%plan_c2r_z, this%fhatz, arr_out(:,:,1:this%physdecomp%zsz(3)))
         arr_out(:,:,1:this%physdecomp%zsz(3)) = arr_out(:,:,1:this%physdecomp%zsz(3))*this%normfactz
         arr_out(:,:,this%physdecomp%zsz(3)+1) = arr_out(:,:,1)
      end if

    end subroutine 
    
    subroutine interp_C2E_spect_cmplx(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)  ), intent(in )  :: arr_in
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)+1), intent(out) :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft(this%plan_c2c_fwd_z_oop, arr_in, this%ctmpz)
         do k = 1,this%spectdecomp%zsz(3)
            do j = 1,this%spectdecomp%zsz(2)
               !$omp simd
               do i = 1,this%spectdecomp%zsz(1)
                  this%ctmpz(i,j,k) = this%ctmpz(i,j,k)*this%C2Eshift(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft(this%plan_c2c_bwd_z_oop, this%ctmpz, arr_out(:,:,1:this%spectdecomp%zsz(3)))
         arr_out(:,:,1:this%spectdecomp%zsz(3)) =  this%normfactz*arr_out(:,:,1:this%spectdecomp%zsz(3))
         arr_out(:,:,this%spectdecomp%zsz(3)+1) =  arr_out(:,:,1)
      end if

    end subroutine 

    subroutine d2dz2_C2C_spect_cmplx(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)), intent(in )  :: arr_in
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)), intent(out) :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft(this%plan_c2c_fwd_z_oop, arr_in, arr_out)
         do k = 1,this%spectdecomp%zsz(3)
            do j = 1,this%spectdecomp%zsz(2)
               !$omp simd
               do i = 1,this%spectdecomp%zsz(1)
                  arr_out(i,j,k) = arr_out(i,j,k)*this%mk3sq(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft(this%plan_c2c_bwd_z_ip,arr_out, arr_out)
         arr_out = arr_out*this%normfactz
      end if 
    end subroutine 
    
    subroutine d2dz2_C2C_spect_real(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)), intent(in )  :: arr_in
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)), intent(out)  :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft_r2c(this%plan_r2c_z, arr_in, this%fhatz)
         do k = 1,this%physdecomp%zsz(3)/2
            do j = 1,this%physdecomp%zsz(2)
               !$omp simd
               do i = 1,this%physdecomp%zsz(1)
                  this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%mk3sq(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft_c2r(this%plan_c2r_z,this%fhatz, arr_out)
         arr_out = arr_out*this%normfactz
      end if 
    end subroutine 
    
    subroutine d2dz2_E2E_spect_cmplx(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)+1), intent(in )  :: arr_in
      complex(rkind), dimension(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)+1), intent(out) :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft(this%plan_c2c_fwd_z_oop, arr_in(:,:,1:this%spectdecomp%zsz(3)), arr_out(:,:,1:this%spectdecomp%zsz(3)))
         do k = 1,this%spectdecomp%zsz(3)
            do j = 1,this%spectdecomp%zsz(2)
               !$omp simd
               do i = 1,this%spectdecomp%zsz(1)
                  arr_out(i,j,k) = arr_out(i,j,k)*this%mk3sq(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft(this%plan_c2c_bwd_z_ip,arr_out(:,:,1:this%spectdecomp%zsz(3)), arr_out(:,:,1:this%spectdecomp%zsz(3)))
         arr_out(:,:,1:this%spectdecomp%zsz(3)) = arr_out(:,:,1:this%spectdecomp%zsz(3))*this%normfactz
         arr_out(:,:,this%spectdecomp%zsz(3)+1) = arr_out(:,:,1)
      end if 
    end subroutine 

    subroutine d2dz2_E2E_spect_real(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)+1), intent(in )  :: arr_in
      real(rkind), dimension(this%physdecomp%zsz(1),this%physdecomp%zsz(2),this%physdecomp%zsz(3)+1), intent(out)  :: arr_out
      integer :: i, j, k

      if (this%init_periodicInZ) then
         call dfftw_execute_dft_r2c(this%plan_r2c_z, arr_in(:,:,1:this%physdecomp%zsz(3)), this%fhatz)
         do k = 1,this%physdecomp%zsz(3)/2
            do j = 1,this%physdecomp%zsz(2)
               !$omp simd
               do i = 1,this%physdecomp%zsz(1)
                  this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%mk3sq(k) 
               end do 
            end do 
         end do
         call dfftw_execute_dft_c2r(this%plan_c2r_z,this%fhatz, arr_out(:,:,1:this%physdecomp%zsz(3)))
         arr_out(:,:,1:this%physdecomp%zsz(3)) = arr_out(:,:,1:this%physdecomp%zsz(3))*this%normfactz
         arr_out(:,:,this%physdecomp%zsz(3)+1) = arr_out(:,:,1)
      end if 
    end subroutine 

    subroutine TestFilter_ip(this, fhat)
        class(spectral),  intent(in)         :: this
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(inout) :: fhat
        integer :: i, j, k

        do k = 1,this%fft_size(3)
            do j = 1,this%fft_size(2)
                !$omp simd 
                do i = 1,this%fft_size(1)
                    fhat(i,j,k) = fhat(i,j,k)*this%GTestFilt(i,j,k)
                end do 
            end do 
        end do 
         
    end subroutine

    subroutine SurfaceFilter_ip(this, fhat)
        class(spectral), intent(in) :: this
        complex(rkind), dimension(this%spectdecomp%zsz(1), this%spectdecomp%zsz(2)), intent(inout) :: fhat
        fhat = fhat * this%GsurfaceFilter
    end subroutine 

    pure subroutine SurfaceFilter_oop(this, fin, fout)
        class(spectral), intent(in) :: this
        complex(rkind), dimension(this%spectdecomp%zsz(1), this%spectdecomp%zsz(2)), intent(in) :: fin
        complex(rkind), dimension(this%spectdecomp%zsz(1), this%spectdecomp%zsz(2)), intent(out) :: fout

        fout = fin * this%GsurfaceFilter
    end subroutine 
    
    subroutine TestFilter_oop(this, fhat,fhatout)
        class(spectral),  intent(in)         :: this
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(in) :: fhat
        complex(rkind), dimension(this%fft_size(1),this%fft_size(2),this%fft_size(3)), intent(out) :: fhatout
        integer :: i, j, k

        do k = 1,this%fft_size(3)
            do j = 1,this%fft_size(2)
                !$omp simd 
                do i = 1,this%fft_size(1)
                    fhatout(i,j,k) = fhat(i,j,k)*this%GTestFilt(i,j,k)
                end do 
            end do 
        end do 
         
    end subroutine

      subroutine init_periodic_inZ_procedures(this, dx, dy, dz)
         use constants, only: imi
         class(spectral),  intent(inout)         :: this
         real(rkind), intent(in) :: dx, dy, dz
         real(rkind), dimension(:,:,:),  allocatable :: rbuffz, rbuffz1
         real(rkind), dimension(:), allocatable :: k3_1d
         complex(rkind), dimension(:,:,:), allocatable :: cbuffz
         real(rkind) :: kdealiasx, kdealiasy, kdealiasz
         integer :: nz, nxT, nyT

         if (.not. this%isinitialized) then
            call gracefulExit("You cannot call INIT_PERIODIC_INZ_PROCEDURES before initializing the spectral derived type",103) 
         end if 

         deallocate(this%Gdealias)
         
         nxT = this%spectdecomp%zsz(1)
         nyT = this%spectdecomp%zsz(2)
         nz = this%spectdecomp%zsz(3)
         if (mod(nz,2) .ne. 0) then
            call gracefulExit("You cannot initialize a periodic_inZ spectral type with an odd values nz",104)
         end if
         
         allocate(this%Gdealias(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)))
         allocate(this%ctmpz(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)))
         allocate(rbuffz(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)))
         allocate(cbuffz(this%spectdecomp%zsz(1),this%spectdecomp%zsz(2),this%spectdecomp%zsz(3)))

         this%Gdealias = one 

         kdealiasx = (this%dealiasFact*pi/dx)
         kdealiasy = (this%dealiasFact*pi/dy)
         kdealiasz = (this%dealiasFact*pi/dz)

         call transpose_y_to_z(this%k1, rbuffz, this%spectdecomp)
         where (abs(rbuffz) >= kdealiasx)    
            this%Gdealias = zero
         end where
      

         call transpose_y_to_z(this%k2, rbuffz, this%spectdecomp)
         where (abs(rbuffz) >= kdealiasy)    
            this%Gdealias = zero
         end where

         call transpose_y_to_z(this%k3, rbuffz, this%spectdecomp)
         where (abs(rbuffz) >= kdealiasz)
            this%Gdealias = zero 
         end where
         
         call dfftw_plan_many_dft(this%plan_c2c_fwd_z_ip, 1, nz,nxT*nyT, this%ctmpz, nz, &
                      nxT*nyT, 1, this%ctmpz, nz,nxT*nyT, 1, FFTW_FORWARD , FFTW_EXHAUSTIVE)   
         
         call dfftw_plan_many_dft(this%plan_c2c_fwd_z_oop, 1, nz,nxT*nyT, cbuffz, nz, &
                      nxT*nyT, 1, this%ctmpz, nz,nxT*nyT, 1, FFTW_FORWARD , FFTW_EXHAUSTIVE)   
         
         call dfftw_plan_many_dft(this%plan_c2c_bwd_z_ip, 1, nz,nxT*nyT, this%ctmpz, nz, &
                      nxT*nyT, 1, this%ctmpz, nz,nxT*nyT, 1, FFTW_BACKWARD , FFTW_EXHAUSTIVE)   
         
         call dfftw_plan_many_dft(this%plan_c2c_bwd_z_oop, 1, nz,nxT*nyT, cbuffz, nz, &
                      nxT*nyT, 1, this%ctmpz, nz,nxT*nyT, 1, FFTW_BACKWARD , FFTW_EXHAUSTIVE)   

         allocate(rbuffz1(this%physdecomp%zsz(1), this%physdecomp%zsz(2), this%physdecomp%zsz(3)      ))
         allocate(this%fhatz(this%physdecomp%zsz(1), this%physdecomp%zsz(2), this%physdecomp%zsz(3)/2 + 1))
                      
         call dfftw_plan_many_dft_r2c(this%plan_r2c_z, 1, this%physdecomp%zsz(3), &
                 this%physdecomp%zsz(1)*this%physdecomp%zsz(2),rbuffz1 , this%physdecomp%zsz(3), &
                 this%physdecomp%zsz(1)*this%physdecomp%zsz(2), 1, this%fhatz, this%physdecomp%zsz(3)/2 + 1, &
                 this%physdecomp%zsz(1)*this%physdecomp%zsz(2), 1, FFTW_EXHAUSTIVE)

         call dfftw_plan_many_dft_c2r(this%plan_c2r_z, 1, this%physdecomp%zsz(3), &
                 this%physdecomp%zsz(1)*this%physdecomp%zsz(2), this%fhatz , this%physdecomp%zsz(3)/2 + 1, &
                 this%physdecomp%zsz(1)*this%physdecomp%zsz(2), 1, rbuffz1, this%physdecomp%zsz(3), &
                 this%physdecomp%zsz(1)*this%physdecomp%zsz(2), 1, FFTW_EXHAUSTIVE)

         allocate(this%k3_C2Eshift(nz))
         allocate(this%k3_E2Cshift(nz))
         allocate(this%C2Eshift(nz))
         allocate(this%E2Cshift(nz))
         allocate(this%mk3sq(nz))
         allocate(k3_1d(nz))
         k3_1d = GetWaveNums(nz,dz) 
         this%mk3sq = -(k3_1d**2)
         this%k3_C2Eshift = imi*k3_1d*exp(-imi*k3_1d*dz/two)
         this%k3_E2Cshift = imi*k3_1d*exp( imi*k3_1d*dz/two)
         this%C2Eshift = exp(-imi*k3_1d*dz/two)
         this%E2Cshift = exp( imi*k3_1d*dz/two)

         deallocate(k3_1d)
         this%normfactZ = one/real(nz,rkind)

         this%fhatz = cmplx(zero,zero,rkind)
         this%init_periodicinZ = .true. 

         deallocate(rbuffz, cbuffz, rbuffz1)
      end subroutine 

    subroutine init(this,pencil, nx_g, ny_g, nz_g, dx, dy, dz, scheme, filt, dimTransform, fixOddball, use2decompFFT, useConsrvD2, createK, exhaustiveFFT, init_periodicInZ, dealiasF) 
        class(spectral),  intent(inout)         :: this
        character(len=1), intent(in)            :: pencil              ! PHYSICAL decomposition direction
        integer,          intent(in)            :: nx_g, ny_g, nz_g    ! Global data size
        real(rkind),      intent(in)            :: dx, dy, dz          ! PHYSICAL grid spacing 
        character(len=*), intent(in)            :: scheme              ! Scheme used for modified wavenumbers
        character(len=*), intent(in)            :: filt                ! Scheme used for modified wavenumbers
        integer,          intent(in), optional  :: dimTransform        ! 2 or 3 (number of periodic directions) - default is 3
        logical,                      optional  :: fixOddball          ! Fix the oddball wavenumber - default is TRUE
        logical,                      optional  :: use2decompFFT       ! Use the 3d fft procedures defined in 2decomp - default is TRUE
        logical,                      optional  :: useConsrvD2         ! Use the conservative form of 2nd derivative - default is TRUE
        logical, intent(in),          optional  :: createK 
        logical, intent(in),          optional  :: exhaustiveFFT 
        logical, intent(in),          optional  :: init_periodicInZ 
        real(rkind), intent(in),      optional  :: dealiasF

        if (present(createK)) then
            this%storeK = createK
        end if 

        if (this%isInitialized) then
            call GracefulExit("You are trying to reinitialize SPECTRAL derived type. This is not allowed", 111)
        end if
        
        if (present(exhaustiveFFT)) then
            useExhaustiveFFT = exhaustiveFFT
        end if

        this%nx_g = nx_g
        this%ny_g = ny_g
        this%nz_g = nz_g

        if (present(dealiasF)) this%dealiasFact = dealiasF
        if (present(fixOddball)) this%fixOddball = fixOddball
        if (present(use2decompFFT)) this%use2decompFFT = use2decompFFT 
        if (present(useConsrvD2)) this%useConsrvD2 = useConsrvD2

        if (present(dimTransform)) then
            if (dimTransform == 2) then
                this%is3dFFT = .false.
                call this%initializeEverything(pencil,nx_g,ny_g,nz_g,dx,dy,dz,scheme,filt,.true.)
            else if (dimTransform == 3) then
                call this%initializeEverything(pencil,nx_g,ny_g,nz_g,dx,dy,dz,scheme,filt)
            else
                call GracefulExit("Incorrect choice for DIMTRANSFORM while initializing SPECTRAL derived type. &
                                   & Available options include 2 and 3", 102)
            end if 
        else
            call this%initializeEverything(pencil,nx_g,ny_g,nz_g,dx,dy,dz,scheme,filt)
        end if 
      
        this%normfact = one/real(nx_g)/real(ny_g)/real(nz_g) 
        this%isInitialized = .true.  

        if (present(init_periodicInZ)) then
            if (init_periodicInZ) then
               call this%init_periodic_inZ_procedures(dx, dy, dz)
            end if
        end if 
    end subroutine 


    subroutine initializeEverything(this,pencil,nx_g,ny_g,nz_g,dx,dy,dz,scheme,filt,nonPeriodic)
        class(spectral),  intent(inout), target :: this
        character(len=1), intent(in)            :: pencil              ! PHYSICAL decomposition direction
        integer,          intent(in)            :: nx_g, ny_g, nz_g    ! Global data size
        real(rkind),      intent(in)            :: dx, dy, dz          ! PHYSICAL grid spacing 
        character(len=*), intent(in)            :: scheme              ! Scheme used for modified wavenumbers
        character(len=*), intent(in)            :: filt                ! Scheme used for modified wavenumbers

        real(rkind), dimension(nx_g) :: k1_1d 
        real(rkind), dimension(ny_g) :: k2_1d 
        real(rkind), dimension(nz_g) :: k3_1d 
        type(decomp_info), allocatable :: spectdecomp
        real(rkind), dimension(:,:,:), allocatable :: tmp1, tmp2
        integer :: i, j, k, rPencil, ierr  
        !real(rkind), dimension(:,:,:), allocatable :: k1four, k2four, k3four
        logical, optional, intent(in) :: nonPeriodic
        logical :: TwoPeriodic = .false. 
        real(rkind) ::kdealiasx, kdealiasy 

        if (present(nonPeriodic)) then
                if (nonPeriodic) TwoPeriodic = .true.
        end if 

        ! STEP 0: Figure out what the input decomposition is going to be  
        select case (pencil)
        case("x")
            rPencil = 1
        case("z")
            rPencil = 3
        case("y")
            call GracefulExit("Y - input is currently unsupported",100) 
        case default
            call GracefulExit("Incorrect pencil direction chosen during initializing SPECTRAL derived type",101)
        end select

        this%rPencil = rPencil

        call message("===========================================================")
        call message("Now Generating the SPECTRAL - Derived Type for the problem.")

        ! STEP 1: Start the 3d FFT engine
        call message ( " ***** Using the FFTW (version 3.x) engine &
                & (found in: src/utilities/fft_3d.F90)  ***** ")
        allocate(this%FT)
        ierr = this%FT%init(nx_g,ny_g,nz_g,pencil, dx, &
            dy,dz, useExhaustiveFFT, .false., .false.)  
        if (ierr == 0) then
            call message (3, "Successfully initialized!")
        else
            call GracefulExit("Couldn't initialize 3d FFT inside SPECTRAL derived type",123)
        end if   
        this%use2decompFFT = .false. 

        ! STEP 2: Allocate wavenumbers and temporary decomp for spectral transposes
        allocate(spectdecomp)
        allocate(this%spectdecomp)
        allocate(this%physdecomp)
        call decomp_info_init(this%nx_g, this%ny_g, this%nz_g, this%physdecomp)
        select case (rPencil)
        case(1)
           call decomp_info_init(nx_g/2+1, ny_g, nz_g, spectdecomp)
           call decomp_info_init(nx_g/2+1, ny_g, nz_g, this%spectdecomp)
           if (TwoPeriodic) then
                this%fft_size(1) = spectdecomp%ysz(1)    
                this%fft_size(2) = spectdecomp%ysz(2)    
                this%fft_size(3) = spectdecomp%ysz(3) 
           else 
                this%fft_size(1) = spectdecomp%zsz(1)    
                this%fft_size(2) = spectdecomp%zsz(2)    
                this%fft_size(3) = spectdecomp%zsz(3) 
           end if    
        case(3)
           call decomp_info_init(nx_g, ny_g, nz_g/2+1, spectdecomp) 
           call decomp_info_init(nx_g, ny_g, nz_g/2+1, this%spectdecomp)
           this%fft_size(1) = spectdecomp%xsz(1)    
           this%fft_size(2) = spectdecomp%xsz(2)    
           this%fft_size(3) = spectdecomp%xsz(3)    
        end select 
       
        if (this%storeK) then 

            if (allocated(this%k1)) deallocate(this%k1)
            allocate (this%k1(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
            if (allocated(this%k2)) deallocate(this%k2)
            allocate (this%k2(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
            if (allocated(this%k3)) deallocate(this%k3)
            allocate (this%k3(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
            if (allocated(this%kabs_sq)) deallocate(this%kabs_sq)
            allocate (this%kabs_sq(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
            if (allocated(this%Gdealias)) deallocate(this%Gdealias)
            allocate (this%Gdealias(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
            if (allocated(this%GTestFilt)) deallocate(this%GTestFilt)
            allocate (this%GTestFilt(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
            
            ! STEP 3: Generate 1d wavenumbers 
            k1_1d = GetWaveNums(nx_g,dx) 
            k2_1d = GetWaveNums(ny_g,dy) 
            k3_1d = GetWaveNums(nz_g,dz) 
            ! flip the sign of the oddball (helps later)
            k1_1d(nx_g/2+1) = -k1_1d(nx_g/2+1)
            k2_1d(ny_g/2+1) = -k2_1d(ny_g/2+1)
            k3_1d(nz_g/2+1) = -k3_1d(nz_g/2+1)


       
            ! STEP 4: Create temporary array for k1 and transpose it to the appropriate dimension
            allocate(tmp1(spectdecomp%xsz(1),spectdecomp%xsz(2),spectdecomp%xsz(3)))
            do k = 1,size(tmp1,3)
                do j = 1,size(tmp1,2)
                    tmp1(1:spectdecomp%xsz(1),j,k) = k1_1d(1:spectdecomp%xsz(1))
                end do 
            end do 
            
            select case (rPencil)
            case(1)
                allocate(tmp2(spectdecomp%ysz(1),spectdecomp%ysz(2),spectdecomp%ysz(3)))
                if (TwoPeriodic) then
                    call transpose_x_to_y(tmp1,this%k1,spectdecomp)
                else
                    call transpose_x_to_y(tmp1,tmp2,spectdecomp)
                    call transpose_y_to_z(tmp2,this%k1,spectdecomp)
                end if 
                deallocate(tmp2)
            case(3)
                this%k1 = tmp1
            end select     
            deallocate(tmp1)

            ! STEP 5: Create temporary array for k2 and transpose it to the appropriate dimension
            allocate(tmp1(spectdecomp%ysz(1),spectdecomp%ysz(2),spectdecomp%ysz(3))) 
            do k = 1,size(tmp1,3)
                do i = 1,size(tmp1,1)
                    tmp1(i,1:spectdecomp%ysz(2),k) = k2_1d(1:spectdecomp%ysz(2))
                end do 
            end do 

            select case (rPencil)
            case (1)
                if (TwoPeriodic) then
                    this%k2 = tmp1
                else
                    call transpose_y_to_z(tmp1,this%k2,spectdecomp)
                end if 
            case(3)
                call transpose_y_to_x(tmp1,this%k2,spectdecomp)
            end select 
            deallocate(tmp1)

            ! STEP 6: Create temporary array for k3 and transpose it to the appropriate dimension
            allocate(tmp1(spectdecomp%zsz(1),spectdecomp%zsz(2),spectdecomp%zsz(3)))
            do j = 1,size(tmp1,2)
                do i = 1,size(tmp1,1)
                    tmp1(i,j,1:spectdecomp%zsz(3)) = k3_1d(1:spectdecomp%zsz(3))
                end do 
            end do 

            select case (rPencil)
            case (1)
                if (TwoPeriodic) then
                    allocate(tmp2(spectdecomp%ysz(1),spectdecomp%ysz(2),spectdecomp%ysz(3)))
                    call transpose_z_to_y(tmp1,tmp2,spectdecomp)
                    this%k3 = tmp2 
                    deallocate(tmp2)
                else
                    this%k3 = tmp1
                end if 
            case (3)
                allocate(tmp2(spectdecomp%ysz(1),spectdecomp%ysz(2),spectdecomp%ysz(3)))
                call transpose_z_to_y(tmp1,tmp2,spectdecomp)
                call transpose_y_to_x(tmp2,this%k3,spectdecomp)
                deallocate(tmp2)
            end select
            deallocate(tmp1)
            deallocate (spectdecomp)
            

            if (TwoPeriodic) then
                    do k = 1,size(this%k1,3)
                        do j = 1,size(this%k1,2)
                            do i = 1,size(this%k1,1)
                                this%kabs_sq(i,j,k) = this%k1(i,j,k)**2 + this%k2(i,j,k)**2
                            end do 
                        end do 
                    end do                 

            else 
                    this%kabs_sq = this%k1**2 + this%k2**2 + this%k3**2 
            end if 
            

            ! STEP 7: Create the dealiasing filter transfer function 
            !select case (filt)
            !case ("2/3rd")
            !    if (TwoPeriodic) then
            !        this%Gdealias = TwoThirdsRule2D(dx,dy,this%k1,this%k2)
           !     else
           !         this%Gdealias = TwoThirdsRule(nx_g,ny_g,nz_g,this%kabs_sq)
           !     end if 
           ! case ("cf90")
           !     allocate(tmp1(size(this%Gdealias,1),size(this%Gdealias,2),size(this%Gdealias,3)))
           !     tmp1 = GetCF90TransferFunction(this%k1,dx)
           !     this%Gdealias = tmp1
           !     tmp1 = GetCF90TransferFunction(this%k2,dy)
           !     this%Gdealias = this%Gdealias*tmp1
           !     tmp1 = GetCF90TransferFunction(this%k3,dz)
           !     if (TwoPeriodic) then
           !         this%Gdealias = this%Gdealias
           !     else
           !         this%Gdealias = this%Gdealias*tmp1
           !     end if 
           !     deallocate(tmp1)
           ! case ("none")
           !     this%Gdealias = one 
           ! case default
           !     call GracefulExit("The dealiasing filter specified is incorrect.",104)
           ! end select
            kdealiasx = ((two/three)*pi/dx)
            kdealiasy = ((two/three)*pi/dy)
            do k = 1,size(this%k1,3)
                do j = 1,size(this%k1,2)
                    do i = 1,size(this%k1,1)
                        if ((abs(this%k1(i,j,k)) < kdealiasx) .and. (abs(this%k2(i,j,k))< kdealiasy)) then
                            this%Gdealias(i,j,k) = one
                        else
                            this%Gdealias(i,j,k) = zero
                        end if
                    end do 
                end do  
            end do 
            call message(1, "Dealiasing Summary:")
            call message(2, "Total non zero:", p_sum(sum(this%Gdealias)))


            kdealiasx = kdealiasx/3.d0
            kdealiasy = kdealiasy/3.d0
            do k = 1,size(this%k1,3)
                do j = 1,size(this%k1,2)
                    do i = 1,size(this%k1,1)
                        if ((abs(this%k1(i,j,k)) < kdealiasx) .and. (abs(this%k2(i,j,k))< kdealiasy)) then
                            this%GTestFilt(i,j,k) = one
                        else
                            this%GTestFilt(i,j,k) = zero
                        end if
                    end do 
                end do  
            end do 
            call message(1, "TestFilter Summary:")
            call message(2, "Total non zero:", p_sum(sum(this%GTestFilt)))


        end if    

        ! STEP 10: Determine the sizes of the fft and ifft input arrays
        this%nx_c = this%fft_size(1); this%ny_c = this%fft_size(2); this%nz_c = this%fft_size(3)
        if (allocated(spectdecomp)) deallocate(spectdecomp)
        allocate(spectdecomp)
        call decomp_info_init(nx_g, ny_g, nz_g, spectdecomp) 
        select case (rPencil)
        case (1)
            this%nx_r = spectdecomp%xsz(1)
            this%ny_r = spectdecomp%xsz(2)
            this%nz_r = spectdecomp%xsz(3)
        case (3)
            this%nx_r = spectdecomp%zsz(1)
            this%ny_r = spectdecomp%zsz(2)
            this%nz_r = spectdecomp%zsz(3)
        end select
        deallocate(spectdecomp)

        ! STEP 11: Fix Oddball wavenumbers
        if (this%StoreK) then
            if (this%fixOddball) then 
                allocate(spectdecomp)
                select case (rPencil) 
                case (1) 
                    call decomp_info_init(nx_g/2+1, ny_g, nz_g, spectdecomp)
                    
                    if (TwoPeriodic) then
                        allocate (tmp1(spectdecomp%zsz(1),spectdecomp%zsz(2),spectdecomp%zsz(3)))
                        ! Fuck k3
                        this%k2(:,ny_g/2+1,:) = zero
                        allocate (tmp2(spectdecomp%xsz(1),spectdecomp%xsz(2),spectdecomp%xsz(3)))
                        call transpose_y_to_x(this%k1,tmp2,spectdecomp)
                        tmp2(nx_g/2+1,:,:) = zero
                        call transpose_x_to_y(tmp2,this%k1,spectdecomp)
                    else
                        this%k3(:,:,nz_g/2+1) = zero
                        allocate (tmp1(spectdecomp%ysz(1),spectdecomp%ysz(2),spectdecomp%ysz(3)))
                        allocate (tmp2(spectdecomp%xsz(1),spectdecomp%xsz(2),spectdecomp%xsz(3)))
                        
                        call transpose_z_to_y(this%k2,tmp1,spectdecomp)
                        tmp1(:,ny_g/2+1,:) = zero
                        call transpose_y_to_z(tmp1,this%k2,spectdecomp)
                        
                        call transpose_z_to_y(this%k1,tmp1,spectdecomp)
                        call transpose_y_to_x(tmp1,tmp2,spectdecomp)
                        tmp2(nx_g/2+1,:,:) = zero
                        call transpose_x_to_y(tmp2,tmp1,spectdecomp)
                        call transpose_y_to_z(tmp1,this%k1,spectdecomp)
                    end if 
                case (3)
                    call decomp_info_init(nx_g, ny_g, nz_g/2+1, spectdecomp) 
                    this%k1(nx_g/2+1,:,:) = zero
                    allocate (tmp1(spectdecomp%ysz(1),spectdecomp%ysz(2),spectdecomp%ysz(3)))
                    allocate (tmp2(spectdecomp%zsz(1),spectdecomp%zsz(2),spectdecomp%zsz(3)))
                    
                    call transpose_x_to_y(this%k2,tmp1,spectdecomp)
                    tmp1(:,ny_g/2+1,:) = zero
                    call transpose_y_to_x(tmp1,this%k2,spectdecomp)

                    call transpose_x_to_y(this%k3,tmp1,spectdecomp)
                    call transpose_y_to_z(tmp1,tmp2,spectdecomp)
                    tmp2(:,:,nz_g/2+1) = zero
                    call transpose_z_to_y(tmp2,tmp1,spectdecomp)
                    call transpose_y_to_x(tmp1,this%k3,spectdecomp)
                end select
                deallocate (tmp1, tmp2)
                deallocate(spectdecomp)
            end if
        end if 
        
        ! STEP 12: Carrying zero wavenumber

        if (TwoPeriodic) then
            do j = 1,size(this%k1,2)
                do i = 1,size(this%k1,1)
                    if ((abs(this%k1(i,j,1))<1.D-13) .and. (abs(this%k2(i,j,1))<1.D-13)) then
                        this%carryingZeroK = .true.
                        this%ZeroK_i = i
                        this%ZeroK_j = j
                        if ((i .ne. 1) .and. (j .ne. 1)) then
                            print*, nrank, i, j
                            call GracefulExit("Catastrophic failure while initializing spectral &
                                    derived type. Unable to isolate k1 = 0 and k2 = 0 wavenumbers.",312)
                        end if
                        !print*,  "Identified ZERO wavenumber on process:", nrank
                        !print*,  "i - index:", i
                        !print*,  "j - index:", j
                    end if 
                end do 
            end do 
        end if 
       

        ! STEP 13: Surface Filter
        allocate(tmp1(this%spectdecomp%zsz(1), this%spectdecomp%zsz(2), this%spectdecomp%zsz(3)))
        allocate(tmp2(this%spectdecomp%zsz(1), this%spectdecomp%zsz(2), this%spectdecomp%zsz(3)))
        call transpose_y_to_z(this%k1,tmp1,this%spectdecomp)
        call transpose_y_to_z(this%k2,tmp2,this%spectdecomp)
        allocate(this%Gsurfacefilter(this%spectdecomp%zsz(1), this%spectdecomp%zsz(2)))
        kdealiasx = (one/three)*pi/dx
        kdealiasy = (one/three)*pi/dy
        do j = 1,size(tmp1,2)
            do i = 1,size(tmp1,1)
                if ((abs(tmp1(i,j,1)) < kdealiasx) .and. (abs(tmp2(i,j,1))< kdealiasy)) then
                    this%GSurfaceFilter(i,j) = one
                else
                    this%GSurfaceFilter(i,j) = zero
                end if
            end do 
        end do 
     
        ! STEP 14: Prep filter for KS  
        allocate(this%GksPrep1(this%spectdecomp%ysz(1), this%spectdecomp%ysz(2), this%spectdecomp%ysz(3)))
        kdealiasx = (one/four)*pi/dx
        kdealiasy = (one/four)*pi/dy
        do k = 1,this%spectdecomp%ysz(3)
            do j = 1,this%spectdecomp%ysz(2)
                do i = 1,this%spectdecomp%ysz(1)
                    if ((abs(this%k1(i,j,k)) < kdealiasx) .and. (abs(this%k2(i,j,k))< kdealiasy)) then
                        this%GksPrep1(i,j,k) = one
                    else
                        this%GksPrep1(i,j,k) = zero
                    end if
                end do 
            end do 
        end do 
        allocate(this%GksPrep2(this%spectdecomp%ysz(1), this%spectdecomp%ysz(2), this%spectdecomp%ysz(3)))
        kdealiasx = (one/eight)*pi/dx
        kdealiasy = (one/eight)*pi/dy
        do k = 1,this%spectdecomp%ysz(3)
            do j = 1,this%spectdecomp%ysz(2)
                do i = 1,this%spectdecomp%ysz(1)
                    if ((abs(this%k1(i,j,k)) < kdealiasx) .and. (abs(this%k2(i,j,k))< kdealiasy)) then
                        this%GksPrep2(i,j,k) = one
                    else
                        this%GksPrep2(i,j,k) = zero
                    end if
                end do 
            end do 
        end do 
        ! STEP 15: Allocate 3/2s rule arrays
        !call this%FT%alloc_upsampledArr(this%arr1Up)
        !call this%FT%alloc_upsampledArr(this%arr2Up)


        ! STEP 16: Print completion message 
        call message("SPECTRAL - Derived Type for the problem generated successfully.")
        call message("===============================================================")
        ! Finished !
    end subroutine
    
    subroutine destroy(this)
        class(spectral), intent(inout) :: this
      
        if (.not. this%isInitialized) then
            call GracefulExit("You are trying to destroy a SPECTRAL derived type before initializing it",110)
        end if 
        if (this%use2decompFFT) then
            call decomp_2d_fft_finalize 
            this%use2decompFFT = .false. 
        else
            call this%FT%destroy
            deallocate(this%FT)
        end if 
        if (allocated(this%Gdealias)) then
                deallocate(this%Gdealias)
        end if
        if (allocated(this%GsurfaceFilter)) deallocate(this%GsurfaceFilter)
        if (allocated(this%k1)) deallocate(this%k1) 
        if (allocated(this%k2)) deallocate(this%k2) 
        if (allocated(this%k3)) deallocate(this%k3) 
        if (allocated(this%kabs_sq)) deallocate(this%kabs_sq) 
        if (allocated(this%k1_der2)) deallocate(this%k1_der2)
        if (allocated(this%k2_der2)) deallocate(this%k2_der2)
        if (allocated(this%k3_der2)) deallocate(this%k3_der2)
        if (allocated(this%one_by_kabs_sq)) deallocate(this%one_by_kabs_sq)
      
        if (allocated(this%arr1Up)) deallocate(this%arr1Up)
        if (allocated(this%arr2Up)) deallocate(this%arr2Up)
        if (allocated(this%spectdecomp)) deallocate(this%spectdecomp)
        this%isInitialized = .false. 
   

    end subroutine 

    subroutine fft1_x2y(this,arr_in,arr_out)  !this%u,this%cbuffC(:,:,:,1))
        !use decomp_2d_fft, only: decomp_2d_fft_3d
        class(spectral), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in) :: arr_in
        complex(rkind), dimension(:,:,:), intent(out) :: arr_out

        call this%FT%fft1_x2y(arr_in,arr_out)

    end subroutine 

    subroutine fft(this,arr_in,arr_out)
        use decomp_2d_fft, only: decomp_2d_fft_3d
        class(spectral), intent(inout) :: this
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(in) :: arr_in
        complex(rkind), dimension(this%nx_c,this%ny_c,this%nz_c), intent(out) :: arr_out

        if (this%is3dFFT) then
            select case (this%rPencil) 
            case (1)
                call this%FT%fft3_x2z(arr_in,arr_out)
            case (3)
                call this%FT%fft3_z2x(arr_in,arr_out)
            end select
        else
            call this%FT%fft2_x2y(arr_in,arr_out)        
        end if 

    end subroutine

    subroutine ifft(this,arr_in,arr_out,setOddball)
        use decomp_2d_fft, only: decomp_2d_fft_3d
        class(spectral), intent(inout) :: this
        complex(rkind), dimension(this%nx_c,this%ny_c,this%nz_c), intent(in) :: arr_in
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(out) :: arr_out 
        logical, intent(in), optional :: setOddball

        if (this%is3dFFT) then
            select case (this%rPencil)
            case (1)
                call this%FT%ifft3_z2x(arr_in,arr_out)
            case (3)
                call this%FT%ifft3_x2z(arr_in,arr_out)
            end select
        else
            if (present(setOddBall)) then
                    call this%FT%ifft2_y2x(arr_in,arr_out,setOddball) 
            else                    
               call this%FT%ifft2_y2x(arr_in,arr_out,.false.)
            end if
        end if 
    end subroutine
   
    subroutine fft_y2z(this, arr_in, arr_out)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%ysz(1), this%spectdecomp%ysz(2), this%spectdecomp%ysz(3)), intent(in)  :: arr_in
      complex(rkind), dimension(this%spectdecomp%zsz(1), this%spectdecomp%zsz(2), this%spectdecomp%zsz(3)), intent(out) :: arr_out

      call transpose_y_to_z(arr_in, arr_out, this%spectdecomp)
      call dfftw_execute_dft(this%plan_c2c_fwd_z_ip, arr_out , arr_out)  

    end subroutine 


    subroutine ifft_z2y(this, arr_in, arr_out) 
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%ysz(1), this%spectdecomp%ysz(2), this%spectdecomp%ysz(3)), intent(in)  :: arr_in
      complex(rkind), dimension(this%spectdecomp%zsz(1), this%spectdecomp%zsz(2), this%spectdecomp%zsz(3)), intent(out) :: arr_out

      call dfftw_execute_dft(this%plan_c2c_bwd_z_oop, arr_in , this%ctmpz)  
      this%ctmpz = this%normfactz*this%ctmpz
      call transpose_z_to_y(this%ctmpz, arr_out, this%spectdecomp)

    end subroutine 

    subroutine take_fftz(this, arr_in)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%ysz(1), this%spectdecomp%ysz(2), this%spectdecomp%ysz(3)), intent(in)  :: arr_in

      call transpose_y_to_z(arr_in, this%ctmpz, this%spectdecomp)
      call dfftw_execute_dft(this%plan_c2c_fwd_z_ip, this%ctmpz , this%ctmpz)  

    end subroutine 
  

    subroutine take_fft1d_z2z_ip(this, arr_inout)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%zsz(1), this%spectdecomp%zsz(2), this%spectdecomp%zsz(3)), intent(inout)  :: arr_inout
      call dfftw_execute_dft(this%plan_c2c_fwd_z_ip, arr_inout , arr_inout)  
    end subroutine 

    subroutine take_ifft1d_z2z_ip(this, arr_inout)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%zsz(1), this%spectdecomp%zsz(2), this%spectdecomp%zsz(3)), intent(inout)  :: arr_inout
      call dfftw_execute_dft(this%plan_c2c_bwd_z_ip, arr_inout , arr_inout)  
      arr_inout = this%normfactz*arr_inout
    end subroutine 


    subroutine take_ifftz(this, arr_out)
      class(spectral), intent(inout) :: this
      complex(rkind), dimension(this%spectdecomp%ysz(1), this%spectdecomp%ysz(2), this%spectdecomp%ysz(3)), intent(out)  :: arr_out

      call dfftw_execute_dft(this%plan_c2c_bwd_z_ip, this%ctmpz , this%ctmpz)  
      this%ctmpz = this%normfactz*this%ctmpz
      call transpose_z_to_y(this%ctmpz, arr_out, this%spectdecomp)

    end subroutine 


    !!!!!!!!!!!!!!!!!!!!!! DEALISED MULTIPLICATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dealiasedMult_oop(this,arr1,arr2,arrout)
        class(spectral), intent(inout) :: this
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(in) :: arr1, arr2
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(out) :: arrOut
       
        if (use3by2rule) then 
            call this%FT%upsample(arr1,this%arr1Up)
            call this%FT%upsample(arr2,this%arr2Up)
            this%arr1Up = this%arr1Up * this%arr2Up
            call this%FT%downsample(this%arr1Up,arrout)
        else
            arrout = arr1*arr2
        end if 


    end subroutine

    subroutine dealiasedMult_ip(this,arr1,arr2)
        class(spectral), intent(inout) :: this
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(in) :: arr2
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(inout) :: arr1
       
        if (use3by2rule) then 
            call this%FT%upsample(arr1,this%arr1Up)
            call this%FT%upsample(arr2,this%arr2Up)
            this%arr1Up = this%arr1Up * this%arr2Up
            call this%FT%downsample(this%arr1Up,arr1)
        else
            arr1 = arr1*arr2
        end if 


    end subroutine

    subroutine dealiasedDiv_oop(this,arrNum,arrDen,arrout)
        class(spectral), intent(inout) :: this
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(in) :: arrNum, arrDen
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(out) :: arrOut
       
        if (use3by2rule) then 
            call this%FT%upsample(arrNum,this%arr1Up)
            call this%FT%upsample(arrDen,this%arr2Up)
            this%arr1Up = this%arr1Up / ( this%arr2Up + 1d-13)
            call this%FT%downsample(this%arr1Up,arrout)
        else
            arrout = arrNum/(arrDen + 1d-13)
        end if 

    end subroutine
    
    subroutine dealiasedDiv_ip(this,arrNum,arrDen)
        class(spectral), intent(inout) :: this
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(in) :: arrDen
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(inout) :: arrNum
       
        if (use3by2rule) then 
            call this%FT%upsample(arrNum,this%arr1Up)
            call this%FT%upsample(arrDen,this%arr2Up)
            this%arr1Up = this%arr1Up / ( this%arr2Up + 1d-13)
            call this%FT%downsample(this%arr1Up,arrNum)
        else
            arrNum = arrNum/(arrDen + 1d-13)
        end if 

    end subroutine
    
    subroutine dealiasedSquare_oop(this,arr1,arrout)
        class(spectral), intent(inout) :: this
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(in) :: arr1
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(out) :: arrOut
      
        if (use3by2rule) then 
            call this%FT%upsample(arr1,this%arr1Up)
            this%arr1Up = this%arr1Up * this%arr1Up
            call this%FT%downsample(this%arr1Up,arrout)
        else
            arrout = arr1*arr1
        end if 

    end subroutine
    
    subroutine dealiasedSquare_ip(this,arr1)
        class(spectral), intent(inout) :: this
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r), intent(inout) :: arr1
      
        if (use3by2rule) then 
            call this%FT%upsample(arr1,this%arr1Up)
            this%arr1Up = this%arr1Up * this%arr1Up
            call this%FT%downsample(this%arr1Up,arr1)
        else
            arr1 = arr1*arr1
        end if 

    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine alloc_r2c_out_Rank3(this,arr_inout)
        class(spectral), intent(in) :: this
        complex(rkind), dimension(:,:,:), allocatable, intent(out) :: arr_inout
        
        if (.not. this%isInitialized) then
            call GracefulExit("You cannot use alloc_r2c_out before initializing the SPECTRAL derived type",112)
        end if 
        allocate(arr_inout(this%fft_size(1),this%fft_size(2),this%fft_size(3)))
    end subroutine
    
    subroutine alloc_r2c_out_Rank4(this,arr_inout,numvars)
        class(spectral), intent(in) :: this
        integer, intent(in) :: numvars
        complex(rkind), dimension(:,:,:,:), allocatable, intent(out) :: arr_inout
        
        if (.not. this%isInitialized) then
            call GracefulExit("You cannot use alloc_r2c_out before initializing the SPECTRAL derived type",112)
        end if 
        allocate(arr_inout(this%fft_size(1),this%fft_size(2),this%fft_size(3),numvars))
    end subroutine
    
    function TwoThirdsRule(Nx,Ny,Nz,kin) result(Tf)
        use constants, only: one, zero, three
        real(rkind), intent(in), dimension(:,:,:) :: kin
        integer, intent(in) :: Nx, Ny, Nz
        real(rkind), dimension(size(kin,1),size(kin,2),size(kin,3)) :: Tf

        real(rkind) :: kd1, kd2, kd3, kdealias
        
        kd1 = real(floor(real(Nx)/three))
        kd2 = real(floor(real(Ny)/three))
        kd3 = real(floor(real(Nz)/three))
        kdealias = min(kd1,kd2,kd3)
        where (kin .ge. kdealias**2) 
            Tf = zero
        elsewhere
            Tf = one
        end where

    end function

    function TwoThirdsRule2D(dx,dy,k1,k2) result(Tf)
        use constants, only: one, pi, zero, two, three
        real(rkind), intent(in) :: dx, dy 
        real(rkind), dimension(:,:,:), intent(in) :: k1, k2
        real(rkind), dimension(size(k1,1),size(k1,2), size(k1,3)) :: Tf
        real(rkind) :: kxd, kyd, kdealias

        kxd = (two/three)*pi/dx
        kyd = (two/three)*pi/dy
        kdealias = min(kxd,kyd)
        Tf = one
        where(abs(k1)>kdealias)
                Tf = zero
        end where
        
        where(abs(k2)>kdealias)
                Tf = zero
        end where


    end function
    

    pure elemental function GetCF90TransferFunction(kin,dx) result(T)
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
        kp = kp/dx 

    end function


   
    pure elemental function GetCD10D2ModWaveNum(kin,dx) result(wpp)
        use cd10stuff, only: alpha10d2, beta10d2, a10d2, b10d2, c10d2
        real(rkind), intent(in) :: kin,dx
        real(rkind) :: w
        real(rkind) :: wpp
        real(rkind) :: a, b, c, alpha, beta

        alpha = alpha10d2
        beta = beta10d2
        a = a10d2*1._rkind
        b = b10d2*4._rkind
        c = c10d2*9._rkind

        w = kin*dx
        wpp = (2._rkind*a*(1._rkind-cos(w)) + (b/2._rkind)*(1._rkind - cos(2._rkind*w)) &
        + (2*c/9)*(1._rkind - cos(3._rkind*w)))/(1._rkind + 2._rkind*alpha*cos(w) + 2*beta*cos(2*w))
  
        wpp = wpp/dx

    end function


    pure elemental function GetCD06ModWaveNum(kin,dx) result(kp)
        use constants
        real(rkind), intent(in) :: kin,dx
        real(rkind) :: k
        real(rkind) :: kp
        real(rkind) :: alpha, beta, a, b, c, d
        real(rkind) :: num, den 

        alpha = one/three
        beta = zero
        a = (two/nine)*(eight - three*alpha)
        b = (one/18._rkind)*(-17._rkind + 57._rkind*alpha)
        c = zero
        d = zero

        k = kin*dx

        den = 1 + 2*alpha*cos(k) + 2*beta*cos(2*k)
        num = a*sin(k) + (b/two)*sin(two*k) + (c/three)*sin(three*k) + (d/four)*sin(four*k)
        kp = num/den 
        
        kp = kp/dx 
    
    end function
    
    pure function GetWaveNums(nx,dx) result(k)
        use constants, only: pi, two
        integer, intent(in) :: nx
        real(rkind), intent(in) :: dx
        real(rkind), dimension(nx) :: k

        integer :: i,dummy

        dummy = nx - MOD(nx,2)

        do i = 1,nx
            k(i) = ( -pi + (i-1)*two*pi/real(dummy,rkind) ) / dx
        end do

        k = ifftshift(k)

    end function
    
    pure function ifftshift(k) result(kshift)

        real(rkind), dimension(:), intent(in) :: k
        real(rkind), dimension(SIZE(k)) :: kshift
        integer :: n

        n = SIZE(k)

        select case ( MOD(n,2) )
        case (0)
            kshift(1:n/2) = k(n/2+1:n)
            kshift(n/2+1:n) = k(1:n/2)
        case (1)
            kshift(1:(n+1)/2) = k((n+1)/2:n)
            kshift((n+1)/2+1:n) = k(1:(n-1)/2)
        end select

    end function
   
    subroutine initPP(this,ffactx, ffacty, dx, dy)
      use reductions, only: p_sum
      use constants, only: pi
      class(spectral), intent(inout) :: this
      integer :: i, j, k
      real(rkind), intent(in) :: ffactx, ffacty, dx, dy
      real(rkind) :: kfilx, kfily
      
      if (this%initPostProcessor) call this%destroyPP()
      
      kfilx = ffactx*pi/dx
      kfily = ffacty*pi/dy

      allocate(this%G_PostProcess(this%spectdecomp%ysz(1),this%spectdecomp%ysz(2),this%spectdecomp%ysz(3)))
      do k = 1,this%spectdecomp%ysz(3)
          do j = 1,this%spectdecomp%ysz(2)
              do i = 1,this%spectdecomp%ysz(1)
                  if ((abs(this%k1(i,j,k)) < kfilx) .and. (abs(this%k2(i,j,k))< kfily)) then
                      this%G_PostProcess(i,j,k) = 1
                  else
                      this%G_PostProcess(i,j,k) = 0
                  end if
              end do 
          end do 
      end do 

      call message(0,"Post-processing filter initialized")
      call message(1,"Size of G:",p_sum(size(this%G_PostProcess)))
      call message(2,"Number of non-zeros:",p_sum(sum(this%G_PostProcess)))
      this%initPostProcessor = .true.
    end subroutine

    subroutine destroyPP(this)
      class(spectral), intent(inout) :: this

      this%initPostProcessor = .false.
      if (allocated(this%G_PostProcess)) deallocate(this%G_PostProcess)

    end subroutine

    subroutine spectralfilter_ip(this,arr)
      class(spectral), intent(in) :: this
      complex(rkind), dimension(this%spectdecomp%ysz(1),this%spectdecomp%ysz(2),this%spectdecomp%ysz(3)), intent(inout) :: arr

      arr = this%G_PostProcess*arr
    end subroutine
end module 
