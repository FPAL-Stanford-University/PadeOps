module fourierGaboMod
   use kind_parameters, only: rkind
   use nufftMod, only: nufft
   use constants, only: imi, one, two, four, three, six, zero, pi, half 
   implicit none
   private
   public :: fourierGabor


   type :: fourierGabor
      private
      real(rkind),    dimension(:,:), allocatable :: wavenums
      complex(rkind), dimension(:,:), allocatable :: velhats
      real(rkind)  :: xcenter, ycenter, zcenter
      real(rkind), dimension(:,:,:), allocatable :: window, f
      integer :: nx, ny, nz
 
      real(rkind), dimension(:), allocatable :: rbuff1, rbuff2, rbuff3
      complex(rkind), dimension(:), allocatable :: cbuff1
      complex(rkind), dimension(:,:,:), allocatable :: cbuff2

      real(rkind) :: dxFine, dyFine, dzFine
      real(rkind) :: xshift, yshift, zshift 
      type(nufft) :: my_nufft
      
      integer :: nmodes
      logical :: isActive = .false. 

      contains
         procedure :: init 
         procedure :: destroy
         procedure :: getField
   end type

contains 

subroutine init(this, kmin, kmax, nk, ntheta, xcenter, ycenter, zcenter, xwidth, ywidth, zwidth, xline, yline, zline, seed, &
                     eps_nufft)
   class(fourierGabor), intent(inout), target :: this
   real(rkind), intent(in) :: kmin, kmax, xcenter, ycenter, zcenter, xwidth, ywidth, zwidth
   integer, intent(in) :: nk, ntheta, seed
   real(rkind), dimension(:), allocatable :: xline, yline, zline
   real(rkind), intent(in) :: eps_nufft
   real(rkind), dimension(:), pointer :: kx, ky, kz
   complex(rkind), dimension(:), pointer :: uhat, vhat, what 

   this%nx = size(xline); this%ny = size(yline); this%nz = size(zline)  
   allocate(this%window(this%nx, this%ny, this%nz))
   call get_window(xline,yline,zline,[xcenter,ycenter,zcenter],[two*xwidth,two*ywidth,two*zwidth], this%window, this%isActive)
   
   if (this%isActive) then
      allocate(this%wavenums(nk*ntheta,3), this%velhats(nk*ntheta,3))
      allocate(this%cbuff1(nk*ntheta))
      allocate(this%cbuff2(this%nx, this%ny, this%nz))
      uhat => this%velhats(:,1); vhat => this%velhats(:,2); what => this%velhats(:,3)
      kx  => this%wavenums(:,1); ky => this%wavenums(:,2);  kz => this%wavenums(:,3)
      
      this%nmodes = nk*ntheta
      call getIsotropicModes(nk,ntheta,kmin,kmax,seed,kx,ky,kz,uhat,vhat,what)

      this%xshift = half*(maxval(xline) + minval(xline))
      this%yshift = half*(maxval(xline) + minval(xline))
      this%zshift = half*(maxval(xline) + minval(xline))

      this%dxFine = xline(2) - xline(1) 
      this%dyFine = yline(2) - yline(1)
      this%dzFine = zline(2) - zline(1)

      call this%my_nufft%init(eps_nufft, this%nx, this%ny, this%nz)
      nullify(uhat, vhat, what, kx, ky, kz)
   else
      deallocate(this%window)
   end if 

end subroutine 

subroutine getField(this, field, fid)
   class(fourierGabor), intent(inout) :: this
   real(rkind), dimension(this%nx, this%ny, this%nz), intent(inout) :: field
   integer, intent(in) :: fid

   if (this%isActive) then
      this%cbuff1 = this%velhats(:,fid)*exp(imi*(this%wavenums(:,1)*this%xshift + &
                    this%wavenums(:,2)*this%yshift + this%wavenums(:,3)*this%zshift))
      
      this%rbuff1 = this%dxFine*this%wavenums(:,1)
      this%rbuff2 = this%dyFine*this%wavenums(:,2)
      this%rbuff3 = this%dzFine*this%wavenums(:,3)

      call this%my_nufft%fft(this%nmodes, this%rbuff1, this%rbuff2, this%rbuff3, this%cbuff1, this%cbuff2) 

      field = field + this%nmodes*this%window*real(this%cbuff2,rkind)
   end if 

end subroutine 

subroutine destroy(this)
   class(fourierGabor), intent(inout) :: this

   if (this%isActive) then
      deallocate(this%rbuff1, this%rbuff2, this%rbuff3)
      deallocate(this%cbuff1, this%cbuff2, this%window)
   end if
   this%isActive = .false. 

end subroutine 


subroutine get_window(xin,yin,zin,center,width, window, HasInfluence)
   implicit none
   real(rkind), intent(in), dimension(:) :: xin, yin, zin
   real(rkind), intent(in), dimension(3) :: center, width
   real(rkind), intent(out), dimension(:,:,:) :: window
   logical, intent(out) :: HasInfluence 

   real(rkind) :: xleft, yleft, zleft, xright, yright, zright
   integer :: i, j, k
   real(rkind), dimension(size(xin)) :: x
   real(rkind), dimension(size(yin)) :: y
   real(rkind), dimension(size(zin)) :: z
   integer :: iterator 

   xleft = center(1) - width(1)/two;
   xright = center(1) + width(1)/two;
   yleft = center(2) - width(2)/two
   yright = center(2) + width(2)/two
   zleft = center(3) - width(3)/two
   zright = center(3) + width(3)/two
   
   x = pi*(xin - xleft)/(xright - xleft);
   y = pi*(yin - yleft)/(yright - yleft)
   z = pi*(zin - zleft)/(zright - zleft)
   
   window = zero
   iterator = 0

   do k = 1,size(window,3)
      if ((z(k) < pi) .and. (z(k) > zero)) then
         do j = 1,size(window,2)
            if ((y(j) < pi) .and. (y(j) > zero)) then
               do i = 1,size(window,1)
                  if ((x(i) < pi) .and. (x(i) > zero)) then
                     window(i,j,k) = sin(x(i))*sin(y(j))*sin(z(k))
                     iterator = iterator + 1
                  end if 
               end do 
            end if 
         end do 
      end if 
   end do
  
   if (iterator == 0) then
      HasInfluence = .false. 
   else
      HasInfluence = .true. 
   end if
end subroutine 

subroutine getIsotropicModes(nk,ntheta,kmin,kmax, seed, kx, ky, kz, uhat, vhat, what)
    use gridtools, only: logspace
    use random, only: uniform_random, gaussian_random  
    implicit none
    integer, intent(in) :: nk, ntheta
    real(rkind), intent(in) :: kmin, kmax
    integer, intent(in) :: seed
    integer :: seed1, seed2, seed3, seed4, seed5 

    real(rkind), dimension(:,:,:), allocatable, target :: randomNums
    real(rkind), dimension(:), allocatable, target :: k_abs, Ek
    real(rkind), dimension(:), allocatable :: k1, k2, k3, amag, bmag
    real(rkind), dimension(:,:), allocatable:: a, b, ua, ub
    integer :: idx, str_idx, end_idx

    real(rkind) :: dk, uabs
    real(rkind), pointer :: mag, zeta(:), theta(:), alpha(:) 
    complex(rkind), dimension(nk*ntheta), intent(out) :: uhat, vhat, what
    real(rkind), dimension(nk*ntheta), intent(out) :: kx, ky, kz

    real(rkind), parameter :: C_vk = 1.4528_rkind

    seed1 = 2134*seed; seed2 = seed1 + 353123; seed3 = seed1 + 42341
    seed4 = seed1 + 12344; seed5 = seed1 + 915412; 

    ! Step 1: Generate shell radii
    allocate(k_abs(nk))
    k_abs = logspace(real(log10(kmin),rkind),real(log10(kmax),rkind),nk)
   
    ! Step 2: Allocate all the other local scope variables
    allocate(Ek(nk))
    allocate(randomNums(ntheta, nk, 4)) 

    ! Step 3: Create the energy spectra for KE and PE
    Ek = C_vk*(k_abs**four)/((one + k_abs**two)**(17.0_rkind/six))

    ! Step 4: Generate the random numbers
    call uniform_random(randomNums(:,:,1),-one,one,seed1)  
    call uniform_random(randomNums(:,:,2:4),zero,two*pi,seed2)  

    ! Step 5: Allocate local quantities for each shell
    allocate(k1(ntheta),k2(ntheta),k3(ntheta))
    allocate(a(ntheta,3),b(ntheta,3), amag(ntheta),bmag(ntheta))
    allocate(ua(ntheta,3),ub(ntheta,3))

    
    ! Step 6: Create the shells
    do idx = 1,nk
        mag => k_abs(idx)
        zeta => randomNums(:,idx,1)
        theta => randomNums(:,idx,2)
        k1 = mag*sqrt(1 - zeta**2)*cos(theta)
        k2 = mag*sqrt(1 - zeta**2)*sin(theta)
        k3 = mag*zeta

        if (idx .eq. 1) then
            dk = (k_abs(2) - k_abs(1))/two
        elseif (idx .eq. nk) then
            dk = (k_abs(nk) - k_abs(nk-1))/two
        else
            dk = (k_abs(idx+1) - k_abs(idx))/two + (k_abs(idx) - k_abs(idx-1))/two
        end if

        uabs = sqrt(two*Ek(idx)*dk/real(ntheta,rkind))
        a(:,1) =  zero
        a(:,2) = -k3
        a(:,3) =  k2
        b(:,1) =  k2**2 + k3**2
        b(:,2) = -k1*k2
        b(:,3) = -k1*k3
        amag = sqrt(a(:,1)**2 + a(:,2)**2 + a(:,3)**2)
        bmag = sqrt(b(:,1)**2 + b(:,2)**2 + b(:,3)**2)
        
        a(:,1) = a(:,1)/amag
        a(:,2) = a(:,2)/amag
        a(:,3) = a(:,3)/amag
        
        b(:,1) = b(:,1)/bmag
        b(:,2) = b(:,2)/bmag
        b(:,3) = b(:,3)/bmag

        alpha => randomNums(:,idx,3)
        ua(:,1) = uabs*(cos(alpha)*a(:,1) + sin(alpha)*b(:,1))
        ua(:,2) = uabs*(cos(alpha)*a(:,2) + sin(alpha)*b(:,2))
        ua(:,3) = uabs*(cos(alpha)*a(:,3) + sin(alpha)*b(:,3))

        alpha => randomNums(:,idx,4)
        ub(:,1) = uabs*(cos(alpha)*a(:,1) + sin(alpha)*b(:,1))
        ub(:,2) = uabs*(cos(alpha)*a(:,2) + sin(alpha)*b(:,2))
        ub(:,3) = uabs*(cos(alpha)*a(:,3) + sin(alpha)*b(:,3))

        str_idx = (idx - 1)*ntheta + 1
        end_idx = str_idx + ntheta - 1
        kx(str_idx:end_idx) = k1
        ky(str_idx:end_idx) = k2
        kz(str_idx:end_idx) = k3
        uhat(str_idx:end_idx) = cmplx(ua(:,1),ub(:,1),rkind) 
        vhat(str_idx:end_idx) = cmplx(ua(:,2),ub(:,2),rkind) 
        what(str_idx:end_idx) = cmplx(ua(:,3),ub(:,3),rkind) 

    end do  

    ! Penultimate Step: Deallocate all local scope variables
    deallocate(k_abs, Ek, ua, ub, amag, bmag, a, b)
    nullify(alpha, theta, zeta)
    deallocate(randomNums)
    
    !call uniform_random(xPos,zero,one,seed3) 
    !call uniform_random(yPos,zero,one,seed4) 
    !call uniform_random(zPos,zero,one,seed5) 

    !xPos = widths(1)*(xPos - half) + centers(1)
    !yPos = widths(2)*(yPos - half) + centers(2)
    !zPos = widths(3)*(zPos - half) + centers(3)

end subroutine





end module 
