module StratifiedShearLayer_IO

    use kind_parameters, only: rkind, clen
    use decomp_2d,       only: decomp_info,nrank,nproc
    use basic_io,        only: read_2d_ascii 
    use exits,           only: message
    implicit none

contains

subroutine read_Domain_info(Lx,Ly,Lz,fname)
   real(rkind), intent(out) :: Lx, Ly, Lz
   character(len=*), intent(in) :: fname
   
   real(rkind), dimension(:,:), allocatable :: data2read

   call read_2d_ascii(data2read,fname) 

   Lx = data2read(1,1)
   Ly = data2read(2,1)
   Lz = data2read(3,1)

   deallocate(data2read)
end subroutine

subroutine get_perturbations(gp, x, y, InitFileTag, InitFileDirectory, u, v, w, T)
   use constants, only: imi, pi   
   type(decomp_info), intent(in) :: gp
   character(len=*), intent(in) :: InitFileTag, InitFileDirectory
   real(rkind), dimension(:,:,:), intent(in) :: x, y
   real(rkind), dimension(:,:,:), intent(out) :: u, v, w, T

   real(rkind), dimension(:), allocatable :: kx, ky, kmode, beta
   real(rkind), dimension(:,:), allocatable :: data2read
   character(len=clen) :: fname 
   integer :: nmodes, nz
   real(rkind), dimension(:,:), allocatable :: uhat_real, uhat_imag, vhat_real, vhat_imag, what_real, what_imag, That_real, That_imag
   complex(rkind), dimension(:,:), allocatable :: uhat, vhat, what, That
   complex(rkind)  :: expfact
   
   integer :: modeID, j, k, ziter, i

   fname = InitFileDirectory(:len_trim(InitFileDirectory))//"/"//trim(InitFileTag)//"_mode_info.dat"
   call read_2d_ascii(data2read,fname)
   nmodes = size(data2read,1)
   allocate(kx(nmodes))
   allocate(kmode(nmodes))
   kmode = data2read(:,1)
   kx = kmode

   deallocate(kmode)
   deallocate(data2read)
   call message(0, "Number of normal modes being used:",nmodes)
   
   nz = gp%zsz(3) ! Global nz
   allocate(uhat_real(nz,nmodes),vhat_real(nz,nmodes),what_real(nz,nmodes),uhat_imag(nz,nmodes),vhat_imag(nz,nmodes),what_imag(nz,nmodes))
   allocate(That_real(nz,nmodes), That_imag(nz,nmodes))

   fname = InitFileDirectory(:len_trim(InitFileDirectory))//"/"//trim(InitFileTag)//"_init_info_imag.dat"
   call read_2d_ascii(data2read,fname)
   uhat_imag = reshape(data2read(:,1),[nz,nmodes])
   vhat_imag = reshape(data2read(:,2),[nz,nmodes])
   what_imag = reshape(data2read(:,3),[nz,nmodes])
   That_imag = reshape(data2read(:,4),[nz,nmodes])
   deallocate(data2read)

   fname = InitFileDirectory(:len_trim(InitFileDirectory))//"/"//trim(InitFileTag)//"_init_info_real.dat"
   call read_2d_ascii(data2read,fname)
   uhat_real = reshape(data2read(:,1),[nz,nmodes])
   vhat_real = reshape(data2read(:,2),[nz,nmodes])
   what_real = reshape(data2read(:,3),[nz,nmodes])
   That_real = reshape(data2read(:,4),[nz,nmodes])
   deallocate(data2read)
   
   allocate(uhat(nz,nmodes), vhat(nz,nmodes), what(nz,nmodes), That(nz,nmodes))

   uhat = uhat_real + imi*uhat_imag
   vhat = vhat_real + imi*vhat_imag
   what = what_real + imi*what_imag
   That = That_real + imi*That_imag
   deallocate(uhat_real,vhat_real,what_real,uhat_imag,vhat_imag,what_imag,That_real,That_imag)
 

   u = 0.d0
   v = 0.d0
   w = 0.d0
   T = 0.d0 

   do modeID = 1,nmodes
      ziter = 1
      do k = gp%xst(3),gp%xen(3)
         do j = 1,gp%xsz(2)
            do i = 1,gp%xsz(1)
               expfact = exp(imi*(kx(modeID)*x(i,1,1)))
               u(i,j,ziter) = u(i,j,ziter) + real(uhat(k, modeID)*expfact , rkind) 
               v(i,j,ziter) = v(i,j,ziter) + real(vhat(k, modeID)*expfact , rkind) 
               w(i,j,ziter) = w(i,j,ziter) + real(what(k, modeID)*expfact , rkind)
               T(i,j,ziter) = T(i,j,ziter) + real(That(k, modeID)*expfact , rkind)
            end do 
         end do 
         ziter = ziter + 1
      end do
   end do 

   deallocate(uhat, vhat, what, That)
   deallocate(kx)

end subroutine 

end module 
