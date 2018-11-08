module MixedConvection_IO 

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

subroutine get_perturbations(gp, x, y, InitFileTag, InitFileDirectory, u, v, w)
   use constants, only: imi, pi   
   type(decomp_info), intent(in) :: gp
   character(len=*), intent(in) :: InitFileTag, InitFileDirectory
   real(rkind), dimension(:,:,:), intent(in) :: x, y
   real(rkind), dimension(:,:,:), intent(out) :: u, v, w

   real(rkind), dimension(:), allocatable :: kx, ky, kmode, beta
   real(rkind), dimension(:,:), allocatable :: data2read
   character(len=clen) :: fname 
   integer :: nmodes, nz
   real(rkind), dimension(:,:), allocatable :: uhat_real, uhat_imag, vhat_real, vhat_imag, what_real, what_imag
   complex(rkind), dimension(:,:), allocatable :: uhat, vhat, what
   complex(rkind)  :: expfact, uhatn, vhatn, whatn
   
   integer :: modeID, j, k, ziter, i

   fname = InitFileDirectory(:len_trim(InitFileDirectory))//"/"//trim(InitFileTag)//"_mode_info.dat"
   call read_2d_ascii(data2read,fname)
   nmodes = size(data2read,1)
   allocate(kx(nmodes), ky(nmodes))
   allocate(beta(nmodes), kmode(nmodes))
   beta  = data2read(:,1)
   kmode = data2read(:,2)
   kx = -kmode*sin(beta*pi/180.d0)
   ky =  kmode*cos(beta*pi/180.d0)

   deallocate(kmode)
   deallocate(data2read)
   call message(0, "Number of normal modes being used:",nmodes)
   
   nz = gp%zsz(3) ! Global nz
   allocate(uhat_real(nz,nmodes),vhat_real(nz,nmodes),what_real(nz,nmodes),uhat_imag(nz,nmodes),vhat_imag(nz,nmodes),what_imag(nz,nmodes))


   fname = InitFileDirectory(:len_trim(InitFileDirectory))//"/"//trim(InitFileTag)//"_init_info_imag.dat"
   call read_2d_ascii(data2read,fname)
   uhat_imag = reshape(data2read(:,1),[nz,nmodes])
   vhat_imag = reshape(data2read(:,2),[nz,nmodes])
   what_imag = reshape(data2read(:,3),[nz,nmodes])
   deallocate(data2read)

   fname = InitFileDirectory(:len_trim(InitFileDirectory))//"/"//trim(InitFileTag)//"_init_info_real.dat"
   call read_2d_ascii(data2read,fname)
   uhat_real = reshape(data2read(:,1),[nz,nmodes])
   vhat_real = reshape(data2read(:,2),[nz,nmodes])
   what_real = reshape(data2read(:,3),[nz,nmodes])
   deallocate(data2read)
   
   allocate(uhat(nz,nmodes), vhat(nz,nmodes), what(nz,nmodes))

   uhat = uhat_real + imi*uhat_imag
   vhat = vhat_real + imi*vhat_imag
   what = what_real + imi*what_imag
   deallocate(uhat_real,vhat_real,what_real,uhat_imag,vhat_imag,what_imag)
 

   u = 0.d0
   v = 0.d0
   w = 0.d0

   do modeID = 1,nmodes
      ziter = 1
      do k = gp%xst(3),gp%xen(3)
         do j = 1,gp%xsz(2)
            do i = 1,gp%xsz(1)
               expfact = exp(imi*(ky(modeID)*y(1,j,1) + kx(modeID)*x(i,1,1)))
               uhatn = uhat(k,modeID)*cos(beta(modeID)*pi/180.d0) - vhat(k,modeID)*sin(beta(modeID)*pi/180.d0)
               vhatn = uhat(k,modeID)*sin(beta(modeID)*pi/180.d0) + vhat(k,modeID)*cos(beta(modeID)*pi/180.d0)
               whatn = what(k,modeID) 
               u(i,j,ziter) = u(i,j,ziter) + real(uhatn*expfact , rkind) 
               v(i,j,ziter) = v(i,j,ziter) + real(vhatn*expfact , rkind) 
               w(i,j,ziter) = w(i,j,ziter) + real(whatn*expfact , rkind)
            end do 
         end do 
         ziter = ziter + 1
      end do
   end do 

   deallocate(beta)
   deallocate(uhat, vhat, what)
   deallocate(kx, ky)

end subroutine 

end module 
