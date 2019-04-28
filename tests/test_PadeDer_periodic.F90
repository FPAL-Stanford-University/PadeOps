program test_padeder_periodic
   use kind_parameters, only: rkind, clen
   use decomp_2d 
   use spectralMod, only: spectral
   use mpi
   use constants, only: two, pi
   use PadeDerOps, only: Pade6stagg
   
   implicit none 

   real(rkind), dimension(:,:,:), allocatable :: x, y, z, xE, yE, zE, fE, fC, tmpC, tmpE, rbuffyC, rbuffzC
   complex(rkind), dimension(:,:,:), allocatable :: fhatC, fhatzC, fhatzE, fhatE, fhatzC2, fhatzE2
   integer :: nx = 4, ny = 4, nz = 32
   type(spectral), target :: spectC, spectE
   type(decomp_info) :: gpC, gpE
   type(decomp_info), pointer :: sp_gpC, sp_gpE
   type(Pade6stagg) :: derZ
   real(rkind) :: dx, dy, dz
   integer, parameter :: scheme = 2
   integer :: i, j, k, ierr

   
   call MPI_Init(ierr)
   call decomp_2d_init(nx, ny, nz, 1, 1)
   call get_decomp_info(gpC)
   call decomp_info_init(nx,ny,nz+1,gpE)
   
   allocate(x(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
   allocate(y(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
   allocate(z(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))

   allocate(xE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))
   allocate(yE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))
   allocate(zE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))

   dx = two*pi/real(nx,rkind)
   dy = two*pi/real(ny,rkind)
   dz = two*pi/real(nz,rkind)
   do k=1,nz
       do j=1,ny
           do i=1,nx
               x(i,j,k) = real(  i , rkind ) * dx
               y(i,j,k) = real(  j , rkind ) * dy
               z(i,j,k) = real(  k , rkind ) * dz + dz/two
           end do
       end do
   end do
   x = x - dx; y = y - dy; z = z - dz 
   zE(:,:,1:nz) = z - dz/2.d0; zE(:,:,nz+1) = zE(:,:,nz) + dz;
   xE(:,:,1:nz) = x; xE(:,:,nz+1) = xE(:,:,nz)
   yE(:,:,1:nz) = y; yE(:,:,nz+1) = yE(:,:,nz)

   allocate(fE(nx, ny, nz+1))
   allocate(fC(nx, ny, nz))
   allocate(tmpC(nx, ny, nz))
   allocate(tmpE(nx, ny, nz+1))


   fC = cos(x)*cos(y)*cos(z)
   fE = cos(xE)*cos(yE)*cos(zE)
   !fC = cos(z)
   !fE = cos(zE)

   call spectC%init("x", nx, ny, nz  , dx, dy, dz, &
                "four", "2/3rd", 2 , fixOddball=.false., exhaustiveFFT=.true., init_periodicInZ=.true.)
   call spectE%init("x", nx, ny, nz+1, dx, dy, dz, &
                "four", "2/3rd", 2 , fixOddball=.false., exhaustiveFFT=.true., init_periodicInZ=.false.)

   sp_gpC => spectC%spectdecomp
   sp_gpE => spectE%spectdecomp
   allocate(fhatC(sp_gpC%ysz(1), sp_gpC%ysz(2), sp_gpC%ysz(3)))
   allocate(fhatE(sp_gpE%ysz(1), sp_gpE%ysz(2), sp_gpE%ysz(3)))
   allocate(fhatzC(sp_gpC%zsz(1), sp_gpC%zsz(2), sp_gpC%zsz(3)))
   allocate(fhatzC2(sp_gpC%zsz(1), sp_gpC%zsz(2), sp_gpC%zsz(3)))
   allocate(fhatzE(sp_gpE%zsz(1), sp_gpE%zsz(2), sp_gpE%zsz(3)))
   allocate(fhatzE2(sp_gpE%zsz(1), sp_gpE%zsz(2), sp_gpE%zsz(3)))
   
   call derZ%init(gpC, sp_gpC, gpE, sp_gpE, dz, scheme, .true., spectC)

   ! ddz_C2C(fC)
   
   tmpC = fC
   call spectC%ddz_C2C_spect(tmpC)
   print*, "max error (ddz_C2C): ", maxval(tmpC - (-cos(x)*cos(y)*sin(z)))

   ! ddz_C2E(fC)
   call spectC%fft(fC, fhatC)
   call transpose_y_to_z(fhatC, fhatzC, sp_gpC)
   call derZ%ddz_C2E(fhatzC, fhatzE, 0, 0)
   call transpose_z_to_y(fhatzE, fhatE, sp_gpE)
   call spectE%ifft(fhatE, tmpE)
   print*, "max error (ddz_C2E_CMPLX):", maxval(abs(tmpE - (-cos(xE)*cos(yE)*sin(zE)) ))
   call derZ%ddz_C2E(fC,tmpE,0,0)
   print*, "max error (ddz_C2E_REAL):", maxval(abs(tmpE - (-cos(xE)*cos(yE)*sin(zE)) ))

   ! interp_C2E(fC)
   call spectC%fft(fC, fhatC)
   call transpose_y_to_z(fhatC, fhatzC, sp_gpC)
   call derZ%interpz_C2E(fhatzC, fhatzE, 0, 0)
   call transpose_z_to_y(fhatzE, fhatE, sp_gpE)
   call spectE%ifft(fhatE, tmpE)
   print*, "max error (interpz_C2E_CMPLX):", maxval(abs(tmpE - (cos(xE)*cos(yE)*cos(zE)) ))
   call derZ%interpz_C2E(fC,tmpE,0,0)
   print*, "max error (interpz_C2E_REAL):", maxval(abs(tmpE - (cos(xE)*cos(yE)*cos(zE)) ))


   ! d2dz2_C2C(fC)
   call spectC%fft(fC, fhatC)
   call transpose_y_to_z(fhatC, fhatzC, sp_gpC)
   call derZ%d2dz2_C2C(fhatzC, fhatzC2, 0, 0)
   call transpose_z_to_y(fhatzC2, fhatC, sp_gpC)
   call spectE%ifft(fhatC, tmpC)
   print*, "max error (d2dz2_C2C_CMPLX):", maxval(abs(tmpC - (-cos(x)*cos(y)*cos(z)) ))
   call derZ%d2dz2_C2C(fC, tmpC,0,0)
   print*, "max error (d2dz2_C2C_REAL):", maxval(abs(tmpC - (-cos(x)*cos(y)*cos(z)) ))



   ! ddz_E2C(fE)
   call spectE%fft(fE, fhatE)
   call transpose_y_to_z(fhatE, fhatzE, sp_gpE)
   call derZ%ddz_E2C(fhatzE, fhatzC, 0, 0)
   call transpose_z_to_y(fhatzC, fhatC, sp_gpC)
   call spectC%ifft(fhatC, tmpC)
   print*, "max error (ddz_E2C_REAL):", maxval(abs(tmpC - (-cos(x)*cos(y)*sin(z)) ))
   call derZ%ddz_E2C(fE, tmpC,0,0)
   print*, "max error (ddz_E2C_REAL):", maxval(abs(tmpC - (-cos(x)*cos(y)*sin(z)) ))

   ! interp_E2C(fE)
   call spectE%fft(fE, fhatE)
   call transpose_y_to_z(fhatE, fhatzE, sp_gpE)
   call derZ%interpz_E2C(fhatzE, fhatzC, 0, 0)
   call transpose_z_to_y(fhatzC, fhatC, sp_gpC)
   call spectC%ifft(fhatC, tmpC)
   print*, "max error (interpz_E2C_CMPLX):", maxval(abs(tmpC - (cos(x)*cos(y)*cos(z)) ))
   call derZ%interpz_E2C(fE, tmpC, 0, 0)
   print*, "max error (interpz_E2C_REAL):", maxval(abs(tmpC - (cos(x)*cos(y)*cos(z)) ))


   ! d2dz2_E2E(fE)
   call spectE%fft(fE, fhatE)
   call transpose_y_to_z(fhatE, fhatzE, sp_gpE)
   call derZ%d2dz2_E2E(fhatzE, fhatzE2, 0, 0)
   call transpose_z_to_y(fhatzE2, fhatE, sp_gpE)
   call spectE%ifft(fhatE, tmpE)
   print*, "max error (d2dz2_E2E_REAL):", maxval(abs(tmpE - (-cos(xE)*cos(yE)*cos(zE)) ))
   call derZ%d2dz2_E2E(fE, tmpE,0,0)
   print*, "max error (d2dz2_E2E_REAL):", maxval(abs(tmpE - (-cos(xE)*cos(yE)*cos(zE)) ))

   call derZ%destroy()
   nullify(sp_gpC, sp_gpE)
   call spectC%destroy()
   call spectE%destroy()

   deallocate(fhatC, fhatzC, fhatzC2)
   deallocate(fhatE, fhatzE, fhatzE2)
   call MPI_Finalize(ierr)

end program 
