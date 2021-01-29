program getWakeCenterline
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two, half, one, zero
   use reductions, only: p_sum
   use mpi
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3, buff4
   real(rkind), dimension(:,:,:), allocatable :: omega1, omega2, omega3, x, y, z
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, fconv, fconv_p, buff1_p 
   real(rkind), dimension(:,:  ), allocatable :: data2write, raddist, fmask
   real(rkind), dimension(:    ), allocatable :: yc, zc
   integer,     dimension(:    ), allocatable :: tidx_arr
   real(rkind) :: dx, dy, dz, Re = 3000.d0, turbDiam = one, turbRadSq
   integer :: nx, ny, nz, RunID, TIDX, tcoun
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   character(len=4)    :: flabel
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0, numer, denom, ycen, zcen, xturb = 0.2_rkind, yturb, zturb, xst, xen
   integer :: ierr, NumericalSchemeVert = 1, i, j, k, jst, jen, kst, ken, shift_y, shift_z, ii2(2), ist, ien
   integer :: tsnap_st = 70000, dtsnaps = 40, nsnaps = 40, max_nyen, max_nzen, min_nyen, min_nzen, jturb, kturb, num_radii_y=6, num_radii_z=3
   logical :: isZPeriodic = .false.

   namelist /INPUT/ InputDir, OutputDir, RunID, nx, ny, nz, Lx, Ly, Lz, Re, NumericalSchemeVert, turbDiam, xturb,  yturb, zturb, tsnap_st, dtsnaps, nsnaps, num_radii_y, num_radii_z
   
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)          
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)

   dx = Lx/real(nx,rkind) 
   dy = Ly/real(ny,rkind) 
   dz = Lz/real(nz,rkind)

   ! Initialize the operator class
   call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NUmericalSchemeVert)

   ! Allocate all the needed memory 
   call ops%allocate3DField(u)
   call ops%allocate3DField(v)
   call ops%allocate3DField(w)
   
   call ops%allocate3DField(x)
   call ops%allocate3DField(y)
   call ops%allocate3DField(z)

   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4)
   
   call ops%allocate3DField(omega1)
   call ops%allocate3DField(omega2)
   call ops%allocate3DField(omega3)

   allocate(yc(nx),zc(nx))
   allocate(fconv_p(ny,nz,nx),fconv(nx,ny,nz))
   allocate(buff1_p(ny,nz,nx))
   allocate(data2write(nx,3))

   call ops%getGrid(x,y,z)
   turbRadSq = turbDiam*turbDiam/4.0d0
   shift_y = floor(turbDiam/two/dy)
   shift_z = floor(turbDiam/two/dz)

   jst = shift_y+1; jen = ny-shift_y
   kst = shift_z+1; ken = nz-shift_z

   xst = xturb-0.5*turbDiam;     ist = minloc(abs(x(:,1,1)-xst),1)
   xen = xturb+20.0*turbDiam;    ien = minloc(abs(x(:,1,1)-xen),1)

   !! trial 1
   !min_nyen = 1
   !min_nzen = 1
   !max_nyen = floor(3.0_rkind*ny/4.0_rkind)
   !max_nzen = floor(3.0_rkind*nz/4.0_rkind)

   ! trial 2
   jturb = minloc(abs(y(1,:,1)-yturb),1)
   kturb = minloc(abs(z(1,1,:)-zturb),1)
   min_nyen = jturb - num_radii_y*shift_y; min_nyen = max(min_nyen, 1)
   max_nyen = jturb + num_radii_y*shift_y; max_nyen = min(max_nyen, floor(3.0*ny/4.0))
   min_nzen = kturb - num_radii_z*shift_z; min_nzen = max(min_nzen, 1)
   max_nzen = kturb + num_radii_z*shift_z; max_nzen = min(max_nzen, floor(3.0*nz/4.0))
   print *, 'y limits = ', min_nyen, max_nyen
   print *, 'z limits = ', min_nzen, max_nzen

   allocate(tidx_arr(nsnaps))
   do tcoun = 1, nsnaps
       tidx_arr(tcoun) = tsnap_st + (tcoun-1) * dtsnaps
   enddo

   do tcoun = 1, size(tidx_arr)

     tidx = tidx_arr(tcoun)
     call message(0, "Reading fields for tid:", TIDX)
     call tic()
     call ops%ReadField3D(u,"uVel",TIDX)
     call ops%ReadField3D(v,"vVel",TIDX)
     call ops%ReadField3D(w,"wVel",TIDX)
  

     fconv_p = zero;     yc = zero;     zc = zero

     ! linear impulse
     !call ops%getCurl(u,v,w,omega1,omega2,omega3,0,0,0,0)
     !buff1 = half*(y*omega3 - z*omega2)

     !do i=1,size(x,1)
     !  numer = p_sum(sum(buff1(i,:,:)*y(i,:,:)))
     !  denom = p_sum(sum(buff1(i,:,:)))
     !  yc(i) = numer/(denom+1.0d-18)
 
     !  numer = p_sum(sum(buff1(i,:,:)*z(i,:,:)))
     !  zc(i) = numer/(denom+1.0d-18)
     !enddo

     ! calculate power
     buff1 = u*u + v*v + w*w
     buff1 = half*buff1*u

     do i = 1, nx
       do k = 1, nz
         do j = 1, ny
           buff1_p(j,k,i) = buff1(i,j,k)
         enddo
       enddo
     enddo

     do i = ist, ien
       do k = 1, nz
          buff1_p(:,k,i) = buff1_p(:,k,i) - sum(buff1_p(:,k,i))/real(ny, rkind)
       enddo
       do j = jst, jen
         do k = kst, ken
             fconv_p(j,k,i) = sum(buff1_p(j-shift_y:j+shift_y,k-shift_z:k+shift_z,i))
         enddo
       enddo
       ii2 = minloc(fconv_p(min_nyen:max_nyen,min_nzen:max_nzen,i))
       ii2(1) = ii2(1) + min_nyen-1;       ii2(2) = ii2(2) + min_nzen-1
       yc(i) = y(i,ii2(1),ii2(2))
       zc(i) = z(i,ii2(1),ii2(2))
       print '(3(i5,1x),2(f12.5,1x))', i, ii2, y(i,ii2(1),ii2(2)), z(i,ii2(1),ii2(2))
     enddo

     do k = 1, nz
       do j = 1, ny
         do i = 1, nx
           fconv(i,j,k) = fconv_p(j,k,i)
         enddo
       enddo
     enddo
     call ops%WriteField3D(fconv, "fcnv", tidx)

     !call ops%WriteField3D(buff2, "pcon", tidx)

     data2write(:,1) = x(:,1,1)
     data2write(:,2) = yc
     data2write(:,3) = zc
     write(flabel,"(A1,I3.3)") "c", tcoun
     !print *, filename
     call ops%WriteASCII_2D(data2write, flabel)
 
   enddo

   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

     ! go to y decomp
     ! cast all local arrarys s.t. u(ny, nzloc, nxloc)
     ! compute buff1_p with local spanwise mean subtracted
     ! compute sum in y. then
 
