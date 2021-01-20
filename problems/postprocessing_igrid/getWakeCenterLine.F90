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
   real(rkind) :: dx, dy, dz, Re = 3000.d0, turbDiam = one, turbRadSq
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile!, filename
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0, numer, denom, ycen, zcen, xturb, yturb, zturb
   integer :: ierr, tsnapshot = 0, NumericalSchemeVert = 1, i, j, k, jst, jen, kst, ken, shift_y, shift_z, ii2(2)
   logical :: isZPeriodic = .false. 

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tsnapshot, nx, ny, nz, Re, NumericalSchemeVert, turbDiam, xturb,  yturb, zturb
   
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


   tidx = tsnapshot
   call message(0, "Reading fields for tid:", TIDX)
   call tic()
   call ops%ReadField3D(u,"uVel",TIDX)
   call ops%ReadField3D(v,"vVel",TIDX)
   call ops%ReadField3D(w,"wVel",TIDX)
  
   allocate(yc(nx),zc(nx))
   !allocate(raddist(ny,nz),fconv(ny,nz))
   allocate(fconv_p(ny,nz,nx),fconv(nx,ny,nz))
   allocate(buff1_p(ny,nz,nx))


   fconv = zero; fconv_p = zero; buff1_p = zero

   ! linear impulse
   call ops%getGrid(x,y,z)
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

   
 
   turbRadSq = turbDiam*turbDiam/4.0d0


   do i = 1, nx
     do k = 1, nz
       do j = 1, ny
         buff1_p(j,k,i) = buff1(i,j,k)
       enddo
     enddo
   enddo

   shift_y = floor(turbDiam/two/dy)
   shift_z = floor(turbDiam/two/dz)

   jst = shift_y+1; jen = ny-shift_y
   kst = shift_z+1; ken = nz-shift_z

   do i = 1, nx
     do k = 1, nz
        buff1_p(:,k,i) = buff1_p(:,k,i) - sum(buff1_p(:,k,i))/real(ny, rkind)
     enddo
     do j = jst, jen
       do k = kst, ken
           fconv_p(j,k,i) = sum(buff1_p(j-shift_y:j+shift_y,k-shift_z:k+shift_z,i))
       enddo
     enddo
     ii2 = minloc(fconv_p)
     yc(i) = y(i,ii2(1),ii2(2))
     zc(i) = z(i,ii2(1),ii2(2))
     print *, i
   enddo

   do k = 1, nz
     do j = 1, ny
       do i = 1, nx
         fconv(i,j,k) = fconv_p(j,k,i)
       enddo
     enddo
   enddo
   call ops%WriteField3D(fconv, "fcnv", tidx)

   !buff2 = -100.0d0
   !kst = minloc(abs(z(1,1,:)-(zturb-two*turbDiam)),dim=1)
   !ken = minloc(abs(z(1,1,:)-(zturb+two*turbDiam)),dim=1)
   !jst = minloc(abs(y(1,:,1)-(yturb-two*turbDiam)),dim=1)
   !jen = minloc(abs(y(1,:,1)-(yturb+two*turbDiam)),dim=1)
   !print *, jst, jen
   !print *, kst, ken
   !print *, zturb, yturb, turbDiam
   !print *, minval(z), maxval(z)
   !print *, minval(y), maxval(y)
 
   !! sharp mask function
   !do k = kst, ken
   !  do j = jst, jen
   !    print *, j, k
   !    ycen = y(1,j,k); zcen = z(1,j,k)
   !    raddist(:,:) = (y(1,:,:) - ycen)**2 + (z(1,:,:) - zcen)**2
   !    where(raddist < turbRadSq)
   !      fmask = -one
   !    elsewhere
   !      fmask = zero
   !    endwhere
   !   
   !    do i = 1, nx 
   !      buff2(i,j,k) = p_sum(sum(buff1(i,:,:)*fmask))
   !    enddo

   !  enddo
   !enddo 

   !call ops%WriteField3D(buff1, "powd", tidx)
   call ops%WriteField3D(buff2, "pcon", tidx)
   

   allocate(data2write(nx,3))
   data2write(:,1) = x(:,1,1)
   data2write(:,2) = yc
   data2write(:,3) = zc
   !write(filename,"(A5,I6.6,A4)") "yczc_", tidx, ".dat"
   !print *, filename
   call ops%WriteASCII_2D(data2write, 'yczc')
 
   call message(0,"Done computing x impulse")
   call ops%WriteField3D(buff1, "impx", tidx)
   call message(0, "x impulse field written.")
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

