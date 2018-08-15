program StratifiedShearLayerProfiles
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3, buff4, buff5, buff6
   real(rkind), dimension(:,:,:), allocatable :: x, y, z, ufluct, vfluct, w, T
   real(rkind) :: dx, dy, dz
   integer :: nx=256, ny=256, nz=256
   type(igrid_ops) :: ops
   character(len=clen) :: inputfile="/"
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   integer :: ierr, NumericalSchemeVert = 1
   integer :: i, j, k
   logical :: isZPeriodic = .false.
   real(rkind) :: cumerror, temp

   call MPI_Init(ierr)

   dx =     Lx/real(nx,rkind)
   dy =     Ly/real(ny,rkind)
   dz = two*Lz/real(nz,rkind)

   ! Initialize the operator class
   call ops%init(nx, ny, nz, dx, dy, dz, inputfile, inputfile, 0, isZPeriodic, NumericalSchemeVert)

   ! Allocate all the needed memory
   call ops%allocate3DField(x)
   call ops%allocate3DField(y)
   call ops%allocate3DField(z)
   call ops%allocate3DField(ufluct)
   call ops%allocate3DField(vfluct)
   call ops%allocate3DField(w)
   call ops%allocate3DField(T)

   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4)
   call ops%allocate3DField(buff5)
   call ops%allocate3DField(buff6)

   ! Initialise grid  
   do k = 1,ops%gp%xsz(3)
      do j = 1,ops%gp%xsz(2)
         do i = 1,ops%gp%xsz(1)
           x(i,j,k) = real(ops%gp%xst(1) - 1 + i - 1, rkind)*dx
           y(i,j,k) = real(ops%gp%xst(2) - 1 + j - 1, rkind)*dy
           z(i,j,k) = real(ops%gp%xst(3) - 1 + k - 1, rkind)*dz
          
           ufluct(i,j,k) = Sin(y(i,j,k)*two*pi/Ly)*Cos(z(i,j,k)*pi/Lz)
           vfluct(i,j,k) = Sin(x(i,j,k)*two*pi/Lx)*Cos(z(i,j,k)*pi/Lz)
           w(i,j,k) = Sin(z(i,j,k)*pi/Lz)
           T(i,j,k) = Cos(z(i,j,k)*pi/Lz)**2.d0
         end do
      end do
   end do




   ! T = Rib*(T - Tref)  ! Rescale Potential temperature to buoyancy variable: b

   call ops%getCurl(ufluct,vfluct,w, buff1,buff2,buff3,1,1,1,1)

   call ops%ComputeF_mvOmega(buff1,buff2,buff3,T,buff4,buff5,buff6)
   
     
  
   ! Check Solution
   cumerror = 0.d0

   do k = 1,ops%gp%xsz(3)
      do j = 1,ops%gp%xsz(2)
         do i = 1,ops%gp%xsz(1)

           temp = -4.d0 * pi**2.d0 * Lz * Cos(pi*z(i,j,k)/Lz) * Sin(2.d0 * pi * y(i,j,k)/Ly) / (Ly*Ly*Lz) 
           cumerror = cumerror + abs(buff4(i,j,k) - temp)

           !call message(0, "value1:", buff4(i,j,k))
           !call message(0, "answer1:", temp)

           temp = -4.d0 * pi**2.d0 * Lz * Cos(pi*z(i,j,k)/Lz) * Sin(2.d0 * pi *x(i,j,k)/Lx) / (Lx*Lx*Lz)
           cumerror = cumerror + abs(buff5(i,j,k) - temp)

           !call message(0, "value2:", buff5(i,j,k))
           !call message(0, "answer2:", temp)

           cumerror = cumerror + abs(buff6(i,j,k))
     
           !call message(0, "value3:", buff6(i,j,k))
           !call message(0, "answer3:", 0)

         end do
      end do
   end do

   call message(0, "Total Error:", cumerror)  

   call ops%destroy()
   call MPI_Finalize(ierr)


end program
