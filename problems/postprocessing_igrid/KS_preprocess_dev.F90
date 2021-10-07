program KS_preprocess_dev
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi
   use mpi
   use timer, only: tic, toc
   use exits, only: message, gracefulExit
   implicit none

   ! Input file variables
   integer :: nx, ny, nz, nxLES, nyLES, nzLES, nxQH, nyQH, nzQH
   integer :: tid_st, tid_en, tid_stride = 1
   real(rkind) :: Lx, Ly, Lz
   real(rkind) :: Lfact = 1.d0, disp_fact = 1.d0
   character(len=clen) :: inputdir, outputdir
   logical :: isZperiodic = .true.
   integer :: NumericalSchemeVert = 2
   integer :: RunID

   ! Other variables
   real(rkind) :: kco1, kco2
   real(rkind), dimension(:,:,:), allocatable :: u, v, w
   real(rkind), dimension(:,:,:), allocatable :: dudx, dudy, dudz
   real(rkind), dimension(:,:,:), allocatable :: dvdx, dvdy, dvdz
   real(rkind), dimension(:,:,:), allocatable :: dwdx, dwdy, dwdz
   type(igrid_ops) :: ops
   real(rkind) :: dx, dy, dz
   real(rkind) :: dxLES, dyLES, dzLES
   real(rkind) :: dxQH, dyQH, dzQH
   integer :: ierr, tid
   character(len=clen) :: inputfile

   namelist /INPUT/ nx, ny, nz, nxLES, nyLES, nzLES, nxQH, nyQH, nzQH, Lx, Ly, &
     Lz, Lfact, disp_fact, inputdir, outputdir, tid_st, tid_en, tid_stride, &
     RunID, isZperiodic, NumericalSchemeVert

   call MPI_Init(ierr)
   call GETARG(1,inputfile)
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)

   dx = Lx/nx
   dy = Ly/ny
   dz = Lz/nz
   dxLES = Lx/nxLES
   dyLES = Ly/nyLES
   dzLES = Lz/nzLES
   dxQH = Lx/nxQH
   dyQH = Ly/nyQH
   dzQH = Lz/nzQH

   ! Initialize the operator classes
   call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, &
     isZPeriodic, NUmericalSchemeVert, FFT3D = .true.)
   
   ! Allocate all the needed memory
   call ops%allocate3DField(u)
   call ops%allocate3DField(v)
   call ops%allocate3DField(w)

   call ops%allocate3DField(dudx)
   call ops%allocate3DField(dudy)
   call ops%allocate3DField(dudz)

   call ops%allocate3DField(dvdx)
   call ops%allocate3DField(dvdy)
   call ops%allocate3DField(dvdz)

   call ops%allocate3DField(dwdx)
   call ops%allocate3DField(dwdy)
   call ops%allocate3DField(dwdz)

   ! Loop over tids
   do tid = tid_st, tid_en, tid_stride
     ! Read velocity data
     call ops%ReadField3D(u,'uVel',tid)
     call ops%ReadField3D(v,'vVel',tid)
     call ops%ReadField3D(w,'wVel',tid)
     
     ! Filter velocity data
     call ops%initFilter(nxLES, nyLES, nzLES, 1)
     call ops%FilterField_inplace(u)
     call ops%FilterField_inplace(v)
     call ops%FilterField_inplace(w)

     call ops%WriteField3D(u,'uVel',tid)
     call ops%WriteField3D(v,'vVel',tid)
     call ops%WriteField3D(w,'wVel',tid)
     
     ! Compute gradients
     call ops%ddx(u,dudx)
     call ops%ddy(u,dudy)
     call ops%ddz(u,dudz)

     call ops%ddx(v,dvdx)
     call ops%ddy(v,dvdy)
     call ops%ddz(v,dvdz)

     call ops%ddx(w,dwdx)
     call ops%ddy(w,dwdy)
     call ops%ddz(w,dwdz)

     call ops%WriteField3D(dudx,'dudx',tid)
     call ops%WriteField3D(dudy,'dudy',tid)
     call ops%WriteField3D(dudz,'dudz',tid)
     
     call ops%WriteField3D(dvdx,'dvdx',tid)
     call ops%WriteField3D(dvdy,'dvdy',tid)
     call ops%WriteField3D(dvdz,'dvdz',tid)
     
     call ops%WriteField3D(dwdx,'dwdx',tid)
     call ops%WriteField3D(dwdy,'dwdy',tid)
     call ops%WriteField3D(dwdz,'dwdz',tid)
     
     call gracefulExit('Finished gradients',ierr)
     
     ! Downsample velocity data
     

     ! further downsample (for QH only)
     ! Save data
   end do
end program KS_preprocess_dev
