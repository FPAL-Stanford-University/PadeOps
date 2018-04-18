program HIT_AD_getRapidSlowPressure
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use timer, only: tic, toc
   use exits, only: message, gracefulExit
   use fof_mod, only: fof  
   implicit none
   real(rkind), dimension(:,:,:), allocatable :: uM, vM, wM 
   real(rkind), dimension(:,:,:), allocatable :: uF, vF, wF
   real(rkind), dimension(:,:,:), allocatable :: u, v, w
   real(rkind), dimension(:,:,:), allocatable :: dudx, dudy, dudz 
   real(rkind), dimension(:,:,:), allocatable :: dvdx, dvdy, dvdz 
   real(rkind), dimension(:,:,:), allocatable :: dwdx, dwdy, dwdz 
   
   real(rkind), dimension(:,:,:), allocatable :: dudxM, dudyM, dudzM 
   real(rkind), dimension(:,:,:), allocatable :: dvdxM, dvdyM, dvdzM 
   real(rkind), dimension(:,:,:), allocatable :: dwdxM, dwdyM, dwdzM 
   
   real(rkind), dimension(:,:,:), allocatable :: R11, R12, R13, R22, R23, R33
   real(rkind), dimension(:,:,:), allocatable :: rbuff1, rbuff2
   real(rkind), dimension(:,:,:), allocatable :: pslow, prapid, pfull 

   complex(rkind), dimension(:,:,:), allocatable :: cbuffzC

   real(rkind) :: dx, dy, dz
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir, FilterInfoDir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 10.d0*pi, Ly = 2.d0*pi, Lz = 2.d0*pi
   integer ::  ierr, tsnapshot, NumericalSchemeVert = 2
   logical :: isZPeriodic = .true. 
   integer :: MeanTIDX 
   integer :: topBC = 0, botBC = 0
   type(fof) :: fof1, fof2

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tsnapshot, nx, ny, nz, MeanTIDX, FilterInfoDir

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
   
   call ops%alloc_cbuffz(cbuffzC)
   call fof1%init(RunID, FilterInfoDir , outputdir, 1, ops%spect, ops%cbuffy1, cbuffzC, isZPeriodic, ops%gp)
   call fof2%init(RunID, FilterInfoDir , outputdir, 2, ops%spect, ops%cbuffy1, cbuffzC, isZPeriodic, ops%gp)

   ! Allocate all the needed memory 
   call ops%allocate3DField(u)
   call ops%allocate3DField(v)
   call ops%allocate3DField(w)
   
   call ops%allocate3DField(uM)
   call ops%allocate3DField(vM)
   call ops%allocate3DField(wM)
   
   call ops%allocate3DField(uF)
   call ops%allocate3DField(vF)
   call ops%allocate3DField(wF)
   
   call ops%allocate3DField(R11)
   call ops%allocate3DField(R12)
   call ops%allocate3DField(R13)
   call ops%allocate3DField(R22)
   call ops%allocate3DField(R23)
   call ops%allocate3DField(R33)
  
   call ops%allocate3DField(dudx)
   call ops%allocate3DField(dudy)
   call ops%allocate3DField(dudz)
   call ops%allocate3DField(dvdx)
   call ops%allocate3DField(dvdy)
   call ops%allocate3DField(dvdz)
   call ops%allocate3DField(dwdx)
   call ops%allocate3DField(dwdy)
   call ops%allocate3DField(dwdz)

   call ops%allocate3DField(dudxM)
   call ops%allocate3DField(dudyM)
   call ops%allocate3DField(dudzM)
   call ops%allocate3DField(dvdxM)
   call ops%allocate3DField(dvdyM)
   call ops%allocate3DField(dvdzM)
   call ops%allocate3DField(dwdxM)
   call ops%allocate3DField(dwdyM)
   call ops%allocate3DField(dwdzM)

   call ops%allocate3DField(rbuff1)
   call ops%allocate3DField(rbuff2)
   call ops%allocate3DField(pslow)
   call ops%allocate3DField(prapid)
   call ops%allocate3DField(pfull)
   
   call ops%ReadField3D(uM,"uMmn",MeanTIDX)
   call ops%ReadField3D(vM,"vMmn",MeanTIDX)
   call ops%ReadField3D(wM,"wMmn",MeanTIDX)

   call ops%ReadField3D(R11,"R11m",MeanTIDX)
   call ops%ReadField3D(R12,"R12m",MeanTIDX)
   call ops%ReadField3D(R13,"R13m",MeanTIDX)
   call ops%ReadField3D(R22,"R22m",MeanTIDX)
   call ops%ReadField3D(R23,"R23m",MeanTIDX)
   call ops%ReadField3D(R33,"R33m",MeanTIDX)

   ! Get mean flow gradients 
   call ops%GetGradient(uM, dudxM, dudyM, dudzM, botBC, topBC)
   call ops%GetGradient(vM, dvdxM, dvdyM, dvdzM, botBC, topBC)
   call ops%GetGradient(wM, dwdxM, dwdyM, dwdzM, botBC, topBC)
   
  
   call ops%initPoissonSolver(dx, dy, dz)

   tidx = tsnapshot
   call message(0, "Reading fields for tid:", TIDX)
   call tic()
   call ops%ReadField3D(u,"uVel",TIDX)
   call ops%ReadField3D(v,"vVel",TIDX)
   call ops%ReadField3D(w,"wVel",TIDX)
   
   uF = u - uM
   vF = v - vM
   wF = w - wM
   
   call ops%GetGradient(uF, dudx, dudy, dudz, botBC, topBC)
   call ops%GetGradient(vF, dvdx, dvdy, dvdz, botBC, topBC)
   call ops%GetGradient(wF, dwdx, dwdy, dwdz, botBC, topBC)

   prapid = -2.d0*(dudxM*dudx + dudyM*dvdx + dudzM*dwdx &
              &  + dvdxM*dudy + dvdyM*dvdy + dvdzM*dwdy & 
              &  + dwdxM*dudz + dwdyM*dvdz + dwdzM*dwdz )  

   call ops%PoissonSolvePeriodic_inplace(prapid,dealiasRHS=.true.)

   rbuff1 = R11 - uF*uF
   call ops%ddx(rbuff1,rbuff2)
   call ops%ddx(rbuff2,rbuff1)
   pslow = rbuff1

   rbuff1 = 2.d0*(R12 - uF*vF)
   call ops%ddx(rbuff1,rbuff2)
   call ops%ddy(rbuff2,rbuff1)
   pslow = pslow + rbuff1

   rbuff1 = 2.d0*(R13 - uF*wF)
   call ops%ddx(rbuff1,rbuff2)
   call ops%ddz(rbuff2,rbuff1, botBC, topBC)
   pslow = pslow + rbuff1

   rbuff1 = R22 - vF*vF
   call ops%ddy(rbuff1,rbuff2)
   call ops%ddy(rbuff2,rbuff1)
   pslow = pslow + rbuff1

   rbuff1 = 2.d0*(R23 - vF*wF)
   call ops%ddy(rbuff1,rbuff2)
   call ops%ddz(rbuff2,rbuff1, botBC, topBC)
   pslow = pslow + rbuff1

   rbuff1 = R33 - wF*wF
   call ops%ddz(rbuff1,rbuff2, botBC, topBC)
   call ops%ddz(rbuff2,rbuff1, botBC, topBC)
   pslow = pslow + rbuff1

   call ops%PoissonSolvePeriodic_inplace(pslow,dealiasRHS=.true.)

   pfull = prapid + pslow
   
   call toc()

   call ops%writeField3D(prapid,"PRap",tidx)
   call ops%writeField3D(pslow,"PSlo",tidx)
   call ops%writeField3D(pfull,"PFlt",tidx)


   ! Try Filtering 
   call fof1%filter_Real2Real(uF,rbuff1)
   call ops%writeField3D(rbuff1,"uSc1",tidx)
   call fof2%filter_Real2Real(uF,rbuff2)
   rbuff2 = rbuff2 - rbuff1
   call ops%writeField3D(rbuff2,"uSc2",tidx)
   rbuff2 = rbuff2 + rbuff1
   rbuff1 = uF - (rbuff2)
   call ops%writeField3D(rbuff1,"uSc3",tidx)
   
   call fof1%filter_Real2Real(vF,rbuff1)
   call ops%writeField3D(rbuff1,"vSc1",tidx)
   call fof2%filter_Real2Real(vF,rbuff2)
   rbuff2 = rbuff2 - rbuff1
   call ops%writeField3D(rbuff2,"vSc2",tidx)
   rbuff2 = rbuff2 + rbuff1
   rbuff1 = vF - (rbuff2)
   call ops%writeField3D(rbuff1,"vSc3",tidx)
   
   call fof1%filter_Real2Real(wF,rbuff1)
   call ops%writeField3D(rbuff1,"wSc1",tidx)
   call fof2%filter_Real2Real(wF,rbuff2)
   rbuff2 = rbuff2 - rbuff1
   call ops%writeField3D(rbuff2,"wSc2",tidx)
   rbuff2 = rbuff2 + rbuff1
   rbuff1 = wF - (rbuff2)
   call ops%writeField3D(rbuff1,"wSc3",tidx)
   
   call fof1%filter_Real2Real(pfull,rbuff1)
   call ops%writeField3D(rbuff1,"pSc1",tidx)
   call fof2%filter_Real2Real(pfull,rbuff2)
   rbuff2 = rbuff2 - rbuff1
   call ops%writeField3D(rbuff2,"pSc2",tidx)
   rbuff2 = rbuff2 + rbuff1
   rbuff1 = pfull - (rbuff2)
   call ops%writeField3D(rbuff1,"pSc3",tidx)
   
   
   
   
   call message(0,"Done writing filtered data")

   call ops%destroy()
   call MPI_Finalize(ierr)           

end program 

