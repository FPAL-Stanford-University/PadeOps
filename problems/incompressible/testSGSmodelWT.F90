program testSGSmodelWT
   use mpi
   use kind_parameters,  only: rkind, clen
   use timer, only: tic, toc
   use sgsmod_igrid, only: sgs_igrid
   use PadeDerOps, only: Pade6stagg
   use spectralMod, only: spectral
   use decomp_2d 
   use decomp_2d_io
   use turbineMod, only: turbineArray

   complex(rkind), dimension(:,:,:), allocatable :: uhatC, vhatC, whatE, uhatE, vhatE, whatC, ThatC
   real(rkind), dimension(:,:,:), allocatable :: uC, vC, wC, uE, vE, wC, fbody
   real(rkind), dimension(:,:,:,:), allocatable :: duidxjE, duidxjC
   complex(rkind), dimension(:,:,:,:), allocatable :: duidxjEhat
   type(sgs_igrid) :: sgs
    type(turbineArray) :: turbArray

   real(rkind), parameter :: Re = 1.d3, Fr = 1.d10
   real(rkind), parameter :: Tsurf = 1.d0, ThetaRef = 1.d0
   real(rkind) :: dx, dy, dz, Lx, Ly, Lz
   real(rkind), dimension(:,:,:,:), allocatable :: mesh
   real(rkind), dimension(:,:,:), allocatable ::  zMeshE
   type(spectral)  :: spectE, spectC
   type(decomp_info) :: gpC, gpE
   type(Pade6stagg) :: PadeDer
   integer :: ierr
   
   character(len=clen) :: inputdir, outputdir 
   integer :: nx, ny, nz
   real(rkind) :: z0init
   namelist /INPUT/ Lx, Ly, Lz, z0init, outputdir, inputdir, nx, ny, nz

   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
   read(unit=ioUnit, NML=PBLINPUT)
   close(ioUnit)    

   !Lx = two*pi; Ly = two*pi; Lz = one

   ! Do MPI stuff
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)         

   ! Do file IO - input file


   call decomp_2d_init(nx, ny, nz, 0, 0)
   call get_decomp_info(gpC)

   ! Allocate memory


   ! Allocate Buffers

   ! Create Mesh


   ! Read in the fields (from restart files)

   ! Initialize spectral


   ! Initialize Padeder


   ! Initialize sgs
   call sgs%init(this, gpC, gpE, spectC, spectE, dx, dy, dz, inputfile, zMeshE, zMeshC, fBody, computeFbody, cbuffyC, cbuffzC, rbuffxC, rbuffyC, rbuffzC, rbuffyE, rbuffzE, Tsurf, ThetaRef, Fr, Re)

   ! Initialize WT
   call turbArray%init(inputFile, gpC, gpE, spectC, spectE, rbuffxC, cbuffyC, cbuffyE, cbuffzC, cbuffzE, mesh, dx, dy, dz) 


   ! Interpolations
   
   
   ! duidxj
   

   ! Get RHS WT

   ! Get RHS constant force
   fbody(:,:,:,1) = fbody(:,:,:,1) + 1.d0 
   
   ! get tau_sgs



   ! deallocate everything











end program 
