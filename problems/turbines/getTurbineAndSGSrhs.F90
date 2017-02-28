#include "getTurbineAndSGSrhs_files/initialize.F90"       

program getTurbineAndSGSrhs
   use kind_parameters, only: rkind, clen
   use mpi
   use decomp_2d
   use spectralMod, only: spectral
   use initprocedures, only: meshgen,readVisualizationFile,writeVisualizationFile
   use turbineMod, only: turbineArray
   use sgsmod, only: sgs
   use PadeDerOps, only: Pade6Stagg

   implicit none
   complex(rkind), parameter :: czero = dcmplx(0.d0, 0.d0)
   type(turbineArray) :: turbArray
   integer :: ierr, nx, ny, nz, tid, rid, ioUnit
   character(len=clen) :: inputfile, inputdir, outputdir
   real(rkind) :: z0, Lx, Ly, Lz, z0init, ncWall, Cs, dx, dy, dz, dt, Pr = 1.d0
   integer :: SGSModelID, wallMType
   logical :: useSGS, useDynamicProcedure, useSGSClipping, useVerticalTfilter, useWallDamping
   type(decomp_info) :: gpC, gpE
   type(decomp_info), pointer :: sp_gpC, sp_gpE
   type(spectral), target :: spectC, spectE
   type(sgs) :: SGSmodel
   real(rkind),    dimension(:,:,:,:), allocatable :: mesh
   real(rkind),    dimension(:,:,:,:), allocatable :: rbuffxC
   complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyC, cbuffyE, cbuffzC, cbuffzE
   type(Pade6Stagg) :: Pade6opZ
   real(rkind), dimension(8) :: inst_horz_avg_turb

   real(rkind), dimension(:,:,:), allocatable :: u, v, w, wC, uE, vE, u_rhsR, v_rhsR, w_rhsR
   complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what, u_rhs, v_rhs, w_rhs, wC_rhs 

   namelist /INPUT/ nx, ny, nz, inputdir, outputdir, tid, rid
   namelist /PBLINPUT/ Lx, Ly, Lz, z0init
   namelist /LES/useSGS, SGSModelID, useDynamicProcedure, useSGSClipping, useVerticalTfilter, useWallDamping, ncWall, Cs
   namelist /WALLMODEL/ wallMType, z0

   ! Start MPI
   call MPI_Init(ierr)
   
   ! Read command line argument
   call GETARG(1,inputfile)          !<-- Get the location of the input file

   ! Read input file
   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=ioUnit, NML=INPUT)
   read(unit=ioUnit, NML=LES)
   read(unit=ioUnit, NML=WALLMODEL)
   read(unit=ioUnit, NML=PBLINPUT)
   close(ioUnit)


   ! Start decomp_2d
   call decomp_2d_init(nx, ny, nz, 0, 0)
   call get_decomp_info(gpC)
   call decomp_info_init(nx, ny, nz+1, gpE)

   ! Allocate mesh 
   allocate(mesh(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3),3))
   dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
   call meshgen(gpC, dx, dy, dz, mesh, inputfile)

   ! Initialize spectral 
   call spectC%init("x", nx, ny, nz  , dx, dy, dz, "four", "2/3rd", 2,.false.)
   call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", "2/3rd", 2,.false.)
   sp_gpC => spectC%spectdecomp; sp_gpE => spectE%spectdecomp

   ! Allocate buffers
   allocate(rbuffxC(   gpC%xsz(1),   gpC%xsz(2),   gpC%xsz(3),2))
   allocate(cbuffyC(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3),2))
   allocate(cbuffyE(sp_gpE%ysz(1),sp_gpE%ysz(2),sp_gpE%ysz(3),2))
   allocate(cbuffzC(sp_gpC%zsz(1),sp_gpC%zsz(2),sp_gpC%zsz(3),2))
   allocate(cbuffzE(sp_gpE%zsz(1),sp_gpE%zsz(2),sp_gpE%zsz(3),2))
   allocate(      u(   gpC%xsz(1),   gpC%xsz(2),   gpC%xsz(3)  ))
   allocate(      v(   gpC%xsz(1),   gpC%xsz(2),   gpC%xsz(3)  ))
   allocate(     wC(   gpC%xsz(1),   gpC%xsz(2),   gpC%xsz(3)  ))
   allocate(      w(   gpE%xsz(1),   gpE%xsz(2),   gpE%xsz(3)  ))
   allocate( u_rhsR(   gpC%xsz(1),   gpC%xsz(2),   gpC%xsz(3)  ))
   allocate( v_rhsR(   gpC%xsz(1),   gpC%xsz(2),   gpC%xsz(3)  ))
   allocate( w_rhsR(   gpC%xsz(1),   gpC%xsz(2),   gpC%xsz(3)  ))

   call spectC%alloc_r2c_out(u_rhs)
   call spectC%alloc_r2c_out(v_rhs)
   call spectE%alloc_r2c_out(w_rhs)
   call spectC%alloc_r2c_out(wC_rhs)
   call spectC%alloc_r2c_out(uhat )
   call spectC%alloc_r2c_out(vhat )
   call spectE%alloc_r2c_out(what )


   ! Initialize turbineMod 
   call turbArray%init(inputFile, gpC, gpE, spectC, spectE, rbuffxC, cbuffyC, cbuffyE, cbuffzC, cbuffzE, mesh, dx, dy, dz) 


   ! Initialize SGS model
   call SGSmodel%init(SGSModelID, spectC, spectE, gpC, gpE, dx, dy, dz, useDynamicProcedure, useSGSclipping, mesh(:,:,:,3), z0, &
       .true., WallMType, useVerticalTfilter, Pr, useWallDamping, nCWall, Cs, .true., .false. )


   ! Initialize the derivative derived type
   call Pade6OpZ%init(gpC,sp_gpC,dz)


   ! Read velocity fields
   call readVisualizationFile(tid, rid, u, v, wC, gpC, inputdir)


   ! Interpolate the fields 


   ! compute duidxj_tensor



   ! Compute wind turbine RHS
   u_rhs = czero; v_rhs = czero; w_rhs = czero
   call turbArray%getForceRHS(dt,u,v,wC,u_rhs,v_rhs,w_rhs,inst_horz_avg_turb)
   print*, shape(w_rhs)
   print*, shape(cbuffzE)
   call transpose_y_to_z(w_rhs, cbuffzE(:,:,:,1), sp_gpE)
   call Pade6OpZ%interpz_E2C(cbuffzE(:,:,:,1), cbuffzC(:,:,:,1),-1,-1)
   call transpose_z_to_y(cbuffzC(:,:,:,1),wC_rhs, sp_gpC)
   print*, 1
   call spectC%ifft(u_rhs , u_rhsR) 
   call spectC%ifft(v_rhs , v_rhsR) 
   print*, 2
   call spectC%ifft(wC_rhs, w_rhsR) 

   print*, 3
   call writeVisualizationFile(tid, rid, u_rhsR, gpC, "xtrb", OutputDir)
   call writeVisualizationFile(tid, rid, v_rhsR, gpC, "ytrb", OutputDir)
   call writeVisualizationFile(tid, rid, w_rhsR, gpC, "ztrb", OutputDir)




! write wind turbine RHS to disk


! COmpute SGS model RHS


! write SGS model RHS to disk


! deallocate memory 
   deallocate(mesh)

! End MPI
   call MPI_Finalize(ierr)           !<-- Terminate MPI 



end program


