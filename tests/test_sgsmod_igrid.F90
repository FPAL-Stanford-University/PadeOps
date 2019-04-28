program test_sgsmod_igrid
   use kind_parameters, only: rkind, clen
   use mpi
   use decomp_2d
   use decomp_2d_io
   use exits, only: message 
   use spectralMod, only: spectral
   use constants, only: pi, two
   use PadeDerOps, only: Pade6stagg
   use sgsmod_igrid, only: sgs_igrid 
   use timer, only: tic, toc

   implicit none
   real   (rkind), dimension(:,:,:), allocatable :: u, v, w, wC, T 
   real   (rkind), dimension(:,:,:,:), allocatable :: duidxjC, duidxjE, rbuffxC, rbuffxE, rbuffyC, rbuffzC, rbuffzE, rbuffyE
   complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyE, cbuffyC, cbuffzE, cbuffzC
   real(rkind),    dimension(:,:,:), allocatable :: fbody_x, fbody_y, fbody_z, dTdxC, dTdyC, dTdzC, dTdzE
   real(rkind),    dimension(:),     allocatable :: zMeshE, zMeshC
   complex(rkind),    dimension(:,:,:), allocatable :: uhatC, vhatC, whatC, ThatC, whatE, urhs, vrhs, wrhs, Trhs
   logical, parameter :: isPeriodic = .false.
   integer, parameter :: scheme = 1

   real(rkind), dimension(:,:,:), pointer :: nuSGS, kappaSGS, tau13, tau23, q1, q2, q3
   real(rkind), dimension(:,:,:,:), pointer :: tauSGS_ij

   integer, parameter :: nx = 128, ny = 128, nz = 128
   real(rkind), parameter :: dx = 9.d0*pi/real(nx,rkind), dy = 9.d0*pi/real(ny,rkind), dz = two*8.d0/real(nz,rkind) 
   type(decomp_info)          :: gp, gpE
   type(decomp_info), pointer :: sp_gp, sp_gpE
   integer :: ierr 
   type(spectral), allocatable, target :: spectC, spectE
   type(Pade6stagg), allocatable :: der
   type(sgs_igrid), allocatable :: sgsmodel
   integer :: botBC_temp = 1
   logical :: computeFbody = .false.  
   character(len=clen) :: inputfile = "/work/04076/tg833754/stampede2/me461/testing_data/inputsgs.dat"
   real(rkind), dimension(:,:,:), allocatable :: dTdxE, dTdyE

   call mpi_init(ierr)
   call decomp_2d_init(nx,ny,nz,0,0)
   call get_decomp_info(gp)
   call decomp_info_init(nx,ny,nz+1, gpE)

   allocate(spectC, spectE)
   call spectC%init("x", nx, ny, nz  , dx, dy, dz, "four", "2/3rd", 2,fixOddball=.false., init_periodicInZ=.false. )
   call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", "2/3rd", 2,fixOddball=.false., init_periodicInZ=.false.)
   sp_gp  => spectC%spectdecomp
   sp_gpE => spectE%spectdecomp

   allocate(u (gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(v (gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(wC(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(w (gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))
   allocate(T (gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(duidxjC(gp %xsz(1), gp %xsz(2), gp %xsz(3),9))
   allocate(duidxjE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3),9))
   
   allocate(dTdxC(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(dTdyC(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(dTdzC(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(dTdzE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))

   allocate(fbody_x(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(fbody_y(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(fbody_z(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(zMeshE (gpE%xsz(3)))
   allocate(zMeshC (gp %xsz(3)))
   
   allocate(rbuffxC(gp %xsz(1), gp %xsz(2), gp %xsz(3),3))
   allocate(rbuffyC(gp %ysz(1), gp %ysz(2), gp %ysz(3),3))
   allocate(rbuffzC(gp %zsz(1), gp %zsz(2), gp %zsz(3),3))
   
   allocate(rbuffxE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3),3))
   allocate(rbuffyE(gpE%ysz(1), gpE%ysz(2), gpE%ysz(3),3))
   allocate(rbuffzE(gpE%zsz(1), gpE%zsz(2), gpE%zsz(3),3))

   allocate(cbuffyC(sp_gp %ysz(1), sp_gp %ysz(2), sp_gp %ysz(3),3))
   allocate(cbuffyE(sp_gpE%ysz(1), sp_gpE%ysz(2), sp_gpE%ysz(3),3))

   allocate(cbuffzC(sp_gp %zsz(1), sp_gp %zsz(2), sp_gp %zsz(3),3))
   allocate(cbuffzE(sp_gpE%zsz(1), sp_gpE%zsz(2), sp_gpE%zsz(3),3))
   
   allocate(uhatC(sp_gp %ysz(1), sp_gp %ysz(2), sp_gp %ysz(3)))
   allocate(vhatC(sp_gp %ysz(1), sp_gp %ysz(2), sp_gp %ysz(3)))
   allocate(whatC(sp_gp %ysz(1), sp_gp %ysz(2), sp_gp %ysz(3)))
   allocate(whatE(sp_gpE%ysz(1), sp_gpE%ysz(2), sp_gpE%ysz(3)))
   allocate(ThatC(sp_gp %ysz(1), sp_gp %ysz(2), sp_gp %ysz(3)))
   
   allocate(urhs(sp_gp %ysz(1), sp_gp %ysz(2), sp_gp %ysz(3)))
   allocate(vrhs(sp_gp %ysz(1), sp_gp %ysz(2), sp_gp %ysz(3)))
   allocate(wrhs(sp_gpE%ysz(1), sp_gpE%ysz(2), sp_gpE%ysz(3)))
   allocate(Trhs(sp_gp %ysz(1), sp_gp %ysz(2), sp_gp %ysz(3)))
   
   call decomp_2d_read_one(1,u ,'/work/04076/tg833754/stampede2/me461/testing_data/RESTART_Run01_u.001263',gp )
   call decomp_2d_read_one(1,v ,'/work/04076/tg833754/stampede2/me461/testing_data/RESTART_Run01_v.001263',gp )
   call decomp_2d_read_one(1,w ,'/work/04076/tg833754/stampede2/me461/testing_data/RESTART_Run01_w.001263',gpE)
   call decomp_2d_read_one(1,T ,'/work/04076/tg833754/stampede2/me461/testing_data/RESTART_Run01_T.001263',gp )
   
   call message(0, "Testing data read successfully.")


   allocate(der)
   call der%init( gp, spectC%spectdecomp, gpE, spectE%spectdecomp, dz, scheme, isPeriodic, spectC)

   allocate(sgsmodel) 
   call message(0, "Now initializing the SGS model.")
   call sgsmodel%init( gp , gpE, spectC, spectE, dx, dy, dz, inputfile, zMeshE, zMeshC, fBody_x, fBody_y, fBody_z, &
                        & computeFbody, der, cbuffyC, cbuffzC,     cbuffyE, cbuffzE, rbuffxC, rbuffyC, rbuffzC,&
                        & rbuffyE, rbuffzE, 100.d0, 100.d0, 0.0d0, 1.d0, 3000.d0, .false., .true., botBC_temp, .false.)

   call message(0, "All operators initiatized successfully.")

   ! STEP 1: Compute uhatC, vhatC, whatC, ThatC
   call spectC%fft(u, uhatC)
   call spectC%fft(v, vhatC)
   call spectE%fft(w, whatE)
   call spectC%fft(T, ThatC)
   call transpose_y_to_z(whatE, cbuffzE(:,:,:,1), sp_gpE)
   call der%interpz_E2C(cbuffzE(:,:,:,1), cbuffzC(:,:,:,1), -1, -1)
   call transpose_z_to_y(cbuffzC(:,:,:,1), whatC, sp_gp )
   call spectC%ifft(whatC, wC)

   allocate(dTdxE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
   allocate(dTdyE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))

   ! STEP 2: Compute dudxC, dudyC, dudzC, dudxE, dudyE, dudzE
   call spectC%mtimes_ik1_oop(uhatC, cbuffyC(:,:,:,1))
   call spectC%ifft(cbuffyC(:,:,:,1), duidxjC(:,:,:,1))
   call spectC%mtimes_ik2_oop(uhatC, cbuffyC(:,:,:,1))
   call spectC%ifft(cbuffyC(:,:,:,1), duidxjC(:,:,:,2))
   call transpose_y_to_z(uhatC, cbuffzC(:,:,:,1), sp_gp)
   call der%ddz_C2E(cbuffzC(:,:,:,1),cbuffzE(:,:,:,1),1,1)
   call transpose_z_to_y(cbuffzE(:,:,:,1),cbuffyE(:,:,:,1),sp_gpE)
   call spectE%ifft(cbuffyE(:,:,:,1),duidxjE(:,:,:,3))
   call der%interpz_C2E(cbuffzC(:,:,:,1),cbuffzE(:,:,:,1),1,1)
   call der%ddz_E2C(cbuffzE(:,:,:,1),cbuffzC(:,:,:,2),1,1)
   call transpose_z_to_y(cbuffzC(:,:,:,2),cbuffyC(:,:,:,1),sp_gp)
   call spectC%ifft(cbuffyC(:,:,:,1),duidxjC(:,:,:,3))
   call transpose_z_to_y(cbuffzE(:,:,:,1),cbuffyE(:,:,:,1),sp_gpE)
   call spectE%mtimes_ik1_oop(cbuffyE(:,:,:,1), cbuffyE(:,:,:,2))
   call spectE%ifft(cbuffyE(:,:,:,2), duidxjE(:,:,:,1))
   call spectE%mtimes_ik2_oop(cbuffyE(:,:,:,1), cbuffyE(:,:,:,2))
   call spectE%ifft(cbuffyE(:,:,:,2), duidxjE(:,:,:,2))

   ! STEP 3: Compute dvdxC, dvdyC, dvdzC, dvdxE, dvdyE, dvdzE
   call spectC%mtimes_ik1_oop(vhatC, cbuffyC(:,:,:,1))
   call spectC%ifft(cbuffyC(:,:,:,1), duidxjC(:,:,:,4))
   call spectC%mtimes_ik2_oop(vhatC, cbuffyC(:,:,:,1))
   call spectC%ifft(cbuffyC(:,:,:,1), duidxjC(:,:,:,5))
   call transpose_y_to_z(vhatC, cbuffzC(:,:,:,1), sp_gp)
   call der%ddz_C2E(cbuffzC(:,:,:,1),cbuffzE(:,:,:,1),1,1)
   call transpose_z_to_y(cbuffzE(:,:,:,1),cbuffyE(:,:,:,1),sp_gpE)
   call spectE%ifft(cbuffyE(:,:,:,1),duidxjE(:,:,:,6))
   call der%interpz_C2E(cbuffzC(:,:,:,1),cbuffzE(:,:,:,1),1,1)
   call der%ddz_E2C(cbuffzE(:,:,:,1),cbuffzC(:,:,:,2),1,1)
   call transpose_z_to_y(cbuffzC(:,:,:,2),cbuffyC(:,:,:,1),sp_gp)
   call spectC%ifft(cbuffyC(:,:,:,1),duidxjC(:,:,:,6))
   call transpose_z_to_y(cbuffzE(:,:,:,1),cbuffyE(:,:,:,1),sp_gpE)
   call spectE%mtimes_ik1_oop(cbuffyE(:,:,:,1), cbuffyE(:,:,:,2))
   call spectE%ifft(cbuffyE(:,:,:,2), duidxjE(:,:,:,4))
   call spectE%mtimes_ik2_oop(cbuffyE(:,:,:,1), cbuffyE(:,:,:,2))
   call spectE%ifft(cbuffyE(:,:,:,2), duidxjE(:,:,:,5))

   ! STEP 4: Compute dwdxC, dwdyC, dwdzC, dwdxE, dwdyE, dwdzE
   call spectC%mtimes_ik1_oop(whatC, cbuffyC(:,:,:,1))
   call spectC%ifft(cbuffyC(:,:,:,1), duidxjC(:,:,:,7))
   call spectC%mtimes_ik2_oop(whatC, cbuffyC(:,:,:,1))
   call spectC%ifft(cbuffyC(:,:,:,1), duidxjC(:,:,:,8))
   call transpose_y_to_z(whatC, cbuffzC(:,:,:,1), sp_gp)
   call der%ddz_C2E(cbuffzC(:,:,:,1),cbuffzE(:,:,:,1),1,1)
   call transpose_z_to_y(cbuffzE(:,:,:,1),cbuffyE(:,:,:,1),sp_gpE)
   call spectE%ifft(cbuffyE(:,:,:,1),duidxjE(:,:,:,9))
   call der%interpz_C2E(cbuffzC(:,:,:,1),cbuffzE(:,:,:,1),1,1)
   call der%ddz_E2C(cbuffzE(:,:,:,1),cbuffzC(:,:,:,2),1,1)
   call transpose_z_to_y(cbuffzC(:,:,:,2),cbuffyC(:,:,:,1),sp_gp)
   call spectC%ifft(cbuffyC(:,:,:,1),duidxjC(:,:,:,9))
   call transpose_z_to_y(cbuffzE(:,:,:,1),cbuffyE(:,:,:,1),sp_gpE)
   call spectE%mtimes_ik1_oop(cbuffyE(:,:,:,1), cbuffyE(:,:,:,2))
   call spectE%ifft(cbuffyE(:,:,:,2), duidxjE(:,:,:,7))
   call spectE%mtimes_ik2_oop(cbuffyE(:,:,:,1), cbuffyE(:,:,:,2))
   call spectE%ifft(cbuffyE(:,:,:,2), duidxjE(:,:,:,8))


   ! STEP 5: Compute dTdxC, dTdyC, dTdzC, dTdzE
   call spectC%mtimes_ik1_oop(ThatC, cbuffyC(:,:,:,1))
   call spectC%ifft(cbuffyC(:,:,:,1),dTdxC)
   call spectC%mtimes_ik2_oop(ThatC, cbuffyC(:,:,:,1))
   call spectC%ifft(cbuffyC(:,:,:,1),dTdyC)
   call transpose_y_to_z(ThatC, cbuffzC(:,:,:,1), sp_gp)
   call der%ddz_C2E(cbuffzC(:,:,:,1),cbuffzE(:,:,:,1),1,1)
   call transpose_z_to_y(cbuffzE(:,:,:,1),cbuffyE(:,:,:,1),sp_gpE)
   call spectE%ifft(cbuffyE(:,:,:,1),dTdzE)
   call der%interpz_C2E(cbuffzC(:,:,:,1),cbuffzE(:,:,:,1),1,1)
   call der%ddz_E2C(cbuffzE(:,:,:,1),cbuffzC(:,:,:,2),1,1)
   call transpose_z_to_y(cbuffzC(:,:,:,2),cbuffyC(:,:,:,1),sp_gp)
   call spectC%ifft(cbuffyC(:,:,:,1),dTdzC)
   call message(0, "Done computing all the gradients needed in the SGS model.")
 
   ! STEP 6: Compute the SGS rhs terms 
   call sgsmodel%link_pointers(nuSGS, tauSGS_ij, tau13, tau23, q1, q2, q3, kappaSGS)
   call sgsmodel%getRHS_SGS(urhs, vrhs, wrhs, duidxjC, duidxjE, uhatC, vhatC, whatC, ThatC, u, v, wC, .true., dTdxC, dTdyC, dTdzC, dTdxE, dTdyE, dTdzE)
   call sgsmodel%getRHS_SGS_Scalar(Trhs, dTdxC, dTdyC, dTdzC, dTdzE, u, v, wC, T, ThatC, duidxjC)
! subroutine getRHS_SGS_Scalar(this, Trhs, dTdxC, dTdyC, dTdzC, dTdzE, u, v, w, T, That, duidxjC, TurbPrandtlNum, Cy, lowbound, highbound)

   call decomp_2d_write_one(1,nuSGS,'/work/04076/tg833754/stampede2/me461/testing_data/nuSGS_verify.out',gp)
   call decomp_2d_write_one(1,kappaSGS,'/work/04076/tg833754/stampede2/me461/testing_data/kappaSGS_verify.out',gp)
   call message(0,"Data successfully written to disk")

   call sgsmodel%destroy()
   deallocate(der)
   deallocate(sgsmodel)
   !deallocate(uhat, vhat, what)
   !deallocate(u, v, w, wChat, wChatz, whatz)
   
   call MPI_Finalize(ierr)

end program 
