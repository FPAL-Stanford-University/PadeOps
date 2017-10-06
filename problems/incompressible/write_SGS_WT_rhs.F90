program write_SGS_wt_rhs
   use mpi
   use constants
   use kind_parameters,  only: rkind, clen
   use timer, only: tic, toc
   use sgsmod_igrid, only: sgs_igrid
   use sgsmod, only: sgs
   use PadeDerOps, only: Pade6stagg
   use spectralMod, only: spectral
   use decomp_2d 
   use decomp_2d_io
   use turbineMod, only: turbineArray
   implicit none

   complex(rkind), dimension(:,:,:), allocatable :: uhatC, vhatC, whatE, uhatE, vhatE, whatC, ThatC,u_rhs,v_rhs,w_rhs
   real(rkind), dimension(:,:,:), allocatable :: uC, vC, wC, uE, vE, wE, fbody_x, fbody_y, fbody_z
   real(rkind), dimension(:,:,:,:), allocatable, target :: duidxjE, duidxjC,fbody,rbuffxC,rbuffyC,rbuffzC,rbuffyE,rbuffzE, duidxjE2
   complex(rkind), dimension(:,:,:,:), allocatable, target :: duidxjEhat,duidxjChat,cbuffyC,cbuffzC,cbuffyE,cbuffzE
   type(sgs_igrid) :: newsgs
   type(sgs) :: sgsold
   type(turbineArray) :: turbArray

   complex(rkind), parameter :: zeroC = dcmplx(0.0D0, 0.0D0)
   real(rkind), parameter :: Re = 1.d10, Fr = 1.d10, Pr = 1.d0
   real(rkind), parameter :: Tsurf = 1.d0, ThetaRef = 1.d0
   real(rkind) :: dx, dy, dz, Lx, Ly, Lz
   real(rkind), dimension(:,:,:,:), allocatable :: mesh
   real(rkind), dimension(:,:,:), allocatable ::  zMeshE, filteredSpeedSq
   type(spectral), target  :: spectE, spectC
   type(decomp_info) :: gpC, gpE
   type(decomp_info), pointer :: sp_gpC, sp_gpE
   type(Pade6stagg) :: Pade6opZ
   integer :: ierr, ix1, ixn, iy1, iyn, iz1, izn, TID, RID
   logical :: computeFbody
   
    integer :: wBC_bottom     = -1, wBC_top     = -1
    integer :: uBC_bottom     =  0, uBC_top     =  1
    integer :: vBC_bottom     =  0, vBC_top     =  1
    !integer :: TBC_bottom     =  1, TBC_top     =  0
    !integer :: WdUdzBC_bottom = -1, WdUdzBC_top =  1
    !integer :: WdVdzBC_bottom = -1, WdVdzBC_top =  1
    !integer :: WdWdzBC_bottom = -1, WdWdzBC_top = -1
    !integer :: WWBC_bottom    =  1, WWBC_top    =  1
    !integer :: UWBC_bottom    = -1, UWBC_top    = -1
    !integer :: VWBC_bottom    = -1, VWBC_top    = -1
    !integer :: WTBC_bottom    = -1, WTBC_top    = -1
    integer :: dUdzBC_bottom  =  0, dUdzBC_top  =  -1
    integer :: dVdzBC_bottom  =  0, dVdzBC_top  =  -1
    integer :: dWdzBC_bottom  =  1, dWdzBC_top  =  1
    !integer :: dTdzBC_bottom  =  -1, dTdzBC_top  =  -1

    integer :: scheme = 1
   character(len=clen) :: inputdir, outputdir, inputFile 
   integer :: nx, ny, nz, ioUnit, i, j, k
   real(rkind) :: z0init, dt, inst_horz_avg_turb(8), tsim, max_nuSGS
   namelist /INPUT/ Lx, Ly, Lz, z0init, outputdir, inputdir, nx, ny, nz, TID, RID


   ! Do MPI stuff
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)         

   ! Do file IO - input file
   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
   read(unit=ioUnit, NML=INPUT)
   close(ioUnit)    

   print*, "1", rid, tid
   
   call decomp_2d_init(nx, ny, nz, 0, 0)
   call get_decomp_info(gpC)

   call decomp_info_init(nx, ny, nz+1, gpE)

   ! Initialize spectral
   dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
    call spectC%init("x", nx, ny, nz, dx, dy, dz, "four", '2/3rd', 2 , .false.)
    call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", '2/3rd', 2 , .false.)

   sp_gpC => spectC%spectdecomp
   sp_gpE => spectE%spectdecomp

   ! Allocate memory
   allocate( mesh(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3),3) )
   allocate( zMeshE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)) )
   allocate( fbody_x(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)) )
   allocate( fbody_y(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)) )
   allocate( fbody_z(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)) )

   allocate( duidxjC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3),9) )
   allocate( duidxjE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3),9) )
   allocate( duidxjEhat(sp_gpE%ysz(1),sp_gpE%ysz(2),sp_gpE%ysz(3),9) )
   allocate( duidxjChat(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3),9) )
   allocate( uhatE(sp_gpE%ysz(1),sp_gpE%ysz(2),sp_gpE%ysz(3)) )
   allocate( vhatE(sp_gpE%ysz(1),sp_gpE%ysz(2),sp_gpE%ysz(3)) )
   allocate( whatE(sp_gpE%ysz(1),sp_gpE%ysz(2),sp_gpE%ysz(3)) )
   allocate( uhatC(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3)) )
   allocate( vhatC(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3)) )
   allocate( whatC(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3)) )
   allocate( ThatC(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3)) )
   allocate( uC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)) )
   allocate( vC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)) )
   allocate( wC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)) )
   allocate( uE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)) )
   allocate( vE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)) )
   allocate( wE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)) )
   allocate( u_rhs(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3))) ! -- what should this size be?
   allocate( v_rhs(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3))) ! -- what should this size be?
   allocate( w_rhs(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3))) ! -- what should this size be?
   allocate(filteredSpeedSq(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
   allocate( duidxjE2(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3),4) )

   ! Allocate Buffers
   allocate(cbuffzE(sp_gpE%zsz(1),sp_gpE%zsz(2),sp_gpE%zsz(3),2)) ! -- what should this size be?
   allocate(cbuffyE(sp_gpE%ysz(1),sp_gpE%ysz(2),sp_gpE%ysz(3),2)) ! -- what should this size be?
   allocate(cbuffyC(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3),2)) ! -- what should this size be?
   allocate(cbuffzC(sp_gpC%zsz(1),sp_gpC%zsz(2),sp_gpC%zsz(3),2)) ! -- what should this size be?
   allocate(rbuffxC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3),2)) ! -- what should this size be?
   allocate(rbuffyC(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3),2)) ! -- what should this size be?
   allocate(rbuffzC(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3),2)) ! -- what should this size be?
   allocate(rbuffyE(gpE%ysz(1),gpE%ysz(2),gpE%ysz(3),2)) ! -- what should this size be?
   allocate(rbuffzE(gpE%zsz(1),gpE%zsz(2),gpE%zsz(3),2)) ! -- what should this size be?

   ! Create Mesh
   ix1 = gpC%xst(1); iy1 = gpC%xst(2); iz1 = gpC%xst(3)
   ixn = gpC%xen(1); iyn = gpC%xen(2); izn = gpC%xen(3)
   do k=1,size(mesh,3)
       do j=1,size(mesh,2)
           do i=1,size(mesh,1)
               mesh(i,j,k,1) = real( ix1 + i - 1, rkind ) * dx
               mesh(i,j,k,2) = real( iy1 + j - 1, rkind ) * dy
               mesh(i,j,k,3) = real( iz1 + k - 1, rkind ) * dz + dz/two
           end do
       end do
   end do
   mesh(:,:,:,1) = mesh(:,:,:,1) - dx; mesh(:,:,:,2) = mesh(:,:,:,2) - dy; mesh(:,:,:,3) = mesh(:,:,:,3) - dz 

   iz1 = gpE%xst(3);  izn = gpE%xen(3)
   do k=1,size(zMeshE,3)
       do j=1,size(zMeshE,2)
           do i=1,size(zMeshE,1)
               zMeshE(i,j,k) = real( iz1 + k - 1, rkind ) * dz 
           end do
       end do
   end do
print *, 'dxdydz = ', (dx*dy*dz)**(2.0d0/3.0d0)

   ! Read in the fields (from restart files)
   computeFbody = .true.
   call readRestartFile(TID, RID)
   call spectC%fft(uC, uhatC)
   call spectC%fft(vC, vhatC)
   call spectE%fft(wE, whatE)

   ! Initialize Padeder
   call Pade6opz%init(gpC, sp_gpC, gpE, sp_gpE, dz, scheme,.false.)

   ! Initialize sgs
   call sgsold%init(1, spectC, spectE, gpC, gpE, dx, dy, dz, .false., .false., mesh(:,:,:,3), z0init,  .true., 1, .false., 0.7d0, .false., 1.0D0, 1.0d0, .true., .false.)
   call newsgs%init(gpC, gpE, spectC, spectE, dx, dy, dz, inputfile, zMeshE(1,1,:), mesh(1,1,:,3), fbody_x, fbody_y, fbody_z, &
      computeFbody, Pade6opZ, cbuffyC, cbuffzC, cbuffyE, cbuffzE, rbuffxC, rbuffyC, rbuffzC, rbuffyE, rbuffzE, Tsurf, ThetaRef, Fr, Re, & 
      Pr, .false., .false.,1)

   ! Initialize WT
   call turbArray%init(inputFile, gpC, gpE, spectC, spectE, cbuffyC, cbuffyE, cbuffzC, cbuffzE, mesh, dx, dy, dz) 
   
   ! Interpolations
   call interp_primitivevars()
   
   ! duidxj
   call compute_duidxj()
   print *, maxval(abs(duidxjC(:,:,:,1)+duidxjC(:,:,:,5) + duidxjC(:,:,:,9)))

   ! Get RHS WT
   u_rhs = zeroC; v_rhs = zeroC; w_rhs = zeroC
   call turbArray%getForceRHS(dt, uC, vC, wC, u_rhs, v_rhs, w_rhs, .true., inst_horz_avg_turb)
   
   call spectC%ifft(u_rhs,fbody_x)
   call spectC%ifft(v_rhs,fbody_y)
   call spectE%ifft(w_rhs,fbody_z)


   call dumpFullField(fbody_x,"xtrb")
   call dumpFullField(fbody_y,"ytrb")

   ! Get RHS constant force
   fbody_x = fbody_x + 1.d0 


   u_rhs = zeroC; v_rhs = zeroC; w_rhs = zeroC
  

   duidxjE2(:,:,:,1) = duidxjE(:,:,:,7) 
   duidxjE2(:,:,:,2) = duidxjE(:,:,:,8) 
   duidxjE2(:,:,:,3) = duidxjE(:,:,:,3) 
   duidxjE2(:,:,:,4) = duidxjE(:,:,:,6) 
   call sgsold%getRHS_SGS_WallM(duidxjC, duidxjE2 , duidxjChat ,& 
                                   u_rhs  , v_rhs          , w_rhs      ,&
                                   uhatC   , vhatC           , whatC      ,&
                                   uC      , vC              , wC         ,&
                                   0.1d0  , 1.0d0            , 0.0d0        ,&
                                   1.0d0  , filteredSpeedSq, 0.0d0,&
                                   max_nuSGS, inst_horz_avg_turb)
   
   call spectC%ifft(u_rhs,fbody_x)
   call spectC%ifft(v_rhs,fbody_y)
   call spectE%ifft(w_rhs,fbody_z)
  
   call dumpFullField(fbody_x,"xsgs")
   call dumpFullField(fbody_y,"ysgs")
   !call writeVisualizationFile(tid, rid, fbody_z, gpE, "zsgs", OutputDir)
   
   ! get tau_sgs
   !call newsgs%getTauSGS(duidxjC, duidxjE, duidxjEhat, uhatE, vhatE, whatE, uhatC, vhatC, ThatC, uC, vC, uE, vE, wE, newTimeStep = .true.)
  





   ! deallocate everything
   !deallocate(wE)
   !deallocate(vE)
   !deallocate(uE)
   !deallocate(wC)
   !deallocate(vC)
   !deallocate(uC)
   !deallocate(ThatC)
   !deallocate(whatC)
   !deallocate(vhatC)
   !deallocate(uhatC)
   !deallocate(whatE)
   !deallocate(vhatE)
   !deallocate(uhatE)
   !deallocate(duidxjEhat)
   !deallocate(duidxjChat)
   !deallocate(duidxjE)
   !deallocate(duidxjC)
   !deallocate(rbuffzE)
   !deallocate(rbuffyE)
   !deallocate(rbuffzC)
   !deallocate(rbuffyC)
   !deallocate(rbuffxC)
   !deallocate(cbuffzC)
   !deallocate(cbuffyC)
   !deallocate(cbuffyE)
   !deallocate(cbuffzE)
   !deallocate(fbody)
   !deallocate(zMeshE)
   !deallocate(mesh)
   !deallocate(w_rhs)
   !deallocate(v_rhs)
   !deallocate(u_rhs)







contains

    subroutine compute_duidxj()
        complex(rkind), dimension(:,:,:), pointer :: ctmpz1, ctmpz2
        complex(rkind), dimension(:,:,:), pointer :: ctmpz3, ctmpz4
        complex(rkind), dimension(:,:,:), pointer :: ctmpy1, ctmpy2
        real(rkind),    dimension(:,:,:), pointer :: dudx, dudy, dudz
        real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind),    dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC 
        real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC
        real(rkind),    dimension(:,:,:), pointer :: dwdzE, dudxE
        real(rkind),    dimension(:,:,:), pointer :: dudyE, dvdxE 
        real(rkind),    dimension(:,:,:), pointer :: dvdyE
        
        complex(rkind), dimension(:,:,:), pointer :: dudxH, dudyH, dudzH 
        complex(rkind), dimension(:,:,:), pointer :: dvdxH, dvdyH, dvdzH
        complex(rkind), dimension(:,:,:), pointer :: dwdxH, dwdyH, dwdzH
        
        complex(rkind), dimension(:,:,:), pointer :: dudxEH, dudyEH, dudzEH 
        complex(rkind), dimension(:,:,:), pointer :: dvdxEH, dvdyEH, dvdzEH
        complex(rkind), dimension(:,:,:), pointer :: dwdxEH, dwdyEH, dwdzEH

        dudx  => duidxjC(:,:,:,1); dudy  => duidxjC(:,:,:,2); dudzC => duidxjC(:,:,:,3) 
        dvdx  => duidxjC(:,:,:,4); dvdy  => duidxjC(:,:,:,5); dvdzC => duidxjC(:,:,:,6) 
        dwdxC => duidxjC(:,:,:,7); dwdyC => duidxjC(:,:,:,8); dwdz  => duidxjC(:,:,:,9) 

        dudxH => duidxjChat(:,:,:,1); dudyH => duidxjChat(:,:,:,2); dudzH => duidxjChat(:,:,:,3) 
        dvdxH => duidxjChat(:,:,:,4); dvdyH => duidxjChat(:,:,:,5); dvdzH => duidxjChat(:,:,:,6) 
        dwdxH => duidxjChat(:,:,:,7); dwdyH => duidxjChat(:,:,:,8); dwdzH => duidxjChat(:,:,:,9) 
        
        dudxEH => duidxjEhat(:,:,:,1); dudyEH => duidxjEhat(:,:,:,2); dudzEH => duidxjEhat(:,:,:,3) 
        dvdxEH => duidxjEhat(:,:,:,4); dvdyEH => duidxjEhat(:,:,:,5); dvdzEH => duidxjEhat(:,:,:,6) 
        dwdxEH => duidxjEhat(:,:,:,7); dwdyEH => duidxjEhat(:,:,:,8); dwdzEH => duidxjEhat(:,:,:,9)
       
        dudxE => duidxjE(:,:,:,1); dudyE => duidxjE(:,:,:,2); dudz  => duidxjE(:,:,:,3)
        dvdxE => duidxjE(:,:,:,4); dvdyE => duidxjE(:,:,:,5); dvdz  => duidxjE(:,:,:,6)
        dwdx  => duidxjE(:,:,:,7); dwdy  => duidxjE(:,:,:,8); dwdzE => duidxjE(:,:,:,9)

        ctmpz1 => cbuffzC(:,:,:,1); ctmpz2 => cbuffzE(:,:,:,1); 
        ctmpz3 => cbuffzC(:,:,:,2); ctmpz4 => cbuffzE(:,:,:,2)
        ctmpy1 => cbuffyC(:,:,:,1); ctmpy2 => cbuffyE(:,:,:,1)

      
        ! dudx
        call spectC%mTimes_ik1_oop(uhatC,dudxH)
        call spectC%ifft(dudxH,dudx)
        call spectE%mTimes_ik1_oop(uhatE,dudxEH)
        call spectE%ifft(dudxEH,dudxE)
         
        ! dudy
        call spectC%mTimes_ik2_oop(uhatC,dudyH)
        call spectC%ifft(dudyH,dudy)
        call spectE%mTimes_ik2_oop(uhatE,dudyEH)
        call spectE%ifft(dudyEH,dudyE)

        ! dvdx 
        call spectC%mTimes_ik1_oop(vhatC,dvdxH)
        call spectC%ifft(dvdxH,dvdx)
        call spectE%mTimes_ik1_oop(vhatE,dvdxEH)
        call spectE%ifft(dvdxEH,dvdxE)

        ! dvdy
        call spectC%mTimes_ik2_oop(vhatC,dvdyH)
        call spectC%ifft(dvdyH,dvdy)
        call spectE%mTimes_ik2_oop(vhatE,dvdyEH)
        call spectE%ifft(dvdyEH,dvdyE)

        ! dwdx
        call spectC%mTimes_ik1_oop(whatC,dwdxH)
        call spectC%ifft(dwdxH,dwdxC)
        call spectE%mTimes_ik1_oop(whatE, dwdxEH)
        call spectE%ifft(dwdxEH,dwdx)

        ! dwdy
        call spectC%mTimes_ik2_oop(whatC,dwdyH)
        call spectC%ifft(dwdyH,dwdyC)
        call spectE%mTimes_ik2_oop(whatE,dwdyEH)
        call spectE%ifft(dwdyEH,dwdy)
       
        ! dwdz
        call transpose_y_to_z(whatE,ctmpz2,sp_gpE)
        call Pade6opZ%ddz_E2C(ctmpz2,ctmpz1,wBC_bottom,wBC_top)
        call transpose_z_to_y(ctmpz1,dwdzH,sp_gpC)
        call spectC%ifft(dwdzH,dwdz)
        call Pade6opZ%interpz_C2E(ctmpz1,ctmpz4,dWdzBC_bottom, dWdzBC_top)
        call transpose_z_to_y(ctmpz4,dwdzEH,sp_gpE)
        call spectE%ifft(dwdzEH,dwdzE)

        !! d2wdz2
        !if(.not. isinviscid) then
        !   call Pade6opZ%d2dz2_E2E(ctmpz2,ctmpz4,wBC_bottom,wBC_top)
        !   call transpose_z_to_y(ctmpz4,d2wdz2hatE,sp_gpE)
        !end if

        ! dudz and d2udz2
        call transpose_y_to_z(uhatC,ctmpz1,sp_gpC)
        call Pade6opZ%ddz_C2E(ctmpz1,ctmpz2,uBC_bottom,uBC_top)
        call transpose_z_to_y(ctmpz2,dudzEH,sp_gpE)
        call spectE%ifft(dudzEH,dudz)
        !if (.not. isinviscid) then
        !       if ((uBC_top == 0) .or. (uBC_bottom == 0)) then
        !          call Pade6opZ%ddz_C2E(ctmpz1,ctmpz4,uBC_bottom,uBC_top)
        !          call Pade6opZ%ddz_E2C(ctmpz4,ctmpz3,dUdzBC_bottom,dUdzBC_top)
        !       else
        !          call Pade6opZ%d2dz2_C2C(ctmpz1,ctmpz3,uBC_bottom,uBC_top)
        !       end if
        !       call transpose_z_to_y(ctmpz3,d2udz2hatC,sp_gpC)
        !end if 
        call Pade6opZ%interpz_E2C(ctmpz2,ctmpz1,dUdzBC_bottom,dUdzBC_top)
        call transpose_z_to_y(ctmpz1,dudzH,sp_gpC)
        call spectC%ifft(dudzH,dudzC)
      
        ! dvdz and d2vdz2
        call transpose_y_to_z(vhatC,ctmpz1,sp_gpC)
        call Pade6opZ%ddz_C2E(ctmpz1,ctmpz2,vBC_bottom,vBC_top)
        call transpose_z_to_y(ctmpz2,dvdzEH,sp_gpE)
        call spectE%ifft(dvdzEH,dvdz)
        !if (.not. isinviscid) then
        !    if ((vBC_top == 0) .or. (vBC_bottom == 0)) then
        !       call Pade6opZ%ddz_C2E(ctmpz1,ctmpz4,vBC_bottom,vBC_top)
        !       call Pade6opZ%ddz_E2C(ctmpz4,ctmpz3,dVdzBC_bottom,dVdzBC_top)
        !    else
        !       call Pade6opZ%d2dz2_C2C(ctmpz1,ctmpz3,vBC_bottom,vBC_top)
        !    end if
        !    call transpose_z_to_y(ctmpz3,d2vdz2hatC,sp_gpC)
        !end if 
        call Pade6opZ%interpz_E2C(ctmpz2,ctmpz1,dVdzBC_bottom,dVdzBC_top)
        call transpose_z_to_y(ctmpz1,dvdzH,sp_gpC)
        call spectC%ifft(dvdzH,dvdzC)

    end subroutine

    subroutine interp_PrimitiveVars()
        complex(rkind), dimension(:,:,:), pointer :: ybuffC, zbuffC, zbuffE
        
        zbuffE => cbuffzE(:,:,:,1)
        zbuffC => cbuffzC(:,:,:,1)
        ybuffC => cbuffyC(:,:,:,1)

        ! Step 1: Interpolate w -> wC
        call transpose_y_to_z(whatE,zbuffE,sp_gpE)
        call Pade6opZ%interpz_E2C(zbuffE,zbuffC,wBC_bottom, wBC_top)
        call transpose_z_to_y(zbuffC,whatC,sp_gpC)
        call spectC%ifft(whatC,wC)

        ! Step 2: Interpolate u -> uE
        call transpose_y_to_z(uhatC,zbuffC,sp_gpC)
        call Pade6opZ%interpz_C2E(zbuffC,zbuffE,uBC_bottom, uBC_top)
        call transpose_z_to_y(zbuffE,uhatE, sp_gpE)
        call spectE%ifft(uhatE, uE)

        ! Step 3: Interpolate v -> vE
        call transpose_y_to_z(vhatC,zbuffC,sp_gpC)
        call Pade6opZ%interpz_C2E(zbuffC,zbuffE,vBC_bottom, vBC_top)
        call transpose_z_to_y(zbuffE,vhatE, sp_gpE)
        call spectE%ifft(vhatE, vE)
        

        !! Step 4: Interpolate T
        !if (isStratified) then
        !    call transpose_y_to_z(That,zbuffC,sp_gpC)
        !    call Pade6opZ%interpz_C2E(zbuffC,zbuffE,TBC_bottom, TBC_top)
        !    if (botBC_Temp == 0) then 
        !        zbuffE(:,:,1) = zero 
        !        if (nrank == 0) then
        !            zbuffE(1,1,1) = Tsurf*real(nx,rkind)*real(ny,rkind)
        !        end if 
        !    end if
        !    call transpose_z_to_y(zbuffE,TEhat,sp_gpE)
        !    call spectE%ifft(TEhat,TE)
        !end if 
    end subroutine

    subroutine dumpFullField(arr,label)
        use decomp_2d_io
        use mpi
        character(len=clen) :: tempname, fname
        real(rkind), dimension(:,:,:), intent(in) :: arr
        character(len=4), intent(in) :: label

        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_",label,"_t",tid,".out"
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,arr,fname, gpC)

    end subroutine

    subroutine readRestartFile(tid, rid)
        use decomp_2d_io
        use mpi
        use exits, only: message
        integer, intent(in) :: tid, rid
        character(len=clen) :: tempname, fname
        integer :: ierr

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_u.",tid
        fname = trim(InputDir)//"/"//trim(tempname)
        call decomp_2d_read_one(1,uC,fname, gpC)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_v.",tid
        fname = trim(InputDir)//"/"//trim(tempname)
        call decomp_2d_read_one(1,vC,fname, gpC)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_w.",tid
        fname = trim(InputDir)//"/"//trim(tempname)
        call decomp_2d_read_one(1,wE,fname, gpE)

        !if (this%isStratified) then
        !    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_T.",tid
        !    fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
        !    call decomp_2d_read_one(1,this%T,fname, this%gpC)
        !end if 

        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",rid, "_info.",tid
        fname = trim(InputDir)//"/"//trim(tempname)
        open(unit=10,file=fname,access='sequential',form='formatted')
        read (10, *)  tsim
        close(10)

        call mpi_barrier(mpi_comm_world, ierr)
        call message("================= RESTART FILE USED ======================")
        call message(0, "Simulation Time at restart:", tsim)
        call message("=================================== ======================")

    end subroutine

end program 

