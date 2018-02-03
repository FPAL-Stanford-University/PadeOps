program getEnergyBudgetTerms
   use mpi
   use constants
   use kind_parameters,  only: rkind, clen
   use timer, only: tic, toc
   use sgsmod_igrid, only: sgs_igrid
   use PadeDerOps, only: Pade6stagg
   use spectralMod, only: spectral
   use decomp_2d 
   use decomp_2d_io
   use turbineMod, only: turbineArray
   implicit none

   complex(rkind), dimension(:,:,:), allocatable :: uhatC, vhatC, whatE, uhatE, vhatE, whatC, ThatC,u_rhs,v_rhs,w_rhs
   real(rkind), dimension(:,:,:), allocatable :: pC, uC, vC, wC, uE, vE, wE, fbody_x, fbody_y, fbody_z, fbody_zC
   real(rkind), dimension(:,:,:,:), allocatable, target :: duidxjE, duidxjC, rbuffxC,rbuffyC,rbuffzC,rbuffyE,rbuffzE, duidxjE2
   complex(rkind), dimension(:,:,:,:), allocatable, target :: duidxjEhat,duidxjChat,cbuffyC,cbuffzC,cbuffyE,cbuffzE
   type(sgs_igrid) :: newsgs
   type(turbineArray), allocatable :: turbArray

   complex(rkind), parameter :: zeroC = dcmplx(0.0D0, 0.0D0)
   real(rkind), parameter :: Re = 1.d10, Fr = 1.d10, Pr = 1.d10
   real(rkind), parameter :: Tsurf = 1.d0, ThetaRef = 1.d0
   real(rkind) :: dx, dy, dz, Lx, Ly, Lz
   real(rkind), dimension(:,:,:,:), allocatable :: mesh
   real(rkind), dimension(:,:,:), allocatable :: zMeshE, filteredSpeedSq, fx_turb_store, fy_turb_store, fz_turb_store 
   real(rkind), dimension(:,:,:), allocatable :: dPdx_store, dPdy_store, dPdz_store, u_store, v_store, w_store, P_store, uu_store, uv_store, uw_store, vv_store, vw_store, ww_store
   real(rkind), dimension(:,:,:), allocatable :: u_bar, v_bar, w_bar, p_bar, uprime_uprime_bar, uprime_vprime_bar, uprime_wprime_bar, vprime_vprime_bar, vprime_wprime_bar, wprime_wprime_bar
   real(rkind), dimension(:,:,:), allocatable :: xtrbprime_uprime_bar, xsgsprime_uprime_bar, ysgsprime_vprime_bar, zsgsprime_wprime_bar
   real(rkind), dimension(:,:,:), allocatable :: K_bar, fx_turb_bar, fy_turb_bar, fz_turb_bar, fx_sgs_bar, fy_sgs_bar, fz_sgs_bar
   real(rkind), dimension(:,:,:), allocatable :: fx_sgs_store, fy_sgs_store, fz_sgs_store, xtrbu_store, xsgsu_store, ysgsv_store, zsgsw_store
   type(spectral), target  :: spectE, spectC
   type(decomp_info) :: gpC, gpE
   type(decomp_info), pointer :: sp_gpC, sp_gpE
   type(Pade6stagg) :: Pade6opZ
   integer :: ierr, ix1, ixn, iy1, iyn, iz1, izn, RID
   logical :: computeFbody 
  
   real(rkind), dimension(:,:,:)  , pointer :: nuSGS, kappaSGS
   real(rkind), dimension(:,:,:)  , pointer :: tau13, tau23
   real(rkind), dimension(:,:,:,:), pointer :: tauSGS_ij
   real(rkind), dimension(:,:,:)  , pointer :: q1, q2, q3
   integer :: wBC_bottom     = -1, wBC_top     = -1
   integer :: uBC_bottom     =  0, uBC_top     =  1
   integer :: vBC_bottom     =  0, vBC_top     =  1
   integer :: dUdzBC_bottom  =  0, dUdzBC_top  =  -1
   integer :: dVdzBC_bottom  =  0, dVdzBC_top  =  -1
   integer :: dWdzBC_bottom  =  1, dWdzBC_top  =  1

   character(len=clen) :: inputdir, outputdir, inputFile 
   integer :: nx, ny, nz, ioUnit, i, j, k, nvis = 0, tid_initial, tid_final, dtid, ind=0
   real(rkind) :: dt, inst_horz_avg_turb(8), tsim
   namelist /INPUT/ Lx, Ly, Lz, outputdir, inputdir, nx, ny, nz, RID, tid_initial, tid_final, dtid
   
   logical :: PeriodicInZ = .false.
   integer :: botWall, topWall, NumericalSchemeVert = 1
   namelist /BCs/ PeriodicInZ, botWall, topWall

   namelist /NUMERICS/ NumericalSchemeVert


   ! Do MPI stuff
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)         

   ! Do file IO - input file
   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
   read(unit=ioUnit, NML=INPUT)
   close(ioUnit)    

   call decomp_2d_init(nx, ny, nz, 0, 0)
   call get_decomp_info(gpC)

   call decomp_info_init(nx, ny, nz+1, gpE)

   ! Initialize spectral
   dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
    call spectC%init("x", nx, ny, nz, dx, dy, dz, "four", '2/3rd', 2 , .false.)
    call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", '2/3rd', 2 , .false.)

   sp_gpC => spectC%spectdecomp
   sp_gpE => spectE%spectdecomp

   call initializeEverything()

   ! DO LOOP START HERE
   do ind = tid_initial, tid_final, dtid 
        call tic()
        
        ! Initialize WT
        allocate(turbArray)
        call turbArray%init(inputFile, gpC, gpE, spectC, spectE, cbuffyC, cbuffyE, cbuffzC, cbuffzE, mesh, dx, dy, dz) 
        
        
        ! READ FIELDS)
        ! need to add read in pressure to readVisualizationFile subroutine
        call readVisualizationFile(ind, RID)
        print *, ind
        call spectC%fft(uC, uhatC)
        call spectC%fft(vC, vhatC)
        call spectE%fft(wE, whatE)
 
        ! PREPROCESS FIELDS
        call interp_primitivevars()
        call compute_duidxj()

        ! WIND TURBINE STUFF
        u_rhs = zeroC; v_rhs = zeroC; w_rhs = zeroC
        call turbArray%getForceRHS(dt, uC, vC, wC, u_rhs, v_rhs, w_rhs, .true., inst_horz_avg_turb)
        call spectC%ifft(u_rhs,fbody_x)
        call spectC%ifft(v_rhs,fbody_y)
        call spectE%ifft(w_rhs,fbody_z)
        
        call transpose_x_to_y(fbody_z,rbuffyE(:,:,:,1),gpE)
        call transpose_y_to_z(rbuffyE(:,:,:,1),rbuffzE(:,:,:,1),gpE)
        call Pade6opz%interpz_E2C(rbuffzE(:,:,:,1),rbuffzC(:,:,:,1),0,0)
        call transpose_z_to_y(rbuffzC(:,:,:,1),rbuffyC(:,:,:,1),gpC)
        call transpose_y_to_x(rbuffyC(:,:,:,1),fbody_zC,gpC)

        fx_turb_store = fx_turb_store + fbody_x 
        fy_turb_store = fy_turb_store + fbody_y
        fz_turb_store = fz_turb_store + fbody_zC 

        xtrbu_store = xtrbu_store + fbody_x*uC       
 
        ! SGS MODEL STUFF
        u_rhs = zeroC; v_rhs = zeroC; w_rhs = zeroC
        !call newsgs%getRHS_SGS(u_rhs, v_rhs, w_rhs, duidxjC, duidxjE, duidxjEhat, uhatE, vhatE, whatE, uhatC, vhatC, ThatC, uC, vC, uE, vE, wE, .true.)
        call newsgs%getRHS_SGS(u_rhs, v_rhs, w_rhs, duidxjC, duidxjE,  uhatC, vhatC, whatC, ThatC, uC, vC, wC, .true.)
        call spectC%ifft(u_rhs,fbody_x)
        call spectC%ifft(v_rhs,fbody_y)
        call spectE%ifft(w_rhs,fbody_z)
        
        call transpose_x_to_y(fbody_z,rbuffyE(:,:,:,1),gpE)
        call transpose_y_to_z(rbuffyE(:,:,:,1),rbuffzE(:,:,:,1),gpE)
        call Pade6opz%interpz_E2C(rbuffzE(:,:,:,1),rbuffzC(:,:,:,1),0,0)
        call transpose_z_to_y(rbuffzC(:,:,:,1),rbuffyC(:,:,:,1),gpC)
        call transpose_y_to_x(rbuffyC(:,:,:,1),fbody_zC,gpC)

        fx_sgs_store = fx_sgs_store + fbody_x 
        fy_sgs_store = fy_sgs_store + fbody_y
        fz_sgs_store = fz_sgs_store + fbody_zC 

        xsgsu_store = xsgsu_store + fbody_x*uC        
        ysgsv_store = ysgsv_store + fbody_y*vC        
        zsgsw_store = zsgsw_store + fbody_zC*wC        

        ! GRADIENT OF  ADVECTION TERM
        ! need to compute avg(ui)avg(uj)
        ! then take derivative d/dx
        u_store = u_store +uC
        v_store = v_store +vC
        w_store = w_store +wC

        ! GRADIENT OF REYNOLDS STRESS
        ! need to compute avg(ui'uj')
        ! then take derivative d/dx
        uu_store = uu_store +(uC*uC)
        uv_store = uv_store +uC*vC
        uw_store = uw_store +uC*wC
        vv_store = vv_store +vC*vC
        vw_store = vw_store +vC*wC
        ww_store = ww_store +wC*wC

        ! GRADIENT OF PRESSURE 
        ! derivative of P: dP/dx
        P_store = P_store + pC

        ! WRAP UP 
        deallocate(turbArray)
        call toc()

        nvis = nvis + 1
        
        ! END DO LOOP
   end do
   
   u_bar = u_store/real(nvis,rkind)
   v_bar = v_store/real(nvis,rkind)
   w_bar = w_store/real(nvis,rkind)
   p_bar = P_store/real(nvis,rkind)
   K_bar = 0.5*(u_bar**2 +v_bar**2 +w_bar**2)
   fx_turb_bar = fx_turb_store/real(nvis,rkind)
   fy_turb_bar = fy_turb_store/real(nvis,rkind)
   fz_turb_bar = fz_turb_store/real(nvis,rkind)
   fx_sgs_bar  = fx_sgs_store/real(nvis,rkind)
   fy_sgs_bar  = fy_sgs_store/real(nvis,rkind)
   fz_sgs_bar  = fz_sgs_store/real(nvis,rkind)

        ! Turbine Force
   call dumpFullField(fx_turb_bar,"xtrb")
   call dumpFullField(fy_turb_bar,"ytrb")
   call dumpFullField(fz_turb_bar,"ztrb")
   
        ! SGS Force
   call dumpFullField(fx_sgs_bar,"xsgs")
   call dumpFullField(fy_sgs_bar,"ysgs")
   call dumpFullField(fz_sgs_bar,"zsgs")
        
        !Advection of Mean Kinetic Energy
   call getddx(K_bar*u_bar,rbuffxC(:,:,:,1))
   call getddy(K_bar*v_bar,rbuffxC(:,:,:,2))
   call getddz(K_bar*w_bar,rbuffxC(:,:,:,3))
   call dumpFullField(rbuffxC(:,:,:,1)+rbuffxC(:,:,:,2)+rbuffxC(:,:,:,3),"mken")

        !Pressure Gradient
   call getddx(p_bar,rbuffxC(:,:,:,1))
   call dumpFullField(rbuffxC(:,:,:,1),"dPdx")
   call getddy(p_bar,rbuffxC(:,:,:,1))
   call dumpFullField(rbuffxC(:,:,:,1),"dPdy")
   call getddz(p_bar,rbuffxC(:,:,:,1))
   call dumpFullField(rbuffxC(:,:,:,1),"dPdz")

     !Mean momentum gradient
   call getddx(u_bar*u_bar,rbuffxC(:,:,:,1))
   call getddy(u_bar*v_bar,rbuffxC(:,:,:,2))
   call getddz(u_bar*w_bar,rbuffxC(:,:,:,3))
   call dumpFullField((rbuffxC(:,:,:,1)+rbuffxC(:,:,:,2)+rbuffxC(:,:,:,3)),"duux")
   
   call getddx(v_bar*u_bar,rbuffxC(:,:,:,1))
   call getddy(v_bar*v_bar,rbuffxC(:,:,:,2))
   call getddz(v_bar*w_bar,rbuffxC(:,:,:,3))
   call dumpFullField((rbuffxC(:,:,:,1)+rbuffxC(:,:,:,2)+rbuffxC(:,:,:,3)),"duuy")
   
   call getddx(w_bar*u_bar,rbuffxC(:,:,:,1))
   call getddy(w_bar*v_bar,rbuffxC(:,:,:,2))
   call getddz(w_bar*w_bar,rbuffxC(:,:,:,3))
   call dumpFullField((rbuffxC(:,:,:,1)+rbuffxC(:,:,:,2)+rbuffxC(:,:,:,3)),"duuz")

   uprime_uprime_bar = uu_store/real(nvis,rkind)-u_bar*u_bar
   uprime_vprime_bar = uv_store/real(nvis,rkind)-u_bar*v_bar
   uprime_wprime_bar = uw_store/real(nvis,rkind)-u_bar*w_bar
   vprime_vprime_bar = vv_store/real(nvis,rkind)-v_bar*v_bar
   vprime_wprime_bar = vw_store/real(nvis,rkind)-v_bar*w_bar
   wprime_wprime_bar = ww_store/real(nvis,rkind)-w_bar*w_bar

      !Body Force Fluctuation terms
   xtrbprime_uprime_bar = xtrbu_store/real(nvis,rkind)-fx_turb_bar*u_bar
   xsgsprime_uprime_bar = xsgsu_store/real(nvis,rkind)-fx_sgs_bar*u_bar
   ysgsprime_vprime_bar = ysgsv_store/real(nvis,rkind)-fy_sgs_bar*v_bar
   zsgsprime_wprime_bar = zsgsw_store/real(nvis,rkind)-fz_sgs_bar*w_bar

   call dumpFullField(xtrbprime_uprime_bar,"xtup")
   call dumpFullField(xsgsprime_uprime_bar,"xsup")
   call dumpFullField(ysgsprime_vprime_bar,"ysvp")
   call dumpFullField(zsgsprime_wprime_bar,"zswp")

        !Advection of Turbulent Kinetic Energy (Turbulent Transport))
   call getddx(u_bar*uprime_uprime_bar+v_bar*uprime_vprime_bar+w_bar*uprime_wprime_bar,rbuffxC(:,:,:,1))
   call getddy(u_bar*uprime_vprime_bar+v_bar*vprime_vprime_bar+w_bar*vprime_wprime_bar,rbuffxC(:,:,:,2))
   call getddz(u_bar*uprime_wprime_bar+v_bar*vprime_wprime_bar+w_bar*wprime_wprime_bar,rbuffxC(:,:,:,3))
   call dumpFullField(rbuffxC(:,:,:,1)+rbuffxC(:,:,:,2)+rbuffxC(:,:,:,3),"tken")

        !Body Forces
   call dumpFullField(u_bar*fx_turb_bar+v_bar*fy_turb_bar+w_bar*fz_turb_bar,"etrb")
   call dumpFullField(u_bar*fx_sgs_bar +v_bar*fy_sgs_bar +w_bar*fz_sgs_bar ,"esgs")
   call dumpFullField(u_bar*1          +v_bar*1          +w_bar*1          ,"egrd")

        ! Pressure Diffusion
   call getddx(u_bar*P_bar,rbuffxC(:,:,:,1))
   call getddy(v_bar*P_bar,rbuffxC(:,:,:,2))
   call getddz(w_bar*P_bar,rbuffxC(:,:,:,3))
   call dumpFullField((rbuffxC(:,:,:,1)+rbuffxC(:,:,:,2)+rbuffxC(:,:,:,3)),"pdif")
        
        !Production
   call getddx(u_bar,rbuffxC(:,:,:,1))
   call getddy(u_bar,rbuffxC(:,:,:,2))
   call getddz(u_bar,rbuffxC(:,:,:,3))
   rbuffxC(:,:,:,4) = uprime_uprime_bar*rbuffxC(:,:,:,1)+uprime_vprime_bar*rbuffxC(:,:,:,2)+uprime_wprime_bar*rbuffxC(:,:,:,3)

   call getddx(v_bar,rbuffxC(:,:,:,1))
   call getddy(v_bar,rbuffxC(:,:,:,2))
   call getddz(v_bar,rbuffxC(:,:,:,3))
   rbuffxC(:,:,:,4) = rbuffxC(:,:,:,4)+uprime_vprime_bar*rbuffxC(:,:,:,1)+vprime_vprime_bar*rbuffxC(:,:,:,2)+vprime_wprime_bar*rbuffxC(:,:,:,3)

   call getddx(w_bar,rbuffxC(:,:,:,1))
   call getddy(w_bar,rbuffxC(:,:,:,2))
   call getddz(w_bar,rbuffxC(:,:,:,3))
   rbuffxC(:,:,:,4) = rbuffxC(:,:,:,4)+uprime_wprime_bar*rbuffxC(:,:,:,1)+vprime_wprime_bar*rbuffxC(:,:,:,2)+wprime_wprime_bar*rbuffxC(:,:,:,3)

   call dumpFullField(rbuffxC(:,:,:,4),"prod")
        
        !Reynolds stress gradient
   call getddx(uprime_uprime_bar,rbuffxC(:,:,:,1))
   call getddy(uprime_vprime_bar,rbuffxC(:,:,:,2))
   call getddz(uprime_wprime_bar,rbuffxC(:,:,:,3))
   call dumpFullField((rbuffxC(:,:,:,1)+rbuffxC(:,:,:,2)+rbuffxC(:,:,:,3)),"reyx")

   call getddx(uprime_vprime_bar,rbuffxC(:,:,:,1))
   call getddy(vprime_vprime_bar,rbuffxC(:,:,:,2))
   call getddz(vprime_wprime_bar,rbuffxC(:,:,:,3))
   call dumpFullField((rbuffxC(:,:,:,1)+rbuffxC(:,:,:,2)+rbuffxC(:,:,:,3)),"reyy")

   call getddx(uprime_wprime_bar,rbuffxC(:,:,:,1))
   call getddy(vprime_wprime_bar,rbuffxC(:,:,:,2))
   call getddz(wprime_wprime_bar,rbuffxC(:,:,:,3))
   call dumpFullField((rbuffxC(:,:,:,1)+rbuffxC(:,:,:,2)+rbuffxC(:,:,:,3)),"reyz")

        !Time averaged velocity and pressure
   call dumpFullField(u_bar,"uVel")
   call dumpFullField(v_bar,"vVel")
   call dumpFullField(w_bar,"wVel")
   call dumpFullField(p_bar,"prss")

        ! Time-averaged Reynolds Stress fields
   call dumpFullField(uprime_uprime_bar,"upup")
   call dumpFullField(uprime_vprime_bar,"upvp")
   call dumpFullField(uprime_wprime_bar,"upwp")
   call dumpFullField(vprime_vprime_bar,"vpvp")
   call dumpFullField(vprime_wprime_bar,"vpwp")
   call dumpFullField(wprime_wprime_bar,"wpwp")

        ! Divergence
   call getddx(u_bar,rbuffxC(:,:,:,1))
   call getddy(v_bar,rbuffxC(:,:,:,2))
   call getddz(w_bar,rbuffxC(:,:,:,3))
   call dumpFullField((rbuffxC(:,:,:,1)+rbuffxC(:,:,:,2)+rbuffxC(:,:,:,3)),"divr")
   
   
   call mpi_barrier(mpi_comm_world, ierr)
   !stop
   call MPI_Finalize(ierr)           !<-- Terminate MPI 


contains
   subroutine initializeEverything()

      ! Allocate memory
      allocate( mesh(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3),3) )
      allocate( zMeshE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)) )
      allocate( fbody_x(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)) )
      allocate( fbody_y(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)) )
      allocate( fbody_zC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)) )
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
      allocate( pC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)) )
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

      allocate(fx_turb_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(fy_turb_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(fz_turb_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(fx_turb_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(fy_turb_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(fz_turb_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      
      allocate(fx_sgs_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(fy_sgs_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(fz_sgs_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(fx_sgs_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(fy_sgs_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(fz_sgs_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))

      allocate(dPdx_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(dPdy_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(dPdz_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(P_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(u_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(v_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(w_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))

      allocate(uu_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(uv_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(uw_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(vv_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(vw_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(ww_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))

      allocate(u_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(v_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(w_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(p_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(K_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))

      allocate(uprime_uprime_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(uprime_vprime_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(uprime_wprime_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(vprime_vprime_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(vprime_wprime_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(wprime_wprime_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))

      allocate(xtrbu_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(xsgsu_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(ysgsv_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(zsgsw_store(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(xtrbprime_uprime_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(xsgsprime_uprime_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(ysgsprime_vprime_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
      allocate(zsgsprime_wprime_bar(gpC%xsz(1),gpC%xsz(2), gpC%xsz(3)))
        
        !initialize all stored quantitites to zero
      fx_turb_store = 0
      fy_turb_store = 0
      fz_turb_store = 0
      
      fx_sgs_store = 0
      fy_sgs_store = 0
      fz_sgs_store = 0

      dPdx_store = 0
      dPdy_store = 0
      dPdz_store = 0
      P_store = 0
      u_store = 0
      v_store = 0
      w_store = 0

      uu_store = 0
      uv_store = 0
      uw_store = 0
      vv_store = 0
      vw_store = 0
      ww_store = 0
     
      xtrbu_store = 0
      xsgsu_store = 0
      ysgsv_store = 0
      zsgsw_store = 0
 
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

      ! Read in the fields (from restart files)
      computeFbody = .true.

      ! Initialize Padeder
      call Pade6opz%init(gpC, sp_gpC, gpE, sp_gpE, dz, NumericalSchemeVert,PeriodicInZ,spectC)

      ! Initialize sgs
      call newsgs%init(gpC, gpE, spectC, spectE, dx, dy, dz, inputfile, zMeshE(1,1,:), mesh(1,1,:,3), fbody_x, fbody_y, &
                      fbody_z, computeFbody, Pade6opZ, cbuffyC, cbuffzC, cbuffyE, cbuffzE, rbuffxC, rbuffyC, rbuffzC, &
                      rbuffyE, rbuffzE, Tsurf, ThetaRef, Fr, Re, .false., .false.,1)
      call newsgs%link_pointers(nuSGS, tauSGS_ij, tau13, tau23, q1, q2, q3, kappaSGS)


!subroutine init(Re, Pr, isInviscid, isStratified, botBC_temp)

   end subroutine



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

        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_",label,"_t",tid_final,".out"
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,arr,fname, gpC)

    end subroutine

!    subroutine readRestartFile(tid, rid)
!        use decomp_2d_io
!       use mpi
!        use exits, only: message
!        integer, intent(in) :: tid, rid
!        character(len=clen) :: tempname, fname
!        integer :: ierr
!
!        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_u.",tid
!        fname = trim(InputDir)//"/"//trim(tempname)
!        call decomp_2d_read_one(1,uC,fname, gpC)
!
!        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_v.",tid
!        fname = trim(InputDir)//"/"//trim(tempname)
!        call decomp_2d_read_one(1,vC,fname, gpC)
!
!        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_w.",tid
!        fname = trim(InputDir)//"/"//trim(tempname)
!        call decomp_2d_read_one(1,wE,fname, gpE)
!
        !if (this%isStratified) then
        !    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_T.",tid
        !    fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
        !    call decomp_2d_read_one(1,this%T,fname, this%gpC)
        !end if 

!        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",rid, "_info.",tid
!        fname = trim(InputDir)//"/"//trim(tempname)
!        open(unit=10,file=fname,access='sequential',form='formatted')
!        read (10, *)  tsim
!        close(10)

!        call mpi_barrier(mpi_comm_world, ierr)
!        call message("================= RESTART FILE USED ======================")
!        call message(0, "Simulation Time at restart:", tsim)
!        call message("=================================== ======================")

!    end subroutine
    
    subroutine readVisualizationFile(tid, rid)
        use decomp_2d_io
        use mpi
        use exits, only: message
        integer, intent(in) :: tid, rid
        character(len=clen) :: tempname, fname
        integer :: ierr
        character(len=4) :: label
        print *, 'inside readVisualizationFile, ind=', tid
        label = "uVel"
        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_",label,"_t",tid,".out"
        fname = trim(InputDir)//"/"//trim(tempname)
        call decomp_2d_read_one(1,uC,fname, gpC)

        !write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_v.",tid
        label = "vVel"
        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_",label,"_t",tid,".out"
        fname = trim(InputDir)//"/"//trim(tempname)
        call decomp_2d_read_one(1,vC,fname, gpC)

        !write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_w.",tid
        label = "wVel"
        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_",label,"_t",tid,".out"
        fname = trim(InputDir)//"/"//trim(tempname)
        call decomp_2d_read_one(1,wC,fname, gpC)
       
        !add in pressure
        label = "prss"
        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_",label,"_t",tid,".out"
        fname = trim(InputDir)//"/"//trim(tempname)
        call decomp_2d_read_one(1,pC,fname, gpC)
        
        !if (this%isStratified) then
        !    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_T.",ind
        !    fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
        !    call decomp_2d_read_one(1,this%T,fname, this%gpC)
        !end if 

        call transpose_x_to_y(wC, rbuffyC(:,:,:,1), gpC)
        call transpose_y_to_z(rbuffyC(:,:,:,1), rbuffzC(:,:,:,1), gpC)
        call Pade6opz%interpz_C2E(rbuffzC(:,:,:,1), rbuffzE(:,:,:,1),-1,-1)
        call transpose_z_to_y(rbuffzE(:,:,:,1),rbuffyE(:,:,:,1), gpE)
        call transpose_y_to_x(rbuffyE(:,:,:,1), wE, gpE)

        call mpi_barrier(mpi_comm_world, ierr)
        call message("================= RESTART FILE USED ======================")
        call message(0, "Simulation Time at restart:", tsim)
        call message("=================================== ======================")

    end subroutine

    subroutine getddx(f,dfdx)
        real(rkind), dimension(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)),  intent(in) :: f
        real(rkind), dimension(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)),  intent(out) :: dfdx
        
        call spectC%fft(f,cbuffyC(:,:,:,1))
        call spectC%mtimes_ik1_ip(cbuffyC(:,:,:,1))
        call spectC%ifft(cbuffyC(:,:,:,1),dfdx)    

    end subroutine

    subroutine getddy(f,dfdy)
        real(rkind), dimension(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)),  intent(in) :: f
        real(rkind), dimension(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)),  intent(out) :: dfdy
        
        call spectC%fft(f,cbuffyC(:,:,:,1))
        call spectC%mtimes_ik2_ip(cbuffyC(:,:,:,1))
        call spectC%ifft(cbuffyC(:,:,:,1),dfdy)    

    end subroutine


    subroutine getddz(f,dfdz)
        real(rkind), dimension(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)),  intent(in) :: f
        real(rkind), dimension(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)),  intent(out) :: dfdz
        
        call transpose_x_to_y(f,rbuffyC(:,:,:,1),gpC)
        call transpose_y_to_z(rbuffyC(:,:,:,1),rbuffzC(:,:,:,1),gpC)
        call Pade6opz%interpz_C2E(rbuffzC(:,:,:,1), rbuffzE(:,:,:,1),0,0)
        call Pade6opz%ddz_E2C(rbuffzE(:,:,:,1),rbuffzC(:,:,:,1),0,0)
        call transpose_z_to_y(rbuffzC(:,:,:,1),rbuffyC(:,:,:,1),gpC)
        call transpose_y_to_x(rbuffyC(:,:,:,1),dfdz,gpC)

    end subroutine

end program 
