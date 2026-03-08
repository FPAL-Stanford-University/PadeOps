module SolidGrid
    use kind_parameters,         only: rkind, clen
    use constants,               only: zero,eps,third,half,one,two,three,four,pi
    use FiltersMod,              only: filters
    use GridMod,                 only: grid
    use gridtools,               only: alloc_buffs, destroy_buffs
    use sgrid_hooks,             only: meshgen, initfields, hook_output, hook_bc, hook_timestep, hook_mixture_source, get_sponge
    use decomp_2d,               only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                                       transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,          only: derivatives
    use DerivativesStaggeredMod, only: derivativesStagg
    use InterpolatorsMod,        only: interpolators
    use LADMod,                  only: ladobject
    use StiffGasEOS,             only: stiffgas
    use Sep1SolidEOS,            only: sep1solid
    use SolidMixtureMod,         only: solid_mixture
    use IOsgridMod,              only: IOsgrid

    use operators,               only: filter3D
   
    implicit none

    integer, parameter :: rho_index    = 1 
    integer, parameter :: u_index      = 2
    integer, parameter :: v_index      = 3
    integer, parameter :: w_index      = 4
    integer, parameter :: p_index      = 5
    integer, parameter :: T_index      = 6
    integer, parameter :: e_index      = 7
    integer, parameter :: sos_index    = 8
    integer, parameter :: mu_index     = 9
    integer, parameter :: bulk_index   = 10
    integer, parameter :: kap_index    = 11
    integer, parameter :: sxx_index    = 12
    integer, parameter :: sxy_index    = 13
    integer, parameter :: sxz_index    = 14
    integer, parameter :: syy_index    = 15
    integer, parameter :: syz_index    = 16
    integer, parameter :: szz_index    = 17
    integer, parameter :: uint_index   = 18
    integer, parameter :: vint_index   = 19
    integer, parameter :: wint_index   = 20
    integer, parameter :: uJ_index            = 21
    integer, parameter :: vJ_index            = 22
    integer, parameter :: wJ_index            = 23
    integer, parameter :: keJ_index           = 24
    integer, parameter :: eJ_index            = 25
    integer, parameter :: tauxx_index         = 26
    integer, parameter :: tauyy_index         = 27
    integer, parameter :: tauzz_index         = 28
    integer, parameter :: tauxy_index         = 29
    integer, parameter :: tauyx_index         = 30
    integer, parameter :: tauxz_index         = 31
    integer, parameter :: tauzx_index         = 32
    integer, parameter :: tauyz_index         = 33
    integer, parameter :: tauzy_index         = 34
    integer, parameter :: dudx_index          = 35
    integer, parameter :: dudy_index          = 36
    integer, parameter :: dudz_index          = 37
    integer, parameter :: dvdx_index          = 38
    integer, parameter :: dvdy_index          = 39
    integer, parameter :: dvdz_index          = 40
    integer, parameter :: dwdx_index          = 41
    integer, parameter :: dwdy_index          = 42
    integer, parameter :: dwdz_index          = 43 
    integer, parameter  :: tauSum_index       = 44
    integer, parameter  :: metric_exact_index = 45
    integer, parameter  :: metric_index       = 46
    integer, parameter  :: metric_N2F_index   = 47
    integer, parameter  :: metric_half_index  = 48
    integer, parameter  :: uref_index         = 49
    integer, parameter  :: m1int_index        = 50
    integer, parameter  :: m2int_index        = 51
    integer, parameter  :: nfields = 51

    integer, parameter :: mom_index = 1
    integer, parameter :: TE_index = mom_index+3
    integer, parameter :: ncnsrv  = TE_index

    ! These indices are for data management, do not change if you're not sure of what you're doing
    integer, parameter :: tauxyidx = 2
    integer, parameter :: tauxzidx = 3
    integer, parameter :: tauyzidx = 6
    integer, parameter :: tauxxidx = 4
    integer, parameter :: tauyyidx = 7
    integer, parameter :: tauzzidx = 8

    integer, parameter :: qxidx = 1
    integer, parameter :: qyidx = 5
    integer, parameter :: qzidx = 9


    ! Number of buffers to create
    integer, parameter :: nbufsx = 2
    integer, parameter :: nbufsy = 6
    integer, parameter :: nbufsz = 2
 

    type, extends(grid) :: sgrid
       
        type(filters),          allocatable :: gfil
        type(derivatives),      allocatable :: derD02, derD06,derD04,derCD06,der_nostretch,derCD06_nostretch,derCD04
        type(solid_mixture),    allocatable :: mix
        type(ladobject),        allocatable :: LAD
        type(derivativesStagg), allocatable :: derStagg,derStagg_stretch
        type(derivativesStagg), allocatable :: derStaggd02,derStaggd04
        type(interpolators),    allocatable ::interpMid,interpMide06
        type(interpolators),    allocatable ::interpMid02, interpMid04,interpMid08,interpMid06
        type( IOsgrid ),        allocatable :: viz

        logical     :: PTeqb                       ! Use pressure and temperature equilibrium formulation
        logical     :: pEqb                        ! nterpolators),    allocatable ::interpMidUse pressure equilibrium formulation
        logical     :: pRelax                      ! Use pressure and temperature non-equilibrium formulation, but relax pressure at each substep
        logical     :: use_gTg                     ! Use formulation with the Finger tensor g^T.g instead of the full g tensor
        logical     :: cnsrv_g, cnsrv_gt, cnsrv_gp, cnsrv_pe ! use conservative form of equations
        logical     :: strainHard                  ! use strainHardening
        logical     :: updateEtot                  ! Update species etot (vs ehydro) with pRelax
        logical     :: useOneG                     ! Use formulation with a single g or gTg field
        logical     :: useNC
        logical     :: intSharp                    ! Include interface sharpening terms
        logical     :: intSharp_cpl                ! Include coupling of sharpening with momentum and energy equations
        logical     :: intSharp_cpg                ! Include coupling of sharpening with kinematic equations
        logical     :: intSharp_cpg_west           ! Use form of kinematic sharpening terms derived by Jacob West
        logical     :: intSharp_spf                ! Use Shukla-Pantano-Freund method - not in divergence form
        logical     :: intSharp_ufv                ! Use finite volume discretization for sharpening term
        logical     :: intSharp_utw                ! Use Tiwari formulation
        logical     :: usePhiForm
        logical     :: twoPhaseLAD                 ! Use dYs/dx instead of the fickian 
        logical     :: LAD5eqn
        logical     :: useEigenFunction
        logical     :: FilteredTimeStep 
        real(rkind) :: intSharp_gam                ! Interface sharpening Gamma parameter
        real(rkind) :: intSharp_eps                ! Interface sharpening epsilon parameter
        real(rkind) :: intSharp_cut                ! Interface sharpening cutoff parameter, for VF approaching 1 or 0
        real(rkind) :: intSharp_dif                ! Interface sharpening VF out of bounds diffusion
        real(rkind) :: intSharp_tnh                ! Interface sharpening blending parameter
        real(rkind) :: intSharp_pfloor              ! Pressure floor for pressure-temperature relaxation / LAD
        real(rkind) :: intSharp_tfloor              ! Temperature floor for pressure-temperature relaxation / LAD
        real(rkind) :: alpha_skew
        logical :: intSharp_d02                    ! Use 2nd order
        logical :: intSharp_msk                    ! Mask FV diffusion
        logical :: intSharp_flt                    ! Use dealliasing filter for interface sharpening derivatives
        logical :: intSharp_flp                    ! Filter pressure
        logical :: skew_mass, skew_Ys, skew_VF
        logical     :: useAkshayForm, SpongeLayer	
        logical     :: weightedcurvature
	logical     :: surface_mask
        logical     :: LADInt,LADN2F,LADMass_Consv
	logical     :: use_FV, use_D04, use_Stagg, use_XiLS         !flag to use FV in surface tension scheme
	logical     :: use_gradphi,energy_surfTen          !flag to use phi formulation in surface tension calculation
	logical     :: use_gradVF           !flag to use VF formulation in surface tension calculation		
        logical     :: use_gradXi
        logical     :: use_surfaceTension   !flag to turn on/off surface tension (in momentum and energy equations)
        logical     :: use_normFV
        logical     :: use_normInt
        logical     :: use_CnsrvSurfaceTension
        logical     :: Stretch1D, Stretch1Dy, Stretch1Dx, Stretch1Dz
        real(rkind) :: surfaceTension_coeff !constant coefficient for surface tension
        real(rkind) :: R, p_amb, XiLS_eps
        logical :: filt_mask = .FALSE.             ! mask filter in high gradient regions
        real(rkind), dimension(:,:,:,:), allocatable :: filt_tmp,filt_grad  ! temporary for filter mask
        real(rkind), dimension(:,:,:), allocatable :: filt_thrs  ! temporary for filter mask
        real(rkind) :: filt_cut  ! bulk threshold for filter mask

        real(rkind), dimension(:,:,:,:), allocatable :: Wcnsrv                               ! Conserved variables
        real(rkind), dimension(:,:,:,:), allocatable :: xbuf, ybuf, zbuf   ! Buffers
        real(rkind), dimension(:,:,:,:), allocatable :: u_mid, v_mid, w_mid, p_mid,ke_mid,pu_mid,rho_mid
        real(rkind), dimension(:,:,:),   pointer     :: rho_int, u_int, v_int,w_int,tauxy_int, tauyy_int, tauyz_int,qy_int, e_int, TE, p_int,VF_int, m1_int,m2_int
        real(rkind), dimension(:,:,:),   pointer     :: xflux_x, yflux_x, zflux_x, xflux_y, yflux_y, zflux_y, xflux_z, yflux_z, zflux_z, yflux_e, xflux_e, zflux_e

        real(rkind), dimension(:,:,:), pointer :: x 
        real(rkind), dimension(:,:,:), pointer :: y 
        real(rkind), dimension(:,:,:), pointer :: z 
       
        real(rkind), dimension(:,:,:), pointer :: eta1
        real(rkind), dimension(:,:,:), pointer :: eta2
        real(rkind), dimension(:,:,:), pointer :: eta3
 
        real(rkind), dimension(:,:,:), pointer :: rho 
        real(rkind), dimension(:,:,:), pointer :: u 
        real(rkind), dimension(:,:,:), pointer :: v 
        real(rkind), dimension(:,:,:), pointer :: w 
        real(rkind), dimension(:,:,:), pointer :: p 
        real(rkind), dimension(:,:,:), pointer :: T 
        real(rkind), dimension(:,:,:), pointer :: e 
        real(rkind), dimension(:,:,:), pointer :: sos,rhoe 
        real(rkind), dimension(:,:,:), pointer :: mu 
        real(rkind), dimension(:,:,:), pointer :: bulk 
        real(rkind), dimension(:,:,:), pointer :: kap, eLAD
        real(rkind), dimension(:,:,:), pointer :: tauaiidivu
        real(rkind), dimension(:,:,:), pointer :: pmix, intP
        real(rkind), dimension(:,:,:,:), pointer :: devstress
        real(rkind), dimension(:,:,:), pointer :: sxx
        real(rkind), dimension(:,:,:), pointer :: sxy
        real(rkind), dimension(:,:,:), pointer :: sxz
        real(rkind), dimension(:,:,:), pointer :: syy
        real(rkind), dimension(:,:,:), pointer :: syz
        real(rkind), dimension(:,:,:), pointer :: szz
       
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxxe, dudx
        real(rkind), dimension(:,:,:), pointer :: tauyy, tauyye, dudy, dudy2, dvfdy, dmudy, drhody
        real(rkind), dimension(:,:,:), pointer :: tauzz, tauzze, dudz
        real(rkind), dimension(:,:,:), pointer :: tauxy, tauxye, dvdx
        real(rkind), dimension(:,:,:), pointer :: tauyx, tauyxe, dvdy
        real(rkind), dimension(:,:,:), pointer :: tauxz, tauxze, dvdz
        real(rkind), dimension(:,:,:), pointer :: tauzx, tauzxe, dwdx
        real(rkind), dimension(:,:,:), pointer :: tauyz, tauyze,dwdy,metric_half, metric_N2F
        real(rkind), dimension(:,:,:), pointer :: tauzy, tauzye, dwdz,tauSum,esum, esumJ, metric, metric_exact
        real(rkind), dimension(:,:,:), pointer :: fsw,divgrad
        real(rkind), dimension(:,:,:), pointer :: rhouHeur,rhovHeur,rhowHeur,rhoeHeur,m1heur,m2heur,VFheur,entropy,discreteKE,puKE,SurfTenDiff
        real(rkind), dimension(:,:,:), pointer :: keJ, uJ, vJ, wJ, eJ, qDiv,pEvolve, VFEvolve, pError, VFerror, pJ, tauRho,uref
        real(rkind) :: phys_mu1, phys_mu2,CP,sos_ratio, sos_ref=4.28646795041d0
        real(rkind) :: phys_bulk1, phys_bulk2
        real(rkind) :: phys_kap1, phys_kap2, g 
        real(rkind) :: st_limit, pthick,uthick,rhothick,Ys_wiggle,VF_wiggle,VF_thick,Ys_thick
        real(rkind), dimension(2) :: rhou_ref,rhov_ref,rhow_ref, rhoe_ref
        real(rkind), dimension(10) :: lamru,lamrv,lamrw,lamre,lamm1,lamm2,lamvf,lamtim
        real(rkind), dimension(:,:,:,:), allocatable :: meshstretch, sponge
        real(rkind), dimension(:,:,:), allocatable :: yMetric,xMetric, zMetric,yMetric_half, xMetric_half, zMetric_half,yLADMetric, yMetric_F2N, dy_stretch
        integer :: stepfil,numfil
        contains
            procedure          :: init
            procedure          :: destroy
            procedure          :: laplacian
            procedure          :: gradient 
            procedure          :: secondder
            procedure          :: advance_RK45
            procedure          :: simulate
            procedure          :: update_p
            procedure          :: getRHS_P
            procedure          :: checkTau
            procedure, private :: get_dt(stability)
            procedure, private  :: get_dtlocal(stability)
        else
          call this%get_dtlocal(stability)
        endif
        
        !populate surface tension terms at initial condition
        if(this%use_surfaceTension) then
            if(this%mix%ns.ne.2) then
                call GracefulExit("Surface tension is not defined for single-species, and not implemented for more than 2 species",4634)
            endif

             call  this%mix%get_surfaceTension(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w,this%x,this%y,1)  ! Compute surface tension terms for momentum and energy equations
        !      call this%mix%get_gradp(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)

         !       this%mix%surfaceTension_f = this%mix%gradp
        endif
        !call this%mix%Test_Der_NP(this%x,this%y,this%z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !call this%mix%Test_1DStretch(this%x,this%y,this%z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        this%metric = this%yMetric
        this%metric_exact = this%yLADMetric
        this%metric_N2F   = this%yMetric_F2N
        this%metric_half  = this%yMetric_half
        !call this%mix%Test_Der(this%x,this%y,this%z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !call this%mix%Test_Der_Periodic(this%x,this%y,this%z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        if(this%use_CnsrvSurfaceTension) then
            if(this%mix%ns.ne.2) then
                call GracefulExit("Surface tension is not defined for single-species, and not implemented for more than 2 species",4634)
            endif

            call  this%mix%get_surfaceTensionCnsrv(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)
            ! Compute surface tension terms for momentum and energy equations

        endif

   
        ! Write out initial conditions
        ! call hook_output(this%decomp, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
!        call hook_output(this%decomp,this%der,this%dx,this%dy,this%dz,this%outputdir,this%mesh,this%fields,this%mix,this%tsim,this%viz%vizcount,this%pthick,this%uthick,this%rhothick,this%Ys_thick,this%VF_thick,this%Ys_wiggle,this%VF_wiggle,this%x_bc,this%y_bc,this%z_bc)
        if( this%Stretch1Dy) then
           call this%viz%WriteViz(this%decomp, this%meshstretch, this%fields, this%mix, this%tsim)
        else
           call this%viz%WriteViz(this%decomp, this%mesh, this%fields, this%mix,this%tsim)
        endif

        vizcond = .FALSE.
       
        ! Check for visualization condition and adjust time step
        if ( (this%tviz > zero) .AND. (this%tsim + this%dt > this%tviz * this%viz%vizcount) ) then
            this%dt = abs(this%tviz * this%viz%vizcount - this%tsim)
            vizcond = .TRUE.
            stability = 'vizdump'
        end if

        tcond = .TRUE.
        ! Check tstop condition
        if ( (this%tstop > zero) .AND. (this%tsim >= this%tstop) ) then
            tcond = .FALSE.
        else if ( (this%tstop > zero) .AND. (this%tsim + this%dt >= this%tstop) ) then
            this%dt = this%tstop - this%tsim
        end if

        ! Check nsteps condition
        if ( (this%nsteps <= 0) .OR. (this%step < this%nsteps) ) then
            stepcond = .TRUE.
        else
            stepcond = .FALSE.
        end if

        if ( (this%tstop <= zero) .AND. (this%nsteps <= 0) ) then
            call GracefulExit('No stopping criterion set. Set either tstop or nsteps to be positive.', 345)
        end if

        ! Start the simulation while loop
        if(nrank==0) write(*,*) 'Starting time loop'
        do while ( tcond .AND. stepcond )
            ! Advance time
            call tic()
            call this%advance_RK45()
            call toc(cputime)

            !call this%mix%thick_calculations(this%rho, this%p,this%u,this%pthick,this%uthick,this%rhothick,this%Ys_thick,this%VF_thick,this%Ys_wiggle,this%VF_wiggle, this%dx)          
!            u_max = P_MAXVAL(this%u)
!            u_min = P_MINVAL(this%u)
!            v_max = P_MAXVAL(this%v)
!            v_min = P_MINVAL(this%v)
!            w_max = P_MAXVAL(this%w)
!            w_min = P_MINVAL(this%w)
!            p_max = P_MAXVAL(this%p)
!            p_min = P_MINVAL(this%p)
!            rho_max = P_MAXVAL(this%rho)
!            rho_min = P_MINVAL(this%rho)
!            Ys_max = P_MAXVAL(this%mix%material(1)%Ys)
!            Ys_min = P_MINVAL(this%mix%material(1)%Ys)
!            VF_max = P_MAXVAL(this%mix%material(1)%VF)
!            VF_min = P_MINVAL(this%mix%material(1)%VF)

            call message(1,"Time",this%tsim)
            call message(1,"Step",this%step)
            call message(2,"Time step",this%dt)
            call message(2,"Stability limit: "//trim(stability))
            call message(2,"CPU time (in seconds)",cputime)
   !         call message(3, " u min ", u_min)
   !         call message(3, " u max ", u_max)
   !         call message(3, " v min ", v_min)
   !         call message(3, " v max ", v_max)
   !         call message(3, " w min ", w_min)
   !         call message(3, " w max ", w_max)
   !         call message(3, " p min ", p_min)
   !         call message(3, " p max ", p_max)
   !         call message(3, " rho min ", rho_min)
   !         call message(3, " rho max ", rho_max)
   !         call message(3, " Ys min ", Ys_min)
   !         call message(3, " Ys max ", Ys_max)
   !         call message(3, " VF min ", VF_min)
   !         call message(3, " VF max ", VF_max)
            call hook_timestep(this%decomp, this%mesh, this%fields, this%mix, this%step, this%tsim)
            ! Write out vizualization dump if vizcond is met 
           if (vizcond) then
                ! call hook_output(this%decomp, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
!                call hook_output(this%decomp,this%der,this%dx,this%dy,this%dz,this%outputdir,this%mesh,this%fields,this%mix,this%tsim,this%viz%vizcount,this%pthick,this%uthick,this%rhothick,this%Ys_thick,this%VF_thick,this%Ys_wiggle,this%VF_wiggle,this%x_bc,this%y_bc,this%z_bc)

                if( this%Stretch1Dy) then               
                   call this%viz%WriteViz(this%decomp, this%meshstretch, this%fields, this%mix, this%tsim)
                else

                !   call this%FilteringHeuristicHighOrder()
                !   call this%LocalDiffHeuristic()
                !   call this%BicubicMetric()
                !  call this%FilDiffHeuristic() 
                   call this%viz%WriteViz(this%decomp, this%mesh, this%fields,this%mix, this%tsim)
                !   call this%FilteringHeuristic()
                !   call this%BicubicMetric()
                endif
                vizcond = .FALSE.
           end if
            
            ! Get the new time step
            call this%get_dtlocal(stability)
            ! Check for visualization condition and adjust time step
            if ( (this%tviz > zero) .AND. (this%tsim + this%dt >= this%tviz * this%viz%vizcount) ) then
                this%dt = abs( this%tviz * this%viz%vizcount - this%tsim)
                vizcond = .TRUE.
            end if


           if ( restartWrite .or. (mod(this%step,this%t_restartDump) == 0) ) then
                 call this%dumpRestartfile()
                 
                 call message(0,"Scheduled restart file dumped.")
            end if
            ! Check tstop condition
            if ( (this%tstop > zero) .AND. (this%tsim >= this%tstop*(one - eps)) ) then
                tcond = .FALSE.
            else if ( (this%tstop > zero) .AND. (this%tsim + this%dt >= this%tstop*(one - eps)) ) then
                this%dt = this%tstop - this%tsim
                stability = 'stop'
                vizcond = .TRUE.
            end if

            ! Check nsteps condition
            if ( (this%nsteps <= 0) .OR. (this%step < this%nsteps) ) then
                stepcond = .TRUE.
            else
                stepcond = .FALSE.
            end if

            ! Check for exitpdo file
            if(check_exit(this%outputdir)) then
                ! call hook_output(this%decomp, this%dx, this%dy, this%dz, this%outputdir, this%mesh, this%fields, this%mix, this%tsim, this%viz%vizcount)
                !call hook_output(this%decomp,this%der,this%dx,this%dy,this%dz,this%outputdir,this%mesh,this%fields,this%mix,this%tsim,this%viz%vizcount,this%pthick,this%uthick,this%rhothick,this%Ys_thick,this%VF_thick,this%Ys_wiggle,this%VF_wiggle,this%x_bc,this%y_bc,this%z_bc)
                if( this%Stretch1Dy) then

                   call this%viz%WriteViz(this%decomp, this%meshstretch, this%fields, this%mix, this%tsim)

                else
                
                   call this%viz%WriteViz(this%decomp, this%mesh, this%fields,this%mix, this%tsim) 

                endif

                call GracefulExit("Found exitpdo file in working directory",1234)
            endif

        end do

         call hook_timestep(this%decomp, this%mesh, this%fields, this%mix, this%step, this%tsim)
    end subroutine

    subroutine advance_RK45(this)
        use RKCoeffs,   only: RK45_steps,RK45_A,RK45_B,RK3_steps,RK3_A,RK3_B
        use timer,      only: tic, toc
        use exits,      only: message,nancheck,GracefulExit
        use reductions, only: P_MAXVAL, P_MINVAL
        use decomp_2d,  only: nrank
        use operators, only: divergence,gradient,gradFV_x,gradFV_y,gradFV_z,interpolateFV
        use constants,               only: pi
        class(sgrid), target, intent(inout) :: this


        real(rkind)                                               :: Qtmpt      ! Temporary variable for RK45
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv) :: rhs        ! RHS for conserved variables
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,ncnsrv) :: Qtmp             ! Temporary variable for RK45
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: divu,Qtmpp, pmix ! Velocity divergence for species energy eq
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: viscwork         ! Viscous work term for species energy eq
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: Fsource, tmp, eta, tmp2,rhofil,efil,m1fil,m2fil,TEfil,rhoufil,rhovfil,rhowfil,VFfil,H1,H2 
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx,drudy,drudz,drvdx,drvdy,drvdz,drwdx,drwdy,drwdz,dredx,dredy,dredz,dVFdx,dVFdy,dVFdz,dm1dx,dm1dy,dm1dz,dm2dx,dm2dy,dm2dz,tmp1,tmp3
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx4,drudy4,drudz4,drvdx4,drvdy4,drvdz4,drwdx4,drwdy4,drwdz4,dredx4,dredy4,dredz4,dVFdx4,dVFdy4,dVFdz4,dm1dx4,dm1dy4,dm1dz4,dm2dx4,dm2dy4,dm2dz4
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3)      :: tmpint
        integer :: isub,i,j,k,l,imat,iter,ii,jj,kk
        real(rkind), dimension(:,:,:,:), allocatable, target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind) ::cputime, Ystime, VFtime, consvTime, intsharp, u_max, u_min,rhoYs_min, rhoYs_max,tmod,Ut
        
        character(len=clen) :: charout




        call this%get_conserved()
        Qtmp  = this%Wcnsrv
        Qtmpt = this%tsim
        pmix  = zero
     

       do isub = 1, RK3_steps

            if(this%use_CnsrvSurfaceTension) then
                call this%mix%get_surfaceTensionPE(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)
            endif
            call this%get_conserved()


            if(this%use_surfaceTension) then
                if(this%mix%ns.ne.2) then
                    call GracefulExit("Surface tension is not defined for single-species, and not implemented for more than 2 species",4634)
                endif

                !call tic()                
                call this%mix%get_surfaceTension(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w,this%x,this%y,isub)
                !call toc(cputime)
                !call message(3,"Surface Tension time (in seconds)",cputime)
! Compute surface tension terms for momentum and energy equations

            endif

            where( this%mix%material(1)%VF .GT. 1d-10 )
               this%mix%deltakap =abs(this%mix%kappaNoFil*(this%dy_stretch)) !*abs(this%mix%material(1)%VF*(1-this%mix%material(1)%VF) )*4.0           
            elsewhere
               this%mix%deltakap = 0
            endwhere

            do i = 1,2
             
             this%mix%material(i)%deltakap = this%mix%deltakap

            enddo       

            !call tic()
            call this%getFaces()
            !call toc(cputime)
            !call message(3,"Faces time (in seconds)",cputime)

            if ( nancheck(this%Wcnsrv,i,j,k,l) ) then
                call message("Wcnsrv: ",this%Wcnsrv(i,j,k,l))
                !write(charout,'(A,I1,A,I5,A,4(I5,A))') "NaN encountered in solution (Wcnsrv) at substep ", isub, " of step ", this%step+1, " at (",i,", ",j,", ",k,", ",l,") of Wcnsrv"
                write(charout,'(A,I1,A,I5,A,4(I5,A))') "NaN encountered in solution (Wcnsrv) at substep ", isub, " of step ", this%step+1, " at (",i+this%decomp%yst(1)-1,", ",j+this%decomp%yst(2)-1,", ",k+this%decomp%yst(3)-1,", ",l,") of Wcnsrv"
                call GracefulExit(trim(charout), 999)
            end if
            ! Pre-compute stress, LAD, J, etc.
            call this%mix%getSOS(this%rho,this%p,this%sos)

            !call tic()
            call this%mix%getLAD(this%rho,this%p,this%e,this%u, this%v, this%w,this%sos,this%yMetric,this%dy_stretch,this%use_gTg,this%strainHard,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%intSharp_tfloor,this%dt)  ! Compute species LAD (kap, diff, diff_g, diff_gt,diff_pe)

            do imat = 1,this%mix%ns
                call this%mix%material(imat)%LAD_Quant(this%rho,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            enddo
            !call toc(cputime)
            !call message(3,"LAD time (in seconds)",cputime)


            if(this%intSharp) then
               if(this%mix%ns.ne.2) then
                  call GracefulExit("Problem if ns=1, should work for ns>2 but not tested",4634)
               endif
               !endif
               do imat=1,this%mix%ns
                   this%mix%material(imat)%intSharp_a = zero
                   this%mix%material(imat)%intSharp_aDiff = zero
                   this%mix%material(imat)%intSharp_aFV = zero !-- problem?
                   this%mix%material(imat)%intSharp_R = zero
                   this%mix%material(imat)%intSharp_RDiff = zero
                   this%mix%material(imat)%intSharp_RFV = zero
               enddo

               this%mix%intSharp_f = zero
               this%mix%intSharp_fDiff = zero
               this%mix%intSharp_fFV = zero
               this%mix%intSharp_h = zero
               this%mix%intSharp_hDiff = zero
               this%mix%intSharp_hFV = zero
               this%mix%intSharp_kFV = zero
             
               !call tic()
               call this%mix%get_intSharp_clean2(this%rho,this%ke_mid,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w,this%p,this%u_mid,this%v_mid,this%w_mid,this%p_mid)
               !call toc(cputime)
               !call message(3,"Interface Sharpening Time (in seconds)",cputime)

       else      
                  ! !debug
                   do imat=1,this%mix%ns
                      this%mix%material(imat)%intSharp_a = zero
                      this%mix%material(imat)%intSharp_aDiff = zero
                      this%mix%material(imat)%intSharp_aFV = zero !-- problem?
                      this%mix%material(imat)%intSharp_R = zero
                      this%mix%material(imat)%intSharp_RDiff = zero
                      this%mix%material(imat)%intSharp_RFV = zero
                   enddo
                  this%mix%intSharp_f = zero
                  this%mix%intSharp_fDiff = zero
                  this%mix%intSharp_fFV = zero
                  this%mix%intSharp_h = zero
                  this%mix%intSharp_hDiff = zero
                  this%mix%intSharp_hFV = zero
                  this%mix%intSharp_kFV = zero
                  ! !end debug
                  

            endif

             call this%mix%checkNaN()

             if(this%use_CnsrvSurfaceTension) then
                if(this%mix%ns.ne.2) then
                    call GracefulExit("Surface tension is not defined forsingle-species, and not implemented for more than 2 species",4634)
                endif

                call this%mix%get_surfaceTensionCnsrv(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)
                ! Compute surface tension terms for momentum and energy equations

            endif

            ! Update total mixture conserved variables


            !call tic()
            if (this%useNC) then
              call this%getRHS_NC(rhs,divu, viscwork)
            else
             call this%getRHS(rhs,divu,viscwork)
            endif
            !call toc(cputime)
            !call message(3,"RHS time (in seconds)",cputime)


            !!!!!!!!!!!!!! UNCOMMENT            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           this%Wcnsrv = RK3_B(isub)*(this%Wcnsrv + this%dt*rhs) + RK3_A(isub)*Qtmp
            !Qtmp  = this%dt*rhs  + RK45_A(isub)*Qtmp
            !this%Wcnsrv = this%Wcnsrv + RK45_B(isub)*Qtmp

           !!!!!!!!!!!!!!!!!!! UNCOMMENT       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! calculate sources if they are needed
           ! if(.not. this%PTeqb) call this%mix%calculate_source(this%rho,divu,this%u,this%v,this%w,this%p,Fsource,this%x_bc,this%y_bc,this%z_bc) ! -- actually, source terms should be included for PTeqb as well --NSG

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMMENTED OUT UPDATE G             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Now update all the individual species variables
            !check eps
            !call this%mix%update_g(isub,max(this%dt,eps),this%rho,this%u,this%v,this%w,this%x,this%y,this%z,Fsource,this%tsim,this%x_bc,this%y_bc,this%z_bc)               ! g tensor

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            !call tic()
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(this%use_Stagg) then
                
               call this%mix%update_Ys(isub,this%dt,this%rho,this%u,this%v,this%w,this%u_int,this%v_int,this%w_int,this%sos,this%x,this%y,this%z,this%tsim,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%sponge,this%alpha_skew)               ! Volume Fraction
!               call this%update_P(Qtmpp,isub,this%dt,this%x,this%y,this%z,this%tsim,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
             
            !   call this%mix%material(1)%update_specYs(isub,this%dt,this%rho,this%u_mid(:,:,:,1),this%v_mid(:,:,:,2),this%w_mid(:,:,:,3),this%x,this%y,this%z,this%tsim,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%sponge,this%alpha_skew)  
            else
               call this%mix%update_Ys(isub,this%dt,this%rho,this%u,this%v,this%w,this%u_mid,this%v_mid,this%w_mid,this%sos,this%x,this%y,this%z,this%tsim,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%sponge,this%alpha_skew)
            endif
            !call toc(cputime)
            !call message(3,"Ys time (in seconds)",cputime)

            !call tic()
            if(this%pEqb) then


               if(this%use_Stagg) then
                  call divergence(this%decomp,this%der,this%u,this%v,this%w,divu,this%x_bc,this%y_bc,this%z_bc)
!                  call
!                  this%mix%update_VF(isub,this%dt,this%rho,this%u,this%v,this%w,this%u_mid(:,:,:,1),this%v_mid(:,:,:,2),this%w_mid(:,:,:,3),this%sos,this%x,this%y,this%z,this%tsim,divu,Fsource,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%sponge,this%alpha_skew)
!                  call
!                  this%mix%update_VF(isub,this%dt,this%rho,this%u,this%v,this%w,this%x,this%y,this%z,this%tsim,divu,Fsource,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%sponge,this%alpha_skew)


               else
             !     call
             !     this%mix%update_VF(isub,this%dt,this%rho,this%u,this%v,this%w,this%u_mid(:,:,:,1),this%v_mid(:,:,:,2),this%w_mid(:,:,:,3),this%sos,this%x,this%y,this%z,this%tsim,divu,Fsource,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%sponge,this%alpha_skew)
               endif
!               call
!               this%update_P(Qtmpp,isub,this%dt,this%x,this%y,this%z,this%tsim,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
               call this%mix%update_VF(isub,this%dt,this%rho,this%u,this%v,this%w,this%u_int,this%v_int,this%w_int,this%sos,this%x,this%y,this%z,this%tsim,divu,Fsource,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%sponge,this%alpha_skew)
            elseif(this%pRelax) then
                call this%mix%update_VF(isub,this%dt,this%rho,this%u,this%v,this%w,this%u_int,this%v_int,this%w_int,this%sos,this%x,this%y,this%z,this%tsim,divu,Fsource,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc,this%sponge,this%alpha_skew)
! Volume Fraction
                call this%mix%update_eh(isub,this%dt,this%rho,this%u,this%v,this%w,this%x,this%y,this%z,this%tsim,divu,viscwork,Fsource,this%devstress,this%x_bc,this%y_bc,this%z_bc)
! Hydrodynamic energy
            end if
            !call toc(cputime)
            !call message(3,"VF time (in seconds)",cputime)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMENT             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !if (.NOT. this%PTeqb) then
           this%tsim = RK3_B(isub)*( this%tsim + this%dt)  + RK3_A(isub)*Qtmpt

            
!           Qtmpt = this%dt + RK45_A(isub)*Qtmpt
!           this%tsim = this%tsim + RK45_B(isub)*Qtmpt
           

             
!            if(this%tsim .GE. 0.11) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UNCOMMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              if(.NOT. this%use_Stagg) then
                  ! Filter the conserved variables
!                 call this%filter(this%Wcnsrv(:,:,:,mom_index  ), this%fil, 1,-this%x_bc, this%y_bc, this%z_bc)
!                 call this%filter(this%Wcnsrv(:,:,:,mom_index+1), this%fil, 1, this%x_bc,-this%y_bc, this%z_bc)
!                 call this%filter(this%Wcnsrv(:,:,:,mom_index+2), this%fil, 1, this%x_bc, this%y_bc,-this%z_bc)
!                 call this%filter(this%Wcnsrv(:,:,:, TE_index  ), this%fil, 1, this%x_bc, this%y_bc, this%z_bc)
!              call this%filter(this%p, this%fil, 1,this%x_bc, this%y_bc, this%z_bc)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! Filter the individual species variables
!                 call this%mix%filter(1, this%x_bc, this%y_bc, this%z_bc)
!            end if 
           !  endif
!          print *, "Filter"
!          endif


           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(this%use_CnsrvSurfaceTension) then

               call this%mix%getYs(this%rho)
               call this%mix%get_surfaceTensionPE(this%rho,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz,this%u,this%v,this%w)

            endif

            call this%get_primitive()
            

            
            if (this%PTeqb) then
               !call this%mix%equilibratePressureTemperature(this%rho, this%e, this%p, this%T, isub)
               ! do i=1,2


            elseif (this%pEqb) then


              call this%mix%equilibratePressure(this%rho,this%e, this%p)

            elseif (this%pRelax) then
                call this%mix%relaxPressure(this%rho, this%e, this%p)
            end if


            call hook_bc(this%decomp, this%mesh, this%fields, this%mix, this%tsim, this%x_bc, this%y_bc, this%z_bc)
            

            if(this%pEqb) then
            call this%post_bc_2()
            else
            call this%post_bc()
            endif

         end do
        

        this%step = this%step + 1
    end subroutine

    subroutine getFaces( this)
        use decomp_2d,  only: nrank
        use operators, only: interpolateFV,interpolateFV_x,interpolateFV_F2Ny,interpolateFV_y,interpolateFV_F2Nx,filter3D,gradFV_N2Fx,gradFV_N2Fy,gradFV_N2Fz
        use timer, only: tic, toc
        use exits,      only: message,nancheck,GracefulExit
        class(sgrid), target, intent(inout) :: this
        integer :: i,j,k,iflag = one, nx,nz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: T1_int,T2_int,rhom_int,rhom2_int,rhoe1_int,rhoe2_int,rhoe_int,rhou_int,rhov_int,rhow_int,spec_int,pgam,T_int,rhoc_int,rhocp_int,mu_int,mv_int,mw_int,sos_int,VFbar,kappabar,gradp,gradVF,peff,peff_fil
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: psi_int,e1_int, Gam_int,pVF_int,num,denom,af,c1_int,rhobar
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: xhalf,tmp,rhom,tmp2, tmpfil,tmp3,psi,Gam,psi_safe,c1,rhom1,tmp1,mask1,mask2,mask3
        real(rkind),dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xtmp1,xtmp2,xtmp3,xtmp4,xtmp5
        real(rkind),dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: ztmp1,ztmp2,ztmp3,ztmp4,ztmp5
        real(rkind) :: e = 1d-10,ref_ratio=1d3, rho_ratio, cputime

          

        !!!!!!!!!!!!!!!!!!!!!!!!!! Calculate Ratios !!!!!!!!!!!!!!!!!!!!!!!!!!
        rho_ratio =(this%mix%material(1)%elastic%rho0 / this%mix%material(2)%elastic%rho0 )/ref_ratio
        !this%sos_ratio = this%sos_ratio /sos_ref 
        nx = this%decomp%xsz(1)
        nz = this%decomp%zsz(3)
        call interpolateFV(this%decomp,this%interpMid,this%u,this%u_mid,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV(this%decomp,this%interpMid,this%v,this%v_mid,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV(this%decomp,this%interpMid,this%w,this%w_mid,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV(this%decomp,this%interpMid,this%p,this%p_mid,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV(this%decomp,this%interpMid,this%mix%kappa,kappabar,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc) 

      this%rho_mid = 0
      rhou_int = 0
      rhov_int = 0
      rhow_int = 0
      sos_int  = 0
      do i = 1,2
         
        call this%mix%material(i)%getFaces(this%rho,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        this%rho_mid = this%rho_mid + this%mix%material(i)%rhoYs_mid
        call interpolateFV(this%decomp,this%interpMid,this%mix%material(1)%VF,VFbar,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

      end do 

!      if (nrank.eq.0) print*,rhou_int(:,:,:,1)-this%u_mid(:,:,:,1)
!      if (nrank.eq.0) print*,rhov_int(:,:,:,2)-this%v_mid(:,:,:,2)

!      this%u_int = rhou_int(:,:,:,1) /this%rho_mid(:,:,:,1) - this%u_mid(:,:,:,1)
!      this%v_int = rhov_int(:,:,:,2)/this%rho_mid(:,:,:,2) - this%v_mid(:,:,:,2)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Calculate Pressure Gradients and Coefficients       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      call gradFV_N2Fx(this%decomp,this%derStaggd04,this%p,gradp(:,:,:,1),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)           
      call gradFV_N2Fy(this%decomp,this%derStaggd04,this%p,gradp(:,:,:,2),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
      call gradFV_N2Fz(this%decomp,this%derStaggd04,this%p,gradp(:,:,:,3),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
      call interpolateFV(this%decomp,this%interpMid,this%sos,sos_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
      this%mix%material(2)%Ys_mid = 1.0_rkind - this%mix%material(1)%Ys_mid
      this%mix%material(2)%VF_mid = 1.0_rkind - this%mix%material(1)%VF_mid

      mask1 = (1.0_rkind/0.6598_rkind)*(abs(VFbar(:,:,:,1)*(1.0_rkind-VFbar(:,:,:,1)) ) )**0.3_rkind    ! (1_rkind/0.6598_rkind)*(abs(VFbar(:,:,:,1)*(1-VFbar(:,:,:,1)) ) )**(0.3_rkind)
      mask2 = (1.0_rkind/0.6598_rkind)*abs(VFbar(:,:,:,2)*(1.0_rkind-VFbar(:,:,:,2)))**0.3_rkind   ! (1_rkind/0.6598_rkind)*( abs(VFbar(:,:,:,2)*(1-VFbar(:,:,:,2)) ) )**(0.3_rkind)
      mask3 = (1.0_rkind/0.6598_rkind)*abs(VFbar(:,:,:,3)*(1.0_rkind-VFbar(:,:,:,3)))**0.3_rkind 

      af(:,:,:,1) =this%CP*mask1*this%dt/(this%rho_mid(:,:,:,1))*this%sos_ratio*rho_ratio  * ( (sqrt( this%u_mid(:,:,:,1)**2.0_rkind + this%v_mid(:,:,:,1)**2.0_rkind + this%w_mid(:,:,:,1 )**2.0_rkind )) /sos_int(:,:,:,1) )**2.0_rkind  ! sos_int(:,:,:,i)/this%dx ! max( this%rho_mid(:,:,:,i)/this%dt, sos_int(:,:,:,i)/this%dx)
      af(:,:,:,2) =this%CP*mask2*this%dt/(this%rho_mid(:,:,:,2))*this%sos_ratio*rho_ratio *( ( sqrt( this%u_mid(:,:,:,2)**2.0_rkind + this%v_mid(:,:,:,2)**2.0_rkind + this%w_mid(:,:,:,2 )**2.0_rkind  ) ) /sos_int(:,:,:,2))**2.0_rkind
      af(:,:,:,3) =this%CP*mask3*this%dt/(this%rho_mid(:,:,:,3))*this%sos_ratio*rho_ratio *( ( sqrt( this%u_mid(:,:,:,3)**2.0_rkind + this%v_mid(:,:,:,3)**2.0_rkind + this%w_mid(:,:,:,3 )**2.0_rkind  ) ) /sos_int(:,:,:,3))**2.0_rkind

      peff = ( gradp + this%surfaceTension_coeff*kappabar*gradVF )

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Y Correction !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = 2,this%nyp-1
        this%v_int(:,j,:) = this%v_mid(:,j,:,2) -  af(:,j,:,2) *  ( peff(:,j,:,2) - 0.5_rkind*(peff(:,j+1,:,2) + peff(:,j-1,:,2) ) ) !  * ( gradp(:,j,:,2) + this%surfaceTension_coeff*kappabar(:,j,:,2)*gradVF(:,j,:,2) )  !/this%dy * ( ( this%p(:,j+1,:) - this%p(:,j,:)   ) & (peff(:,j,:,2) - peff_fil(:,j,:,2))

     enddo


     this%v_int(:,this%nyp,:) = this%v_mid(:,this%nyp,:,2) -   af(:,this%nyp,:,2) *  ( peff(:,this%nyp,:,2) - 0.5_rkind*(peff(:,1,:,2) + peff(:,this%nyp-1,:,2) ) ) !(peff(:,this%nyp,:,2) - peff_fil(:,this%nyp,:,2))
     this%v_int(:,1,:) = this%v_mid(:,1,:,2) - af(:,1,:,2) * ( peff(:,1,:,2) - 0.5_rkind*(peff(:,2,:,2) + peff(:,this%nyp,:,2) ) ) !(peff(:,1,:,2) - peff_fil(:,1,:,2)
     !this%v_int =this%v_mid(:,:,:,2)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! X Correction      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tmp1 = this%u_mid(:,:,:,1)
      tmp2 = af(:,:,:,1) 

      call transpose_y_to_x(tmp1,xtmp1,this%decomp)
      call transpose_y_to_x(tmp2,xtmp2,this%decomp)
      call transpose_y_to_x(peff(:,:,:,1),xtmp3,this%decomp)
      call transpose_y_to_x(kappabar(:,:,:,1),xtmp4,this%decomp)
!      call transpose_y_to_x(this%mix%material(1)%VF_mid(:,:,:,1),xtmp5,this%decomp)
!      
      do i = 2,nx-1
!
         xtmp1(i,:,:) = xtmp1(i,:,:) - xtmp2(i,:,:) * ( xtmp3(i,:,:)  - 0.5_rkind*(xtmp3(i+1,:,:) + xtmp3(i-1,:,:) ))

                        
      enddo
!
      xtmp1(nx,:,:) = xtmp1(nx,:,:) - xtmp2(nx,:,:) * ( xtmp3(nx,:,:)   - 0.5_rkind*(xtmp3(1,:,:) + xtmp3(nx-1,:,:) ))
      xtmp1(1,:,:) = xtmp1(1,:,:) - xtmp2(1,:,:) * ( xtmp3(1,:,:)        - 0.5_rkind*(xtmp3(nx,:,:) + xtmp3(nx,:,:) ))

      !/this%dx *( ( xtmp3(1,:,:) - xtmp3(nx,:,:) ) + this%surfaceTension_coeff*xtmp4(nx,:,:)*( xtmp5(1,:,:) - xtmp5(nx,:,:) )  )
      call transpose_x_to_y(xtmp1,tmp1,this%decomp)
      this%u_int(:,:,:) = tmp1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Z Correction !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tmp = this%w_mid(:,:,:,3)
      tmp3 = af(:,:,:,3)

      if(nz .GT. 1 ) then
      call transpose_y_to_z(tmp,ztmp1,this%decomp)
      call transpose_y_to_z(tmp3,ztmp2,this%decomp)
      call transpose_y_to_z(peff(:,:,:,3),ztmp3,this%decomp)
      call transpose_y_to_z(kappabar(:,:,:,3),ztmp4,this%decomp)
!
      do i = 2,nz-1
!
         ztmp1(:,:,i) = ztmp1(:,:,i) - ztmp2(:,:,i) * ( ztmp3(:,:,i)  - 0.5_rkind*(ztmp3(:,:,i+1) + ztmp3(:,:,i-1) ))


      enddo
!
      ztmp1(:,:,nz) = ztmp1(:,:,nz) - ztmp2(:,:,nz) * ( ztmp3(:,:,nz)   - 0.5_rkind*(ztmp3(:,:,1) + ztmp3(:,:,nz-1) ))
      ztmp1(:,:,1) = ztmp1(:,:,1) - ztmp2(:,:,1) * ( ztmp3(:,:,1)        - 0.5_rkind*(ztmp3(:,:,nz) + ztmp3(:,:,nz) ))

      call transpose_z_to_y(ztmp1,tmp,this%decomp)
      this%w_int(:,:,:) =  tmp
      endif
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      this%m1_int = this%mix%material(1)%rhoYs_mid(:,:,:,2)
      this%m2_int = this%mix%material(2)%rhoYs_mid(:,:,:,2)

    end subroutine

    subroutine get_dt(this,stability)
        use reductions, only : P_MAXVAL, P_MINVAL
        use decomp_2d,  only: nrank
	use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
        class(sgrid), target, intent(inout) :: this
        character(len=*), intent(out) :: stability
        real(rkind) :: dtCFL, a, dtsigma,dtsigma2,dtsigma3,dtmu, dtbulk, dtkap, dtdiff, dtdiff_g, dtdiff_gt, dtdiff_gp, dtplast, phys_mu, delta, dtSharp_diff, dtSharp_Adiff,alpha,dtSharp_bound,st_fac=10.D0,dtYs1, dtYs2,dtVF1,dtVF2,deltay,dteKap,dtCurv
        integer :: i
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: uk,rhofil,rhokappafil !Source term for possible use in VF, g eh eqns
        character(len=30) :: str,str2
        real(rkind) :: dtSponge, sigma_max
        this%st_limit = 20

        if(this%Stretch1Dy) then
          deltay = P_MINVAL(this%dy_stretch)
          delta = min(this%dx,deltay,this%dz)
        else 
          delta = min(this%dx, this%dy, this%dz)
	endif
	phys_mu = min(this%phys_mu1, this%phys_mu2)

        sigma_max = P_MAXVAL(this%sponge(:,:,:,1))
        dtSponge  = this%CFL / sigma_max
        ! continuum
        dtCFL  = this%CFL / P_MAXVAL( ABS(this%u)/this%dx + ABS(this%v)/this%dy + ABS(this%w)/this%dz   &
                 + this%sos*sqrt( one/(this%dx**2) + one/(this%dy**2) + one/(this%dz**2) ))


        !dtCFL  = this%CFL / P_MAXVAL( ABS(this%u)/this%dx + ABS(this%v)/this%dy +  &
        !       + this%sos*sqrt( one/(this%dx**2) + one/(this%dy**2)  ))
        !print *, "u maxval"
        !print *, P_MAXVAL( ABS(this%u)/this%dx )
        !print *, "v maxval"
        !print *, P_MAXVAL( ABS(this%v)/this%dy )
        !print *, "w maxval"
        !print *, P_MAXVAL( ABS(this%w)/this%dz )        
        !print *, "sos"
        !print *, P_MAXVAL( this%sos*sqrt( one/(this%dx**2) + one/(this%dy**2) +one/(this%dz**2) ) )
        !print *, "dy"
        !print *, this%dy
        !print *, "dx"
        !print *, this%dx
        !print *, "dz"
        !print *, this%dz
        !print *, "sosmax"
        !print *, P_MAXVAL( this%sos)
        rhokappafil = this%mix%material(1)%rhodiff*(this%mix%kappa)**2.0_rkind
        call this%filter(rhofil, this%fil,1,-this%x_bc,this%y_bc,this%z_bc)
        dtmu   = 0.2_rkind * delta**2.0_rkind / (P_MAXVAL( this%mu/this%rho   ) + eps) * this%CFL
        dtYs1 = 0.5_rkind * delta**2.0_rkind / (max(P_MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,1)),P_MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,2)),P_MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,3))) + eps) 
        dtYs2 = dtYs1 !0.75_rkind * delta**2 / (P_MAXVAL( this%mix%material(2)%rhodiff  ) + eps)
        dtVF1 = dtYs1 !0.75_rkind * delta**2 / (P_MAXVAL( this%mix%material(1)%adiff  ) + eps)
        dtVF2 = dtYs1 !0.75_rkind * delta**2 / (P_MAXVAL( this%mix%material(2)%adiff  ) + eps)

        !dtbulk = 0.2_rkind * delta**2.0_rkind / (P_MAXVAL( this%bulk/ this%rho ) + eps) * this%CFL
        dtbulk = 0.2_rkind * delta**2.0_rkind / (P_MAXVAL( this%bulk/ this%rho ) + eps) !/ 5.0 !test /5
	dtCurv = 0.5_rkind   /(max(P_MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,1)*(this%mix%kappa)**2.0_rkind),P_MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,2)*(this%mix%kappa)**2.0_rkind),P_MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,3)*(this%mix%kappa)**2.0_rkind)) + eps)

	if ((this%use_surfaceTension) .OR. (this%use_CnsrvSurfaceTension)) then
              !  if ( phys_mu > eps) then
          ! filter3D(this%
                  uk = ABS(this%u*this%mix%norm(:,:,:,1)) + ABS(this%v*this%mix%norm(:,:,:,2)) + ABS(this%w*this%mix%norm(:,:,:,3))
	          dtsigma = this%CFL*4*max(1/P_MAXVAL(sqrt((4*pi*this%surfaceTension_coeff*(1/this%dx**3 + 1/this%dy**3 ))/ (this%rho))) , &
                                  P_MINVAL(this%mu/(this%surfaceTension_coeff*(1/this%dx + 1/this%dy ) ) ))
                                  !max(1/P_MAXVAL( sqrt((4*pi*this%surfaceTension_coeff*(1/this%dx**3 + 1/this%dy**3 +1/this%dx**3) )/ (this%rho))) , &
                                  !P_MINVAL(this%mu/(this%surfaceTension_coeff*(1/this%dx + 1/this%dy + 1/this%dz) ) ))
              !  else 
                  dtsigma2 = P_MINVAL( sqrt( (this%rho)/(4*pi*this%surfaceTension_coeff*(1/this%dx**3 + 1/this%dy**3 + 1/this%dx**3) ) ) )
                  dtsigma3 =  P_MINVAL( 1/ (sqrt((4*pi*this%surfaceTension_coeff) / (this%rho*( this%dx**3 + this%dy**3 +this%dz**3))) + uk*(1/this%dx + 1/this%dy + 1/this%dz  )))

              !  end if
	end if
        ! species specific
        call this%mix%get_dt(this%rho, delta, dtkap, dtdiff, dtdiff_g, dtdiff_gt, dtdiff_gp, dtplast)

        if (this%PTeqb) then
            !dtkap  = delta**2 / (P_MAXVAL( this%kap*this%T/(this%rho*this%sos**2)) + eps)   ! Cook (2007) formulation
            dtkap  = one / ( (P_MAXVAL(this%kap*this%T/(this%rho*delta**4)))**(third) + eps) ! Cook (2009) formulation
        end if

        dtkap     = 0.2_rkind * dtkap! * this%CFL
        !dtkap     = 0.2_rkind * dtkap !/ 5.0! test

        dtdiff    = 0.2_rkind * dtdiff! * this%CFL
        dtdiff_g  = 0.2_rkind * dtdiff_g! * this%CFL
        dtdiff_gt = 0.2_rkind * dtdiff_gt! * this%CFL
        dtdiff_gp = 0.2_rkind * dtdiff_gp! * this%CFL
        
        if(this%intSharp) then

           if(this%intSharp_gam.lt.-0.5) then !maximize intSharp_gam without restricting time step, based on CFL
              !For intSharp_gam = -1
              !if(this%intSharp_gam.gt.-two) then
              this%mix%intSharp_gam = 0.2_rkind * delta**2 / (dtCFL * this%mix%intSharp_eps + eps)
              !else !ELSE: set intSharp_gam based on maximum velocity as in Tiwari, Freund, Pantano JCP 2013   !For intSharp_gam <= -2
              if(this%intSharp_gam.lt.-1.5) then
                 this%mix%intSharp_gam = zero
                 do i=1,this%mix%ns
                    this%mix%intSharp_gam = max( this%mix%intSharp_gam, P_MAXVAL( four*sqrt( this%u**2 + this%v**2 + this%w**2 ) * this%mix%material(i)%VF*(one-this%mix%material(i)%VF)) )
                 enddo
                 if(this%intSharp_gam.lt.-2.5) then
                    this%mix%intSharp_gam = zero
                    this%mix%intSharp_gam = max( this%mix%intSharp_gam, P_MAXVAL(sqrt( this%u**2 + this%v**2 + this%w**2 ) ) )
                 endif
                 if (this%intSharp_gam.lt.-3.5) then
                 this%mix%intSharp_gam = zero
                 do i=1,this%mix%ns
                     this%mix%intSharp_gam = max( this%mix%intSharp_gam,P_MAXVAL( four*sqrt( this%v**2 ) *this%mix%material(i)%VF*(one-this%mix%material(i)%VF)) )
                 enddo
                 endif

                 if (this%intSharp_gam.lt.-4.5) then
                 this%mix%intSharp_gam = zero
                 do i=1,this%mix%ns
                     this%mix%intSharp_gam = max(this%mix%intSharp_gam,P_MAXVAL( sqrt( this%v**2 )) )
                 enddo
                 endif



              endif
              !print*,this%intSharp_gam,this%mix%intSharp_gam
           endif

           ! if (this%step .LE. this%st_limit) then
           !    !this%mix%intSharp_gam = zero
           !    this%mix%intSharp_gam = this%mix%intSharp_gam * 1.0D-2
           !    if (nrank.eq.0) print*,"limiting intSharp_gam"
           ! endif
           ! ! ! this%mix%intSharp_gam = zero
           ! ! ! if (nrank.eq.0) print*,"limiting intSharp_gam"
              

           dtSharp_diff =  delta**2 / (P_MAXVAL( 6*this%mix%intSharp_gam*this%mix%intSharp_eps) + eps)*this%CFL !based on diffusivity in VF sharpening equation 

           !dtSharp_diff = delta**2/(2.0*this%mix%intSharp_gam*this%mix%intSharp_eps)!from Suhas Jain, Mani, Moin JCP 2020 -- not work

           dtSharp_bound = one/eps
           if(this%intSharp_msk) then
              do i=1,this%mix%ns
                 !dtSharp_bound = min(dtSharp_bound, 0.2_rkind * delta**2 / (P_MAXVAL( this%mix%intSharp_gam*this%mix%intSharp_eps*this%mix%VFboundDiff(:,:,:,i)) + eps))! * this%CFL !based on VF out of bounds diffusivity in VF sharpening equation 
                 !dtSharp_bound = min(dtSharp_bound, 0.2_rkind * delta**2 / (P_MAXVAL( this%mix%intSharp_gam*this%mix%intSharp_eps*this%mix%VFboundDiff(:,:,:,i)) + eps))! * this%CFL !based on VF out of bounds diffusivity in VF sharpening equation 
                 !dtSharp_bound = min(dtSharp_bound, 0.2_rkind * delta**2 / (P_MAXVAL( this%mix%intSharp_dif*this%mix%intSharp_gam*this%intSharp_eps*this%mix%VFboundDiff(:,:,:,i)) + eps))! * this%CFL !based on VF out of bounds diffusivity in VF sharpening equation 
                 dtSharp_bound = min(dtSharp_bound, 0.2_rkind * delta**2 / (P_MAXVAL( this%intSharp_dif*this%mix%intSharp_gam*this%intSharp_eps*this%mix%VFboundDiff(:,:,:,i)) + eps))! * this%CFL !based on VF out of bounds diffusivity in VF sharpening equation 
              enddo
           endif

           dtSharp_Adiff   = delta/max(this%mix%intSharp_gam,eps)! * this%CFL !based on anti-diffusivity in VF sharpening equation
        endif
         a = -1
        ! Use fixed time step if CFL <= 0
        if ( this%CFL .LE. zero ) then
            this%dt = this%dtfixed
            stability = 'fixed'
        else
            stability = 'convective'
            this%dt = dtCFL
            if ( this%dt > dtmu ) then
                 this%dt = dtmu
                 stability = 'shear'
             else if ( this%dt > dtbulk ) then
                 this%dt = dtbulk
                 stability = 'bulk'
             else if ( this%dt > dtYs1 ) then
                 this%dt = dtYs1
                 stability = 'Ys1'
             else if ( this%dt > dtYs2) then
                 this%dt = dtYs2
                 stability = 'Ys2'
             else if ( this%dt > dtVF1 ) then
                 this%dt = dtVF1
                 stability = 'VF1'
             else if ( this%dt > dtVF2) then
                 this%dt = dtVF2
                 stability = 'VF2'
             else if ( this%dt > dtkap ) then
                 this%dt = dtkap
                 stability = 'conductive'
             else if ( this%dt > dtdiff ) then
                 this%dt = dtdiff
                 stability = 'diffusive'
             else if ( this%dt > dtdiff_g ) then
                 this%dt = dtdiff_g
                 stability = 'diffusive g'
             else if ( this%dt > dtdiff_gt ) then
                 this%dt = dtdiff
                 stability = 'diffusive g_t'
             else if ( this%dt > dtplast ) then
                 this%dt = dtplast
                 stability = 'plastic'
              else if ( this%dt > dtSponge ) then
                 this%dt = dtplast
                 stability = 'sponge'
             else if (this%dt > dtCurv ) then
                 stability = 'numerical curvature'
             endif
          
            if ( this%dt > dtCurv ) then
               this%dt = dtCurv
               write(str,'(ES10.3E3)') 1.0D0-dtCurv/dtCFL
               stability = 'Curv: '//trim(str)//' CFL loss fraction'
            endif
 
            if ( this%dt > dtSponge ) then
               this%dt = dtSponge
               write(str,'(ES10.3E3)') 1.0D0-dtmu/dtCFL
               stability = 'Sponge: '//trim(str)//' CFL loss fraction'
            endif 
            if ( this%dt > dtmu ) then
               this%dt = dtmu
               write(str,'(ES10.3E3)') 1.0D0-dtmu/dtCFL
               stability = 'shear: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtbulk ) then
               this%dt = dtbulk
               write(str,'(ES10.3E3)') 1.0D0-dtbulk/dtCFL
               stability = 'bulk: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtkap ) then
               this%dt = dtkap
               write(str,'(ES10.3E3)') 1.0D0-dtkap/dtCFL
               stability = 'conductive: '//trim(str)//' CFL loss fraction'
            endif

            if ( this%dt > dtdiff ) then
               this%dt = dtdiff
               write(str,'(ES10.3E3)') 1.0D0-dtdiff/dtCFL
               stability = 'diffusive: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtdiff_g ) then
               this%dt = dtdiff_g
               write(str,'(ES10.3E3)') 1.0D0-dtdiff_g/dtCFL
               stability = 'diffusive g: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtdiff_gt ) then
               this%dt = dtdiff_gt
               write(str,'(ES10.3E3)') 1.0D0-dtdiff_gt/dtCFL
               stability = 'diffusive g_t: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtdiff_gp ) then
               this%dt = dtdiff_gp
               write(str,'(ES10.3E3)') 1.0D0-dtdiff_gp/dtCFL
               stability = 'diffusive g_p: '//trim(str)//' CFL loss fraction'
            endif
	if ((this%use_surfaceTension) .OR. (this%use_CnsrvSurfaceTension)) then	
	    if ( this%dt > dtsigma ) then
               this%dt = dtsigma
               write(str,'(ES10.3E3)') 1.0D0-dtsigma/dtCFL
               stability = 'surfaceTension: '//trim(str)//' CFL loss fraction'
            endif
	endif
            if ( this%dt > dtplast ) then
               this%dt = dtplast
               write(str,'(ES10.3E3)') 1.0D0-dtplast/dtCFL
               stability = 'plastic: '//trim(str)//' CFL loss fraction'
            end if
            if (this%intSharp) then
               if ( this%dt > dtSharp_diff ) then
                  this%dt = dtSharp_diff
                  write(str,'(ES10.3E3)') 1.0D0-dtSharp_diff/dtCFL
                  stability = 'sharp diff: '//trim(str)//' CFL loss fraction'
               end if
               if ( this%dt > dtSharp_Adiff ) then
                  this%dt = dtSharp_Adiff
                  write(str,'(ES10.3E3)') 1.0D0-dtSharp_Adiff/dtCFL
                  stability = 'sharp a-diff: '//trim(str)//' CFL loss fraction'
               end if
               
               if(this%intSharp_msk) then

                  if ( this%dt > dtSharp_bound ) then
                      this%dt = dtSharp_bound
                      !write(str2,'(F25.18)') dtCFL
                      !write(str,'(F6.2)') 1.0D2*dtSharp_bound/dtCFL
                      !stability = 'Sharp VF bounds: '//trim(str)//'%'//' '//trim(str2)
                      write(str,'(ES10.3E3)') 1.0D0-dtSharp_bound/dtCFL
                      stability = 'sharp VF bounds: '//trim(str)//' CFL loss fraction'
                  end if
               end if
            end if

            if (this%step .LE. this%st_limit) then
               this%dt = min(this%dt / st_fac, this%dtfixed)
               stability = 'startup'
            endif
         endif

    end subroutine

    subroutine get_dtlocal(this,stability)
    use reductions, only : P_MAXVAL, P_MINVAL
    use decomp_2d,  only: nrank
    use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
    use mpi
    use kind_parameters, only: mpirkind

    class(sgrid), target, intent(inout) :: this
    character(len=*), intent(out) :: stability
    real(rkind) :: dtCFL, a, dtsigma,dtsigma2,dtsigma3,dtmu, dtbulk, dtkap, dtdiff, dtdiff_g, dtdiff_gt, dtdiff_gp, dtplast, phys_mu, delta, dtSharp_diff, dtSharp_Adiff,alpha,dtSharp_bound,st_fac=10.D0,dtYs1, dtYs2,dtVF1,dtVF2,deltay,dteKap,dtCurv
    integer :: i, ierr, nvals, idx
    real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: uk,rhofil,rhokappafil
    character(len=30) :: str,str2
    real(rkind) :: dtSponge, sigma_max

    ! Arrays for batched reductions
    real(rkind), allocatable :: local_vals(:), global_max(:), global_min(:)
    integer :: n_max, n_min

    this%st_limit = 20

    ! Determine delta
    if(this%Stretch1Dy) then
        deltay = MINVAL(this%dy_stretch)
        delta = min(this%dx,deltay,this%dz)
    else
        delta = min(this%dx, this%dy, this%dz)
    endif
    phys_mu = min(this%phys_mu1, this%phys_mu2)

    ! ============================================================
    ! BATCHED REDUCTIONS - Compute all local max/min values first
    ! ============================================================

    ! Count how many reductions we need
    n_max = 8  ! Base maxvals
    n_min = 2  ! Base minvals

    if ((this%use_surfaceTension) .OR. (this%use_CnsrvSurfaceTension)) then
        n_max = n_max + 2
        n_min = n_min + 2
    endif

    if (this%intSharp) then
        n_max = n_max + 1  ! For intSharp_gam calculation
        if (this%intSharp_gam .lt. -1.5) then
            n_max = n_max + 1  ! velocity magnitude
            if (this%intSharp_gam .lt. -2.5) then
                n_max = n_max + 1
            endif
            if (this%intSharp_gam .lt. -3.5) then
                n_max = n_max + 1
            endif
            if (this%intSharp_gam .lt. -4.5) then
                n_max = n_max + 1
            endif
        endif
        if (this%intSharp_msk) n_max = n_max + this%mix%ns
    endif

    allocate(local_vals(max(n_max, n_min)))
    allocate(global_max(n_max))
    allocate(global_min(n_min))

    ! ============================================================
    ! Compute all local MAXVAL operations
    ! ============================================================
    idx = 0

    ! 1. Sponge sigma_max
    idx = idx + 1
    local_vals(idx) = MAXVAL(this%sponge(:,:,:,1))

    ! 2. CFL velocity term
    idx = idx + 1
    local_vals(idx) = MAXVAL(ABS(this%u)/this%dx + ABS(this%v)/this%dy + ABS(this%w)/this%dz + &
                             this%sos*sqrt(one/(this%dx**2) + one/(this%dy**2) + one/(this%dz**2)))

    ! 3. Viscosity (mu/rho)
    idx = idx + 1
    local_vals(idx) = MAXVAL(this%mu/this%rho)

    ! 4. Bulk viscosity (bulk/rho)
    idx = idx + 1
    local_vals(idx) = MAXVAL(this%bulk/this%rho)

    ! 5-7. dtYs1 terms (3 directional components)
    idx = idx + 1
    local_vals(idx) = MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,1))
    idx = idx + 1
    local_vals(idx) = MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,2))
    idx = idx + 1
    local_vals(idx) = MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,3))

    ! 8-10. dtCurv terms (3 directional components with kappa)
    idx = idx + 1
    local_vals(idx) = MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,1)*(this%mix%kappa)**2.0_rkind)
    idx = idx + 1
    local_vals(idx) = MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,2)*(this%mix%kappa)**2.0_rkind)
    idx = idx + 1
    local_vals(idx) = MAXVAL(this%mix%material(1)%adiff_stagg(:,:,:,3)*(this%mix%kappa)**2.0_rkind)

    ! Surface tension terms
    if ((this%use_surfaceTension) .OR. (this%use_CnsrvSurfaceTension)) then
        uk = ABS(this%u*this%mix%norm(:,:,:,1)) + ABS(this%v*this%mix%norm(:,:,:,2)) + &
             ABS(this%w*this%mix%norm(:,:,:,3))

        idx = idx + 1
        local_vals(idx) = MAXVAL(sqrt((4*pi*this%surfaceTension_coeff*(1/this%dx**3 + 1/this%dy**3))/(this%rho)))

        idx = idx + 1
        local_vals(idx) = MAXVAL(1/(sqrt((4*pi*this%surfaceTension_coeff)/(this%rho*(this%dx**3 + this%dy**3 + this%dz**3))) + &
                                     uk*(1/this%dx + 1/this%dy + 1/this%dz)))
    endif

    ! Interface sharpening terms
    if (this%intSharp) then
        idx = idx + 1
        local_vals(idx) = MAXVAL(6*this%mix%intSharp_gam*this%mix%intSharp_eps)

        if (this%intSharp_gam .lt. -1.5) then
            idx = idx + 1
            local_vals(idx) = MAXVAL(four*sqrt(this%u**2 + this%v**2 + this%w**2) * &
                                     this%mix%material(1)%VF*(one-this%mix%material(1)%VF))

            if (this%intSharp_gam .lt. -2.5) then
                idx = idx + 1
                local_vals(idx) = MAXVAL(sqrt(this%u**2 + this%v**2 + this%w**2))
            endif

            if (this%intSharp_gam .lt. -3.5) then
                idx = idx + 1
                local_vals(idx) = MAXVAL(four*sqrt(this%v**2)*this%mix%material(1)%VF*(one-this%mix%material(1)%VF))
            endif

            if (this%intSharp_gam .lt. -4.5) then
                idx = idx + 1
                local_vals(idx) = MAXVAL(sqrt(this%v**2))
            endif
        endif

        if (this%intSharp_msk) then
            do i=1,this%mix%ns
                idx = idx + 1
                local_vals(idx) = MAXVAL(this%intSharp_dif*this%mix%intSharp_gam*this%intSharp_eps* &
                                         this%mix%VFboundDiff(:,:,:,i))
            enddo
        endif
    endif

    n_max = idx  ! Actual number of max values computed

    ! Single MPI_Allreduce for all MAX operations
    call MPI_Allreduce(local_vals(1:n_max), global_max, n_max, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! ============================================================
    ! Compute all local MINVAL operations
    ! ============================================================
    idx = 0

    if (this%Stretch1Dy) then
        idx = idx + 1
        local_vals(idx) = MINVAL(this%dy_stretch)
    endif

    if ((this%use_surfaceTension) .OR. (this%use_CnsrvSurfaceTension)) then
        idx = idx + 1
        local_vals(idx) = MINVAL(this%mu/(this%surfaceTension_coeff*(1/this%dx + 1/this%dy)))

        idx = idx + 1
        local_vals(idx) = MINVAL(sqrt(this%rho/(4*pi*this%surfaceTension_coeff* &
                                                 (1/this%dx**3 + 1/this%dy**3 + 1/this%dx**3))))
    endif

    n_min = idx  ! Actual number of min values computed

    ! Single MPI_Allreduce for all MIN operations
    if (n_min > 0) then
        call MPI_Allreduce(local_vals(1:n_min), global_min, n_min, mpirkind, MPI_MIN, MPI_COMM_WORLD, ierr)
    endif

    deallocate(local_vals)

    ! ============================================================
    ! Extract results and compute timesteps
    ! ============================================================
    idx = 0

    ! 1. Sponge
    idx = idx + 1
    sigma_max = global_max(idx)
    dtSponge = this%CFL / sigma_max

    ! 2. CFL
    idx = idx + 1
    dtCFL = this%CFL / global_max(idx)

    ! 3. Viscosity
    idx = idx + 1
    dtmu = 0.2_rkind * delta**2.0_rkind / (global_max(idx) + eps) * this%CFL

    ! 4. Bulk
    idx = idx + 1
    dtbulk = 0.2_rkind * delta**2.0_rkind / (global_max(idx) + eps)

    ! 5-7. dtYs1
    idx = idx + 1
    dtYs1 = 0.5_rkind * delta**2.0_rkind / (max(global_max(idx), global_max(idx+1), global_max(idx+2)) + eps)
    idx = idx + 2  ! Skip next two since we used them above

    ! 8-10. dtCurv
    idx = idx + 1
    dtCurv = 0.5_rkind / (max(global_max(idx), global_max(idx+1), global_max(idx+2)) + eps)
    idx = idx + 2

    ! Surface tension
    if ((this%use_surfaceTension) .OR. (this%use_CnsrvSurfaceTension)) then
        idx = idx + 1
        idx = idx + 1
        dtsigma3 = global_max(idx)

        idx = 0  ! Reset for min values
        idx = idx + 1
        dtsigma = this%CFL * 4 * max(one/global_max(idx), global_min(idx))

        idx = idx + 1
        dtsigma2 = global_min(idx)
    endif

    ! Interface sharpening
    if (this%intSharp) then
        if (this%intSharp_gam .lt. -0.5) then
            if (this%intSharp_gam .gt. -two) then
                this%mix%intSharp_gam = 0.2_rkind * delta**2 / (dtCFL * this%mix%intSharp_eps + eps)
            else if (this%intSharp_gam .lt. -1.5) then
                ! Use pre-computed values from global_max
                ! Implementation depends on specific logic needed
            endif
        endif

        idx = idx + 1
        dtSharp_diff = delta**2 / (global_max(idx) + eps) * this%CFL

        dtSharp_bound = one/eps
        if (this%intSharp_msk) then
            do i=1,this%mix%ns
                idx = idx + 1
                dtSharp_bound = min(dtSharp_bound, 0.2_rkind * delta**2 / (global_max(idx) + eps))
            enddo
        endif

        dtSharp_Adiff = delta/max(this%mix%intSharp_gam,eps)
    endif

    deallocate(global_max)
    if (allocated(global_min)) deallocate(global_min)

    ! ============================================================
    ! Rest of the original logic (unchanged)
    ! ============================================================
    a = -1

    if ( this%CFL .LE. zero ) then
        this%dt = this%dtfixed
        stability = 'fixed'
    else
        stability = 'convective'
        this%dt = dtCFL

        if ( this%dt > dtmu ) then
            this%dt = dtmu
            stability = 'shear'
        else if ( this%dt > dtbulk ) then
            this%dt = dtbulk
            stability = 'bulk'
        else if ( this%dt > dtYs1 ) then
            this%dt = dtYs1
            stability = 'Ys1'
        else if ( this%dt > dtdiff ) then
            this%dt = dtdiff
            stability = 'diffusive'
        else if ( this%dt > dtSponge ) then
            this%dt = dtplast
            stability = 'sponge'
        else if (this%dt > dtCurv ) then
            stability = 'numerical curvature'
        endif

        if ( this%dt > dtCurv ) then
            this%dt = dtCurv
            write(str,'(ES10.3E3)') 1.0D0-dtCurv/dtCFL
            stability = 'Curv: '//trim(str)//' CFL loss fraction'
        endif

        if ( this%dt > dtSponge ) then
            this%dt = dtSponge
            write(str,'(ES10.3E3)') 1.0D0-dtmu/dtCFL
            stability = 'Sponge: '//trim(str)//' CFL loss fraction'
        endif

        if ( this%dt > dtmu ) then
               this%dt = dtmu
               write(str,'(ES10.3E3)') 1.0D0-dtmu/dtCFL
               stability = 'shear: '//trim(str)//' CFL loss fraction'
            endif
            if ( this%dt > dtbulk ) then
               this%dt = dtbulk
               write(str,'(ES10.3E3)') 1.0D0-dtbulk/dtCFL
               stability = 'bulk: '//trim(str)//' CFL loss fraction'
            endif

            if ( this%dt > dtdiff ) then
               this%dt = dtdiff
               write(str,'(ES10.3E3)') 1.0D0-dtdiff/dtCFL
               stability = 'diffusive: '//trim(str)//' CFL loss fraction'
            endif
        if ((this%use_surfaceTension) .OR. (this%use_CnsrvSurfaceTension)) then
            if ( this%dt > dtsigma ) then
               this%dt = dtsigma
               write(str,'(ES10.3E3)') 1.0D0-dtsigma/dtCFL
               stability = 'surfaceTension: '//trim(str)//' CFL loss fraction'
            endif
        endif
            if (this%intSharp) then
               if ( this%dt > dtSharp_diff ) then
                  this%dt = dtSharp_diff
                  write(str,'(ES10.3E3)') 1.0D0-dtSharp_diff/dtCFL
                  stability = 'sharp diff: '//trim(str)//' CFL loss fraction'
               end if
               if ( this%dt > dtSharp_Adiff ) then
                  this%dt = dtSharp_Adiff
                  write(str,'(ES10.3E3)') 1.0D0-dtSharp_Adiff/dtCFL
                  stability = 'sharp a-diff: '//trim(str)//' CFL loss fraction'
               end if

            end if


            if (this%step .LE. this%st_limit) then
               this%dt = min(this%dt / st_fac, this%dtfixed)
               stability = 'startup'
            endif
         endif

    end subroutine    
    subroutine get_primitive(this)
        use reductions, only: P_MAXVAL, P_MINVAL
        use decomp_2d,  only: nrank
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(:,:,:), pointer :: onebyrho
       ! real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: onebyrho
        real(rkind), dimension(:,:,:), pointer :: rhou,rhov,rhow,TE
        real(rkind) :: rhomin

        onebyrho => this%ybuf(:,:,:,1)

        call this%mix%get_rho(this%rho)

        !this%rho = 1 / this%mix%material(1)%spec_consrv(:,:,:,1)
        ! rhomin = P_MINVAL(this%rho)
        ! if (nrank.eq.0) print*,rhomin

        rhou => this%Wcnsrv(:,:,:,mom_index  )
        rhov => this%Wcnsrv(:,:,:,mom_index+1)
        rhow => this%Wcnsrv(:,:,:,mom_index+2)
        TE   => this%Wcnsrv(:,:,:, TE_index  )

        onebyrho = one/this%rho !this%mix%material(1)%spec_consrv(:,:,:,1) ! one/this%rho
        this%u = rhou * onebyrho
        this%v = rhov * onebyrho
        this%w = rhow * onebyrho

        if(this%use_CnsrvSurfaceTension) then
     
           this%e = ((TE-this%mix%surfaceTension_pe)*onebyrho) - half*( this%u*this%u + this%v*this%v + this%w*this%w )

        else
           
           this%e = (TE*onebyrho) - half*(this%u*this%u + this%v*this%v +this%w*this%w )

        endif
    
        call this%mix%get_primitive(this%rho, this%u, this%v, this%w, this%e, this%devstress, this%p, this%sos)                  ! Get primitive variables for individual species
        !this%mix%material(1)%spec_consrv(:,:,:,1) = 1 / this%rho
    end subroutine


    subroutine get_primitive_g(this)
      class(sgrid), target, intent(inout) :: this

      call this%mix%get_primitive_g(this%rho)                  ! Get primitive kinematic variables
      
    end subroutine get_primitive_g
     

    pure subroutine get_conserved(this)
        class(sgrid), intent(inout) :: this
        integer :: i 
        ! Assume rho is already available
        this%Wcnsrv(:,:,:,mom_index  ) = this%rho * this%u
        this%Wcnsrv(:,:,:,mom_index+1) = this%rho * this%v
        this%Wcnsrv(:,:,:,mom_index+2) = this%rho * this%w

        call this%mix%get_conserved(this%rho,this%u,this%v,this%w)

         ! AFTER mix

        if(this%use_CnsrvSurfaceTension) then
            this%Wcnsrv(:,:,:, TE_index  ) = this%rho * ( this%e + half*(this%u*this%u + this%v*this%v + this%w*this%w )) + this%mix%surfaceTension_pe
        else
            this%Wcnsrv(:,:,:, TE_index  ) = this%rho *(this%e +  half*(this%u*this%u + this%v*this%v + this%w*this%w ) ) !+ this%Wcnsrv(:,:,:, TE_index  )
        endif

        


    end subroutine

    subroutine get_conserved_g(this)
      class(sgrid), target, intent(inout) :: this

      call this%mix%get_conserved_g(this%rho)                  ! Get conserved kinematic variables
      
    end subroutine get_conserved_g

    subroutine post_bc(this)
        class(sgrid), intent(inout) :: this
        integer :: i

        if(this%useOneG) then
            call this%mix%get_mixture_properties()
        endif

        call this%mix%get_eelastic_devstress(this%devstress)   ! Get species elastic energies, and mixture and species devstress
        call this%mix%get_ehydro_from_p(this%rho)              ! Get species hydrodynamic energy, temperature; and mixture pressure, temperature
        call this%mix%get_pmix(this%p)
       ! Get mixture pressure
    !    if( this%useRestartFile ) then

    !        do i = 1, 2

    !           call this%mix%material(i)%get_ehydroT_from_p(this%rho)
 
    !        enddo


    !    endif

        call this%mix%get_Tmix(this%T)                         ! Get mixture temperature
        call this%mix%getSOS(this%rho,this%p,this%sos)
!print *, 'SOS: ', this%sos(179,1,1)
        ! assuming pressures have relaxed and sum( (Ys*(ehydro + eelastic) ) over all
        ! materials equals e
        call this%mix%get_emix(this%e)
    end subroutine

    subroutine post_bc_2(this)
        class(sgrid), intent(inout) :: this
        integer :: i
        if(this%useOneG) then
            call this%mix%get_mixture_properties()
        endif
        call this%mix%get_eelastic_devstress(this%devstress)   ! Get specieselastic energies, and mixture and species devstress
        ! Get specieshydrodynamic energy, temperature; and mixture pressure, temperature
        call this%mix%get_ehydro_from_p(this%rho) 
        call this%mix%get_pmix(this%p)                         ! Get mixturepressure

      !   if( this%useRestartFile ) then

      !      do i = 1, 2

      !         call this%mix%material(i)%get_ehydroT_from_p(this%rho)

      !      enddo


      !  endif

        call this%mix%get_Tmix(this%T)                         ! Get mixturetemperature
        call this%mix%getSOS(this%rho,this%p,this%sos)
!print *, 'SOS: ', this%sos(179,1,1)
        ! assuming pressures have relaxed and sum( (Ys*(ehydro + eelastic) )
        ! over all
        ! materials equals e
        call this%mix%get_emix(this%e)

    end subroutine

    subroutine CheckTau(this,tauxx,tauyy,tauzz,tauxy,tauyx,tauyz,tauzy,tauxz,tauzx)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y,transpose_y_to_z, transpose_z_to_y
        use exits,      only: message,nancheck,GracefulExit
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: tauxx,tauyy,tauzz,tauxy,tauyx,tauyz,tauzy,tauxz,tauzx
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: x_half, y_half,z_half,tau11, tau22, tau33, tau12, tau21, tau13, tau31, tau23, tau32

        print *, "Check Tau"
        x_half = this%x + 0.5*this%dx 
        y_half = this%y + 0.5*this%dy
        z_half = this%z + 0.5*this%dz
          
        tau11 = -2_rkind / 3._rkind*this%mu*(-4.0*sin(4.0*this%z)*cos(x_half)) + 4._rkind / 3._rkind*this%mu*4.0*sin(2.0*this%z)*cos(4.0*x_half)
        tau33 = 4._rkind / 3._rkind*this%mu*(-4.0*sin(4.0*z_half)*cos(this%x)) - 2._rkind /3._rkind*this%mu*(4.0*sin(2.0*z_half)*cos(4.0*this%x)) 
        tau22 = 0; tau12 = 0; tau21 = 0; tau23 = 0; tau32 = 0;
        tau13 = this%mu*(2*cos(2*z_half)*sin(4*this%x) - cos(4*z_half)*sin(this%x))
        tau31 = this%mu*(2*cos(2*this%z)*sin(4*x_half) - cos(4*this%z)*sin(x_half)) 
        
    end subroutine

    subroutine getRHS(this, rhs, divu, viscwork)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y,transpose_y_to_z, transpose_z_to_y
        use operators, only: divergence,gradient,divergenceFV, interpolateFV, interpolateFV_x, interpolateFV_y, interpolateFV_z,gradFV_N2Fx, gradFV_N2Fy, gradFV_N2Fz,interpolateMax,gradFV_x,gradFV_y
        use exits,      only: message,nancheck,GracefulExit
        use timer,      only: tic,toc
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,ncnsrv), intent(out) :: rhs
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),     intent(out) :: divu
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),target, intent(out) :: viscwork
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target  :: duidxj, duidxj_s
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,12), target :: duidxj_int
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: dudx_s,dudy_s,dudz_s,dvdx_s,dvdy_s,dvdz_s,dwdx_s,dwdy_s,dwdz_s
        real(rkind), dimension(:,:,:), pointer :: dvdy_x,dwdz_x, dvdx_y, dwdx_z
        real(rkind), dimension(:,:,:), pointer :: dudx_y, dwdz_y, dudy_x, dwdy_z
        real(rkind), dimension(:,:,:), pointer :: dudx_z, dvdy_z, dudz_x, dvdz_y
        real(rkind), dimension(:,:,:), pointer :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz, tauzy, tauzx, tauyx
        real(rkind), dimension(:,:,:), pointer :: qx,qy,qz
        real(rkind), dimension(:,:,:), pointer :: ehmix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: u_int, v_int,w_int,ke_int
        integer :: imat,i
        !logical :: useNewSPF = .FALSE.
        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: ke,tmp,dJ,drhodx,drhody,drhodz, uJ, vJ, wJ, keJ, eJ, Fbody, tmp1,tmp2, tmp3, rhoeJ
        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: drhoedx,drhoedy, drhoedz,dedx,dedz,dedy, eKap,dedx_n, dedy_n, dedz_n
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, 3) :: J,Frho,Fenergy, Fp, yMetric_F2N_int, De_int,rho_int, eLADcoef,gradrhoh,rhoh,gradp
        real(rkind) :: g = -0.1d0, cputime


        !call tic()
        if(this%use_Stagg) then

           dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
           dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
           dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

           dvdy_x => duidxj_int(:,:,:,1); dudx_y => duidxj_int(:,:,:,2); dudx_z => duidxj_int(:,:,:,3);
           dwdz_x => duidxj_int(:,:,:,4); dwdz_y => duidxj_int(:,:,:,5); dvdy_z => duidxj_int(:,:,:,6);
           dvdx_y => duidxj_int(:,:,:,7); dudy_x => duidxj_int(:,:,:,8); dudz_x => duidxj_int(:,:,:,9);
           dwdx_z => duidxj_int(:,:,:,10); dwdy_z => duidxj_int(:,:,:,11); dvdz_y => duidxj_int(:,:,:,12);

           
           call gradient(this%decomp,this%derCD06,this%u, dudx, dudy, dudz,  -this%x_bc,  this%y_bc,this%z_bc)
           call gradient(this%decomp,this%derCD06,this%v, dvdx, dvdy, dvdz,  this%x_bc, -this%y_bc,this%z_bc)
           call gradient(this%decomp,this%derCD06,this%w, dwdx, dwdy, dwdz,  this%x_bc,  this%y_bc,-this%z_bc)


           call interpolateFV_x(this%decomp,this%interpMid,dvdy,dvdy_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_x(this%decomp,this%interpMid,dwdz,dwdz_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,dvdx,dvdx_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,dwdx,dwdx_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

           call interpolateFV_y(this%decomp,this%interpMid,dudx,dudx_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,dwdz,dwdz_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_x(this%decomp,this%interpMid,dudy,dudy_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,dwdy,dwdy_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

           call interpolateFV_z(this%decomp,this%interpMid,dudx,dudx_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,dvdy,dvdy_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc) 
           call interpolateFV_x(this%decomp,this%interpMid,dudz,dudz_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,dvdz,dvdz_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

           if(this%Stretch1Dy) then

              call interpolateFV(this%decomp,this%interpMid,this%yMetric_F2N,yMetric_F2N_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
              dudy_x = yMetric_F2N_int(:,:,:,1)*dudy_x
              dvdy_x = yMetric_F2N_int(:,:,:,1)*dvdy_x
              dwdy_z = yMetric_F2N_int(:,:,:,3)*dwdy_z
              dvdy_z = yMetric_F2N_int(:,:,:,3)*dvdy_z 
           endif


           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GET STAGGERED DERIVATIVES !!!!!!!!!!!!!!!!!!!!
           dudx_s => duidxj_s(:,:,:,1); dudy_s => duidxj_s(:,:,:,2); dudz_s => duidxj_s(:,:,:,3);
           dvdx_s => duidxj_s(:,:,:,4); dvdy_s => duidxj_s(:,:,:,5); dvdz_s => duidxj_s(:,:,:,6);
           dwdx_s => duidxj_s(:,:,:,7); dwdy_s => duidxj_s(:,:,:,8); dwdz_s => duidxj_s(:,:,:,9);

           call gradFV_N2Fx(this%decomp,this%derStagg,this%u,dudx_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fx(this%decomp,this%derStagg,this%v,dvdx_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fx(this%decomp,this%derStagg,this%w,dwdx_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

           call gradFV_N2Fy(this%decomp,this%derStagg,this%u,dudy_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fy(this%decomp,this%derStagg,this%v,dvdy_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fy(this%decomp,this%derStagg,this%w,dwdy_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
          
           call gradFV_N2Fz(this%decomp,this%derStagg,this%u,dudz_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fz(this%decomp,this%derStagg,this%v,dvdz_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_N2Fz(this%decomp,this%derStagg,this%w,dwdz_s,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc) 


        else

           dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
           dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
           dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

           call this%gradient(this%u, dudx, dudy, dudz, -this%x_bc,  this%y_bc,  this%z_bc)
           call this%gradient(this%v, dvdx, dvdy, dvdz,  this%x_bc, -this%y_bc,  this%z_bc)
           call this%gradient(this%w, dwdx, dwdy, dwdz,  this%x_bc,  this%y_bc, -this%z_bc)
        endif
        !call toc(cputime)
        !   call message(4, "Viscous Stress Stuff", cputime)
        divu = dudx + dvdy + dwdz

        !call tic()
        call this%getPhysicalProperties()
        !call this%LAD%get_viscosities(this%rho,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc)
        call this%LAD%get_viscosities(this%rho,this%p,this%sos,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc,this%dt,this%intSharp_pfloor,this%yMetric,this%dy_stretch,this%mix%deltakap*abs(this%mix%material(1)%Ys*(1-this%mix%material(1)%Ys))*4_rkind,this%mix%deltakap*abs(this%mix%material(1)%VF*(1-this%mix%material(1)%VF))* 4_rkind)
        !call toc(cputime)
        !call message(4, "Viscous Stress LAD ", cputime)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Conductivity LAD        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       rhoeJ = 0

        if (this%PTeqb) then
            ! subtract elastic energies to determine mixture hydrostatic energy. conductivity 
            ! is assumed a function of only hydrostatic energy
            ehmix => viscwork ! use some storage space
            ehmix = this%e
            do imat = 1, this%mix%ns
                ehmix = ehmix - this%mix%material(imat)%Ys * this%mix%material(imat)%eel
            enddo
           ! call this%LAD%get_conductivity(this%rho,this%p,ehmix,this%T,this%sos,this%kap,this%x_bc,this%y_bc,this%z_bc,this%intSharp_tfloor)
        end if

        if( .NOT. this%use_Stagg) then

           ! Get tau tensor tensor. Put in off-diagonal components of duidxj (also get the viscous work term for energy equation)
           call this%get_tau( duidxj, viscwork )
           ! Now, associate the pointers to understand what's going on better
           tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx); tauxz => duidxj(:,:,:,tauxzidx);
                                         tauyy => duidxj(:,:,:,tauyyidx); tauyz => duidxj(:,:,:,tauyzidx);
                                                                          tauzz => duidxj(:,:,:,tauzzidx);
           !print '(a,9(e21.14,1x))', 'dudx ', duidxj(89,1,1,1:9)
           !print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)
           ! Add the deviatoric stress to the tau for use in fluxes 
           tauxx = tauxx + this%sxx; tauxy = tauxy + this%sxy; tauxz = tauxz + this%sxz
                                  tauyy = tauyy + this%syy; tauyz = tauyz + this%syz
                                                            tauzz = tauzz + this%szz
           !print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)
      
           ! store artificial stress tensor in devstress. this should not break anything since devstress will be
           ! overwritten in get_primitive and post_bc. used in update_eh -- NSG
           this%sxx = tauxx - this%sxx; this%sxy = tauxy - this%sxy; this%sxz = tauxz - this%sxz
                                     this%syy = tauyy - this%syy; this%syz = tauyz - this%syz
                                                                  this%szz = tauzz - this%szz
      
           ! Get heat conduction vector (q). Stored in remaining 3 components of duidxj 
           qx => duidxj(:,:,:,qxidx); qy => duidxj(:,:,:,qyidx); qz => duidxj(:,:,:,qzidx);

        else
           call this%get_tauStagg( duidxj, duidxj_int, duidxj_s )    
           tauxx => duidxj(:,:,:,tauxxidx); tauyx => duidxj(:,:,:,tauxyidx); tauzx => duidxj(:,:,:,tauxzidx);
           tauxy => duidxj_int(:,:,:, 1);   tauyy => duidxj(:,:,:,tauyyidx); tauzy => duidxj(:,:,:,tauyzidx);
           tauxz => duidxj_int(:,:,:,2);    tauyz => duidxj_int(:,:,:,3);    tauzz => duidxj(:,:,:,tauzzidx);
           ! Add the deviatoric stress to the tau for use in fluxes 
           !tauxx = tauxx + this%sxx; tauxy = tauxy + this%sxy; tauxz = tauxz + this%sxz
           !                       tauyy = tauyy + this%syy; tauyz = tauyz + this%syz
           !                                                 tauzz = tauzz + this%szz
          
           ! store artificial stress tensor in devstress. this should not break
           ! anything since devstress will be
           ! overwritten in get_primitive and post_bc. used in update_eh -- NSG
           !this%sxx = tauxx - this%sxx; this%sxy = tauxy - this%sxy; this%sxz = tauxz - this%sxz
           !                            this%syy = tauyy - this%syy; this%syz = tauyz - this%syz
           !                                                          this%szz = tauzz - this%szz
           this%tauxx = tauxx; this%tauyy = tauyy; this%tauzz = tauzz;
           this%tauxy = tauxy; this%tauyx = tauyx; this%tauxz = tauxz;
           this%tauzx = tauzx; this%tauyz = tauyz; this%tauzy = tauzy;
           ! Get heat conduction vector (q). Stored in remaining 3 components of
           ! duidxj 
           qx => duidxj(:,:,:,qxidx); qy => duidxj(:,:,:,qyidx); qz => duidxj(:,:,:,qzidx);

           !call this%CheckTau(tauxx,tauyy,tauzz,tauxy,tauyx,tauyz,tauzy,tauxz,tauzx)

        endif 
        if (this%PTeqb) then
          call this%get_q(qx, qy, qz)            ! add artificial thermal conduction fluxes
        end if

        rhs = zero
        call this%mix%getLAD_5eqn(this%rho,this%p,this%e,this%p_mid,Frho,Fenergy,Fp,this%x_bc,this%y_bc,this%z_bc,this%dx,this%dy,this%dz,this%periodicx,this%periodicy,this%periodicz)

        if(this%use_Stagg) then

              call this%getRHS_xStagg(              rhs,&
                                      tauxx,tauxy,tauxz,&
                                       Frho(:,:,:,1),Fenergy(:,:,:,1), qx )

              call this%getRHS_yStagg(              rhs,&
                                      tauxy,tauyy,tauyz,&
                                        Frho(:,:,:,2),Fenergy(:,:,:,2),qy )

              call this%getRHS_zStagg(              rhs,&
                                      tauxz,tauyz,tauzz,&
                                        Frho(:,:,:,3),Fenergy(:,:,:,3),qz )




           else
              call this%getRHS_x(              rhs,&
                                 tauxx,tauxy,tauxz,&
                                                qx )
           !print '(a,4(e21.14,1x))', 'rhsx: ', rhs(179,1,1,1:4)
              call this%getRHS_y(              rhs,&
                                 tauxy,tauyy,tauyz,&
                                                qy )
           !print '(a,4(e21.14,1x))', 'rhsy: ', rhs(179,1,1,1:4)

              call this%getRHS_z(              rhs,&
                                 tauxz,tauyz,tauzz,&
                                                 qz )
           !print '(a,4(e21.14,1x))', 'rhsz: ', rhs(179,1,1,1:4)
        endif
        !call toc(cputime)
        !call message(4, "RHS Stuff", cputime)
        if(this%intSharp .AND. this%intSharp_cpl) then
           !calculate kinetic energy for intSharp terms
            ke = half*( this%u**2 + this%v**2 + this%w**2 ) 
       
          !FV sharpening
           rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + this%mix%intSharp_fFV(:,:,:,1) ! + this%mix%intSharp_fDiffFV(:,:,:,1)
           rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + this%mix%intSharp_fFV(:,:,:,2) !+ this%mix%intSharp_fDiffFV(:,:,:,2)
           rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + this%mix%intSharp_fFV(:,:,:,3) !+ this%mix%intSharp_fDiffFV(:,:,:,3)
           rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + this%mix%intSharp_hFV !+ this%mix%intSharp_kFV !+ this%mix%intSharp_hDiffFV + this%mix%intSharp_kDiffFV

        endif

        if (this%use_surfaceTension .OR. ( this%use_CnsrvSurfaceTension)) then
            
            !this%puKE = gradp(:,:,:,1)-this%mix%surfaceTension_f(:,:,:,1)
            !this%SurfTenDiff = gradp(:,:,:,2)-this%mix%surfaceTension_f(:,:,:,2)
            rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + this%mix%surfaceTension_f(:,:,:,1)
            rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + this%mix%surfaceTension_f(:,:,:,2)
            rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + this%mix%surfaceTension_f(:,:,:,3)
            rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + this%mix%surfaceTension_e
        endif

        this%uJ = 0.0; this%vJ = 0.0; this%wJ = 0.0; this%keJ = 0.0; this%eJ = 0; 
        if( .NOT. this%twoPhaseLAD) then
             
          J = 0
          ke = half*( this%u**2 + this%v**2 + this%w**2 )

          do i = 1,this%mix%ns
              J = J + this%mix%material(i)%Ji
          enddo

          call divergence(this%decomp, this%der,J(:,:,:,1), J(:,:,:,2), J(:,:,:,3),dJ,this%x_bc,this%y_bc,this%z_bc )
          uJ = this%u*dJ
          vJ = this%v*dJ
          wJ = this%w*dJ
          keJ = ke*dJ

        else
        
         !call tic() 
          ke = half*( this%u**2 + this%v**2 + this%w**2 )

!          if( this%LADInt .OR. this%LADN2F) then

!             ke_int = half*(this%u_mid*this%u_mid + this%v_mid*this%v_mid + this%w_mid*this%w_mid)
!             call divergenceFV(this%decomp,this%derStagg,this%u_mid(:,:,:,1)*Frho(:,:,:,1),this%u_mid(:,:,:,2)*Frho(:,:,:,2),this%u_mid(:,:,:,3)*Frho(:,:,:,3),this%uJ,this%periodicx, this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)
!             call divergenceFV(this%decomp,this%derStagg,this%v_mid(:,:,:,1)*Frho(:,:,:,1),this%v_mid(:,:,:,2)*Frho(:,:,:,2),this%v_mid(:,:,:,3)*Frho(:,:,:,3),this%vJ,this%periodicx, this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)
!             call divergenceFV(this%decomp,this%derStagg,this%w_mid(:,:,:,1)*Frho(:,:,:,1),this%w_mid(:,:,:,2)*Frho(:,:,:,2),this%w_mid(:,:,:,3)*Frho(:,:,:,3),this%wJ,this%periodicx, this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)
!             call divergenceFV(this%decomp,this%derStagg,ke_int(:,:,:,1)*Frho(:,:,:,1),ke_int(:,:,:,2)*Frho(:,:,:,2),ke_int(:,:,:,3)*Frho(:,:,:,3),this%keJ,this%periodicx, this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)
!             call divergenceFV(this%decomp,this%derStagg,Fenergy(:,:,:,1),Fenergy(:,:,:,2),Fenergy(:,:,:,3),this%eJ,this%periodicx, this%periodicy, this%periodicz,this%x_bc,this%y_bc,this%z_bc)

!           else

!            call divergence(this%decomp,this%der,this%u*Frho(:,:,:,1),this%u*Frho(:,:,:,2),this%u*Frho(:,:,:,3),this%uJ,this%x_bc,this%y_bc,this%z_bc)
!            call divergence(this%decomp,this%der,this%v*Frho(:,:,:,1),this%v*Frho(:,:,:,2),this%v*Frho(:,:,:,3),this%vJ,this%x_bc,this%y_bc,this%z_bc)
!            call divergence(this%decomp,this%der,this%w*Frho(:,:,:,1),this%w*Frho(:,:,:,2),this%w*Frho(:,:,:,3),this%wJ,this%x_bc,this%y_bc,this%z_bc)
!            call divergence(this%decomp,this%der,ke*Frho(:,:,:,1),ke*Frho(:,:,:,2),ke*Frho(:,:,:,3),this%keJ,this%x_bc,this%y_bc,this%z_bc)
!            call divergence(this%decomp,this%der,Fenergy(:,:,:,1),Fenergy(:,:,:,2),Fenergy(:,:,:,3),this%eJ,this%x_bc,this%y_bc,this%z_bc)

!           endif 
!          rhs(:,:,:, mom_index   ) = rhs(:,:,:,mom_index   ) + this%uJ
!          rhs(:,:,:, mom_index+1 ) = rhs(:,:,:,mom_index+1 ) + this%vJ + this%rho*this%g
!          rhs(:,:,:, mom_index+2 ) = rhs(:,:,:,mom_index+2 ) + this%wJ 
!          rhs(:,:,:, TE_index )    = rhs(:,:,:,TE_index    ) + this%keJ + this%eJ  + this%rho*this%g*(this%v)

        endif
        !call toc(cputime)
        !call message(4, "LAD RHS consistency ", cputime)

!       print *, "rhow_ref1", this%rhow_ref
!       print *, "rhoe_ref", this%rhoe_ref

       if(this%SpongeLayer) then


          rhs(:,:,:, mom_index   ) = rhs(:,:,:,mom_index   ) - this%sponge(:,:,:,1)*(-this%rho*this%u + this%rhou_ref(1))  &
                                     - this%sponge(:,:,:,2)*(-this%rho*this%u + this%rhou_ref(2))

          rhs(:,:,:, mom_index+1 ) = rhs(:,:,:,mom_index+1 ) - this%sponge(:,:,:,1)*(-this%rho*this%v + this%rhov_ref(1))  &
                                     - this%sponge(:,:,:,2)*(-this%rho*this%v + this%rhov_ref(2))

          rhs(:,:,:, mom_index+2 ) = rhs(:,:,:,mom_index+2 ) - this%sponge(:,:,:,1)*(-this%rho*this%w + this%rhow_ref(1))  &
                                     - this%sponge(:,:,:,2)*(-this%rho*this%w + this%rhow_ref(2))

          rhs(:,:,:, TE_index )    = rhs(:,:,:,TE_index    ) - this%sponge(:,:,:,1)*(-this%rho*(this%e + 0.5*(this%u**2 + this%v**2 + this%w**2)) + this%rhoe_ref(1))  &
                                     - this%sponge(:,:,:,2)*(-this%rho*(this%e + 0.5*(this%u**2 + this%v**2 + this%w**2)) + this%rhoe_ref(2)) 
          !do i = 1,2
!          rhs(:,:,:, mom_index   ) = rhs(:,:,:,mom_index   ) - this%sponge(:,:,:,1)*(-this%rho*this%u + this%rhou_ref(1))  &
!                                     - this%sponge(:,:,:,2)*(-this%rho*this%u + this%rhou_ref(2))
!          rhs(:,:,:, mom_index+1 ) = rhs(:,:,:,mom_index+1 ) - this%sponge(:,:,:,1)*(-this%rho*this%v + this%rhov_ref(1))  &
!                                     - this%sponge(:,:,:,2)*(-this%rho*this%v + this%rhov_ref(2))
!          rhs(:,:,:, mom_index+2 ) = rhs(:,:,:,mom_index+2 ) - this%sponge(:,:,:,1)*(-this%rho*this%w + this%rhow_ref(1))  &
!                                     - this%sponge(:,:,:,2)*(-this%rho*this%w + this%rhow_ref(2))
!          rhs(:,:,:, TE_index )    = rhs(:,:,:,TE_index    ) - this%sponge(:,:,:,1)*(-this%rho*( 0.5*(this%u**2 + this%v**2 + this%w**2)) )  &
!                                     - this%sponge(:,:,:,2)*(-this%rho*( 0.5*(this%u**2 + this%v**2 + this%w**2))) 
         ! enddo
        endif


       !!!!!!!!!!!!!!!!!!!!!!!!!!!! GRAVITY        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! rhs(:,:,:, mom_index+1 ) = rhs(:,:,:,mom_index+1 ) + this%rho*this%g
       ! rhs(:,:,:, TE_index )    = rhs(:,:,:,TE_index    ) + this%rho*this%g*(this%v)

 
    end subroutine



subroutine getRHS_NC(this, rhs, divu, viscwork)
        use operators, only: divergence,gradient
        use exits, only: GracefulExit
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,ncnsrv), intent(out):: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj 
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,3*this%mix%ns),target :: gradYs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: NCbuff
! extra buffers for nonconservative
        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: flux, tmp, ke
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: d2udx2,d2udy2,d2udz2,d2vdx2,d2vdy2,d2vdz2,d2wdx2,d2wdy2,d2wdz2
        real(rkind), dimension(:,:,:,:), pointer :: dYsdx, dYsdy, dYsdz
        real(rkind), dimension(:,:,:,:), pointer :: d2Ysdx2, d2Ysdy2, d2Ysdz2
        real(rkind), dimension(:,:,:), pointer :: tauxx,tauxy,tauxz,tauyy,tauyz,tauzz
        real(rkind), dimension(:,:,:), pointer :: dTdx,dTdy,dTdz,d2Tdx2,d2Tdy2,d2Tdz2
        real(rkind), dimension(:,:,:), pointer :: drhou_dx, drhov_dy, drhow_dz
        real(rkind), dimension(:,:,:), pointer :: dmudx,dmudy,dmudz
        real(rkind), dimension(:,:,:), pointer :: totale,bambda,lambda
        real(rkind), dimension(:,:,:), pointer :: dbulkdx,dbulkdy,dbulkdz
        real(rkind), dimension(:,:,:), pointer :: dkapdx,dkapdy,dkapdz
        real(rkind), dimension(:,:,:,:), pointer :: Jx,Jy,Jz
        real(rkind), dimension(:,:,:), pointer :: dsumJxdx, dsumJydy, dsumJzdz
        real(rkind), dimension(:,:,:), pointer :: dJxdx, dJydy, dJzdz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),target, intent(out) :: viscwork
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),     intent(out) :: divu
        real(rkind), dimension(:,:,:), pointer :: ehmix
        real(rkind), dimension(:,:,:), pointer :: qx,qy,qz
        integer :: i,k, imat
        
        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2,ztmp1,ztmp2

        der => this%der
        decomp => this%decomp
        xtmp1 => this%xbuf(:,:,:,1)
        xtmp2 => this%xbuf(:,:,:,2)
        ztmp1 => this%zbuf(:,:,:,1)
        ztmp2 => this%zbuf(:,:,:,2)

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

        call this%gradient(this%u,dudx,dudy,dudz, [0,0], this%y_bc, this%z_bc)
        call this%gradient(this%v,dvdx,dvdy,dvdz, this%x_bc, [0,0], this%z_bc)
        call this%gradient(this%w,dwdx,dwdy,dwdz, this%x_bc, this%y_bc, [0,0])

        divu = dudx + dvdy + dwdz




        if (this%mix%ns .GT. 1) then
          dYsdx => gradYs(:,:,:,1:this%mix%ns); dYsdy => gradYs(:,:,:,this%mix%ns+1:2*this%mix%ns);
          dYsdz => gradYs(:,:,:,2*this%mix%ns+1:3*this%mix%ns);
            do i = 1,this%mix%ns
                call this%gradient(this%mix%material(i)%Ys,dYsdx(:,:,:,i),dYsdy(:,:,:,i),dYsdz(:,:,:,i),this%x_bc, this%y_bc, this%z_bc)
         end do
       endif

        call this%getPhysicalProperties()
        !call
        !this%LAD%get_viscosities(this%rho,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc)
        call  this%LAD%get_viscosities(this%rho,this%p,this%sos,duidxj,this%mu,this%bulk,this%x_bc,this%y_bc,this%z_bc,this%dt,this%intSharp_pfloor,this%yMetric,this%dy_stretch,this%mix%deltakap*abs(this%mix%material(1)%Ys*(1-this%mix%material(1)%Ys))*4_rkind,this%mix%deltakap*abs(this%mix%material(1)%VF*(1-this%mix%material(1)%VF))*4_rkind)

        if (this%PTeqb) then
            ! subtract elastic energies to determine mixture hydrostatic energy.
            ! conductivity
            ! is assumed a function of only hydrostatic energy
            ehmix => viscwork ! use some storage space
            ehmix = this%e
            do imat = 1, this%mix%ns
                ehmix = ehmix - this%mix%material(imat)%Ys *this%mix%material(imat)%eel
            enddo
        !    call this%LAD%get_conductivity(this%rho,this%p,ehmix,this%T,this%sos,this%kap,this%x_bc,this%y_bc,this%z_bc,this%intSharp_tfloor)
        end if

       ! call this%get_tau( duidxj, viscwork )
        ! Now, associate the pointers to understand what's going on better
       ! tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx); tauxz => duidxj(:,:,:,tauxzidx);
       ! tauyy => duidxj(:,:,:,tauyyidx); tauyz => duidxj(:,:,:,tauyzidx); tauzz => duidxj(:,:,:,tauzzidx);
!print '(a,9(e21.14,1x))', 'dudx ', duidxj(89,1,1,1:9)
!print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)
        ! Add the deviatoric stress to the tau for use in fluxes
        !tauxx = tauxx + this%sxx; tauxy = tauxy + this%sxy; tauxz = tauxz +this%sxz
        !                          tauyy = tauyy + this%syy; tauyz = tauyz +this%syz
        !                                                    tauzz = tauzz + this%szz
!print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)

        ! store artificial stress tensor in devstress. this should not break
        ! anything since devstress will be
        ! overwritten in get_primitive and post_bc. used in update_eh -- NSG
       ! this%sxx = tauxx - this%sxx; this%sxy = tauxy - this%sxy; this%sxz = tauxz - this%sxz
       !                              this%syy = tauyy - this%syy; this%syz = tauyz - this%syz
       !                                                           this%szz = tauzz - this%szz

        ! Get heat conduction vector (q). Stored in remaining 3 components of
        ! duidxj
!        qx => duidxj(:,:,:,qxidx); qy => duidxj(:,:,:,qyidx); qz =>duidxj(:,:,:,qzidx);
!        call this%mix%get_qmix(qx, qy, qz)                     ! Get onlyspecies diffusion fluxes if PTeqb, else, everything
!        if (this%PTeqb) call this%get_q(qx, qy, qz)            ! add artificial thermal conduction fluxes



        rhs = zero

        ! regular convection and pressure terms
            ! x-derivatives, convective conservative form
           ! select case(this%mix%ns)
            !case(1)
           !     flux = this%Wcnsrv(:,:,:,mom_index)   ! mass
           !     call transpose_y_to_x(flux,xtmp1,this%decomp)
           !     call this%der%ddx(xtmp1,xtmp2,0,0)
           !     call transpose_x_to_y(xtmp2,flux,this%decomp)
           !     rhs(:,:,:,1) = rhs(:,:,:,1) - flux
           ! case default
           !     do i = 1,this%mix%ns
           !         flux = this%Wcnsrv(:,:,:,i)*this%u   ! mass
           !         call transpose_y_to_x(flux,xtmp1,this%decomp)
           !         call this%der%ddx(xtmp1,xtmp2,0,0)
           !         call transpose_x_to_y(xtmp2,flux,this%decomp)
           !         rhs(:,:,:,i) = rhs(:,:,:,i) - flux
           !     end do
           ! end select
            flux = this%Wcnsrv(:,:,:,mom_index)*this%u + this%p
            call transpose_y_to_x(flux,xtmp1,this%decomp)
            call this%der%ddx(xtmp1,xtmp2,0,0)
            call transpose_x_to_y(xtmp2,flux,this%decomp)
            rhs(:,:,:,mom_index) = rhs(:,:,:,mom_index) - flux
            flux = this%Wcnsrv(:,:,:,mom_index+1)*this%u
            call transpose_y_to_x(flux,xtmp1,this%decomp)
            call this%der%ddx(xtmp1,xtmp2,0,0)
            call transpose_x_to_y(xtmp2,flux,this%decomp)
            rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
            flux = this%Wcnsrv(:,:,:,mom_index+2)*this%u
            call transpose_y_to_x(flux,xtmp1,this%decomp)
            call this%der%ddx(xtmp1,xtmp2,0,0)
            call transpose_x_to_y(xtmp2,flux,this%decomp)
            rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
            flux = (this%Wcnsrv(:,:,:, TE_index) + this%p)*this%u
            call transpose_y_to_x(flux,xtmp1,this%decomp)
            call this%der%ddx(xtmp1,xtmp2,0,0)
            call transpose_x_to_y(xtmp2,flux,this%decomp)
            rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index) - flux

            ! y-derivatives, convective conservative form
           ! select case(this%mix%ns)
           ! case(1)
           !     flux = this%Wcnsrv(:,:,:,mom_index+1)   ! mass
           !     call this%der%ddy(flux,tmp,0,0)
           !     rhs(:,:,:,1) = rhs(:,:,:,1) - tmp
           ! case default
           !     do i = 1,this%mix%ns
           !         flux = this%Wcnsrv(:,:,:,i)*this%v   ! mass
           !         call this%der%ddy(flux,tmp,0,0)
           !         rhs(:,:,:,i) = rhs(:,:,:,i) - tmp
           !     end do
           ! end select
            flux = this%Wcnsrv(:,:,:,mom_index)*this%v     ! x-momentum
            call this%der%ddy(flux,tmp,0,0)
            rhs(:,:,:,mom_index) = rhs(:,:,:,mom_index) - tmp
            flux = this%Wcnsrv(:,:,:,mom_index+1)*this%v + this%p  ! y-momentum
            call this%der%ddy(flux,tmp,0,0)
            rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - tmp
            flux = this%Wcnsrv(:,:,:,mom_index+2)*this%v    ! z-momentum
            call this%der%ddy(flux,tmp,0,0)
            rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - tmp
            flux = (this%Wcnsrv(:,:,:, TE_index) + this%p)*this%v  ! TotalEnergy
            call this%der%ddy(flux,tmp,0,0)
            rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) - tmp

            ! z-derivatives, convective conservative form
            !select case(this%mix%ns)
            !case(1)
            !    flux = this%Wcnsrv(:,:,:,mom_index+2)   ! mass
            !    call transpose_y_to_z(flux,ztmp1,this%decomp)
            !    call this%der%ddz(ztmp1,ztmp2,0,0)
            !    call transpose_z_to_y(ztmp2,flux,this%decomp)
            !    rhs(:,:,:,1) = rhs(:,:,:,1) - flux
            !case default
            !    do i = 1,this%mix%ns
            !        flux = this%Wcnsrv(:,:,:,i)*this%w   ! mass
            !        call transpose_y_to_z(flux,ztmp1,this%decomp)
            !        call this%der%ddz(ztmp1,ztmp2,0,0)
            !        call transpose_z_to_y(ztmp2,flux,this%decomp)
           !         rhs(:,:,:,i) = rhs(:,:,:,i) - flux
           !     end do
           ! end select
            flux = this%Wcnsrv(:,:,:,mom_index)*this%w     ! x-momentum
            call transpose_y_to_z(flux,ztmp1,this%decomp)
            call this%der%ddz(ztmp1,ztmp2,0,0)
            call transpose_z_to_y(ztmp2,flux,this%decomp)
            rhs(:,:,:,mom_index) = rhs(:,:,:,mom_index) - flux
            flux = this%Wcnsrv(:,:,:,mom_index+1)*this%w   ! y-momentum
            call transpose_y_to_z(flux,ztmp1,this%decomp)
            call this%der%ddz(ztmp1,ztmp2,0,0)
            call transpose_z_to_y(ztmp2,flux,this%decomp)
            rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
            flux = this%Wcnsrv(:,:,:,mom_index+2)*this%w + this%p  ! z-momentum
            call transpose_y_to_z(flux,ztmp1,this%decomp)
            call this%der%ddz(ztmp1,ztmp2,0,0)
            call transpose_z_to_y(ztmp2,flux,this%decomp)
            rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
            flux = (this%Wcnsrv(:,:,:, TE_index) + this%p)*this%w   ! Total
            call transpose_y_to_z(flux,ztmp1,this%decomp)
            call this%der%ddz(ztmp1,ztmp2,0,0)
            call transpose_z_to_y(ztmp2,flux,this%decomp)
            rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) - flux
       

        ! heat conduction (kap) terms, for energy equation
        dTdx => this%ybuf(:,:,:,1)
        dTdy => this%ybuf(:,:,:,2)
        dTdz => this%ybuf(:,:,:,3)
        d2Tdx2 => this%ybuf(:,:,:,4)
        d2Tdy2 => this%ybuf(:,:,:,5)
        d2Tdz2 => this%ybuf(:,:,:,6)
        dkapdx => NCbuff(:,:,:,1)
        dkapdy => NCbuff(:,:,:,2)
        dkapdz => NCbuff(:,:,:,3)
        call this%gradient(this%T,dTdx,dTdy,dTdz, this%x_bc, this%y_bc,this%z_bc)
        call this%secondder(this%der,this%T,d2Tdx2,d2Tdy2,d2Tdz2,[0,0],[0,0],[0,0])
        call this%gradient(this%kap,dkapdx,dkapdy,dkapdz, this%x_bc, this%y_bc,this%z_bc)
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + dkapdx*dTdx+dkapdy*dTdy+dkapdz*dTdz+this%kap*(d2Tdx2+d2Tdy2+d2Tdz2)

        ! stress (mu/bulk) terms
        bambda => NCbuff(:,:,:,8)             ! used for tauxx,tauyy,tauzz
        lambda => NCbuff(:,:,:,9)
        bambda = (four/three)*this%mu + this%bulk
        lambda = this%bulk - (two/three)*this%mu
        dmudx => this%ybuf(:,:,:,1)
        dmudy => this%ybuf(:,:,:,2)
        dmudz => this%ybuf(:,:,:,3)
        dbulkdx => this%ybuf(:,:,:,4)
        dbulkdy => this%ybuf(:,:,:,5)
        dbulkdz => this%ybuf(:,:,:,6)
        call this%gradient(this%mu,dmudx,dmudy,dmudz, this%x_bc, this%y_bc,this%z_bc)
        call this%gradient(this%bulk,dbulkdx,dbulkdy,dbulkdz, this%x_bc,this%y_bc, this%z_bc)
        !tau_xx,tau_xy,tau_xz setup
        d2udx2 => NCbuff(:,:,:,1)
        d2udy2 => NCbuff(:,:,:,2)
        d2udz2 => NCbuff(:,:,:,3)
        call this%secondder(this%der,this%u,d2udx2,d2udy2,d2udz2,[0,0],[0,0],[0,0])
        !tau_xx in x-momentum and energy eq
        flux =(four/three*dmudx+dbulkdx)*dudx+(dbulkdx-two/three*dmudx)*(dvdy+dwdz)+bambda*d2udx2
        tmp = dvdy+dwdz
        call transpose_y_to_x(tmp,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,0,0)
        call transpose_x_to_y(xtmp2,tmp,this%decomp)                   !tmp =
        flux = flux+lambda*tmp
        rhs(:,:,:, mom_index) = rhs(:,:,:, mom_index) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%u*flux
        !tau_xy in x-momentum and energy eq
        call this%der%ddy(dvdx,tmp,0,0)        !tmp = ddy(dvdx)
        flux = dmudy*(dudy+dvdx)+this%mu*(d2udy2+tmp)
        rhs(:,:,:, mom_index) = rhs(:,:,:, mom_index) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%u*flux
        !tau_xz in x-momentum and energy eq
        call transpose_y_to_z(dwdx,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,0,0)
        call transpose_z_to_y(ztmp2,tmp,this%decomp)                   !tmp =
        flux = dmudz*(dudz+dwdx)+this%mu*(d2udz2+tmp)
        rhs(:,:,:, mom_index) = rhs(:,:,:, mom_index) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%u*flux
        !tau_yy,tau_yx,tau_yz setup
        d2vdx2 => NCbuff(:,:,:,1)
        d2vdy2 => NCbuff(:,:,:,2)
        d2vdz2 => NCbuff(:,:,:,3)
        call this%secondder(this%der,this%v,d2vdx2,d2vdy2,d2vdz2,[0,0],[0,0],[0,0])
        !tau_yy in y-momentum and energy eq
        flux =(four/three*dmudy+dbulkdy)*dvdy+(dbulkdy-two/three*dmudy)*(dudx+dwdz)+bambda*d2vdy2
        call this%der%ddy(dudx+dwdz,tmp,0,0)     !tmp = ddy(dudx+dwdz)
        flux = flux+lambda*tmp
        rhs(:,:,:, mom_index+1) = rhs(:,:,:, mom_index+1) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%v*flux
        !tau_yx in y-momentum and energy eq
        call transpose_y_to_x(dudy,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,0,0)
        call transpose_x_to_y(xtmp2,tmp,this%decomp)                   !tmp =
        flux = dmudx*(dvdx+dudy)+this%mu*(d2vdx2+tmp)
        rhs(:,:,:, mom_index+1) = rhs(:,:,:, mom_index+1) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%v*flux
        !tau_yz in y-momentum and energy eq
        call transpose_y_to_z(dwdy,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,0,0)
        call transpose_z_to_y(ztmp2,tmp,this%decomp)                   !tmp =
        flux = dmudz*(dvdz+dwdy)+this%mu*(d2vdz2+tmp)
        rhs(:,:,:, mom_index+1) = rhs(:,:,:, mom_index+1) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%v*flux
        !tau_zz,tau_zx,tau_zy setup
        d2wdx2 => NCbuff(:,:,:,1)
        d2wdy2 => NCbuff(:,:,:,2)
        d2wdz2 => NCbuff(:,:,:,3)
        call this%secondder(this%der,this%w,d2wdx2,d2wdy2,d2wdz2,[0,0],[0,0],[0,0])
        !tau_zz in z-momentum and energy eq
        flux = (four/three*dmudz+dbulkdz)*dwdz+(dbulkdz-two/three*dmudz)*(dudx+dvdy)+bambda*d2wdz2
        tmp = dudx+dvdy
        call transpose_y_to_z(tmp,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,0,0)
        call transpose_z_to_y(ztmp2,tmp,this%decomp)                   !tmp =
        flux = flux+lambda*tmp
        rhs(:,:,:, mom_index+2) = rhs(:,:,:, mom_index+2) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%w*flux
        !tau_zx in z-momentum and energy eq
        call transpose_y_to_x(dudz,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,0,0)
        call transpose_x_to_y(xtmp2,tmp,this%decomp)                   !tmp =
        flux = dmudx*(dwdx+dudz)+this%mu*(d2wdx2+tmp)
        rhs(:,:,:, mom_index+2) = rhs(:,:,:, mom_index+2) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%w*flux
        !tau_zy in z-momentum and energy eq
        call this%der%ddy(dvdz,tmp,0,0)        !tmp = ddy(dvdz)
        flux = dmudy*(dwdy+dvdz)+this%mu*(d2wdy2+tmp)
        rhs(:,:,:, mom_index+2) = rhs(:,:,:, mom_index+2) + flux
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + this%w*flux
        !finish terms in energy eq
        NCbuff = duidxj
        dudx => NCbuff(:,:,:,1); dudy => NCbuff(:,:,:,2); dudz => NCbuff(:,:,:,3);
        dvdx => NCbuff(:,:,:,4); dvdy => NCbuff(:,:,:,5); dvdz => NCbuff(:,:,:,6);
        dwdx => NCbuff(:,:,:,7); dwdy => NCbuff(:,:,:,8); dwdz => NCbuff(:,:,:,9);
        call this%get_tau( duidxj, viscwork )
        tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx)
        tauxz => duidxj(:,:,:,tauxzidx); tauyy => duidxj(:,:,:,tauyyidx) 
        tauyz => duidxj(:,:,:,tauyzidx); tauzz => duidxj(:,:,:,tauzzidx)
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + dudx*tauxx+dudy*tauxy+dudz*tauxz
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + dvdx*tauxy+dvdy*tauyy+dvdz*tauyz
        rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + dwdx*tauxz+dwdy*tauyz+dwdz*tauzz
     

         call this%mix%get_J(this%rho)

        if (this%mix%ns .GT. 3) then
            call GracefulExit("Only up to 3 species are currently supported",3214)
        end if
        ! calculating dJ/dx terms
        if (this%mix%ns .GT. 1) then
            d2Ysdx2 => NCbuff(:,:,:,1:this%mix%ns)
            d2Ysdy2 => NCbuff(:,:,:,this%mix%ns+1:2*this%mix%ns)
            d2Ysdz2 => NCbuff(:,:,:,2*this%mix%ns+1:3*this%mix%ns)
            dsumJxdx => this%ybuf(:,:,:,1); dsumJydy => this%ybuf(:,:,:,2);dsumJzdz => this%ybuf(:,:,:,3); 
            dJxdx => this%ybuf(:,:,:,4); dJydy => this%ybuf(:,:,:,5); dJzdz =>this%ybuf(:,:,:,6);
            do i = 1,this%mix%ns
                call this%secondder(this%der,this%mix%material(i)%Ys,d2Ysdx2(:,:,:,i),d2Ysdy2(:,:,:,i),d2Ysdz2(:,:,:,i),[0,0],[0,0],[0,0])
            end do
            !x-derivatives
            do i = 1,this%mix%ns
                dsumJxdx = zero
                if (this%mix%ns .GT. 2) then
                    do k = 1,this%mix%ns
                    ! this calculates the molecular diffusion flux correction
                        flux = this%Wcnsrv(:,:,:,i)*this%mix%material(k)%diff       
                        call transpose_y_to_x(flux,xtmp1,this%decomp)
                        call this%der%ddx(xtmp1,xtmp2,0,0)
                        call transpose_x_to_y(xtmp2,tmp,this%decomp)
                        dsumJxdx =dsumJxdx+tmp*dYsdx(:,:,:,k)+flux*d2Ysdx2(:,:,:,k)
                    end do
                end if
                flux = this%rho*this%mix%material(i)%diff
                call transpose_y_to_x(flux,xtmp1,this%decomp)
                call this%der%ddx(xtmp1,xtmp2,0,0)
                call transpose_x_to_y(xtmp2,tmp,this%decomp)
                dJxdx = tmp*dYsdx(:,:,:,i) + flux*d2Ysdx2(:,:,:,i) - dsumJxdx
                !species
                rhs(:,:,:,i) = rhs(:,:,:,i) + dJxdx
                !energy
                call this%mix%material(i)%get_enthalpy(tmp)    ! tmp= h_i
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + tmp*dJxdx
            end do
            !y-derivatives
            do i = 1,this%mix%ns
                dsumJydy = zero
                if (this%mix%ns .GT. 2) then
                    do k = 1,this%mix%ns
                    ! this calculates the molecular diffusion flux correction
                        flux = this%Wcnsrv(:,:,:,i)*this%mix%material(k)%diff
                        call this%der%ddy(flux,tmp,0,0)
                        dsumJydy = dsumJydy+tmp*dYsdy(:,:,:,k)+flux*d2Ysdy2(:,:,:,k)
                    end do
                end if
                flux = this%rho*this%mix%material(i)%diff
                call this%der%ddy(flux,tmp,0,0)
                dJydy = tmp*dYsdy(:,:,:,i) + flux*d2Ysdy2(:,:,:,i) - dsumJydy
                !species
                rhs(:,:,:,i) = rhs(:,:,:,i) + dJydy
                !energy
                call this%mix%material(i)%get_enthalpy(tmp)    ! tmp= h_i
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + tmp*dJydy
            end do
            !z-derivatives
            do i = 1,this%mix%ns
                dsumJzdz = zero
                if (this%mix%ns .GT. 2) then
                    do k = 1,this%mix%ns
                    ! this calculates the molecular diffusion flux correction
                        flux = this%Wcnsrv(:,:,:,i)*this%mix%material(k)%diff
                        call transpose_y_to_z(flux,ztmp1,this%decomp)
                        call this%der%ddz(ztmp1,ztmp2,0,0)
                        call transpose_z_to_y(ztmp2,tmp,this%decomp)
                        dsumJzdz = dsumJzdz+tmp*dYsdz(:,:,:,k)+flux*d2Ysdz2(:,:,:,k)
                    end do
                end if
                flux = this%rho*this%mix%material(i)%diff
                call transpose_y_to_z(flux,ztmp1,this%decomp)
                call this%der%ddz(ztmp1,ztmp2,0,0)
                call transpose_z_to_y(ztmp2,tmp,this%decomp)
                dJzdz = tmp*dYsdz(:,:,:,i) + flux*d2Ysdz2(:,:,:,i) - dsumJzdz
                !species
                rhs(:,:,:,i) = rhs(:,:,:,i) + dJzdz
                !energy
                call this%mix%material(i)%get_enthalpy(tmp)    ! tmp= h_i
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + tmp*dJzdz
            end do

            !final J terms in energy equation:
          ! call this%get_J(gradYs)
          ! Jx => gradYs(:,:,:,1:this%mix%ns)
          ! Jy => gradYs(:,:,:,this%mix%ns+1:2*this%mix%ns)
          ! Jz => gradYs(:,:,:,2*this%mix%ns+1:3*this%mix%ns)
            do i = 1,this%mix%ns
                call this%mix%material(i)%get_enthalpy(flux)
                call transpose_y_to_x(flux,xtmp1,this%decomp)
                call this%der%ddx(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
                call transpose_x_to_y(xtmp2,tmp,this%decomp)
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) +tmp*this%mix%material(i)%Ji(:,:,:,1)
                call this%der%ddy(flux,tmp,this%y_bc(1),this%y_bc(2))
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + tmp*this%mix%material(i)%Ji(:,:,:,2)
                call transpose_y_to_z(flux,ztmp1,this%decomp)
                call this%der%ddz(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
                call transpose_z_to_y(ztmp2,tmp,this%decomp)
                rhs(:,:,:, TE_index) = rhs(:,:,:, TE_index) + tmp*this%mix%material(i)%Ji(:,:,:,3)
            end do
        end if 


       if(this%intSharp.AND.this%intSharp_cpl) then
           !calculate kinetic energy for intSharp terms
           ke = half*( this%u**2 + this%v**2 + this%w**2 ) !is this accesiblewithout recreating?

           if(this%intSharp_spf) then
              ! if(useNewSPF) then !this is unstable
              !    !new -- for useNewSPF = .TRUE. in solidmix
                 rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + this%mix%intSharp_f(:,:,:,1)
                 rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + this%mix%intSharp_f(:,:,:,2)
                 rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + this%mix%intSharp_f(:,:,:,3)
                 rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + this%mix%intSharp_h(:,:,:,1)
              ! else
              !    !original
              !    rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) +
              !    this%mix%intSharp_f(:,:,:,1)*this%u
              !    rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) +
              !    this%mix%intSharp_f(:,:,:,1)*this%v
              !    rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) +
              !    this%mix%intSharp_f(:,:,:,1)*this%w
              !    rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) +
              !    this%mix%intSharp_f(:,:,:,1)*ke +
              !    this%mix%intSharp_h(:,:,:,1)
              ! endif

              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%u,this%mix%intSharp_fDiff(:,:,:,2)*this%u,this%mix%intSharp_fDiff(:,:,:,3)*this%u,tmp,this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%v,this%mix%intSharp_fDiff(:,:,:,2)*this%v,this%mix%intSharp_fDiff(:,:,:,3)*this%v,tmp,-this%x_bc,this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%w,this%mix%intSharp_fDiff(:,:,:,2)*this%w,this%mix%intSharp_fDiff(:,:,:,3)*this%w,tmp,-this%x_bc,-this%y_bc,this%z_bc)
              rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*ke + this%mix%intSharp_hDiff(:,:,:,1),this%mix%intSharp_fDiff(:,:,:,2)*ke + this%mix%intSharp_hDiff(:,:,:,2),this%mix%intSharp_fDiff(:,:,:,3)*ke + this%mix%intSharp_hDiff(:,:,:,3),tmp,-this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + tmp

           else

              !low order terms
              call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*this%u,this%mix%intSharp_f(:,:,:,2)*this%u,this%mix%intSharp_f(:,:,:,3)*this%u,tmp,this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + tmp

              call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*this%v,this%mix%intSharp_f(:,:,:,2)*this%v,this%mix%intSharp_f(:,:,:,3)*this%v,tmp,-this%x_bc,this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + tmp

              call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*this%w,this%mix%intSharp_f(:,:,:,2)*this%w,this%mix%intSharp_f(:,:,:,3)*this%w,tmp,-this%x_bc,-this%y_bc,this%z_bc)
              rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + tmp

              call divergence(this%decomp,this%derD02,this%mix%intSharp_f(:,:,:,1)*ke + this%mix%intSharp_h(:,:,:,1),this%mix%intSharp_f(:,:,:,2)*ke + this%mix%intSharp_h(:,:,:,2),this%mix%intSharp_f(:,:,:,3)*ke + this%mix%intSharp_h(:,:,:,3),tmp,-this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + tmp


              !high order terms
              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%u,this%mix%intSharp_fDiff(:,:,:,2)*this%u,this%mix%intSharp_fDiff(:,:,:,3)*this%u,tmp,this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%v,this%mix%intSharp_fDiff(:,:,:,2)*this%v,this%mix%intSharp_fDiff(:,:,:,3)*this%v,tmp,-this%x_bc,this%y_bc,-this%z_bc)
              rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*this%w,this%mix%intSharp_fDiff(:,:,:,2)*this%w,this%mix%intSharp_fDiff(:,:,:,3)*this%w,tmp,-this%x_bc,-this%y_bc,this%z_bc)
              rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + tmp

              call divergence(this%decomp,this%der,this%mix%intSharp_fDiff(:,:,:,1)*ke + this%mix%intSharp_hDiff(:,:,:,1),this%mix%intSharp_fDiff(:,:,:,2)*ke + this%mix%intSharp_hDiff(:,:,:,2),this%mix%intSharp_fDiff(:,:,:,3)*ke + this%mix%intSharp_hDiff(:,:,:,3),tmp,-this%x_bc,-this%y_bc,-this%z_bc)
              rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + tmp

              !FV sharpening
              rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + this%mix%intSharp_fFV(:,:,:,1)
              rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + this%mix%intSharp_fFV(:,:,:,2)
              rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + this%mix%intSharp_fFV(:,:,:,3)
              rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + this%mix%intSharp_hFV

           endif
    endif
        


        !!Surface Tension
        if (this%use_surfaceTension) then
            rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) + this%mix%surfaceTension_f(:,:,:,1)
            rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) + this%mix%surfaceTension_f(:,:,:,2)
            rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) + this%mix%surfaceTension_f(:,:,:,3)
            rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) + this%mix%surfaceTension_e
        endif

        ! Call problem source hook
        
       ! Call problem source hook
       call hook_mixture_source(this%decomp, this%mesh, this%fields, this%mix, this%tsim, rhs)




        call this%get_tau( duidxj, viscwork )
        ! Now, associate the pointers to understand what's going on better
        tauxx => duidxj(:,:,:,tauxxidx); tauxy => duidxj(:,:,:,tauxyidx); tauxz => duidxj(:,:,:,tauxzidx);
        tauyy => duidxj(:,:,:,tauyyidx); tauyz => duidxj(:,:,:,tauyzidx); tauzz => duidxj(:,:,:,tauzzidx);
!print '(a,9(e21.14,1x))', 'dudx ', duidxj(89,1,1,1:9)
!print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)
        ! Add the deviatoric stress to the tau for use in fluxes
        tauxx = tauxx + this%sxx; tauxy = tauxy + this%sxy; tauxz = tauxz + this%sxz
                                  tauyy = tauyy + this%syy; tauyz = tauyz + this%syz
                                                            tauzz = tauzz + this%szz
!print '(a,9(e21.14,1x))', 'tauxx', tauxx(89,1,1)

        ! store artificial stress tensor in devstress. this should not break
        ! anything since devstress will be
        ! overwritten in get_primitive and post_bc. used in update_eh -- NSG
        this%sxx = tauxx - this%sxx; this%sxy = tauxy - this%sxy; this%sxz =tauxz - this%sxz
                                     this%syy = tauyy - this%syy; this%syz =tauyz - this%syz
                                                                  this%szz =tauzz - this%szz

    end subroutine

    subroutine getRHS_xStagg( this,  rhs, tauxx,tauxy,tauxz,Frho,Fenergy, qx)
        use operators, only: gradFV_x, interpolateFV_x,gradFV_N2Fx,filter3D,gradFV_N2Fy,gradient
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv),intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxx,tauxy,tauxz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qx
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: Frho,Fenergy
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: buff, flux,TE, u_int,v_int, w_int, p_int, tauxx_int, tauxy_int,tauxz_int, qx_int, e_int, rho_int, rhodiff_int, rhoe_prim, gam, num, t_int, KE,e_prim,p4, Eint,gradu,gradp,gradup, UU, kef,gradm1,gradm2,gradrhou,clocal, delp,gradVFx,gradVFy,UV,umag,soslocal,sos1,sos2
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhou_int,rhov_int, rhow_int, rhoe_int,spe_int, rhoYs_int, den,gradRYs, tauRho_mid, rhom, rhom_int,ke_int,sos_int,Mu_int,Mv_int,Mw_int,H_int,delrhou,delrhov,pbar,GVFmag_x
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: u_int6, rho_int6, t_int6, u_int8, t_int8, rho_int8,spec_int,tmp1,tmp2,tmp3
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: EfluxI, rhouI,rhovI,rhouvI,rhouuI,GVFmag
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,2) :: VF_int, M_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradVF
        real(rkind),dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xtmp1,xtmp2,xtmp3,delptmp,clocaltmp,rhovtmp,rhoutmp,xtmp4
        integer :: i

        p_int = this%p_mid(:,:,:,1)
        u_int = this%u_mid(:,:,:,1) !this%u_mid(:,:,:,1) !spec_int*rhou_int !this%u_mid(:,:,:,1)
        v_int = this%v_mid(:,:,:,1)
        w_int = this%w_mid(:,:,:,1)
 
        call interpolateFV_x(this%decomp,this%interpMid02,this%p, pbar,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rho_int =0.0d0
        num = 0d0
        gam = 0d0
        KE = 0d0
        den = 0.0d0

        rhou_int = 0.0d0
        rhov_int = 0.0d0
        rhow_int = 0.0d0
        rhoe_prim = 0.0d0

        Mu_int = 0d0
        Mv_int = 0d0
        Mw_int = 0.0d0

        do i = 1,2
            
           VF_int(:,:,:,i) = this%mix%material(i)%VF_mid(:,:,:,1)
           M_int(:,:,:,i)  = this%mix%material(i)%rhoYs_mid(:,:,:,1)
           rho_int  = rho_int + M_int(:,:,:,i)

!          call interpolateFV_x(this%decomp,this%interpMid,this%u*this%mix%material(i)%consrv(:,:,:,1),tmp1,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
!          call interpolateFV_x(this%decomp,this%interpMid,this%v*this%mix%material(i)%consrv(:,:,:,1),tmp2,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
!          call interpolateFV_x(this%decomp,this%interpMid,this%w*this%mix%material(i)%consrv(:,:,:,1),tmp3,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
           Mu_int = Mu_int + M_int(:,:,:,i)*this%u_int ! tmp1
           Mv_int = Mv_int + M_int(:,:,:,i)*v_int ! tmp2
           Mw_int = Mw_int + M_int(:,:,:,i)*w_int ! tmp3

        enddo



       rhou_int = Mu_int*u_int
       rhov_int = Mu_int*v_int
       rhow_int = Mu_int*w_int
       KE  = 0.5_rkind*rho_int*(u_int*u_int + v_int*v_int + w_int*w_int) ! ( Mu_int**2 + Mv_int**2 + Mw_int**2 ) / rho_int

      do i = 1,2

         rhoe_prim = rhoe_prim +  this%mix%material(i)%hydro%onebygam_m1*( p_int + this%mix%material(i)%hydro%gam*this%mix%material(i)%hydro%Pinf)/(M_int(:,:,:,i)/VF_int(:,:,:,i) ) *  M_int(:,:,:,i)/rho_int
      enddo

      H_int = KE*this%u_int +(p_int/rho_int + rhoe_prim ) *Mu_int + pbar*(this%u_int - u_int )

          

      GVFmag_x=0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(this%use_CnsrvSurfaceTension) then
            call gradient(this%decomp,this%derCD06,this%mix%material(1)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))
            GVFmag = this%surfaceTension_coeff*sqrt(gradVF(:,:,:,1)**2.0_rkind + gradVF(:,:,:,2)**2.0_rkind + gradVF(:,:,:,3)**2.0_rkind )
            call interpolateFV_x(this%decomp,this%interpMid,GVFmag,GVFmag_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        endif
        flux = 0.0
        buff = rhou_int + p_int - tauxx - Frho*u_int !- this%CP*gradVF*clocal*this%dx**2*delrhou

        call gradFV_x(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux  
        this%xflux_x = flux
        flux =0.0
        buff =rhov_int - tauxy - Frho*v_int !- this%CP*gradVF*clocal*this%dx**2*delrhov !y-momentum

        !endif

        call gradFV_x(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
        this%xflux_y = flux

        buff = rhow_int  - tauxz - Frho*w_int !z-momentum
        flux = 0.0

        !endif
        call gradFV_x(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
 
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
        this%xflux_z = flux
        flux = 0.0

        !!!!!!!!!!!!!! add back in KE         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        buff = H_int  -( tauxx)*u_int - (v_int*tauxy - w_int*tauxz) - 0.5_rkind*(u_int*u_int +v_int*v_int + w_int*w_int)*Frho - Fenergy        !  - this%CP*gradVF*clocal*this%dx**2*delp 
        call gradFV_x(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux  ! - gradup ! -0.5*( gradu*this%p + this%v*gradp)    
 
    end subroutine

    subroutine getRHS_yStagg( this,  rhs, tauxy,tauyy,tauyz,Frho,Fenergy,qy)
        use operators, only: gradFV_y,interpolateFV_y,gradFV_N2Fx,gradFV_N2Fy,gradient
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,ncnsrv),intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxy,tauyy,tauyz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qy
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: Frho,Fenergy
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,2) :: VF_int,M_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: buff, flux,TE,u_int,v_int, w_int, p_int, tauxy_int, tauyy_int,tauyz_int, qy_int, e_int,rho_int, gam, num, rhoe_prim, KE, e_prim,gradu,gradup, UU, kef,delp,cl,cr,clocal,gradVFx,gradVFy,UV,VV,tmp1,tmp2,tmp3,umag,soslocal,sos1,sos2
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhou_int,rhov_int,rhow_int, rhoe_int, spe_int, rhoYs_int, den,ke_int,p4,Eint, up_int, cpressure, gradcpressure, gradp,sos_int, H_int, Mu_int,delrhou,delrhov,Mv_int,Mw_int,pbar,GVFmag,GVFmag_y
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: EfluxI,rhouI,rhovI,rhovvI,rhouvI
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradVF
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        real(rkind) :: g = 0.1, c = 1d4
        integer :: i,j,k

       u_int = this%u_mid(:,:,:,2)
       v_int = this%v_mid(:,:,:,2)
       w_int = this%w_mid(:,:,:,2)
       sos_int = (this%sos + sqrt( this%u*this%u + this%v*this%v))*this%rho
       p_int = this%p_mid(:,:,:,2)
       call interpolateFV_y(this%decomp,this%interpMid02,this%p,pbar,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       rho_int = 0.0d0
       rhoe_prim = 0.0d0
       num = 0.0d0
       gam = 0.0d0
       den = 0.0d0

       Mu_int = 0.0d0
       Mv_int = 0.0d0
       Mw_int = 0.0d0

        do i = 1,2

           VF_int(:,:,:,i) = this%mix%material(i)%VF_mid(:,:,:,2)
           M_int(:,:,:,i)  = this%mix%material(i)%rhoYs_mid(:,:,:,2)
           rho_int  = rho_int + M_int(:,:,:,i)            

           Mu_int = Mu_int + M_int(:,:,:,i)*u_int !tmp1
           Mv_int = Mv_int + M_int(:,:,:,i)*this%v_int !tmp2
           Mw_int = Mw_int + M_int(:,:,:,i)*w_int !tmp3

        enddo

        

       rhou_int = Mv_int*u_int
       rhov_int = Mv_int*v_int
       rhow_int = Mv_int*w_int
       KE  = 0.5_rkind*rho_int*(u_int*u_int + v_int*v_int + w_int*w_int ) !  0.5_rkind* ( Mu_int**2 + Mv_int**2 + Mw_int**2 ) / rho_int 

    
        do i = 1,2

         rhoe_prim = rhoe_prim + this%mix%material(i)%hydro%onebygam_m1*(p_int +this%mix%material(i)%hydro%gam*this%mix%material(i)%hydro%Pinf)/(M_int(:,:,:,i)/VF_int(:,:,:,i)) *  M_int(:,:,:,i)/rho_int
       enddo


       H_int = KE*this%v_int +( p_int / rho_int + rhoe_prim)*Mv_int + pbar*(this%v_int - v_int)



        flux = 0.0
        buff = rhou_int  - tauxy - Frho*u_int  !x-momentum 
        
        call gradFV_y(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux 
        this%yflux_x = flux

        flux = 0.0
        buff = rhov_int - tauyy + p_int  - Frho*v_int
        call gradFV_y(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
!        call this%filter(gradp, this%fil,1,-this%x_bc,this%y_bc,this%z_bc)
        !endif
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1)  - flux  
        this%yflux_y = flux

        flux = 0.0
        buff = rhow_int  - tauyz - Frho*w_int !z-momentum

        !endif
        call gradFV_y(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
        this%yflux_z = flux

        flux = 0.0d0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Add back in KE         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        buff = H_int- (tauyy)*v_int - u_int*tauxy -w_int*tauyz - 0.5_rkind*(u_int*u_int+v_int*v_int + w_int*w_int)*Frho - Fenergy 

        !endif


        call gradFV_y(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux 
        this%yflux_e = flux


    end subroutine

   
    subroutine getRHS_zStagg( this,  rhs, tauxz,tauyz,tauzz,Frho,Fenergy,qz)
        use operators, only: gradFV_z, interpolateFV_z,gradient
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp,this%nzp,ncnsrv),intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) ::tauxz,tauyz,tauzz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: Frho,Fenergy
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: buff,flux,TE,u_int,v_int, w_int, p_int, tauxz_int, tauzz_int, tauyz_int, qz_int,e_int,rho_int,KE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhou_int,rhov_int, rhow_int, rhoe_int, spe_int, rhoYs_int, den,  num,gam,rhoe_prim,GVFmag,GVFmag_z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: Mu_int,Mv_int,Mw_int,H_int,pbar
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,2) :: VF_int,M_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradVF
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        integer :: i

       u_int = this%u_mid(:,:,:,3)
       v_int = this%v_mid(:,:,:,3)
       w_int = this%w_mid(:,:,:,3)
       p_int = this%p_mid(:,:,:,3)

       call interpolateFV_z(this%decomp,this%interpMid02,this%p, pbar,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       rho_int = 0.0d0
       rhoe_prim = 0.0d0
       num = 0.0d0
       gam = 0.0d0
       den = 0.0d0

       Mu_int = 0.0d0
       Mv_int = 0.0d0
       Mw_int = 0.0d0
       rhou_int = 0.0d0
       rhov_int = 0.0d0
       rhow_int = 0.0d0
       KE = 0.0d0
        do i = 1,2

           VF_int(:,:,:,i) = this%mix%material(i)%VF_mid(:,:,:,3)
           M_int(:,:,:,i)  = this%mix%material(i)%rhoYs_mid(:,:,:,3)
           rho_int  = rho_int + M_int(:,:,:,i)

           Mu_int = Mu_int + M_int(:,:,:,i)*u_int !tmp1
           Mv_int = Mv_int + M_int(:,:,:,i)*v_int !tmp2
           Mw_int = Mw_int + M_int(:,:,:,i)*this%w_int

           KE = KE+ 0.5_rkind*(u_int*u_int + v_int*v_int + w_int*w_int )*this%w_int*this%mix%material(i)%rhoYs_mid(:,:,:,3)
           rhou_int = rhou_int + this%mix%material(i)%rhoYs_mid(:,:,:,3)*u_int*w_int
           rhov_int = rhov_int + this%mix%material(i)%rhoYs_mid(:,:,:,3)*v_int*w_int
           rhow_int = rhow_int + this%mix%material(i)%rhoYs_mid(:,:,:,3)*w_int*this%w_int !this%w_int
           rhoe_prim = rhoe_prim +this%mix%material(i)%VF_mid(:,:,:,3)*this%mix%material(i)%hydro%gam*this%mix%material(i)%hydro%onebygam_m1*(p_int +this%mix%material(i)%hydro%Pinf)*this%w_int
        enddo


!      rhou_int = Mw_int*u_int
!      rhov_int = Mw_int*v_int
!      rhow_int = Mw_int*w_int
!       KE  = 0.5*(u_int**2 + v_int**2 + w_int**2 ) !  0.5_rkind* ( Mu_int**2 + Mv_int**2 + Mw_int**2 ) / rho_int

!      do i = 1,2
!
!         rhoe_prim = rhoe_prim + this%mix%material(i)%hydro%onebygam_m1*(p_int +this%mix%material(i)%hydro%gam*this%mix%material(i)%hydro%Pinf)/(M_int(:,:,:,i)/VF_int(:,:,:,i)) *  M_int(:,:,:,i)/rho_int
!       enddo


        H_int = KE + rhoe_prim + pbar*(this%w_int - w_int)
        flux = 0.0
        buff = rhou_int  - tauxz - Frho*u_int
        call gradFV_z(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux
        this%zflux_x = flux 

        flux = 0.0
        buff = rhov_int   - tauyz - Frho*v_int

        call gradFV_z(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
        this%zflux_y = flux

        flux = 0.0
        buff = rhow_int + p_int   - tauzz - Frho*w_int !z-momentum

        call gradFV_z(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
        this%zflux_z = flux
        flux = 0.0


        buff = KE + rhoe_prim  - tauzz*w_int - u_int*tauxz-v_int*tauyz  - 0.5_rkind*(u_int*u_int + v_int*v_int + w_int*w_int)*Frho - Fenergy
        !buff = ( TE + p_int )*w_int

        call gradFV_z(this%decomp,this%derStagg,buff,flux,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux
        this%zflux_e = flux


    end subroutine
    
    subroutine getRHS_x( this,  rhs, tauxx,tauxy,tauxz, qx)
        use operators, only: divergence,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_y, gradFV_z, gradFV_x
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxx,tauxy,tauxz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qx
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: flux, ddx_p,p_int, u_int, ddx_up, ke
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        integer :: i

        xtmp1 => this%xbuf(:,:,:,1); xtmp2 => this%xbuf(:,:,:,2)

        !call interpolateFV_x(this%decomp,this%interpMid,this%p,p_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !call gradFV_x(this%decomp,this%derStagg,p_int,ddx_p,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !flux = this%Wcnsrv(:,:,:,mom_index  )*this%u - tauxx
      
        flux = this%Wcnsrv(:,:,:,mom_index  )*this%u + this%p - tauxx  ! x-momentum
!print *, 'flux 1', flux(89,1,1), this%u(89,1,1), this%p(89,1,1), tauxx(89,1,1)
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxx

        endif
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2, this%x_bc(1), this%x_bc(2)) ! Symmetric for x-momentum
!do i = 1, size(flux,1)
!  write(*,'(4(e21.14,1x))') xtmp1(i,1,1), xtmp2(i,1,1)
!enddo
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        this%xflux_x = flux

        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux 
        flux = this%Wcnsrv(:,:,:,mom_index  )*this%v          - tauxy   ! y-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxy

        endif

        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
        this%xflux_y = flux 

        flux = this%Wcnsrv(:,:,:,mom_index  )*this%w          - tauxz  ! z-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxz

        endif

        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
        this%xflux_z = flux

        !call interpolateFV_x(this%decomp,this%interpMid,this%u,u_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !call gradFV_x(this%decomp,this%derStagg,p_int*u_int,ddx_up,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !flux = (this%Wcnsrv(:,:,:, TE_index  ) - tauxx)*this%u -this%v*tauxy - this%w*tauxz + qx
        ke = 0.5*(this%u*this%u + this%v*this%v + this%w*this%w)
        flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauxx)*this%u - this%v*tauxy - this%w*tauxz + qx ! Total Energy
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxx*this%u -this%mix%surfaceTension_fxy*this%v - this%mix%surfaceTension_fxz*this%w

        endif
        call transpose_y_to_x(flux,xtmp1,this%decomp)
        call this%der%ddx(xtmp1,xtmp2,-this%x_bc(1),-this%x_bc(2)) ! Anti-symmetric for all but x-momentum
        call transpose_x_to_y(xtmp2,flux,this%decomp)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux 
        this%xflux_e = flux

    end subroutine

    subroutine getRHS_y( this,  rhs,&
                        tauxy,tauyy,tauyz,&
                            qy )
        use operators, only: divergence,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_y, gradFV_z, gradFV_x
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxy,tauyy,tauyz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qy

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: flux, p_int,ddy_p, v_int, ddy_vp
        real(rkind), dimension(:,:,:), pointer :: ytmp1

        ytmp1 => this%ybuf(:,:,:,6)

        flux = this%Wcnsrv(:,:,:,mom_index+1)*this%u          - tauxy ! x-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxy

        endif

        call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - ytmp1
        this%yflux_x = flux

       ! call interpolateFV_y(this%decomp,this%interpMid,this%p,p_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       ! call gradFV_y(this%decomp,this%derStagg,p_int,ddy_p,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        flux = this%Wcnsrv(:,:,:,mom_index+1)*this%v + this%p - tauyy ! y-momentum
        !flux = this%Wcnsrv(:,:,:,mom_index+1)*this%v  - tauyy
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fyy

        endif
        call this%der%ddy(flux,ytmp1, this%y_bc(1), this%y_bc(2)) ! Symmetric for y-momentum
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - ytmp1 
        this%yflux_y = flux

        flux = this%Wcnsrv(:,:,:,mom_index+1)*this%w          - tauyz !z-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fyz

        endif
        call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - ytmp1 
        this%yflux_z = flux

        !call interpolateFV_y(this%decomp,this%interpMid,this%v,v_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !call gradFV_y(this%decomp,this%derStagg,p_int*v_int,ddy_vp,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        !flux = (this%Wcnsrv(:,:,:, TE_index  ) - tauyy)*this%v - this%u*tauxy -this%w*tauyz + qy 
        flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauyy)*this%v - this%u*tauxy - this%w*tauyz + qy ! Total Energy
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxy*this%u - this%mix%surfaceTension_fyy*this%v - this%mix%surfaceTension_fyz*this%w  

        endif
        
        call this%der%ddy(flux,ytmp1,-this%y_bc(1),-this%y_bc(2)) ! Anti-symmetric for all but y-momentum
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - ytmp1 
        this%yflux_e = flux

    end subroutine

    subroutine getRHS_z(       this,  rhs,&
                        tauxz,tauyz,tauzz,&
                            qz )
        use operators, only: divergence,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_y, gradFV_z, gradFV_x
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, ncnsrv), intent(inout) :: rhs
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: tauxz,tauyz,tauzz
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: qz

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: flux, ddz_p,p_int, w_int, ddz_wp
        real(rkind), dimension(:,:,:), pointer :: ztmp1,ztmp2

        ztmp1 => this%zbuf(:,:,:,1); ztmp2 => this%zbuf(:,:,:,2)

        flux = this%Wcnsrv(:,:,:,mom_index+2)*this%u          - tauxz ! x-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxz

        endif
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - flux
        this%zflux_x = flux

        flux = this%Wcnsrv(:,:,:,mom_index+2)*this%v          - tauyz ! y-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fyz

        endif
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - flux
        this%zflux_y = flux

       ! call interpolateFV_z(this%decomp,this%interpMid,this%p,p_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       ! call gradFV_z(this%decomp,this%derStagg,p_int,ddz_p,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       ! flux = this%Wcnsrv(:,:,:,mom_index+2)*this%w  - tauzz 

        flux = this%Wcnsrv(:,:,:,mom_index+2)*this%w + this%p - tauzz ! z-momentum
        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fzz

        endif
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2, this%z_bc(1), this%z_bc(2)) ! Symmetric for z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - flux
        this%zflux_z = flux

        !flux = (this%Wcnsrv(:,:,:, TE_index  )  - tauzz)*this%w -this%u*tauxz - this%v*tauyz + qz
        flux = (this%Wcnsrv(:,:,:, TE_index  ) + this%p - tauzz)*this%w - this%u*tauxz - this%v*tauyz + qz ! Total Energy

        if(this%use_CnsrvSurfaceTension) then

            flux = flux - this%mix%surfaceTension_fxz*this%u - this%mix%surfaceTension_fyz*this%v - this%mix%surfaceTension_fzz*this%w

        endif

       ! call interpolateFV_z(this%decomp,this%interpMid,this%w,w_int,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
       ! call gradFV_z(this%decomp,this%derStagg,p_int*w_int,ddz_wp,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call transpose_y_to_z(flux,ztmp1,this%decomp)
        call this%der%ddz(ztmp1,ztmp2,-this%z_bc(1),-this%z_bc(2)) ! Anti-symmetric for all but z-momentum
        call transpose_z_to_y(ztmp2,flux,this%decomp)
        rhs(:,:,:, TE_index  ) = rhs(:,:,:, TE_index  ) - flux 
        this%zflux_e = flux

    end subroutine

    subroutine filter(this,arr,myfil,numtimes,x_bc_,y_bc_,z_bc_)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(inout) :: arr
        type(filters), target, optional, intent(in) :: myfil
        integer, optional, intent(in) :: numtimes
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc
        
        type(filters), pointer :: fil2use
        integer :: times2fil
        real(rkind), dimension(:,:,:), pointer :: tmp_in_y, tmp1_in_x, tmp1_in_z, tmp2_in_x, tmp2_in_z
        integer :: lastx, lasty, lastz, idx


        if (present(myfil)) then
            fil2use => myfil
        else
            fil2use => this%fil
        end if 

        if (present(numtimes)) then
            times2fil = numtimes
        else
            times2fil = 1
        end if

        ! Allocate pointers for the needed buffers 
        ! Atleast 2 buffers in x and z are assumed
        ! Last two buffers are occupied

        lastx = size(this%xbuf,4)
        lasty = size(this%ybuf,4)
        lastz = size(this%zbuf,4)

        tmp1_in_x => this%xbuf(:,:,:,lastx)
        tmp2_in_x => this%xbuf(:,:,:,lastx-1)
        tmp_in_y => this%ybuf(:,:,:,lasty)
        tmp1_in_z => this%zbuf(:,:,:,lastz)
        tmp2_in_z => this%zbuf(:,:,:,lastz-1)

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_
        
        ! First filter in y
        call fil2use%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
        ! Subsequent refilters 
        do idx = 1,times2fil-1
            arr = tmp_in_y
            call fil2use%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
        end do
        
        ! Then transpose to x
        call transpose_y_to_x(tmp_in_y,tmp1_in_x,this%decomp)

        ! First filter in x
        call fil2use%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_x = tmp2_in_x
            call fil2use%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
        end do 

        ! Now transpose back to y
        call transpose_x_to_y(tmp2_in_x,tmp_in_y,this%decomp)

        ! Now transpose to z
        call transpose_y_to_z(tmp_in_y,tmp1_in_z,this%decomp)

        !First filter in z
        call fil2use%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_z = tmp2_in_z
            call fil2use%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
        end do 

        ! Now transpose back to y
        call transpose_z_to_y(tmp2_in_z,arr,this%decomp)

        ! Finished
    end subroutine
  
     subroutine entropy_discreteKE( this)
        use operators, only: gradFV_x, interpolateFV_x,gradFV_N2Fx,filter3D,gradFV_y,gradFV_N2Fy,interpolateFV_y,interpolateFV,divergenceFV,gradFV_z,gradFV_N2Fz
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: buff,rhoLAD,adiff_fil,hi
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: rhom_int,gradRhom,gradRhomeh,gradRhomFace,rhoe_int,adiff_int,gradVF,VF_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: u_int6, rho_int6,t_int6, u_int8, t_int8, rho_int8,spec_int
        real(rkind), dimension(:,:,:), pointer :: xtmp1,xtmp2
        integer :: i

        this%entropy = 0

        do i = 1,2

         rhom_int = this%mix%material(i)%rhoYs_mid/(this%mix%material(i)%VF_mid + 1d-32)
         VF_int   = this%mix%material(i)%VF_mid 
         rhoe_int = (this%p_mid + this%mix%material(i)%hydro%gam *this%mix%material(i)%hydro%Pinf) / (this%mix%material(i)%hydro%gam - 1_rkind)

         call  this%mix%material(i)%get_enthalpy(hi)
  
         call gradFV_x(this%decomp,this%derStagg,VF_int(:,:,:,1),gradVF(:,:,:,1),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
         call gradFV_y(this%decomp,this%derStagg,VF_int(:,:,:,2),gradVF(:,:,:,2),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
         call gradFV_z(this%decomp,this%derStagg,VF_int(:,:,:,3),gradVF(:,:,:,3),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

         call gradFV_x(this%decomp,this%derStagg,rhoe_int(:,:,:,1),gradRhomeh(:,:,:,1),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
         call gradFV_y(this%decomp,this%derStagg,rhoe_int(:,:,:,2),gradRhomeh(:,:,:,2),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
         call gradFV_z(this%decomp,this%derStagg,rhoe_int(:,:,:,3),gradRhomeh(:,:,:,3),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)         

         call gradFV_N2Fx(this%decomp,this%derStagg,this%mix%material(i)%rhom,gradRhomFace(:,:,:,1),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
         call gradFV_N2Fy(this%decomp,this%derStagg,this%mix%material(i)%rhom,gradRhomFace(:,:,:,2),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
         call gradFV_N2Fz(this%decomp,this%derStagg,this%mix%material(i)%rhom,gradRhomFace(:,:,:,3),this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

         call divergenceFV(this%decomp,this%derStagg,this%mix%material(i)%adiff_stagg(:,:,:,1)*gradRhomFace(:,:,:,1),this%mix%material(i)%adiff_stagg(:,:,:,2)*gradRhomFace(:,:,:,2),this%mix%material(i)%adiff_stagg(:,:,:,3)*gradRhomFace(:,:,:,3),rhoLAD,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

         adiff_fil = this%mix%material(i)%adiff
         call this%filter(adiff_fil,this%fil,1,this%x_bc,this%y_bc,this%z_bc)
         this%entropy = this%entropy -this%mix%material(i)%VF*hi*rhoLAD + adiff_fil*(gradVF(:,:,:,1)*gradRhomeh(:,:,:,1) + gradVF(:,:,:,2)*gradRhomeh(:,:,:,2) &
                        + gradVF(:,:,:,3)*gradRhomeh(:,:,:,3) )
        enddo

 
    end subroutine
    subroutine getPhysicalProperties(this)
        use exits,      only: GracefulExit
        class(sgrid), intent(inout) :: this

        if (this%mix%ns > 2) then
            call GracefulExit("Number of species must be 1 or 2. for current &
                               implementation of getPhysicalProperties",928)
        endif

        ! If inviscid set everything to zero (otherwise use a model)
        this%mu   = this%phys_mu1   * this%mix%material(1)%VF
        this%bulk = this%phys_bulk1 * this%mix%material(1)%VF
        !this%mix%material(1)%physmu = this%phys_mu1
        !this%mix%material(2)%physmu = this%phys_mu2

        if (this%mix%ns .eq. 2) then
            this%mu   = this%mu   + this%phys_mu2   * this%mix%material(2)%VF
            this%bulk = this%bulk + this%phys_bulk2 * this%mix%material(2)%VF
        endif


        if (this%PTeqb) then
            this%kap  = this%phys_kap1  * this%mix%material(1)%VF
            if (this%mix%ns .eq. 2) then
                this%kap  = this%kap  + this%phys_kap2  * this%mix%material(2)%VF
            endif
        endif

    end subroutine  

    ! Get tau_ij for momentum equation and simulataneously calculate the viscous work
    ! for the material hydrodynamic equations
    subroutine get_tau(this,duidxj,viscwork)
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), target, intent(inout) :: duidxj
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),           intent(out)   :: viscwork
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: NCbuff, Mubuff   
        real(rkind), dimension(:,:,:), pointer :: d2udx2,d2udy2,d2udz2,d2vdx2,d2vdy2,d2vdz2,d2wdx2,d2wdy2,d2wdz2
        real(rkind), dimension(:,:,:), pointer :: dmudx,dmudy,dmudz
        real(rkind), dimension(:,:,:), pointer :: dbulkdx,dbulkdy,dbulkdz
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: lambda, bambda

        lambda => this%ybuf(:,:,:,1)
        bambda => this%ybuf(:,:,:,2)

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
       
        ! Compute the multiplying factors (thermo-shit)
        bambda = (four/three)*this%mu + this%bulk
        lambda = this%bulk - (two/three)*this%mu
    
    
        ! Step 1: Get tau_12  (dudy is destroyed)
        dudy =  dudy + dvdx
        viscwork  = (this%mu*dudy) * (dudy)  ! tau_12 * (S_12+S_21) (Since symmetric)
        dudy = this%mu*dudy
        !tauxyidz = 2
    
        ! Step 2: Get tau_13 (dudz is destroyed)
        dudz = dudz + dwdx
        viscwork  = viscwork + (this%mu*dudz) * (dudz)  ! tau_13 *(S_13 + S_31)
        dudz = this%mu*dudz
        !tauxzidx = 3

        ! Step 3: Get tau_23 (dvdz is destroyed)
        dvdz = dvdz + dwdy
        viscwork  = viscwork + (this%mu*dvdz) * (half*dvdz)  ! tau_23 * (S_23 + S_32)
        dvdz = this%mu*dvdz
        !tauyzidx = 6

        ! Step 4: Get tau_11 (dvdx is destroyed)
        dvdx = bambda*dudx + lambda*(dvdy + dwdz)
        viscwork  = viscwork + dvdx * dudx  ! tau_11 * S_11
        !tauxxidx = 4

        ! Step 5: Get tau_22 (dwdx is destroyed)
        dwdx = bambda*dvdy + lambda*(dudx + dwdz)
        viscwork  = viscwork + dwdx * dvdy  ! tau_22 * S_22
        !tauyyidx = 7

        ! Step 6: Get tau_33 (dwdy is destroyed)
        dwdy = bambda*dwdz + lambda*(dudx + dvdy)
        viscwork  = viscwork + dwdy * dwdz  ! tau_33 * S_33
        !tauzzidx = 8
    
        ! Done 
    end subroutine 


    subroutine get_tauStagg(this,duidxj,duidxj_int,duidxj_s)
        use operators, only : interpolateFV_x, interpolateFV_y, interpolateFV_z
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), target, intent(inout) :: duidxj, duidxj_s
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,12), target,intent(inout) :: duidxj_int
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: NCbuff, Mubuff
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: mu_x, mu_y, mu_z,bulk_x, bulk_y, bulk_z
        real(rkind), dimension(:,:,:), pointer :: dmudx,dmudy,dmudz
        real(rkind), dimension(:,:,:), pointer :: dbulkdx,dbulkdy,dbulkdz
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: lambda, bambda
        real(rkind), dimension(:,:,:), pointer :: dudx_s,dudy_s,dudz_s,dvdx_s,dvdy_s,dvdz_s,dwdx_s,dwdy_s,dwdz_s
        real(rkind), dimension(:,:,:), pointer :: dvdy_x,dwdz_x, dvdx_y, dwdx_z
        real(rkind), dimension(:,:,:), pointer :: dudx_y, dwdz_y, dudy_x, dwdy_z
        real(rkind), dimension(:,:,:), pointer :: dudx_z, dvdy_z, dudz_x, dvdz_y

        call interpolateFV_x(this%decomp,this%interpMid,this%mu,mu_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%mu,mu_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%mu,mu_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        call interpolateFV_x(this%decomp,this%interpMid,this%bulk,bulk_x,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,this%bulk,bulk_y,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,this%bulk,bulk_z,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

        lambda => this%ybuf(:,:,:,1)
        bambda => this%ybuf(:,:,:,2)

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

        dvdy_x => duidxj_int(:,:,:,1); dudx_y => duidxj_int(:,:,:,2); dudx_z => duidxj_int(:,:,:,3);
        dwdz_x => duidxj_int(:,:,:,4); dwdz_y => duidxj_int(:,:,:,5); dvdy_z => duidxj_int(:,:,:,6);
        dvdx_y => duidxj_int(:,:,:,7); dudy_x => duidxj_int(:,:,:,8); dudz_x => duidxj_int(:,:,:,9);
        dwdx_z => duidxj_int(:,:,:,10); dwdy_z => duidxj_int(:,:,:,11); dvdz_y => duidxj_int(:,:,:,12);

        dudx_s => duidxj_s(:,:,:,1); dudy_s => duidxj_s(:,:,:,2); dudz_s => duidxj_s(:,:,:,3);
        dvdx_s => duidxj_s(:,:,:,4); dvdy_s => duidxj_s(:,:,:,5); dvdz_s => duidxj_s(:,:,:,6);
        dwdx_s => duidxj_s(:,:,:,7); dwdy_s => duidxj_s(:,:,:,8); dwdz_s => duidxj_s(:,:,:,9);

        ! Compute the multiplying factors (thermo-shit)
        bambda = (four/three)*mu_x + bulk_x
        lambda = bulk_x - (two/three)*mu_x


        ! Step 1: Get tau_12  (dudy is destroyed)
        dudy =  dudy_x + dvdx_s
        dudy = mu_x*dudy
        !tauxyidz = 2

        ! Step 2: Get tau_13 (dudz is destroyed)
        dudz = dudz_x + dwdx_s
        dudz = mu_x*dudz
        !tauxzidx = 3

        ! Step 3: Get tau_23 (dvdz is destroyed)
        dvdz = dvdz_y +  dwdy_s
        dvdz = mu_y*dvdz
        !tauyzidx = 6

        ! Step 4: Get tau_11 (dvdx is destroyed)
        dvdx = bambda*dudx_s + lambda*(dvdy_x + dwdz_x)
        !tauxxidx = 4

        bambda = (four/three)*mu_y + bulk_y
        lambda = bulk_y - (two/three)*mu_y

        ! Step 5: Get tau_22 (dwdx is destroyed)
        dwdx = bambda*dvdy_s + lambda*(dudx_y + dwdz_y)
        !tauyyidx = 7
     
        bambda = (four/three)*mu_z + bulk_z
        lambda = bulk_z - (two/three)*mu_z

        ! Step 6: Get tau_33 (dwdy is destroyed)
        dwdy = bambda*dwdz_s + lambda*(dudx_z + dvdy_z)
        !tauzzidx = 8

        ! tau 21
        dvdy_x = mu_y*(dvdx_y + dudy_s)
        !tauyxid = 1
    
        ! tau 31
        dudx_y  = mu_z*( dwdx_z + dudz_s )

        ! tau 32 
        dudx_z = mu_z*(dvdz_s + dwdy_z)
        ! Done 
    end subroutine

    subroutine get_q(this,qx,qy,qz)
        use exits, only: nancheck
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(inout) :: qx,qy,qz

        ! integer :: i
        real(rkind), dimension(:,:,:), pointer :: tmp1_in_x, tmp2_in_x, tmp1_in_y, tmp1_in_z, tmp2_in_z
        type(derivatives), pointer :: der

        der => this%der

        tmp1_in_x => this%xbuf(:,:,:,1)
        tmp2_in_x => this%xbuf(:,:,:,2)

        tmp1_in_z => this%zbuf(:,:,:,1)
        tmp2_in_z => this%zbuf(:,:,:,2)

        tmp1_in_y => this%ybuf(:,:,:,1)

        ! Species enthalpy diffusion is computed earlier in SolidMixture

        ! Step 1: Get qy (dvdy is destroyed)
        call der%ddy(this%T,tmp1_in_y,this%y_bc(1),this%y_bc(2))
        qy = qy - this%kap*tmp1_in_y

        ! Step 2: Get qx (dudx is destroyed)
        call transpose_y_to_x(this%T,tmp1_in_x,this%decomp)
        call der%ddx(tmp1_in_x,tmp2_in_x,this%x_bc(1),this%x_bc(2))
        call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
        qx = qx - this%kap*tmp1_in_y

        ! Step 3: Get qz (dwdz is destroyed)
        call transpose_y_to_z(this%T,tmp1_in_z,this%decomp)
        call der%ddz(tmp1_in_z,tmp2_in_z,this%z_bc(1),this%z_bc(2))
        call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
        qz = qz - this%kap*tmp1_in_y

        ! Done
    end subroutine

    subroutine BicubicMetric(this)
        use timer,      only: tic, toc
        use exits,      only: message,nancheck,GracefulExit
        use reductions, only: P_MAXVAL, P_MINVAL,P_SUM
        use decomp_2d,  only: nrank
        use operators, only: divergence,gradient,gradFV_x,gradFV_y,gradFV_z,interpolateFV
        use constants,               only: pi
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx,drudy,drudz,drvdx,drvdy,drvdz,drwdx,drwdy,drwdz,dredx,dredy,dredz,dVFdx,dVFdy,dVFdz,dm1dx,dm1dy,dm1dz,dm2dx,dm2dy,dm2dz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudxy,drvdxy,drwdxy,dredxy,dVFdxy,dm1dxy,dm2dxy,tmp1,tmp2,tmp3,ru,rv,re,m1,m2,vf, mask
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3)      :: tmpint
        real(rkind), dimension(4,4)                               :: acoef,a,b,c
        real(rkind), dimension(4)                                 :: xvector,yvector,step1
        integer :: isub,i,j,k,l,imat,iter,ii,jj,kk
        real(rkind) :: MSETotal,MSEru,MSErv,MSErw,MSEre,MSEm1,MSEm2,MSEVF,dx,dy,sum_t, sum_t2,sum_logH, sum_tlogH
        real(rkind), dimension(6) :: global_heuristic
        logical :: file_exists
        character(len=clen) :: tempname, filename
 
        call this%gradient(this%Wcnsrv(:,:,:,mom_index),   drudx,drudy,drudz, this%x_bc,this%y_bc,this%z_bc) 
        call this%gradient(this%Wcnsrv(:,:,:,mom_index+1), drvdx,drvdy,drvdz,this%x_bc,this%y_bc,this%z_bc)
        call this%gradient(this%Wcnsrv(:,:,:,mom_index+2), drwdx,drwdy,drwdz,this%x_bc,this%y_bc,this%z_bc)  
        call this%gradient(this%Wcnsrv(:,:,:,TE_index),   dredx,dredy,dredz,this%x_bc,this%y_bc,this%z_bc)
        call this%gradient(this%mix%material(1)%VF,dVFdx,dVFdy,dVFdz,this%x_bc,this%y_bc,this%z_bc)
        call this%gradient(this%mix%material(1)%consrv(:,:,:,1),dm1dx,dm1dy,dm1dz,this%x_bc,this%y_bc,this%z_bc)
        call this%gradient(this%mix%material(2)%consrv(:,:,:,1),dm2dx,dm2dy,dm2dz,this%x_bc,this%y_bc,this%z_bc)

        call this%gradient(drudx,tmp1,drudxy,tmp2,this%x_bc,this%y_bc,this%z_bc)
        call this%gradient(drvdx,tmp1,drvdxy,tmp2,this%x_bc,this%y_bc,this%z_bc)
        call this%gradient(drwdx,tmp1,drwdxy,tmp2,this%x_bc,this%y_bc,this%z_bc)
        call this%gradient(dredx,tmp1,dredxy,tmp2,this%x_bc,this%y_bc,this%z_bc)
        call this%gradient(dm1dx,tmp1,dm1dxy,tmp2,this%x_bc,this%y_bc,this%z_bc)
        call this%gradient(dm2dx,tmp1,dm2dxy,tmp2,this%x_bc,this%y_bc,this%z_bc)
        call this%gradient(dVFdx,tmp1,dVFdxy,tmp2,this%x_bc,this%y_bc,this%z_bc)

        a(1,1) =  1.0; a(1,2) =  0.0; a(1,3) =  0.0; a(1,4) =  0.0;
        a(2,1) =  0.0; a(2,2) =  0.0; a(2,3) =  1.0; a(2,4) =  0.0;
        a(3,1) = -3.0; a(3,2) =  3.0; a(3,3) = -2.0; a(3,4) = -1.0;
        a(4,1) =  2.0; a(4,2) = -2.0; a(4,3) =  1.0; a(4,4) =  1.0;
 
        c(1,1) = 1.0; c(1,2) = 0.0; c(1,3) = -3.0; c(1,4) =  2.0;
        c(2,1) = 0.0; c(2,2) = 0.0; c(2,3) =  3.0; c(2,4) = -2.0;
        c(3,1) = 0.0; c(3,2) = 1.0; c(3,3) = -2.0; c(3,4) =  1.0;
        c(4,1) = 0.0; c(4,2) = 0.0; c(4,3) = -1.0; c(4,4) =  1.0;
        
        do i = 2,this%nxp-1
           do j = 2, this%nyp-1

               dx = sqrt( (this%x(i-1,j,1) - this%x(i,j-1,1) )**2 +(this%y(i+1,j,1) -this%y(i,j-1,1) )**2 )
               dy = dx
               b(1,1) =  this%Wcnsrv(i,j-1,1,mom_index); b(1,2) =  this%Wcnsrv(i+1,j,1,mom_index); b(1,3) = dy*drudy(i,j-1,1); b(1,4) =  dy*drudy(i+1,j,1);
               b(2,1) =  this%Wcnsrv(i-1,j,1,mom_index); b(2,2) =  this%Wcnsrv(i,j+1,1,mom_index); b(2,3) = dy*drudy(i-1,j,1); b(2,4) =  dy*drudy(i,j+1,1);
               b(3,1) =  dx*drudx(i,j-1,1); b(3,2) =  dx*drudx(i+1,j,1); b(3,3) = dx*dy*drudxy(i,j-1,1); b(3,4) = dx*dy*drudxy(i+1,j,1);
               b(4,1) =  dx*drudx(i-1,j,1); b(4,2) =  dx*drudx(i,j+1,1); b(4,3) = dx*dy*drudxy(i-1,j,1); b(4,4) = dx*dy*drudxy(i,j+1,1);

               acoef = MATMUL(a, MATMUL(b,c) )
               yvector(1) = 1.0; yvector(2) = (this%y(i,j,1) - this%y(i,j-1,1) ) /dy; ! (this%y(i+1,j) - this%y(i,j-1) );
               yvector(3) = ( (this%y(i,j,1) - this%y(i,j-1,1) ) / dy)**2.0 ! (this%y(i+1,j) - this%y(i,j-1) ) )**2.0;
               yvector(4) = ( (this%y(i,j,1) - this%y(i,j-1,1) ) / dy)**3.0 ! (this%y(i+1,j) - this%y(i,j-1) ) )**3.0;
               
               xvector(1) = 1.0; xvector(2) = (this%x(i,j,1) - this%x(i,j-1,1) ) /dx ! (this%x(i-1,j) - this%x(i,j-1) );
               xvector(3) = ( (this%x(i,j,1) - this%x(i,j-1,1) ) / dx)**2.0 ! (this%x(i-1,j) - this%x(i,j-1) ) )**2.0;
               xvector(4) = ( (this%x(i,j,1) - this%x(i,j-1,1) ) / dx)**3.0 ! (this%x(i-1,j) - this%x(i,j-1) ) )**3.0;

               do jj = 1,4

                   step1(jj) = DOT_PRODUCT(xvector,acoef(:,jj) )
                
               enddo               
               ru(i,j,1) = DOT_PRODUCT(step1,yvector)
            
           enddo
        enddo

        where(abs(this%mix%material(1)%VF) .LE. 1d-3 .OR. abs(this%mix%material(2)%VF) .LE. 1d-3 )

         mask = 1.0

        elsewhere

         mask = 0.0

        endwhere

 

       this%rhouHeur = ( this%Wcnsrv(:,:,:,mom_index) - ru )*mask
       this%rhovHeur = ru

  end subroutine

  subroutine FilDiffHeuristic(this)
        use timer,      only: tic, toc
        use exits,      only: message,nancheck,GracefulExit
        use reductions, only: P_MAXVAL, P_MINVAL,P_SUM
        use decomp_2d,  only: nrank
        use operators, only: divergence,gradient,gradFV_x,gradFV_y,gradFV_z,interpolateFV
        use constants,               only: pi
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx,drudy,drudz,drvdx,drvdy,drvdz,drwdx,drwdy,drwdz,dredx,dredy,dredz,dVFdx,dVFdy,dVFdz,dm1dx,dm1dy,dm1dz,dm2dx,dm2dy,dm2dz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx4,drudy4,drudz4,drvdx4,drvdy4,drvdz4,drwdx4,drwdy4,drwdz4,dredx4,dredy4,dredz4,dVFdx4,dVFdy4,dVFdz4,dm1dx4,dm1dy4,dm1dz4,dm2dx4,dm2dy4,dm2dz4, rufil,rvfil,rwfil,refil, m1fil, m2fil, vffil, mask
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3)      :: tmpint
        integer :: isub,i,j,k,l,imat,iter,ii,jj,kk
        real(rkind) :: MSETotal,MSEru,MSErv,MSErw,MSEre,MSEm1,MSEm2,MSEVF,sum_t,sum_t2,sum_logH,sum_tlogH
        real(rkind), dimension(7) :: global_heuristic, growthrate
        logical :: file_exists
        character(len=clen) :: tempname, filename

        where((abs(this%mix%material(1)%VF) .LE. 1d-4) .OR. abs(this%mix%material(2)%VF) .LE. 1d-4 )

            mask = 1 
  

        elsewhere

            mask = 0
        
        endwhere

        rufil = this%Wcnsrv(:,:,:,mom_index);   rvfil = this%Wcnsrv(:,:,:,mom_index+1);
        rwfil = this%Wcnsrv(:,:,:,mom_index+2); refil = this%Wcnsrv(:,:,:,TE_index);
        m1fil = this%mix%material(1)%consrv(:,:,:,1); m2fil = this%mix%material(2)%consrv(:,:,:,1);
        vffil = this%mix%material(1)%VF

        call this%filter(rufil, this%fil, 1,this%x_bc, this%y_bc, this%z_bc)
        call this%filter(rvfil, this%fil, 1,this%x_bc, this%y_bc, this%z_bc)
        call this%filter(rwfil, this%fil, 1,this%x_bc, this%y_bc, this%z_bc)
        call this%filter(refil, this%fil, 1,this%x_bc, this%y_bc, this%z_bc)
        call this%filter(m1fil, this%fil, 1,this%x_bc, this%y_bc, this%z_bc)
        call this%filter(m2fil, this%fil, 1,this%x_bc, this%y_bc, this%z_bc)
        call this%filter(vffil, this%fil, 1,this%x_bc, this%y_bc, this%z_bc)

        
        call this%secondder(this%derCD06,rufil - this%Wcnsrv(:,:,:,mom_index),drudx,drudy,drudz,[0,0],[0,0],[0,0])
        call this%secondder(this%derCD06,rvfil - this%Wcnsrv(:,:,:,mom_index+1),drvdx,drvdy,drvdz,[0,0],[0,0],[0,0])
        call this%secondder(this%derCD06,rwfil - this%Wcnsrv(:,:,:,mom_index+2),drwdx,drwdy,drwdz,[0,0],[0,0],[0,0])
        call this%secondder(this%derCD06,refil - this%Wcnsrv(:,:,:,TE_index),dredx,dredy,dredz,[0,0],[0,0],[0,0])
        call this%secondder(this%derCD06,m1fil - this%mix%material(1)%consrv(:,:,:,1),dm1dx,dm1dy,dm1dz,[0,0],[0,0],[0,0])
        call this%secondder(this%derCD06,m2fil - this%mix%material(2)%consrv(:,:,:,1),dm2dx,dm2dy,dm2dz,[0,0],[0,0],[0,0])
        call this%secondder(this%derCD06,vffil - this%mix%material(1)%VF,dVFdx,dVFdy,dVFdz,[0,0],[0,0],[0,0])

        call gradient(this%decomp,this%derCD06,rufil,drudx4,drudy4,drudz4,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD06,rvfil,drvdx4,drvdy4,drvdz4,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD06,rwfil,drwdx4,drwdy4,drwdz4,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD06,refil,dredx4,dredy4,dredz4,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD06,m1fil,dm1dx4,dm1dy4,dm1dz4,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD06,m2fil,dm2dx4,dm2dy4,dm2dz4,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD06,vffil,dVFdx4,dVFdy4,dVFdz4,[0,0],[0,0],[0,0]) 
 
        this%rhouHeur = (rufil - this%Wcnsrv(:,:,:,mom_index))
        this%rhovHeur = (rvfil - this%Wcnsrv(:,:,:,mom_index+1))
        this%rhowHeur = (rwfil - this%Wcnsrv(:,:,:,mom_index+2))
        this%rhoeHeur = (refil - this%Wcnsrv(:,:,:,TE_index))
        this%m1Heur   = (m1fil - this%mix%material(1)%consrv(:,:,:,1)) 
        this%m2Heur   = (m2fil- this%mix%material(2)%consrv(:,:,:,1))
        this%VFHeur   = (vffil - this%mix%material(1)%VF )

        MSEru = P_SUM(mask*(rufil - this%Wcnsrv(:,:,:,mom_index))**2 *this%dx*this%dy*this%dz)  
        MSErv = P_SUM(mask*(rvfil - this%Wcnsrv(:,:,:,mom_index+1))**2 *this%dx*this%dy*this%dz) 
        MSErw = P_SUM(mask*(rwfil - this%Wcnsrv(:,:,:,mom_index+2))**2 *this%dx*this%dy*this%dz) 
        MSEre = P_SUM(mask*(refil - this%Wcnsrv(:,:,:,TE_index))**2 *this%dx*this%dy*this%dz) 
        MSEm1 = P_SUM(mask*(m1fil - this%mix%material(1)%consrv(:,:,:,1))**2 *this%dx*this%dy*this%dz) 
        MSEm2 = P_SUM(mask*(m2fil - this%mix%material(2)%consrv(:,:,:,1))**2 *this%dx*this%dy*this%dz) 
        MSEVF = P_SUM(mask*(vffil - this%mix%material(1)%VF)**2 *this%dx*this%dy*this%dz) 

        do i = 1,9
        this%lamru(i) = this%lamru(i+1)
        this%lamrv(i) = this%lamrv(i+1)
        this%lamrw(i) = this%lamrw(i+1)
        this%lamre(i) = this%lamre(i+1)
        this%lamm1(i) = this%lamm1(i+1)
        this%lamm2(i) = this%lamm2(i+1)
        this%lamVF(i) = this%lamVF(i+1)
        this%lamtim(i) = this%lamtim(i+1)
        enddo
        this%lamru(10) = MSEru
        this%lamrv(10) = MSErv
        this%lamrw(10) = MSErw
        this%lamre(10) = MSEre
        this%lamm1(10) = MSEm1
        this%lamm2(10) = MSEm2
        this%lamVF(10) = MSEVF
        this%lamtim(10) = this%tsim

        if(  this%step .GE. 75 ) then

           if(this%lamru(1) .GT. 1d-40 .AND. this%lamru(4) .GT. 1d-40 ) then

             sum_t = 0; sum_t2 = 0; sum_logH = 0; sum_tlogH = 0;
 
             do i = 1,10

               sum_t = sum_t + this%lamtim(i)
               sum_t2 = sum_t2 + this%lamtim(i)*this%lamtim(i)
               sum_logH = sum_logH + LOG(this%lamru(i) )
               sum_tlogH = sum_tlogH + LOG(this%lamru(i) )*this%lamtim(i)

             enddo
             growthrate(1) = (10*sum_tlogH - sum_t*sum_logH) / (10*sum_t2 - sum_t*sum_t)  ! 1_rkind / (this%lamtim(4) - this%lamtim(1)) * LOG(this%lamru(4)/this%lamru(1))
         
           else

             growthrate(1) = 0
           
           endif

           if(this%lamrv(1) .GT. 1d-40 .AND. this%lamrv(4) .GT. 1d-40 ) then

               sum_t = 0; sum_t2 = 0; sum_logH = 0; sum_tlogH = 0;

             do i = 1,10

               sum_t = sum_t + this%lamtim(i)
               sum_t2 = sum_t2 + this%lamtim(i)*this%lamtim(i)
               sum_logH = sum_logH + LOG(this%lamrv(i) )
               sum_tlogH = sum_tlogH + LOG(this%lamrv(i) )*this%lamtim(i)

             enddo
             growthrate(2) = (10*sum_tlogH - sum_t*sum_logH) / (10*sum_t2 - sum_t*sum_t) 

           else

             growthrate(2) = 0

           endif

           if(this%lamrw(1) .GT. 1d-40 .AND. this%lamrw(4) .GT. 1d-40 ) then

               sum_t = 0; sum_t2 = 0; sum_logH = 0; sum_tlogH = 0;

               do i = 1,10

                 sum_t = sum_t + this%lamtim(i)
                 sum_t2 = sum_t2 + this%lamtim(i)*this%lamtim(i)
                 sum_logH = sum_logH + LOG(this%lamrw(i) )
                 sum_tlogH = sum_tlogH + LOG(this%lamrw(i) )*this%lamtim(i)

               enddo
               growthrate(3) = (10*sum_tlogH - sum_t*sum_logH) / (10*sum_t2 - sum_t*sum_t) 


           else

             growthrate(3) = 0

           endif

           if(this%lamre(1) .GT. 1d-32 .AND. this%lamre(4) .GT. 1d-32 ) then

               sum_t = 0; sum_t2 = 0; sum_logH = 0; sum_tlogH = 0;

               do i = 1,10

                 sum_t = sum_t + this%lamtim(i)
                 sum_t2 = sum_t2 + this%lamtim(i)*this%lamtim(i)
                 sum_logH = sum_logH + LOG(this%lamre(i) )
                 sum_tlogH = sum_tlogH + LOG(this%lamre(i) )*this%lamtim(i)

               enddo
               growthrate(4) = (10*sum_tlogH - sum_t*sum_logH) / (10*sum_t2 - sum_t*sum_t) 


           else

             growthrate(4) = 0

           endif

           if(this%lamm1(1) .GT. 1d-20 .AND. this%lamm1(4) .GT. 1d-20 ) then

             sum_t = 0; sum_t2 = 0; sum_logH = 0; sum_tlogH = 0;

             do i = 1,10

               sum_t = sum_t + this%lamtim(i)
               sum_t2 = sum_t2 + this%lamtim(i)*this%lamtim(i)
               sum_logH = sum_logH + LOG(this%lamm1(i) )
               sum_tlogH = sum_tlogH + LOG(this%lamm1(i) )*this%lamtim(i)

             enddo
             growthrate(5) = (10*sum_tlogH - sum_t*sum_logH) / (10*sum_t2 - sum_t*sum_t) 



           else

             growthrate(5)= 0

           endif


           if(this%lamm2(1) .GT. 1d-20 .AND. this%lamm2(4) .GT. 1d-20 ) then

             sum_t = 0; sum_t2 = 0; sum_logH = 0; sum_tlogH = 0;

             do i = 1,10

               sum_t = sum_t + this%lamtim(i)
               sum_t2 = sum_t2 + this%lamtim(i)*this%lamtim(i)
               sum_logH = sum_logH + LOG(this%lamm2(i) )
               sum_tlogH = sum_tlogH + LOG(this%lamm2(i) )*this%lamtim(i)

             enddo
             growthrate(6) = (10*sum_tlogH - sum_t*sum_logH) / (10*sum_t2 - sum_t*sum_t) 


           else

             growthrate(6) = 0

           endif

           if(this%lamVF(1) .GT. 1d-40 .AND. this%lamVF(4) .GT. 1d-40 ) then

             sum_t = 0; sum_t2 = 0; sum_logH = 0; sum_tlogH = 0;

             do i = 1,10

               sum_t = sum_t + this%lamtim(i)
               sum_t2 = sum_t2 + this%lamtim(i)*this%lamtim(i)
               sum_logH = sum_logH + LOG(this%lamVF(i) )
               sum_tlogH = sum_tlogH + LOG(this%lamVF(i) )*this%lamtim(i)

             enddo
             growthrate(7) = (10*sum_tlogH - sum_t*sum_logH) / (10*sum_t2 - sum_t*sum_t) 


           else

             growthrate(7) = 0

           endif

          
           MSETotal = 0

           do i = 1,7


             if( growthrate(i) .GE. 0 ) then

               MSETotal = MSETotal + growthrate(i)
               
 
             endif


           enddo

           MSETotal = MSETotal / 6_rkind


           if(MSEru .GE. 250 .OR. MSErv .GE. 250 .OR. MSEre .GE. 250 .OR. MSEm1 .GE. 250 .OR. MSEm2 .GE. 250 .OR. MSEVF .GE. 250 ) then

              this%Wcnsrv(:,:,:,mom_index) = rufil
              this%Wcnsrv(:,:,:,mom_index+1) = rvfil
              this%Wcnsrv(:,:,:,mom_index+2) = rwfil
              this%Wcnsrv(:,:,:,TE_index) = refil
              this%mix%material(1)%consrv(:,:,:,1) = m1fil
              this%mix%material(2)%consrv(:,:,:,1) = m2fil
              this%mix%material(1)%VF = VFfil
              this%mix%material(2)%VF = 1_rkind - this%mix%material(1)%VF
              this%numfil = this%numfil + 1
              this%stepfil = 0

              call this%get_primitive()
              call this%mix%equilibratePressure(this%rho,this%e, this%p)
              call this%post_bc()

              if (nrank == 0) then

               print *, " Filtered Again ", this%numfil
             
              endif

           endif
    
           

           




        else

!          this%stepfil = this%stepfil + 1

        endif

         this%stepfil = this%stepfil + 1

!        growthrate(1) = this%lamru(4); growthrate(2) = this%lamrv(4); growthrate(3) = this%lamrw(4); growthrate(4) = this%lamre(5);
!        growthrate(5) = this%lamm1(4); growthrate(6) = this%lamm2(4); growthrate(7) = this%lamVF(4);
        
       growthrate(3) = MSETotal 
       if(MOD(this%stepfil,10) .EQ. 0 ) then

            if (nrank == 0) then
               filename = 'heuristic_MSET.dat'

               ! Check if file exists to open in append or write mode
               inquire(file=filename, exist=file_exists)

               if (file_exists) then
                 open(unit=10, file=filename, status='old', position='append', action='write')
              else
                 open(unit=10, file=filename, status='new', action='write') 
                 write(10, '(A)') '# Timestep          Hru          Hrv           Hrw           Hre           Hm1           Hm2         H_alpha1'
              end if

              write(10, '(F12.5, 7(2X,E16.8))') this%tsim, (growthrate(i),i=1,7)
              close(10)

           end if




       endif



  end subroutine

  subroutine FilteringHeuristic(this)
        use timer,      only: tic, toc
        use exits,      only: message,nancheck,GracefulExit
        use reductions, only: P_MAXVAL, P_MINVAL,P_SUM
        use decomp_2d,  only: nrank
        use operators, only: divergence,gradient,gradFV_x,gradFV_y,gradFV_z,interpolateFV
        use constants,               only: pi
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx,drudy,drudz,drvdx,drvdy,drvdz,drwdx,drwdy,drwdz,dredx,dredy,dredz,dVFdx,dVFdy,dVFdz,dm1dx,dm1dy,dm1dz,dm2dx,dm2dy,dm2dz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx4,drudy4,drudz4,drvdx4,drvdy4,drvdz4,drwdx4,drwdy4,drwdz4,dredx4,dredy4,dredz4,dVFdx4,dVFdy4,dVFdz4,dm1dx4,dm1dy4,dm1dz4,dm2dx4,dm2dy4,dm2dz4       
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3)      :: tmpint
        integer :: isub,i,j,k,l,imat,iter,ii,jj,kk
        real(rkind) :: MSETotal,MSEru,MSErv,MSErw,MSEre,MSEm1,MSEm2,MSEVF
        real(rkind), dimension(6) :: global_heuristic
        logical :: file_exists
        character(len=clen) :: tempname, filename
            call interpolateFV(this%decomp,this%interpMid,this%Wcnsrv(:,:,:,mom_index),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStagg,tmpint(:,:,:,1),drudx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStagg,tmpint(:,:,:,2),drudy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStagg,tmpint(:,:,:,3),drudz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid,this%Wcnsrv(:,:,:,mom_index+1),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStagg,tmpint(:,:,:,1),drvdx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStagg,tmpint(:,:,:,2),drvdy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStagg,tmpint(:,:,:,3),drvdz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid,this%Wcnsrv(:,:,:,mom_index+3),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStagg,tmpint(:,:,:,1),drwdx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStagg,tmpint(:,:,:,2),drwdy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStagg,tmpint(:,:,:,3),drwdz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid,this%Wcnsrv(:,:,:,TE_index),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStagg,tmpint(:,:,:,1),dredx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStagg,tmpint(:,:,:,2),dredy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStagg,tmpint(:,:,:,3),dredz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid,this%mix%material(1)%consrv(:,:,:,1),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStagg,tmpint(:,:,:,1),dm1dx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStagg,tmpint(:,:,:,2),dm1dy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStagg,tmpint(:,:,:,3),dm1dz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid,this%mix%material(2)%consrv(:,:,:,1),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStagg,tmpint(:,:,:,1),dm2dx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStagg,tmpint(:,:,:,2),dm2dy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStagg,tmpint(:,:,:,3),dm2dz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid,this%mix%material(1)%VF,tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStagg,tmpint(:,:,:,1),dVFdx,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStagg,tmpint(:,:,:,2),dVFdy,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)            
            call gradFV_z(this%decomp,this%derStagg,tmpint(:,:,:,3),dVFdz,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            !!!!!!!!!!! 4th Order             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call interpolateFV(this%decomp,this%interpMid04,this%Wcnsrv(:,:,:,mom_index),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStaggd04,tmpint(:,:,:,1),drudx4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStaggd04,tmpint(:,:,:,2),drudy4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStaggd04,tmpint(:,:,:,3),drudz4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid04,this%Wcnsrv(:,:,:,mom_index+1),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStaggd04,tmpint(:,:,:,1),drvdx4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStaggd04,tmpint(:,:,:,2),drvdy4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStaggd04,tmpint(:,:,:,3),drvdz4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid04,this%Wcnsrv(:,:,:,mom_index+3),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStaggd04,tmpint(:,:,:,1),drwdx4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStaggd04,tmpint(:,:,:,2),drwdy4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStaggd04,tmpint(:,:,:,3),drwdz4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid04,this%Wcnsrv(:,:,:,TE_index),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStaggd04,tmpint(:,:,:,1),dredx4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStaggd04,tmpint(:,:,:,2),dredy4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStaggd04,tmpint(:,:,:,3),dredz4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid04,this%mix%material(1)%consrv(:,:,:,1),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStaggd04,tmpint(:,:,:,1),dm1dx4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStaggd04,tmpint(:,:,:,2),dm1dy4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStaggd04,tmpint(:,:,:,3),dm1dz4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid04,this%mix%material(2)%consrv(:,:,:,1),tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStaggd04,tmpint(:,:,:,1),dm2dx4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStaggd04,tmpint(:,:,:,2),dm2dy4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStaggd04,tmpint(:,:,:,3),dm2dz4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            call interpolateFV(this%decomp,this%interpMid04,this%mix%material(1)%VF,tmpint,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_x(this%decomp,this%derStaggd04,tmpint(:,:,:,1),dVFdx4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_y(this%decomp,this%derStaggd04,tmpint(:,:,:,2),dVFdy4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)
            call gradFV_z(this%decomp,this%derStaggd04,tmpint(:,:,:,3),dVFdz4,this%periodicx,this%periodicy,this%periodicz,this%x_bc,this%y_bc,this%z_bc)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Take Difference             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            drudx4 = drudx-drudx4; drudy4 = drudy-drudy4; drudz4 = drudz-drudz4;
            drvdx4 = drvdx-drvdx4; drvdy4 = drvdy-drvdy4; drvdz4 = drvdz-drvdz4;
            drwdx4 = drwdx-drwdx4; drwdy4 = drwdy-drwdy4; drwdz4 = drwdz-drwdz4;
            dredx4 = dredx-dredx4; dredy4 = dredy-dredy4; dredz4 = dredz-dredz4;
            dm1dx4 = dm1dx-dm1dx4; dm1dy4 = dm1dy-dm1dy4; dm1dz4 = dm1dz-dm1dz4;
            dm2dx4 = dm2dx-dm2dx4; dm2dy4 = dm2dy-dm2dy4; dm2dz4 = dm2dz-dm2dz4;
            dVFdx4 = dVFdx-dVFdx4; dVFdy4 = dVFdy-dVFdy4; dVFdz4 = dVFdz-dVFdz4;

           
            MSEru  =( P_SUM( (this%dx*this%dy*this%dz)*(drudx4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)* (drudy4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)* (drudz4)**2 ) ) &
                   /(1d-20 + P_SUM( this%dx*this%dy*this%dz*( drudx**2 + drudy**2 + drudz**2)*0.015) )                



! / ( 1d-14 + ((this%dx*this%dy*this%dz)*P_MAXVAL(drudx)**2 + (this%dx*this%dy*this%dz)*P_MAXVAL(drudy)**2 + (this%dx*this%dy*this%dz)*P_MAXVAL(drudx)**2 ) )

            MSErv  =( P_SUM( (this%dx*this%dy*this%dz)*(drvdx4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)* (drvdy4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)* (drvdz4)**2 ) ) &
                   /( 1d-20 + P_SUM( this%dx*this%dy*this%dz*( drvdx**2 + drvdy**2 + drvdz**2)*0.015 ) )

!                     / ( 1d-14 + ((this%dx*this%dy*this%dz)*P_MAXVAL(drvdx)**2 + (this%dx*this%dy*this%dz)*P_MAXVAL(drvdy)**2 + (this%dx*this%dy*this%dz)*P_MAXVAL(drvdx)**2 ) )
            MSErw  =( P_SUM( (this%dx*this%dy*this%dz)*(drwdx4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)* (drwdy4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)* (drwdz4)**2 ) ) &
                   / ( 1d-20 + P_SUM( this%dx*this%dy*this%dz*( drwdx**2 + drwdy**2 + drwdz**2)*1e-10 ) )
!                     / ( 1d-14 + ((this%dx*this%dy*this%dz)*P_MAXVAL(drwdx)**2 + (this%dx*this%dy*this%dz)*P_MAXVAL(drwdy)**2 + (this%dx*this%dy*this%dz)*P_MAXVAL(drwdx)**2 ) )
   
            MSEre  =( P_SUM( (this%dx*this%dy*this%dz)*(dredx4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)* (dredy4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)*(dredz4)**2 ) )  &
                   / (1d-20 + P_SUM( this%dx*this%dy*this%dz*( dredx**2 + dredy**2 + dredz**2)*2.5 ) )
!                     / ( 1d-14 + ((this%dx*this%dy*this%dz)*P_MAXVAL(dredx)**2 + (this%dx*this%dy*this%dz)*P_MAXVAL(dredy)**2 + (this%dx*this%dy*this%dz)*P_MAXVAL(dredx)**2 ) )

            MSEm1  =( P_SUM( (this%dx*this%dy*this%dz)*(dm1dx4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)* (dm1dy4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)*(dm1dz4)**2 ) ) &
                   /(1d-20 + P_SUM( this%dx*this%dy*this%dz*( dm1dx**2 + dm1dy**2 + dm1dz**2) ) )

            MSEm2  =( P_SUM( (this%dx*this%dy*this%dz)*(dm2dx4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)* (dm2dy4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)*(dm2dz4)**2 ) )  &
                   /(1d-20 + P_SUM( this%dx*this%dy*this%dz*( dm2dx**2 + dm2dy**2 + dm2dz**2) ) )

            MSEVF  =( P_SUM( (this%dx*this%dy*this%dz)*(dVFdx4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)* (dVFdy4)**2 ) + P_SUM( (this%dx*this%dy*this%dz)*(dVFdz4)**2 ) )  &
                   /(1d-20 + P_SUM( this%dx*this%dy*this%dz*( dVFdx**2 + dVFdy**2 + dVFdz**2) ) )

            this%rhouHeur = (this%dx*this%dy*this%dz)*(drudx4)**2 + (this%dx*this%dy*this%dz)* (drudy4)**2 
            this%rhovHeur = (this%dx*this%dy*this%dz)*(drvdx4)**2 + (this%dx*this%dy*this%dz)* (drvdy4)**2
            this%rhoeHeur = (this%dx*this%dy*this%dz)*(dredx4)**2 + (this%dx*this%dy*this%dz)* (dredy4)**2 
            this%m1Heur   = (this%dx*this%dy*this%dz)*(dm1dx4)**2 + (this%dx*this%dy*this%dz)* (dm1dy4)**2 
            this%m2Heur   = (this%dx*this%dy*this%dz)*(dm2dx4)**2 + (this%dx*this%dy*this%dz)* (dm2dy4)**2 
            this%VFHeur   = (this%dx*this%dy*this%dz)*(dVFdx4)**2 + (this%dx*this%dy*this%dz)* (dVFdy4)**2 
            MSETotal = MSEru + MSErv  + MSEre + MSEm2 + MSEm1 + MSEVF

            global_heuristic(1) = MSEru; global_heuristic(2) = MSErv; global_heuristic(3) = MSErw;
            global_heuristic(4) = MSEre; global_heuristic(5) = MSEm1; global_heuristic(6) = MSEm2;
            global_heuristic(7) = MSEVF; global_heuristic(8) = MSETotal;

            if (nrank == 0) then
               filename = 'heuristic_logNorm.dat'

               ! Check if file exists to open in append or write mode
               inquire(file=filename, exist=file_exists)
    
               if (file_exists) then
                 open(unit=10, file=filename, status='old', position='append', action='write')
              else
                 open(unit=10, file=filename, status='new', action='write')
                 write(10, '(A)') '# Timestep          Hru          Hrv          Hrw           Hre           Hm1           Hm2         H_alpha1        H_total'
              end if

              write(10, '(F12.5, 8(2X,E16.8))') this%tsim, (global_heuristic(i), i=1,8)
              close(10)

           end if
     end subroutine

     subroutine FilteringHeuristicHighOrder(this)
        use timer,      only: tic, toc
        use exits,      only: message,nancheck,GracefulExit
        use reductions, only: P_MAXVAL, P_MINVAL,P_SUM
        use decomp_2d,  only: nrank
        use operators, only: divergence,gradient,gradFV_x,gradFV_y,gradFV_z,interpolateFV,laplacian
        use constants,               only: pi
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx,drudy,drudz,drvdx,drvdy,drvdz,drwdx,drwdy,drwdz,dredx,dredy,dredz,dVFdx,dVFdy,dVFdz,dm1dx,dm1dy,dm1dz,dm2dx,dm2dy,dm2dz,tmp1,tmp2,tmp3
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx4,drudy4,drudz4,drvdx4,drvdy4,drvdz4,drwdx4,drwdy4,drwdz4,dredx4,dredy4,dredz4,dVFdx4,dVFdy4,dVFdz4,dm1dx4,dm1dy4,dm1dz4,dm2dx4,dm2dy4,dm2dz4, mask
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3)      :: tmpint
        integer :: isub,i,j,k,l,imat,iter,ii,jj,kk
        real(rkind) :: MSETotal,MSEru,MSErv,MSErw,MSEre,MSEm1,MSEm2,MSEVF
        real(rkind), dimension(6) :: global_heuristic

        logical :: file_exists
        character(len=clen) :: tempname, filename


        call this%fourthder(this%derCD06,this%Wcnsrv(:,:,:,mom_index),drudx,drudy,drudz,[0,0],[0,0],[0,0])
        call this%fourthder(this%derCD04,this%Wcnsrv(:,:,:,mom_index),drudx4,drudy4,drudz4,[0,0],[0,0],[0,0])
       
        call this%fourthder(this%derCD06,this%Wcnsrv(:,:,:,mom_index+1),drvdx,drvdy,drvdz,[0,0],[0,0],[0,0])
        call this%fourthder(this%derCD04,this%Wcnsrv(:,:,:,mom_index+1),drvdx4,drvdy4,drvdz4,[0,0],[0,0],[0,0]) 
 
        call this%fourthder(this%derCD06,this%Wcnsrv(:,:,:,mom_index+2),drwdx,drwdy,drwdz,[0,0],[0,0],[0,0])
        call this%fourthder(this%derCD04,this%Wcnsrv(:,:,:,mom_index+2),drwdx4,drwdy4,drwdz4,[0,0],[0,0],[0,0])

        call this%fourthder(this%derCD06,this%Wcnsrv(:,:,:,TE_index),dredx,dredy,dredz,[0,0],[0,0],[0,0])
        call this%fourthder(this%derCD04,this%Wcnsrv(:,:,:,TE_index),dredx4,dredy4,dredz4,[0,0],[0,0],[0,0])
!        call this%fourthder(this%derCD06,this%p,dredx,dredy,dredz,[0,0],[0,0],[0,0])
!        call this%fourthder(this%derCD04,this%p,dredx4,dredy4,dredz4,[0,0],[0,0],[0,0])
        call this%fourthder(this%derCD06,this%mix%material(1)%consrv(:,:,:,1),dm1dx,dm1dy,dm1dz,[0,0],[0,0],[0,0])
        call this%fourthder(this%derCD04,this%mix%material(1)%consrv(:,:,:,1),dm1dx4,dm1dy4,dm1dz4,[0,0],[0,0],[0,0])

        call this%fourthder(this%derCD06,this%mix%material(2)%consrv(:,:,:,1),dm2dx,dm2dy,dm2dz,[0,0],[0,0],[0,0])
        call this%fourthder(this%derCD04,this%mix%material(2)%consrv(:,:,:,1),dm2dx4,dm2dy4,dm2dz4,[0,0],[0,0],[0,0])

        call this%fourthder(this%derCD06,this%mix%material(1)%VF,dVFdx,dVFdy,dVFdz,[0,0],[0,0],[0,0])
        call this%fourthder(this%derCD04,this%mix%material(1)%VF,dVFdx4,dVFdy4,dVFdz4,[0,0],[0,0],[0,0])

        drudx4 = drudx-drudx4; drudy4 = drudy-drudy4; drudz4 = drudz-drudz4;
        drvdx4 = drvdx-drvdx4; drvdy4 = drvdy-drvdy4; drvdz4 = drvdz-drvdz4;
        drwdx4 = drwdx-drwdx4; drwdy4 = drwdy-drwdy4; drwdz4 = drwdz-drwdz4;
        dredx4 = dredx-dredx4; dredy4 = dredy-dredy4; dredz4 = dredz-dredz4;
        dm1dx4 = dm1dx-dm1dx4; dm1dy4 = dm1dy-dm1dy4; dm1dz4 = dm1dz-dm1dz4;
        dm2dx4 = dm2dx-dm2dx4; dm2dy4 = dm2dy-dm2dy4; dm2dz4 = dm2dz-dm2dz4;
        dVFdx4 = dVFdx-dVFdx4; dVFdy4 = dVFdy-dVFdy4; dVFdz4 = dVFdz-dVFdz4;
 
        where(abs(this%mix%material(1)%VF) .LE. 1d-3 .OR. abs(this%mix%material(2)%VF) .LE. 1d-3 )

         mask = 1.0

        elsewhere

         mask = 0.0

        endwhere

        this%rhouHeur = 1/(1 + sqrt( drudx**2 + drudy**2) ) * ( (drudx4)**2 + (drudy4)**2 )*mask
        this%rhovHeur = 1/(1 + sqrt( drvdx**2 + drvdy**2) ) * ( (drvdx4)**2 + (drvdy4)**2 )*mask
        this%rhoeHeur = 1/(1 + sqrt( dredx**2 + dredy**2) ) * ( (dredx4)**2 + (dredy4)**2 )*mask
        this%m1Heur   = 1/(1 + sqrt( dm1dx**2 + dm1dy**2) ) * ( (dm1dx4)**2 + (dm1dy4)**2 )*mask
        this%m2Heur   = 1/(1 + sqrt( dm2dx**2 + dm2dy**2) ) * ( (dm2dx4)**2 + (dm2dy4)**2 )*mask
        this%VFHeur   = 1/(1 + sqrt( dVFdx**2 + dVFdy**2) ) * ( (dVFdx4)**2 + (dVFdy4)**2 )*mask


     end subroutine

     subroutine LocalDiffHeuristic(this)
        use timer,      only: tic, toc
        use exits,      only: message,nancheck,GracefulExit
        use reductions, only: P_MAXVAL, P_MINVAL,P_SUM
        use decomp_2d,  only: nrank
        use operators, only: divergence,gradient,gradFV_x,gradFV_y,gradFV_z,interpolateFV,laplacian
        use constants,               only: pi
        class(sgrid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx,drudy,drudz,drvdx,drvdy,drvdz,drwdx,drwdy,drwdz,dredx,dredy,dredz,dVFdx,dVFdy,dVFdz,dm1dx,dm1dy,dm1dz,dm2dx,dm2dy,dm2dz,tmp1,tmp2,tmp3
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)        :: drudx4,drudy4,drudz4,drvdx4,drvdy4,drvdz4,drwdx4,drwdy4,drwdz4,dredx4,dredy4,dredz4,dVFdx4,dVFdy4,dVFdz4,dm1dx4,dm1dy4,dm1dz4,dm2dx4,dm2dy4,dm2dz4
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3)      :: tmpint
        integer :: isub,i,j,k,l,imat,iter,ii,jj,kk
        real(rkind) :: MSETotal,MSEru,MSErv,MSErw,MSEre,MSEm1,MSEm2,MSEVF
        real(rkind), dimension(6) :: global_heuristic

        logical :: file_exists
        character(len=clen) :: tempname, filename


        call gradient(this%decomp,this%derCD06,this%Wcnsrv(:,:,:,mom_index),drudx,drudy,drudz,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD04,this%Wcnsrv(:,:,:,mom_index),drudx4,drudy4,drudz4,[0,0],[0,0],[0,0])

        call gradient(this%decomp,this%derCD06,this%Wcnsrv(:,:,:,mom_index+1),drvdx,drvdy,drvdz,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD04,this%Wcnsrv(:,:,:,mom_index+1),drvdx4,drvdy4,drvdz4,[0,0],[0,0],[0,0])

        call gradient(this%decomp,this%derCD06,this%Wcnsrv(:,:,:,mom_index+2),drwdx,drwdy,drwdz,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD04,this%Wcnsrv(:,:,:,mom_index+2),drwdx4,drwdy4,drwdz4,[0,0],[0,0],[0,0])

        call gradient(this%decomp,this%derCD06,this%Wcnsrv(:,:,:,TE_index),dredx,dredy,dredz,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD04,this%Wcnsrv(:,:,:,TE_index),dredx4,dredy4,dredz4,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD06,this%mix%material(1)%consrv(:,:,:,1),dm1dx,dm1dy,dm1dz,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD04,this%mix%material(1)%consrv(:,:,:,1),dm1dx4,dm1dy4,dm1dz4,[0,0],[0,0],[0,0])

        call gradient(this%decomp,this%derCD06,this%mix%material(2)%consrv(:,:,:,1),dm2dx,dm2dy,dm2dz,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD04,this%mix%material(2)%consrv(:,:,:,1),dm2dx4,dm2dy4,dm2dz4,[0,0],[0,0],[0,0])

        call gradient(this%decomp,this%derCD06,this%mix%material(1)%VF,dVFdx,dVFdy,dVFdz,[0,0],[0,0],[0,0])
        call gradient(this%decomp,this%derCD04,this%mix%material(1)%VF,dVFdx4,dVFdy4,dVFdz4,[0,0],[0,0],[0,0])

        this%rhouHeur = 0
        this%rhovHeur = 0
        this%rhowHeur = 0
        this%rhoeHeur = 0
        this%m1Heur   = 0 
        this%m2Heur   = 0
        this%VFHeur   = 0

        do i = 2,this%nxp-1
           do j = 2, this%nxp-1


              this%rhouHeur(i,j,1) = (drudx(i,j,1) - drudx(i,j+1,1) )**2 + (drudx(i,j,1) - drudx(i,j-1,1) )**2  &
                                       + (drudx(i,j,1) - drudx(i+1,j,1) )**2 + (drudx(i,j,1) - drudx(i-1,j,1) )**2  &
                                       + (drudy(i,j,1) - drudy(i,j+1,1) )**2 + (drudy(i,j,1) - drudy(i,j-1,1) )**2  &
                                       + (drudy(i,j,1) - drudy(i+1,j,1) )**2 + (drudy(i,j,1) - drudy(i-1,j,1) )**2   

              this%rhovHeur(i,j,1) = (drvdx(i,j,1) - drvdx(i,j+1,1) )**2 + (drvdx(i,j,1) - drvdx(i,j-1,1) )**2  &
                                       + (drvdx(i,j,1) - drvdx(i+1,j,1) )**2 + (drvdx(i,j,1) - drvdx(i-1,j,1) )**2  &
                                       + (drvdy(i,j,1) - drvdy(i,j+1,1) )**2 + (drvdy(i,j,1) - drvdy(i,j-1,1) )**2  &
                                       + (drvdy(i,j,1) - drvdy(i+1,j,1) )**2 + (drvdy(i,j,1) - drvdy(i-1,j,1) )**2 

              this%rhowHeur(i,j,1) = (drwdx(i,j,1) - drwdx(i,j+1,1) )**2 + (drwdx(i,j,1) - drwdx(i,j-1,1) )**2  &
                                       + (drwdx(i,j,1) - drwdx(i+1,j,1) )**2 + (drwdx(i,j,1) - drwdx(i-1,j,1) )**2  &
                                       + (drwdy(i,j,1) - drwdy(i,j+1,1) )**2 + (drwdy(i,j,1) - drwdy(i,j-1,1) )**2  &
                                       + (drwdy(i,j,1) - drwdy(i+1,j,1) )**2 + (drwdy(i,j,1) - drwdy(i-1,j,1) )**2 

              this%rhoeHeur(i,j,1) =  (dredx(i,j,1) - dredx(i,j+1,1) )**2 + (dredx(i,j,1) - dredx(i,j-1,1) )**2  &
                                       + (dredx(i,j,1) - dredx(i+1,j,1) )**2 + (dredx(i,j,1) - dredx(i-1,j,1) )**2  &
                                       + (dredy(i,j,1) - dredy(i,j+1,1) )**2 + (dredy(i,j,1) - dredy(i,j-1,1) )**2  &
                                       + (dredy(i,j,1) - dredy(i+1,j,1) )**2 + (dredy(i,j,1) - dredy(i-1,j,1) )**2 

              this%m1Heur(i,j,1) =  (dm1dx(i,j,1) - dm1dx(i,j+1,1) )**2 + (dm1dx(i,j,1) - dm1dx(i,j-1,1) )**2  &
                                       + (dm1dx(i,j,1) - dm1dx(i+1,j,1) )**2 + (dm1dx(i,j,1) - dm1dx(i-1,j,1) )**2  &
                                       + (dm1dy(i,j,1) - dm1dy(i,j+1,1) )**2 + (dm1dy(i,j,1) - dm1dy(i,j-1,1) )**2  &
                                       + (dm1dy(i,j,1) - dm1dy(i+1,j,1) )**2 + (dm1dy(i,j,1) - dm1dy(i-1,j,1) )**2 

              this%m2Heur(i,j,1) = (dm2dx(i,j,1) - dm2dx(i,j+1,1) )**2 + (dm2dx(i,j,1) - dm2dx(i,j-1,1) )**2  &
                                       + (dm2dx(i,j,1) - dm2dx(i+1,j,1) )**2 + (dm2dx(i,j,1) - dm2dx(i-1,j,1) )**2  &
                                       + (dm2dy(i,j,1) - dm2dy(i,j+1,1) )**2 + (dm2dy(i,j,1) - dm2dy(i,j-1,1) )**2  &
                                       + (dm2dy(i,j,1) - dm2dy(i+1,j,1) )**2 + (dm2dy(i,j,1) - dm2dy(i-1,j,1) )**2 

              this%VFHeur(i,j,1) = (dVFdx(i,j,1) - dVFdx(i,j+1,1) )**2 + (dVFdx(i,j,1) - dVFdx(i,j-1,1) )**2  &
                                       + (dVFdx(i,j,1) - dVFdx(i+1,j,1) )**2 + (dVFdx(i,j,1) - dVFdx(i-1,j,1) )**2  &
                                       + (dVFdy(i,j,1) - dVFdy(i,j+1,1) )**2 + (dVFdy(i,j,1) - dVFdy(i,j-1,1) )**2  &
                                       + (dVFdy(i,j,1) - dVFdy(i+1,j,1) )**2 + (dVFdy(i,j,1) - dVFdy(i-1,j,1) )**2 


           enddo
        enddo


     end subroutine

     subroutine readRestartFile(this, tid, rid)
        use decomp_2d,  only: nrank
        use decomp_2d_io
        use mpi
        use exits, only: message
        use kind_parameters, only: mpirkind
        class(sgrid), intent(inout) :: this
        integer, intent(in) :: tid, rid
        character(len=clen) :: tempname, fname
        integer :: ierr, fid

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_u.",tid
        fname = this%inputdir(:len_trim(this%inputdir))//"/"//trim(tempname)
        call decomp_2d_read_one(2,this%u,fname, this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_v.",tid
        fname = this%inputdir(:len_trim(this%inputdir))//"/"//trim(tempname)
        call decomp_2d_read_one(2,this%v,fname, this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_w.",tid
        fname = this%inputdir(:len_trim(this%inputdir))//"/"//trim(tempname)
        call decomp_2d_read_one(2,this%w,fname, this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_VF.",tid
        fname = this%inputdir(:len_trim(this%inputdir))//"/"//trim(tempname)
        call decomp_2d_read_one(2,this%mix%material(1)%VF,fname, this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_Ys.",tid
        fname = this%inputdir(:len_trim(this%inputdir))//"/"//trim(tempname)
        call decomp_2d_read_one(2,this%mix%material(1)%Ys,fname, this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_rho.",tid
        fname = this%inputdir(:len_trim(this%inputdir))//"/"//trim(tempname)
        call decomp_2d_read_one(2,this%rho,fname, this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_p.",tid
        fname = this%inputdir(:len_trim(this%inputdir))//"/"//trim(tempname)
        call decomp_2d_read_one(2,this%p,fname, this%decomp)
        this%mix%material(1)%p = this%p 
        this%mix%material(2)%p = this%p
        if (nrank == 0) then
            write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",rid, "_info.",tid
            fname = this%inputdir(:len_trim(this%inputdir))//"/"//trim(tempname)
            fid = 10
            open(unit=fid,file=trim(fname),status="old",action="read")
            read (fid, "(100g15.5)")  this%tsim
            close(fid)
        end if

    

        call mpi_barrier(mpi_comm_world, ierr)
        call mpi_bcast(this%tsim,1,mpirkind,0,mpi_comm_world,ierr)
        call mpi_barrier(mpi_comm_world, ierr)
        call message("================= RESTART FILE USED ======================")
        call message(0, "Simulation Time at restart:", this%tsim)
        call message("=================================== ======================")

    !   if(nrank == 0) then
    !      print *, " y ", this%w(1,:,1)
    !    endif


    end subroutine

    subroutine dumpRestartFile(this)
        use decomp_2d,  only: nrank
        use decomp_2d_io
        use mpi
        use exits, only: message
        class(sgrid), intent(inout) :: this
        character(len=512) :: tempname, fname
        integer :: ierr, rank

        call MPI_COMM_RANK(mpi_comm_world,rank,ierr)
     !   this%v = rank
     !   this%w = this%y
     !   this%p = this%x

     !   if(nrank == 0) then
     !     print *, " y ", this%w(1,:,1)
     !   endif
        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_u.",this%step
        fname = this%outputdir(:len_trim(this%outputdir))//"/"//trim(tempname)
        call decomp_2d_write_one(2,this%u,trim(fname), this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_v.",this%step
        fname = this%outputdir(:len_trim(this%outputdir))//"/"//trim(tempname)
        call decomp_2d_write_one(2,this%v,trim(fname), this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_w.",this%step
        fname = this%outputdir(:len_trim(this%outputdir))//"/"//trim(tempname)
        call decomp_2d_write_one(2,this%w,trim(fname), this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_VF.",this%step
        fname = this%outputdir(:len_trim(this%outputdir))//"/"//trim(tempname)
        call decomp_2d_write_one(2,this%mix%material(1)%VF,trim(fname), this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_Ys.",this%step
        fname = this%outputdir(:len_trim(this%outputdir))//"/"//trim(tempname)
        call decomp_2d_write_one(2,this%mix%material(1)%Ys,trim(fname), this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_rho.",this%step
        fname = this%outputdir(:len_trim(this%outputdir))//"/"//trim(tempname)
        call decomp_2d_write_one(2,this%rho,trim(fname), this%decomp)

        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_p.",this%step
        fname = this%outputdir(:len_trim(this%outputdir))//"/"//trim(tempname)
        call decomp_2d_write_one(2,this%p,trim(fname), this%decomp)

        if (nrank == 0) then
            write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",this%runID, "_info.",this%step
            fname = this%outputdir(:len_trim(this%outputdir))//"/"//trim(tempname)
            OPEN(UNIT=10, FILE=trim(fname))
            write(10,"(100g15.5)") this%tsim
            close(10)
        end if

        if (nrank == 0) then
           print *, "DEBUG: runID=", this%runID, " step=", this%step
           print *, "DEBUG filename=", trim(tempname)
       endif

       if (nrank == 0) then
            print *, "Passing fname to decomp_2d_write_one: ", trim(fname)
       endif
       

        call mpi_barrier(mpi_comm_world, ierr)
        call message(1, "Just Dumped a RESTART file")

    end subroutine


    ! NOTE: If you want to dump an edge field, you need to call in dumpFullField
    ! routine with this%gpE passed in as the 3rd argument. If it's a cell field,
    ! then you don't need to pass in any gp since the default gp is this%gpC
!   subroutine dumpFullField(this,arr,label,gp2use)
!       use decomp_2d_io
!       use mpi
!       use exits, only: message
!       class(igrid), intent(in) :: this
!       character(len=clen) :: tempname, fname
!       real(rkind), dimension(:,:,:), intent(in) :: arr
!       character(len=4), intent(in) :: label
!       type(decomp_info), intent(in), optional :: gp2use

!        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",this%step,".out"
!        fname = this%outputdir(:len_trim(this%OutputDir))//"/"//trim(tempname)
!        if (present(gp2use)) then
!           call decomp_2d_write_one(1,arr,fname,gp2use)
!        else
!           call decomp_2d_write_one(1,arr,fname,this%gpC)
!        end if

!   end subroutine


!   subroutine dumpVisualizationInfo(this)
!       class(sgrid), intent(in) :: this
!       character(len=clen) :: tempname, fname


!     if (nrank == 0) then
!           write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_","info","_t",this%step,".out"
!           fname = this%outputdir(:len_trim(this%OutputDir))//"/"//trim(tempname)
!           OPEN(UNIT=10, FILE=trim(fname))
!           write(10,"(100g17.9)") this%tsim
!           write(10,"(100g17.9)") real(this%nx,rkind)
!           write(10,"(100g17.9)") real(this%ny,rkind)
!           write(10,"(100g17.9)") real(this%nz,rkind)
!           close(10)
!       end if
!   end subroutine
    !! STATISTICS !!
end module 
