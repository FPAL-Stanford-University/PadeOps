module ShearLayer4Mode_data
    use kind_parameters,  only: rkind
    use constants,        only: one,two,eight,three,six,sixth,zero, pi
    use FiltersMod,       only: filters
    use DerivativesMod,   only: derivatives
    use DerivativesStaggeredMod, only: derivativesStagg
    use InterpolatorsMod,        only: interpolators
    use mpi 
    implicit none

    real(rkind) :: p_infty = one, Rgas = one, gamma = 1.4_rkind, mu = 0._rkind, rho_0 = one, p_amb = 0.1_rkind
    real(rkind) :: p_infty_2 = one, Rgas_2 = one, gamma_2 = 1.4_rkind, mu_2 = 0._rkind, rho_0_2 = one, eta_det_ge = one,eta_det_ge_2 = one, eta_det_gp = one,eta_det_gp_2 = one, eta_det_gt = one,eta_det_gt_2 = one,diff_c_ge = one,diff_c_ge_2 = one, diff_c_gp = one,diff_c_gp_2 = one, diff_c_gt = one,diff_c_gt_2 = one
    real(rkind) :: minVF = 0.2_rkind, thick = 0.01, p_ten = one
    logical     :: sharp = .FALSE.
    real(rkind) :: p1,p2,rho1,rho2,u1,u2,g11_1,g11_2,grho1,grho2,a1,a2
    real(rkind) :: rho1_2,rho2_2,u1_2,u2_2,g11_1_2,g11_2_2,grho1_2,grho2_2,a1_2,a2_2
    real(rkind) :: rhoL, rhoR, YsL, YsR, VFL, VFR, vL, vR, uL, uR
    real(rkind) :: yield = 0, yield2 = 0, eta0k = 0.4_rkind
    real(rkind) :: melt_t = one, melt_c = one, melt_t2 = one, melt_c2 = one
    real(rkind) :: kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e, alpha3, alpha4,alpha2, alpha
    real(rkind) :: kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2, v_disturb, epsP = 0, epsRho = 0
    real(rkind) :: v0=zero, v0_2=zero, tau0=1d-14, tau0_2=1d-14, Nrho = 1, U0 = zero, m = 1, p_mu = 1, p_mu2 = 1, epsilonk = 0,    run_id
    integer     :: kos_sh,kos_sh2,pointy, pointx, phase_seed
    logical     :: explPlast = .FALSE., explPlast2 = .FALSE.
    logical     :: plastic = .FALSE., plastic2 = .FALSE.
    real(rkind) :: Ly = 1.0, Lx =4*pi,interface_init = 10d-3, kwave = 4.0_rkind, ksize = 10d0, etasize = 0.5d0,delta_d=0.0125D0,delta = 0.0125D0, delta_rho = 0.0125D0 , Lz=2*pi
    real(rkind) :: U_ref, Rho_ref, P_ref, delta_ref
    character(len=1024) :: base_dir, folder_path
    character(len=30) :: temp_alpha_str, temp_beta_str
    integer, parameter :: MAX_MODES = 20
    real(rkind), dimension(MAX_MODES) :: alpha_modes=0.0, beta_modes=0.0, phase_modes=0.0
    integer :: num_modes=0,bnum_modes=0
    real(rkind) :: alpha_dim, beta_dim
    type(filters) :: mygfil

        !TODO: delete all kos stuff and clean up in general

contains

SUBROUTINE fnumden(pf,fparams,iparams,num,den)

  use constants, only: half,third,twothird
  IMPLICIT NONE
  REAL(rkind), INTENT(IN) :: pf
  REAL(rkind), INTENT(IN), DIMENSION(:) :: fparams
  INTEGER, INTENT(IN), DIMENSION(:) :: iparams
  REAL(rkind), INTENT(OUT) :: num,den

  INTEGER :: i, im
  REAL(rkind) :: fac, rho1, u1, p1, p2, rho0, gam, pinf, mus, gm1, gp1, g11, frho1, grho1, Arho1, Brho1, frho2, grho2, fprho2, gprho2

  real(rkind), parameter :: eleventhird = real(11.D0/3.D0,rkind), seventhird = real(7.D0/3.D0,rkind), sixth = real(one/6.D0,rkind), &
                            sevensixth = real(7.D0/6.D0,rkind), eightthird = real(8.D0/3.D0,rkind), onetwone = real(121.D0, rkind), &
                            thirstysix = real(36.D0,rkind), fortynine = real(49.D0,rkind), eighteen = real(18.D0,rkind),            &
                            eighteenth = one/eighteen, ninth = real(one/9.0D0,rkind), threefourth = real(3.D0/4.D0,rkind),          &
                            eleventwelfth = real(11.D0/12.D0,rkind), fourthird = real(4.D0/3.D0,rkind)

  ! if (iparams(1)==PRESSRELAX) then
  !   num = -one; den = zero;
  !   i = iparams(2)
  !   !do im = 1, NUMMAT
  !   !  fac = vfm(i,im)/MAT_GAM(im)/(MAT_PINF(im)+pf)
  !   !  num = num + fac*(psph(i,im)+MAT_GAM(im)*(MAT_PINF(im)+pf)-pf)
  !   !  den = den - fac*(psph(i,im)+MAT_PINF(im))/(pf+MAT_PINF(im))
  !   !enddo
  ! elseif(iparams(1)==SOLIDSTATSHOCK) then
    rho1 = fparams(1); u1 = fparams(2); p1 = fparams(3)
    rho0 = fparams(4); gam = fparams(5); pinf = fparams(6); mus = fparams(7);
    p2 = fparams(8)

    gm1 = gam-one; gp1 = gam+one

    g11 = rho1/rho0    ! g11_1
    grho1 = g11**eleventhird - g11**(-third) - g11**seventhird + g11**third
    frho1 = (eleventwelfth*g11**eleventhird - sixth*g11**(-third) - sevensixth*g11**seventhird -third*g11**third + threefourth*g11)

    Arho1 = (half*u1**two + gam/gm1*(p1+pinf) + mus*frho1)/rho1
    Brho1 = gam/gm1*(p1+pinf+rho1*u1**two+twothird*mus*grho1)

    g11 = pf/rho0     ! g11_2
    grho2 = g11**eleventhird - g11**(-third) - g11**seventhird + g11**third
    frho2 = (eleventwelfth*g11**eleventhird - sixth*g11**(-third) - sevensixth*g11**seventhird -third*g11**third + threefourth*g11)
    gprho2 = one/rho0*(eleventhird*g11**eightthird + third*g11**(-fourthird) &
                     - seventhird*g11**fourthird + third*g11**(-twothird))
    fprho2 = one/rho0*(onetwone/thirstysix*g11**eightthird - fortynine/eighteen*g11**fourthird &
                     + eighteenth*g11**(-fourthird) - ninth*g11**(-twothird) + threefourth)

    !! based on uL
    !num = Arho1 - (Brho1 - half*gp1/gm1*(rho1*u1)**two/pf + mus*frho2 - twothird*gam/gm1*mus*grho2)/pf
    !den = - one/pf*(half*gp1/gm1*(rho1*u1/pf)**two - twothird*mus*gam/gm1*gprho2 + mus*fprho2) &
    !      + one/pf**two*(-half*gp1/gm1*(rho1*u1)**two/pf - twothird*mus*gam/gm1*grho2 + mus*frho2 + Brho1)

    ! based on pR
    num = (gam/gm1*(p1+pinf) + mus*frho1)/rho1 - half*(one/rho1+one/pf)*(p1-p2+twothird*mus*(grho1-grho2)) - mus*frho2/pf - gam/gm1/pf*(p2+pinf)
    den = half/pf**two*(p1-p2+twothird*mus*(grho1-grho2)+two*gam/gm1*(p2+pinf)) + mus*(third*(one/rho1+one/pf)*gprho2 + frho2/pf**two - fprho2/pf)

  ! endif

END SUBROUTINE fnumden

SUBROUTINE rootfind_nr_1d(pf,fparams,iparams)

  IMPLICIT NONE
  REAL(rkind), INTENT(INOUT) :: pf
  REAL(rkind), INTENT(IN), DIMENSION(:) :: fparams
  INTEGER, INTENT(IN), DIMENSION(:) :: iparams

  INTEGER :: ii, itmax = 1000
  REAL(rkind) :: tol = 1.0d-8
  REAL(rkind) :: dpf, num, den, den_conv

  !pfinitguess = pf
  do ii = 1, itmax
    call fnumden(pf,fparams,iparams,num,den)
    if(dabs(den)>1.0d-12) then
      dpf = num/den
    else
      write(*,*) 'den very small, please check.', num, num/den
      stop
    endif
    pf = pf - dpf
    ! check for convergence
    if(dabs(pf)>1.0d-12) then
      den_conv = dabs(pf)
    else
      den_conv = one
    endif
    if(dabs(dpf)/den_conv<1.0d-8) exit
  enddo
  if(ii==itmax+1) then
    write(*,*) 'Newtons method for pf did not converge. Check details.', iparams(1)
  endif

END SUBROUTINE rootfind_nr_1d

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: one, half, pi
    use decomp_2d,        only: decomp_info
    use exits,            only: warning

    use ShearLayer4Mode_data

    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,1)x[0,1)x[0,1) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nx,rkind)
        dy = Ly/real(ny-1,rkind)
        dz = Lz/real(nz,rkind)

        if(abs(dx-dy)>1.0d-13) then
          call warning("dx not equal to dy")
        endif

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1   + i - 1, rkind ) * dx - 2.0*pi
                    y(i,j,k) = real( iy1 - 1  + j - 1, rkind ) * dy  - 0.5d0
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz   - pi
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,der,derStagg,interpMid,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz,periodicx,periodicy,periodicz, x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
    use SolidGrid,        only: u_index,v_index,w_index,rho_index, uref_index,p_index
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: grady,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_x,gradFV_y, gradFV_z
    use DerivativesMod,   only: derivatives
    use DerivativesStaggeredMod, only: derivativesStagg
    use InterpolatorsMod,        only: interpolators
    use reductions,       only: P_SUM, P_MEAN, P_MAXVAL, P_MINVAL
    use ShearLayer4Mode_data
    use mpi

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(derivatives),               intent(in)    :: der
    type(derivativesStagg),          intent(in)    :: derStagg
    type(interpolators),             intent(in)    :: interpMid
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(inout) :: tstop, dt, tviz
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc

    logical :: periodicx,periodicy,periodicz
    integer :: ioUnit,i,j,k,q,nx,ny,nz,ierr,l
    integer, dimension(2) :: iparams
    real(rkind) :: a0, a0_2,dx1
    logical :: adjustRgas = .TRUE.   ! If true, Rgas is used, Rgas2 adjusted to ensure p-T equilibrium
    logical :: adjustPamb = .FALSE.   ! If true, p_amb is adjusted to ensure p-T equilibrium


    real(rkind), allocatable :: theta_fix(:,:), theta_rand(:,:)
    real(rkind) :: mag, magmax
    integer :: j0

    real(rkind) :: Lr, STRETCH_RATIO = 6.0d0
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: eta,tmp,y_stretched
    ! --- Variables for Eigenfunction Initialization ---

    ! Arrays to hold the read-in, pre-interpolated eigenfunctions
    ! Size is now 'pointy', which is the final grid size ny
    real(rkind), dimension(:,:,:), allocatable :: &
        phi_r, phi_i, u_r, u_i, v_r, v_i, w_r, w_i, &
        p_r, p_i, m1_r, m1_i, m2_r, m2_i

    real(rkind), dimension(:), allocatable :: u_base
    real(rkind) :: arg, u_perturb, v_perturb, w_perturb, p_perturb, phi_perturb, m1_perturb, m2_perturb,current_phase,kmag
    real(rkind) :: y_min = -3.5d0, y_max = 6d0, alphar, betar, gammar, y_lo_anchor, y_hi_anchor
    real(rkind), allocatable :: amp_modes(:,:), phase_modes_final(:,:)
    real(rkind) :: umax_loc, umax_glob, target_amp, epsilonk_sol
    real(rkind), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:), pp_fld(:,:,:)


    namelist /PROBINPUT/ p_infty, Rgas, gamma, mu, rho_0, p_amb, thick, minVF, &
                          p_infty_2, Rgas_2, gamma_2, mu_2, rho_0_2, &
                          interface_init, delta, pointy, epsilonk, v0, v0_2, &
                          base_dir, num_modes, alpha_modes, beta_modes, phase_modes, &
                          U_ref, Rho_ref, P_ref, delta_ref,p_amb,Nrho,bnum_modes,run_id

     ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    ! =========================================================================
    ! 1. READ INPUTS
    ! =========================================================================
    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w =>fields(:,:,:,w_index), uref =>fields(:,:,:,uref_index),rho => fields(:,:,:,rho_index),p=>fields(:,:,:,p_index), x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )


        nx = size(mesh,1); ny = size(mesh,2); nz = size(mesh,3)
        !Ensure temperature equilibrium at start
        if(adjustRgas) Rgas_2 = Rgas * (p_amb+p_infty_2)/(p_amb+p_infty)*rho_0/rho_0_2

        ! speed of sound
        a0   = sqrt((gamma*(p_amb+p_infty) )/rho_0)
        a0_2   = sqrt((gamma_2*(p_amb+p_infty_2))/rho_0_2)

        print *, " a0 " , a0
        print *, " a0_2 ", a0_2
        print *, " gamma ", gamma
        print *, " gamma_2 " , gamma_2
        call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield,1.0D-10,eta_det_ge,eta_det_gp,eta_det_gt,diff_c_ge,diff_c_gp,diff_c_gt,melt_t,melt_c,kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e,kos_sh,nx,ny,nz)) !mca: see Sep1SolidEOS.F90 "init"
        call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,1.0D-10,eta_det_ge_2,eta_det_gp_2,eta_det_gt_2,diff_c_ge_2,diff_c_gp_2,diff_c_gt_2,melt_t2,melt_c2,kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2,kos_sh2,nx,ny,nz))

        ! set logicals for plasticity
        mix%material(1)%plast = plastic ; mix%material(1)%explPlast = explPlast
        mix%material(2)%plast = plastic2; mix%material(2)%explPlast = explPlast2

        delta_rho = Nrho * 0.048d0 * 0.275d0
        !delta_rho = Nrho * dx * 0.275d0 !converts from Nrho to approximate thickness of erf profile
!        eta = atanh(2.0d0*y /(1d0 + 1d0/STRETCH_RATIO))
!        Lr    = 14.0d0/(eta(1,ny,1) - eta(1,1,1))
!        eta = Lr*eta

         y_stretched = atanh(2.0d0 * y / (1.0d0 + 1.0d0/STRETCH_RATIO))
         eta = y_min + (y_stretched - minval(y_stretched)) / &
               (maxval(y_stretched) - minval(y_stretched)) * (y_max - y_min)


!         y_lo_anchor = atanh( 2.0d0 * ( 1.0d0/real(ny-1,rkind) - 0.5d0 ) / (1.0d0 + 1.0d0/STRETCH_RATIO) )
!         y_hi_anchor = atanh( 2.0d0 * ( real(ny-2,rkind)/real(ny-1,rkind) - 0.5d0 ) / (1.0d0 + 1.0d0/STRETCH_RATIO) )

!         y_stretched = atanh(2.0d0 * y / (1.0d0 + 1.0d0/STRETCH_RATIO))
!         eta = y_min + (y_stretched - y_lo_anchor) / (y_hi_anchor - y_lo_anchor) * (y_max - y_min)
!        if (ny /= pointy) then
!            call GracefulExit("Grid size mismatch in initfields.", 1)
!        end if

        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=PROBINPUT)
        close(ioUnit)
       
        ! Set default reference values if not provided in input file
        if (U_ref <= 0.0) U_ref = abs(v0_2) ! Assuming v0_2 is the gas free-stream velocity
        if (Rho_ref <= 0.0) Rho_ref = rho_0_2 ! Assuming material 2 is the gas
        if (P_ref <= 0.0) P_ref = p_amb
        if (delta_ref <= 0.0) delta_ref = delta ! Assuming delta is the reference length scale

        if (nrank == 0) then
            print *, '--- LST Initialization Parameters ---'
            write(*,'(A, I3)') 'Number of modes to apply: ', num_modes
            do q = 1, num_modes
              do l = 1, bnum_modes
                write(*, '(A,I3,A,F8.3,A,F8.3,A,F8.3)') 'Mode ', q, ': alpha=', alpha_modes(q), &
                                                       ', beta=', beta_modes(l)
                 enddo
            end do
            print *, '-------------------------------------'
        end if

    ! =========================================================================
    ! 2. SET BASE STATE
    ! =========================================================================
        ! ... (This part is unchanged, set base u, v, w, VF, rho, Ys etc.) ...
        !where(eta .ge. 0)
        !   u = v0_2*erf(eta/delta)
        !elsewhere(eta .lt. 0 )
        !   u = v0*erf(eta/delta)
        !endwhere
        !uref = u
        v = 0.0
        w = 0.0
        mix%material(1)%p = p_amb
        tmp = (half ) * ( one - erf( (eta)/(delta_rho) ) )
        !set mixture Volume fraction
        mix%material(1)%VF = minVF + (one-two*minVF)*tmp
        mix%material(2)%VF = 1.0_rkind - mix%material(1)%VF
        rho = rho_0*mix%material(1)%VF + rho_0_2*mix%material(2)%VF
        mix%material(1)%Ys = mix%material(1)%VF * rho_0 / rho
        mix%material(2)%Ys = 1.0_rkind - mix%material(1)%Ys
        mix%material(2)%p  = mix%material(1)%p

        ! =========================================================================
    ! 3. COMPUTE RANDOM PHASES (seeded by run_id, generated on rank 0, broadcast)
    ! =========================================================================

            ! Allocate arrays to hold the read-in data for all modes
    allocate( phi_r(pointy, MAX_MODES,MAX_MODES), phi_i(pointy, MAX_MODES,MAX_MODES), u_r(pointy, MAX_MODES,MAX_MODES),u_i(pointy, MAX_MODES,MAX_MODES), &
    v_r(pointy, MAX_MODES,MAX_MODES), v_i(pointy, MAX_MODES,MAX_MODES), w_r(pointy, MAX_MODES,MAX_MODES), w_i(pointy,MAX_MODES,MAX_MODES), &
    p_r(pointy, MAX_MODES,MAX_MODES), p_i(pointy, MAX_MODES,MAX_MODES), m1_r(pointy, MAX_MODES,MAX_MODES),m1_i(pointy, MAX_MODES,MAX_MODES), &
    m2_r(pointy, MAX_MODES,MAX_MODES), m2_i(pointy, MAX_MODES,MAX_MODES))
    allocate(  u_base(pointy) )
    
    allocate(theta_fix(num_modes, bnum_modes))
    allocate(theta_rand(num_modes, bnum_modes))
    theta_fix  = 0.0_rkind
    theta_rand = 0.0_rkind

    ! --- Rank 0 reads all data from files ---
    if (nrank == 0) then
        do q = 1, num_modes
            do l = 1, bnum_modes

                write(temp_alpha_str, '(F20.2)') alpha_modes(q)
                write(temp_beta_str,  '(F20.2)') beta_modes(l)
                folder_path = 'alpha_' // trim(adjustl(temp_alpha_str)) // &
                              '_beta_' // trim(adjustl(temp_beta_str)) // '/'
                print *, 'Reading from: ', trim(folder_path)

                open(unit=22, file=trim(folder_path)//'VF_R.txt',  status='old'); read(22,*) phi_r(:,q,l);  close(22)
                open(unit=23, file=trim(folder_path)//'VF_I.txt',  status='old'); read(23,*) phi_i(:,q,l);  close(23)
                open(unit=24, file=trim(folder_path)//'u_R.txt',   status='old'); read(24,*) u_r(:,q,l);    close(24)
                open(unit=25, file=trim(folder_path)//'u_I.txt',   status='old'); read(25,*) u_i(:,q,l);    close(25)
                open(unit=26, file=trim(folder_path)//'v_R.txt',   status='old'); read(26,*) v_r(:,q,l);    close(26)
                open(unit=27, file=trim(folder_path)//'v_I.txt',   status='old'); read(27,*) v_i(:,q,l);    close(27)
                open(unit=28, file=trim(folder_path)//'w_R.txt',   status='old'); read(28,*) w_r(:,q,l);    close(28)
                open(unit=29, file=trim(folder_path)//'w_I.txt',   status='old'); read(29,*) w_i(:,q,l);    close(29)
                open(unit=30, file=trim(folder_path)//'p_R.txt',   status='old'); read(30,*) p_r(:,q,l);    close(30)
                open(unit=31, file=trim(folder_path)//'p_I.txt',   status='old'); read(31,*) p_i(:,q,l);    close(31)
                open(unit=32, file=trim(folder_path)//'m1_R.txt',  status='old'); read(32,*) m1_r(:,q,l);   close(32)
                open(unit=33, file=trim(folder_path)//'m1_I.txt',  status='old'); read(33,*) m1_i(:,q,l);   close(33)
                open(unit=34, file=trim(folder_path)//'m2_R.txt',  status='old'); read(34,*) m2_r(:,q,l);   close(34)
                open(unit=35, file=trim(folder_path)//'m2_I.txt',  status='old'); read(35,*) m2_i(:,q,l);   close(35)
                open(unit=36, file=trim(folder_path)//'Ubase.txt', status='old'); read(36,*) u_base;         close(36)
            end do
        end do

        ! --- Generate reproducible random phases seeded by run_id (rank 0 only) ---
        ! Matches C++: std::mt19937_64 gen(run_id + 0x9E3779B97F4A7C15)
        ! We use a simple LCG seeded the same way for portability.
        ! Seed = run_id XOR-mixed constant (matching the C++ seed)
        block
            integer(kind=8) :: seed64, a_lcg, c_lcg, state
            real(rkind) :: raw
            integer :: qq, ll
            ! LCG parameters (same quality as a quick hash; reproducible)
            a_lcg = 6364136223846793005_8   ! Knuth MMIX multiplier
            c_lcg = 1442695040888963407_8   ! Knuth MMIX increment
            ! Seed: match C++ -> run_id + 0x9E3779B97F4A7C15
            seed64 = int(run_id, 8) + int(z'9E3779B97F4A7C15', 8)
            state  = seed64
            do qq = 1, num_modes
                do ll = 1, bnum_modes
                    ! Advance LCG
                    state = a_lcg * state + c_lcg
                    ! Map to [0, 2*pi): take upper 32 bits, normalize
                    raw = real(iand(int(ishft(state, -32), 8), int(z'FFFFFFFF', 8)), rkind) &
                          / 4294967296.0_rkind
                    theta_rand(qq,ll) = two * pi * raw
                end do
            end do


        end block

    end if  ! nrank == 0

    ! --- Broadcast all eigenfunction data ---
    call MPI_Bcast(phi_r,  pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(phi_i,  pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(u_r,    pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(u_i,    pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(v_r,    pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(v_i,    pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(w_r,    pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(w_i,    pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p_r,    pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p_i,    pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(m1_r,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(m1_i,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(m2_r,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(m2_i,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(u_base, pointy,                     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(theta_rand, num_modes*bnum_modes,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! =========================================================================
    ! 4. COMPUTE theta_fix (phase-lock to eigenfunction peak) + combined phase
    ! =========================================================================
    ! Reuse theta_fix to store the FINAL phase = theta_fix + theta_rand
    ! Store per-mode amplitude (k^-1.5) and final phase in 1D arrays
    ! =========================================================================
    ! 4. COMPUTE theta_fix + von Karman amplitude + combined phase
    ! =========================================================================
    ! von Karman p=1 spectrum: E(kappa) = r/(1+r^2)^{4/3}, r = kappa/kappa_0
    !   Low-k  : E ~ kappa^{+1}  (gentle rise, low-k modes alive for symmetry break)
    !   Peak   : at kappa = kappa_target = 6  (most unstable mode)
    !   High-k : E ~ kappa^{-5/3} (Kolmogorov inertial range)
    !
    ! kappa_0 is shifted from kappa_target because the analytic peak of
    ! r/(1+r^2)^{4/3} occurs at r=sqrt(3/5), not r=1:
    !   kappa_0 = kappa_target / sqrt(3/5) = 6 / 0.7746 = 7.7460
    !
    ! Amplitude: A_ql = sqrt(E(kappa_ql))   [normalised below so sum(A^2)=1]
    ! epsilonk_sol (Step 7) then sets the true 2.5% perturbation level.
    ! =========================================================================

    allocate(amp_modes(num_modes, bnum_modes))
    allocate(phase_modes_final(num_modes, bnum_modes))

    ! =========================================================================
    ! 4. PHYSICALLY MOTIVATED SPECTRAL AMPLITUDE
    ! =========================================================================
    ! Spectrum: Saffman (s=2) -- physically appropriate for shear layers
    !
    !   E(kappa) = (r^2) / (1 + r^2)^{11/6},   r = kappa/kappa_0
    !
    ! where kappa_0 is chosen so the peak of E coincides with the
    ! most unstable wavenumber kappa_delta = 2*pi/delta.
    !
    ! Peak of r^s/(1+r^2)^{(s+5/3)/2} is at r = sqrt(s / (5/3)):
    !   s=2:  r_peak = sqrt(2/(5/3)) = sqrt(6/5) = 1.095
    !   =>    kappa_0 = kappa_peak / r_peak = kappa_delta / sqrt(6/5)
    !
    ! With delta=1:  kappa_delta = 2*pi ~ 6.28
    !                kappa_0     = 6.28 / 1.095 = 5.74
    !
    ! NO FREE PARAMETERS -- everything set by delta.
    !
    ! Low-k  (kappa << kappa_0): E ~ kappa^2   (Saffman, momentum-conserving)
    ! High-k (kappa >> kappa_0): E ~ kappa^{-5/3} (Kolmogorov inertial range)
    !
    ! Reference: Saffman (1967), J. Fluid Mech. 27, 581
    !            Pope (2000), Turbulent Flows, Ch. 6
    ! =========================================================================
    ! =========================================================================
    ! 4. SAFFMAN SPECTRUM WITH SHELL-COUNT CORRECTION
    ! =========================================================================
    ! Physical basis:
    !   E(kappa) = Saffman s=2 spectrum, peak at kappa_delta = 2*sqrt(pi)/delta
    !   (vorticity thickness = sqrt(pi)*delta for erf profile)
    !
    !   Amplitude per mode corrected for shell occupancy:
    !   A_ql = sqrt( E(kappa_ql) / N_shell(kappa_ql) )
    !
    !   This ensures sum_{modes in shell} A^2 = E(kappa)*dkappa
    !   i.e. the discrete spectrum matches the continuous Saffman spectrum
    !   regardless of how many modes fall in each wavenumber band.
    !
    !   Shell width dkappa = 1 (unit shells in wavenumber space)
    ! =========================================================================

    ! =========================================================================
    ! 4. SAFFMAN SPECTRUM WITH LOW-N: FLOOR FOR SYMMETRY BREAKING
    ! =========================================================================
    ! Physical basis:
    !   E(kappa) = max( E_Saffman(kappa), E_floor )
    !
    !   E_Saffman: Saffman s=2 spectrum, peak at kappa_dw = 2*pi/delta_w
    !              where delta_w = sqrt(pi)*delta_sh is vorticity thickness.
    !              E ~ kappa^2 below peak (momentum conserving)
    !              E ~ kappa^{-5/3} above peak (Kolmogorov)
    !
    !   E_floor = f_floor^2 * E_Saffman(kappa_peak)
    !              Minimum energy at any wavenumber.
    !              Ensures low-k modes have at least f_floor fraction
    !              of the peak-mode amplitude at the interface.
    !              Physical interpretation: background broadband excitation.
    !
    !   f_floor = 0.5: low-k modes at least 50% of peak amplitude
    !             This guarantees visual symmetry breaking regardless
    !             of how many high-k modes are present.
    !
    !   Amplitude: A_ql = sqrt( E_seeded(kappa_ql) / N_shell )
    !   Normalisation: sum(A^2 * N_shell) = 1  (unit total spectral energy)
    !   epsilonk_sol sets the true 2.5% perturbation level.
    ! =========================================================================
    block
        real(rkind), parameter :: s        = 2.0_rkind
        real(rkind), parameter :: delta_sh = 1.0_rkind
        real(rkind), parameter :: delta_w  = sqrt(pi) * delta_sh
        real(rkind), parameter :: kappa_dw = two * pi / delta_w
        real(rkind), parameter :: r_peak   = sqrt(s * 3.0_rkind/5.0_rkind)
        real(rkind), parameter :: kappa_0  = kappa_dw / r_peak
        real(rkind), parameter :: expn     = (s + 5.0_rkind/3.0_rkind)/2.0_rkind

        ! f_floor: minimum amplitude as fraction of peak-mode amplitude.
        ! Floor is applied ONLY for kappa < kappa_dw (below spectral peak).
        ! Above kappa_dw: pure Saffman -5/3 rolloff, no floor.
        ! This prevents artificial boosting of high-k modes that cause ringing
        ! at high density ratio where eigenfunctions are narrow at interface.
        real(rkind), parameter :: f_floor  = 0.7_rkind

        real(rkind), parameter :: dkappa  = 1.0_rkind
        integer,     parameter :: nshells = 25
        integer   :: shell_count(nshells)
        integer   :: qq, ll, jj, ishell
        real(rkind) :: r_vk, Ek_vk, E_peak, E_floor_val
        real(rkind) :: amp_sumsq, amp_norm
        real(rkind) :: A_mean, A_min, A_max, A_sum
        integer :: N_in_shell

        E_peak      = r_peak**s / (one + r_peak**2)**expn
        E_floor_val = f_floor**2 * E_peak

        if (nrank == 0) then
            write(*,'(A)')       ''
            write(*,'(A)')       '================================================'
            write(*,'(A)')       ' Saffman s=2 + ONE-SIDED symmetry-breaking floor'
            write(*,'(A)')       '================================================'
            write(*,'(A,F8.4)') '  delta_sh   (shear thickness)  = ', delta_sh
            write(*,'(A,F8.4)') '  delta_w    (vort  thickness)  = ', delta_w
            write(*,'(A,F8.4)') '  kappa_dw   (peak wavenumber)  = ', kappa_dw
            write(*,'(A,F8.4)') '  kappa_0    (formula param)    = ', kappa_0
            write(*,'(A,F8.4)') '  E_peak     (spectrum max)     = ', E_peak
            write(*,'(A,F8.4)') '  f_floor    (low-k fraction)   = ', f_floor
            write(*,'(A,F8.4)') '  E_floor    (minimum E)        = ', E_floor_val
            write(*,'(A)')       '  Spectrum shape:'
            write(*,'(A)')       '    kappa <  kappa_dw: E = max(E_Saffman, E_floor)'
            write(*,'(A)')       '    kappa >= kappa_dw: E = E_Saffman  (pure -5/3)'
            write(*,'(A)')       '    -> high-k modes NOT boosted -> no ringing'
            write(*,'(A)')       '================================================'
        end if

        shell_count = 0
        do qq = 1, num_modes
            do ll = 1, bnum_modes
                kmag   = sqrt(alpha_modes(qq)**2 + beta_modes(ll)**2)
                ishell = max(1, min(nshells, int(kmag/dkappa) + 1))
                shell_count(ishell) = shell_count(ishell) + 1
            end do
        end do

        amp_sumsq = 0.0_rkind

        do qq = 1, num_modes
            do ll = 1, bnum_modes

                magmax = -1.0_rkind
                j0 = 1
                do jj = 1, pointy
                    mag = sqrt(u_r(jj,qq,ll)**2 + u_i(jj,qq,ll)**2)
                    if (mag > magmax) then
                        magmax = mag
                        j0     = jj
                    end if
                end do
                if (magmax < 1.0e-14_rkind) then
                    theta_fix(qq,ll) = 0.0_rkind
                else
                    theta_fix(qq,ll) = atan2(u_i(j0,qq,ll), u_r(j0,qq,ll))
                end if
                phase_modes_final(qq,ll) = theta_fix(qq,ll) + theta_rand(qq,ll)


                kmag   = sqrt(alpha_modes(qq)**2 + beta_modes(ll)**2)
                ishell = max(1, min(nshells, int(kmag/dkappa) + 1))

                if (kmag > 1.0e-12_rkind .and. shell_count(ishell) > 0) then

                    r_vk  = kmag / kappa_0
                    Ek_vk = r_vk**s / (one + r_vk**2)**expn   ! pure Saffman

                    ! OLD: Ek_vk = max(Ek_vk, E_floor_val)
                    !      (floor applied to ALL kappa -> boosts high-k -> ringing)
                    !
                    ! NEW: floor only below spectral peak kappa_dw
                    !      Above kappa_dw: pure Saffman, no artificial boost
                    if (kmag < kappa_dw) then
                        Ek_vk = max(Ek_vk, E_floor_val)
                    end if

                    amp_modes(qq,ll) = sqrt(Ek_vk / real(shell_count(ishell),rkind))
                else
                    amp_modes(qq,ll) = 0.0_rkind
                end if

                amp_sumsq = amp_sumsq + amp_modes(qq,ll)**2

            end do
        end do

        if (amp_sumsq > 1.0e-30_rkind) then
            amp_norm  = one / sqrt(amp_sumsq)
            amp_modes = amp_modes * amp_norm
        end if

        if (nrank == 0) then
            write(*,'(A)') ''
            write(*,'(A)') '  Per-mode amplitudes:'
            write(*,'(A)') '   q   l  alpha   beta   kappa  side      Ek_raw' // &
                           '   floored?  amp_norm'
            do qq = 1, num_modes
                do ll = 1, bnum_modes
                    kmag  = sqrt(alpha_modes(qq)**2 + beta_modes(ll)**2)
                    r_vk  = kmag / kappa_0
                    Ek_vk = r_vk**s / (one + r_vk**2)**expn   ! raw Saffman
                    write(*,'(2I4, 4F8.3, A8, L9, F10.5)') qq, ll,    &
                        alpha_modes(qq), beta_modes(ll), kmag, Ek_vk,  &
                        merge(' low-k  ', ' high-k ', kmag < kappa_dw),&
                        (kmag < kappa_dw .and. Ek_vk < E_floor_val),   &
                        amp_modes(qq,ll)
                end do
            end do

            write(*,'(A)') ''
            write(*,'(A)') '  Shell spectrum:'
            write(*,'(A)') '  shell  kappa_c   N    A_mean    A_min' // &
                           '    A_max   floored?'
            do ishell = 1, nshells
                if (shell_count(ishell) == 0) cycle
                A_sum      = 0.0_rkind
                A_min      = huge(one)
                A_max      = -huge(one)
                N_in_shell = 0
                do qq = 1, num_modes
                    do ll = 1, bnum_modes
                        kmag = sqrt(alpha_modes(qq)**2 + beta_modes(ll)**2)
                        if (max(1,min(nshells,int(kmag/dkappa)+1)) == ishell) then
                            A_sum      = A_sum + amp_modes(qq,ll)
                            A_min      = min(A_min, amp_modes(qq,ll))
                            A_max      = max(A_max, amp_modes(qq,ll))
                            N_in_shell = N_in_shell + 1
                        end if
                    end do
                end do
                A_mean = A_sum / max(N_in_shell, 1)
                ! shell centre kappa
                kmag = (ishell - 0.5_rkind) * dkappa
                write(*,'(I6, F8.3, I5, 3F9.5, A)') ishell,           &
                    (real(ishell,rkind) - 0.5_rkind)*dkappa,           &
                    shell_count(ishell),                                &
                    A_mean, A_min, A_max,                               &
                    merge('  FLOOR', '       ', kmag < kappa_dw .and.  &
                          A_min*A_min*real(shell_count(ishell),rkind)   &
                          < E_floor_val*1.01_rkind)
            end do

            write(*,'(A)') ''
            write(*,'(A,F7.3,A)') '  Modes below kappa_dw=', kappa_dw, &
                                   ':  floor ACTIVE  (low-k symmetry breaking)'
            write(*,'(A,F7.3,A)') '  Modes above kappa_dw=', kappa_dw, &
                                   ':  pure Saffman  (no boost, no ringing)'

        end if

        if (nrank == 0) then
           write(*,'(A)') '  q   l  theta_fix   theta_rand   phase_final   j0    magmax'
           do qq = 1, num_modes
               do ll = 1, bnum_modes
                   write(*,'(2I4, 3F12.6, I6, ES12.4)') qq, ll, &
                       theta_fix(qq,ll), theta_rand(qq,ll), &
                       phase_modes_final(qq,ll), j0, magmax
               end do
           end do
       end if
     end block
!    block
!        real(rkind), parameter :: kappa_target = 1.5_rkind
!        real(rkind), parameter :: kappa_0_vk   = kappa_target / sqrt(3.0_rkind/5.0_rkind)
!        ! = 7.74596...  puts E peak exactly at kappa=6
!        real(rkind) :: r_vk, Ek_vk, amp_sumsq, amp_norm
!        integer :: qq, ll, jj
!
!        amp_sumsq = 0.0_rkind
!
!        do qq = 1, num_modes
!            do ll = 1, bnum_modes
!
!                magmax = -1.0_rkind
!                j0 = 1
!                do jj = 1, pointy
!                    mag = sqrt(u_r(jj,qq,ll)**2 + u_i(jj,qq,ll)**2)
!                    if (mag > magmax) then
!                        magmax = mag
!                        j0 = jj
!                    end if
!                end do
!                if (magmax < 1.0e-14_rkind) then
!                    theta_fix(qq,ll) = 0.0_rkind
!                else
!                    theta_fix(qq,ll) = atan2(u_i(j0,qq,ll), u_r(j0,qq,ll))
!                end if
!
!                phase_modes_final(qq,ll) = theta_fix(qq,ll) + theta_rand(qq,ll)
!
!                kmag = sqrt(alpha_modes(qq)**2 + beta_modes(ll)**2)
!                if (kmag > 1.0e-12_rkind) then
!                    r_vk           = kmag / kappa_0_vk
!                    Ek_vk          = r_vk / (one + r_vk*r_vk)**(4.0_rkind/3.0_rkind)
!                    amp_modes(qq,ll) = sqrt(Ek_vk)
!                else
!                    amp_modes(qq,ll) = 0.0_rkind
!                end if
!
!                amp_sumsq = amp_sumsq + amp_modes(qq,ll)**2
!
!            end do
!        end do
!
!        if (amp_sumsq > 1.0e-30_rkind) then
!            amp_norm = one / sqrt(amp_sumsq)
!            amp_modes = amp_modes * amp_norm
!        end if
!
!        if (nrank == 0) then
!            write(*,'(A,F8.4,A,F8.4)') &
!                'von Karman: kappa_0 = ', kappa_0_vk, &
!                '  peak at kappa = ', kappa_target
!            write(*,'(A,F10.6)') 'sum(amp^2) after norm = ', sum(amp_modes**2)
!            write(*,'(A)') '  q    l   alpha   beta    kappa     amp_norm'
!            do qq = 1, num_modes
!                do ll = 1, bnum_modes
!                    kmag = sqrt(alpha_modes(qq)**2 + beta_modes(ll)**2)
!                    write(*,'(2I4, 4F9.3)') qq, ll, &
!                        alpha_modes(qq), beta_modes(ll), kmag, amp_modes(qq,ll)
!                end do
!            end do
!        end if
!
!    end block
!
    ! =========================================================================
    ! 5. SET BASE FLOW
    ! =========================================================================
    do j = 1, ny
        uref(:,j,:) = u_base(j) * (v0_2 - v0)
        u(:,j,:)    = u_base(j) * (v0_2 - v0)
    end do

    ! =========================================================================
    ! 6. ACCUMULATE PERTURBATIONS WITH epsilonk = 1  (into u then find max|u'|)
    !    Store raw perturbation sums in temporary arrays, then rescale.
    !    We reuse the existing field 'tmp' for |u'| tracking.
    ! =========================================================================
    ! Declare scratch perturbation accumulators
    allocate(up(nx,ny,nz), vp(nx,ny,nz), wp(nx,ny,nz), pp_fld(nx,ny,nz))
    up = 0.0_rkind
    vp = 0.0_rkind
    wp = 0.0_rkind
    pp_fld = 0.0_rkind

    do q = 1, num_modes
        do l = 1, bnum_modes
            alpha_dim = alpha_modes(q) / delta_ref
            beta_dim  = beta_modes(l)  / delta_ref
            current_phase = phase_modes_final(q,l)

            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        arg = alpha_dim * x(i,j,k) + beta_dim * z(i,j,k) - current_phase
                        up(i,j,k)     = up(i,j,k)     + amp_modes(q,l) * U_ref * &
                                        ( u_r(j,q,l)*cos(arg) - u_i(j,q,l)*sin(arg) )
                        vp(i,j,k)     = vp(i,j,k)     + amp_modes(q,l) * U_ref * &
                                        ( v_r(j,q,l)*cos(arg) - v_i(j,q,l)*sin(arg) )
                        wp(i,j,k)     = wp(i,j,k)     + amp_modes(q,l) * U_ref * &
                                        ( w_r(j,q,l)*cos(arg) - w_i(j,q,l)*sin(arg) )
                        pp_fld(i,j,k) = pp_fld(i,j,k) + amp_modes(q,l) * P_ref * &
                                        ( p_r(j,q,l)*cos(arg) - p_i(j,q,l)*sin(arg) )
                    end do
                end do
            end do
        end do
    end do

    ! =========================================================================
    ! 7. SOLVE FOR epsilonk: max|u'| -> 2.5% of (v0_2 - v0)
    ! =========================================================================

    umax_loc = maxval(abs(up))
    call MPI_Allreduce(umax_loc, umax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    target_amp   = 0.025_rkind * abs(v0_2 - v0)
    epsilonk_sol = 0.0_rkind
    if (umax_glob > 1.0e-30_rkind) epsilonk_sol = target_amp / umax_glob

    if (nrank == 0) then
        write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
            'epsilonk solved = ', epsilonk_sol, &
            '  (umax_glob = ',    umax_glob,    &
            ', target = ',        target_amp, ')'
    end if

    ! =========================================================================
    ! 8. ASSEMBLE FINAL FIELDS: base flow + epsilonk_sol * perturbation
    ! =========================================================================
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                u(i,j,k) = u(i,j,k) + epsilonk_sol * up(i,j,k)
                v(i,j,k) = v(i,j,k) + epsilonk_sol * vp(i,j,k)
                w(i,j,k) = w(i,j,k) + epsilonk_sol * wp(i,j,k)
                mix%material(1)%p(i,j,k) = mix%material(1)%p(i,j,k) + epsilonk_sol * pp_fld(i,j,k)
            end do
        end do
    end do

    deallocate(up, vp, wp, pp_fld)
    deallocate(amp_modes, phase_modes_final)
    deallocate(theta_fix, theta_rand)

    ! =========================================================================
    ! 3. READ AND APPLY EIGENFUNCTION PERTURBATIONS
    ! =========================================================================
        ! Allocate arrays to hold the read-in data for all modes
    !    allocate( phi_r(pointy, MAX_MODES,MAX_MODES), phi_i(pointy, MAX_MODES,MAX_MODES), u_r(pointy, MAX_MODES,MAX_MODES),u_i(pointy, MAX_MODES,MAX_MODES), &
    !              v_r(pointy, MAX_MODES,MAX_MODES), v_i(pointy, MAX_MODES,MAX_MODES), w_r(pointy, MAX_MODES,MAX_MODES), w_i(pointy,MAX_MODES,MAX_MODES), &
    !              p_r(pointy, MAX_MODES,MAX_MODES), p_i(pointy, MAX_MODES,MAX_MODES), m1_r(pointy, MAX_MODES,MAX_MODES),m1_i(pointy, MAX_MODES,MAX_MODES), &
    !              m2_r(pointy, MAX_MODES,MAX_MODES), m2_i(pointy, MAX_MODES,MAX_MODES))
    !    allocate(  u_base(pointy) )
    !    
    !    allocate(theta_fix(MAX_MODES,MAX_MODES))
    !    allocate(theta_rand(MAX_MODES,MAX_MODES))
    !    theta_fix  = 0.0_rkind
    !    theta_rand = 0.0_rkind

    !    ! --- Rank 0 reads all data from files ---
    !    if (nrank == 0) then
    !        do q = 1, num_modes
    !          do l = 1,bnum_modes


    !          ! Write to temporary strings with plenty of space
    !          write(temp_alpha_str, '(F20.2)') alpha_modes(q)
    !          write(temp_beta_str, '(F20.2)') beta_modes(l)

    !          ! Assemble the final string, trimming all padding
    !          folder_path = 'alpha_' // trim(adjustl(temp_alpha_str)) // &
    !          '_beta_' // trim(adjustl(temp_beta_str)) // '/'
    !          ! Construct the directory path for the current mode
    !          !  write(folder_path, '("alpha_", F0.2, "_beta_", F0.2, "/")') &
    !          !                  alpha_modes(q), beta_modes(l)
    !            print *, 'Reading from: ', trim(folder_path)

    !            ! Open, read, and close each file
    !            open(unit=22, file=trim(folder_path)//'VF_R.txt', status='old'); read(22,*) phi_r(:,q,l); close(22)
    !            open(unit=23, file=trim(folder_path)//'VF_I.txt', status='old'); read(23,*) phi_i(:,q,l); close(23)
    !            open(unit=24, file=trim(folder_path)//'u_R.txt', status='old');  read(24,*) u_r(:,q,l);   close(24)
    !            open(unit=25, file=trim(folder_path)//'u_I.txt', status='old');  read(25,*) u_i(:,q,l);   close(25)
    !            open(unit=26, file=trim(folder_path)//'v_R.txt', status='old');  read(26,*) v_r(:,q,l);   close(26)
    !            open(unit=27, file=trim(folder_path)//'v_I.txt', status='old');  read(27,*) v_i(:,q,l);   close(27)
    !            open(unit=28, file=trim(folder_path)//'w_R.txt', status='old');  read(28,*) w_r(:,q,l);   close(28)
    !            open(unit=29, file=trim(folder_path)//'w_I.txt', status='old');  read(29,*) w_i(:,q,l);   close(29)
    !            open(unit=30, file=trim(folder_path)//'p_R.txt', status='old');  read(30,*) p_r(:,q,l);   close(30)
    !            open(unit=31, file=trim(folder_path)//'p_I.txt', status='old');  read(31,*) p_i(:,q,l);   close(31)
    !            open(unit=32, file=trim(folder_path)//'m1_R.txt', status='old'); read(32,*) m1_r(:,q,l);  close(32)
    !            open(unit=33, file=trim(folder_path)//'m1_I.txt', status='old'); read(33,*) m1_i(:,q,l);  close(33)
    !            open(unit=34, file=trim(folder_path)//'m2_R.txt', status='old'); read(34,*) m2_r(:,q,l);  close(34)
    !            open(unit=35, file=trim(folder_path)//'m2_I.txt', status='old'); read(35,*) m2_i(:,q,l);  close(35)
    !            open(unit=36, file=trim(folder_path)//'Ubase.txt', status='old'); read(36,*) u_base;  close(36)
    !        end do
    !        end do
    !    end if

    !    ! --- Broadcast the data from Rank 0 to all other processes ---
    !    call MPI_Bcast(phi_r, pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(phi_i, pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(u_r,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(u_i,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(v_r,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(v_i,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(w_r,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(w_i,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(p_r,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(p_i,   pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(m1_r,  pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(m1_i,  pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(m2_r,  pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !    call MPI_Bcast(m2_i,  pointy*MAX_MODES*MAX_MODES, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    !    call MPI_Bcast(u_base, pointy, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    !    do q = 1, num_modes
    !        do l = 1, bnum_modes
    !      
    !          magmax = -1.0_rkind
    !          j0 = 1
    !          do j = 1, pointy
    !            mag = sqrt(u_r(j,q,l)**2 + u_i(j,q,l)**2)
    !            if (mag > magmax) then
    !              magmax = mag
    !              j0 = j
    !            end if
    !          end do
    !      
    !          if (magmax < 1.0e-14_rkind) then
    !            theta_fix(q,l) = 0.0_rkind   ! fallback: mode has ~zero u everywhere
    !          else
    !            theta_fix(q,l) = atan2(u_i(j0,q,l), u_r(j0,q,l))
    !          end if
    !      
    !        end do
    !    end do


    !    do j = 1, ny

    !       uref(:,j,:) = u_base(j)*(v0_2 - v0)
    !       u(:,j,:)    = u_base(j)*(v0_2 - v0)

    !    enddo

    !    ! --- Loop over modes and add perturbations to the 3D field ---
    !    do q = 1, num_modes
    !      do l = 1, bnum_modes
    !        alpha_dim = alpha_modes(q) / delta_ref
    !        beta_dim  = beta_modes(l) / delta_ref
    !        kmag      = sqrt( alpha_dim*alpha_dim + beta_dim*beta_dim)
    !        alphar = 0.6180339887498949_rkind   ! golden ratio conjugate
    !        betar  = 0.4142135623730950_rkind   ! sqrt(2)-1
    !        gammar = 0.7320508075688772_rkind   ! sqrt(3)-1

    !        current_phase = theta_fix(q,l) + two*pi*fract( alphar*real(run_id,rkind) + &
    !                                          betar *real(q,rkind)      + &
    !                                          gammar*real(l,rkind) )
    !        !current_phase = 2.0_rkind * pi * sin( real(q)*(1.27d0) + real(l)*3.81d0)**2d0
    !        do k = 1, nz
    !            do j = 1, ny
    !                do i = 1, nx
    !                    arg = alpha_dim * x(i,j,k) + beta_dim * z(i,j,k) - current_phase

    !                    ! Construct dimensional perturbations using pre-interpolated eigenfunction values at y-index 'j'
    !                    u_perturb = kmag**(-3.0_rkind / 2.0_rkind) * U_ref * ( u_r(j,q,l)*cos(arg) - u_i(j,q,l)*sin(arg) )
    !                    v_perturb = kmag**(-3.0_rkind / 2.0_rkind) * U_ref * ( v_r(j,q,l)*cos(arg) - v_i(j,q,l)*sin(arg) )
    !                    w_perturb = kmag**(-3.0_rkind / 2.0_rkind) * U_ref * ( w_r(j,q,l)*cos(arg) - w_i(j,q,l)*sin(arg) )
    !                    p_perturb = kmag**(-3.0_rkind / 2.0_rkind) * P_ref * ( p_r(j,q,l)*cos(arg) - p_i(j,q,l)*sin(arg) )
    !    !                phi_perturb =       ( phi_r(j,q,l)*cos(arg) - phi_i(j,q,l)*sin(arg) )
    !    !                m1_perturb = Rho_ref * ( m1_r(j,q,l)*cos(arg) - m1_i(j,q,l)*sin(arg) )
    !    !                m2_perturb = Rho_ref * ( m2_r(j,q,l)*cos(arg) - m2_i(j,q,l)*sin(arg) )

    !                    ! Apply the perturbations, scaled by the single amplitude 'epsilonk'
    !                    u(i,j,k) = u(i,j,k) + epsilonk * u_perturb
    !                    v(i,j,k) = v(i,j,k) + epsilonk * v_perturb
    !                    w(i,j,k) = w(i,j,k) + epsilonk * w_perturb
    !                    mix%material(1)%p(i,j,k) = mix%material(1)%p(i,j,k) + epsilonk * p_perturb
    !     !               mix%material(1)%VF(i,j,k) = mix%material(1)%VF(i,j,k) + epsilonk * phi_perturb
    !     !               rho(i,j,k) = rho(i,j,k) + epsilonk * (m1_perturb + m2_perturb)
    !                end do
    !            end do
    !        end do
    !    end do
    !    end do

        ! --- Finalize mixture properties after all perturbations are added ---
        !mix%material(1)%VF = max(minVF, min(1.0_rkind - minVF, mix%material(1)%VF))
        mix%material(2)%VF = 1.0_rkind - mix%material(1)%VF
        !mix%material(1)%Ys = mix%material(1)%VF * rho_0 / rho
        !mix%material(2)%Ys = max(0.0_rkind, 1.0_rkind - mix%material(1)%Ys) ! Ensure Ys sums to 1 and is non-negative
        mix%material(2)%p  = mix%material(1)%p
        p = mix%material(1)%p
        deallocate(phi_r, phi_i, u_r, u_i, v_r, v_i, w_r, w_i, p_r, p_i, m1_r, m1_i, m2_r, m2_i)
        deallocate(u_base)
    ! =========================================================================
    ! 4. CLEANUP AND BOUNDARY CONDITIONS
    ! =========================================================================
            ! Set initial values of g (inverse deformation gradient)
        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one

        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

        !Stuff for boundary conditions
        rhoL = rho(1,1,1)
        rhoR = rho(1,decomp%ysz(2),1)
        uL = uref(1,1,1)
        uR = uref(1,decomp%ysz(2),1)
        vL = 0 !v(1,1,1)
        vR = 0 !v(1,decomp%ysz(2),1)
        YsL  = mix%material(1)%Ys(1,1,1)
        YsR  = mix%material(1)%Ys(1,decomp%ysz(2),1)
        VFL  = mix%material(1)%VF(1,1,1)
        VFR  = mix%material(1)%VF(1,decomp%ysz(2),1)


    end associate

    contains

  real(rkind) function fract(x)
    real(rkind), intent(in) :: x
    fract = x - floor(x)
  end function fract
  real(rkind) function phase_hash(q,l,seed)
       use kind_parameters, only: rkind
       use constants, only: pi, two
       implicit none
       integer, intent(in) :: q,l,seed
       real(rkind) :: t, u

       t = 12.9898_rkind*real(q,rkind) + 78.233_rkind*real(l,rkind) + 37.719_rkind*real(seed,rkind)
       u = sin(t) * 43758.5453_rkind
       u = u - floor(u)          ! fract(u) inline
       phase_hash = two*pi*u
  end function phase_hash
end subroutine initfields

subroutine get_sponge(decomp,dx,dy,dz,mesh,fields,mix,rhou,rhov,rhow,rhoe,sponge,mask)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
    use SolidGrid,        only: u_index,v_index,w_index,rho_index,e_index, uref_index
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit
    use SolidMixtureMod,  only: solid_mixture
    use ShearLayer4Mode_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    type(solid_mixture),             intent(inout) :: mix
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout):: sponge
    real(rkind), dimension(:,:,:), intent(inout):: mask
    real(rkind), dimension(2), intent(inout) :: rhou, rhov,rhow,rhoe
    integer :: ioUnit,i,iy
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp,dum, eta, eta2, yphys,y_stretched
    real(rkind) :: fac, Lr, STRETCH_RATIO = 6.0d0,int_KE
    integer, dimension(2) :: iparams
    real(rkind) :: a0, a0_2, sigma1, sigma2, y_min = -3.5d0, y_max = 6d0, y_lo_anchor, y_hi_anchor
    integer :: nx,ny,nz,k,ix,j
    integer :: ierr, rank,fh, filesize, chunksize, offset, offset2,totalproc
    integer, allocatable :: data(:), recvbuf(:) 


        associate(u => fields(:,:,:,u_index), v => fields(:,:,:,v_index),w => fields(:,:,:,w_index),uref => fields(:,:,:,uref_index), rho => fields(:,:,:,rho_index), e => fields(:,:,:,e_index), x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )


        
        nx = size(mesh,1); ny = size(mesh,2); nz = size(mesh,3)
!        yphys = atanh(2.0d0*y /(1d0 + 1d0/STRETCH_RATIO))
!        Lr    = 14.0d0/(yphys(1,ny,1) - yphys(1,1,1))
!        yphys = Lr*yphys
 
        y_stretched = atanh(2.0d0 * y / (1.0d0 + 1.0d0/STRETCH_RATIO))
        yphys = y_min + (y_stretched - minval(y_stretched)) / &
               (maxval(y_stretched) - minval(y_stretched)) * (y_max - y_min)

        !y_lo_anchor = atanh( 2.0d0 * ( 1.0d0/real(ny-1,rkind) - 0.5d0 ) / (1.0d0 + 1.0d0/STRETCH_RATIO) )
        !y_hi_anchor = atanh( 2.0d0 * ( real(ny-2,rkind)/real(ny-1,rkind) - 0.5d0 ) / (1.0d0 + 1.0d0/STRETCH_RATIO) )

        !y_stretched = atanh(2.0d0 * y / (1.0d0 + 1.0d0/STRETCH_RATIO))
        !yphys = y_min + (y_stretched - y_lo_anchor) / (y_hi_anchor - y_lo_anchor) * (y_max - y_min)
        sigma1 = -1000d0 ! -2400 ! -80000

        where(yphys .LE. -2.0d0)
           sponge(:,:,:,1) = sigma1*( (yphys + 2.0d0)/0.5d0)**2.0d0 * 2.0d0
           mask  = 1.0
        elsewhere
           sponge(:,:,:,1) = 0d0
           mask = 0.0
        endwhere

        where(yphys .GE. 5.5d0)
           sponge(:,:,:,2) = sigma1*( (yphys- 5.5d0)/0.5d0)**2.0d0 
           mask = 1.0
        elsewhere
           sponge(:,:,:,2) = 0
           mask = 0.0
        endwhere

        rhou(1) = rho(1,1,1)*v0!uref(1,1,1)
        rhou(2) = rho(1,ny,1)*v0_2 !uref(1,ny,1)
        rhov(1) = 0 !-1.060981230880199d-5 !rho(1,1,1)*v(1,1,1)
        rhov(2) = 0 !2.175685479370164d-08
        rhow(1) = 0 !rho(1,1,1)*w(1,1,1)
        rhow(2) = 0 !rho(1,ny,1)*w(1,ny,1)
        rhoe(1) = rho(1,1,1)*(e(1,1,1) + 0.5d0*(v0**2d0)) ! 828.903*(3.4899086 + 0.5*(v0**2)) !1d3*(581967.7419+ 0.5*(v0**2)) !1.d0*(103.176 + 0.5*(v0**2))
        rhoe(2) = rho(1,ny,1)*(e(1,ny,1) + 0.5d0*(v0_2**2d0)) !1d0*(1.785714 + 0.5*(v0_2**2)) !1d0*(250000 + 0.5*(v0_2**2))
        do i = 1,2
          mix%material(i)%VF_ref(1) = mix%material(i)%VF(1,1,1)
          mix%material(i)%VF_ref(2) = mix%material(i)%VF(1,ny,1) 
          mix%material(i)%Ys_ref(1) = mix%material(i)%Ys(1,1,1)*rho(1,1,1)
          mix%material(i)%Ys_ref(2) = mix%material(i)%Ys(1,ny,1)*rho(1,ny,1)
        enddo

    !    print *, "rhou ", rhou
    !    print *, "rhov ", rhov
    !    print *, "rhow ", rhow
    !    print *, "rhoe ", rhoe
    !    print *, "VF ", mix%material(1)%VF_ref
    !    print *, "Ys ", mix%material(1)%Ys_ref
        end associate

end subroutine

subroutine initparam_restart(decomp,der,derStagg,interpMid,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
    use SolidGrid,        only: u_index,v_index,w_index,rho_index, uref_index
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: grady,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_x,gradFV_y, gradFV_z
    use DerivativesMod,   only: derivatives  
    use DerivativesStaggeredMod, only: derivativesStagg
    use InterpolatorsMod,        only: interpolators
    use reductions,       only: P_SUM, P_MEAN, P_MAXVAL, P_MINVAL 
    use ShearLayer4Mode_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(derivatives),               intent(in)    :: der
    type(derivativesStagg),          intent(in)    :: derStagg
    type(interpolators),             intent(in)    :: interpMid
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(inout) :: tstop, dt, tviz
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc

    logical :: periodicx,periodicy,periodicz
    integer :: ioUnit,i,iy
    real(rkind), dimension(8) :: fparams
    real(rkind), dimension(4) :: alphai, phase
    real(rkind) :: fac, Lr, STRETCH_RATIO = 6.0, int_KE
    integer, dimension(2) :: iparams
    real(rkind) :: a0, a0_2,dx1
    logical :: adjustRgas = .TRUE.   ! If true, Rgas is used, Rgas2 adjusted to ensure p-T equilibrium
    logical :: adjustPamb = .FALSE.   ! If true, p_amb is adjusted to ensure p-T equilibrium

    integer :: nx,ny,nz,k,ix,j
    integer :: ierr, rank,fh, filesize, chunksize, offset, offset2,totalproc
    integer, allocatable :: data(:), recvbuf(:)
    character(len=50) :: filename,phiIname,phiRname, DphiIname, DphiRname
    character(len=50) :: rhoRname, rhoIname, pIname,pRname
    !nteger(kind=MPI_OFFSET_KIND) :: disp
    !nteger(kind=MPI_STATUS_SIZE) :: status(MPI_STATUS_SIZE)
    logical :: flag

    ! Initialize MPI
    !all MPI_Init(ierr)

    namelist /PROBINPUT/ p_infty, Rgas, gamma, mu, rho_0, p_amb, thick, minVF, &
                          p_infty_2, Rgas_2, gamma_2, mu_2, rho_0_2, &
                          interface_init, delta, pointy, epsilonk, v0, v0_2, &
                          base_dir, num_modes, alpha_modes, beta_modes, phase_modes, &
                          U_ref, Rho_ref, P_ref, delta_ref,p_amb,Nrho,bnum_modes, run_id

     ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)



    call mygfil%init(                        decomp, &
                     .FALSE.,     .TRUE.,    .TRUE., &
                  "gaussian", "gaussian", "gaussian" )


    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w =>fields(:,:,:,w_index), uref => fields(:,:,:,uref_index),rho => fields(:,:,:,rho_index), x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        if (mix%ns /= 2) then
            call GracefulExit("Number of species must be 2 for this problem.Check the input file.",928)
        end if

        nx = decomp%ysz(1)
        ny = decomp%ysz(2)
        nz = decomp%ysz(3)

        print *, nx
        print *, ny
        print *, nz    
        !Ensure temperature equilibrium at start
        if(adjustRgas) Rgas_2 = Rgas * (p_amb+p_infty_2)/(p_amb+p_infty)*rho_0/rho_0_2
        ! speed of sound
        a0   = sqrt((gamma*(p_amb+p_infty) )/rho_0)
        a0_2   = sqrt((gamma_2*(p_amb+p_infty_2)  )/rho_0_2)

       !Stuff for boundary conditions
        rhoL = rho_0
        rhoR = rho_0_2
        uL = v0
        uR = v0_2
        vL = 0 !v(1,1,1)
        vR = 0 !v(1,decomp%ysz(2),1)
        YsL  = mix%material(1)%Ys(1,1,1)
        YsR  = mix%material(1)%Ys(1,decomp%ysz(2),1)
        VFL  = mix%material(1)%VF(1,1,1)
        VFR  = mix%material(1)%VF(1,decomp%ysz(2),1)

        ! write material properties
!        if (nrank == 0) then
!            print *, '---Material 1---'
!            write(*,'(3(a,e12.5))') 'rho_0 = ', rho_0, ', gam  = ', gamma, ',p_infty = ', p_infty
!            write(*,'(3(a,e12.5))') 'shearMod    = ', mu,    ', Rgas = ', Rgas,', SOS = ', a0
!            write(*,'(3(a,e12.5))') 'tau0  = ', tau0
!            print *, '---Material 2---'
!            write(*,'(3(a,e12.5))') 'rho_0 = ', rho_0_2, ', gam  = ', gamma_2,', p_infty = ', p_infty_2
!            write(*,'(3(a,e12.5))') 'shearMod    = ', mu_2,    ', Rgas = ',Rgas_2, ', SOS = ', a0_2
!            write(*,'(3(a,e12.5))') 'tau0  = ', tau0_2
!            write(*,*) 'p_amb = ', p_amb
!        end if

      !  stiffgas(gamma  ,Rgas  ,p_infty )
      !  call mix%set_material_restart(1,stiffgas(gamma, Rgas, p_infty))

!mca: see Sep1SolidEOS.F90 "init"
       ! sep1solid(rho_0_2,mu_2,yield2,1.0D-10,eta_det_ge_2,eta_det_gp_2,eta_det_gt_2,diff_c_ge_2,diff_c_gp_2,diff_c_gt_2,melt_t2,melt_c2,kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2,kos_sh2,nx,ny,nz)
      !  call mix%set_material_restart(2,stiffgas(gamma_2,Rgas_2,p_infty_2))
        
        call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty),sep1solid(rho_0  ,mu ,yield,1.0D-10,eta_det_ge,eta_det_gp,eta_det_gt,diff_c_ge,diff_c_gp,diff_c_gt,melt_t,melt_c,kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e,kos_sh,nx,ny,nz))
!mca: see Sep1SolidEOS.F90 "init"
        call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,1.0D-10,eta_det_ge_2,eta_det_gp_2,eta_det_gt_2,diff_c_ge_2,diff_c_gp_2,diff_c_gt_2,melt_t2,melt_c2,kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2,kos_sh2,nx,ny,nz))


        ! set logicals for plasticity
        mix%material(1)%plast = plastic ; mix%material(1)%explPlast = explPlast
        mix%material(2)%plast = plastic2; mix%material(2)%explPlast = explPlast2

        ! Set initial values of g (inverse deformation gradient)
        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one

        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

        !Stuff for boundary conditions
        
    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount,pthick,rhothick,uthick,Ysthick,VFthick,Ys_wiggle,VF_wiggle,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,eps,half,one,two,pi,four,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                sxx_index,syy_index,szz_index,sxy_index,sxz_index,syz_index,sos_index, uref_index
    use decomp_2d,        only: decomp_info, nrank
    use DerivativesMod,   only: derivatives
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: curl
    use reductions,       only: P_SUM, P_MEAN, P_MAXVAL, P_MINVAL

    use ShearLayer4Mode_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der   
    real(rkind),                     intent(in) :: dx,dy,dz,tsim,uthick,rhothick,pthick,Ysthick,VFthick,Ys_wiggle,VF_wiggle
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix
    integer, dimension(2),           intent(in) :: x_bc, y_bc, z_bc
    integer                                     :: outputunit=229

    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),3) :: vort
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)  ) :: tmp
    real(rkind), dimension(decomp%ysz(1)) :: Ys1_mean,Ys2_mean
    real(rkind) :: vort_pos, vort_neg, mixwidth, Al_mass, xspike, xbubbl, xspike_proc, xbubbl_proc
    real(rkind) :: YsGrowth,VFGrowth, VFmin_proc, VFmax_proc, VFmin, VFmax
    character(len=clen) :: outputfile, str
    integer :: i, j, k

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 sxx  => fields(:,:,:, sxx_index), syy => fields(:,:,:,syy_index), &
                 szz  => fields(:,:,:, szz_index), sxy => fields(:,:,:,sxy_index), &
                 sxz  => fields(:,:,:, sxz_index), syz => fields(:,:,:,syz_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3),       &
                 sos  => fields(:,:,:,sos_index), uref =>fields(:,:,:,uref_index ))

       write(str,'(I4.4,A,ES7.1E2,A,ES7.1E2)') nrank, "_", minVF, "_", rho_0_2/rho_0
       
       if (mix%use_gTg) then
           str = trim(str)//'_gTg'
       else
           str = trim(str)//'_g'
       end if

       if (decomp%ysz(2) == 1) then
           write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ShearLayer4Mode_"//trim(str)//"_", vizcount, ".dat"

           open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
           write(outputunit,'(4ES27.16E3)') tsim, minVF, thick, rho_0_2/rho_0
           do i=1,decomp%ysz(1)
               write(outputunit,'(23ES27.16E3)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                              mix%material(1)%p (i,1,1), mix%material(2)%p (i,1,1), &
                                              mix%material(1)%Ys(i,1,1), mix%material(2)%Ys(i,1,1), &
                                              mix%material(1)%VF(i,1,1), mix%material(2)%VF(i,1,1), &
                                              mix%material(1)%eh(i,1,1), mix%material(2)%eh(i,1,1), &
                                              mix%material(1)%T (i,1,1), mix%material(2)%T (i,1,1), &
                                              mix%material(1)%g11(i,1,1), mix%material(2)%g11(i,1,1), &
                                              mu(i,1,1), bulk(i,1,1), mix%material(1)%kap(i,1,1), mix%material(2)%kap(i,1,1), &
                                              mix%material(1)%diff(i,1,1), mix%material(2)%diff(i,1,1)
           end do
           close(outputunit)
       end if

       call curl(decomp, der, u, v, w, vort, x_bc, y_bc, z_bc)
       
       tmp = zero
       where (vort(:,:,:,3) .GE. zero)
           tmp = vort(:,:,:,3)
       end where
       vort_pos = P_MEAN(tmp)*six*one
       
       tmp = zero
       where (vort(:,:,:,3) .LE. zero)
           tmp = vort(:,:,:,3)
       end where
       vort_neg = P_MEAN(tmp)*six*one

       Ys1_mean = SUM(mix%material(1)%Ys(:,:,1),2) / real(decomp%ysz(2),rkind)
       Ys2_mean = SUM(mix%material(2)%Ys(:,:,1),2) / real(decomp%ysz(2),rkind)

       Ys1_mean = four * Ys1_mean * Ys2_mean
       mixwidth = P_SUM(Ys1_mean) * dx

       Al_mass = P_MEAN(rho*mix%material(2)%Ys)*six*one

       xspike_proc = 0
       xbubbl_proc = Ly
       VFmin_proc  = Ly 
       VFmax_proc  = 0
       do j = 1,decomp%ysz(2)
           do i = 1,decomp%ysz(1)
               if (mix%material(1)%Ys(i,j,1) .GE. half) xspike_proc = max(xspike_proc,y(i,j,1))
               if (mix%material(1)%Ys(i,j,1) .LE. half) xbubbl_proc = min(xbubbl_proc,y(i,j,1))
               if (mix%material(1)%VF(i,j,1) .GE. half) VFmax_proc  = max(VFmax_proc,y(i,j,1))
               if (mix%material(1)%VF(i,j,1) .LE. half) VFmin_proc  = min(VFmin_proc,y(i,j,1))
           end do
       end do
       xspike = P_MAXVAL(xspike_proc)
       xbubbl = P_MINVAL(xbubbl_proc)
       VFmax  = P_MAXVAL(VFmax_proc)
       VFmin  = P_MINVAL(VFmin_proc)
       YsGrowth = xspike - xbubbl
       VFGrowth  = VFmax - VFmin   
       write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ShearLayer4Mode_statistics.dat"

       if (vizcount == 0) then
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
           write(outputunit,'(7A27)') 'tsim', 'YsAmp', 'VFAmp'
       else
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED', action='WRITE', status='OLD', position='APPEND')
       end if
       
       if(nrank .eq. 0) then
         write(outputunit,'(7ES27.16E3)') tsim, YsGrowth, VFGrowth

       endif
       close(outputunit)
 
       write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ShearLayer_LADstatistics.dat"

       if (vizcount == 0) then
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED',status='REPLACE')
           write(outputunit,'(7A27)') 'tsim', 'Ys_thick', 'VF_thick', 'pthick','uthick', 'Ys_wiggle', 'VF_wiggle', "rhothick"
       else
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED',action='WRITE', status='OLD', position='APPEND')
       end if

       if(nrank .eq. 0) then
         write(outputunit,'(7ES27.16E3)') tsim, Ysthick, VFthick, pthick, uthick,Ys_wiggle, VF_wiggle, rhothick
       endif
       close(outputunit)


    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,uref_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: filter3D

    use ShearLayer4Mode_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc
    
    integer :: nx,ny, i, j
    real(rkind) :: dy, yspng, tspng, yspngR, yspngL, Lr, STRETCH_RATIO = 6.0
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum, dumL, dumR, yphys

    ! Coefficients matching C++ homoNeumannBC/dirichletBC (two_phase_shear_layer.cpp)
    real(rkind), parameter :: nb_ainv = -24.0_rkind/23.0_rkind, nb_bneg = -7.0_rkind/8.0_rkind, &
                               nb_cneg = -1.0_rkind/8.0_rkind,  nb_dneg =  1.0_rkind/24.0_rkind
    real(rkind), parameter :: db_ainv =  16.0_rkind/5.0_rkind,  db_bneg = -15.0_rkind/16.0_rkind, &
                               db_cneg =   5.0_rkind/16.0_rkind, db_dneg = -1.0_rkind/16.0_rkind

    nx = decomp%ysz(1)
    ny = decomp%ysz(2)

    !print *, "ny", ny
    
    mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
    mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
    mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one

    mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
    mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
    mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 uref => fields(:,:,:,uref_index),                                  &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )


        if(decomp%yst(2)==1) then
         if(y_bc(1)==0) then
             ! Matches C++ setBoundaryConditionsAtCollocationPoints, eta_lo:
             !   u,w,p,rho -> homoNeumannBC(f(j+1),f(j+2),f(j+3))
             !   v         -> dirichletBC(0.0, ...)
             !   VF -> phi -> dirichletBC(VFL, ...) ; Ys -> m1/m2 -> dirichletBC(YsL, ...)
             rho( :,1,:) = nb_ainv*( nb_bneg*rho( :,2,:) + nb_cneg*rho( :,3,:) + nb_dneg*rho( :,4,:) )
             u  ( :,1,:) = nb_ainv*( nb_bneg*u  ( :,2,:) + nb_cneg*u  ( :,3,:) + nb_dneg*u  ( :,4,:) )
             v  ( :,1,:) = db_ainv*( db_bneg*v  ( :,2,:) + db_cneg*v  ( :,3,:) + db_dneg*v  ( :,4,:) + zero )
             w  ( :,1,:) = nb_ainv*( nb_bneg*w  ( :,2,:) + nb_cneg*w  ( :,3,:) + nb_dneg*w  ( :,4,:) )

             mix%material(1)%p(:,1,:) = nb_ainv*( nb_bneg*mix%material(1)%p(:,2,:) + nb_cneg*mix%material(1)%p(:,3,:) + nb_dneg*mix%material(1)%p(:,4,:) )
             mix%material(2)%p(:,1,:) = nb_ainv*( nb_bneg*mix%material(2)%p(:,2,:) + nb_cneg*mix%material(2)%p(:,3,:) + nb_dneg*mix%material(2)%p(:,4,:) )
             p(:,1,:)                 = nb_ainv*( nb_bneg*p(:,2,:) + nb_cneg*p(:,3,:) + nb_dneg*p(:,4,:) )

             mix%material(1)%VF ( :,1,:) = db_ainv*( db_bneg*mix%material(1)%VF ( :,2,:) + db_cneg*mix%material(1)%VF ( :,3,:) + db_dneg*mix%material(1)%VF ( :,4,:) + VFL )
             mix%material(2)%VF ( :,1,:) = one - mix%material(1)%VF ( :,1,:)
             mix%material(1)%Ys ( :,1,:) = db_ainv*( db_bneg*mix%material(1)%Ys ( :,2,:) + db_cneg*mix%material(1)%Ys ( :,3,:) + db_dneg*mix%material(1)%Ys ( :,4,:) + YsL )
             mix%material(2)%Ys ( :,1,:) = one - mix%material(1)%Ys ( :,1,:)
         end if
       endif

       if(decomp%yen(2)==decomp%ysz(2)) then
         if(y_bc(2)==0) then
             ! Matches C++ setBoundaryConditionsAtCollocationPoints, eta_hi
             rho( :,ny,:) = nb_ainv*( nb_bneg*rho( :,ny-1,:) + nb_cneg*rho( :,ny-2,:) + nb_dneg*rho( :,ny-3,:) )
             u  ( :,ny,:) = nb_ainv*( nb_bneg*u  ( :,ny-1,:) + nb_cneg*u  ( :,ny-2,:) + nb_dneg*u  ( :,ny-3,:) )
             v  ( :,ny,:) = db_ainv*( db_bneg*v  ( :,ny-1,:) + db_cneg*v  ( :,ny-2,:) + db_dneg*v  ( :,ny-3,:) + zero )
             w  ( :,ny,:) = nb_ainv*( nb_bneg*w  ( :,ny-1,:) + nb_cneg*w  ( :,ny-2,:) + nb_dneg*w  ( :,ny-3,:) )

             mix%material(1)%p(:,ny,:) = nb_ainv*( nb_bneg*mix%material(1)%p(:,ny-1,:) + nb_cneg*mix%material(1)%p(:,ny-2,:) + nb_dneg*mix%material(1)%p(:,ny-3,:) )
             mix%material(2)%p(:,ny,:) = nb_ainv*( nb_bneg*mix%material(2)%p(:,ny-1,:) + nb_cneg*mix%material(2)%p(:,ny-2,:) + nb_dneg*mix%material(2)%p(:,ny-3,:) )
             p(:,ny,:)                 = nb_ainv*( nb_bneg*p(:,ny-1,:) + nb_cneg*p(:,ny-2,:) + nb_dneg*p(:,ny-3,:) )

             mix%material(1)%VF ( :,ny,:) = db_ainv*( db_bneg*mix%material(1)%VF ( :,ny-1,:) + db_cneg*mix%material(1)%VF ( :,ny-2,:) + db_dneg*mix%material(1)%VF ( :,ny-3,:) + VFR )
             mix%material(2)%VF ( :,ny,:) = one - mix%material(1)%VF ( :,ny,:)
             mix%material(1)%Ys ( :,ny,:) = db_ainv*( db_bneg*mix%material(1)%Ys ( :,ny-1,:) + db_cneg*mix%material(1)%Ys ( :,ny-2,:) + db_dneg*mix%material(1)%Ys ( :,ny-3,:) + YsR )
             mix%material(2)%Ys ( :,ny,:) = one - mix%material(1)%Ys ( :,ny,:)
         end if
       endif

     
        
        
  ! apply sponge at left and right boundaries to damp outgoing waves
        yphys = atanh(2.0*y /(1 + 1/STRETCH_RATIO))
        Lr    = 24D0
        yphys = Lr*yphys

        yspngL = -0.85 !250
        yspngR = 0.85
        tspng = 0.05
        dumL = half*(one - tanh( (y-yspngL)/(tspng) ))
        dumR = half*(one + tanh( (y-yspngR)/(tspng) ))
        dum  = dumL+dumR

       ! do i=1,4
       !     tmp = u
       !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
       !     u = u + dum*(tmp - u)

       !     tmp = v
       !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
       !     v = v + dum*(tmp - v)

       !     tmp = w
       !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
       !     w = w + dum*(tmp - w)

       !     tmp = e
       !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
       !     e = e + dum*(tmp - e)

       !     tmp = rho
       !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
       !     rho = rho + dum*(tmp - rho)

       !     tmp = mix%material(1)%p
       !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
       !     mix%material(1)%p = mix%material(1)%p + dum*(tmp - mix%material(1)%p)

       !     tmp = mix%material(2)%p
       !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
       !     mix%material(2)%p = mix%material(2)%p + dum*(tmp - mix%material(2)%p)

            ! TODO: delete tmp = mix%material(1)%pe
            ! TODO: delete call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            ! TODO: delete mix%material(1)%pe = mix%material(1)%pe + dum*(tmp - mix%material(1)%pe)

            ! TODO: delete tmp = mix%material(2)%pe
            ! TODO: delete call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            ! TODO: delete mix%material(2)%pe = mix%material(2)%pe + dum*(tmp - mix%material(2)%pe)

       !     do j = 1,9
       !         tmp = mix%material(1)%g(:,:,:,j)
       !         call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
       !         mix%material(1)%g(:,:,:,j) = mix%material(1)%g(:,:,:,j) + dum*(tmp - mix%material(1)%g(:,:,:,j))

       !         tmp = mix%material(2)%g(:,:,:,j)
       !         call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
       !         mix%material(2)%g(:,:,:,j) = mix%material(2)%g(:,:,:,j) + dum*(tmp - mix%material(2)%g(:,:,:,j))

                ! TODO: delete tmp = mix%material(1)%g_t(:,:,:,j)
                ! TODO: delete call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                ! TODO: delete mix%material(1)%g_t(:,:,:,j) = mix%material(1)%g_t(:,:,:,j) + dum*(tmp - mix%material(1)%g_t(:,:,:,j))

                ! TODO: delete tmp = mix%material(2)%g_t(:,:,:,j)
                ! TODO: delete call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                ! TODO: delete mix%material(2)%g_t(:,:,:,j) = mix%material(2)%g_t(:,:,:,j) + dum*(tmp - mix%material(2)%g_t(:,:,:,j))

                ! TODO: delete tmp = mix%material(1)%g_p(:,:,:,j)
                ! TODO: delete call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                ! TODO: delete mix%material(1)%g_p(:,:,:,j) = mix%material(1)%g_p(:,:,:,j) + dum*(tmp - mix%material(1)%g_p(:,:,:,j))

                ! TODO: delete tmp = mix%material(2)%g_p(:,:,:,j)
                ! TODO: delete call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                ! TODO: delete mix%material(2)%g_p(:,:,:,j) = mix%material(2)%g_p(:,:,:,j) + dum*(tmp - mix%material(2)%g_p(:,:,:,j))
       !     end do

            !mca add for stability

            !tmp = T
            !call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            !T = T + dum*(tmp - T)

        !    tmp = mix%material(1)%T
        !    call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
        !    mix%material(1)%T = mix%material(1)%T + dum*(tmp - mix%material(1)%T)

        !    tmp = mix%material(2)%T
        !    call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
        !    mix%material(2)%T = mix%material(2)%T + dum*(tmp - mix%material(2)%T)

        !    tmp = mix%material(1)%Ys
        !    call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
        !    mix%material(1)%Ys = mix%material(1)%Ys + dum*(tmp - mix%material(1)%Ys)

        !    tmp = mix%material(2)%Ys
        !    call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
        !    mix%material(2)%Ys = mix%material(2)%Ys + dum*(tmp - mix%material(2)%Ys)

        !end do

    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use SolidMixtureMod,  only: solid_mixture

    use ShearLayer4Mode_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix
    integer                                     :: imin, ind(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        ! ! determine interface velocity
        ! ind = minloc(abs(mix%material(1)%VF(:,1,1)-0.5d0))
        ! imin = ind(1)
        ! !vfdiff = mix%material(1)%VF(imin,1,1) - half
        ! !do i=1,size(mix%material(1)%VF(:,1,1))
        ! !  vfdiffloc = mix%material(1)%VF(i,1,1) - half
        ! !  if
        ! write(975,*) tsim, x(imin,1,1), u(imin,1,1)

    end associate
end subroutine

subroutine hook_mixture_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,&
                                mom_index,TE_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use ShearLayer4Mode_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    type(solid_mixture),             intent(in)    :: mix

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
    end associate
end subroutine

subroutine hook_material_g_source(decomp,hydro,elastic,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid

    use ShearLayer4Mode_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

end subroutine

subroutine hook_material_mass_source(decomp,hydro,elastic,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid

    use ShearLayer4Mode_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine

subroutine hook_material_energy_source(decomp,hydro,elastic,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid

    use ShearLayer4Mode_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine

subroutine hook_material_VF_source(decomp,hydro,elastic,x,y,z,tsim,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid

    use ShearLayer4Mode_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
