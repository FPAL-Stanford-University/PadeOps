module ShearLayerLinear_data
    use kind_parameters,  only: rkind
    use constants,        only: one,two,eight,three,six,sixth,zero, pi
    use FiltersMod,       only: filters
    use DerivativesMod,   only: derivatives
    use DerivativesStaggeredMod, only: derivativesStagg
    use InterpolatorsMod,        only: interpolators
    use mpi 
    implicit none

    real(rkind) :: p_infty = one, Rgas = one, gamma = 1.4_rkind, mu = 10._rkind, rho_0 = one, p_amb = 0.1_rkind
    real(rkind) :: p_infty_2 = one, Rgas_2 = one, gamma_2 = 1.4_rkind, mu_2 = 10._rkind, rho_0_2 = one, eta_det_ge = one,eta_det_ge_2 = one, eta_det_gp = one,eta_det_gp_2 = one, eta_det_gt = one,eta_det_gt_2 = one,diff_c_ge = one,diff_c_ge_2 = one, diff_c_gp = one,diff_c_gp_2 = one, diff_c_gt = one,diff_c_gt_2 = one
    real(rkind) :: minVF = 0.2_rkind, thick = 0.01, p_ten = one
    logical     :: sharp = .FALSE.
    real(rkind) :: p1,p2,rho1,rho2,u1,u2,g11_1,g11_2,grho1,grho2,a1,a2
    real(rkind) :: rho1_2,rho2_2,u1_2,u2_2,g11_1_2,g11_2_2,grho1_2,grho2_2,a1_2,a2_2
    real(rkind) :: rhoL, rhoR, YsL, YsR, VFL, VFR, vL, vR, uL, uR
    real(rkind) :: yield = one, yield2 = one, eta0k = 0.4_rkind
    real(rkind) :: melt_t = one, melt_c = one, melt_t2 = one, melt_c2 = one
    real(rkind) :: kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e, alpha3, alpha4,alpha2, alpha
    real(rkind) :: kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2, v_disturb
    real(rkind) :: v0=zero, v0_2=zero, tau0=1d-14, tau0_2=1d-14, Nrho = 1, U0 = zero, m = 1, p_mu = 1, p_mu2 = 1, epsilonk = 0
    integer     :: kos_sh,kos_sh2,pointy, pointx
    logical     :: explPlast = .FALSE., explPlast2 = .FALSE.
    logical     :: plastic = .FALSE., plastic2 = .FALSE.
    real(rkind) :: Ly = 1.0, Lx = 2*pi, interface_init = 10d-3, kwave = 4.0_rkind, ksize = 10d0, etasize = 0.5d0, delta_d = 0.0125D0, delta = 0.0125D0, delta_rho = 0.0125D0 

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

    use ShearLayerLinear_data

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
        dz = dx

        if(abs(dx-dy)>1.0d-13) then
          call warning("dx not equal to dy")
        endif

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1   + i - 1, rkind ) * dx -pi
                    y(i,j,k) = real( iy1 - 1  + j - 1, rkind ) * dy - 0.5
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,der,derStagg,interpMid,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz,periodicx,periodicy,periodicz, x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
    use SolidGrid,        only: u_index,v_index,w_index,rho_index, uref_index
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: grady, divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_x, gradFV_y, gradFV_z
    use DerivativesMod,   only: derivatives  
    use DerivativesStaggeredMod, only: derivativesStagg
    use InterpolatorsMod,        only: interpolators
    use reductions,       only: P_SUM, P_MEAN, P_MAXVAL, P_MINVAL 
    use ShearLayerLinear_data

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
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum, eta, eta2, yphys, u_perturb, KE
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2)) :: v_perturb
    real(rkind), dimension(decomp%ysz(2),4) :: phi_i,phi_r, p_perturbr, p_perturbi, Dphi_i, Dphi_r
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2),decomp%ysz(3),4) :: phi_i3, phi_r3,phi_i_int, phi_r_int
    real(rkind), dimension(8) :: fparams
    real(rkind), dimension(4) :: alphai
    real(rkind) :: fac, Lr, STRETCH_RATIO = 10.0, int_KE
    integer, dimension(2) :: iparams
    real(rkind) :: a0, a0_2
    logical :: adjustRgas = .TRUE.   ! If true, Rgas is used, Rgas2 adjusted to ensure p-T equilibrium
    logical :: adjustPamb = .FALSE.   ! If true, p_amb is adjusted to ensure p-T equilibrium

    integer :: nx,ny,nz,k,ix,j
    integer :: ierr, rank,fh, filesize, chunksize, offset, offset2,totalproc
    integer, allocatable :: data(:), recvbuf(:)
    character(len=12) :: filename
    !nteger(kind=MPI_OFFSET_KIND) :: disp
    !nteger(kind=MPI_STATUS_SIZE) :: status(MPI_STATUS_SIZE)
    logical :: flag

    ! Initialize MPI
    !all MPI_Init(ierr)
    call MPI_Comm_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_SIZE(MPI_COMM_WORLD, totalproc, ierr)
    nx = size(mesh,1); ny = size(mesh,2); nz = size(mesh,3)

    namelist /PROBINPUT/  p_infty, Rgas, gamma, mu, rho_0, p_amb, thick, minVF,  &
                          p_infty_2, Rgas_2, gamma_2, mu_2, rho_0_2, plastic, explPlast, yield,   &
                          plastic2, explPlast2, yield2, interface_init, kwave,delta, delta_d, delta_rho, &
                          melt_t, melt_c, melt_t2, melt_c2, &
                          kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e,kos_sh, &
                          kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2,kos_sh2, &
                          eta_det_ge,eta_det_ge_2,eta_det_gp,eta_det_gp_2,eta_det_gt,eta_det_gt_2, &
                          diff_c_ge,diff_c_ge_2,diff_c_gp,diff_c_gp_2,diff_c_gt,diff_c_gt_2,alpha, &
                          v0, alpha4, v_disturb, alpha3, alpha2,v0_2, tau0, tau0_2, eta0k, ksize, etasize, p_mu, p_mu2, Nrho,pointy, pointx, epsilonk        
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    ! Initialize mygfil
    call mygfil%init(                        decomp, &
                     .FALSE.,     .TRUE.,    .TRUE., &
                  "gaussian", "gaussian", "gaussian" )

    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index), &
                 uref => fields(:,:,:,uref_index),                                                    &
                 rho => fields(:,:,:,rho_index), x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        if (mix%ns /= 2) then
            call GracefulExit("Number of species must be 2 for this problem. Check the input file.",928)
        end if
        
        !Ensure temperature equilibrium at start
        if(adjustRgas) Rgas_2 = Rgas * (p_amb+p_infty_2)/(p_amb+p_infty)*rho_0/rho_0_2

        ! speed of sound
        a0   = sqrt((gamma*(p_amb+p_infty) + 4.0d0/3.0d0*mu)/rho_0)
        a0_2   = sqrt((gamma*(p_amb+p_infty_2) + 4.0d0/3.0d0*mu_2)/rho_0_2)

        ! write material properties
        if (nrank == 0) then
            print *, '---Material 1---'
            write(*,'(3(a,e12.5))') 'rho_0 = ', rho_0, ', gam  = ', gamma, ', p_infty = ', p_infty
            write(*,'(3(a,e12.5))') 'shearMod    = ', mu,    ', Rgas = ', Rgas, ', SOS = ', a0
            write(*,'(3(a,e12.5))') 'tau0  = ', tau0
            print *, '---Material 2---'
            write(*,'(3(a,e12.5))') 'rho_0 = ', rho_0_2, ', gam  = ', gamma_2, ', p_infty = ', p_infty_2
            write(*,'(3(a,e12.5))') 'shearMod    = ', mu_2,    ', Rgas = ', Rgas_2, ', SOS = ', a0_2
            write(*,'(3(a,e12.5))') 'tau0  = ', tau0_2
            write(*,*) 'p_amb = ', p_amb
        end if

        ! Set materials
        !call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield,tau0))
        !TODO: delete call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield,tau0,eta_det_ge,eta_det_gp,eta_det_gt,diff_c_ge,diff_c_gp,diff_c_gt,melt_t,melt_c,kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e,kos_sh,nx,ny,nz)) !mca: see Sep1SolidEOS.F90 "init"
        !call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,tau0_2))
        !TODO: delete call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,tau0_2,eta_det_ge_2,eta_det_gp_2,eta_det_gt_2,diff_c_ge_2,diff_c_gp_2,diff_c_gt_2,melt_t2,melt_c2,kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2,kos_sh2,nx,ny,nz))
        call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield,1.0D-10,eta_det_ge,eta_det_gp,eta_det_gt,diff_c_ge,diff_c_gp,diff_c_gt,melt_t,melt_c,kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e,kos_sh,nx,ny,nz)) !mca: see Sep1SolidEOS.F90 "init"
        call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,1.0D-10,eta_det_ge_2,eta_det_gp_2,eta_det_gt_2,diff_c_ge_2,diff_c_gp_2,diff_c_gt_2,melt_t2,melt_c2,kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2,kos_sh2,nx,ny,nz))

        ! set logicals for plasticity
        mix%material(1)%plast = plastic ; mix%material(1)%explPlast = explPlast
        mix%material(2)%plast = plastic2; mix%material(2)%explPlast = explPlast2

        ! Set up smearing function for VF based on interface location and thickness
        !tmp = half * ( one - erf( (x-(interface_init+eta0k/(2.0_rkind*pi*kwave)*sin(2.0_rkind*kwave*pi*y)))/(thick*dx) ) )

        delta_rho = Nrho * dx * 0.275d0 !converts from Nrho to approximate thickness of erf profile
	!delta_rho = Nrho*0.275d0
        !eta =(y-interface_init) !-epsilonk*sin(alpha*x))
        yphys = atanh(2.0*y /(1 + 1/STRETCH_RATIO))
        Lr    = 12.0/(yphys(1,ny,1) - yphys(1,1,1))
        yphys = Lr*yphys
        

	tmp = (half ) * ( one - erf( (yphys)/(delta_rho) ) )

	!set mixture Volume fraction
	mix%material(1)%VF = minVF + (one-two*minVF)*tmp
	mix%material(2)%VF = 1 - mix%material(1)%VF

        !Set density profile and mass fraction based on volume fraction
	rho = rho_0*mix%material(1)%VF + rho_0_2*mix%material(2)%VF
	mix%material(1)%Ys = mix%material(1)%VF * rho_0 / rho
	mix%material(2)%Ys = one - mix%material(1)%Ys ! Enforce sum to unity



       ! eta0k = etasize*delta
       ! kwave = ksize*delta
        eta2 = yphys !eta !-delta_rho

        

        where(eta2 .ge. 0)
           u = v0_2*erf(eta2/delta)
        elsewhere(eta2 .lt. 0 )
           u = v0*erf(eta2/delta)
        endwhere
           uref = u
           v = 0
           w = 0

        
        open (unit=8, file="Phi_I.txt", status='old', action='read' )
        open (unit=12, file="Phi_R.txt", status='old', action='read' )
        open (unit=16, file="DPhi_I.txt", status='old', action='read' )
        open (unit=18, file="DPhi_R.txt", status='old', action='read' )

!        open (unit=20, file="Phi_I_2.txt", status='old', action='read' )
!        open (unit=22, file="Phi_R_2.txt", status='old', action='read' )
!        open (unit=24, file="DPhi_I_2.txt", status='old', action='read' )
!        open (unit=26, file="DPhi_R_2.txt", status='old', action='read' )

!        open (unit=28, file="Phi_I_3.txt", status='old', action='read' )
!        open (unit=30, file="Phi_R_3.txt", status='old', action='read' )
!        open (unit=32, file="DPhi_I_3.txt", status='old', action='read' )
!        open (unit=34, file="DPhi_R_3.txt", status='old', action='read' )

!        open (unit=36, file="Phi_I_4.txt", status='old', action='read' )
!        open (unit=38, file="Phi_R_4.txt", status='old', action='read' )
!        open (unit=40, file="DPhi_I_4.txt", status='old', action='read' )
!        open (unit=42, file="DPhi_R_4.txt", status='old', action='read' )

!        open (unit=44, file="Phi_I_5.txt", status='old', action='read' )
!        open (unit=46, file="Phi_R_5.txt", status='old', action='read' )
!        open (unit=48, file="DPhi_I_5.txt", status='old', action='read' )
!        open (unit=50, file="DPhi_R_5.txt", status='old', action='read' )

!        open (unit=52, file="Phi_I_6.txt", status='old', action='read' )
!        open (unit=54, file="Phi_R_6.txt", status='old', action='read' )
!        open (unit=56, file="DPhi_I_6.txt", status='old', action='read' )
!        open (unit=58, file="DPhi_R_6.txt", status='old', action='read' )

!        open (unit=60, file="Phi_I_7.txt", status='old', action='read' )
!        open (unit=62, file="Phi_R_7.txt", status='old', action='read' )
!        open (unit=64, file="DPhi_I_7.txt", status='old', action='read' )
!        open (unit=66, file="DPhi_R_7.txt", status='old', action='read' )

!        do i = 1,4

!           alphai(i) = i*0.5

!        enddo
 
        alphai(1) = 1.0                     
        do i = 1,pointy
          k = pointy - i + 1
          read(8,*) phi_i(i,1)
          read(12,*) phi_r(i,1)
          read(16,*) Dphi_i(i,1)
          read(18,*) Dphi_r(i,1)

         ! read(20,*) phi_i(i,2)
         ! read(22,*) phi_r(i,2)
         ! read(24,*) Dphi_i(i,2)
         ! read(26,*) Dphi_r(i,2)

         ! read(28,*) phi_i(i,3)
         ! read(30,*) phi_r(i,3)
         ! read(32,*) Dphi_i(i,3)
         ! read(34,*) Dphi_r(i,3)


         ! read(36,*) phi_i(i,4)
         ! read(38,*) phi_r(i,4)
         ! read(40,*) Dphi_i(i,4)
         ! read(42,*) Dphi_r(i,4)

          !read(44,*) phi_i(i,5)
          !read(46,*) phi_r(i,5)
          !read(48,*) Dphi_i(i,5)
          !read(50,*) Dphi_r(i,5)
           
          !read(52,*) phi_i(i,6)
          !read(54,*) phi_r(i,6)
          !read(56,*) Dphi_i(i,6)
          !read(58,*) Dphi_r(i,6)

          !read(60,*) phi_i(i,7)
          !read(62,*) phi_r(i,7)
          !read(64,*) Dphi_i(i,7)
          !read(66,*) Dphi_r(i,7)


        end do

        mix%material(1)%p  = p_amb

    !    do i = 1,2
    !       call interpolateFV_y(decomp,interpMid,phi_i3(:,:,:,i),phi_i_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
 
    !       call interpolateFV_y(decomp,interpMid,phi_r3(:,:,:,i),phi_r_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

    !       call gradFV_y(decomp,derStagg,phi_i_int,Dphi_i(:,:,:,i),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
    !       call gradFV_y(decomp,derStagg,phi_r_int,Dphi_r(:,:,:,i),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
    !    enddo
        v = 0
        u_perturb = 0
        j = 1
        do i = 1,pointy
          k = pointy+1 - i
        ! do j = 1,4
            v(:,i,:) = v(:,i,:) + alphai(j)*(phi_i(i,j)*cos(alphai(j)*x(:,i,:)) + phi_r(i,j)*sin(alphai(j)*x(:,i,:)) )
            u_perturb(:,i,:)  = u_perturb(:,i,:) - Dphi_i(i,j)*sin(alphai(j)*x(:,i,:)) + Dphi_r(i,j)*cos(alphai(j)*x(:,i,:))
        ! enddo

 
        enddo
       
        
        KE = 0.5*(v**2 + u_perturb**2)
        int_KE = P_SUM(KE/(nx*ny*nz))
        v = epsilonk*v /sqrt(int_KE)
        u = u + epsilonk*u_perturb /sqrt(int_KE)
        mix%material(2)%p  = mix%material(1)%p 


        ! Set initial values of g (inverse deformation gradient)
        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one

        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

        !Stuff for boundary conditions
        rhoL = 0.99d-1 !rho(1,1,1)
        rhoR = 1d-1 !rho(1,decomp%ysz(2),1)
        uL = -1*v0 !u(1,1,1)
        uR = v0_2 !u(1,decomp%ysz(2),1)
        vL = 0 !v(1,1,1)
        vR = 0 !v(1,decomp%ysz(2),1)
        YsL  = mix%material(1)%Ys(1,1,1)
        YsR  = mix%material(1)%Ys(1,decomp%ysz(2),1)
        VFL  = mix%material(1)%VF(1,1,1)
        VFR  = mix%material(1)%VF(1,decomp%ysz(2),1)
 
        
    end associate

end subroutine

subroutine get_sponge(decomp,dx,dy,dz,mesh,fields,mix,rhou,rhov,rhow,rhoe,sponge)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
    use SolidGrid,        only: u_index,v_index,w_index,rho_index,e_index
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit
    use SolidMixtureMod,  only: solid_mixture
    use ShearLayerLinear_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    type(solid_mixture),             intent(inout) :: mix
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout):: sponge
    real(rkind), dimension(2), intent(inout) :: rhou, rhov,rhow,rhoe
    integer :: ioUnit,i,iy
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp,dum, eta, eta2, yphys
    real(rkind) :: fac, Lr, STRETCH_RATIO = 10.0, int_KE
    integer, dimension(2) :: iparams
    real(rkind) :: a0, a0_2, sigma1, sigma2
    integer :: nx,ny,nz,k,ix,j
    integer :: ierr, rank,fh, filesize, chunksize, offset, offset2,totalproc
    integer, allocatable :: data(:), recvbuf(:) 


        associate(u => fields(:,:,:,u_index), v => fields(:,:,:,v_index),w => fields(:,:,:,w_index), rho => fields(:,:,:,rho_index), e => fields(:,:,:,e_index), x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )


        
        nx = size(mesh,1); ny = size(mesh,2); nz = size(mesh,3)
        yphys = atanh(2.0*y /(1 + 1/STRETCH_RATIO))
        Lr    = 12.0/(yphys(1,ny,1) - yphys(1,1,1))
        yphys = Lr*yphys


        sigma1 = -2400 !43028

        where(yphys .LE. -4.5)
           sponge(:,:,:,1) = sigma1*( (yphys + 4.5)/1.50)**2.0
        elsewhere
           sponge(:,:,:,1) = 0
        endwhere

        where(yphys .GE. 4.5)
           sponge(:,:,:,2) = sigma1*( (yphys- 4.5 )/1.50)**2.0
        elsewhere
           sponge(:,:,:,2) = 0
        endwhere

        rhou(1) = -0.099*v0 !rho(1,1,1)*u(1,1,1)
        rhou(2) = 1d-1*v0_2 !rho(1,ny,1)*u(1,ny,1)
        rhov(1) = 0 !-1.060981230880199d-5 !rho(1,1,1)*v(1,1,1)
        rhov(2) = 0 !2.175685479370164d-08
        rhow(1) = rho(1,1,1)*w(1,1,1)
        rhow(2) = rho(1,ny,1)*w(1,ny,1)
        rhoe(1) = 0.099*(92.840166369578 + 0.5*(v0**2)) !1.d0*(103.176 + 0.5*(v0**2))
        rhoe(2) = 1d-1*(25 + 0.5*(v0_2**2))

        do i = 1,2
          mix%material(i)%VF_ref(1) = mix%material(i)%VF(1,1,1)
          mix%material(i)%VF_ref(2) = mix%material(i)%VF(1,ny,1) 
          mix%material(i)%Ys_ref(1) = mix%material(i)%Ys(1,1,1)*rho(1,1,1)
          mix%material(i)%Ys_ref(2) = mix%material(i)%Ys(1,ny,1)*rho(1,ny,1)
        enddo

        print *, "rhou ", rhou
        print *, "rhov ", rhov
        print *, "rhow ", rhow
        print *, "rhoe ", rhoe
        print *, "VF ", mix%material(1)%VF_ref
        print *, "Ys ", mix%material(1)%Ys_ref
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

    use ShearLayerLinear_data

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
           write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ShearLayerLinear_"//trim(str)//"_", vizcount, ".dat"

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
       write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ShearLayerLinear_statistics.dat"

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

    use ShearLayerLinear_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc
    
    integer :: nx,ny, i, j
    real(rkind) :: dy, yspng, tspng, yspngR, yspngL, Lr, STRETCH_RATIO = 10.0
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum, dumL, dumR, yphys
    
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
              rho( :,1,:) = rhoL
              u  ( :,1,:) = uL
              v  ( :,1,:) = vL !v(:,2,:)
              w  ( :,1,:) = zero
              mix%material(1)%p(:,1,:) =  1.0 !mix%material(1)%p(:,2,:)
              mix%material(2)%p(:,1,:) =  1.0 !mix%material(2)%p(:,2,:)
 
              mix%material(1)%VF ( :,1,:) = VFL
              mix%material(2)%VF ( :,1,:) = one - VFL
              mix%material(1)%Ys ( :,1,:) = YsL
              mix%material(2)%Ys ( :,1,:) = one - YsL
          end if
        endif

          
        !print *, "decompyen", decomp%yen(2)
        !print *, "decompyst", decomp%yst(2)
        !!print *, "rhoL", rhoL
        !print *, "rhoR", rhoR
        !print *, "vl", vl
        !print *, "vr", vr
        !print *, "Ys", YsR
        !print *, "YsL", YsL
        !print *, "VFL", VFL
        !print *, "VFR", VFR 

        if(decomp%yen(2)==decomp%ysz(2)) then
          if(y_bc(2)==0) then
              rho( :,ny,:) = rhoR
              u  ( :,ny,:) = uR
              v  ( :,ny,:) = vR !v(:,ny-1,:)
              w  ( :,ny,:) = zero
              mix%material(1)%p(:,ny,:) =  1.0 !mix%material(1)%p(:,ny-1,:)
              mix%material(2)%p(:,ny,:) =  1.0 !mix%material(2)%p(:,ny-1,:)
              
              mix%material(1)%VF ( :,ny,:) = VFR
              mix%material(2)%VF ( :,ny,:) = one - VFR
              mix%material(1)%Ys ( :,ny,:) = YsR
              mix%material(2)%Ys ( :,ny,:) = one - YsR
          end if
        endif
     
        
        
  ! apply sponge at left and right boundaries to damp outgoing waves
        yphys = atanh(2.0*y /(1 + 1/STRETCH_RATIO))
        Lr    = 12.0/(yphys(1,ny,1) - yphys(1,1,1))
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

    use ShearLayerLinear_data

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

    use ShearLayerLinear_data

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

    use ShearLayerLinear_data

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

    use ShearLayerLinear_data

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

    use ShearLayerLinear_data

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

    use ShearLayerLinear_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
