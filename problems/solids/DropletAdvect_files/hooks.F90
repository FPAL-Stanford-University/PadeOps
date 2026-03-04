module DropAdvect_data
    use kind_parameters,  only: rkind
    use constants,        only: one,two,eight,three,six,sixth,zero
    use FiltersMod,       only: filters
    use mpi
    implicit none

    real(rkind) :: p_infty = one, Rgas = one, gamma = 1.4_rkind, mu = 10._rkind, rho_0 = one, p_amb = 0.1_rkind
    real(rkind) :: p_infty_2 = one, Rgas_2 = one, gamma_2 = 1.4_rkind, mu_2 = 10._rkind, rho_0_2 = one, eta_det_ge = one,eta_det_ge_2 = one, eta_det_gp = one,eta_det_gp_2 = one, eta_det_gt = one,eta_det_gt_2 = one,diff_c_ge = one,diff_c_ge_2 = one, diff_c_gp = one,diff_c_gp_2 = one, diff_c_gt = one,diff_c_gt_2 = one
    real(rkind) :: minVF = 0.2_rkind, thick = one
    real(rkind) :: rhoRatio = one, pRatio = two, p_ten = one
    logical     :: sharp = .FALSE.
    real(rkind) :: p1,p2,rho1,rho2,u1,u2,g11_1,g11_2,grho1,grho2,a1,a2
    real(rkind) :: rho1_2,rho2_2,u1_2,u2_2,g11_1_2,g11_2_2,grho1_2,grho2_2,a1_2,a2_2
    real(rkind) :: rhoL, rhoR, YsL, YsR, VFL, VFR
    real(rkind) :: yield = one, yield2 = one, eta0k = 0.4_rkind, a_ratio = 1.0_rkind
    real(rkind) :: melt_t = one, melt_c = one, melt_t2 = one, melt_c2 = one
    real(rkind) :: kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e
    real(rkind) :: kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2
    integer     :: kos_sh,kos_sh2
    logical     :: explPlast = .FALSE., explPlast2 = .FALSE.
    logical     :: plastic = .FALSE., plastic2 = .FALSE.
    real(rkind) :: Ly = 1, Lx = 1, interface_init = 0.75_rkind,Tp = 4d0, shock_init = 0.6_rkind, kwave = 4.0_rkind,  v0 = 1d0, v0_2 = 1d0, tau0 =1d0, R = 1
    integer     :: pointx = 1, pointy = 1
    real(rkind) :: tau0_2=1d0, Nvel=1d0, minYs,etasize=1d0, ksize =1d0, delta_rho = 1d0, Nrho = 1d0, delta = 1d0


    type(filters) :: mygfil

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
    use constants,        only: one, half
    use decomp_2d,        only: decomp_info
    use exits,            only: warning

    use DropAdvect_data

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
        dy = Ly/real(ny,rkind)
        dz = dx

        if(abs(dx-dy)>1.0d-13) then
          call warning("dx not equal to dy")
        endif

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1  + i - 1, rkind ) * dx -2.0  ! x \in (-2,4]
                    y(i,j,k) = real( iy1  + j - 1, rkind ) * dy - 2.0
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,der,derStagg,interpMid,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
    use SolidGrid,        only: u_index,v_index,w_index,rho_index,uref_index,p_index
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
    use DropAdvect_data

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

     real(rkind),allocatable, dimension(:,:) :: p,Fx,Fy
    integer :: ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp2,tmp, dum, eta, yphys
    real(rkind), dimension(8) :: fparams
    real(rkind) :: fac,Lr, STRETCH_RATIO = 2.0
    integer, dimension(2) :: iparams
    real(rkind) :: a0, a0_2, vc
    logical :: adjustRgas = .TRUE.   ! If true, Rgas is used, Rgas2 adjusted to ensure p-T equilibrium
    logical :: adjustPamb = .FALSE.   ! If true, p_amb is adjusted to ensure p-T equilibrium	
    integer :: nx,ny,nz,ix
    integer :: ierr, rank,fh, filesize, chunksize, offset, offset2,totalproc
    integer, allocatable :: data(:), recvbuf(:)
    character(len=12) :: filename
    !nteger(kind=MPI_OFFSET_KIND) :: disp
    !nteger(kind=MPI_STATUS_SIZE) :: status(MPI_STATUS_SIZE)
    logical :: flag
    integer :: j, k, ios, iunit
    character(len=256) :: infile
    integer :: xs, xe    ! local x bounds from 2decomp
    real(rkind), allocatable :: row(:)
    integer :: ys, ye
    integer :: ix_global, local_ix
    integer :: nx_global, ny_global, nz_local
    integer :: ix_start, nx_local,ix_end

    ! Initialize MPI
    !all MPI_Init(ierr)
!    call MPI_Comm_RANK(MPI_COMM_WORLD, rank, ierr)
!    call MPI_Comm_SIZE(MPI_COMM_WORLD, totalproc, ierr)

    namelist /PROBINPUT/  p_infty,p_ten, Rgas, gamma, mu, rho_0, p_amb, thick, minVF, rhoRatio, pRatio, &
                          p_infty_2, Rgas_2, gamma_2, mu_2, rho_0_2, plastic, explPlast, yield,   &
                          plastic2, explPlast2, yield2, interface_init, kwave, a_ratio,&
                          melt_t, melt_c, melt_t2, melt_c2,&
                          kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e,kos_sh, &
                          kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2,kos_sh2, &
                          eta_det_ge,eta_det_ge_2,eta_det_gp,eta_det_gp_2,eta_det_gt,eta_det_gt_2, &
                          diff_c_ge,diff_c_ge_2,diff_c_gp,diff_c_gp_2,diff_c_gt,diff_c_gt_2, v0, v0_2, tau0, &
			tau0_2, eta0k, Nvel, etasize, ksize, delta_rho, Nrho, delta, Tp, thick, pointx, pointy

    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)
    allocate(Fx(pointx, pointy))
    ! Initialize mygfil
    call mygfil%init(                        decomp, &
                     .FALSE.,     .TRUE.,    .TRUE., &
                  "gaussian", "gaussian", "gaussian" )

    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index), &
                 p => fields(:,:,:,p_index), rho => fields(:,:,:,rho_index), x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )        
        if (mix%ns /= 2) then
            call GracefulExit("Number of species must be 2 for this problem. Check the input file.",928)
        end if

     nx = size(mesh,1); ny = size(mesh,2); nz = size(mesh,3)
     !   if(rhoRatio > 0) then
         ! if rhoRatio is positive, only rho_0 is different. Rgas is set such
          ! that Temperature equilibrium condition is satisfied
        !  gamma_2 = gamma; Rgas_2 = Rgas/rhoRatio; p_infty_2 = p_infty; 
        !  rho_0_2 = rho_0*rhoRatio; mu_2 = mu
       ! else
          ! if rhoRatio is negative, all quantities except Rgas need to be
          ! specified in input file. Rgas is then set such
          ! that Temperature equilibrium condition is satisfied
          if(adjustRgas) Rgas_2 = Rgas * (p_amb+p_infty_2)/(p_amb+p_infty)*rho_0/rho_0_2

          ! determine p_amb that guarantees T equilibrium
        !  if(adjustPamb) then
        !    fac = Rgas_2*rho_0_2/Rgas/rho_0
        !    p_amb = (fac*p_infty - p_infty_2)/(one - fac)
        !  endif
       ! endif
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
        ! call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield,1.0D-10)) !mca: see Sep1SolidEOS.F90 "init"
        ! call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,1.0D-10))
        call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield,1.0D-10,eta_det_ge,eta_det_gp,eta_det_gt,diff_c_ge,diff_c_gp,diff_c_gt,melt_t,melt_c,kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e,kos_sh,nx,ny,nz)) !mca: see Sep1SolidEOS.F90 "init"
        call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,1.0D-10,eta_det_ge_2,eta_det_gp_2,eta_det_gt_2,diff_c_ge_2,diff_c_gp_2,diff_c_gt_2,melt_t2,melt_c2,kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2,kos_sh2,nx,ny,nz))


        ! set logicals for plasticity
        mix%material(1)%plast = plastic ; mix%material(1)%explPlast = explPlast
        mix%material(2)%plast = plastic2; mix%material(2)%explPlast = explPlast2



        ! speed of sound
        a1 = sqrt((gamma*(p_amb+p_ten+p_infty) + 4.0d0/3.0d0*mu)/rho1)
        a2 = sqrt((gamma_2*(p_amb+p_infty_2) + 4.0d0/3.0d0*mu_2)/rho2)
        tmp = half * ( one - erf((625.0_rkind/7921.0_rkind - (x-0.5)*(x-0.5) - (y-0.5)*(y-0.5))/(thick*dx) ) )

        

        !where( (eta .le. (R)**2  )
        mix%material(1)%VF = (one-2*minVF)*tmp+minVF
        mix%material(1)%VF = one - mix%material(2)%VF
	!endwhere 
               
	!Set density profile and mass fraction based on volume fraction
	rho = rho_0*mix%material(1)%VF + rho_0_2*mix%material(2)%VF
	mix%material(2)%Ys = mix%material(2)%VF * rho_0_2 / rho
	mix%material(1)%Ys = one - mix%material(2)%Ys ! Enforce sum to unity

        minYs = minVF*rho_0


	u = v0
	v = 0
	w = 0


        !set mixture pressure (uniform)
	mix%material(1)%p =  p_amb
	mix%material(2)%p = mix%material(1)%p
        p = mix%material(2)%p
       !mix%surfaceTension_f(:,:,1,1) = Fx
       !mix%surfaceTension_f(:,:,1,2) = Fy
        ! Set initial values of g (inverse deformation gradient)
        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one

        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
	mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one
end associate
end subroutine
subroutine get_sponge(decomp,dx,dy,dz,mesh,fields,mix,rhou,rhov,rhow,rhoe,sponge)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
    use SolidGrid,        only: u_index,v_index,w_index,rho_index,e_index
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit
    use SolidMixtureMod,  only: solid_mixture
    use MultiphaseAdvection_data

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
    real(rkind) :: fac, Lr, STRETCH_RATIO = 5.0, int_KE
    integer, dimension(2) :: iparams
    real(rkind) :: a0, a0_2, sigma1, sigma2
    integer :: nx,ny,nz,k,ix,j
    integer :: ierr, rank,fh, filesize, chunksize, offset, offset2,totalproc
    integer, allocatable :: data(:), recvbuf(:)

    sponge = 0.0
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
    use operators,        only: grady,divergenceFV,interpolateFV,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_x,gradFV_y,gradFV_z
    use DerivativesMod,   only: derivatives
    use DerivativesStaggeredMod, only: derivativesStagg
    use InterpolatorsMod,        only: interpolators
    use reductions,       only: P_SUM, P_MEAN, P_MAXVAL, P_MINVAL
    use MultiphaseAdvection_data

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
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp,dum, eta, eta2, yphys, u_perturb, KE
    real(rkind) :: fac, Lr, STRETCH_RATIO = 5.0, int_KE
    integer, dimension(2) :: iparams
    real(rkind) :: a0, a0_2,dx1
    logical :: adjustRgas = .TRUE.   ! If true, Rgas is used, Rgas2 adjusted toensure p-T equilibrium
    logical :: adjustPamb = .FALSE.   ! If true, p_amb is adjusted to ensure p-T equilibrium

    integer :: nx,ny,nz,k,ix,j
    integer :: ierr, rank,fh, filesize, chunksize, offset, offset2,totalproc
    integer, allocatable :: data(:), recvbuf(:)
    character(len=50) :: filename,phiIname,phiRname, DphiIname, DphiRname
    character(len=50) :: rhoRname, rhoIname, pIname,pRname
    !nteger(kind=MPI_OFFSET_KIND) :: disp
    !nteger(kind=MPI_STATUS_SIZE) :: status(MPI_STATUS_SIZE)
    logical :: flag

end subroutine

subroutine  hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: filter3D

    use DropAdvect_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc
    
    integer :: ny,nx, i, j
    real(rkind) :: dx, xspng, tspng, xspngL, xspngR,STRETCH_RATIO = 2.0, Lr
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp,dumL, dumR,  dum, yphys
    
    nx = decomp%ysz(1)
    ny = decomp%ysz(2)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        !!! Hack to stop liquid's g from blowing up

        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one

        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

    !    if(decomp%yst(1)==1) then
    !      if(x_bc(1)==0) then
    !          rho( 1,:,:) = rhoL
    !          u  ( 1,:,:) = zero
    !          v  ( 1,:,:) = 0
    !          w  ( 1,:,:) = zero
    !          mix%material(1)%p(1,:,:) = p_amb
    !          mix%material(2)%p(1,:,:) = p_amb
!        !      
    !          mix%material(1)%VF ( 1,:,:) = minVF
    !          mix%material(2)%VF ( 1,:,:) = one - minVF
    !          mix%material(1)%Ys ( 1,:,:) = minYs
    !          mix%material(2)%Ys ( 1,:,:) = one - minYs
    !      end if
    !    endif

    !    if(decomp%yen(1)==decomp%xsz(1)) then
    !      if(x_bc(2)==0) then
    !          rho( nx,:,:) = rhoR
    !          u  ( nx,:,:) = zero
    !          v  ( nx,:,:) = 0
    !          w  ( nx,:,:) = zero
    !          mix%material(1)%p(nx,:,:) = p_amb
    !          mix%material(2)%p(nx,:,:) = p_amb

    !          mix%material(1)%VF ( nx,:,:) = minVF
    !          mix%material(2)%VF ( nx,:,:) = one - minVF
    !          mix%material(1)%Ys ( nx,:,:) = minYs
    !          mix%material(2)%Ys ( nx,:,:) = one - minYs
    !      end if
    !    endif

    !    if(decomp%yst(2)==1) then
    !      if(y_bc(1)==0) then
    !          rho( :,1,:) = rhoL*(one-minYs)
    !          u  ( :,1,:) = zero
    !          v  ( :,1,:) = 0
    !          w  ( :,1,:) = zero
    !          mix%material(1)%p(:,1,:) = p_amb
    !          mix%material(2)%p(:,1,:) = p_amb
        !
    !          mix%material(1)%VF ( :,1,:) = minVF
    !          mix%material(2)%VF ( :,1,:) = one - minVF
    !          mix%material(1)%Ys ( :,1,:) = minYs
    !          mix%material(2)%Ys ( :,1,:) = one - minYs
    !      end if
    !    end if

    !   if(decomp%yen(2)==decomp%ysz(2)) then
    !      if(y_bc(2)==0) then
    !          rho( :,ny,:) = rhoR*(one-minYs)
    !          u  ( :,ny,:) = zero
    !          v  ( :,ny,:) = zero
    !          w  ( :,ny,:) = zero
    !          mix%material(1)%p(:,ny,:) = p_amb
    !          mix%material(2)%p(:,ny,:) = p_amb
       !
    !          mix%material(1)%VF ( :,ny,:) = minVF
    !          mix%material(2)%VF ( :,ny,:) = one - minVF
    !          mix%material(1)%Ys ( :,ny,:) = minYs
    !          mix%material(2)%Ys ( :,ny,:) = one - minYs
    !      end if
    !    endif




        !xspng = -two + half
        !tspng = 0.2_rkind
        !dx = x(2,1,1) - x(1,1,1)
        !dum = half*(one - tanh( (x-xspng)/(tspng) ))

        yphys = atanh(2.0*y /(1 + 1/STRETCH_RATIO))
        Lr    = Lx/(yphys(1,ny,1) - yphys(1,1,1))
        yphys = y !Lr*yphys

        xspngL = -2.0  + 0.125
        xspngR =  2.0 - 0.125 
        dx = x(2,1,1) - x(1,1,1)
        tspng = 0.006
        dumL = half*(one - tanh((yphys - xspngL)/(tspng) )) + half*(one - tanh((x - xspngL)/(tspng) ))
        dumR = half*(one + tanh((yphys-xspngR)/(tspng) )) + half*(one + tanh((x - xspngR)/(tspng) ))      
        dum  = dumL+dumR


     !   do i=1,4
     !       tmp = u
     !       call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
     !       u = u + dum*(tmp - u)

     !       tmp = v
     !       call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
     !       v = v + dum*(tmp - v)

     !       tmp = w
     !       call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
     !       w = w + dum*(tmp - w)

     !      tmp = e
     !       call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
     !       e = e + dum*(tmp - e)

     !       tmp = rho
     !       call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
     !       rho = rho + dum*(tmp - rho)

     !       tmp = mix%material(1)%p
     !       call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
     !       mix%material(1)%p = mix%material(1)%p + dum*(tmp - mix%material(1)%p)

     !       tmp = mix%material(2)%p
     !       call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
     !       mix%material(2)%p = mix%material(2)%p + dum*(tmp - mix%material(2)%p)

!            tmp = mix%material(1)%pe
!            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
!            mix%material(1)%pe = mix%material(1)%pe + dum*(tmp - mix%material(1)%pe)

!            tmp = mix%material(2)%pe
!            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
!            mix%material(2)%pe = mix%material(2)%pe + dum*(tmp - mix%material(2)%pe)

           ! do j = 1,9
           !     tmp = mix%material(1)%g(:,:,:,j)
           !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
           !     mix%material(1)%g(:,:,:,j) = mix%material(1)%g(:,:,:,j) + dum*(tmp - mix%material(1)%g(:,:,:,j))

           !     tmp = mix%material(2)%g(:,:,:,j)
           !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
          !      mix%material(2)%g(:,:,:,j) = mix%material(2)%g(:,:,:,j) + dum*(tmp - mix%material(2)%g(:,:,:,j))

           !     tmp = mix%material(1)%g_t(:,:,:,j)
           !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
           !     mix%material(1)%g_t(:,:,:,j) = mix%material(1)%g_t(:,:,:,j) + dum*(tmp - mix%material(1)%g_t(:,:,:,j))

            !    tmp = mix%material(2)%g_t(:,:,:,j)
            !    call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            !    mix%material(2)%g_t(:,:,:,j) = mix%material(2)%g_t(:,:,:,j) + dum*(tmp - mix%material(2)%g_t(:,:,:,j))

             !   tmp = mix%material(1)%g_p(:,:,:,j)
             !   call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
             !   mix%material(1)%g_p(:,:,:,j) = mix%material(1)%g_p(:,:,:,j) + dum*(tmp - mix%material(1)%g_p(:,:,:,j))

              !  tmp = mix%material(2)%g_p(:,:,:,j)
              !  call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
              !  mix%material(2)%g_p(:,:,:,j) = mix%material(2)%g_p(:,:,:,j) + dum*(tmp - mix%material(2)%g_p(:,:,:,j))
            !end do

            !mca add for stability

!            tmp = T
!            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
!            T = T + dum*(tmp - T)

!            tmp = mix%material(1)%T
!            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
!            mix%material(1)%T = mix%material(1)%T + dum*(tmp - mix%material(1)%T)

!            tmp = mix%material(2)%T
!            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
!            mix%material(2)%T = mix%material(2)%T + dum*(tmp - mix%material(2)%T)

      !      tmp = mix%material(1)%Ys
      !     call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
      !      mix%material(1)%Ys = mix%material(1)%Ys + dum*(tmp - mix%material(1)%Ys)

      !      tmp = mix%material(2)%Ys
      !      call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
      !      mix%material(2)%Ys = mix%material(2)%Ys + dum*(tmp - mix%material(2)%Ys)

      !      tmp = mix%material(1)%VF
      !      call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
      !      mix%material(1)%VF = mix%material(1)%VF + dum*(tmp - mix%material(1)%VF)

      !      tmp = mix%material(2)%VF
      !      call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
      !      mix%material(2)%VF = mix%material(2)%VF + dum*(tmp - mix%material(2)%VF)


      !  end do

       ! if(decomp%yen(1)==decomp%xsz(1)) then
          !if(x_bc(2)==0) then
            !rho(nx,:,:) = rhoR ! rho(nx-1,:,:)
            !u  (nx,:,:) = zero ! zero
            !v  (nx,:,:) = zero ! v(nx-1,:,:)
            !w  (nx,:,:) = zero ! w(nx-1,:,:)
           ! mix%material(1)%p  (nx,:,:) = p1 ! mix%material(1)%p(nx-1,:,:)
            !mix%material(2)%p  (nx,:,:) = p1 ! mix%material(2)%p(nx-1,:,:)
            
            !mix%material(1)%g11(nx,:,:) = one;  mix%material(1)%g12(nx,:,:) = zero; mix%material(1)%g13(nx,:,:) = zero
            !mix%material(1)%g21(nx,:,:) = zero; mix%material(1)%g22(nx,:,:) = one;  mix%material(1)%g23(nx,:,:) = zero
            !mix%material(1)%g31(nx,:,:) = zero; mix%material(1)%g32(nx,:,:) = zero; mix%material(1)%g33(nx,:,:) = one
  
           ! mix%material(2)%g11(nx,:,:) = one;  mix%material(2)%g12(nx,:,:) = zero; mix%material(2)%g13(nx,:,:) = zero
           ! mix%material(2)%g21(nx,:,:) = zero; mix%material(2)%g22(nx,:,:) = one;  mix%material(2)%g23(nx,:,:) = zero
           ! mix%material(2)%g31(nx,:,:) = zero; mix%material(2)%g32(nx,:,:) = zero; mix%material(2)%g33(nx,:,:) = one

            !mix%material(1)%gt11(nx,:,:) = one;  mix%material(1)%gt12(nx,:,:) = zero; mix%material(1)%gt13(nx,:,:) = zero
            !mix%material(1)%gt21(nx,:,:) = zero; mix%material(1)%gt22(nx,:,:) = one;  mix%material(1)%gt23(nx,:,:) = zero
            !mix%material(1)%gt31(nx,:,:) = zero; mix%material(1)%gt32(nx,:,:) = zero; mix%material(1)%gt33(nx,:,:) = one
            
            !mix%material(1)%gp11(nx,:,:) = one;  mix%material(1)%gp12(nx,:,:) = zero; mix%material(1)%gp13(nx,:,:) = zero
            !mix%material(1)%gp21(nx,:,:) = zero; mix%material(1)%gp22(nx,:,:) = one;  mix%material(1)%gp23(nx,:,:) = zero
            !mix%material(1)%gp31(nx,:,:) = zero; mix%material(1)%gp32(nx,:,:) = zero; mix%material(1)%gp33(nx,:,:) = one
            
            !mix%material(1)%pe = zero


           ! mix%material(2)%gt11(nx,:,:) = one;  mix%material(2)%gt12(nx,:,:) = zero; mix%material(2)%gt13(nx,:,:) = zero
           ! mix%material(2)%gt21(nx,:,:) = zero; mix%material(2)%gt22(nx,:,:) = one;  mix%material(2)%gt23(nx,:,:) = zero
           ! mix%material(2)%gt31(nx,:,:) = zero; mix%material(2)%gt32(nx,:,:) = zero; mix%material(2)%gt33(nx,:,:) = one
            
          !  mix%material(2)%gp11(nx,:,:) = one;  mix%material(2)%gp12(nx,:,:) = zero; mix%material(2)%gp13(nx,:,:) = zero
          !  mix%material(2)%gp21(nx,:,:) = zero; mix%material(2)%gp22(nx,:,:) = one;  mix%material(2)%gp23(nx,:,:) = zero
         !   mix%material(2)%gp31(nx,:,:) = zero; mix%material(2)%gp32(nx,:,:) = zero; mix%material(2)%gp33(nx,:,:) = one
            
        !    mix%material(2)%pe = zero
            
            ! mix%material(1)%Ys (nx,:,:) = YsR
            ! mix%material(2)%Ys (nx,:,:) = one - YsR
            
       !     mix%material(1)%VF (nx,:,:) = VFR
       !     mix%material(2)%VF (nx,:,:) = one - VFR
      !   endif
     ! endif

    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use SolidMixtureMod,  only: solid_mixture
    use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps

    use DropAdvect_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(in) :: mix
    integer                                     :: imin, ind(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
    end associate
end subroutine

subroutine hook_mixture_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,&
                                mom_index,TE_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use DropAdvect_data

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

    use DropAdvect_data

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

    use DropAdvect_data

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

    use DropAdvect_data

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

    use DropAdvect_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
