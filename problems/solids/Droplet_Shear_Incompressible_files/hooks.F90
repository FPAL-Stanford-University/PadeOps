module Droplet_Shear_Incompressible_data
    use kind_parameters,  only: rkind
    use constants,        only: one,two,eight,three,six,sixth,zero
    use FiltersMod,       only: filters
    implicit none

    real(rkind) :: p_infty = one, Rgas = one, gamma = 1.4_rkind, mu = 10._rkind, rho_0 = one, p_amb = 0.1_rkind
    real(rkind) :: p_infty_2 = one, Rgas_2 = one, gamma_2 = 1.4_rkind, mu_2 = 10._rkind, rho_0_2 = one, eta_det_ge = one,eta_det_ge_2 = one, eta_det_gp = one,eta_det_gp_2 = one, eta_det_gt = one,eta_det_gt_2 = one,diff_c_ge = one,diff_c_ge_2 = one, diff_c_gp = one,diff_c_gp_2 = one, diff_c_gt = one,diff_c_gt_2 = one
    real(rkind) :: minVF = 0.2_rkind, thick = one
    real(rkind) :: rhoRatio = one, pRatio = two
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
    real(rkind) :: Ly = one, Lx = one, interface_init = 0.75_rkind,Tp = 4d0, shock_init = 0.6_rkind, kwave = 4.0_rkind,  v0 = 1d0, v0_2 = 1d0, tau0 =1d0, R = 0.15
    real(rkind) :: tau0_2=1d0, Nvel=1d0, etasize=1d0, ksize =1d0, delta_rho = 1d0, Nrho = 1d0, delta = 1d0


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

    use Droplet_Shear_Incompressible_data

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
                    x(i,j,k) = real( ix1     + i - 1, rkind ) * dx   ! x \in (-2,4]
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven,pi,eps
    use SolidGrid,        only: u_index,v_index,w_index, rho_index
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    use SolidMixtureMod,  only: solid_mixture
    
    use Droplet_Shear_Incompressible_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(inout) :: tstop, dt, tviz
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum, eta
    real(rkind), dimension(8) :: fparams
    real(rkind) :: fac
    integer, dimension(2) :: iparams
    real(rkind) :: a0, a0_2, vc
	logical :: adjustRgas = .TRUE.   ! If true, Rgas is used, Rgas2 adjusted to ensure p-T equilibrium
    logical :: adjustPamb = .FALSE.   ! If true, p_amb is adjusted to ensure p-T equilibrium
	
    integer :: nx,ny,nz
    nx = size(mesh,1); ny = size(mesh,2); nz = size(mesh,3)

    namelist /PROBINPUT/  p_infty, Rgas, gamma, mu, rho_0, p_amb, thick, minVF, rhoRatio, pRatio, &
                          p_infty_2, Rgas_2, gamma_2, mu_2, rho_0_2, plastic, explPlast, yield,   &
                          plastic2, explPlast2, yield2, interface_init, kwave, a_ratio,&
                          melt_t, melt_c, melt_t2, melt_c2, &
                          kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e,kos_sh, &
                          kos_b2,kos_t2,kos_h2,kos_g2,kos_m2,kos_q2,kos_f2,kos_alpha2,kos_beta2,kos_e2,kos_sh2, &
                          eta_det_ge,eta_det_ge_2,eta_det_gp,eta_det_gp_2,eta_det_gt,eta_det_gt_2, &
                          diff_c_ge,diff_c_ge_2,diff_c_gp,diff_c_gp_2,diff_c_gt,diff_c_gt_2, v0, v0_2, tau0, &
			tau0_2, eta0k, Nvel, etasize, ksize, delta_rho, Nrho, delta, Tp, thick

    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    ! Initialize mygfil
    call mygfil%init(                        decomp, &
                     .FALSE.,     .TRUE.,    .TRUE., &
                  "gaussian", "gaussian", "gaussian" )

    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index), &
                 rho => fields(:,:,:,rho_index), x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )        
        if (mix%ns /= 2) then
            call GracefulExit("Number of species must be 2 for this problem. Check the input file.",928)
        end if
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
        a1 = sqrt((gamma*(p1+p_infty) + 4.0d0/3.0d0*mu)/rho1)
        a2 = sqrt((gamma*(p2+p_infty) + 4.0d0/3.0d0*mu)/rho2)

        
 	        ! Set up smearing function for VF based on interface location and thickness
       ! tmp = half * ( one - erf( (x-(interface_init+eta0k/(2.0_rkind*pi*kwave)*sin(2.0_rkind*kwave*pi*y)))/(thick*dx) ) )
	!	delta_rho = Nrho * dx * 0.275d0 !converts from Nrho to approximate thickness of erf profile
	!delta_rho = Nrho*0.275d0
	eta = (x - 0.5)**2 + (y - 0.75)**2

	tmp = (half-minVF)  * ( one + tanh( (eta-(R**2))/(thick*dx) ) )
	
	!set mixture Volume fraction
	!eta = (x - 0.5)**2 + (y - 0.75)**2
	!where( (eta .le. (R)**2  )  
	mix%material(1)%VF = tmp
	mix%material(2)%VF = one - mix%material(1)%VF
	!elsewhere(eta .ge. (R)**2 )	
	!mix%material(1)%VF = 0
	!mix%material(2)%VF = one - mix%material(1)%VF
	!endwhere 
               
	!Set density profile and mass fraction based on volume fraction
	rho = rho_0*mix%material(1)%VF + rho_0_2*mix%material(2)%VF
	mix%material(1)%Ys = mix%material(1)%VF * rho_0 / rho
	mix%material(2)%Ys = one - mix%material(1)%Ys ! Enforce sum to unity

        !set velocities

	!delta = Nvel*dx
	!eta0k = etasize*delta
	!kwave = ksize*delta
        !eta = (x-interface_init)/(half*delta)     
        !eta =(x-interface_init-eta0k*sin(2*pi*y/kwave))/(0.5*delta)
       	!vc  = (rho_0*v0 + rho_0_2*v0_2)/(rho_0+rho_0_2)
        !where(eta .ge. 1.0d0)
       	!	v = v0_2
        !elsewhere((eta .le. 1d0) .and. (eta .ge. 0d0))
	!       v = (vc + -(3/2)*(vc - v0_2)*eta + -(1/2)*(v0_2-vc)*eta*eta*eta)
        !elsewhere( (eta .le. 0d0) .and. (eta .ge. -1d0) )
	!	 v = (vc + -(3/2)*(v0 - vc)*eta + -(1/2)*(vc-v0)*eta*eta*eta)
	!elsewhere(eta .le. -1.0d0)
        !        v = v0
	!endwhere
        !u   = zero
        !w   = zero

        !set mixture pressure (uniform)
        mix%material(1)%p  = p_amb
        mix%material(2)%p  = mix%material(1)%p


 

     ! Set initial values of g (inverse deformation gradient)
        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one

        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
	mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one





        !gt should be same as g
        mix%material(1)%gt11 = one;  mix%material(1)%gt12 = zero; mix%material(1)%gt13 = zero
        mix%material(1)%gt21 = zero; mix%material(1)%gt22 = one;  mix%material(1)%gt23 = zero
        mix%material(1)%gt31 = zero; mix%material(1)%gt32 = zero; mix%material(1)%gt33 = one
        
        
        mix%material(1)%gp11 = one;  mix%material(1)%gp12 = zero; mix%material(1)%gp13 = zero
        mix%material(1)%gp21 = zero; mix%material(1)%gp22 = one;  mix%material(1)%gp23 = zero
        mix%material(1)%gp31 = zero; mix%material(1)%gp32 = zero; mix%material(1)%gp33 = one
        


        mix%material(2)%gt11 = one;  mix%material(2)%gt12 = zero; mix%material(2)%gt13 = zero
        mix%material(2)%gt21 = zero; mix%material(2)%gt22 = one;  mix%material(2)%gt23 = zero
        mix%material(2)%gt31 = zero; mix%material(2)%gt32 = zero; mix%material(2)%gt33 = one
        
        mix%material(2)%gt11 = mix%material(1)%gt11
        
        mix%material(2)%gp11 = one;  mix%material(2)%gp12 = zero; mix%material(2)%gp13 = zero
        mix%material(2)%gp21 = zero; mix%material(2)%gp22 = one;  mix%material(2)%gp23 = zero
        mix%material(2)%gp31 = zero; mix%material(2)%gp32 = zero; mix%material(2)%gp33 = one
        
    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,eps,half,one,two,pi,four,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                sxx_index,syy_index,szz_index,sxy_index,sxz_index,syz_index,sos_index
    use decomp_2d,        only: decomp_info, nrank
    use DerivativesMod,   only: derivatives
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: curl
    use reductions,       only: P_SUM, P_MEAN, P_MAXVAL, P_MINVAL

    use Droplet_Shear_Incompressible_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der   
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix
    integer, dimension(2),           intent(in) :: x_bc, y_bc, z_bc
    integer                                     :: outputunit=229

    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),3) :: vort
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)  ) :: tmp
    real(rkind), dimension(decomp%ysz(1)) :: Ys1_mean, Ys2_mean
    real(rkind) :: vort_pos, vort_neg, mixwidth, Al_mass, xspike, xbubbl, xspike_proc, xbubbl_proc, ep, epp, pe, det_e, det_t, det_p, curl_e, curl_t, curl_p, rcount=0.0,eh,eel,ek,eh2,eel2,pe2
    character(len=clen) :: outputfile, str,str2
    integer :: i, j, k,nfcns

    !integer :: nx,ny,nz
    !nx = size(mesh,1); ny = size(mesh,2); nz = size(mesh,3)
    !print*,nx,ny,nz,decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 sxx  => fields(:,:,:, sxx_index), syy => fields(:,:,:,syy_index), &
                 szz  => fields(:,:,:, szz_index), sxy => fields(:,:,:,sxy_index), &
                 sxz  => fields(:,:,:, sxz_index), syz => fields(:,:,:,syz_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3),       &
                 sos  => fields(:,:,:,sos_index) )

       if (rhoRatio > 0) then
           write(str,'(I4.4,A,ES7.1E2,A,ES7.1E2)') nrank, "_", minVF, "_", rhoRatio
       else
           write(str,'(I4.4,A,ES7.1E2,A,ES7.1E2)') nrank, "_", minVF, "_", rho_0_2/rho_0
       end if
       
       if (mix%use_gTg) then
           str = trim(str)//'_gTg'
       else
           str = trim(str)//'_g'
       end if

       ! if (decomp%ysz(2) == 1) then
       !     write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Droplet_Shear_Incompressible_"//trim(str)//"_", vizcount, ".dat"

       !     open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
       !     write(outputunit,'(4ES27.16E3)') tsim, minVF, thick, rhoRatio
       !     do i=1,decomp%ysz(1)
       !         write(outputunit,'(23ES27.16E3)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
       !                                        mix%material(1)%p (i,1,1), mix%material(2)%p (i,1,1), &
       !                                        mix%material(1)%Ys(i,1,1), mix%material(2)%Ys(i,1,1), &
       !                                        mix%material(1)%VF(i,1,1), mix%material(2)%VF(i,1,1), &
       !                                        mix%material(1)%eh(i,1,1), mix%material(2)%eh(i,1,1), &
       !                                        mix%material(1)%T (i,1,1), mix%material(2)%T (i,1,1), &
       !                                        mix%material(1)%g11(i,1,1), mix%material(2)%g11(i,1,1), &
       !                                        mu(i,1,1), bulk(i,1,1), mix%material(1)%kap(i,1,1), mix%material(2)%kap(i,1,1), &
       !                                        mix%material(1)%diff(i,1,1), mix%material(2)%diff(i,1,1)
       !     end do
       !     close(outputunit)
       ! end if

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

       xspike_proc = -two
       xbubbl_proc = four
       do j = 1,decomp%ysz(2)
           do i = 1,decomp%ysz(1)
               if (mix%material(1)%Ys(i,j,1) .GE. half) xspike_proc = max(xspike_proc,x(i,j,1))
               if (mix%material(1)%Ys(i,j,1) .LE. half) xbubbl_proc = min(xbubbl_proc,x(i,j,1))
           end do
       end do
       xspike = P_MAXVAL(xspike_proc)
       xbubbl = P_MINVAL(xbubbl_proc)


       ep = P_MEAN(mix%material(1)%e_p)*six*one
       epp = P_MEAN(mix%material(1)%e_pp)*six*one
       !pe = P_MEAN(mix%material(1)%pe/rho)*six*one

       !mca add
       pe = P_MEAN(mix%material(1)%pe)*six*one
       eh = P_MEAN(mix%material(1)%eh)*six*one
       eel = P_MEAN(mix%material(1)%eel)*six*one
       ek = P_MEAN(0.5*rho*(u**2+v**2+w**2))*six*one

       pe2 = P_MEAN(mix%material(2)%pe)*six*one
       eh2 = P_MEAN(mix%material(2)%eh)*six*one
       eel2 = P_MEAN(mix%material(2)%eel)*six*one



       det_e = P_MEAN(mix%material(1)%det_e)*six*one
       det_t = P_MEAN(mix%material(1)%det_t)*six*one
       det_p = P_MEAN(mix%material(1)%det_p)*six*one
       curl_e = P_MEAN(mix%material(1)%curl_e)*six*one
       curl_t = P_MEAN(mix%material(1)%curl_t)*six*one
       curl_p = P_MEAN(mix%material(1)%curl_p)*six*one

          
       write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Droplet_Shear_Incompressible_statistics.dat"

       if (vizcount == 0) then
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
           write(outputunit,'(22A27)') 'tsim', 'mixwidth', 'vort_pos', 'vort_neg', 'Al_mass', 'xspike', 'xbubbl', 'm1_ep', 'm1_epp', 'm1_pe', 'm1_det_e', 'm1_det_t', 'm1_det_p', 'm1_curl_e', 'm1_curl_t', 'm1_curl_p', 'ek', 'm1_eh', 'm1_eel', 'm2_eh', 'm2_eel', 'm2_pe'
       else
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED', action='WRITE', status='OLD', position='APPEND')
       end if
       !write(outputunit,'(7ES27.16E3)') tsim, mixwidth, vort_pos, vort_neg, Al_mass, xspike, xbubbl
       !write(outputunit,'(10ES27.16E3)') tsim, mixwidth, vort_pos, vort_neg, Al_mass, xspike, xbubbl, ep, epp, pe
       if (nrank.eq.0) then
          !write(outputunit,'(12ES27.16E3)') real(nrank,rkind),rcount,tsim, mixwidth, vort_pos, vort_neg, Al_mass, xspike, xbubbl, ep, epp, pe
          !write(outputunit,'(10ES27.16E3)') tsim, mixwidth, vort_pos, vort_neg, Al_mass, xspike, xbubbl, ep, epp, pe
          write(outputunit,'(22F40.24)') tsim, mixwidth, vort_pos, vort_neg, Al_mass, xspike, xbubbl, ep, epp, pe, det_e, det_t, det_p, curl_e, curl_t, curl_p, ek, eh, eel, eh2, eel2, pe2
       endif
       rcount = rcount + 1.0
       close(outputunit)

       ! write tec
        ! write(outputfile,'(4A)') trim(outputdir),"/tec_MultSpecShock_"//trim(str),".dat"
        ! if(vizcount==0) then
        !   open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='replace')
        !   write(outputunit,'(350a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p", &
        !                              "sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar", &
        !                              "p-1","Ys-1","VF-1","eh-1","T-1","g11-1","g12-1","g13-1","g21-1","g22-1","g23-1","g31-1","g32-1","g33-1","Dstar-1","kap-1","rhom-1",&
        !                              "p-2","Ys-2","VF-2","eh-2","T-2","g11-2","g12-2","g13-2","g21-2","g22-2","g23-2","g31-2","g32-2","g33-2","Dstar-2","kap-2","rhom-2",&
        !                              "yield-1","mu-1","yield-2","mu-2","ep-1","ep-2","gt11-1","gt12-1","gt13-1","gt21-1","gt22-1","gt23-1","gt31-1","gt32-1","gt33-1","gt11-2","gt12-2","gt13-2","gt21-2","gt22-2","gt23-2","gt31-2","gt32-2","gt33-2"'
        !   write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
        !   write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
        !   do k=1,decomp%ysz(3)
        !    do j=1,decomp%ysz(2)
        !     do i=1,decomp%ysz(1)
        !         write(outputunit,'(56ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &                      ! continuum (9)
        !                                         sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), &      ! continuum (9)
        !                                         mix%material(1)% p(i,j,k),  mix%material(1)% Ys(i,j,k), mix%material(1)% VF(i,j,k), mix%material(1)% eh(i,j,k), &  ! material 1 (14)
        !                                         mix%material(1)% T(i,j,k),  mix%material(1)%g11(i,j,k), mix%material(1)%g12(i,j,k), mix%material(1)%g13(i,j,k), &  ! material 1 
        !                                         mix%material(1)%g21(i,j,k), mix%material(1)%g22(i,j,k), mix%material(1)%g23(i,j,k), mix%material(1)%g31(i,j,k), &  ! material 1 
        !                                         mix%material(1)%g32(i,j,k), mix%material(1)%g33(i,j,k), mix%material(1)%diff(i,j,k), mix%material(1)%kap(i,j,k), mix%material(1)%rhom(i,j,k), & ! material 1 
        !                                         mix%material(2)% p(i,j,k),  mix%material(2)% Ys(i,j,k), mix%material(2)% VF(i,j,k), mix%material(2)% eh(i,j,k), &  ! material 2 (14)
        !                                         mix%material(2)% T(i,j,k),  mix%material(2)%g11(i,j,k), mix%material(2)%g12(i,j,k), mix%material(2)%g13(i,j,k), &  ! material 2
        !                                         mix%material(2)%g21(i,j,k), mix%material(2)%g22(i,j,k), mix%material(2)%g23(i,j,k), mix%material(2)%g31(i,j,k), &  ! material 2
        !                                         mix%material(2)%g32(i,j,k), mix%material(2)%g33(i,j,k), mix%material(2)%diff(i,j,k), mix%material(2)%kap(i,j,k), mix%material(2)%rhom(i,j,k), &    ! material 2
        !                                         mix%material(1)%elastic%yield(i,j,k),mix%material(1)%elastic%mu(i,j,k), mix%material(2)%elastic%yield(i,j,k),mix%material(2)%elastic%mu(i,j,k)! , mix%material(2)%e_p(i,j,k),mix%material(2)%e_p(i,j,k), &
        !                                         ! mix%material(1)%gt11(i,j,k), mix%material(1)%gt12(i,j,k), mix%material(1)%gt13(i,j,k), mix%material(1)%gt21(i,j,k), mix%material(1)%gt22(i,j,k), mix%material(1)%gt23(i,j,k), mix%material(1)%gt31(i,j,k), mix%material(1)%gt32(i,j,k), mix%material(1)%gt33(i,j,k), &
        !                                         ! mix%material(2)%gt11(i,j,k), mix%material(2)%gt12(i,j,k), mix%material(2)%gt13(i,j,k), mix%material(2)%gt21(i,j,k), mix%material(2)%gt22(i,j,k), mix%material(2)%gt23(i,j,k), mix%material(2)%gt31(i,j,k), mix%material(2)%gt32(i,j,k), mix%material(2)%gt33(i,j,k)
        !     end do
        !    end do
        !   end do
        !   close(outputunit)
        ! else
        !   open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='old', action='write', position='append')
        !   write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
        !   write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
        !   write(outputunit,'(a)') ' VARSHARELIST=([1, 2, 3]=1)'
        !   do k=1,decomp%ysz(3)
        !    do j=1,decomp%ysz(2)
        !     do i=1,decomp%ysz(1)
        !         write(outputunit,'(53ES26.16)') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &                                                    ! continuum (6)
        !                                         sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), &      ! continuum (9)
        !                                         mix%material(1)% p(i,j,k),  mix%material(1)% Ys(i,j,k), mix%material(1)% VF(i,j,k), mix%material(1)% eh(i,j,k), &  ! material 1 (14)
        !                                         mix%material(1)% T(i,j,k),  mix%material(1)%g11(i,j,k), mix%material(1)%g12(i,j,k), mix%material(1)%g13(i,j,k), &  ! material 1 
        !                                         mix%material(1)%g21(i,j,k), mix%material(1)%g22(i,j,k), mix%material(1)%g23(i,j,k), mix%material(1)%g31(i,j,k), &  ! material 1 
        !                                         mix%material(1)%g32(i,j,k), mix%material(1)%g33(i,j,k), mix%material(1)%diff(i,j,k), mix%material(1)%kap(i,j,k), mix%material(1)%rhom(i,j,k),&  ! material 1 
        !                                         mix%material(2)% p(i,j,k),  mix%material(2)% Ys(i,j,k), mix%material(2)% VF(i,j,k), mix%material(2)% eh(i,j,k), &  ! material 2 (14)
        !                                         mix%material(2)% T(i,j,k),  mix%material(2)%g11(i,j,k), mix%material(2)%g12(i,j,k), mix%material(2)%g13(i,j,k), &  ! material 2
        !                                         mix%material(2)%g21(i,j,k), mix%material(2)%g22(i,j,k), mix%material(2)%g23(i,j,k), mix%material(2)%g31(i,j,k), &  ! material 2
        !                                         mix%material(2)%g32(i,j,k), mix%material(2)%g33(i,j,k), mix%material(2)%diff(i,j,k), mix%material(2)%kap(i,j,k), mix%material(2)%rhom(i,j,k),&    ! material 2
        !                                         mix%material(1)%elastic%yield(i,j,k),mix%material(1)%elastic%mu(i,j,k), mix%material(2)%elastic%yield(i,j,k),mix%material(2)%elastic%mu(i,j,k)! , mix%material(2)%e_p(i,j,k),mix%material(2)%e_p(i,j,k), &
        !                                         ! mix%material(1)%gt11(i,j,k), mix%material(1)%gt12(i,j,k), mix%material(1)%gt13(i,j,k), mix%material(1)%gt21(i,j,k), mix%material(1)%gt22(i,j,k), mix%material(1)%gt23(i,j,k), mix%material(1)%gt31(i,j,k), mix%material(1)%gt32(i,j,k), mix%material(1)%gt33(i,j,k), &
        !                                         ! mix%material(2)%gt11(i,j,k), mix%material(2)%gt12(i,j,k), mix%material(2)%gt13(i,j,k), mix%material(2)%gt21(i,j,k), mix%material(2)%gt22(i,j,k), mix%material(2)%gt23(i,j,k), mix%material(2)%gt31(i,j,k), mix%material(2)%gt32(i,j,k), mix%material(2)%gt33(i,j,k)
        !     end do
        !    end do
        !   end do
        !   close(outputunit)
        ! endif

        !write plot3D
       !print*, nrank, decomp%xsz(1),decomp%xsz(2),decomp%xsz(3), decomp%ysz(1),decomp%ysz(2),decomp%ysz(3), decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)
       !print*,nrank,decomp%yst(1), decomp%yst(2), decomp%yst(3), decomp%yen(1), decomp%yen(2), decomp%yen(3)
       !print*, nrank,decomp%xsz(1),decomp%ysz(2),decomp%zsz(3),decomp%yst(1),decomp%yen(1),decomp%yst(2),decomp%yen(2),decomp%yst(3),decomp%yen(3) !global nx,ny,nz -- ist,ien,jst,jen,kst,ken


        !write(str,'(I4.4)') nrank
        write(str,*) nrank
        str = adjustl(str)

        if(vizcount==0) then
           write(outputfile,'(4A)') trim(outputdir),"/P3D_MSS_",trim(str),".grd"
           open(unit=outputunit, file=trim(outputfile), form='UNFORMATTED', status='replace')
           write(outputunit) decomp%xsz(1),decomp%ysz(2),decomp%zsz(3),decomp%yst(1),decomp%yen(1),decomp%yst(2),decomp%yen(2),decomp%yst(3),decomp%yen(3) !global nx,ny,nz -- ist,ien,jst,jen,kst,ken
           write(outputunit) decomp%ysz(1),decomp%ysz(2),decomp%ysz(3) !local nx,ny,nz
           write(outputunit) x,y,z
           close(outputunit)
        endif

        if(.false.) then !write P3D
           nfcns=15 + 50*2
           write(str2,*) vizcount
           str2 = adjustl(str2)
           !write(outputfile,'(4A,I4.4,A)') trim(outputdir),"/P3D_MSS_",trim(str),"_",vizcount,".fcn"
           write(outputfile,'(6A)') trim(outputdir),"/P3D_MSS_",trim(str),"_",trim(str2),".fcn"
           open(unit=outputunit, file=trim(outputfile), form='UNFORMATTED', status='replace')
           write(outputunit) decomp%xsz(1),decomp%ysz(2),decomp%zsz(3),decomp%yst(1),decomp%yen(1),decomp%yst(2),decomp%yen(2),decomp%yst(3),decomp%yen(3) !global nx,ny,nz -- ist,ien,jst,jen,kst,ken
           write(outputunit) decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),nfcns !local nx,ny,nz, nfunctions
           write(outputunit) rho,u,v,w,e,p,sxx,sxy,sxz,syy,syz,szz,mu,bulk,kap, & !15
                !material 1
                mix%material(1)%p, mix%material(1)%T, mix%material(1)%Ys, mix%material(1)%VF, mix%material(1)%eh, mix%material(1)%eel, &
                mix%material(1)%elastic%mu, mix%material(1)%elastic%yield, mix%material(1)%e_p, mix%material(1)%e_pp, mix%material(1)%pe, & !11
                mix%material(1)%g11, mix%material(1)%g12, mix%material(1)%g13, & !27
                mix%material(1)%g21, mix%material(1)%g22, mix%material(1)%g23, &
                mix%material(1)%g31, mix%material(1)%g32, mix%material(1)%g33, &
                mix%material(1)%gt11, mix%material(1)%gt12, mix%material(1)%gt13, &
                mix%material(1)%gt21, mix%material(1)%gt22, mix%material(1)%gt23, &
                mix%material(1)%gt31, mix%material(1)%gt32, mix%material(1)%gt33, &
                mix%material(1)%gp11, mix%material(1)%gp12, mix%material(1)%gp13, &
                mix%material(1)%gp21, mix%material(1)%gp22, mix%material(1)%gp23, &
                mix%material(1)%gp31, mix%material(1)%gp32, mix%material(1)%gp33, &
                mix%material(1)%kap, mix%material(1)%diff, mix%material(1)%diff_g, mix%material(1)%diff_gt, mix%material(1)%diff_gp, mix%material(1)%diff_pe,& !12
                mix%material(1)%det_e, mix%material(1)%det_p, mix%material(1)%det_t, mix%material(1)%curl_e, mix%material(1)%curl_p, mix%material(1)%curl_t, &
                !material 2
                mix%material(2)%p, mix%material(2)%T, mix%material(2)%Ys, mix%material(2)%VF, mix%material(2)%eh, mix%material(2)%eel, &
                mix%material(2)%elastic%mu, mix%material(2)%elastic%yield, mix%material(2)%e_p, mix%material(2)%e_pp, mix%material(2)%pe, & !11
                mix%material(2)%g11, mix%material(2)%g12, mix%material(2)%g13, & !27
                mix%material(2)%g21, mix%material(2)%g22, mix%material(2)%g23, &
                mix%material(2)%g31, mix%material(2)%g32, mix%material(2)%g33, &
                mix%material(2)%gt11, mix%material(2)%gt12, mix%material(2)%gt13, &
                mix%material(2)%gt21, mix%material(2)%gt22, mix%material(2)%gt23, &
                mix%material(2)%gt31, mix%material(2)%gt32, mix%material(2)%gt33, &
                mix%material(2)%gp11, mix%material(2)%gp12, mix%material(2)%gp13, &
                mix%material(2)%gp21, mix%material(2)%gp22, mix%material(2)%gp23, &
                mix%material(2)%gp31, mix%material(2)%gp32, mix%material(2)%gp33, &
                mix%material(2)%kap, mix%material(2)%diff, mix%material(2)%diff_g, mix%material(2)%diff_gt, mix%material(2)%diff_gp, mix%material(2)%diff_pe,& !12
                mix%material(2)%det_e, mix%material(2)%det_p, mix%material(2)%det_t, mix%material(2)%curl_e, mix%material(2)%curl_p, mix%material(2)%curl_t
           close(outputunit)
        endif

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: filter3D

    use Droplet_Shear_Incompressible_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc
    
    integer :: nx, ny, i, j
    real(rkind) :: dx, dy, xspng, tspng, xspng2, yspng, yspng2
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum
    
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
	
	if(decomp%yst(1)==1) then
          if(x_bc(1)==0) then
        !      rho( 1,:,:) = rhoL
        !      u  ( 1,:,:) = zero
        !      v  ( 1,:,:) = vL
        !      w  ( 1,:,:) = zero
        !      mix%material(1)%p(1,:,:) = p_amb
        !      mix%material(2)%p(1,:,:) = p_amb
        !      
              mix%material(2)%VF ( 1,:,:) = minVF
              mix%material(1)%VF ( 1,:,:) = one - minVF
              mix%material(2)%Ys ( 1,:,:) = minVF
              mix%material(1)%Ys ( 1,:,:) = one - minVF
          end if
        endif

        if(decomp%yen(1)==decomp%xsz(1)) then
          if(x_bc(2)==0) then
        !      rho( nx,:,:) = rhoR
        !      u  ( nx,:,:) = zero
        !      v  ( nx,:,:) = vR
        !      w  ( nx,:,:) = zero
        !      mix%material(1)%p(nx,:,:) = p_amb
        !      mix%material(2)%p(nx,:,:) = p_amb
        !      
              mix%material(2)%VF ( nx,:,:) = minVF
              mix%material(1)%VF ( nx,:,:) = one - minVF
              mix%material(2)%Ys ( nx,:,:) = minVF
              mix%material(1)%Ys ( nx,:,:) = one - minVF
          end if
        endif
        if(decomp%yst(2)==1) then
          if(y_bc(1)==0) then
        !      rho( 1,:,:) = rhoL
        !      u  ( 1,:,:) = zero
        !      v  ( 1,:,:) = vL
        !      w  ( 1,:,:) = zero
        !      mix%material(1)%p(1,:,:) = p_amb
        !      mix%material(2)%p(1,:,:) = p_amb
        !
              mix%material(2)%VF ( :,1,:) = minVF
              mix%material(1)%VF ( :,1,:) = one - minVF
              mix%material(2)%Ys ( :,1,:) = minVF
              mix%material(1)%Ys ( :,1,:) = one - minVF
          end if
        endif

      if(decomp%yen(2)==decomp%ysz(2)) then
          if(y_bc(2)==0) then
        !      rho( nx,:,:) = rhoR
        !      u  ( nx,:,:) = zero
        !      v  ( nx,:,:) = vR
        !      w  ( nx,:,:) = zero
        !      mix%material(1)%p(nx,:,:) = p_amb
        !      mix%material(2)%p(nx,:,:) = p_amb
        !
              mix%material(2)%VF ( :,ny,:) = minVF
              mix%material(1)%VF ( :,ny,:) = one - minVF
              mix%material(2)%Ys ( :,ny,:) = minVF
              mix%material(1)%Ys ( :,ny,:) = one - minVF
          end if
        endif
        !if(decomp%yst(1)==1) then
         !if(x_bc(1)==0) then
          !    ! rho( 1,:,:) = rhoL
           !   ! u  ( 1,:,:) = (u2-u1)
           !   v  ( 1,:,:) = zero
           !   w  ( 1,:,:) = zero
           !   do i=1,5
           !       mix%material(1)%p( i,:,:) = mix%material(1)%p(6,:,:)
           !       mix%material(2)%p( i,:,:) = mix%material(2)%p(6,:,:)
           !   end do
              
              ! mix%material(1)%g11( 1,:,:) = rho2/rho_0; mix%material(1)%g12( 1,:,:) = zero; mix%material(1)%g13( 1,:,:) = zero
              ! mix%material(1)%g21( 1,:,:) = zero; mix%material(1)%g22( 1,:,:) = one;  mix%material(1)%g23( 1,:,:) = zero
              ! mix%material(1)%g31( 1,:,:) = zero; mix%material(1)%g32( 1,:,:) = zero; mix%material(1)%g33( 1,:,:) = one
  
              ! mix%material(2)%g11( 1,:,:) = rho2/rho_0;  mix%material(2)%g12( 1,:,:) = zero; mix%material(2)%g13( 1,:,:) = zero
              ! mix%material(2)%g21( 1,:,:) = zero; mix%material(2)%g22( 1,:,:) = one;  mix%material(2)%g23( 1,:,:) = zero
              ! mix%material(2)%g31( 1,:,:) = zero; mix%material(2)%g32( 1,:,:) = zero; mix%material(2)%g33( 1,:,:) = one
              
              ! mix%material(1)%Ys ( 1,:,:) = YsL
              ! mix%material(2)%Ys ( 1,:,:) = one - YsL
  
            !  mix%material(1)%VF ( 1,:,:) = VFL
            !  mix%material(2)%VF ( 1,:,:) = one - VFL
         ! end if
        !endif

        !xspng = 0
	!xspng2 = 1
	!yspng = 0
	!yspng2 = 1
        !tspng = 0.1_rkind
       !dx = x(2,1,1) - x(1,1,1)
       !dum = half*(one - tanh( (x-xspng)/(tspng) ))+half*( tanh( (x-xspng2)/(tspng) ))+half*(one - tanh( (y-yspng)/(tspng) ))+half*(tanh( (y-yspng2)/(tspng) ))

        !do i=1,4
        !    tmp = u
        !    call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
        !    u = u + dum*(tmp - u)

        !    tmp = v
        !    call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
        !    v = v + dum*(tmp - v)

        !   tmp = w
        !   call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
        !   w = w + dum*(tmp - w)

        !    tmp = e
        !    call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
         !   e = e + dum*(tmp - e)

          !  tmp = rho
          !  call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
          !  rho = rho + dum*(tmp - rho)

          !  tmp = mix%material(1)%p
          !  call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
          !  mix%material(1)%p = mix%material(1)%p + dum*(tmp - mix%material(1)%p)

          !  tmp = mix%material(2)%p
          !  call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
          !  mix%material(2)%p = mix%material(2)%p + dum*(tmp - mix%material(2)%p)

          !  tmp = mix%material(1)%pe
          !  call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
          !  mix%material(1)%pe = mix%material(1)%pe + dum*(tmp - mix%material(1)%pe)

          !  tmp = mix%material(2)%pe
          !  call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
          !  mix%material(2)%pe = mix%material(2)%pe + dum*(tmp - mix%material(2)%pe)

          !  do j = 1,9
          !      tmp = mix%material(1)%g(:,:,:,j)
          !      call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
          !      mix%material(1)%g(:,:,:,j) = mix%material(1)%g(:,:,:,j) + dum*(tmp - mix%material(1)%g(:,:,:,j))

          !      tmp = mix%material(2)%g(:,:,:,j)
          !      call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
          !      mix%material(2)%g(:,:,:,j) = mix%material(2)%g(:,:,:,j) + dum*(tmp - mix%material(2)%g(:,:,:,j))

             !  tmp = mix%material(1)%g_t(:,:,:,j)
             !  call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
             !  mix%material(1)%g_t(:,:,:,j) = mix%material(1)%g_t(:,:,:,j) + dum*(tmp - mix%material(1)%g_t(:,:,:,j))

             !   tmp = mix%material(2)%g_t(:,:,:,j)
             !   call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
             !   mix%material(2)%g_t(:,:,:,j) = mix%material(2)%g_t(:,:,:,j) + dum*(tmp - mix%material(2)%g_t(:,:,:,j))

             !   tmp = mix%material(1)%g_p(:,:,:,j)
             !   call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
             !   mix%material(1)%g_p(:,:,:,j) = mix%material(1)%g_p(:,:,:,j) + dum*(tmp - mix%material(1)%g_p(:,:,:,j))

              !  tmp = mix%material(2)%g_p(:,:,:,j)
              !  call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
              !  mix%material(2)%g_p(:,:,:,j) = mix%material(2)%g_p(:,:,:,j) + dum*(tmp - mix%material(2)%g_p(:,:,:,j))
            !end do

            !mca add for stability

            !tmp = T
            !call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            !T = T + dum*(tmp - T)

            !tmp = mix%material(1)%T
            !call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            !mix%material(1)%T = mix%material(1)%T + dum*(tmp - mix%material(1)%T)

            !tmp = mix%material(2)%T
            !call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            !mix%material(2)%T = mix%material(2)%T + dum*(tmp - mix%material(2)%T)

            !tmp = mix%material(1)%Ys
            !call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            !mix%material(1)%Ys = mix%material(1)%Ys + dum*(tmp - mix%material(1)%Ys)

            !tmp = mix%material(2)%Ys
        !    call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
        !    mix%material(2)%Ys = mix%material(2)%Ys + dum*(tmp - mix%material(2)%Ys)

       ! end do

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

    use Droplet_Shear_Incompressible_data

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

         ! set u and v
       ! IF (tsim .LT. 4) THEN
	u = -  ((sin(pi*x))**two)*(sin(two*pi*y))*(cos(pi*tsim/Tp))
	v =  ((sin(pi*y))**two)*(sin(two*pi*x))*(cos(pi*tsim/Tp))
	w = 0
	p = p_amb
!	T = 6010
!	e = 7768
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

    use Droplet_Shear_Incompressible_data

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

    use Droplet_Shear_Incompressible_data

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

    use Droplet_Shear_Incompressible_data

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

    use Droplet_Shear_Incompressible_data

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

    use Droplet_Shear_Incompressible_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
