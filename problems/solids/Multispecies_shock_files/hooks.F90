module Multispecies_shock_data
    use kind_parameters,  only: rkind
    use constants,        only: one,two,eight,three,six
    use FiltersMod,       only: filters
    implicit none

    real(rkind) :: p_infty = one, Rgas = one, gamma = 1.4_rkind, mu = 10._rkind, rho_0 = one, p_amb = 0.1_rkind
    real(rkind) :: p_infty_2 = one, Rgas_2 = one, gamma_2 = 1.4_rkind, mu_2 = 10._rkind, rho_0_2 = one
    real(rkind) :: minVF = 0.2_rkind, thick = one
    real(rkind) :: rhoRatio = one, pRatio = two
    logical     :: sharp = .FALSE.
    real(rkind) :: p1,p2,rho1,rho2,u1,u2,g11_1,g11_2,grho1,grho2,a1,a2
    real(rkind) :: rho1_2,rho2_2,u1_2,u2_2,g11_1_2,g11_2_2,grho1_2,grho2_2,a1_2,a2_2
    real(rkind) :: rhoL, rhoR, YsL, YsR, VFL, VFR
    real(rkind) :: yield = one, yield2 = one, eta0k = 0.4_rkind
    logical     :: explPlast = .FALSE., explPlast2 = .FALSE.
    logical     :: plastic = .FALSE., plastic2 = .FALSE.
    real(rkind) :: Ly = one, Lx = six, interface_init = 0.75_rkind, shock_init = 0.6_rkind, kwave = 4.0_rkind
    logical     :: sliding = .false.

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

    use Multispecies_shock_data

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
                    x(i,j,k) = real( ix1     + i - 1, rkind ) * dx - two  ! x \in (-2,4]
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven,pi
    use SolidGrid,        only: u_index,v_index,w_index
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    use SolidMixtureMod,  only: solid_mixture
    
    use Multispecies_shock_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(inout) :: tstop, dt, tviz
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum
    real(rkind), dimension(8) :: fparams
    integer, dimension(2) :: iparams

    namelist /PROBINPUT/  p_infty, Rgas, gamma, mu, rho_0, p_amb, thick, minVF, rhoRatio, pRatio, &
                          p_infty_2, Rgas_2, gamma_2, mu_2, rho_0_2, plastic, explPlast, yield,   &
                          plastic2, explPlast2, yield2, interface_init, kwave, sliding
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    ! Initialize mygfil
    call mygfil%init(                        decomp, &
                     .FALSE.,     .TRUE.,    .TRUE., &
                  "gaussian", "gaussian", "gaussian" )

    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        if (mix%ns /= 2) then
            call GracefulExit("Number of species must be 2 for this problem. Check the input file.",928)
        end if

        if(rhoRatio > 0) then
          ! if rhoRatio is positive, only rho_0 is different. Rgas is set such
          ! that Temperature equilibrium condition is satisfied
          gamma_2 = gamma; Rgas_2 = Rgas/rhoRatio; p_infty_2 = p_infty; 
          rho_0_2 = rho_0*rhoRatio; mu_2 = mu
        else
          ! if rhoRatio is negative, all quantities except Rgas need to be
          ! specified in input file. Rgas is then set such
          ! that Temperature equilibrium condition is satisfied
          Rgas_2 = Rgas * (p_amb+p_infty_2)/(p_amb+p_infty)*rho_0/rho_0_2
        endif

        ! write material properties
        if (nrank == 0) then
            print *, '---Material 1---'
            write(*,'(3(a,e12.5))') 'rho_0 = ', rho_0, ', gam  = ', gamma, ', p_infty = ', p_infty
            write(*,'(3(a,e12.5))') 'mu    = ', mu,    ', Rgas = ', Rgas
            print *, '---Material 2---'
            write(*,'(3(a,e12.5))') 'rho_0 = ', rho_0_2, ', gam  = ', gamma_2, ', p_infty = ', p_infty_2
            write(*,'(3(a,e12.5))') 'mu    = ', mu_2,    ', Rgas = ', Rgas_2
        end if

        ! Set materials
        call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield,1.0D-10))
        call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,1.0D-10))

        ! set logicals for plasticity
        mix%material(1)%plast = plastic; mix%material(1)%explPlast = explPlast
        mix%material(2)%plast = plastic2; mix%material(2)%explPlast = explPlast2

        ! Set logicals for sliding
        mix%material(1)%sliding = sliding
        mix%material(2)%sliding = sliding

        ! determine jump conditions for material 1
        p1 = p_amb
        rho1 = rho_0
        p2 = p1 * pRatio
        u1 = zero
        fparams(1) = rho1; fparams(2) = u1; fparams(3) = p1
        ! fparams(4) = rho_0; fparams(5) = gamma; fparams(6) = p_infty; fparams(7) = zero;
        fparams(4) = rho1; fparams(5) = gamma; fparams(6) = p_infty; fparams(7) = mu;
        fparams(8) = p2
        rho2 = rho1*min(one + p1/p_infty, one) ! Init guess
        call rootfind_nr_1d(rho2,fparams,iparams)
        ! print *, 'After root finding: ', rho2
        g11_1 = fparams(1)/fparams(4);   grho1 = g11_1**real(11.D0/3.D0,rkind) - g11_1**(-third) - g11_1**(seven*third) + g11_1**third
        g11_2 = rho2/fparams(4);         grho2 = g11_2**real(11.D0/3.D0,rkind) - g11_2**(-third) - g11_2**(seven*third) + g11_2**third
        ! u2 = -sqrt(rho1/rho2/(rho1-rho2)*(p1-p2+twothird*zero*(grho1-grho2)))
        u2 = -sqrt(rho1/rho2/(rho1-rho2)*(p1-p2+twothird*mu*(grho1-grho2)))
        u1 = rho2*u2/rho1

        ! speed of sound
        a1 = sqrt((gamma*(p1+p_infty) + 4.0d0/3.0d0*mu)/rho1)
        a2 = sqrt((gamma*(p2+p_infty) + 4.0d0/3.0d0*mu)/rho2)

        if (nrank == 0) then
            print*, '----Shock Initialization-----'
            print*, "Mass flux: ", rho1*u1, rho2*u2
            print*, "Momentum flux: ", rho1*u1*u1+p1, rho2*u2*u2+p2
            print*, "g flux: ", g11_1*u1, g11_2*u2
            
            print*, "rho1, rho2 = ", rho1, rho2
            print*, "u1, u2 = ", u1, u2
            print*, "p1, p2 = ", p1, p2
            print*, "a1, a2 = ", a1, a2
            print*, "M1, M2 = ", u1/a1, u2/a2
            print*, "p_infty = ", p_infty
            print*, "rhoRatio = ", rhoRatio
        end if

        ! determine jump conditions for material 2
        p1 = p_amb
        rho1_2 = rho_0_2
        p2 = p1 * pRatio
        u1_2 = zero
        fparams(1) = rho1_2; fparams(2) = u1_2; fparams(3) = p1
        ! fparams(4) = rho_0; fparams(5) = gamma; fparams(6) = p_infty; fparams(7) = zero;
        fparams(4) = rho1_2; fparams(5) = gamma_2; fparams(6) = p_infty_2; fparams(7) = mu_2;
        fparams(8) = p2
        rho2_2 = rho1_2*min(one + p1/p_infty_2, one) ! Init guess
        call rootfind_nr_1d(rho2_2,fparams,iparams)
        ! print *, 'After root finding: ', rho2_2
        g11_1_2 = fparams(1)/fparams(4);   grho1_2 = g11_1_2**real(11.D0/3.D0,rkind) - g11_1_2**(-third) - g11_1_2**(seven*third) + g11_1_2**third
        g11_2_2 = rho2_2/fparams(4);         grho2_2 = g11_2_2**real(11.D0/3.D0,rkind) - g11_2_2**(-third) - g11_2_2**(seven*third) + g11_2_2**third
        u2_2 = -sqrt(rho1_2/rho2_2/(rho1_2-rho2_2)*(p1-p2+twothird*mu_2*(grho1_2-grho2_2)))
        u1_2 = rho2_2*u2_2/rho1_2

        ! speed of sound
        a1_2 = sqrt((gamma_2*(p1+p_infty_2) + 4.0d0/3.0d0*mu_2)/rho1_2)
        a2_2 = sqrt((gamma_2*(p2+p_infty_2) + 4.0d0/3.0d0*mu_2)/rho2_2)

        if (nrank == 0) then
            print*, '----Shock Initialization-----'
            print*, "Mass flux: ", rho1_2*u1_2, rho2_2*u2_2
            print*, "Momentum flux: ", rho1_2*u1_2*u1_2+p1, rho2_2*u2_2*u2_2+p2
            print*, "g flux: ", g11_1_2*u1_2, g11_2_2*u2_2
            
            print*, "rho1, rho2 = ", rho1_2, rho2_2
            print*, "u1, u2 = ", u1_2, u2_2
            print*, "p1, p2 = ", p1, p2
            print*, "a1, a2 = ", a1_2, a2_2
            print*, "M1, M2 = ", u1_2/a1_2, u2_2/a2_2
            print*, "p_infty = ", p_infty_2
            print*, "rhoRatio = ", rhoRatio
        end if

        ! Get mixture momentum (put in u1 and u2)
        u1 = (one-minVF)*rho1*u1 + minVF*rho1_2*u1_2
        u2 = (one-minVF)*rho2*u2 + minVF*rho2_2*u2_2
        
        ! Get mixture density
        rho1 = (one-minVF)*rho1 + minVF*rho1_2
        rho2 = (one-minVF)*rho2 + minVF*rho2_2

        ! Get mixture velocity
        u1 = u1 / rho1
        u2 = u2 / rho2

        shock_init = interface_init - 1.0_rkind  ! (10*thick) grid points away from the interface
        dum = half * ( one - erf( (x-shock_init)/(two*dx) ) )

        u   = (u2-u1)*dum
        v   = zero
        w   = zero

        tmp = half * ( one - erf( (x-(interface_init+eta0k/(2.0_rkind*pi*kwave)*sin(2.0_rkind*kwave*pi*y)))/(thick*dx) ) )

        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one
        
        mix%material(1)%g11 = (rho2*dum + rho1*(one-dum))/rho_0
        if (mix%use_gTg) then
            mix%material(1)%g11 = mix%material(1)%g11**2
        end if

        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

        mix%material(2)%g11 = mix%material(1)%g11

        mix%material(1)%p  = p2*dum + p1*(one-dum)
        mix%material(2)%p  = mix%material(1)%p

        mix%material(1)%VF = minVF + (one-two*minVF)*tmp
        mix%material(2)%VF = one - mix%material(1)%VF

        if (mix%use_gTg) then
            tmp = rho_0*mix%material(1)%VF*sqrt(mix%material(1)%g11) + rho_0_2*sqrt(mix%material(2)%g11)*(one-mix%material(1)%VF) ! Mixture density
        else
            tmp = rho_0*mix%material(1)%VF*mix%material(1)%g11 + rho_0_2*mix%material(2)%g11*(one-mix%material(1)%VF) ! Mixture density
        end if
        mix%material(1)%Ys = mix%material(1)%VF * rho_0 / tmp
        mix%material(2)%Ys = one - mix%material(1)%Ys ! Enforce sum to unity

        rhoL = tmp(1,1,1)
        rhoR = tmp(decomp%ysz(1),1,1)
        YsL  = mix%material(1)%Ys(1,1,1)
        YsR  = mix%material(1)%Ys(decomp%ysz(1),1,1)
        VFL  = mix%material(1)%VF(1,1,1)
        VFR  = mix%material(1)%VF(decomp%ysz(1),1,1)

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

    use Multispecies_shock_data

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
    real(rkind) :: vort_pos, vort_neg, mixwidth, Al_mass, xspike, xbubbl, xspike_proc, xbubbl_proc
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

       if (decomp%ysz(2) == 1) then
           write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Multispecies_shock_"//trim(str)//"_", vizcount, ".dat"

           open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
           write(outputunit,'(4ES27.16E3)') tsim, minVF, thick, rhoRatio
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
           
       write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Multispecies_shock_statistics.dat"

       if (vizcount == 0) then
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
           write(outputunit,'(7A27)') 'tsim', 'mixwidth', 'vort_pos', 'vort_neg', 'Al_mass', 'xspike', 'xbubbl'
       else
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED', action='WRITE', status='OLD', position='APPEND')
       end if
       write(outputunit,'(7ES27.16E3)') tsim, mixwidth, vort_pos, vort_neg, Al_mass, xspike, xbubbl
       close(outputunit)

        ! write(outputfile,'(4A)') trim(outputdir),"/tec_MultSpecShock_"//trim(str),".dat"
        ! if(vizcount==0) then
        !   open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='replace')
        !   write(outputunit,'(350a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p", &
        !                              "sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar", &
        !                              "p-1","Ys-1","VF-1","eh-1","T-1","g11-1","g12-1","g13-1","g21-1","g22-1","g23-1","g31-1","g32-1","g33-1","Dstar-1","kap-1",&
        !                              "p-2","Ys-2","VF-2","eh-2","T-2","g11-2","g12-2","g13-2","g21-2","g22-2","g23-2","g31-2","g32-2","g33-2","Dstar-2","kap-2"'
        !   write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
        !   write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
        !   do k=1,decomp%ysz(3)
        !    do j=1,decomp%ysz(2)
        !     do i=1,decomp%ysz(1)
        !         write(outputunit,'(50ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &                      ! continuum (9)
        !                                         sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), &      ! continuum (9)
        !                                         mix%material(1)% p(i,j,k),  mix%material(1)% Ys(i,j,k), mix%material(1)% VF(i,j,k), mix%material(1)% eh(i,j,k), &  ! material 1 (14)
        !                                         mix%material(1)% T(i,j,k),  mix%material(1)%g11(i,j,k), mix%material(1)%g12(i,j,k), mix%material(1)%g13(i,j,k), &  ! material 1 
        !                                         mix%material(1)%g21(i,j,k), mix%material(1)%g22(i,j,k), mix%material(1)%g23(i,j,k), mix%material(1)%g31(i,j,k), &  ! material 1 
        !                                         mix%material(1)%g32(i,j,k), mix%material(1)%g33(i,j,k), mix%material(1)%diff(i,j,k), mix%material(1)%kap(i,j,k),&  ! material 1 
        !                                         mix%material(2)% p(i,j,k),  mix%material(2)% Ys(i,j,k), mix%material(2)% VF(i,j,k), mix%material(2)% eh(i,j,k), &  ! material 2 (14)
        !                                         mix%material(2)% T(i,j,k),  mix%material(2)%g11(i,j,k), mix%material(2)%g12(i,j,k), mix%material(2)%g13(i,j,k), &  ! material 2
        !                                         mix%material(2)%g21(i,j,k), mix%material(2)%g22(i,j,k), mix%material(2)%g23(i,j,k), mix%material(2)%g31(i,j,k), &  ! material 2
        !                                         mix%material(2)%g32(i,j,k), mix%material(2)%g33(i,j,k), mix%material(2)%diff(i,j,k), mix%material(2)%kap(i,j,k)    ! material 2
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
        !         write(outputunit,'(47ES26.16)') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &                                                    ! continuum (6)
        !                                         sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), &      ! continuum (9)
        !                                         mix%material(1)% p(i,j,k),  mix%material(1)% Ys(i,j,k), mix%material(1)% VF(i,j,k), mix%material(1)% eh(i,j,k), &  ! material 1 (14)
        !                                         mix%material(1)% T(i,j,k),  mix%material(1)%g11(i,j,k), mix%material(1)%g12(i,j,k), mix%material(1)%g13(i,j,k), &  ! material 1 
        !                                         mix%material(1)%g21(i,j,k), mix%material(1)%g22(i,j,k), mix%material(1)%g23(i,j,k), mix%material(1)%g31(i,j,k), &  ! material 1 
        !                                         mix%material(1)%g32(i,j,k), mix%material(1)%g33(i,j,k), mix%material(1)%diff(i,j,k), mix%material(1)%kap(i,j,k),&  ! material 1 
        !                                         mix%material(2)% p(i,j,k),  mix%material(2)% Ys(i,j,k), mix%material(2)% VF(i,j,k), mix%material(2)% eh(i,j,k), &  ! material 2 (14)
        !                                         mix%material(2)% T(i,j,k),  mix%material(2)%g11(i,j,k), mix%material(2)%g12(i,j,k), mix%material(2)%g13(i,j,k), &  ! material 2
        !                                         mix%material(2)%g21(i,j,k), mix%material(2)%g22(i,j,k), mix%material(2)%g23(i,j,k), mix%material(2)%g31(i,j,k), &  ! material 2
        !                                         mix%material(2)%g32(i,j,k), mix%material(2)%g33(i,j,k), mix%material(2)%diff(i,j,k), mix%material(2)%kap(i,j,k)    ! material 2
        !     end do
        !    end do
        !   end do
        !   close(outputunit)
        ! endif


    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: filter3D

    use Multispecies_shock_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc
    
    integer :: nx, i, j
    real(rkind) :: dx, xspng, tspng
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum
    
    nx = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        !!! Hack to stop liquid's g from blowing up
        !mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        !mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        !mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

        if(decomp%yst(1)==1) then
          if(x_bc(1)==0) then
              ! rho( 1,:,:) = rhoL
              ! u  ( 1,:,:) = (u2-u1)
              v  ( 1,:,:) = zero
              w  ( 1,:,:) = zero
              do i=1,5
                  mix%material(1)%p( i,:,:) = mix%material(1)%p(6,:,:)
                  mix%material(2)%p( i,:,:) = mix%material(2)%p(6,:,:)
              end do
              
              ! mix%material(1)%g11( 1,:,:) = rho2/rho_0; mix%material(1)%g12( 1,:,:) = zero; mix%material(1)%g13( 1,:,:) = zero
              ! mix%material(1)%g21( 1,:,:) = zero; mix%material(1)%g22( 1,:,:) = one;  mix%material(1)%g23( 1,:,:) = zero
              ! mix%material(1)%g31( 1,:,:) = zero; mix%material(1)%g32( 1,:,:) = zero; mix%material(1)%g33( 1,:,:) = one
  
              ! mix%material(2)%g11( 1,:,:) = rho2/rho_0;  mix%material(2)%g12( 1,:,:) = zero; mix%material(2)%g13( 1,:,:) = zero
              ! mix%material(2)%g21( 1,:,:) = zero; mix%material(2)%g22( 1,:,:) = one;  mix%material(2)%g23( 1,:,:) = zero
              ! mix%material(2)%g31( 1,:,:) = zero; mix%material(2)%g32( 1,:,:) = zero; mix%material(2)%g33( 1,:,:) = one
              
              ! mix%material(1)%Ys ( 1,:,:) = YsL
              ! mix%material(2)%Ys ( 1,:,:) = one - YsL
  
              mix%material(1)%VF ( 1,:,:) = VFL
              mix%material(2)%VF ( 1,:,:) = one - VFL
          end if
        endif

        xspng = -two + half
        tspng = 0.2_rkind
        dx = x(2,1,1) - x(1,1,1)
        dum = half*(one - tanh( (x-xspng)/(tspng) ))

        do i=1,4
            tmp = u
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            u = u + dum*(tmp - u)

            tmp = v
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            v = v + dum*(tmp - v)

            tmp = w
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            w = w + dum*(tmp - w)

            tmp = e
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            e = e + dum*(tmp - e)

            tmp = rho
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            rho = rho + dum*(tmp - rho)

            tmp = mix%material(1)%p
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            mix%material(1)%p = mix%material(1)%p + dum*(tmp - mix%material(1)%p)

            tmp = mix%material(2)%p
            call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            mix%material(2)%p = mix%material(2)%p + dum*(tmp - mix%material(2)%p)

            do j = 1,9
                tmp = mix%material(1)%g(:,:,:,j)
                call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                mix%material(1)%g(:,:,:,j) = mix%material(1)%g(:,:,:,j) + dum*(tmp - mix%material(1)%g(:,:,:,j))

                tmp = mix%material(2)%g(:,:,:,j)
                call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                mix%material(2)%g(:,:,:,j) = mix%material(2)%g(:,:,:,j) + dum*(tmp - mix%material(2)%g(:,:,:,j))
            end do
        end do

        if(decomp%yen(1)==decomp%xsz(1)) then
          if(x_bc(2)==0) then
            rho(nx,:,:) = rhoR ! rho(nx-1,:,:)
            u  (nx,:,:) = zero ! zero
            v  (nx,:,:) = zero ! v(nx-1,:,:)
            w  (nx,:,:) = zero ! w(nx-1,:,:)
            mix%material(1)%p  (nx,:,:) = p1 ! mix%material(1)%p(nx-1,:,:)
            mix%material(2)%p  (nx,:,:) = p1 ! mix%material(2)%p(nx-1,:,:)
            
            mix%material(1)%g11(nx,:,:) = one;  mix%material(1)%g12(nx,:,:) = zero; mix%material(1)%g13(nx,:,:) = zero
            mix%material(1)%g21(nx,:,:) = zero; mix%material(1)%g22(nx,:,:) = one;  mix%material(1)%g23(nx,:,:) = zero
            mix%material(1)%g31(nx,:,:) = zero; mix%material(1)%g32(nx,:,:) = zero; mix%material(1)%g33(nx,:,:) = one
  
            mix%material(2)%g11(nx,:,:) = one;  mix%material(2)%g12(nx,:,:) = zero; mix%material(2)%g13(nx,:,:) = zero
            mix%material(2)%g21(nx,:,:) = zero; mix%material(2)%g22(nx,:,:) = one;  mix%material(2)%g23(nx,:,:) = zero
            mix%material(2)%g31(nx,:,:) = zero; mix%material(2)%g32(nx,:,:) = zero; mix%material(2)%g33(nx,:,:) = one
  
            ! mix%material(1)%Ys (nx,:,:) = YsR
            ! mix%material(2)%Ys (nx,:,:) = one - YsR
  
            mix%material(1)%VF (nx,:,:) = VFR
            mix%material(2)%VF (nx,:,:) = one - VFR
          endif
        endif

    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use SolidMixtureMod,  only: solid_mixture

    use Multispecies_shock_data

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

    use Multispecies_shock_data

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

    use Multispecies_shock_data

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

    use Multispecies_shock_data

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

    use Multispecies_shock_data

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

    use Multispecies_shock_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
