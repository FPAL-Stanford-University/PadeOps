module Multispecies_contact_data
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
    real(rkind) :: Ly = one, Lx = 1d0, interface_init = 5d-1, shock_init = 0.6_rkind, kwave = 4.0_rkind
    logical     :: sliding = .false.
    logical     :: ignore_gij = .false.
    real(rkind) :: uAdjust = 0d0    
    logical     :: pAdjust = .false.

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

    use Multispecies_contact_data

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
                    x(i,j,k) = real( ix1     + i - 1, rkind ) * dx 
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
    
    use Multispecies_contact_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(inout) :: tstop, dt, tviz
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit, ind
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum
    real(rkind), dimension(8) :: fparams
    real(rkind) :: fac
    integer, dimension(2) :: iparams
    logical :: adjustRgas = .TRUE.   ! If true, Rgas is used, Rgas2 adjusted to ensure p-T equilibrium

    namelist /PROBINPUT/  p_infty, Rgas, gamma, mu, rho_0, p_amb, thick, minVF, rhoRatio, pRatio,   &
                          p_infty_2, Rgas_2, gamma_2, mu_2, rho_0_2, plastic, explPlast, yield,     &
                          plastic2, explPlast2, yield2, interface_init, kwave, sliding, ignore_gij, &
                          uAdjust, pAdjust
    
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

        !Adjust gas constant for gas 2 to allow pressure equilibrium and satisfy
        !the provided densities
        Rgas_2 = Rgas * (p_amb+p_infty_2)/(p_amb+p_infty)*rho_0/rho_0_2

        ! Set materials
        call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield,1.0D-10))
        call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,1.0D-10))

        ! set logicals for plasticity
        mix%material(1)%plast = plastic; mix%material(1)%explPlast = explPlast
        mix%material(2)%plast = plastic2; mix%material(2)%explPlast = explPlast2

        ! Set logicals for sliding
        mix%material(1)%sliding = sliding
        mix%material(2)%sliding = sliding

        ! Set logical for whether or not to ignore g_ij (useful for fluids only
        ! problems)
        mix%ignore_gij = ignore_gij

        !pressure and velocity: stationary material interface
        p1 = p_amb
        p2 = p_amb
        rho1 = rho_0
        rho2 =  rho_0_2

        ! speed of sound
        a1 = sqrt((gamma*(p1+p_infty) + 4.0d0/3.0d0*mu)/rho1)
        a2 = sqrt((gamma_2*(p2+p_infty) + 4.0d0/3.0d0*mu)/rho2)

        ! write material properties
        if (nrank == 0) then
            print *, '---Simulating Stationary Contact Discontinuity---'
            print *, '---Material 1---'
            write(*,'(3(a,e12.5))') 'rho_0 = ', rho_0, ', gam  = ', gamma, ', p_infty = ', p_infty
            write(*,'(3(a,e12.5))') 'mu    = ', mu,    ', Rgas = ', Rgas,  ', a_1 = ', a1
            print *, '---Material 2---'
            write(*,'(3(a,e12.5))') 'rho_0 = ', rho_0_2, ', gam  = ', gamma_2, ', p_infty = ', p_infty_2
            write(*,'(3(a,e12.5))') 'mu    = ', mu_2,    ', Rgas = ', Rgas_2,  ', a_2 = ', a2
            write(*,*) 'p_amb = ', p_amb
        end if

        !Set mixture velocity
        u   = zero
        v   = zero
        w   = zero
        
        !set mixture pressure
        mix%material(1)%p  = p_amb
        mix%material(2)%p  = mix%material(1)%p

        ! Set up smearing function based on interface location and thickness
        tmp = half * ( one - erf( (x-(interface_init))/(thick*dx) ) )

        ! Set initial values of g (inverse deformation gradient)
        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one

        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

        if (mix%use_gTg) then
            do ind=1,2
                mix%material(ind)%g11 = mix%material(ind)%g11**2
                mix%material(ind)%g22 = mix%material(ind)%g22**2
                mix%material(ind)%g33 = mix%material(ind)%g33**2
            enddo
        end if
        
        !set mixture Volume fraction
        mix%material(1)%VF = minVF + (one-two*minVF)*tmp
        mix%material(2)%VF = one - mix%material(1)%VF

        if (mix%use_gTg) then
            tmp = rho_0*mix%material(1)%VF*sqrt(mix%material(1)%g11) + rho_0_2*sqrt(mix%material(2)%g11)*(one-mix%material(1)%VF) ! Mixture density
        else
            tmp = rho_0*mix%material(1)%VF*mix%material(1)%g11 + rho_0_2*mix%material(2)%g11*(one-mix%material(1)%VF) ! Mixture density
        end if

        mix%material(1)%Ys = mix%material(1)%VF * rho_0 / tmp
        mix%material(2)%Ys = one - mix%material(1)%Ys ! Enforce sum to unity
        
        !Set density profile based on volume fraction
        rho = tmp

        !experiment to reduce waves
        u   = -uAdjust*exp(-((x-interface_init)/(2*dx))**2)
        if (pAdjust .eq. .true.) then
            print *, "Adjusting Pressure Based on Imposed Velocity"
            mix%material(1)%p  = p_amb - rho*u**2
            mix%material(2)%p  = mix%material(1)%p
        endif

        !hard coded IC coming from soln. at t100
        !print *, (/1,2,4/)

!        u(:,1,1) = (/0e0, & !-1.24098151e-01, &
!           0e0, & !-1.60681051e-01, &
!           0e0, & !-2.23376287e-01, &
!           0e0, & !-2.91349548e-01, &
!           0e0, & !-3.60240148e-01, &
!           0e0, & !-4.29942094e-01, &
!           0e0, & !-4.99443083e-01, &
!           0e0, & !-5.68314078e-01, &
!           0e0, & !-6.36364755e-01, &
!           0e0, & !-7.03529712e-01, &
!           0e0, & !-7.69789266e-01, &
!           0e0, & !-8.35126138e-01, &
!           0e0, & !-8.99503046e-01, &
!           0e0, & !-9.62837642e-01, &
!           0e0, & !-1.02499770e+00, &
!           0e0, & !-1.08549608e+00, &
!           0e0, & !-1.14518056e+00, &
!           0e0, & !-1.19495418e+00, &
!           0e0, & !-1.26354166e+00, &
!           0e0, & !-1.27002151e+00, &
!           0e0, & !-1.32905396e+00, &
!           0e0, & !-1.39744124e+00, &
!           0e0, & !-1.24800972e+00, &
!           0e0, & !-1.47811819e+00, &
!           0e0, & !-1.26400729e+00, &
!           0e0, & !-1.44623051e+00, &
!           0e0, & !-1.30219796e+00, &
!           0e0, & !-1.29082454e+00, &
!           0e0, & !-1.47204878e+00, &
!           0e0, & !-1.13128717e+00, &
!           0e0, & !-1.63741428e+00, &
!           0e0, & !-1.10471622e+00, &
!           0e0, & !-1.65263110e+00, &
!           0e0, & !-1.25998636e+00, &
!           0e0, & !-1.41699415e+00, &
!           0e0, & !-1.54085666e+00, &
!           0e0, & !-1.08577164e+00, &
!           0e0, & !-1.96077483e+00, &
!           0e0, & !-4.45233068e-01, &
!           0e0, & !-2.98174090e+00, &
! -0.0005423111058375747 , & !      1.4040664540077083e-06 , & !   -2.81665214e-01, &
! 0.0036666130103198626  , & !      1.9556124878148797e-05 , & !   -1.45134154e+00, &
! 0.0281897990309554     , & !      0.0002723817070212334  , & !   -7.09225034e+00, &
! -0.30942266469594554   , & !      0.0037937880920979493  , & !    1.28226223e+00, &
! -3.199547019723409     , & !      0.05284062727564551    , & !   -3.70422932e+00, &
! -0.9704649754654571    , & !      -0.11439764981046495   , & !   -2.33274332e+01, &
! -9.701925685599942     , & !      -3.607673304755435     , & !   -2.09292423e+00, &
! -8.376032686683951     , & !      -16.203732818359406    , & !   -3.55687990e+01, &
! -72.5688228934471      , & !      -0.7333733321812328    , & !   -1.78421420e+02, &
! -144.77855922262165    , & !      -168.3260081638686     , & !   -1.18163160e+02, &
! -36.04835059767455     , & !      -29.565214469758214    , & !   -4.52063600e+01, &
! -8.908730457410657     , & !      0.17718252331215306    , & !   -2.08077000e+01, &
! 0.8986787247627515     , & !      -3.1847021478089332    , & !   -5.32752285e+00, &
! -2.015548472609517     , & !      -1.4085465073118417    , & !    1.90188679e+00, &
! -0.8028801325064869    , & !      -0.8227110292529738    , & !   -3.86040908e+00, &
! -0.7755585695655322    , & !      -0.26164807247639255   , & !   -8.86481715e-01, &
! -1.5254365793019313    , & !      -0.012002383142710462  , & !   -2.47696696e-01, &
! -0.10099425530514686   , & !      -0.0008617322698835296 , & !   -2.05501776e+00, &
! -0.011600059366519203  , & !      -6.186958784601842e-05 , & !   -4.57531445e-01, &
! 0.0017157071429469404  , & !      -4.442036121917162e-06 , & !   -1.49005101e+00, &
! 0.00022315956373795826 , & !      -3.1892381165386317e-07, & !   -7.55814569e-01, &
!           0e0, & !-9.14355738e-01, &
!           0e0, & !-8.07239566e-01, &
!           0e0, & !-4.85288546e-01, &
!           0e0, & !-1.87708901e+00, &
!           0e0, & !-2.26704589e+00, &
!           0e0, & !-3.06942004e+00, &
!           0e0, & !-2.66793991e+00, &
!           0e0, & !-2.22793401e+00, &
!           0e0, & !-3.05998642e+00, &
!           0e0, & !-2.35498298e+00, &
!           0e0, & !-1.92121336e+00, &
!           0e0, & !-2.84125460e+00, &
!           0e0, & !-2.22502263e+00, &
!           0e0, & !-8.90615711e-01, &
!           0e0, & !-1.59917413e+00, &
!           0e0, & !-2.95165808e+00, &
!           0e0, & !-3.93521263e+00, &
!           0e0, & !-4.90725708e+00, &
!           0e0, & !-5.38741023e+00, &
!           0e0, & !-5.52825692e+00, &
!           0e0, & !-5.50643988e+00, &
!           0e0, & !-5.37547034e+00, &
!           0e0, & !-5.19215112e+00, &
!           0e0, & !-4.97558916e+00, &
!           0e0, & !-4.73921412e+00, &
!           0e0, & !-4.48915634e+00, &
!           0e0, & !-4.22880217e+00, &
!           0e0, & !-3.96000120e+00, &
!           0e0, & !-3.68385778e+00, &
!           0e0, & !-3.40109792e+00, &
!           0e0, & !-3.11224881e+00, &
!           0e0, & !-2.81773563e+00, &
!           0e0, & !-2.51798222e+00, &
!           0e0, & !-2.21359646e+00, &
!           0e0, & !-1.90571342e+00, &
!           0e0, & !-1.59649064e+00, &
!           0e0, & !-1.29063030e+00, &
!           0e0, & !-9.88748396e-01, &
!           0e0, & !-7.10310551e-01, &
!           0e0/) !-5.47959915e-01/)

!        mix%material(1)%p(:,1,1)  = (/1.00000164e+10, &
!            1.00000212e+10, &
!            1.00000295e+10, &
!            1.00000384e+10, &
!            1.00000476e+10, &
!            1.00000568e+10, &
!            1.00000660e+10, &
!            1.00000751e+10, &
!            1.00000841e+10, &
!            1.00000929e+10, &
!            1.00001017e+10, &
!            1.00001104e+10, &
!            1.00001189e+10, &
!            1.00001273e+10, &
!            1.00001355e+10, &
!            1.00001435e+10, &
!            1.00001512e+10, &
!            1.00001591e+10, &
!            1.00001638e+10, &
!            1.00001739e+10, &
!            1.00001739e+10, &
!            1.00001711e+10, &
!            1.00001906e+10, &
!            1.00001666e+10, &
!            1.00001965e+10, &
!            1.00001721e+10, &
!            1.00001818e+10, &
!            1.00001839e+10, &
!            1.00001662e+10, &
!            1.00002048e+10, &
!            1.00001537e+10, &
!            1.00002231e+10, &
!            1.00001590e+10, &
!            1.00002142e+10, &
!            1.00001873e+10, &
!            1.00001747e+10, &
!            1.00002337e+10, &
!            1.00001310e+10, &
!            1.00002727e+10, &
!            1.00001159e+10, &
!            1.00002609e+10, &
!            1.00001668e+10, &
!            1.00002123e+10, &
!            1.00002140e+10, &
!            1.00001789e+10, &
!            1.00003120e+10, &
!            9.99983439e+09, &
!            1.00003888e+10, &
!            9.99987680e+09, &
!            1.00000128e+10, &
!            9.99996967e+09, &
!            9.99920434e+09, &
!            9.99998389e+09, &
!            9.99929312e+09, &
!            9.99976732e+09, &
!            9.99978959e+09, &
!            9.99946724e+09, &
!            1.00000436e+10, &
!            9.99951970e+09, &
!            9.99999283e+09, &
!            9.99984041e+09, &
!            9.99982311e+09, &
!            1.00001346e+10, &
!            9.99981467e+09, &
!            9.99968760e+09, &
!            9.99916871e+09, &
!            9.99901781e+09, &
!            9.99914440e+09, &
!            9.99922165e+09, &
!            9.99902942e+09, &
!            9.99914539e+09, &
!            9.99944442e+09, &
!            9.99898659e+09, &
!            9.99924819e+09, &
!            9.99981502e+09, &
!            9.99946458e+09, &
!            9.99896184e+09, &
!            9.99851905e+09, &
!            9.99812280e+09, &
!            9.99792437e+09, &
!            9.99785678e+09, &
!            9.99786258e+09, &
!            9.99790904e+09, &
!            9.99797867e+09, &
!            9.99806140e+09, &
!            9.99815233e+09, &
!            9.99824886e+09, &
!            9.99834961e+09, &
!            9.99845382e+09, &
!            9.99856105e+09, &
!            9.99867098e+09, &
!            9.99878342e+09, &
!            9.99889817e+09, &
!            9.99901508e+09, &
!            9.99913389e+09, &
!            9.99925415e+09, &
!            9.99937501e+09, &
!            9.99949461e+09, &
!            9.99961268e+09, &
!            9.99972164e+09, &
!            9.99978520e+09/)
!
!        mix%material(2)%p  = mix%material(1)%p

!        rho(:,1,1) = (/1e0, & !/1.00000856, &
!            1e0, & !1.00000842, &
!            1e0, & !1.00000818, &
!            1e0, & !1.0000079 , &
!            1e0, & !1.00000762, &
!            1e0, & !1.00000733, &
!            1e0, & !1.00000703, &
!            1e0, & !1.00000673, &
!            1e0, & !1.00000642, &
!            1e0, & !1.0000061 , &
!            1e0, & !1.00000578, &
!            1e0, & !1.00000545, &
!            1e0, & !1.00000512, &
!            1e0, & !1.00000478, &
!            1e0, & !1.00000443, &
!            1e0, & !1.00000407, &
!            1e0, & !1.00000369, &
!            1e0, & !1.00000333, &
!            1e0, & !1.00000281, &
!            1e0, & !1.00000248, &
!            1e0, & !1.00000197, &
!            1e0, & !1.0000007 , &
!            1e0, & !1.00000165, &
!            1e0, & !0.99999974, &
!            1e0, & !1.00000061, &
!            1e0, & !1.00000342, &
!            1e0, & !1.00000269, &
!            1e0, & !1.00000105, &
!            1e0, & !1.0000163 , &
!            1e0, & !1.00001635, &
!            1e0, & !0.99999967, &
!            1e0, & !1.00003614, &
!            1e0, & !1.00001184, &
!            1e0, & !0.99998766, &
!            1e0, & !1.00006354, &
!            1e0, & !0.99998838, &
!            1e0, & !0.99996692, &
!            1e0, & !1.00015279, &
!            1e0, & !0.9998769 , &
!            1e0, & !1.00003191, &
! 1.0000010806781179 , &! 0.9999999972020789  , & !              1.00069556, &
! 0.9999926934403438 , &! 0.9999999610299772  , & !              0.99965814, &
! 0.999943825419525  , &! 0.9999994572175506  , & !              1.0019863 , &
! 1.0006683553736462 , &! 0.9999924400150931  , & !              1.00435386, &
! 0.9997699497247661 , &! 0.9998947030544845  , & !              0.9999569 , &
! 1.0048228689773138 , &! 1.0001112807577506  , & !              1.0278853 , &
! 1.0028057499446328 , &! 1.0043092264685467  , & !              1.06400067, &
! 1.0622740010079168 , &! 1.0084897237459853  , & !              1.07160941, &
! 1.119055668042871  , &! 1.0514484418953114  , & !              1.64664682, &
! 3.0939554605473396 , &! 2.5044621877892363  , & !              3.78458856, &
! 7.763402009762208  , &! 8.4121440839748     , & !              7.02678696, &
! 9.973353513996502  , &! 10.002704431264823  , & !              9.38100227, &
! 9.947107564322877  , &! 10.009821347254428  , & !             10.0157938 , &
! 10.016928454355028 , &! 9.991582375638847   , & !              9.9540873 , &
! 9.992517858758465  , &! 10.003243687195367  , & !              9.98327059, &
! 10.003006576598471 , &! 9.998928244822345   , & !             10.01866755, &
! 9.99883522131865   , &! 9.999933129432621   , & !              9.99592562, &
! 9.99952889922817   , &! 9.999986845080555   , & !              9.9947804 , &
! 9.999935069272812  , &! 9.999990701690257   , & !             10.00387959, &
! 9.999999272435979  , &! 9.999990978582346   , & !             10.0000333 , &
! 9.999992075983888  , &! 9.999990998462303   , & !              9.99781284, &
!            10e0, & ! 10.00092743, &
!            10e0, & ! 10.00040658, &
!            10e0, & !  9.99910922, &
!            10e0, & ! 10.00001912, &
!            10e0, & !  9.99969741, &
!            10e0, & !  9.99910613, &
!            10e0, & !  9.99950317, &
!            10e0, & !  9.99962105, &
!            10e0, & !  9.99927314, &
!            10e0, & !  9.99945261, &
!            10e0, & !9.99965541, &
!            10e0, & !9.99927624, &
!            10e0, & !9.99942983, &
!            10e0, & !9.99978054, &
!            10e0, & !9.99955852, &
!            10e0, & !9.99925863, &
!            10e0, & !9.99899902, &
!            10e0, & !9.99876234, &
!            10e0, & !9.99864963, &
!            10e0, & !9.99861247, &
!            10e0, & !9.99862047, &
!            10e0, & !9.99865285, &
!            10e0, & !9.99869898, &
!            10e0, & !9.99875303, &
!            10e0, & !9.998812  , &
!            10e0, & !9.99887433, &
!            10e0, & !9.99893919, &
!            10e0, & !9.99900612, &
!            10e0, & !9.99907486, &
!            10e0, & !9.99914524, &
!            10e0, & !9.99921711, &
!            10e0, & !9.99929037, &
!            10e0, & !9.99936493, &
!            10e0, & !9.99944063, &
!            10e0, & !9.99951719, &
!            10e0, & !9.99959408, &
!            10e0, & !9.99967013, &
!            10e0, & !9.99974519, &
!            10e0, & !9.99981442, &
!            10e0/) !9.99985478/)

!        mix%material(1)%Ys(:,1,1) = (/1e0, & !/9.99990000e-01, &
!            1e0, & !9.99990000e-01, &
!            1e0, & !9.99990002e-01, &
!            1e0, & !9.99989990e-01, &
!            1e0, & !9.99990017e-01, &
!            1e0, & !9.99989993e-01, &
!            1e0, & !9.99989985e-01, &
!            1e0, & !9.99990032e-01, &
!            1e0, & !9.99989994e-01, &
!            1e0, & !9.99989938e-01, &
!            1e0, & !9.99990146e-01, &
!            1e0, & !9.99989860e-01, &
!            1e0, & !9.99990028e-01, &
!            1e0, & !9.99990198e-01, &
!            1e0, & !9.99989727e-01, &
!            1e0, & !9.99990161e-01, &
!            1e0, & !9.99990299e-01, &
!            1e0, & !9.99989473e-01, &
!            1e0, & !9.99990447e-01, &
!            1e0, & !9.99990501e-01, &
!            1e0, & !9.99988634e-01, &
!            1e0, & !9.99991641e-01, &
!            1e0, & !9.99989860e-01, &
!            1e0, & !9.99987896e-01, &
!            1e0, & !9.99993460e-01, &
!            1e0, & !9.99988149e-01, &
!            1e0, & !9.99987797e-01, &
!            1e0, & !9.99996141e-01, &
!            1e0, & !9.99984053e-01, &
!            1e0, & !9.99988451e-01, &
!            1e0, & !1.00000457e+00, &
!            1e0, & !9.99970340e-01, &
!            1e0, & !9.99992439e-01, &
!            1e0, & !1.00002041e+00, &
!            1e0, & !9.99938146e-01, &
!            1e0, & !1.00002129e+00, &
!            1e0, & !1.00004130e+00, &
!            1e0, & !9.99836588e-01, &
!            1e0, & !1.00015795e+00, &
!            1e0, & !9.99953861e-01, &
! 0.9999987982939283     , & !    1.0000000031112681      , & !      9.99240790e-01, &
! 1.0000081248402792     , & !    1.0000000433343774      , & !      1.00041351e+00, &
! 1.000062465717875      , & !    1.0000006035700684      , & !      9.97781846e-01, &
! 0.9992572672639375     , & !    1.0000084066472834      , & !      9.95182663e-01, &
! 1.0002514892866297     , & !    1.0001170894244724      , & !      1.00018508e+00, &
! 0.9946877346221406     , & !    0.9998789004086028      , & !      9.69901419e-01, &
! 0.9970992258424044     , & !    0.9952549533697088      , & !      9.32913160e-01, &
! 0.9328955510836218     , & !    0.9910554378906572      , & !      9.25764592e-01, &
! 0.866436731449663      , & !    0.95346002794527        , & !      5.64049791e-01, &
! 0.257713428143775      , & !    0.34960996208147144     , & !      1.82657456e-01, &
! 0.0330953218392212     , & !    0.021962393539650692    , & !      4.70213997e-02, &
! 8.830429889683507e-05  , & !    7.938357829088379e-05   , & !      7.30179206e-03, &
! 0.0005618599546638198  , & !    -0.00010125741054467109 , & !     -1.76373669e-04, &
! -0.00018368327074924748, & !    9.404748865575044e-05   , & !      5.13860224e-04, &
! 8.363302916203748e-05  , & !    -3.5996621750210675e-05 , & !      1.79165793e-04, &
! -3.3451306895741054e-05, & !    1.1914267629997799e-05  , & !     -2.12344567e-04, &
! 1.2932304117953728e-05 , & !    7.433858820372318e-07   , & !      4.35742384e-05, &
! 5.237346793848187e-06  , & !    1.4619310763507283e-07  , & !      5.72252457e-05, &
! 7.218193816258095e-07  , & !    1.0331659892097916e-07  , & !     -4.84934644e-05, &
! 8.029887964960176e-09  , & !    1.0023820436998543e-07  , & !     -2.24083799e-07, &
! 8.80376577025564e-08   , & !    1.0001718560655235e-07  , & !      2.29246022e-05, &
!           1e-7, & !-1.25527457e-05, &
!           1e-7, & !-3.77755926e-06, &
!           1e-7, & ! 8.45389696e-06, &
!           1e-7, & !-2.87746796e-06, &
!           1e-7, & !-2.36910622e-06, &
!           1e-7, & ! 3.20309751e-06, &
!           1e-7, & !-5.78705785e-07, &
!           1e-7, & !-1.16114912e-06, &
!           1e-7, & ! 1.38880572e-06, &
!           1e-7, & !-1.59306741e-07, &
!           1e-7, & !-3.96464310e-07, &
!           1e-7, & ! 5.93487070e-07, &
!           1e-7, & ! 1.45386163e-08, &
!           1e-7, & !-1.11891145e-07, &
!           1e-7, & ! 3.08050123e-07, &
!           1e-7, & ! 6.51473807e-08, &
!           1e-7, & ! 1.10939415e-09, &
!           1e-7, & ! 2.08705220e-07, &
!           1e-7, & ! 6.17564536e-08, &
!           1e-7, & ! 7.48711761e-08, &
!           1e-7, & ! 1.40310574e-07, &
!           1e-7, & ! 8.24214819e-08, &
!           1e-7, & ! 9.19828521e-08, &
!           1e-7, & ! 1.16259370e-07, &
!           1e-7, & ! 9.26287460e-08, &
!           1e-7, & ! 9.49851705e-08, &
!           1e-7, & ! 1.10835207e-07, &
!           1e-7, & ! 9.17014358e-08, &
!           1e-7, & ! 1.02290863e-07, &
!           1e-7, & ! 1.02098243e-07, &
!           1e-7, & ! 9.70644107e-08, &
!           1e-7, & ! 1.01483240e-07, &
!           1e-7, & ! 1.00067253e-07, &
!           1e-7, & ! 9.95576803e-08, &
!           1e-7, & ! 9.97612522e-08, &
!           1e-7, & ! 1.01138618e-07, &
!           1e-7, & ! 9.84560110e-08, &
!           1e-7, & ! 1.01277153e-07, &
!           1e-7, & ! 9.93669452e-08, &
!           1e-7/)  ! 1.00000090e-07/)

!        mix%material(2)%Ys = one - mix%material(1)%Ys ! Enforce sum to unity
        
!        mix%material(1)%VF(:,1,1) = (/1e0, & !/9.99999000e-01, &
!            1e0, & !9.99999000e-01, &
!            1e0, & !9.99999000e-01, &
!            1e0, & !9.99998999e-01, &
!            1e0, & !9.99999002e-01, &
!            1e0, & !9.99998999e-01, &
!            1e0, & !9.99998998e-01, &
!            1e0, & !9.99999003e-01, &
!            1e0, & !9.99998999e-01, &
!            1e0, & !9.99998994e-01, &
!            1e0, & !9.99999015e-01, &
!            1e0, & !9.99998986e-01, &
!            1e0, & !9.99999003e-01, &
!            1e0, & !9.99999020e-01, &
!            1e0, & !9.99998973e-01, &
!            1e0, & !9.99999016e-01, &
!            1e0, & !9.99999030e-01, &
!            1e0, & !9.99998947e-01, &
!            1e0, & !9.99999045e-01, &
!            1e0, & !9.99999050e-01, &
!            1e0, & !9.99998863e-01, &
!            1e0, & !9.99999164e-01, &
!            1e0, & !9.99998986e-01, &
!            1e0, & !9.99998790e-01, &
!            1e0, & !9.99999346e-01, &
!            1e0, & !9.99998815e-01, &
!            1e0, & !9.99998780e-01, &
!            1e0, & !9.99999614e-01, &
!            1e0, & !9.99998405e-01, &
!            1e0, & !9.99998845e-01, &
!            1e0, & !1.00000046e+00, &
!            1e0, & !9.99997034e-01, &
!            1e0, & !9.99999244e-01, &
!            1e0, & !1.00000204e+00, &
!            1e0, & !9.99993814e-01, &
!            1e0, & !1.00000213e+00, &
!            1e0, & !1.00000413e+00, &
!            1e0, & !9.99983656e-01, &
!            1e0, & !1.00001579e+00, &
!            1e0, & !9.99995386e-01, &
!   0.9999998799246534      , & !   1.0000000003108804     , &!    9.99924027e-01, &
!   1.0000008118399617      , & !   1.0000000043300028     , &!    1.00004134e+00, &
!   1.0000062416200528      , & !   1.0000000603091612     , &!    9.99777741e-01, &
!   0.9999257382918172      , & !   1.0000008399983233     , &!    9.99516169e-01, &
!   1.0000255611416926      , & !   1.0000116996606132     , &!    1.00001850e+00, &
!   0.9994641256691874      , & !   0.9999876354713613     , &!    9.96906339e-01, &
!   0.9996882500061519      , & !   0.9995211970590505     , &!    9.92860230e-01, &
!   0.9930806665546761      , & !   0.9990566973615571     , &!    9.92044969e-01, &
!   0.9867715924396808      , & !   0.9942835064560767     , &!    9.28255751e-01, &
!   0.7673382821614069      , & !   0.8328375346900851     , &!    6.90859239e-01, &
!   0.24851088780419894     , & !   0.17642843511391124    , &!    3.30393790e-01, &
!   0.0029607206670556815   , & !   -0.0003004923627577169 , &!    6.85153566e-02, &
!   0.0058769372974581385   , & !   -0.0010912608060475568 , &!   -1.76654083e-03, &
!   -0.0018809393727808536  , & !   0.0009352915956828852  , &!    5.11494693e-03, &
!   0.000831349026837101    , & !   -0.000360409688374345  , &!    1.78877355e-03, &
!   -0.00033406406649659736 , & !   0.00011908390862866926 , &!   -2.12751156e-03, &
!   0.00012941985348348894  , & !   7.430063041936673e-06  , &!    4.35571567e-04, &
!   5.234453020313774e-05   , & !   1.461657716358542e-06  , &!    5.71957883e-04, &
!   7.2145252433780735e-06  , & !   1.0331455274634642e-06 , &!   -4.85146382e-04, &
!   8.084044692049673e-08   , & !   1.002379739229046e-06  , &!   -2.24084251e-06, &
!   8.804462349636195e-07   , & !   1.000170855228319e-06  , &!    2.29198733e-04, &
!          1e-6, & !-1.25541640e-04, &
!          1e-6, & !-3.77768769e-05, &
!          1e-6, & ! 8.45325379e-05, &
!          1e-6, & !-2.87754248e-05, &
!          1e-6, & !-2.36915674e-05, &
!          1e-6, & ! 3.20300518e-05, &
!          1e-6, & !-5.78708799e-06, &
!          1e-6, & !-1.16116126e-05, &
!          1e-6, & ! 1.38878837e-05, &
!          1e-6, & !-1.59306969e-06, &
!          1e-6, & ! -3.96465724e-06, &
!          1e-6, & !  5.93483900e-06, &
!          1e-6, & !  1.45386144e-07, &
!          1e-6, & ! -1.11891257e-06, &
!          1e-6, & !  3.08049269e-06, &
!          1e-6, & !  6.51473425e-07, &
!          1e-6, & !  1.10939414e-08, &
!          1e-6, & !  2.08704828e-06, &
!          1e-6, & !  6.17564193e-07, &
!          1e-6, & !  7.48711257e-07, &
!          1e-6, & !  1.40310396e-06, &
!          1e-6, & !  8.24214208e-07, &
!          1e-6, & !  9.19827759e-07, &
!          1e-6, & !  1.16259248e-06, &
!          1e-6, & !  9.26286688e-07, &
!          1e-6, & !  9.49850893e-07, &
!          1e-6, & !  1.10835096e-06, &
!          1e-6, & !  9.17013601e-07, &
!          1e-6, & !  1.02290769e-06, &
!          1e-6, & !  1.02098150e-06, &
!          1e-6, & !  9.70643259e-07, &
!          1e-6, & !  1.01483147e-06, &
!          1e-6, & !  1.00067162e-06, &
!          1e-6, & !  9.95575911e-07, &
!          1e-6, & !  9.97611626e-07, &
!          1e-6, & !  1.01138526e-06, &
!          1e-6, & !  9.84559237e-07, &
!          1e-6, & !  1.01277061e-06, &
!          1e-6, & !  9.93668564e-07, &
!          1e-6/)  !  1.00000000e-06/)

!        mix%material(2)%VF = one - mix%material(1)%VF ! Enforce sum to unity

        !Set "left and right" quantities
        YsL  = mix%material(1)%Ys(1,1,1)
        YsR  = mix%material(1)%Ys(decomp%ysz(1),1,1)
        VFL  = mix%material(1)%VF(1,1,1)
        VFR  = mix%material(1)%VF(decomp%ysz(1),1,1)
        rhoL = rho(1,1,1)
        rhoR = rho(decomp%ysz(1),1,1)


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

    use Multispecies_contact_data

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
           write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Multispecies_contact_"//trim(str)//"_", vizcount, ".dat"

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
           
       write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Multispecies_contact_statistics.dat"

       if (vizcount == 0) then
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
           write(outputunit,'(7A27)') 'tsim', 'mixwidth', 'vort_pos', 'vort_neg', 'Al_mass', 'xspike', 'xbubbl'
       else
           open(unit=outputunit, file=trim(outputfile), form='FORMATTED', action='WRITE', status='OLD', position='APPEND')
       end if
       write(outputunit,'(7ES27.16E3)') tsim, mixwidth, vort_pos, vort_neg, Al_mass, xspike, xbubbl
       close(outputunit)

        write(outputfile,'(4A)') trim(outputdir),"/tec_MultSpecShock_"//trim(str),".dat"
        if(vizcount==0) then
          open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='replace')
          write(outputunit,'(350a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p", &
                                     "sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar", &
                                     "p-1","Ys-1","VF-1","eh-1","T-1","g11-1","g12-1","g13-1","g21-1","g22-1","g23-1","g31-1","g32-1","g33-1","Dstar-1","kap-1","rhom-1",&
                                     "p-2","Ys-2","VF-2","eh-2","T-2","g11-2","g12-2","g13-2","g21-2","g22-2","g23-2","g31-2","g32-2","g33-2","Dstar-2","kap-2","rhom-2"'
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(52ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &                      ! continuum (9)
                                                sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), &      ! continuum (9)
                                                mix%material(1)% p(i,j,k),  mix%material(1)% Ys(i,j,k), mix%material(1)% VF(i,j,k), mix%material(1)% eh(i,j,k), &  ! material 1 (14)
                                                mix%material(1)% T(i,j,k),  mix%material(1)%g11(i,j,k), mix%material(1)%g12(i,j,k), mix%material(1)%g13(i,j,k), &  ! material 1 
                                                mix%material(1)%g21(i,j,k), mix%material(1)%g22(i,j,k), mix%material(1)%g23(i,j,k), mix%material(1)%g31(i,j,k), &  ! material 1 
                                                mix%material(1)%g32(i,j,k), mix%material(1)%g33(i,j,k), mix%material(1)%diff(i,j,k), mix%material(1)%kap(i,j,k), mix%material(1)%rhom(i,j,k), & ! material 1 
                                                mix%material(2)% p(i,j,k),  mix%material(2)% Ys(i,j,k), mix%material(2)% VF(i,j,k), mix%material(2)% eh(i,j,k), &  ! material 2 (14)
                                                mix%material(2)% T(i,j,k),  mix%material(2)%g11(i,j,k), mix%material(2)%g12(i,j,k), mix%material(2)%g13(i,j,k), &  ! material 2
                                                mix%material(2)%g21(i,j,k), mix%material(2)%g22(i,j,k), mix%material(2)%g23(i,j,k), mix%material(2)%g31(i,j,k), &  ! material 2
                                                mix%material(2)%g32(i,j,k), mix%material(2)%g33(i,j,k), mix%material(2)%diff(i,j,k), mix%material(2)%kap(i,j,k), mix%material(2)%rhom(i,j,k)    ! material 2
            end do
           end do
          end do
          close(outputunit)
        else
          open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='old', action='write', position='append')
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          write(outputunit,'(a)') ' VARSHARELIST=([1, 2, 3]=1)'
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(49ES26.16)') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &                                                    ! continuum (6)
                                                sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), &      ! continuum (9)
                                                mix%material(1)% p(i,j,k),  mix%material(1)% Ys(i,j,k), mix%material(1)% VF(i,j,k), mix%material(1)% eh(i,j,k), &  ! material 1 (14)
                                                mix%material(1)% T(i,j,k),  mix%material(1)%g11(i,j,k), mix%material(1)%g12(i,j,k), mix%material(1)%g13(i,j,k), &  ! material 1 
                                                mix%material(1)%g21(i,j,k), mix%material(1)%g22(i,j,k), mix%material(1)%g23(i,j,k), mix%material(1)%g31(i,j,k), &  ! material 1 
                                                mix%material(1)%g32(i,j,k), mix%material(1)%g33(i,j,k), mix%material(1)%diff(i,j,k), mix%material(1)%kap(i,j,k), mix%material(1)%rhom(i,j,k),&  ! material 1 
                                                mix%material(2)% p(i,j,k),  mix%material(2)% Ys(i,j,k), mix%material(2)% VF(i,j,k), mix%material(2)% eh(i,j,k), &  ! material 2 (14)
                                                mix%material(2)% T(i,j,k),  mix%material(2)%g11(i,j,k), mix%material(2)%g12(i,j,k), mix%material(2)%g13(i,j,k), &  ! material 2
                                                mix%material(2)%g21(i,j,k), mix%material(2)%g22(i,j,k), mix%material(2)%g23(i,j,k), mix%material(2)%g31(i,j,k), &  ! material 2
                                                mix%material(2)%g32(i,j,k), mix%material(2)%g33(i,j,k), mix%material(2)%diff(i,j,k), mix%material(2)%kap(i,j,k), mix%material(2)%rhom(i,j,k)    ! material 2
            end do
           end do
          end do
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

    use Multispecies_contact_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc
    
    integer :: nx, i, j
    real(rkind) :: dx, xspng, tspng, xspngL, xspngR
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum, dumL, dumR
    
    nx = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if(decomp%yst(1)==1) then
          if(x_bc(1)==0) then
              rho( 1,:,:) = rhoL
              u  ( 1,:,:) = zero
              v  ( 1,:,:) = zero
              w  ( 1,:,:) = zero
              mix%material(1)%p  (1,:,:) = p1 ! mix%material(1)%p(nx-1,:,:)
              mix%material(2)%p  (1,:,:) = p1 ! mix%material(2)%p(nx-1,:,:)
              
              mix%material(1)%g11( 1,:,:) = one; mix%material(1)%g12( 1,:,:) = zero; mix%material(1)%g13( 1,:,:) = zero
              mix%material(1)%g21( 1,:,:) = zero; mix%material(1)%g22( 1,:,:) = one;  mix%material(1)%g23( 1,:,:) = zero
              mix%material(1)%g31( 1,:,:) = zero; mix%material(1)%g32( 1,:,:) = zero; mix%material(1)%g33( 1,:,:) = one
  
              mix%material(2)%g11( 1,:,:) = one;  mix%material(2)%g12( 1,:,:) = zero; mix%material(2)%g13( 1,:,:) = zero
              mix%material(2)%g21( 1,:,:) = zero; mix%material(2)%g22( 1,:,:) = one;  mix%material(2)%g23( 1,:,:) = zero
              mix%material(2)%g31( 1,:,:) = zero; mix%material(2)%g32( 1,:,:) = zero; mix%material(2)%g33( 1,:,:) = one
              
              mix%material(1)%Ys ( 1,:,:) = YsL
              mix%material(2)%Ys ( 1,:,:) = one - YsL
  
              mix%material(1)%VF ( 1,:,:) = VFL       !Dirichlet BC: set left volume fraction to VFL
              mix%material(2)%VF ( 1,:,:) = one - VFL !Dirichlet BC: set left volume fraction to 1- VFL
          end if
        endif
        
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
  
            mix%material(1)%Ys (nx,:,:) = YsR
            mix%material(2)%Ys (nx,:,:) = one - YsR
  
            mix%material(1)%VF (nx,:,:) = VFR
            mix%material(2)%VF (nx,:,:) = one - VFR
          endif
        endif

        ! apply sponge at left and right boundaries to damp outgoing waves
        xspngL = 0.15d0
        xspngR = 0.85d0
        tspng = 0.03d0
        dx = x(2,1,1) - x(1,1,1)
        dumL = half*(one - tanh( (x-xspngL)/(tspng) ))
        dumR = half*(one + tanh( (x-xspngR)/(tspng) ))
        dum  = dumL+dumR

        do i=1,10
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


    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use constants,        only: zero, one
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use SolidMixtureMod,  only: solid_mixture

    use Multispecies_contact_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix
    integer                                     :: imin, ind(1), i
    !logical                                     :: ignore_gij = .true.    

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

        !Use this if the simulation only involves fluids. In this case, g_ij is
        !meaningless and can only mess up the simulation. Therefore, reset it to
        !1 at each time step so it stays well-behaved.
        if (mix%ignore_gij) then
            !enforce gij = delta_ij (applies regardless of g or gTg formulation)
            do i=1, mix%ns
                mix%material(i)%g11 = one;  mix%material(i)%g12 = zero; mix%material(i)%g13 = zero
                mix%material(i)%g21 = zero; mix%material(i)%g22 = one;  mix%material(i)%g23 = zero
                mix%material(i)%g31 = zero; mix%material(i)%g32 = zero; mix%material(i)%g33 = one
            enddo
        endif

    end associate
end subroutine

subroutine hook_mixture_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,&
                                mom_index,TE_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use Multispecies_contact_data

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

    use Multispecies_contact_data

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

    use Multispecies_contact_data

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

    use Multispecies_contact_data

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

    use Multispecies_contact_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
