module ShockEntropy_data
    use kind_parameters,  only: rkind
    use constants,        only: one,third,half,twothird,two,three,four,seven
    implicit none
    
    real(rkind) :: pRatio = real(1.96D3,rkind)
    real(rkind) :: thick = one
    real(rkind) :: p1 = real(1.D5,rkind), p2, rho1, rho2, u1, u2, g11_1, g11_2
    real(rkind) :: Cbeta
    real(rkind) :: shmod
    real(rkind) :: xs = real(-9.50D0,rkind)
    real(rkind) :: xe = real(-8.85D0,rkind)
    real(rkind) :: rho_0
    real(rkind) :: ximp = real(-3.85D0,rkind)
    real(rkind) :: uimpact = real(300.0D0,rkind)
contains

SUBROUTINE fnumden(pf,fparams,iparams,num,den)

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
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use ShockEntropy_data

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
    
    ! Create mesh from [0,2*pi)x[0,2*pi)x[0,2*pi) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = real(20.D0,rkind)/real(nx-1,rkind)
        dy = dx
        dz = dx

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1 + i - 1, rkind ) * dx - real(10.D0,rkind)
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,eostype,eosparams,rho0,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,eps,third,half,one,two,pi
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,&
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info
    
    use ShockEntropy_data

    implicit none
    character(len=*),                                               intent(in)    :: inputfile
    type(decomp_info),                                              intent(in)    :: decomp
    real(rkind),                                                    intent(in)    :: dx,dy,dz
    integer,                                                        intent(in)    :: eostype
    real(rkind), dimension(:),                                      intent(inout) :: eosparams
    real(rkind),                                          optional, intent(inout) :: rho0, tstop, dt, tviz
    real(rkind), dimension(:,:,:,:),     intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp
    real(rkind) :: p_star, rhoRatio, tfactor, h1, grho1, grho2
    integer, dimension(2) :: iparams
    real(rkind), dimension(8) :: fparams
    integer :: nx
    real(rkind) :: mu, gam, PInf, yield, tau0

    namelist /PROBINPUT/  pRatio, p1, thick, Cbeta, xs, xe, ximp, uimpact
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    if(eostype==1) then
        gam   = eosparams(1);                           PInf = eosparams(3);
        mu = eosparams(4);   yield = eosparams(5);   tau0 = eosparams(6);
    else
    endif

    associate( rho => fields(:,:,:,rho_index),   u => fields(:,:,:,  u_index), &
                 v => fields(:,:,:,  v_index),   w => fields(:,:,:,  w_index), &
                 p => fields(:,:,:,  p_index),   T => fields(:,:,:,  T_index), &
                 e => fields(:,:,:,  e_index), g11 => fields(:,:,:,g11_index), &
               g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), & 
               g23 => fields(:,:,:,g23_index), g31 => fields(:,:,:,g31_index), & 
               g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        p_star = p1
        PInf = PInf / p_star
        mu = mu / p_star
        yield = yield / p_star
        p1 = p1 / p_star
        p2 = pRatio*p1
        rho0 = one
        rho1 = one

        fparams(1) = rho1; fparams(2) = u1; fparams(3) = p1
        fparams(4) = rho0; fparams(5) = gam; fparams(6) = PInf; fparams(7) = mu;
        fparams(8) = p2
        rho2 = rho1*min(one + p1/PInf, one) ! Init guess
        call rootfind_nr_1d(rho2,fparams,iparams)
        write(*,*) 'After root finding: ', rho2
        g11_1 = fparams(1)/fparams(4);   grho1 = g11_1**real(11.D0/3.D0,rkind) - g11_1**(-third) - g11_1**(seven*third) + g11_1**third
        g11_2 = rho2/fparams(4);         grho2 = g11_2**real(11.D0/3.D0,rkind) - g11_2**(-third) - g11_2**(seven*third) + g11_2**third
        u2 = -sqrt(rho1/rho2/(rho1-rho2)*(p1-p2+twothird*mu*(grho1-grho2)))
        u1 = rho2*u2/rho1

        print*, "Mass flux: ", rho1*u1, rho2*u2
        print*, "Momentum flux: ", rho1*u1*u1+p1, rho2*u2*u2+p2
        print*, "g flux: ", g11_1*u1, g11_2*u2
        print*, "M1 = ", u1/sqrt((gam*(p1+Pinf)+4.0_rkind/3.0_rkind*mu)/rho1)
        print*, "M2 = ", u2/sqrt((gam*(p2+Pinf)+4.0_rkind/3.0_rkind*mu)/rho2)

        ! initialize shock at xs
        tmp = half * ( one + erf( (x-xs)/(thick*dx) ) )

        u   = (one-tmp)*  (u2-u1) ! + tmp*  u1
        rho = (one-tmp)*rho2 + tmp*rho1
        v   = zero
        w   = zero
        p   = (one-tmp)*  p2 + tmp*  p1

        ! set location of impact if this is to be used instead of normal shock
        tmp = half * ( one + erf( (x-ximp)/(thick*dx) ) )
        u = (one-tmp)*two*uimpact

        ! add vorticity/entropy fluctuations starting from xe
        tmp = half * ( one + erf( (x-xe)/(eps*dx) ) )

        rho = rho*(one-tmp) + tmp*exp( -0.01_rkind * sin(13._rkind*(x-xe)))
        
        rho1 = rho(decomp%yen(1),1,1)
        u1   = u  (decomp%yen(1),1,1)
        u2   = u  (            1,1,1)

        g11 = rho/rho0; g12 = zero; g13 = zero
        g21 = zero;     g22 = one;  g23 = zero
        g31 = zero;     g32 = zero; g33 = one

        ! set pressure fluctutations to get constant sig11 - naturally reduces
        ! to constant pressure for gases and liquids since shear modulus is zero
        p = p - tmp*twothird*mu*(g11**(-third)*(g11**4-one) - g11**third*(g11**2-one))
  
        ! Get rho compatible with det(g) and rho0
        tmp = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        rho = rho0 * tmp

        rho_0 = rho0

        ! tmp = sqrt( (gam*(p+pInf) + four/three * mu)/rho )    ! Speed of sound
        ! tfactor = one / minval(tmp)
        ! tstop = tstop * tfactor
        ! tviz = tviz * tfactor

        ! write out solution of nonlinear problem, used to initialize the
        ! simulation
        print*, "rho1 = ", rho1
        print*, "rho2 = ", rho2
        print*, "u1 = ", u1
        print*, "u2 = ", u2
        print*, "p1 = ", p1
        print*, "p2 = ", p2
        print*, "a1 = ", sqrt((gam*(p1+Pinf)+4.0_rkind/3.0_rkind*mu)/rho1)
        print*, "a2 = ", sqrt((gam*(p2+Pinf)+4.0_rkind/3.0_rkind*mu)/rho2)
        print*, "Pinf = ", PInf

        print*, "rho diff: ", (rho(1,1,1) - rho2)/rho2

        print*, "tviz = ", tviz
        print*, "tstop = ", tstop

        ! store state variables at boundaries for use in hooks_bc
        nx = decomp%ysz(1)
        rho1 = rho(nx,1,1); u1 = u(nx,1,1); p1 = p(nx,1,1)
        rho2 = rho( 1,1,1); u2 = u( 1,1,1); p2 = p( 1,1,1)

    end associate

end subroutine

subroutine hook_output(decomp,der,fil,dx,dy,dz,outputdir,mesh,fields,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index, &
                                sxx_index, sxy_index, sxz_index, syy_index, syz_index, szz_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use FiltersMod,       only: filters
    use operators,        only: curl

    use ShockEntropy_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(filters),                   intent(in) :: fil
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    integer, dimension(2),           intent(in) :: x_bc, y_bc, z_bc
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile, str
    integer :: i,j,k
    real(rkind), allocatable, dimension(:,:,:,:) :: curlg
    real(rkind), allocatable, dimension(:,:,:) :: detg


    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => fields(:,:,:,g11_index), g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index),    & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), g23 => fields(:,:,:,g23_index),    &
               g31 => fields(:,:,:,g31_index), g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index),    & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3),                                      &
               sxx => fields(:,:,:, sxx_index), sxy => fields(:,:,:, sxy_index), sxz => fields(:,:,:, sxz_index), &
               syy => fields(:,:,:, syy_index), syz => fields(:,:,:, syz_index), szz => fields(:,:,:, szz_index)  )

        write(str,'(ES7.1E2,A,ES7.1E2)') pRatio, "_", Cbeta
        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ShockEntropy_"//trim(str)//"_", vizcount, ".dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        do i=1,decomp%ysz(1)
            write(outputunit,'(10ES27.16E3)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                           g11(i,1,1), g21(i,1,1), mu(i,1,1), bulk(i,1,1), kap(i,1,1)
        
        end do
        close(outputunit)

        ! do post-processing stuff
        allocate(curlg(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),9))
        allocate(detg(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))

        call curl(decomp, der, g11, g12, g13, curlg(:,:,:,1:3),-x_bc, y_bc, z_bc)
        call curl(decomp, der, g21, g22, g23, curlg(:,:,:,4:6), x_bc,-y_bc, z_bc)
        call curl(decomp, der, g31, g32, g33, curlg(:,:,:,7:9), x_bc, y_bc,-z_bc)


        detg = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)

        !! ------debug block------
        !! check if curlg13 (actually (curl(g^T))_31) is being computed correctly
        !! Get dv/dx
        !call transpose_y_to_x( g12, xtmp, decomp)
        !call der%ddx(xtmp,xdum,-x_bc(1),-x_bc(2))
        !call transpose_x_to_y(xdum, curlg(:,:,:,1) )

        !! Get du/dy
        !call der%ddy( g11, curlg(:,:,:,2), y_bc(1), y_bc(2) )

        !curlg(:,:,:,4) = curlg(:,:,:,1) - curlg(:,:,:,2) ! dv/dx - du/dy
        !write(*,*) 'curlg 2 methods diff: ', maxval(curlg(:,:,:,4)-curlg(:,:,:,3)), minval(curlg(:,:,:,4)-curlg(:,:,:,3))

        !! check commutation of filter and derivative
        !!filter curlg1 -> curlg3; curlg3 contains filter(ddx(g12))
        !curlg(:,:,:,3) = curlg(:,:,:,1)
        !call filter(decomp,curlg(:,:,:,3),fil,1,x_bc,y_bc,z_bc)

        !!filter curlg2 -> curlg4; curlg4 contains filter(ddy(g11))
        !curlg(:,:,:,4) = curlg(:,:,:,2)
        !call filter(decomp,curlg(:,:,:,4),fil,1,x_bc,y_bc,z_bc)

        !!filter g12 -> curlg9;
        !curlg(:,:,:,9) = g12
        !call filter(decomp,curlg(:,:,:,9),fil,1,x_bc,y_bc,z_bc)
        !! ddx curlg9 -> curlg5; curlg5 contains ddx(filter(g12))
        !call transpose_y_to_x( curlg(:,:,:,9), xtmp, decomp)
        !call der%ddx(xtmp,xdum,-x_bc(1),-x_bc(2))
        !call transpose_x_to_y(xdum, curlg(:,:,:,5) )

        !!filter g11 -> curlg9; ddy curlg9 -> curlg6; curlg6 contains ddy(filter(g11))
        !curlg(:,:,:,9) = g11
        !call filter(decomp,curlg(:,:,:,9),fil,1,x_bc,y_bc,z_bc)
        !call der%ddy( curlg(:,:,:,9), curlg(:,:,:,6), y_bc(1), y_bc(2) )

        !write(*,*) 'Filter commutation: F(ddx(g12)), ddx(F(g12)): ', maxval(curlg(:,:,:,3)-curlg(:,:,:,5)), minval(curlg(:,:,:,3)-curlg(:,:,:,5))
        !write(*,*) 'Filter commutation: F(ddy(g11)), ddy(F(g11)): ', maxval(curlg(:,:,:,4)-curlg(:,:,:,6)), minval(curlg(:,:,:,4)-curlg(:,:,:,6))
        !! ------debug block------

        write(outputfile,'(2A)') trim(outputdir),"/tec_ShEn_"//trim(str)//".dat"
        if(vizcount==0) then

          open(unit=outputunit, file=trim(outputfile), form='FORMATTED',STATUS='unknown')
          write(outputunit,'(290a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p","g11","g12","g13","g21","g22","g23","g31","g32","g33","sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar","curlg11","curlg12","curlg13","curlg21","curlg22","curlg23","curlg31","curlg32","curlg33","conterr"'
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES27.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(37ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
                                               g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
                                               sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), curlg(i,j,k,1:9), &
                                               rho(i,j,k)-rho_0*detg(i,j,k)
          
            end do
           end do
          end do
          close(outputunit)
        else
          open(unit=outputunit, file=trim(outputfile), form='FORMATTED',STATUS='old',ACTION='write',POSITION='append')
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          write(outputunit,'(a)') ' VARSHARELIST=([1, 2, 3]=1)'
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(34ES26.16)') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
                                               g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
                                               sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), curlg(i,j,k,1:9), &
                                               rho(i,j,k)-rho_0*detg(i,j,k)
          
            end do
           end do
          end do
          close(outputunit)
        endif

        deallocate(curlg,detg)
    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info

    use ShockEntropy_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind), dimension(2),       intent(in)    :: x_bc, y_bc, z_bc

    integer :: nx

    nx = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => fields(:,:,:,g11_index), g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), g23 => fields(:,:,:,g23_index), &
               g31 => fields(:,:,:,g31_index), g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        rho( 1,:,:) = rho2
        u  ( 1,:,:) = u2
        v  ( 1,:,:) = zero
        w  ( 1,:,:) = zero
        p  ( 1,:,:) = p2
        
        g11( 1,:,:) = rho2; g12( 1,:,:) = zero; g13( 1,:,:) = zero
        g21( 1,:,:) = zero; g22( 1,:,:) = one;  g23( 1,:,:) = zero
        g31( 1,:,:) = zero; g32( 1,:,:) = zero; g33( 1,:,:) = one

        rho(nx,:,:) = rho1
        u  (nx,:,:) = u1
        v  (nx,:,:) = zero
        w  (nx,:,:) = zero
        p  (nx,:,:) = p1
        
        g11(nx,:,:) = rho1; g12(nx,:,:) = zero; g13(nx,:,:) = zero
        g21(nx,:,:) = zero; g22(nx,:,:) = one;  g23(nx,:,:) = zero
        g31(nx,:,:) = zero; g32(nx,:,:) = zero; g33(nx,:,:) = one

    end associate
end subroutine

subroutine hook_timestep(decomp,der,mesh,fields,step,tsim,dt,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind, clen
    use constants,        only: zero, half, one, two
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use DerivativesMod,   only: derivatives

    use ShockEntropy_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind),                     intent(in) :: dt
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    integer,     dimension(2),       intent(in) :: x_bc, y_bc, z_bc

    integer :: nx, istart, iend
    real(rkind), dimension(decomp%ysz(1)) :: p_exact, dpdx
    real(rkind) :: dx, sthick, mwa
    integer :: iounit = 229
    character(len=clen) :: outputfile

    nx = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = x(2,1,1) - x(1,1,1)

        p_exact = p2
        where (x(:,1,1) .GT. half)
            p_exact = p1
        end where

        dpdx = zero
        dpdx(2:nx-1) = ( p(3:nx,1,1)-p(1:nx-2,1,1) ) / (two*dx)
        sthick = abs(p2-p1)/maxval(dx*abs(dpdx))

        istart = 1
        do while ( x(istart,1,1) .LT. half-sthick*dx )
            istart = istart + 1
        end do

        iend = nx
        do while ( x(iend,1,1) .GT. half+sthick*dx )
            iend = iend - 1
        end do

        ! mwa = maxval(abs(p(1:istart,1,1)-p_exact(1:istart)))/abs(p2-p1)
        ! mwa = max(mwa,maxval(abs(p(1:istart,1,1)-p_exact(1:istart)))/abs(p2-p1))

        mwa = maxval(p(:,1,1)-p2)
        mwa = max(mwa,maxval(p1-p(:,1,1)))
        mwa = mwa / abs(p2-p1)

        call message(2,"Shock thickness", sthick)
        call message(2,"Maximum Wiggles Amplitude", mwa)

        write(outputfile,'(A,ES8.2E2,A)') "ShockEntropy_stats_", Cbeta,".dat"
        if (step == 1) then
            open(unit=iounit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
            write(iounit,'(3A26)') "Time", "Shock thickness", "MWA"
        else
            open(unit=iounit, file=trim(outputfile), form='FORMATTED', position='APPEND', status='OLD')
        end if
        write(iounit,'(3ES26.16)') tsim, sthick, mwa
        close(iounit)

    end associate
end subroutine

subroutine hook_source(decomp,mesh,fields,tsim,rhs,rhsg)
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhsg

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
    end associate
end subroutine

subroutine hook_finalize()

end subroutine

