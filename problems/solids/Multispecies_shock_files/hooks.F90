module Multispecies_shock_data
    use kind_parameters,  only: rkind
    use constants,        only: one,two,eight
    implicit none

    real(rkind) :: p_infty = one, Rgas = one, gamma = 1.4_rkind, mu = 10._rkind, rho_0 = one, p_amb = 0.1_rkind
    real(rkind) :: minVF = 0.2_rkind, thick = one
    real(rkind) :: rhoRatio = one, pRatio = two
    logical     :: sharp = .FALSE.
    real(rkind) :: p1,p2,rho1,rho2,u1,u2,g11_1,g11_2,grho1,grho2
    real(rkind) :: rhoL, rhoR, YsL, YsR, VFL, VFR

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
    use constants,        only: one
    use decomp_2d,        only: decomp_info

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

        dx = one/real(nx-1,rkind)
        dy = dx
        dz = dx

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1 + i - 1, rkind ) * dx
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,twothird,one,two,seven
    use SolidGrid,        only: u_index,v_index,w_index
    use decomp_2d,        only: decomp_info
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

    namelist /PROBINPUT/  p_infty, Rgas, gamma, mu, rho_0, p_amb, thick, minVF, rhoRatio, pRatio
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        if (mix%ns /= 2) then
            call GracefulExit("Number of species must be 2 for this problem. Check the input file.",928)
        end if

        ! Set materials (both are the same material here)
        ! call mix%set_material(1,stiffgas(gamma,Rgas,p_infty),sep1solid(rho_0,zero,1.0D30,1.0D-10))
        call mix%set_material(1,stiffgas(gamma,Rgas,p_infty),sep1solid(rho_0,mu,1.0D30,1.0D-10))
        call mix%set_material(2,stiffgas(gamma,Rgas/rhoRatio,p_infty),sep1solid(rhoRatio*rho_0,mu,1.0D30,1.0D-10))

        p1 = p_amb
        rho1 = rho_0
        p2 = p1 * pRatio
        rho2 = rho_0
        u1 = zero
        fparams(1) = rho1; fparams(2) = u1; fparams(3) = p1
        ! fparams(4) = rho_0; fparams(5) = gamma; fparams(6) = p_infty; fparams(7) = zero;
        fparams(4) = rho_0; fparams(5) = gamma; fparams(6) = p_infty; fparams(7) = mu;
        fparams(8) = p2
        rho2 = rho1*min(one + p1/p_infty, one) ! Init guess
        call rootfind_nr_1d(rho2,fparams,iparams)
        write(*,*) 'After root finding: ', rho2
        g11_1 = fparams(1)/fparams(4);   grho1 = g11_1**real(11.D0/3.D0,rkind) - g11_1**(-third) - g11_1**(seven*third) + g11_1**third
        g11_2 = rho2/fparams(4);         grho2 = g11_2**real(11.D0/3.D0,rkind) - g11_2**(-third) - g11_2**(seven*third) + g11_2**third
        ! u2 = -sqrt(rho1/rho2/(rho1-rho2)*(p1-p2+twothird*zero*(grho1-grho2)))
        u2 = -sqrt(rho1/rho2/(rho1-rho2)*(p1-p2+twothird*mu*(grho1-grho2)))
        u1 = rho2*u2/rho1

        print*, "Mass flux: ", rho1*u1, rho2*u2
        print*, "Momentum flux: ", rho1*u1*u1+p1, rho2*u2*u2+p2
        print*, "g flux: ", g11_1*u1, g11_2*u2
        
        print*, "rho1 = ", rho1
        print*, "rho2 = ", rho2
        print*, "u1 = ", u1
        print*, "u2 = ", u2
        print*, "p1 = ", p1
        print*, "p2 = ", p2
        print*, "p_infty = ", p_infty
        print*, "rhoRatio = ", rhoRatio

        dum = half * ( one - erf( (x-0.25_rkind)/(dx) ) )

        u   = (u2-u1)*dum
        v   = zero
        w   = zero

        tmp = half * ( one - erf( (x-half)/(thick*dx) ) )

        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one
        
        mix%material(1)%g11 = (rho2*dum + rho1*(one-dum))/rho1

        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

        mix%material(2)%g11 = mix%material(1)%g11

        mix%material(1)%p  = p2*dum + p1*(one-dum)
        mix%material(2)%p  = mix%material(1)%p

        mix%material(1)%VF = minVF + (one-two*minVF)*tmp
        mix%material(2)%VF = one - mix%material(1)%VF

        tmp = rho_0*(mix%material(1)%VF*mix%material(1)%g11 + rhoRatio*mix%material(2)%g11*(one-mix%material(1)%VF)) ! Mixture density
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

subroutine hook_output(decomp,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,eps,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use Multispecies_shock_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile, str
    integer :: i

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (sharp) then
            write(str,'(I4.4,A,ES7.1E2,A,ES7.1E2,A)') decomp%ysz(1), "_", minVF, "_", rhoRatio, "_sharp"
        else
            write(str,'(I4.4,A,ES7.1E2,A,ES7.1E2,A)') decomp%ysz(1), "_", minVF, "_", rhoRatio, "_smooth"
        end if

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

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use Multispecies_shock_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    integer :: nx

    nx = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        rho( 1,:,:) = rhoL
        u  ( 1,:,:) = (u2-u1)
        v  ( 1,:,:) = zero
        w  ( 1,:,:) = zero
        mix%material(1)%p  ( 1,:,:) = p2
        mix%material(2)%p  ( 1,:,:) = p2
        
        mix%material(1)%g11( 1,:,:) = rho2/rho_0; mix%material(1)%g12( 1,:,:) = zero; mix%material(1)%g13( 1,:,:) = zero
        mix%material(1)%g21( 1,:,:) = zero; mix%material(1)%g22( 1,:,:) = one;  mix%material(1)%g23( 1,:,:) = zero
        mix%material(1)%g31( 1,:,:) = zero; mix%material(1)%g32( 1,:,:) = zero; mix%material(1)%g33( 1,:,:) = one

        mix%material(2)%g11( 1,:,:) = rho2/rho_0;  mix%material(2)%g12( 1,:,:) = zero; mix%material(2)%g13( 1,:,:) = zero
        mix%material(2)%g21( 1,:,:) = zero; mix%material(2)%g22( 1,:,:) = one;  mix%material(2)%g23( 1,:,:) = zero
        mix%material(2)%g31( 1,:,:) = zero; mix%material(2)%g32( 1,:,:) = zero; mix%material(2)%g33( 1,:,:) = one
        
        ! mix%material(1)%Ys ( 1,:,:) = YsL
        ! mix%material(2)%Ys ( 1,:,:) = one - YsL

        mix%material(1)%VF ( 1,:,:) = VFL
        mix%material(2)%VF ( 1,:,:) = one - VFL

        rho(nx,:,:) = rhoR
        u  (nx,:,:) = zero
        v  (nx,:,:) = zero
        w  (nx,:,:) = zero
        mix%material(1)%p  (nx,:,:) = p1
        mix%material(2)%p  (nx,:,:) = p1
        
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
