module EntrAdvec_data
    use kind_parameters,  only: rkind
    use constants,        only: one,third,half,twothird,two,three,four,seven
    implicit none
   
    include "fftw3.f"
 
    logical     :: broadband = .false.
    integer     :: k0 = 12
    real(rkind) :: delta = 1.0d-2
    real(rkind) :: Lx = one, rho_0
!contains
!
!SUBROUTINE fnumden(pf,fparams,iparams,num,den)
!
!  IMPLICIT NONE
!  REAL(rkind), INTENT(IN) :: pf
!  REAL(rkind), INTENT(IN), DIMENSION(:) :: fparams
!  INTEGER, INTENT(IN), DIMENSION(:) :: iparams
!  REAL(rkind), INTENT(OUT) :: num,den
!
!  INTEGER :: i, im
!  REAL(rkind) :: fac, rho1, u1, p1, p2, rho0, gam, pinf, mus, gm1, gp1, g11, frho1, grho1, Arho1, Brho1, frho2, grho2, fprho2, gprho2
!
!  real(rkind), parameter :: eleventhird = real(11.D0/3.D0,rkind), seventhird = real(7.D0/3.D0,rkind), sixth = real(one/6.D0,rkind), &
!                            sevensixth = real(7.D0/6.D0,rkind), eightthird = real(8.D0/3.D0,rkind), onetwone = real(121.D0, rkind), &
!                            thirstysix = real(36.D0,rkind), fortynine = real(49.D0,rkind), eighteen = real(18.D0,rkind),            &
!                            eighteenth = one/eighteen, ninth = real(one/9.0D0,rkind), threefourth = real(3.D0/4.D0,rkind),          &
!                            eleventwelfth = real(11.D0/12.D0,rkind), fourthird = real(4.D0/3.D0,rkind)
!
!  ! if (iparams(1)==PRESSRELAX) then
!  !   num = -one; den = zero;
!  !   i = iparams(2)
!  !   !do im = 1, NUMMAT
!  !   !  fac = vfm(i,im)/MAT_GAM(im)/(MAT_PINF(im)+pf)
!  !   !  num = num + fac*(psph(i,im)+MAT_GAM(im)*(MAT_PINF(im)+pf)-pf)
!  !   !  den = den - fac*(psph(i,im)+MAT_PINF(im))/(pf+MAT_PINF(im))
!  !   !enddo
!  ! elseif(iparams(1)==SOLIDSTATSHOCK) then
!    rho1 = fparams(1); u1 = fparams(2); p1 = fparams(3)
!    rho0 = fparams(4); gam = fparams(5); pinf = fparams(6); mus = fparams(7);
!    p2 = fparams(8)
!
!    gm1 = gam-one; gp1 = gam+one
!
!    g11 = rho1/rho0    ! g11_1
!    grho1 = g11**eleventhird - g11**(-third) - g11**seventhird + g11**third
!    frho1 = (eleventwelfth*g11**eleventhird - sixth*g11**(-third) - sevensixth*g11**seventhird -third*g11**third + threefourth*g11)
!
!    Arho1 = (half*u1**two + gam/gm1*(p1+pinf) + mus*frho1)/rho1
!    Brho1 = gam/gm1*(p1+pinf+rho1*u1**two+twothird*mus*grho1)
!
!    g11 = pf/rho0     ! g11_2
!    grho2 = g11**eleventhird - g11**(-third) - g11**seventhird + g11**third
!    frho2 = (eleventwelfth*g11**eleventhird - sixth*g11**(-third) - sevensixth*g11**seventhird -third*g11**third + threefourth*g11)
!    gprho2 = one/rho0*(eleventhird*g11**eightthird + third*g11**(-fourthird) &
!                     - seventhird*g11**fourthird + third*g11**(-twothird))
!    fprho2 = one/rho0*(onetwone/thirstysix*g11**eightthird - fortynine/eighteen*g11**fourthird &
!                     + eighteenth*g11**(-fourthird) - ninth*g11**(-twothird) + threefourth)
!
!    !! based on uL
!    !num = Arho1 - (Brho1 - half*gp1/gm1*(rho1*u1)**two/pf + mus*frho2 - twothird*gam/gm1*mus*grho2)/pf
!    !den = - one/pf*(half*gp1/gm1*(rho1*u1/pf)**two - twothird*mus*gam/gm1*gprho2 + mus*fprho2) &
!    !      + one/pf**two*(-half*gp1/gm1*(rho1*u1)**two/pf - twothird*mus*gam/gm1*grho2 + mus*frho2 + Brho1)
!
!    ! based on pR
!    num = (gam/gm1*(p1+pinf) + mus*frho1)/rho1 - half*(one/rho1+one/pf)*(p1-p2+twothird*mus*(grho1-grho2)) - mus*frho2/pf - gam/gm1/pf*(p2+pinf)
!    den = half/pf**two*(p1-p2+twothird*mus*(grho1-grho2)+two*gam/gm1*(p2+pinf)) + mus*(third*(one/rho1+one/pf)*gprho2 + frho2/pf**two - fprho2/pf)
!
!  ! endif
!
!END SUBROUTINE fnumden
!
!SUBROUTINE rootfind_nr_1d(pf,fparams,iparams)
!
!  IMPLICIT NONE
!  REAL(rkind), INTENT(INOUT) :: pf
!  REAL(rkind), INTENT(IN), DIMENSION(:) :: fparams
!  INTEGER, INTENT(IN), DIMENSION(:) :: iparams
!
!  INTEGER :: ii, itmax = 1000
!  REAL(rkind) :: tol = 1.0d-8
!  REAL(rkind) :: dpf, num, den, den_conv
!
!  !pfinitguess = pf
!  do ii = 1, itmax
!    call fnumden(pf,fparams,iparams,num,den)
!    if(dabs(den)>1.0d-12) then
!      dpf = num/den
!    else
!      write(*,*) 'den very small, please check.', num, num/den
!      stop
!    endif
!    pf = pf - dpf
!    ! check for convergence
!    if(dabs(pf)>1.0d-12) then
!      den_conv = dabs(pf)
!    else
!      den_conv = one
!    endif
!    if(dabs(dpf)/den_conv<1.0d-8) exit
!  enddo
!  if(ii==itmax+1) then
!    write(*,*) 'Newtons method for pf did not converge. Check details.', iparams(1)
!  endif
!
!END SUBROUTINE rootfind_nr_1d

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use EntrAdvec_data

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

        dx = Lx/real(nx,rkind)
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

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,eostype,eosparams,rho0,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,eps,third,half,one,two,pi
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,&
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info
    use random,           only: uniform_random
    
    use EntrAdvec_data

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
    integer :: nx, k
    real(rkind) :: mu, gam, PInf, yield, tau0, ek

    namelist /PROBINPUT/  broadband, k0, delta
    
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
        
        u   = one
        v   = zero
        w   = zero
        p   = one
        if(broadband) then
          rho = one
          call uniform_random(tmp,zero,one,321341)
          do k = 1, size(rho,1)/2
              ek = (real(k,rkind)/real(k0,rkind))**4*exp(-two*(real(k,rkind)/real(k0,rkind))**2)
              rho = rho + delta*sqrt(ek)*sin(two*pi*k*(x+tmp(k,1,1)))
          enddo
        else
          rho = one+half*sin(two*pi*x)
        endif

        !rho1 = rho(decomp%yen(1),1,1)
        !u1   = u  (decomp%yen(1),1,1)
        !u2   = u  (            1,1,1)

        g11 = rho/rho0; g12 = zero; g13 = zero
        g21 = zero;     g22 = one;  g23 = zero
        g31 = zero;     g32 = zero; g33 = one
        print *, 'g11:', maxval(g11)

        ! set pressure fluctutations to get constant sig11 - naturally reduces
        ! to constant pressure for gases and liquids since shear modulus is zero
        p = p - tmp*twothird*mu*(g11**(-third)*(g11**4-one) - g11**third*(g11**2-one))
  
        ! Get rho compatible with det(g) and rho0
        tmp = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        rho = rho0 * tmp

        rho_0 = rho0

        print*, "tviz = ", tviz
        print*, "tstop = ", tstop
        print*, "rho max = ", maxval(rho)
        print*, "g11 max = ", maxval(g11)
        print*, "g12 max = ", maxval(g12)
        print*, "g13 max = ", maxval(g13)
        print*, "g21 max = ", maxval(g21)
        print*, "g22 max = ", maxval(g22)
        print*, "g23 max = ", maxval(g23)
        print*, "g31 max = ", maxval(g31)
        print*, "g32 max = ", maxval(g32)
        print*, "g33 max = ", maxval(g33)
        print*, "detgmax = ", maxval(tmp)
        print*, "rho0    = ", rho0

        ! store state variables at boundaries for use in hooks_bc
        nx = decomp%ysz(1)
        !rho1 = rho(nx,1,1); u1 = u(nx,1,1); p1 = p(nx,1,1)
        !rho2 = rho( 1,1,1); u2 = u( 1,1,1); p2 = p( 1,1,1)

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
    use random,           only: uniform_random
    !use MKL_DFTI

    use EntrAdvec_data

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
    integer :: i,j,k, plan
    real(rkind), allocatable, dimension(:,:,:,:) :: curlg
    real(rkind), allocatable, dimension(:,:,:) :: detg, rhoex, xx
    real(rkind), allocatable, dimension(:) :: rhoin, spec
    complex(rkind), allocatable, dimension(:) :: rhoout
    real(rkind) :: trem, l2err, ek
    !type(dfti_descriptor), pointer :: my_desc1_handle


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

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/EntrAdvec_", vizcount, ".dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        do i=1,decomp%ysz(1)
            write(outputunit,'(10ES27.16E3)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                           g11(i,1,1), g21(i,1,1), mu(i,1,1), bulk(i,1,1), kap(i,1,1)
        
        end do
        close(outputunit)

        ! do post-processing stuff
        allocate(curlg(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),9))
        allocate(detg(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
        allocate(rhoex(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
        allocate(xx(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))
        allocate(rhoin(size(rho,1)), rhoout(size(rho,1)/2+1),spec(size(rho,1)))

        call curl(decomp, der, g11, g12, g13, curlg(:,:,:,1:3),-x_bc, y_bc, z_bc)
        call curl(decomp, der, g21, g22, g23, curlg(:,:,:,4:6), x_bc,-y_bc, z_bc)
        call curl(decomp, der, g31, g32, g33, curlg(:,:,:,7:9), x_bc, y_bc,-z_bc)


        ! compute exact solution error
        trem = mod(tsim, (one/Lx)) ! one/Lx is one flow-through time
        xx = x - trem
        where (xx<0)
          xx = xx + Lx
        endwhere
        if(broadband) then
            rhoex = one
            call uniform_random(detg,zero,one,321341)
            do k = 1, size(rho,1)/2
                ek = (real(k,rkind)/real(k0,rkind))**4*exp(-two*(real(k,rkind)/real(k0,rkind))**2)
                rhoex = rho + delta*sqrt(ek)*sin(two*pi*k*(xx+detg(k,1,1)))
            enddo
        else
            rhoex = one + half*sin(two*pi*xx)
        endif
        l2err = sqrt(sum((rho-rhoex)**2)/real(size(rhoex,1),rkind))
        write(*,*) 'l2error = ', l2err

        ! compute spectrum of rho
        !plan = DftiCreateDescriptor(my_desc1_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, size(rho,1))
        !plan = DftiCommitDescriptor(my_desc1_handle)
        !plan = DftiCopmuteForward(my_desc1_handle,rho(:,1,1), detg)
        !plan = DftiFreeDescriptor(my_desc1_handle)
        call dfftw_plan_dft_r2c_1d(plan, size(rho,1), rhoin, rhoout, FFTW_ESTIMATE)
        rhoin = rho(:,1,1) - one
        call dfftw_execute_dft_r2c(plan, rhoin, rhoout)
        call dfftw_destroy_plan(plan)
        spec(1:size(rho,1)/2+1) = two*abs(rhoout)**2/real(size(rho,1),rkind)
        spec(size(rho,1)/2+2:size(rho,1)) = zero

        detg = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)


        write(outputfile,'(2A,I5.5,A)') trim(outputdir),"/tec_EntrAdvec",decomp%ysz(1),".dat"
        if(vizcount==0) then

          open(unit=outputunit, file=trim(outputfile), form='FORMATTED',STATUS='unknown')
          write(outputunit,'(300a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p","g11","g12","g13","g21","g22","g23","g31","g32","g33","sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar","curlg11","curlg12","curlg13","curlg21","curlg22","curlg23","curlg31","curlg32","curlg33","conterr","rhoex","specrho"'
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES27.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(39ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
                                               g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
                                               sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), curlg(i,j,k,1:9), &
                                               rho(i,j,k)-rho_0*detg(i,j,k), rhoex(i,j,k), spec(i)
          
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
                write(outputunit,'(36ES26.16)') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
                                               g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
                                               sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), curlg(i,j,k,1:9), &
                                               rho(i,j,k)-rho_0*detg(i,j,k), rhoex(i,j,k), spec(i)
          
            end do
           end do
          end do
          close(outputunit)
        endif

        deallocate(curlg,detg,rhoex,xx,rhoin,rhoout,spec)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info

    use EntrAdvec_data

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

    use EntrAdvec_data

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

        !dx = x(2,1,1) - x(1,1,1)

        !p_exact = p2
        !where (x(:,1,1) .GT. half)
        !    p_exact = p1
        !end where

        !dpdx = zero
        !dpdx(2:nx-1) = ( p(3:nx,1,1)-p(1:nx-2,1,1) ) / (two*dx)
        !sthick = abs(p2-p1)/maxval(dx*abs(dpdx))

        !istart = 1
        !do while ( x(istart,1,1) .LT. half-sthick*dx )
        !    istart = istart + 1
        !end do

        !iend = nx
        !do while ( x(iend,1,1) .GT. half+sthick*dx )
        !    iend = iend - 1
        !end do

        !! mwa = maxval(abs(p(1:istart,1,1)-p_exact(1:istart)))/abs(p2-p1)
        !! mwa = max(mwa,maxval(abs(p(1:istart,1,1)-p_exact(1:istart)))/abs(p2-p1))

        !mwa = maxval(p(:,1,1)-p2)
        !mwa = max(mwa,maxval(p1-p(:,1,1)))
        !mwa = mwa / abs(p2-p1)

        !call message(2,"Shock thickness", sthick)
        !call message(2,"Maximum Wiggles Amplitude", mwa)

        !write(outputfile,'(A,ES8.2E2,A)') "ShockEntropy_stats_", Cbeta,".dat"
        !if (step == 1) then
        !    open(unit=iounit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
        !    write(iounit,'(3A26)') "Time", "Shock thickness", "MWA"
        !else
        !    open(unit=iounit, file=trim(outputfile), form='FORMATTED', position='APPEND', status='OLD')
        !end if
        !write(iounit,'(3ES26.16)') tsim, sthick, mwa
        !close(iounit)

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

