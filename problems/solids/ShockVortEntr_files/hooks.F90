module ShockVortEntr_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight,third,half,two,twothird,pi
    use FiltersMod,       only: filters
    implicit none
    
    real(rkind) :: pinit   = real(1.0D5,rkind)
    real(rkind) :: rho_0
    real(rkind) :: p1 = real(1.D5,rkind), pRatio = 2.0_rkind, thick = one, xs = one
    real(rkind) :: rho1, rho2, u1, u2, g11_1, g11_2, p2
    real(rkind) :: Av = 0.025_rkind, Ae = 0.025_rkind, psiang = 45.0_rkind*pi/180.0_rkind, k2wv = 1.0_rkind, k1wv = 1.0_rkind, sinpsi, cospsi
    type(filters) :: mygfil

    ! statistics
    real(rkind), allocatable, dimension(:,:,:,:) :: vort
    real(rkind), allocatable, dimension(:,:,:)   :: stats
    real(rkind), allocatable, dimension(:,:)     :: umean
    real(rkind)                                  :: tavg, tstatstart = 30.0_rkind
    logical                                      :: istatstart = .FALSE., istatwritestart = .FALSE.

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
    use constants,        only: half,one,pi
    use decomp_2d,        only: decomp_info

    use ShockVortEntr_data

    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    real(rkind) :: xa, xb, yc, yd

    xa = 0.0D0; xb = 8.0D0*pi+1.0D0
    yc = -pi; yd = pi

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,2*pi)x[0,2*pi)x[0,2*pi) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = (xb-xa)/real(nx-1,rkind)
        dy = (yd-yc)/real(ny,rkind)
        !dx = real(1.0d0,rkind)/real(nx-1,rkind)
        !dy = dx
        dz = dx

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = xa + real( ix1 - 1 + i - 1, rkind ) * dx
                    y(i,j,k) = yc + real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,eostype,eosparams,rho0,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,one,two,pi,eight,seven
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,&
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info, nrank
    
    use ShockVortEntr_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:),       intent(in)    :: eosparams
    integer,                         intent(in)    :: eostype
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),           optional, intent(inout) :: rho0, tstop, dt, tviz

    integer :: ioUnit, i, j, k, nx
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, entr
    real(rkind) :: stp_x,bkstp_x,bkstp_y,phiang,rad,stp_r1,bkstp_r2,regfrac,dxdy
    real(rkind) :: p_star, grho1, grho2, entr1, entr2
    integer, dimension(2) :: iparams
    real(rkind), dimension(8) :: fparams
    real(rkind) :: gam, PInf, mu, yield, tau0
    real(rkind), dimension(decomp%ysz(1)) :: mask

    namelist /PROBINPUT/  pRatio, p1, thick, xs, Av, Ae, psiang, k1wv, k2wv, tstatstart
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    if(eostype==1) then
        gam   = eosparams(1);                           PInf = eosparams(3);
        mu = eosparams(4);   yield = eosparams(5);   tau0 = eosparams(6);
    else
    endif

    ! Initialize mygfil
    call mygfil%init(                        decomp, &
                     .FALSE.,     .TRUE.,    .TRUE., &
                  "gaussian", "gaussian", "gaussian" )

    ! Initialize arrays for statistics
    allocate(vort(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),3))
    allocate(stats(decomp%ysz(1),decomp%ysz(3),4))
    allocate(umean(decomp%ysz(1),3))
    tavg = zero
    stats = zero

    rho_0 = rho0
    associate( rho => fields(:,:,:,rho_index),   u => fields(:,:,:,  u_index), &
                 v => fields(:,:,:,  v_index),   w => fields(:,:,:,  w_index), &
                 p => fields(:,:,:,  p_index),   T => fields(:,:,:,  T_index), &
                 e => fields(:,:,:,  e_index), g11 => fields(:,:,:,g11_index), &
               g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), & 
               g23 => fields(:,:,:,g23_index), g31 => fields(:,:,:,g31_index), & 
               g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        

        ! solve non-linear problem for determining shock jump
        p_star = one;         
        PInf = PInf / p_star;    mu = mu / p_star;    yield = yield / p_star
        p1 = p1 / p_star;    p2 = pRatio*p1;
        rho0 = rho_0;          rho1 = rho_0

        fparams(1) = rho1; fparams(2) = u1; fparams(3) = p1
        fparams(4) = rho0; fparams(5) = gam; fparams(6) = PInf; fparams(7) = mu;
        fparams(8) = p2
        if(PInf<1.0d-10) then
            rho2 = rho1 ! Init guess
        else
            rho2 = rho1*min(one + p1/PInf, one) ! Init guess
        endif
        call rootfind_nr_1d(rho2,fparams,iparams)
        if(nrank==0) write(*,*) 'After root finding: ', rho2
        g11_1 = fparams(1)/fparams(4);   grho1 = g11_1**real(11.D0/3.D0,rkind) - g11_1**(-third) - g11_1**(seven*third) + g11_1**third
        g11_2 = rho2/fparams(4);         grho2 = g11_2**real(11.D0/3.D0,rkind) - g11_2**(-third) - g11_2**(seven*third) + g11_2**third
        u2 = sqrt(rho1/rho2/(rho1-rho2)*(p1-p2+twothird*mu*(grho1-grho2)))
        u1 = rho2*u2/rho1

        if(nrank==0) then
          ! write out solution of nonlinear problem, used to initialize the
          ! simulation
          print*, "Mass flux: ", rho1*u1, rho2*u2
          print*, "Momentum flux: ", rho1*u1*u1+p1, rho2*u2*u2+p2
          print*, "g flux: ", g11_1*u1, g11_2*u2
          print*, "rho1 = ", rho1
          print*, "rho2 = ", rho2
          print*, "u1 = ", u1
          print*, "u2 = ", u2
          print*, "p1 = ", p1
          print*, "p2 = ", p2
          print*, "a1 = ", sqrt((gam*(p1+Pinf)+4.0_rkind/3.0_rkind*mu)/rho1)
          print*, "a2 = ", sqrt((gam*(p2+Pinf)+4.0_rkind/3.0_rkind*mu)/rho2)
          print*, "M1 = ", u1/sqrt((gam*(p1+Pinf)+4.0_rkind/3.0_rkind*mu)/rho1)
          print*, "M2 = ", u2/sqrt((gam*(p2+Pinf)+4.0_rkind/3.0_rkind*mu)/rho2)
          print*, "Pinf = ", PInf
        endif

        ! initialize shock at xs
        tmp = half * ( one + erf( (x-xs)/(thick*dx) ) )

        u   = (one-tmp)*u1   + tmp*u2
        rho = (one-tmp)*rho1 + tmp*rho2
        v   = zero
        w   = zero

        ! determine p as a function of entropy
        entr1 = log((p1/p1)**(one/gam)*rho1/rho1) ! zero
        entr2 = log((p2/p1)**(one/gam)*rho1/rho2) ! non-zero
        print *, 'entr1 = ', entr1
        print *, 'entr2 = ', entr2
        entr   = (one-tmp)*entr1   + tmp*entr2
        p = p1*(rho/rho1*exp(entr))**gam

        !rho1 = rho(decomp%yen(1),1,1)
        !u1   = u  (decomp%yen(1),1,1)
        !u2   = u  (            1,1,1)

        g11 = rho/rho0; g12 = zero; g13 = zero
        g21 = zero;     g22 = one;  g23 = zero
        g31 = zero;     g32 = zero; g33 = one

        ! Get rho compatible with det(g) and rho0
        tmp = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        rho = rho0 * tmp

        rho_0 = rho0

        if(nrank==0) then
          print*, "tviz = ", tviz
          print*, "tstop = ", tstop
        endif

        ! store state variables at boundaries for use in hooks_bc
        nx = decomp%ysz(1)
        rho2 = rho(nx,1,1); u2 = u(nx,1,1); p2 = p(nx,1,1); g11_2 = rho2/rho0
        rho1 = rho( 1,1,1); u1 = u( 1,1,1); p1 = p( 1,1,1); g11_1 = rho2/rho0

        ! store base state for computation of fluctuating ke
        umean(:,1) = u(:,1,1);     umean(:,2) = v(:,1,1);     umean(:,3) = w(:,1,1)


        ! Add vorticity/entropy fluctuations
        psiang = psiang*pi/180.0_rkind
        sinpsi = sin(psiang); cospsi = cos(psiang)
        if(tan(psiang)<1.0d-6) then
            k1wv = 1.0d0; k2wv = 0.0d0
        elseif(tan(psiang)>1.0d6) then
            k1wv = 0.0d0; k2wv = 1.0d0
        else
            k1wv = k2wv/tan(psiang)
        endif
        if(nrank==0) write(*,*) '(k1, k2): ', k1wv, k2wv

        mask = 0.0d0
        where(x(:,1,1) .lt. xs-one)
          mask = 1.0d0
        endwhere

        do k=1,decomp%ysz(3)
         do j=1,decomp%ysz(2)
          do i=1,decomp%ysz(1)

             !rho(i,j,k) = rho(i,j,k) + mask(i)*(rho1*exp(-Ae* sin(k1wv*x(i,j,k)))  -rho1)
             rho(i,j,k) = rho(i,j,k) + rho1 * Ae *          cos(k1wv*x(i,j,k) + k2wv*y(i,j,k))
             u(i,j,k)   = u(i,j,k)   +   u1 * Av * sinpsi * cos(k1wv*x(i,j,k) + k2wv*y(i,j,k))
             v(i,j,k)   =            -   u1 * Av * cospsi * cos(k1wv*x(i,j,k) + k2wv*y(i,j,k))
          end do 
         end do 
        end do 

        !!g11 = one;  g12 = zero; g13 = zero
        !!g21 = zero; g22 = one;  g23 = zero
        !!g31 = zero; g32 = zero; g33 = one

        !!! Get rho compatible with det(g) and rho0
        !!tmp = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        !!rho = rho0 * tmp
        g11 = rho/rho0

    end associate

end subroutine

subroutine filter(decomp,myfil,arr,numtimes,x_bc_,y_bc_,z_bc_)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use FiltersMod,       only: filters

    implicit none
    type(decomp_info),                                                    intent(in) :: decomp
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), intent(inout) :: arr
    type(filters),                                                        intent(in) :: myfil
    integer,                                                              intent(in) :: numtimes
    integer, dimension(2),                                                intent(in) :: x_bc_, y_bc_, z_bc_

    integer :: times2fil
    integer, dimension(2) :: x_bc, y_bc, z_bc
    real(rkind), dimension(:,:,:), pointer :: tmp_in_y, tmp1_in_x, tmp1_in_z, tmp2_in_x, tmp2_in_z
    real(rkind), dimension(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3),2), target :: xbuf
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),1), target :: ybuf
    real(rkind), dimension(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3),2), target :: zbuf
    integer :: lastx, lasty, lastz, idx


    !if (present(myfil)) then
    !    fil2use => myfil
    !else
    !    !fil2use => this%fil
    !end if 

    !if (present(numtimes)) then
        times2fil = numtimes
    !else
    !    times2fil = 1
    !end if

    ! Allocate pointers for the needed buffers 
    ! Atleast 2 buffers in x and z are assumed
    ! Last two buffers are occupied

    !lastx = size(this%xbuf,4)
    !lasty = size(this%ybuf,4)
    !lastz = size(this%zbuf,4)

    tmp1_in_x => xbuf(:,:,:,1)
    tmp2_in_x => xbuf(:,:,:,2)
    tmp_in_y => ybuf(:,:,:,1)
    tmp1_in_z => zbuf(:,:,:,1)
    tmp2_in_z => zbuf(:,:,:,2)
   
    x_bc = x_bc_
    y_bc = y_bc_
    z_bc = z_bc_
   

    !! Order y -> x 
    !  ! First filter in y
    !  call myfil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
    !  ! Subsequent refilters 
    !  do idx = 1,times2fil-1
    !      arr = tmp_in_y
    !      call myfil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
    !  end do
    !  
    !  ! Then transpose to x
    !  call transpose_y_to_x(tmp_in_y,tmp1_in_x,decomp)

    !  ! First filter in x
    !  call myfil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
    !  ! Subsequent refilters
    !  do idx = 1,times2fil-1
    !      tmp1_in_x = tmp2_in_x
    !      call myfil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
    !  end do 

    !  ! Now transpose back to y
    !  call transpose_x_to_y(tmp2_in_x,tmp_in_y,decomp)
    !! ---------Order y -> x 

    ! Order x -> y 
      ! Transpose to x
      call transpose_y_to_x(arr,tmp1_in_x,decomp)

      ! First filter in x
      call myfil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
      ! Subsequent refilters
      do idx = 1,times2fil-1
          tmp1_in_x = tmp2_in_x
          call myfil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
      end do 

      ! Now transpose back to y
      call transpose_x_to_y(tmp2_in_x,arr,decomp)

      ! First filter in y
      call myfil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
      ! Subsequent refilters 
      do idx = 1,times2fil-1
          arr = tmp_in_y
          call myfil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
      end do
    ! ---------Order x -> y 
      

    ! Now transpose to z
    call transpose_y_to_z(tmp_in_y,tmp1_in_z,decomp)

    !First filter in z
    call myfil%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
    ! Subsequent refilters
    do idx = 1,times2fil-1
        tmp1_in_z = tmp2_in_z
        call myfil%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
    end do 

    ! Now transpose back to y
    call transpose_z_to_y(tmp2_in_z,arr,decomp)

    ! Finished
end subroutine

subroutine hook_output(decomp,der,fil,dx,dy,dz,outputdir,mesh,fields,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index, &
                                sxx_index,sxy_index,sxz_index,syy_index,syz_index,szz_index
    use decomp_2d,        only: decomp_info, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y, nrank
    use DerivativesMod,   only: derivatives
    use FiltersMod,       only: filters
    use operators,        only: curl

    use ShockVortEntr_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(derivatives),               intent(in) :: der
    type(filters),                   intent(in) :: fil
    integer, dimension(2),       intent(in)    :: x_bc, y_bc, z_bc
    integer                                     :: outputunit=229

    real(rkind), allocatable, dimension(:,:,:,:) :: curlg
    real(rkind), allocatable, dimension(:,:,:) :: detg

    real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) ::xtmp,xdum

    character(len=clen) :: outputfile, velstr
    integer :: i,j,k

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => fields(:,:,:,g11_index), g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), g23 => fields(:,:,:,g23_index), &
               g31 => fields(:,:,:,g31_index), g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3),       &
               sxx => fields(:,:,:, sxx_index), sxy => fields(:,:,:, sxy_index), sxz => fields(:,:,:, sxz_index), &
               syy => fields(:,:,:, syy_index), syz => fields(:,:,:, syz_index), szz => fields(:,:,:, szz_index)  )


        ! do post-processing stuff
        allocate(curlg(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),9))
        allocate(detg(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)))

        ! compute vorticity
        call curl(decomp, der, u,   v,   w,   vort, x_bc, y_bc, z_bc)

        ! compute curl of g
        call curl(decomp, der, g11, g12, g13, curlg(:,:,:,1:3),-x_bc, y_bc, z_bc)
        call curl(decomp, der, g21, g22, g23, curlg(:,:,:,4:6), x_bc,-y_bc, z_bc)
        call curl(decomp, der, g31, g32, g33, curlg(:,:,:,7:9), x_bc, y_bc,-z_bc)


        write(velstr,'(I3.3)') nrank
        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/vort_"//trim(velstr)//"_", vizcount, ".dat"

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

        if(vizcount==0) then
          write(outputfile,'(2A)') trim(outputdir),"/tec_vort_"//trim(velstr)//".dat"
          open(unit=outputunit, file=trim(outputfile), form='FORMATTED',STATUS='unknown')
          write(outputunit,'(290a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p","g11","g12","g13","g21","g22","g23","g31","g32","g33","sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar","curlg11","curlg12","curlg13","curlg21","curlg22","curlg23","curlg31","curlg32","curlg33","conterr","vort"'
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES27.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(38ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
                                               g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
                                               sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), curlg(i,j,k,1:9), &
                                               rho(i,j,k)-rho_0*detg(i,j,k), vort(i,j,k,3)
          
            end do
           end do
          end do
          close(outputunit)
        else
          write(outputfile,'(2A)') trim(outputdir),"/tec_vort_"//trim(velstr)//".dat"
          open(unit=outputunit, file=trim(outputfile), form='FORMATTED', STATUS='old', ACTION='write', POSITION='append')
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          write(outputunit,'(a)') ' VARSHARELIST=([1, 2, 3]=1)'
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(35ES26.16)') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
                                               g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
                                               sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), curlg(i,j,k,1:9), &
                                               rho(i,j,k)-rho_0*detg(i,j,k), vort(i,j,k,3)
          
            end do
           end do
          end do
          close(outputunit)
        endif

        ! write out statistics
        if(tsim > tstatstart) then
          write(outputfile,'(2A)') trim(outputdir),"/tec_stats_"//trim(velstr)//".dat"
          if(.not. istatwritestart) then
            istatwritestart = .TRUE.
            open(unit=outputunit, file=trim(outputfile), form='FORMATTED',STATUS='unknown')
            write(outputunit,'(100a)') 'VARIABLES="x","z","vort_sq","ke","vort_fluc_sq","tke"'!,"vort_sq_norm","ke_norm","vort_fluc_sq_norm","tke_norm"'
            write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
            write(outputunit,'(a,ES27.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
            do k=1,decomp%ysz(3)
              do i=1,decomp%ysz(1)
                  !print *, i, k
                  write(outputunit,'(38ES26.16)') x(i,1,k), z(i,1,k), stats(i,k,1)/(tavg+1.0D-32), stats(i,k,2)/(tavg+1.0D-32), (stats(i,k,1)-stats(i,k,3)*stats(i,k,3))/(tavg+1.0D-32), stats(i,k,4)/(tavg+1.0d-32)!, stats(i,k,1)/(stats(1,k,1)+1.0D-32), stats(i,k,2)/(stats(1,k,2)+1.0D-32), (stats(i,k,1)-stats(i,k,3)*stats(i,k,3))/(stats(1,k,1)-stats(1,k,3)*stats(1,k,3)+1.0D-32), stats(i,k,4)/stats(1,k,4)
              enddo
            enddo
          else
            open(unit=outputunit, file=trim(outputfile), form='FORMATTED',STATUS='old',ACTION='write',POSITION='append')
            write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
            write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
            write(outputunit,'(a)') ' VARSHARELIST=([1, 2]=1)'
            do k=1,decomp%ysz(3)
              do i=1,decomp%ysz(1)
                  write(outputunit,'(38ES26.16)') stats(i,k,1)/(tavg+1.0D-32), stats(i,k,2)/(tavg+1.0D-32), (stats(i,k,1)-stats(i,k,3)*stats(i,k,3))/(tavg+1.0D-32), stats(i,k,4)/(tavg+1.0d-32)!, stats(i,k,1)/(stats(1,k,1)+1.0D-32), stats(i,k,2)/(stats(1,k,2)+1.0D-32), (stats(i,k,1)-stats(i,k,3)*stats(i,k,3))/(stats(1,k,1)-stats(1,k,3)*stats(1,k,3)+1.0D-32), stats(i,k,4)/stats(1,k,4)
              enddo
            enddo
          endif
          close(outputunit)
          if(nrank==0) write(*,*) "Done writing statistics to "//trim(outputfile)
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

    use ShockVortEntr_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),       intent(in)    :: x_bc, y_bc, z_bc

    integer :: nx, ny, ax, i, j, k
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, dum
    real(rkind) :: dx, xspng, tspng

    nx = decomp%xsz(1); ny = decomp%ysz(2)
    ax = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => fields(:,:,:,g11_index), g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), g23 => fields(:,:,:,g23_index), &
               g31 => fields(:,:,:,g31_index), g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        ! left x
        if(decomp%yst(1)==1) then
         if(x_bc(1)==0) then
            do k=1,decomp%ysz(3)
             do j=1,decomp%ysz(2)
               !rho(1,j,k) = rho1 * exp(-Ae*sin(k1wv*u1*tsim))
               rho(1,j,k) = rho1 + rho1 * Ae *          cos(k2wv*y(1,j,k) - k1wv*u1*tsim)
                 u(1,j,k) =   u1 +   u1 * Av * sinpsi * cos(k2wv*y(1,j,k) - k1wv*u1*tsim)
                 v(1,j,k) =      -   u1 * Av * cospsi * cos(k2wv*y(1,j,k) - k1wv*u1*tsim)
                 p(1,:,:) = p1
             end do 
            end do 

            g11( 1,:,:) = g11_1;  g12( 1,:,:) = zero; g13( 1,:,:) = zero
            g21( 1,:,:) = zero;   g22( 1,:,:) = one;  g23( 1,:,:) = zero
            g31( 1,:,:) = zero;   g32( 1,:,:) = zero; g33( 1,:,:) = one

         endif
        endif
        
        ! right x
        xspng = 8.0*pi+0.5_rkind
        tspng = 0.5_rkind
        dx = x(2,1,1)-x(1,1,1)
        dum = half*(one + tanh((x-xspng)/tspng))
        if(decomp%yen(1)==nx) then
         if(x_bc(2)==0) then
            rho(ax,:,:) = rho(ax-1,:,:)
            u  (ax,:,:) = u(ax-1,:,:)
            v  (ax,:,:) = zero
            w  (ax,:,:) = zero
            p  (ax,:,:) = p(ax-1,:,:)
            
            g11(ax,:,:) = g11_2;  g12(ax,:,:) = zero; g13(ax,:,:) = zero
            g21(ax,:,:) = zero;   g22(ax,:,:) = one;  g23(ax,:,:) = zero
            g31(ax,:,:) = zero;   g32(ax,:,:) = zero; g33(ax,:,:) = one

         endif
        endif
        

        do i=1,4
            tmp = u
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            u = u + dum*(tmp - u)

            tmp = v
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            v = v + dum*(tmp - v)

            tmp = w
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            w = w + dum*(tmp - w)

            tmp = e
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            e = e + dum*(tmp - e)

            tmp = rho
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            rho = rho + dum*(tmp - rho)

            tmp = p
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            p = p + dum*(tmp - p)

            tmp = g11(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g11(:,:,:) = g11(:,:,:) + dum*(tmp - g11(:,:,:))

            tmp = g12(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g12(:,:,:) = g12(:,:,:) + dum*(tmp - g12(:,:,:))

            tmp = g13(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g13(:,:,:) = g13(:,:,:) + dum*(tmp - g13(:,:,:))

            tmp = g21(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g21(:,:,:) = g21(:,:,:) + dum*(tmp - g21(:,:,:))

            tmp = g22(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g22(:,:,:) = g22(:,:,:) + dum*(tmp - g22(:,:,:))

            tmp = g23(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g23(:,:,:) = g23(:,:,:) + dum*(tmp - g23(:,:,:))

            tmp = g31(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g31(:,:,:) = g31(:,:,:) + dum*(tmp - g31(:,:,:))

            tmp = g32(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g32(:,:,:) = g32(:,:,:) + dum*(tmp - g32(:,:,:))

            tmp = g33(:,:,:)
            call filter(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
            g33(:,:,:) = g33(:,:,:) + dum*(tmp - g33(:,:,:))
        end do

    end associate
end subroutine

subroutine hook_timestep(decomp,der,mesh,fields,step,tsim,dt,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use SolidGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use DerivativesMod,   only: derivatives
    use operators,        only: curl

    use ShockVortEntr_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind),                     intent(in) :: dt
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    integer,     dimension(2),       intent(in) :: x_bc, y_bc, z_bc

    integer :: i, k, nx
    real(rkind) :: sthick, mwa, dx
    real(rkind), dimension(decomp%ysz(1)) :: dpdx

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))

        !! determine shock statistics
        !dx = x(2,1,1) - x(1,1,1)
        !dpdx = 0.0d0
        !nx = size(dpdx,1)
        !!dpdx(2:nx-1) = ( p(3:nx,1,1)-p(1:nx-2,1,1) ) / (two*dx)
        !call der%ddx(p, dpdx, x_bc(1), x_bc(2))
        !sthick = abs(p2-p1)/maxval(dx*abs(dpdx))

        !mwa = maxval(p(:,1,1)-p2)
        !mwa = max(mwa,maxval(p1-p(:,1,1)))
        !mwa = mwa / abs(p2-p1)

        !call message(2,"StatNormalShockess", sthick)
        !call message(2,"Maximum Wiggles Amplitude", mwa)

        

        !k = 1
        if(tsim > tstatstart) then
          !write(*,*) '------'
          !write(*,*) 'vort3: ', maxval(vort(:,:,:,3)), minval(vort(:,:,:,3))
          !write(*,*) 'u    : ', maxval(u(:,:,:)), minval(u(:,:,:))
          !write(*,*) 'v    : ', maxval(v(:,:,:)), minval(v(:,:,:))
          !write(*,*) 'w    : ', maxval(w(:,:,:)), minval(w(:,:,:))
          !write(*,*) 'ke   : ', sum(u(:,:,:)*u(:,:,:) + v(:,:,:)*v(:,:,:) + w(:,:,:)*w(:,:,:))
          if(.not. istatstart) then
            if(nrank==0) write(*,*) "Started collecting statistics at t = ", tsim
            istatstart = .TRUE.
          endif
          ! compute vorticity
          call curl(decomp, der, u,   v,   w,   vort, x_bc, y_bc, z_bc)
          do k = 1, decomp%ysz(3)
            do i = 1, decomp%ysz(1)
                stats(i,k,1) = stats(i,k,1) + dt * sum(vort(i,:,k,3)*vort(i,:,k,3))/real(decomp%ysz(2),rkind)
                stats(i,k,2) = stats(i,k,2) + dt * sum(u(i,:,k)*u(i,:,k) + v(i,:,k)*v(i,:,k) + w(i,:,k)*w(i,:,k))/real(decomp%ysz(2),rkind)
                stats(i,k,3) = stats(i,k,3) + dt * sum(vort(i,:,k,3))/real(decomp%ysz(2),rkind)
                stats(i,k,4) = stats(i,k,4) + dt * sum((u(i,:,k)-umean(i,1))**2 + (v(i,:,k)-umean(i,2))**2 + (w(i,:,k)-umean(i,3))**2)/real(decomp%ysz(2),rkind)
            enddo
          enddo
          tavg = tavg + dt
          !write(*,*) 'dtavg: ', dt, tavg
          !write(*,*) 'stat1: ', maxval(stats(:,:,1)), minval(stats(:,:,1))
          !write(*,*) 'stat2: ', maxval(stats(:,:,2)), minval(stats(:,:,2))
          !write(*,*) '------'
        endif

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

    use ShockVortEntr_data

    deallocate(vort,stats,umean)

end subroutine

