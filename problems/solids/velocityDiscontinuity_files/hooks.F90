!module velocityDiscontinuity_data

    !real(rkind) :: rho_0 = one, tau_0 = 1d-10, u_L = 0d0, u_R = 0d0, v_L = 0d0, v_R =0d0
    !real(rkind) :: rho_L =one, rho_R = one, p_L = one, p_R = one
    !real(rkind) :: ge11_L, ge11_R, ge22_L, ge22_R
    !logical     :: ignore_gij = .false.

!contains

!!!!!!!!!!!!!!!!!!!!!!!!!!

module velocityDiscontinuity_data
    use kind_parameters,  only: rkind
    use constants,        only: one,two,eight,three,six, zero
    use FiltersMod,       only: filters
    implicit none

    real(rkind) :: p_infty = one, Rgas = one, gamma = 1.4_rkind, mu = 10._rkind

    !real(rkind) :: rho_0 = one, p_0 = 0.1_rkind, g_0 = 1d0, tau_0 = 1d-10, u_0 = 0d0, v_0 = 0d0
    real(rkind) :: rho_0 = one, tau_0 = 1d-10, u_L = 0d0, u_R = 0d0, v_L = 0d0, v_R =0d0
    real(rkind) :: rho_L =one, rho_R = one, p_L = one, p_R = one
    real(rkind) :: ge11_L=one, ge11_R=one, ge22_L=one, ge22_R=one
    real(rkind) :: ge21_L=zero, ge21_R=zero, ge12_L=zero, ge12_R=zero
    real(rkind) :: gp11_L=one, gp11_R=one, gp22_L=one, gp22_R=one

    real(rkind) :: interface_init = 0.5, thick=0.01
    real(rkind) :: yield = one 
    real(rkind) :: kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e
    real(rkind) :: melt_t = one, melt_c = one
    integer     :: kos_sh
    logical     :: explPlast = .true.
    logical     :: plastic = .true.
    real(rkind) :: Ly = one, Lx = 1.0d0 
    real(rkind) :: eta_det_ge = one, eta_det_gp = one, eta_det_gt = one, diff_c_ge = one
    real(rkind) :: diff_c_gp = one, diff_c_gt = one
    
    !logical     :: ignore_gij = .false.

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

    use velocityDiscontinuity_data

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
    
    use velocityDiscontinuity_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(inout) :: tstop, dt, tviz
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit, ind
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp_01, dum
    real(rkind) :: a_L, g11_L, g12_L, g13_L, g21_L, g22_L, g23_L, g31_L, g32_L, g33_L, detg_L, detgp_L
    real(rkind) :: a_R, g11_R, g12_R, g13_R, g21_R, g22_R, g23_R, g31_R, g32_R, g33_R, detg_R, detgp_R
    real(rkind) ::             gp12_L, gp13_L, gp21_L,     gp23_L, gp31_L, gp32_L, gp33_L
    real(rkind) ::             gp12_R, gp13_R, gp21_R,     gp23_R, gp31_R, gp32_R, gp33_R

    integer :: nx,ny,nz
    nx = size(mesh,1); ny = size(mesh,2); nz = size(mesh,3)

    namelist /PROBINPUT/  p_infty, Rgas, gamma, mu, rho_0, tau_0, &
                          plastic, explPlast, yield, &
                          interface_init, thick, &
                          p_L, p_R, rho_L, rho_R, ge11_L, ge22_L, ge11_R, ge22_R, & 
                          u_L, U_R, v_L, v_R, ge21_L, ge21_R, ge12_L, ge12_R, gp11_L, gp11_R, gp22_L, gp22_R, &
                          melt_t, melt_c, kos_b, kos_t, kos_h, kos_g, kos_m, &
                          kos_q, kos_f, kos_alpha, kos_beta, kos_e, kos_sh, &    
                          eta_det_ge, eta_det_gp, eta_det_gt, diff_c_ge, &
                          diff_c_gp, diff_c_gt 

    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    !Check # of species
    if (mix%ns /= 1) then
        call GracefulExit("Number of species must be 1 for this problem. Check the input file.",928)
    end if

    ! Initialize mygfil
    call mygfil%init(                        decomp, &
                     .false.,     .TRUE.,    .TRUE., &
                  "gaussian", "gaussian", "gaussian" )

    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index), &
                 rho => fields(:,:,:,rho_index), x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        ! Set materials (including gamma, Pinf, gamma, , rho0, mu, yieldm and tau0)
        !call mix%set_material(1,stiffgas(gamma ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield, tau_0))
        call mix%set_material(1,stiffgas(gamma  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield,tau_0,eta_det_ge,eta_det_gp,eta_det_gt,diff_c_ge,diff_c_gp,diff_c_gt,melt_t,melt_c,kos_b,kos_t,kos_h,kos_g,kos_m,kos_q,kos_f,kos_alpha,kos_beta,kos_e,kos_sh,nx,ny,nz)) !mca: see Sep1SolidEOS.F90 "init"

        ! set logicals for plasticity
        mix%material(1)%plast = plastic; mix%material(1)%explPlast = explPlast

        ! Set logical for whether or not to ignore g_ij (useful for fluids only
        ! problems)
        !mix%ignore_gij = ignore_gij

        ! speed of sound
        a_L = sqrt((gamma*(p_L+p_infty) + 4.0d0/3.0d0*mu)/rho_L)
        a_R = sqrt((gamma*(p_R+p_infty) + 4.0d0/3.0d0*mu)/rho_R)
        
        ! write material properties
        if (nrank == 0) then
            print *, '---Simulating Impact ---'
            write(*,'(3(a,e12.5))') 'rho_L = ', rho_L, ', gam  = ', gamma, ', p_L = ', p_L
            write(*,'(3(a,e12.5))') 'rho_R = ', rho_R, ', gam  = ', gamma, ', p_R = ', p_R
            write(*,'(3(a,e12.5))') 'mu    = ', mu,    ', Rgas = ', Rgas,  ', rho_0 = ', rho_0
            write(*,'(3(a,e12.5))') 'yield = ', yield, ', tau_0 = ', tau_0,  ', p_infty = ', p_infty
        end if

        !Set mixture velocity
        tmp_01 = 0.5d0 * (erf( (x-(interface_init))/(thick*dx) ) + 1.0d0) !goes from 0 -> 1
        u   = ( u_R - u_L ) * tmp_01 + u_L
        v   = ( v_R - v_L ) * tmp_01 + v_L
        w   = zero

        !Set gij tensor
        g11_L = ge11_L;  g12_L = ge12_L; g13_L = zero
        g21_L = ge21_L;  g22_L = ge22_L; g23_L = zero
        g31_L = zero;    g32_L = zero;  
        
        g11_R = ge11_R;  g12_R = ge12_R;   g13_R = zero
        g21_R = ge21_R;  g22_R = ge22_R; g23_R = zero
        g31_R = zero;    g32_R = zero;  

        !Set gij tensor
        gp11_L = gp11_L;  gp12_L = zero;   gp13_L = zero
        gp21_L = zero;    gp22_L = gp22_L; gp23_L = zero
        gp31_L = zero;    gp32_L = zero;  
        
        gp11_R = gp11_R;  gp12_R = zero;   gp13_R = zero
        gp21_R = zero;    gp22_R = gp22_R; gp23_R = zero
        gp31_R = zero;    gp32_R = zero;  


        !set mixture pressure
        mix%material(1)%p = ( p_R - p_L ) * tmp_01 + p_L
        
        ! Make rho compatible with det(g) and rho0 by adjusting ge33 and gp33
        g33_L = (rho_L/rho_0 - g13_L*(g21_L*g32_L-g31_L*g22_L) + (g11_L*g23_L*g32_L - g12_L*g31_L*g23_L)) / (g11_L*g22_L - g12_L*g21_L)
        g33_R = (rho_R/rho_0 - g13_R*(g21_R*g32_R-g31_R*g22_R) + (g11_R*g23_R*g32_R - g12_R*g31_R*g23_R)) / (g11_R*g22_R - g12_R*g21_R)

        gp33_L = (1.0d0 - gp13_L*(gp21_L*gp32_L-gp31_L*gp22_L) + (gp11_L*gp23_L*gp32_L - gp12_L*gp31_L*gp23_L)) / (gp11_L*gp22_L - gp12_L*gp21_L)
        gp33_R = (1.0d0 - gp13_R*(gp21_R*gp32_R-gp31_R*gp22_R) + (gp11_R*gp23_R*gp32_R - gp12_R*gp31_R*gp23_R)) / (gp11_R*gp22_R - gp12_R*gp21_R)

        detg_L  = g11_L*(g22_L*g33_L-g23_L*g32_L) - g12_L*(g21_L*g33_L-g31_L*g23_L) + g13_L*(g21_L*g32_L-g31_L*g22_L)
        detg_R  = g11_R*(g22_R*g33_R-g23_R*g32_R) - g12_R*(g21_R*g33_R-g31_R*g23_R) + g13_R*(g21_R*g32_R-g31_R*g22_R)
        detgp_L  = gp11_L*(gp22_L*gp33_L-gp23_L*gp32_L) - gp12_L*(gp21_L*gp33_L-gp31_L*gp23_L) + gp13_L*(gp21_L*gp32_L-gp31_L*gp22_L)
        detgp_R  = gp11_R*(gp22_R*gp33_R-gp23_R*gp32_R) - gp12_R*(gp21_R*gp33_R-gp31_R*gp23_R) + gp13_R*(gp21_R*gp32_R-gp31_R*gp22_R)

        if (abs(rho_0 * detg_L - rho_L) > 1.0d-6 ) then
            print *, 'rho_L = ', rho_L, ', det_L = ', rho_0*detg_L
            call GracefulExit("Determinant of ge_L is not compatible with rho_L and rho_0. Please Double-check.",928)
        endif
        if (abs(detgp_L - one) > 1.0d-6 ) then
            print *, 'detgp_L = ', detgp_L
            call GracefulExit("Determinant of gp_L is not equal to one. Please Double-check.",928)
        endif
        if (abs(rho_0 * detg_R - rho_R) > 1.0d-6 ) then
            print *, 'rho_R = ', rho_R, ', det_R = ', rho_0*detg_R
            call GracefulExit("Determinant of ge_R is not compatible with rho_R and rho_0. Please Double-check.",928)
        endif
        if (abs(detgp_R - one) > 1.0d-6 ) then
            print *, 'detgp_R = ', detgp_R
            call GracefulExit("Determinant of gp_R is not equal to one. Please Double-check.",928)
        endif

        rho = ( rho_R - rho_L ) * tmp_01 + rho_L

        ! Set initial values of ge (called g) and gp (inverse deformation gradient)
        mix%material(1)%g11 = ( g11_R - g11_L ) * tmp_01 + g11_L; mix%material(1)%g12 = ( g12_R - g12_L ) * tmp_01 + g12_L; mix%material(1)%g13 = ( g13_R - g13_L ) * tmp_01 + g13_L
        mix%material(1)%g21 = ( g21_R - g21_L ) * tmp_01 + g21_L; mix%material(1)%g22 = ( g22_R - g22_L ) * tmp_01 + g22_L; mix%material(1)%g23 = ( g23_R - g23_L ) * tmp_01 + g23_L
        mix%material(1)%g31 = ( g31_R - g31_L ) * tmp_01 + g31_L; mix%material(1)%g32 = ( g32_R - g32_L ) * tmp_01 + g32_L; mix%material(1)%g33 = ( g33_R - g33_L ) * tmp_01 + g33_L

        IF (plastic) THEN
          mix%material(1)%gp11 = ( gp11_R - gp11_L ) * tmp_01 + gp11_L
          mix%material(1)%gp12 = zero
          mix%material(1)%gp13 = zero
          mix%material(1)%gp21 = zero 
          mix%material(1)%gp22 = ( gp22_R - gp22_L ) * tmp_01 + gp22_L
          mix%material(1)%gp23 = zero
          mix%material(1)%gp31 = zero 
          mix%material(1)%gp32 = zero
          mix%material(1)%gp33 = ( gp22_R - gp22_L ) * tmp_01 + gp22_L

          !This is assuming that gp is diagonal
          mix%material(1)%gt11 = mix%material(1)%gp11 * mix%material(1)%g11  
          mix%material(1)%gt12 = mix%material(1)%gp22 * mix%material(1)%g12 
          mix%material(1)%gt13 = mix%material(1)%gp33 * mix%material(1)%g13 
          mix%material(1)%gt21 = mix%material(1)%gp11 * mix%material(1)%g21 
          mix%material(1)%gt22 = mix%material(1)%gp22 * mix%material(1)%g22 
          mix%material(1)%gt23 = mix%material(1)%gp33 * mix%material(1)%g23 
          mix%material(1)%gt31 = mix%material(1)%gp11 * mix%material(1)%g31 
          mix%material(1)%gt32 = mix%material(1)%gp22 * mix%material(1)%g32 
          mix%material(1)%gt33 = mix%material(1)%gp33 * mix%material(1)%g33 
        ELSE
          !deformation elastic only
          mix%material(1)%gp11 = one;  mix%material(1)%gp12 = zero; mix%material(1)%gp13 = zero
          mix%material(1)%gp21 = zero; mix%material(1)%gp22 = one;  mix%material(1)%gp23 = zero
          mix%material(1)%gp31 = zero; mix%material(1)%gp32 = zero; mix%material(1)%gp33 = one

          !gt should be same as g
          mix%material(1)%gt11 = mix%material(1)%g11
          mix%material(1)%gt12 = mix%material(1)%g12
          mix%material(1)%gt13 = mix%material(1)%g13
          mix%material(1)%gt21 = mix%material(1)%g21
          mix%material(1)%gt22 = mix%material(1)%g22
          mix%material(1)%gt23 = mix%material(1)%g23
          mix%material(1)%gt31 = mix%material(1)%g31
          mix%material(1)%gt32 = mix%material(1)%g32
          mix%material(1)%gt33 = mix%material(1)%g33
        ENDIF

        if ( (mix%use_gTg) .AND. (kos_sh .eq. 2) ) then
            !assuming ge-Gp formulation
            mix%material(1)%gp11 = mix%material(1)%gp11**2;  
            mix%material(1)%gp22 = mix%material(1)%gp22**2;  
            mix%material(1)%gp33 = mix%material(1)%gp33**2;  
        else
            !call GracefulExit("This test problems assume ge-Gp formulation, which required use_gTg=.TRUE. and kos_sh=2",928)
        endif

        !set mixture Volume fraction and Mass Fraction
        mix%material(1)%VF = one
        mix%material(1)%Ys = one;
        
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

    use velocityDiscontinuity_data

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
    real(rkind) :: vort_pos, vort_neg, Al_mass, xspike, xbubbl, xspike_proc, xbubbl_proc
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

        
       write(str,'(I4.4,A,ES7.1E2,A,ES7.1E2)') nrank, "_", 0d0, "_", 1d0
       
       if (mix%use_gTg) then
           str = trim(str)//'_gTg'
       else
           str = trim(str)//'_g'
       end if

       if (decomp%ysz(2) == 1) then
           write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/velocityDiscontinuity_"//trim(str)//"_", vizcount, ".dat"

           open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
           write(outputunit,'(4ES27.16E3)') tsim, 0d0, 0d0, 1d0
           do i=1,decomp%ysz(1)
               write(outputunit,'(15ES27.16E3)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                              mix%material(1)%p (i,1,1),  &
                                              mix%material(1)%Ys(i,1,1),  &
                                              mix%material(1)%VF(i,1,1),  &
                                              mix%material(1)%eh(i,1,1),  &
                                              mix%material(1)%T (i,1,1),  &
                                              mix%material(1)%g11(i,1,1), &
                                              mu(i,1,1), bulk(i,1,1), mix%material(1)%kap(i,1,1), &
                                              mix%material(1)%diff(i,1,1)
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

       Al_mass = 0d0

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
       write(outputunit,'(7ES27.16E3)') tsim, 0d0, vort_pos, vort_neg, Al_mass, xspike, xbubbl
       close(outputunit)

        write(outputfile,'(4A)') trim(outputdir),"/tec_MultSpecShock_"//trim(str),".dat"
        if(vizcount==0) then
          open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='replace')
          write(outputunit,'(350a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p", &
                                     "sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar", &
                                     "p-1","Ys-1","VF-1","eh-1","T-1","g11-1","g12-1","g13-1","g21-1","g22-1","g23-1","g31-1","g32-1","g33-1","Dstar-1","kap-1","rhom-1"'
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                write(outputunit,'(35ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &                      ! continuum (9)
                                                sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), &      ! continuum (9)
                                                mix%material(1)% p(i,j,k),  mix%material(1)% Ys(i,j,k), mix%material(1)% VF(i,j,k), mix%material(1)% eh(i,j,k), &  ! material 1 (14)
                                                mix%material(1)% T(i,j,k),  mix%material(1)%g11(i,j,k), mix%material(1)%g12(i,j,k), mix%material(1)%g13(i,j,k), &  ! material 1 
                                                mix%material(1)%g21(i,j,k), mix%material(1)%g22(i,j,k), mix%material(1)%g23(i,j,k), mix%material(1)%g31(i,j,k), &  ! material 1 
                                                mix%material(1)%g32(i,j,k), mix%material(1)%g33(i,j,k), mix%material(1)%diff(i,j,k), mix%material(1)%kap(i,j,k),& 
                                                mix%material(1)%rhom(i,j,k) 
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
                write(outputunit,'(32ES26.16)') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &                                                    ! continuum (6)
                                                sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k), &      ! continuum (9)
                                                mix%material(1)% p(i,j,k),  mix%material(1)% Ys(i,j,k), mix%material(1)% VF(i,j,k), mix%material(1)% eh(i,j,k), &  ! material 1 (14)
                                                mix%material(1)% T(i,j,k),  mix%material(1)%g11(i,j,k), mix%material(1)%g12(i,j,k), mix%material(1)%g13(i,j,k), &  ! material 1 
                                                mix%material(1)%g21(i,j,k), mix%material(1)%g22(i,j,k), mix%material(1)%g23(i,j,k), mix%material(1)%g31(i,j,k), &  ! material 1 
                                                mix%material(1)%g32(i,j,k), mix%material(1)%g33(i,j,k), mix%material(1)%diff(i,j,k), mix%material(1)%kap(i,j,k),&
                                                mix%material(1)%rhom(i,j,k)  ! material 1 
            end do
           end do
          end do
          close(outputunit)
        endif

    end associate

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Adapted from Multispecies_contact/hooks.f90!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: filter3D

    use velocityDiscontinuity_data

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

      !overwrite plastic entropy equation
      mix%material(1)%pe = zero

        if(decomp%yst(1)==1) then
          if(x_bc(1)==0) then
              rho( 1,:,:) = rho_L
              u  ( 1,:,:) = u_L
              v  ( 1,:,:) = v_L
              w  ( 1,:,:) = zero
              mix%material(1)%p  (1,:,:) = p_L ! mix%material(1)%p(nx-1,:,:)
              
              mix%material(1)%g11( 1,:,:)  = ge11_L; mix%material(1)%g12( 1,:,:) = ge12_L; mix%material(1)%g13( 1,:,:) = zero
              mix%material(1)%g21( 1,:,:)  = ge21_L; mix%material(1)%g22( 1,:,:) = ge22_L; mix%material(1)%g23( 1,:,:) = zero
              mix%material(1)%g31( 1,:,:)  = zero;   mix%material(1)%g32( 1,:,:) = zero;   mix%material(1)%g33( 1,:,:) = ge22_L

              mix%material(1)%gp11( 1,:,:)  = gp11_L**2; mix%material(1)%gp12( 1,:,:) = zero;      mix%material(1)%gp13( 1,:,:) = zero
              mix%material(1)%gp21( 1,:,:)  = zero;      mix%material(1)%gp22( 1,:,:) = gp22_L**2; mix%material(1)%gp23( 1,:,:) = zero
              mix%material(1)%gp31( 1,:,:)  = zero;      mix%material(1)%gp32( 1,:,:) = zero;      mix%material(1)%gp33( 1,:,:) = gp22_L**2

              mix%material(1)%gt11( 1,:,:) = ge11_L*gp11_L; mix%material(1)%gt12( 1,:,:) = ge12_L*gp22_L; mix%material(1)%gt13( 1,:,:) = zero
              mix%material(1)%gt21( 1,:,:) = ge21_L*gp11_L; mix%material(1)%gt22( 1,:,:) = ge22_L*gp22_L; mix%material(1)%gt23( 1,:,:) = zero
              mix%material(1)%gt31( 1,:,:) = zero;          mix%material(1)%gt32( 1,:,:) = zero;          mix%material(1)%gt33( 1,:,:) = ge22_L*gp22_L
  
              mix%material(1)%Ys ( 1,:,:) = one
  
              mix%material(1)%VF ( 1,:,:) = one       !Dirichlet BC: set left volume fraction to VFL
          end if
        endif
        
        if(decomp%yen(1)==decomp%xsz(1)) then
          if(x_bc(2)==0) then
            rho(nx,:,:) = rho_R ! rho(nx-1,:,:)
            u  (nx,:,:) = u_R ! zero
            v  (nx,:,:) = v_R ! v(nx-1,:,:)
            w  (nx,:,:) = zero ! w(nx-1,:,:)
            mix%material(1)%p  (nx,:,:) = p_R ! mix%material(1)%p(nx-1,:,:)
            
            mix%material(1)%g11(nx,:,:) = ge11_R;  mix%material(1)%g12(nx,:,:) = ge12_R;   mix%material(1)%g13(nx,:,:) = zero
            mix%material(1)%g21(nx,:,:) = ge21_R;  mix%material(1)%g22(nx,:,:) = ge22_R; mix%material(1)%g23(nx,:,:) = zero
            mix%material(1)%g31(nx,:,:) = zero;    mix%material(1)%g32(nx,:,:) = zero;   mix%material(1)%g33(nx,:,:) = ge22_R
  
            mix%material(1)%gp11(nx,:,:) = gp11_R**2; mix%material(1)%gp12(nx,:,:) = zero;      mix%material(1)%gp13(nx,:,:) = zero
            mix%material(1)%gp21(nx,:,:) = zero;      mix%material(1)%gp22(nx,:,:) = gp22_R**2; mix%material(1)%gp23(nx,:,:) = zero
            mix%material(1)%gp31(nx,:,:) = zero;      mix%material(1)%gp32(nx,:,:) = zero;      mix%material(1)%gp33(nx,:,:) = gp22_R**2
  
            mix%material(1)%gt11(nx,:,:) = ge11_R*gp11_R; mix%material(1)%gt12(nx,:,:) = ge12_R*gp22_R; mix%material(1)%gt13(nx,:,:) = zero
            mix%material(1)%gt21(nx,:,:) = ge21_R*gp11_R; mix%material(1)%gt22(nx,:,:) = ge22_R*gp22_R; mix%material(1)%gt23(nx,:,:) = zero
            mix%material(1)%gt31(nx,:,:) = zero;          mix%material(1)%gt32(nx,:,:) = zero;          mix%material(1)%gt33(nx,:,:) = ge22_R*gp22_R
  
            mix%material(1)%Ys (nx,:,:) = one
  
            mix%material(1)%VF (nx,:,:) = one
          endif
        endif

        ! apply sponge at left and right boundaries to damp outgoing waves
        xspngL = 0.15d0*Lx
        xspngR = 0.85d0*Lx
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

            do j = 1,9
                tmp = mix%material(1)%g(:,:,:,j)
                call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                mix%material(1)%g(:,:,:,j) = mix%material(1)%g(:,:,:,j) + dum*(tmp - mix%material(1)%g(:,:,:,j))

                tmp = mix%material(1)%g_t(:,:,:,j)
                call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                mix%material(1)%g_t(:,:,:,j) = mix%material(1)%g_t(:,:,:,j) + dum*(tmp - mix%material(1)%g_t(:,:,:,j))

                tmp = mix%material(1)%g_p(:,:,:,j)
                call filter3D(decomp,mygfil,tmp,1,x_bc,y_bc,z_bc)
                mix%material(1)%g_p(:,:,:,j) = mix%material(1)%g_p(:,:,:,j) + dum*(tmp - mix%material(1)%g_p(:,:,:,j))

            end do
        end do


    end associate
end subroutine

!!!!Adapted from Multispecies_contact/hooks.f90!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use constants,        only: zero, one
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use SolidMixtureMod,  only: solid_mixture

    use velocityDiscontinuity_data

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
        ! if (mix%ignore_gij) then
        !     !enforce gij = delta_ij (applies regardless of g or gTg formulation)
        !     do i=1, mix%ns
        !         mix%material(i)%g11 = one;  mix%material(i)%g12 = zero; mix%material(i)%g13 = zero
        !         mix%material(i)%g21 = zero; mix%material(i)%g22 = one;  mix%material(i)%g23 = zero
        !         mix%material(i)%g31 = zero; mix%material(i)%g32 = zero; mix%material(i)%g33 = one
        !     enddo
        ! endif

    end associate
end subroutine

subroutine hook_mixture_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,&
                                mom_index,TE_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use velocityDiscontinuity_data

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

    use velocityDiscontinuity_data

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

    use velocityDiscontinuity_data

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

    use velocityDiscontinuity_data

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

    use velocityDiscontinuity_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
