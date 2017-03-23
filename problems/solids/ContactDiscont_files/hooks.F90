module ContactDiscont_data
    use kind_parameters,  only: rkind
    use constants,        only: one,third,half,twothird,two,three,four,seven
    implicit none
    
    real(rkind) :: kparam = real(0.1D0,rkind)
    real(rkind) :: thick = one
    real(rkind) :: p1 = real(1.D5,rkind), p2, rho1, rho2, u1, u2, g11_1, g11_2
    real(rkind) :: rho_0
    real(rkind) :: Ckap

    real(rkind), parameter :: eleventhird = real(11.D0/3.D0,rkind), seventhird = real(7.D0/3.D0,rkind)

contains

SUBROUTINE fnumden(rho1,fparams,iparams,num,den)

  IMPLICIT NONE
  REAL(rkind), INTENT(IN) :: rho1
  REAL(rkind), INTENT(IN), DIMENSION(:) :: fparams
  INTEGER, INTENT(IN), DIMENSION(:) :: iparams
  REAL(rkind), INTENT(OUT) :: num,den

  REAL(rkind) :: g11, grho1, gprho1, grho2, fac, p2, rho0, gam, pinf, mus, onebygam

  real(rkind), parameter :: eightthird = real(8.D0/3.D0,rkind), fourthird = real(4.D0/3.D0,rkind), &
                            eleven = 11.0_rkind

    write(*,*) '-1-'
    rho0 = fparams(1); gam = fparams(2); pinf = fparams(3); mus = fparams(4);
    rho2 = fparams(5); p2  = fparams(6); kparam = fparams(7)

    onebygam = one/gam
    write(*,*) '-2-', rho0

    g11 = rho1/rho0    ! g11_1
    write(*,*) '--', rho1
    grho1 = twothird*(g11**eleventhird - g11**(-third) - g11**seventhird + g11**third)
    write(*,*) '--', g11, grho1
    gprho1 = third*twothird/rho0*(eleven*g11**eightthird - seven*g11**fourthird + g11**(-fourthird) + g11**(-twothird))
    write(*,*) '-3-'

    g11 = rho2/rho0     ! g11_2
    grho2 = twothird*(g11**eleventhird - g11**(-third) - g11**seventhird + g11**third)
    write(*,*) '-4-'

    fac = (one + mus*(grho2 - grho1)/(p2 + pinf))
    write(*,*) '-5-'
    
    num =  fac**onebygam*rho2/rho1 - kparam
    den = -fac**onebygam*rho2/rho1 * (one/rho1 + mus*onebygam*gprho1/((p2 + pinf)*fac))

    write(*,'(a,8(e19.12,1x))') '--', rho0, rho1, rho2, grho1, grho2, fac, num, den

END SUBROUTINE fnumden

SUBROUTINE rootfind_nr_1d(pf,fparams,iparams)

  IMPLICIT NONE
  REAL(rkind), INTENT(INOUT) :: pf
  REAL(rkind), INTENT(IN), DIMENSION(:) :: fparams
  INTEGER, INTENT(IN), DIMENSION(:) :: iparams

  INTEGER :: ii, itmax = 1000
  REAL(rkind) :: tol = 1.0d-12
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
    if(dabs(pf)>tol) then
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

    use ContactDiscont_data

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

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,eostype,eosparams,rho0,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,eps,third,half,one,two,pi
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,&
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info
    
    use ContactDiscont_data

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
    integer, dimension(2) :: iparams
    real(rkind), dimension(8) :: fparams
    integer :: nx
    real(rkind) :: mu, gam, PInf, yield, tau0, grho1, grho2

    namelist /PROBINPUT/  kparam, p2, thick, Ckap
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    write(*,*) 'eostype = ', eostype
    write(*,*) 'eosparams = ', eosparams
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

        rho2 = rho0
        
        fparams(1) = rho0; fparams(2) = gam; fparams(3) = PInf; fparams(4) = mu;
        fparams(5) = rho2; fparams(6) = p2;  fparams(7) = kparam;
        rho1 = rho2/kparam ! Init guess
        write(*,*) 'Before root finding: ', fparams(1:7)
        call rootfind_nr_1d(rho1,fparams,iparams)
        write(*,*) 'After root finding: ', rho1

        grho1 = twothird*((rho1/rho0)**eleventhird - (rho1/rho0)**(-third) - (rho1/rho0)**seventhird + (rho1/rho0)**third)
        grho2 = twothird*((rho2/rho0)**eleventhird - (rho2/rho0)**(-third) - (rho2/rho0)**seventhird + (rho2/rho0)**third)
        p1 = p2 + mu*(grho2-grho1)

        ! initialize contact discontinuity at 0.5
        tmp = half * ( one + erf( (x-half)/(thick*dx) ) )

        u   = zero
        v   = zero
        w   = zero
        rho = (one-tmp)*rho1 + tmp*rho2
        p   = (one-tmp)*  p1 + tmp*  p2

        !rho1 = rho(decomp%yen(1),1,1)
        !u1   = u  (decomp%yen(1),1,1)
        !u2   = u  (            1,1,1)

        g11 = rho/rho0; g12 = zero; g13 = zero
        g21 = zero;     g22 = one;  g23 = zero
        g31 = zero;     g32 = zero; g33 = one

        rho_0 = rho0

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

    use ContactDiscont_data

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

        write(str,'(ES7.1E2,A,ES7.1E2)') kparam, "_", Ckap
        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ContactDiscont_"//trim(str)//"_", vizcount, ".dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        do i=1,decomp%ysz(1)
            write(outputunit,'(10ES27.16E3)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                           g11(i,1,1), g21(i,1,1), mu(i,1,1), bulk(i,1,1), kap(i,1,1)
        
        end do
        close(outputunit)

        ! write(outputfile,'(2A)') trim(outputdir),"/tec_CD_"//trim(str)//".dat"
        ! if(vizcount==0) then

        !   open(unit=outputunit, file=trim(outputfile), form='FORMATTED',STATUS='unknown')
        !   write(outputunit,'(290a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p","g11","g12","g13","g21","g22","g23","g31","g32","g33","sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar"'
        !   write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
        !   write(outputunit,'(a,ES27.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
        !   do k=1,decomp%ysz(3)
        !    do j=1,decomp%ysz(2)
        !     do i=1,decomp%ysz(1)
        !         write(outputunit,'(37ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
        !                                        g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
        !                                        sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k)
        !   
        !     end do
        !    end do
        !   end do
        !   close(outputunit)
        ! else
        !   open(unit=outputunit, file=trim(outputfile), form='FORMATTED',STATUS='old',ACTION='write',POSITION='append')
        !   write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
        !   write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
        !   write(outputunit,'(a)') ' VARSHARELIST=([1, 2, 3]=1)'
        !   do k=1,decomp%ysz(3)
        !    do j=1,decomp%ysz(2)
        !     do i=1,decomp%ysz(1)
        !         write(outputunit,'(34ES26.16)') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), e(i,j,k), p(i,j,k), &
        !                                        g11(i,j,k), g12(i,j,k), g13(i,j,k), g21(i,j,k), g22(i,j,k), g23(i,j,k), g31(i,j,k), g32(i,j,k), g33(i,j,k), &
        !                                        sxx(i,j,k), sxy(i,j,k), sxz(i,j,k), syy(i,j,k), syz(i,j,k), szz(i,j,k), mu(i,j,k), bulk(i,j,k), kap(i,j,k)
        !   
        !     end do
        !    end do
        !   end do
        !   close(outputunit)
        ! endif
    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info

    use ContactDiscont_data

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
        
        g11( 1,:,:) = rho2/rho_0; g12( 1,:,:) = zero; g13( 1,:,:) = zero
        g21( 1,:,:) = zero; g22( 1,:,:) = one;  g23( 1,:,:) = zero
        g31( 1,:,:) = zero; g32( 1,:,:) = zero; g33( 1,:,:) = one

        rho(nx,:,:) = rho1
        u  (nx,:,:) = u1
        v  (nx,:,:) = zero
        w  (nx,:,:) = zero
        p  (nx,:,:) = p1
        
        g11(nx,:,:) = rho1/rho_0; g12(nx,:,:) = zero; g13(nx,:,:) = zero
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

    use ContactDiscont_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind),                     intent(in) :: dt
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    integer,     dimension(2),       intent(in) :: x_bc, y_bc, z_bc

    integer :: nx, istart, iend, i
    real(rkind), dimension(decomp%ysz(1)) :: dedx, drdx
    real(rkind) :: dx, sthick, mwa, e_lo, e_hi, sthick_r, mwa_r, rho_hi, rho_lo
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

        e_hi = max(e(1,1,1), e(nx,1,1))
        e_lo = min(e(1,1,1), e(nx,1,1))
        

        dedx = zero
        dedx(2:nx-1) = ( e(3:nx,1,1)-e(1:nx-2,1,1) ) / (two*dx)
        sthick = (e_hi-e_lo)/max(maxval(dx*abs(dedx)), 1.0d-32)

        mwa = maxval(e(:,1,1)-e_hi)
        mwa = max(mwa, maxval(e_lo-e(:,1,1)))
        mwa = mwa / max(e_hi-e_lo, 1.0D-32)

        call message(2,"Shock thickness", sthick)
        call message(2,"Maximum Wiggles Amplitude", mwa)

        ! ----- rho oscillations -----

          rho_hi = max(rho(1,1,1), rho(nx,1,1))
          rho_lo = min(rho(1,1,1), rho(nx,1,1))

          drdx = zero
          drdx(2:nx-1) = ( rho(3:nx,1,1)-rho(1:nx-2,1,1) ) / (two*dx)
          sthick_r = abs(rho2-rho1)/max(maxval(dx*abs(drdx)), 1.0d-32)

          mwa_r = maxval(rho(:,1,1)-rho_hi)
          mwa_r = max(mwa_r, maxval(rho_lo-rho(:,1,1)))
          mwa_r = mwa_r / max(rho_hi-rho_lo, 1.0D-32)
        ! ----- rho oscillations -----

        write(outputfile,'(A)') "ContactDiscont_stats.dat"
        if (step == 1) then
            open(unit=iounit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
            write(iounit,'(5A26)') '"Time", "E thickness", "MWA", "R thickness", "R MWA"'
        else
            open(unit=iounit, file=trim(outputfile), form='FORMATTED', position='APPEND', status='OLD')
        end if
        write(iounit,'(5ES26.16)') tsim, sthick, mwa,sthick_r,mwa_r
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

