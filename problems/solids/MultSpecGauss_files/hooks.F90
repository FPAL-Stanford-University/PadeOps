module MultSpecGauss_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight,two
    implicit none

    real(rkind) :: p_infty   = one, Rgas   = one, gamma   = 1.4_rkind, mu   = 10._rkind, rho_0   = one, p_amb = 0.1_rkind
    real(rkind) :: p_infty_2 = one, Rgas_2 = one, gamma_2 = 1.4_rkind, mu_2 = 10._rkind, rho_0_2 = one, p_amb_2 = 0.1_rkind
    real(rkind) :: K0 = one, alp = one, CV = one, T0 = one
    real(rkind) :: K0_2 = one, alp_2 = one, CV_2 = one, T0_2 = one
    real(rkind) :: minVF = 0.2_rkind, thick = one, rho_ratio = two, thickness = 0.001_rkind, yield = 0.1_rkind, yield_2 = 0.1_rkind, interfloc = 0.5_rkind, interfwidth = 0.4_rkind, tau0 = 1.0d-14, tau0_2 = 1.0d-14
    logical     :: plastic = .false., plastic_2 = .false.
    real(rkind) :: rhomax, rhomin, umax, umin, vmax, vmin, wmax, wmin, pmax, pmin, Tmax, Tmin, vfmax, vfmin
    real(rkind) :: Ys1max, Ys1min, eh1max, eh1min, g1max, g1min, Ys2max, Ys2min, eh2max, eh2min, g2max, g2min
    real(rkind) :: sigma_0 = one, lamL, lamL_2, cL, cL_2, cL_mix, cL_mix_2, crefl, ctran

end module

subroutine meshgen(decomp, dx, dy, dz, mesh, xcentered)
    use kind_parameters,  only: rkind
    use constants,        only: one, half, zero
    use decomp_2d,        only: decomp_info

    use MultSpecGauss_data

    implicit none

    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    logical,                         intent(in)    :: xcentered

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    real(rkind) :: xfst

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,1)x[0,1)x[0,1) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = one/real(nx,rkind)
        dy = dx
        dz = dx

        if(xcentered) then
           xfst = half*dx
        else
           xfst = zero
        endif

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1 + i - 1, rkind ) * dx + xfst
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,half,one,two,three,four,pi
    use SolidGrid,        only: u_index,v_index,w_index
    use decomp_2d,        only: decomp_info
    use exits,            only: GracefulExit
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    use SolidMixtureMod,  only: solid_mixture
    
    use MultSpecGauss_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(inout) :: tstop, dt, tviz
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit, indwv, indwv2, iprob = 1
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp, sig11, eps11
    real(rkind) :: p_star, rho_star, c_star
    logical :: SOSmodel

    namelist /PROBINPUT/  p_infty, Rgas, gamma, mu, rho_0, p_amb, thick, minVF, rho_ratio, &
                          p_infty_2, Rgas_2, gamma_2, mu_2, rho_0_2, p_amb_2, sigma_0, SOSmodel, thickness, &
                          plastic, plastic_2, yield, yield_2, tau0, tau0_2, iprob, &
                          K0, alp, CV, T0, K0_2, alp_2, CV_2, T0_2, interfloc, interfwidth
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        if (mix%ns /= 2) then
            call GracefulExit("Number of species must be 2 for this problem. Check the input file.",928)
        end if

        mix%SOSmodel = SOSmodel

        !if(hydroeostype==StiffGasEOS) then
            p_star = gamma*p_infty + (four/three)*mu
            rho_star = rho_0
            c_star = sqrt(p_star/rho_star)
            sigma_0 = sigma_0 / p_star
            tstop = tstop * c_star
            dt = dt * c_star
            tviz = tviz * c_star
            tau0 = tau0 * c_star
            tau0_2 = tau0_2 * c_star
            !dtprob = dt
            !print*, "p_star = ", p_star
            !print*, "rho_star = ", rho_star
            !print*, "c_star = ", c_star
            !print*, "sigma_0 = ", sigma_0
            !print*, "tstop = ", tstop
            !print*, "dt = ", dt
            !print*, "tviz = ", tviz

            ! Non-dimensionalize problem parameters
            rho_0 = rho_0 / rho_star;     rho_0_2 = rho_0_2 / rho_star
            mu = mu / p_star;             mu_2 = mu_2 / p_star
            p_infty = p_infty / p_star;   p_infty_2 = p_infty_2 / p_star

            ! Set materials
            if(rho_ratio > zero) then
              ! if rho_ratio > 0, (both have same gamma here)
              p_amb_2 = p_amb
              p_infty_2 = p_infty!rho_ratio*p_infty + p_amb/gamma*(rho_ratio-one)
              Rgas_2 = Rgas/rho_ratio*(p_amb_2+p_infty_2)/(p_amb+p_infty)
              mu_2 = mu
              rho_0_2 = rho_ratio*rho_0
              gamma_2 = gamma
            else
              Rgas_2 = Rgas*rho_0/rho_0_2*(p_amb_2+p_infty_2)/(p_amb+p_infty)
            endif
            call mix%set_material(1,stiffgas(gamma,  Rgas,  p_infty  ), sep1solid(  rho_0,mu,  yield,  tau0))
            call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2), sep1solid(rho_0_2,mu_2,yield_2,tau0_2))

        !!elseif(hydroeostype==GRhydroeos) then
        !    call mix%set_material(1,stiffgas(rho_0,   K0,   alp,   gamma  , CV,   T0  ), sep1solid(rho_0  ,mu  ,yield,1.0D-14))
        !    call mix%set_material(2,stiffgas(rho_0_2, K0_2, alp_2, gamma_2, CV_2, T0_2), sep1solid(rho_0_2,mu_2,yield_2,1.0D-14))
        !!endif

        mix%material(1)%plast = plastic
        mix%material(2)%plast = plastic_2

        if(thick<zero) then
          tmp = half * ( erf( (x-interfloc)/(thickness) ) - erf( (x-interfloc-interfwidth)/(thickness) ) )
        else
          tmp = half * ( erf( (x-interfloc)/(thick*dx) ) - erf( (x-interfloc-interfwidth)/(thick*dx) ) )
          !tmp = half * ( one + erf( (x-0.5_rkind)/(thick*dx) ) )
        endif
        !write(*,*) 'thickness = ', thickness
        !write(*,*) tmp(48:52,1,1)

        ! Set material 1 properties
        mix%material(1)%g11 = one;  mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero; mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero; mix%material(1)%g32 = zero; mix%material(1)%g33 = one

        mix%material(1)%p  = p_amb
        mix%material(1)%VF = minVF + (one-two*minVF)*tmp
        tmp = rho_0*(mix%material(1)%VF + (one-mix%material(1)%VF)*rho_0_2)      ! mixture density
        mix%material(1)%Ys = mix%material(1)%VF*rho_0/tmp

        ! Set material 2 properties
        mix%material(2)%g11 = one;  mix%material(2)%g12 = zero; mix%material(2)%g13 = zero
        mix%material(2)%g21 = zero; mix%material(2)%g22 = one;  mix%material(2)%g23 = zero
        mix%material(2)%g31 = zero; mix%material(2)%g32 = zero; mix%material(2)%g33 = one

        mix%material(2)%p  = p_amb_2
        mix%material(2)%VF = one - mix%material(1)%VF ! Enforce sum to unity
        mix%material(2)%Ys = one - mix%material(1)%Ys ! Enforce sum to unity

        !if(iprob == 1) then

          ! set linear wave speeds
          lamL = gamma * p_infty - (two/three)*mu ! hydroeos = stiffgas
          !lamL = rho_0*K0 - (two/three)*mu !-gam**2*rho_0*CV*T0 (?) ! GRhydroeos
          cL = sqrt( (lamL + two*mu)/ rho_0 )

          lamL_2 = gamma_2 * p_infty_2 - (two/three)*mu_2 ! hydroeos = stiffgas
          !lamL_2 = rho_0_2*K0_2 - (two/three)*mu_2 !-gam_2**2*rho_0_2*CV_2*T0_2 (?) ! GRhydroeos
          cL_2 = sqrt( (lamL_2 + two*mu_2)/ rho_0_2 )

          ! mixture wave speed
          indwv = nint(0.35_rkind/1.0_rkind*decomp%ysz(1)); write(*,*) 'indwv = ', indwv
          if(mix%SOSmodel) then
              cL_mix = sqrt(one / tmp(indwv,1,1) / (mix%material(1)%VF(indwv,1,1)/rho_0/cL**2 + &
                                                    mix%material(2)%VF(indwv,1,1)/rho_0_2/cL_2**2) )
          else
              cL_mix = sqrt(mix%material(1)%Ys(indwv,1,1)*cL**2 + mix%material(2)%Ys(indwv,1,1)*cL_2**2)
          endif  
          write(*,*) 'eqbmodel = ', mix%SOSmodel
          write(*,*) 'wave speeds', cL, cL_2, cL_mix
          write(*,*) 'Ys', mix%material(1)%Ys(indwv,1,1), mix%material(2)%Ys(indwv,1,1)
          write(*,*) 'VF', mix%material(1)%VF(indwv,1,1), mix%material(2)%VF(indwv,1,1)
          ! mixture wave speed
          indwv2 = nint(0.80_rkind/1.0_rkind*decomp%ysz(1)); write(*,*) 'indwv2 = ', indwv2
          if(mix%SOSmodel) then
              cL_mix_2 = sqrt(one / tmp(indwv2,1,1) / (mix%material(1)%VF(indwv2,1,1)/rho_0/cL**2 + &
                                                    mix%material(2)%VF(indwv2,1,1)/rho_0_2/cL_2**2) )
          else
              cL_mix_2 = sqrt(mix%material(1)%Ys(indwv2,1,1)*cL**2 + mix%material(2)%Ys(indwv2,1,1)*cL_2**2)
          endif  
          ! also compute reflection and transmission coefficients
          crefl = ((tmp(indwv2,1,1)*cL_mix_2)/(tmp(indwv,1,1)*cL_mix) - one) / ((tmp(indwv2,1,1)*cL_mix_2)/(tmp(indwv,1,1)*cL_mix) + one)
          ctran = ((tmp(indwv2,1,1)*cL_mix_2)/(tmp(indwv,1,1)*cL_mix) * two) / ((tmp(indwv2,1,1)*cL_mix_2)/(tmp(indwv,1,1)*cL_mix) + one)
          write(*,*) 'cr, ct', crefl, ctran, ctran-crefl

          ! correct fields for initial wave
          sig11 = -sigma_0 * exp( -((x - cL_mix*zero - 0.35_rkind)/(0.03_rkind))**2 )
          eps11 = sig11 / (lamL + two*mu)

          u   = - (cL_mix * sig11) / (lamL + two*mu)
          v   = zero
          w   = zero

          mix%material(1)%g11 = one / (one + eps11)
          mix%material(1)%p = p_infty*(mix%material(1)%g11**gamma - one) ! hydroeos = stiffgas
          mix%material(1)%T = T0*mix%material(1)%g11**gamma
          !mix%material(1)%p = -(lamL + two/three*mu)*eps11     !GRhydroeos
          !mix%material(1)%p = rho_0*mix%material(1)%g11*(K0/alp*mix%material(1)%g11**alp*(mix%material(1)%g11**alp-one))

          mix%material(2)%g11 = one / (one + eps11)
          mix%material(2)%p = p_infty_2*(mix%material(2)%g11**gamma_2 - one) ! hydroeos = stiffgas
          mix%material(2)%T = T0_2*mix%material(2)%g11**gamma_2
          !mix%material(2)%p = mix%material(1)%p  !-(lamL_2 + two/three*mu_2)*eps11     !GRhydroeos
          !mix%material(2)%p = rho_0_2*mix%material(2)%g11*(K0_2/alp_2*mix%material(2)%g11**alp_2*(mix%material(2)%g11**alp_2-one))

        !else
        if(iprob==2) then

          u = zero

          tmp = 1.0d0/(25.d0)*exp(-((x-0.4d0)/0.0625d0)**2/2.0d0)
          print *, maxval(tmp), minval(tmp)
          print *, '--', tmp(1,1,1), tmp(101,1,1)
          mix%material(1)%g11 = 1.1d0
          print *, 1
          mix%material(1)%g22 = 1.1d0*(1.0d0+9.0d0*tmp)
          print *, 2, maxval(mix%material(1)%g22), minval(mix%material(1)%g22)
          mix%material(1)%g33 = 1.1d0/(1.0d0+9.0d0*tmp)
          print *, 3, maxval(mix%material(1)%g33), minval(mix%material(1)%g33)

          mix%material(2)%g11 = mix%material(1)%g11
          mix%material(2)%g22 = mix%material(1)%g22
          mix%material(2)%g33 = mix%material(1)%g33

          mix%material(1)%p  = p_amb
          mix%material(2)%p  = p_amb_2

        endif

    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,eps,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                sxx_index,syy_index,szz_index,sxy_index,sxz_index,syz_index,sos_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use SolidMixtureMod,  only: solid_mixture

    use MultSpecGauss_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix
    integer, dimension(2), intent(in), optional :: x_bc, y_bc, z_bc
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile, str
    integer :: i, j, k
    real(rkind) :: pinc, prefl, ptra

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

        write(str,'(I5.5)') decomp%ysz(1)
        write(outputfile,'(2A,I5.5,A)') trim(outputdir),"/MultSpecGauss_"//trim(str)//"_", vizcount, ".dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        write(outputunit,'(1ES27.16E3)') tsim
        do i=1,decomp%ysz(1)
            write(outputunit,'(24ES27.16E3)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), sos(i,1,1), &
                                           mix%material(1)%p (i,1,1), mix%material(2)%p (i,1,1), &
                                           mix%material(1)%Ys(i,1,1), mix%material(2)%Ys(i,1,1), &
                                           mix%material(1)%VF(i,1,1), mix%material(2)%VF(i,1,1), &
                                           mix%material(1)%eh(i,1,1), mix%material(2)%eh(i,1,1), &
                                           mix%material(1)%T (i,1,1), mix%material(2)%T (i,1,1), &
                                           mix%material(1)%g11(i,1,1), mix%material(2)%g11(i,1,1), &
                                           mu(i,1,1), bulk(i,1,1), kap(i,1,1), kap(i,1,1), &
                                           mix%material(1)%diff(i,1,1), mix%material(2)%diff(i,1,1)
        end do
        close(outputunit)

        write(outputfile,'(5A)') trim(outputdir),"/tec_MultSpecGauss_"//trim(str),".dat"
        if(vizcount==0) then
          open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='replace')
          write(outputunit,'(330a)') 'VARIABLES="x","y","z","rho","u","v","w","e","p", &
                                     "sig11","sig12","sig13","sig22","sig23","sig33","mustar","betstar","kapstar", &
                                     "p-1","Ys-1","VF-1","eh-1","T-1","g11-1","g12-1","g13-1","g21-1","g22-1","g23-1","g31-1","g32-1","g33-1","Dstar-1",&
                                     "p-2","Ys-2","VF-2","eh-2","T-2","g11-2","g12-2","g13-2","g21-2","g22-2","g23-2","g31-2","g32-2","g33-2","Dstar-2","rhom-1","rhom-2","eel-1","eel-2"'
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
                                                mix%material(1)%g32(i,j,k), mix%material(1)%g33(i,j,k), mix%material(1)%diff(i,j,k), &                             ! material 1 
                                                mix%material(2)% p(i,j,k),  mix%material(2)% Ys(i,j,k), mix%material(2)% VF(i,j,k), mix%material(2)% eh(i,j,k), &  ! material 2 (14)
                                                mix%material(2)% T(i,j,k),  mix%material(2)%g11(i,j,k), mix%material(2)%g12(i,j,k), mix%material(2)%g13(i,j,k), &  ! material 2
                                                mix%material(2)%g21(i,j,k), mix%material(2)%g22(i,j,k), mix%material(2)%g23(i,j,k), mix%material(2)%g31(i,j,k), &  ! material 2
                                                mix%material(2)%g32(i,j,k), mix%material(2)%g33(i,j,k), mix%material(2)%diff(i,j,k), &                             ! material 2
                                                mix%material(1)%rhom(i,j,k), mix%material(2)%rhom(i,j,k), &
                                                mix%material(1)%eel(i,j,k),  mix%material(2)%eel(i,j,k)
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
                                                mix%material(1)%g32(i,j,k), mix%material(1)%g33(i,j,k), mix%material(1)%diff(i,j,k), &                             ! material 1 
                                                mix%material(2)% p(i,j,k),  mix%material(2)% Ys(i,j,k), mix%material(2)% VF(i,j,k), mix%material(2)% eh(i,j,k), &  ! material 2 (14)
                                                mix%material(2)% T(i,j,k),  mix%material(2)%g11(i,j,k), mix%material(2)%g12(i,j,k), mix%material(2)%g13(i,j,k), &  ! material 2
                                                mix%material(2)%g21(i,j,k), mix%material(2)%g22(i,j,k), mix%material(2)%g23(i,j,k), mix%material(2)%g31(i,j,k), &  ! material 2
                                                mix%material(2)%g32(i,j,k), mix%material(2)%g33(i,j,k), mix%material(2)%diff(i,j,k), &                             ! material 2
                                                mix%material(1)%rhom(i,j,k), mix%material(2)%rhom(i,j,k), &
                                                mix%material(1)%eel(i,j,k),  mix%material(2)%eel(i,j,k)
            end do
           end do
          end do
          close(outputunit)
        endif

        write(outputfile,'(5A)') trim(outputdir),"/tec_ExactGauss_"//trim(str),".dat"
        if(vizcount==0) then
          open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='replace')
          write(outputunit,'(330a)') 'VARIABLES="x","y","z","pincident","preflected","ptransmitted"'
          write(outputunit,'(6(a,i7),a)') 'ZONE I=', decomp%ysz(1), ' J=', decomp%ysz(2), ' K=', decomp%ysz(3), ' ZONETYPE=ORDERED'
          write(outputunit,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', tsim
          do k=1,decomp%ysz(3)
           do j=1,decomp%ysz(2)
            do i=1,decomp%ysz(1)
                pinc =  sigma_0 * exp( -((x(i,j,k) - cL_mix*tsim - 0.35_rkind)/(0.03_rkind))**2 )
                prefl = sigma_0 * exp( -((-x(i,j,k) - cL_mix*tsim - 0.35_rkind + half + half)/(0.03_rkind))**2 ) * crefl
                ptra =  sigma_0 * exp( -(cL_mix*(x(i,j,k)/cL_mix_2 - tsim - 0.35_rkind/cL_mix + half/cL_mix - half/cL_mix_2)/(0.03_rkind))**2 ) * ctran
                write(outputunit,'(6ES26.16)') x(i,j,k), y(i,j,k), z(i,j,k), pinc, prefl, ptra
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
                pinc =  sigma_0 * exp( -((x(i,j,k) - cL_mix*tsim - 0.35_rkind)/(0.03_rkind))**2 )
                prefl =  sigma_0 * exp( -((-x(i,j,k) - cL_mix*tsim - 0.35_rkind + half + half)/(0.03_rkind))**2 ) * crefl
                ptra =  sigma_0 * exp( -(cL_mix*(x(i,j,k)/cL_mix_2 - tsim - 0.35_rkind/cL_mix + half/cL_mix - half/cL_mix_2)/(0.03_rkind))**2 ) * ctran
                write(outputunit,'(3ES26.16)') pinc, prefl, ptra
            end do
           end do
          end do
          close(outputunit)
        endif

        if(vizcount==0) then
          ! set max and min for later use
          rhomax = maxval(rho); rhomin = minval(rho); umax = half; umin = half; vmax = zero; vmin = zero; wmax = zero; wmin = zero
          vfmax = one-two*minVF; vfmin = minVF; pmax = max(p_amb_2, p_amb); pmin = min(p_amb_2, p_amb);
          Ys1max = maxval(mix%material(1)%Ys); Ys1min = minval(mix%material(1)%Ys);  g1max = maxval(mix%material(1)%g11); g1min = minval(mix%material(1)%g11); 
          Ys2max = maxval(mix%material(2)%Ys); Ys2min = minval(mix%material(2)%Ys);  g2max = maxval(mix%material(2)%g11); g2min = minval(mix%material(2)%g11);

          Tmax = max(maxval(mix%material(1)%T), maxval(mix%material(1)%T)); Tmin = min(minval(mix%material(1)%T), minval(mix%material(1)%T));
          eh1max = maxval(mix%material(1)%eh); eh1min = minval(mix%material(1)%eh);  eh2max = maxval(mix%material(2)%eh); eh2min = minval(mix%material(2)%eh);

          ! open file
          open(unit=11, file='timeevol.dat', status='replace')
          close(11)
        endif
 
    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use MultSpecGauss_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    integer, dimension(2), intent(in), optional    :: x_bc, y_bc, z_bc

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use SolidMixtureMod,  only: solid_mixture

    use MultSpecGauss_data

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


      if(mod(step,40)==0) then
        open(unit=11, file='timeevol.dat', status='old',action='write',position='append')
        write(11,'(45ES26.16)') tsim, maxval(rho)-rhomax, minval(rho)-rhomin, maxval(u)-umax, minval(u)-umin, maxval(v)-vmax, minval(v)-vmin, & 
                                      maxval(w)-wmax, minval(w)-wmin, maxval(p)-pmax, minval(p)-pmin, maxval(T)-Tmax, minval(T)-Tmin,     &
                                      maxval(mix%material(1)%VF)-vfmax, minval(mix%material(1)%VF)-vfmin, maxval(mix%material(1)%Ys)-Ys1max,  minval(mix%material(1)%Ys)-Ys1min, &
                                      maxval(mix%material(1)%eh)-eh1max, minval(mix%material(1)%eh)-eh1min, maxval(mix%material(1)%g11)-g1max, minval(mix%material(1)%g11)-g1min,&
                                      maxval(mix%material(2)%VF)-vfmax, minval(mix%material(2)%VF)-vfmin, maxval(mix%material(2)%Ys)-Ys2max,  minval(mix%material(2)%Ys)-Ys2min, &
                                      maxval(mix%material(2)%eh)-eh2max, minval(mix%material(2)%eh)-eh2min, maxval(mix%material(2)%g11)-g2max, minval(mix%material(2)%g11)-g2min
        close(11)
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

    use MultSpecGauss_data

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

    use MultSpecGauss_data

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

    use MultSpecGauss_data

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

    use MultSpecGauss_data

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

    use MultSpecGauss_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
