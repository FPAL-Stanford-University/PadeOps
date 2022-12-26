module StiffgasShock_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none
    
    real(rkind) :: pRatio = real(4.5,rkind)
    real(rkind) :: pinfbyp1 = real(1.0D3,rkind)
    real(rkind) :: thick = one
    real(rkind) :: p1, p2, rho1, rho2, u1, u2
    real(rkind) :: Cbeta

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use StiffgasShock_data

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

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    
    use StiffgasShock_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz, tstop, dt, tviz
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    
  !real(rkind),            optional, intent(inout) :: rho0, mu, gam, PInf, tstop, dt, tviz, yield, tau0

    integer :: ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp
    real(rkind) :: rhoRatio, tfactor, h1, gam = 1.4d0, Rgas = one, p_infty = 1.0d0
    real(rkind) :: rho_0 = one, mu = one, yield = one

    namelist /PROBINPUT/  pRatio, pinfbyp1, thick, Cbeta, rhoRatio, &
                          gam, Rgas, p_infty, rho_0, mu, yield
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    associate( rho => fields(:,:,:,rho_index),   u => fields(:,:,:,  u_index), &
                 v => fields(:,:,:,  v_index),   w => fields(:,:,:,  w_index), &
                 p => fields(:,:,:,  p_index),   T => fields(:,:,:,  T_index), &
                 e => fields(:,:,:,  e_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        tmp = half * ( one + erf( (x-half)/(thick*dx) ) )

        !rhoRatio = ( (gam+one)*pRatio + (two*gam)*pinfbyp1 + (gam-one) ) / ( (gam-one)*pRatio + (two*gam)*pinfbyp1 + (gam+one) )

        rho1 = one
        rho2 = rho1 * rhoRatio
        p1 = one !(rho1 / gam) / (one + pinfbyp1)   ! Make speed of sound unity
        p2 = pRatio * p1
        u1 = zero  !-sqrt(rhoRatio * (p2-p1) / (rho2-rho1))
        u2 = zero  !rho1 * u1 / rho2

        !PInf = pinfbyp1 * p1

        !h1 = gam*(p1+PInf)/((gam-one)*rho1) + half*u1*u1

        ! rho = (one-tmp)*rho2 + tmp*rho1
        u   = (one-tmp)*  u2  + tmp*  u1
        p   = (one-tmp)*  p2  + tmp*  p1
        rho = (one-tmp)* rho2 + tmp*  rho1

        v   = zero
        w   = zero

        ! p = ( h1 - half*u*u ) * ((gam-one)*rho) / gam - PInf
        ! print*, "p diff: ", (p2 - p(1,1,1))/p2
        ! p2 = p(1,1,1)
        ! print*, "p diff: ", (p1 - p(decomp%ysz(1),1,1))/p1
        ! p1 = p(decomp%ysz(1),1,1)

        ! Set materials
        call mix%set_material(1,stiffgas(gam  ,Rgas  ,p_infty  ),sep1solid(rho_0  ,mu  ,yield, 1.0D-10))
        !call mix%set_material(2,stiffgas(gamma_2,Rgas_2,p_infty_2),sep1solid(rho_0_2,mu_2,yield2,1.0D-10))

        ! set logicals for plasticity
        mix%material(1)%plast = .false.;  mix%material(1)%explPlast = .false.
        !mix%material(2)%plast = plastic2; mix%material(2)%explPlast = explPlast2

        ! Set logicals for sliding
        mix%material(1)%sliding = .false.
        !mix%material(2)%sliding = sliding

        mix%material(1)%g11 = rho/rho_0; mix%material(1)%g12 = zero; mix%material(1)%g13 = zero
        mix%material(1)%g21 = zero;      mix%material(1)%g22 = one;  mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero;      mix%material(1)%g32 = zero; mix%material(1)%g33 = one

        mix%material(1)%p  = p
        mix%material(1)%VF = one
        mix%material(1)%Ys = one
        !! Get rho compatible with det(g) and rho0
        !tmp = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        !rho = rho0 * tmp

        !tmp = sqrt( gam*(p+pInf)/rho )    ! Speed of sound
        !tfactor = one / minval(tmp)
        !tstop = tstop * tfactor
        !tviz = tviz * tfactor

        print*, "rho1 = ", rho1
        print*, "rho2 = ", rho2
        print*, "u1 = ", u1
        print*, "u2 = ", u2
        print*, "p1 = ", p1
        print*, "p2 = ", p2
        !print*, "Pinf = ", PInf

        print*, "rho diff: ", (rho(1,1,1) - rho2)/rho2
    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                sxx_index,syy_index,szz_index,sxy_index,sxz_index,syz_index,sos_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use SolidMixtureMod,  only: solid_mixture

    use StiffgasShock_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der   
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    character(len=*),                intent(in) :: outputdir
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix
    integer,                         intent(in) :: vizcount
    integer, dimension(2),           intent(in) :: x_bc, y_bc, z_bc
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile, str
    integer :: i

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

        write(str,'(ES7.1E2,A1,ES7.1E2)') pRatio, "_", pinfbyp1
        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/StiffgasShock_"//trim(str)//"_", vizcount, ".dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        do i=1,decomp%ysz(1)
            write(outputunit,'(10ES26.16)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                 mix%material(1)%g11(i,1,1), mix%material(1)%g21(i,1,1), mu(i,1,1), bulk(i,1,1), kap(i,1,1)
        
        end do
        close(outputunit)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture
    use operators,        only: filter3D

    use StiffgasShock_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc

    integer :: nx

    nx = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        rho( 1,:,:) = rho2
        u  ( 1,:,:) = u2
        v  ( 1,:,:) = zero
        w  ( 1,:,:) = zero
        p  ( 1,:,:) = p2
        
        rho(nx,:,:) = rho1
        u  (nx,:,:) = u1
        v  (nx,:,:) = zero
        w  (nx,:,:) = zero
        p  (nx,:,:) = p1
        
    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind, clen
    use constants,        only: zero, half, one, two
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use SolidMixtureMod,  only: solid_mixture

    use StiffgasShock_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),             intent(in) :: mix

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

        !!dx = x(2,1,1) - x(1,1,1)

        !!p_exact = p2
        !!where (x(:,1,1) .GT. half)
        !!    p_exact = p1
        !!end where

        !!dpdx = zero
        !!dpdx(2:nx-1) = ( p(3:nx,1,1)-p(1:nx-2,1,1) ) / (two*dx)
        !!sthick = abs(p2-p1)/maxval(dx*abs(dpdx))

        !!istart = 1
        !!do while ( x(istart,1,1) .LT. half-sthick*dx )
        !!    istart = istart + 1
        !!end do

        !!iend = nx
        !!do while ( x(iend,1,1) .GT. half+sthick*dx )
        !!    iend = iend - 1
        !!end do

        !!! mwa = maxval(abs(p(1:istart,1,1)-p_exact(1:istart)))/abs(p2-p1)
        !!! mwa = max(mwa,maxval(abs(p(1:istart,1,1)-p_exact(1:istart)))/abs(p2-p1))

        !!mwa = maxval(p(:,1,1)-p2)
        !!mwa = max(mwa,maxval(p1-p(:,1,1)))
        !!mwa = mwa / abs(p2-p1)

        !!call message(2,"Shock thickness", sthick)
        !!call message(2,"Maximum Wiggles Amplitude", mwa)

        !!write(outputfile,'(A,ES8.2E2,A)') "StiffgasShock_stats_", Cbeta,".dat"
        !!if (step == 1) then
        !!    open(unit=iounit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
        !!    write(iounit,'(3A26)') "Time", "Shock thickness", "MWA"
        !!else
        !!    open(unit=iounit, file=trim(outputfile), form='FORMATTED', position='APPEND', status='OLD')
        !!end if
        !!write(iounit,'(3ES26.16)') tsim, sthick, mwa
        !!close(iounit)

    end associate
end subroutine

!!subroutine hook_source(decomp,mesh,fields,tsim,rhs,rhsg)
!!    use kind_parameters,  only: rkind
!!    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
!!    use decomp_2d,        only: decomp_info
!!    type(decomp_info),               intent(in)    :: decomp
!!    real(rkind),                     intent(in)    :: tsim
!!    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
!!    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
!!    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
!!    real(rkind), dimension(:,:,:,:), intent(inout) :: rhsg
!!
!!    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
!!                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
!!                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
!!                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
!!                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
!!                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
!!    end associate
!!end subroutine

subroutine hook_mixture_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,&
                                mom_index,TE_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use StiffgasShock_data

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

    use StiffgasShock_data

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

    use StiffgasShock_data

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

    use StiffgasShock_data

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

    use StiffgasShock_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
