module impact_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none
    
    real(rkind) :: uimpact = 100._rkind
    real(rkind) :: pinit   = real(1.0D5,rkind), Tinit
    real(rkind) :: rho0, gam, Rgas, PInf, mu, yield, tau0, K0, Cv, T0, B0, alpha, beta

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use impact_data

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
    use kind_parameters,        only: rkind
    use constants,              only: zero,half,one
    use exits,                  only: GracefulExit
    use SolidGrid,              only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index
    use decomp_2d,              only: decomp_info
    use SolidMixtureMod,        only: solid_mixture
    use AbstractEOSMod,         only: abstracteos
    use Sep1SolidEOSMod,        only: sep1solideos
    use GodunovRomenskiiEOSMod, only: godromeos
    
    use impact_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: dx,dy,dz
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(inout) :: tstop, dt, tviz

    ! class(abstracteos), allocatable :: eos
    type(sep1solideos), allocatable :: eos
    type(godromeos   ), allocatable :: eos2
    integer :: eostype, ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp

    namelist /PROBINPUT/  eostype, uimpact, &
                          gam, Rgas, PInf, rho0, mu, yield, tau0, &
                          K0, Cv, T0, B0, alpha, beta
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    ! Need to set mixture density and velocities and set material eos, g, VF and T (and material energy if not using PTeqb)
    associate( rho => fields(:,:,:,rho_index),   u => fields(:,:,:,  u_index), &
                 v => fields(:,:,:,  v_index),   w => fields(:,:,:,  w_index), &
                 e => fields(:,:,:,  e_index),                                 &
                 p => mix%material(1)%p,         T => mix%material(1)%T,       &
               g11 => mix%material(1)%g11, g12 => mix%material(1)%g12, g13 => mix%material(1)%g13, & 
               g21 => mix%material(1)%g21, g22 => mix%material(1)%g22, g23 => mix%material(1)%g23, &
               g31 => mix%material(1)%g31, g32 => mix%material(1)%g32, g33 => mix%material(1)%g33, & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
     
        select case(eostype)
        case(1)
            eos = sep1solideos(gam,Rgas,PInf,rho0,mu,yield,tau0,mix%usegTg)
            call mix%set_material(1, eos)
        case(2)
            eos2 = godromeos(rho0,K0,Cv,T0,B0,alpha,beta,gam,mix%usegTg)
            call mix%set_material(1, eos2)
        case default
            call GracefulExit("Only eostype = 1 (Sep1Solid) or 2 (GodunovRomenskii) are supported.",456)
        end select
        !allocate( eos, source=sep1solideos(gam,Rgas,PInf,rho0,mu,yield,tau0,mix%usegTg) )


        tmp = tanh( (x-half)/dx )

        u   = -uimpact*tmp
        v   = zero
        w   = zero
        p   = pinit

        g11 = one;  g12 = zero; g13 = zero
        g21 = zero; g22 = one;  g23 = zero
        g31 = zero; g32 = zero; g33 = one

        ! Get rho compatible with det(g) and rho0
        tmp = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        rho = rho0 * tmp
      
        select case(eostype)
        case(1)
            ! Get temperature by first getting hydro energy
            call eos%hydro%get_e_from_p(rho,p,e)
            call eos%hydro%get_T(e,T,rho)
        case(2)
            T = T0
        end select
        Tinit = T(1,1,1)

        ! Set volume fraction to one since only 1 material is present
        mix%material(1)%VF = one

        if(allocated(eos )) deallocate(eos )
        if(allocated(eos2)) deallocate(eos2)
    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use SolidMixtureMod,  only: solid_mixture

    use impact_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(solid_mixture),             intent(in) :: mix
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    integer, dimension(2),           intent(in) :: x_bc,y_bc,z_bc

    integer                                     :: outputunit=229
    character(len=clen) :: outputfile, velstr
    integer :: i

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 e => fields(:,:,:,  e_index),                                     &
                 p => mix%material(1)%p,           T => mix%material(1)%T,         &
                 mu  => fields(:,:,:, mu_index), bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => mix%material(1)%g11, g12 => mix%material(1)%g12, g13 => mix%material(1)%g13, & 
               g21 => mix%material(1)%g21, g22 => mix%material(1)%g22, g23 => mix%material(1)%g23, &
               g31 => mix%material(1)%g31, g32 => mix%material(1)%g32, g33 => mix%material(1)%g33, & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        write(velstr,'(I3.3)') int(uimpact)
        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/impact_"//trim(velstr)//"_", vizcount, ".dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        do i=1,decomp%ysz(1)
            write(outputunit,'(10ES26.16)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                           g11(i,1,1), g21(i,1,1), mu(i,1,1), bulk(i,1,1), kap(i,1,1)
        
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

    use impact_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix
    real(rkind),                     intent(in)    :: tsim
    integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc

    integer :: nx

    nx = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 e => fields(:,:,:,  e_index),                                     &
                 p => mix%material(1)%p,           T => mix%material(1)%T,         &
                 mu  => fields(:,:,:, mu_index), bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => mix%material(1)%g11, g12 => mix%material(1)%g12, g13 => mix%material(1)%g13, & 
               g21 => mix%material(1)%g21, g22 => mix%material(1)%g22, g23 => mix%material(1)%g23, &
               g31 => mix%material(1)%g31, g32 => mix%material(1)%g32, g33 => mix%material(1)%g33, & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        rho( 1,:,:) = rho0
        u  ( 1,:,:) = uimpact
        v  ( 1,:,:) = zero
        w  ( 1,:,:) = zero
        p  ( 1,:,:) = pinit
        T  ( 1,:,:) = Tinit
        
        g11( 1,:,:) = one;  g12( 1,:,:) = zero; g13( 1,:,:) = zero
        g21( 1,:,:) = zero; g22( 1,:,:) = one;  g23( 1,:,:) = zero
        g31( 1,:,:) = zero; g32( 1,:,:) = zero; g33( 1,:,:) = one

        rho(nx,:,:) = rho0
        u  (nx,:,:) = -uimpact
        v  (nx,:,:) = zero
        w  (nx,:,:) = zero
        p  (nx,:,:) = pinit
        T  (nx,:,:) = Tinit
        
        g11(nx,:,:) = one;  g12(nx,:,:) = zero; g13(nx,:,:) = zero
        g21(nx,:,:) = zero; g22(nx,:,:) = one;  g23(nx,:,:) = zero
        g31(nx,:,:) = zero; g32(nx,:,:) = zero; g33(nx,:,:) = one

    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture
    use exits,            only: message
    use reductions,       only: P_MAXVAL

    use impact_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(solid_mixture),             intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 e => fields(:,:,:,  e_index),                                     &
                 p => mix%material(1)%p,           T => mix%material(1)%T,         &
                 mu  => fields(:,:,:, mu_index), bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => mix%material(1)%g11, g12 => mix%material(1)%g12, g13 => mix%material(1)%g13, & 
               g21 => mix%material(1)%g21, g22 => mix%material(1)%g22, g23 => mix%material(1)%g23, &
               g31 => mix%material(1)%g31, g32 => mix%material(1)%g32, g33 => mix%material(1)%g33, & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))

    end associate
end subroutine

subroutine hook_mixture_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(solid_mixture),             intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
end subroutine

subroutine hook_material_g_source(decomp,eos,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info
    use AbstractEOSMod,   only: abstracteos

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    class(abstracteos),              intent(in)    :: eos
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
end subroutine

subroutine hook_material_mass_source(decomp,eos,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info
    use AbstractEOSMod,   only: abstracteos

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    class(abstracteos),              intent(in)    :: eos
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs
end subroutine

subroutine hook_material_VF_source(decomp,eos,x,y,z,tsim,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info
    use AbstractEOSMod,   only: abstracteos

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    class(abstracteos),              intent(in)    :: eos
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs
end subroutine

subroutine hook_material_energy_source(decomp,eos,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info
    use AbstractEOSMod,   only: abstracteos

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    class(abstracteos),              intent(in)    :: eos
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs
end subroutine
