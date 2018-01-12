module species_shocktube_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none

    real(rkind), dimension(2) :: gam  = [1.4_rkind, 1.6_rkind]
    real(rkind), dimension(2) :: Rgas = [one, one]
    real(rkind) :: rhoRatio = 0.125_rkind
    real(rkind) :: p0 = one, p1 = 0.1_rkind
    real(rkind) :: thick = one
    real(rkind) :: Cdiff, CY

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use species_shocktube_data

    implicit none

    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Need to set x, y and z as well as  dx, dy and dz
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = one/real(nx-1,rkind)
        dy = dx
        dz = dx

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1 + i - 1, rkind ) * dx - half
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use IdealGasEOS,      only: idealgas
    use exits,            only: GracefulExit
    
    use species_shocktube_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    integer :: iounit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp

    namelist /PROBINPUT/  thick, Cdiff, CY
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),                    &
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),                    &
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),                    &
                 e => fields(:,:,:,  e_index), Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (mix%ns /= 2) call GracefulExit("This problem can only be run with 2 species. Check your input file.",4562)

        ! First set the materials
        Rgas(1) = p0   ! Set Rgas(1) = p0/(rho0*T0) based on rho0 = 1, T0 = 1
        Rgas(2) = p1 / (rhoRatio)  ! Set Rgas(2) = p1/(rho1*T1) based on rho1 = rhoRatio, T1 = 1
        call mix%set_material( 1, idealgas(gam(1),Rgas(1)) )
        call mix%set_material( 2, idealgas(gam(2),Rgas(2)) )

        ! Set the massfractions (must sum to unity)
        Ys(:,:,:,2) = half*( one + erf( (x)/(thick*dx) ) )
        Ys(:,:,:,1) = one - Ys(:,:,:,2)

        print*, "Max Ys = ", maxval(Ys(:,:,:,2)), ", Min Ys = ", minval(Ys(:,:,:,2))

        ! Use pressure and temperature equilibrium to set density
        p = p0 * Ys(:,:,:,1) + p1 * Ys(:,:,:,2)
        T = one
        call mix%update(Ys)
        rho = p / (mix%Rgas * T)

        print*, "Max rho = ", maxval(rho), ", Min rho = ", minval(rho)

        u   = zero
        v   = zero
        w   = zero

    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture

    use species_shocktube_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(mixture),                   intent(in) :: mix
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile
    integer :: i

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/species_shocktube_", vizcount, ".dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        do i=1,decomp%ysz(1)
            write(outputunit,'(12ES26.16)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                           mu(i,1,1), bulk(i,1,1), kap(i,1,1), Ys(i,1,1,1), Ys(i,1,1,2), &
                                           diff(i,1,1,1), diff(i,1,1,2)
        
        end do
        close(outputunit)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture

    use species_shocktube_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,two
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use exits,            only: message
    use reductions,       only: P_MAXVAL,P_MINVAL

    use species_shocktube_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(mixture),                   intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields

    real(rkind) :: dx, Ythick, oob
    real(rkind), dimension(decomp%ysz(1)) :: dY
    integer :: nx
    integer :: iounit = 229
    character(len=clen) :: outputfile

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        nx = decomp%ysz(1)
        dx = x(2,1,1) - x(1,1,1)

        dY = zero
        dY(2:nx-1) = (Ys(3:nx,1,1,1)-Ys(1:nx-2,1,1,1))/(two*dx)
        Ythick = one/maxval(dx*abs(dY))

        oob = maxval( half*(abs(Ys(:,1,1,1))-one + abs(Ys(:,1,1,1)-one)) )

        write(outputfile,'(A,ES8.2E2,A,ES8.2E2,A)') "species_shocktube_stats_", Cdiff, "_", CY, ".dat"
        if (step == 1) then
            open(unit=iounit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
            write(iounit,'(3A26)') "Time", "Shock thickness", "MWA"
        else
            open(unit=iounit, file=trim(outputfile), form='FORMATTED', position='APPEND', status='OLD')
        end if
        write(iounit,'(3ES26.16)') tsim, Ythick, oob
        close(iounit)

        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Minimum bulk viscosity",P_MINVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message(2,"Maximum diffusivity",P_MAXVAL(diff))
        call message(2,"Minimum diffusivity",P_MINVAL(diff))

    end associate
end subroutine

subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters, only: rkind
    use decomp_2d,       only: decomp_info
    use MixtureEOSMod,    only: mixture

    use species_shocktube_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

end subroutine
