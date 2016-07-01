module reconnection_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none

    real(rkind) :: Rgas, gamma, p_infty, rho_0
    real(rkind) :: muL, lamL, cL
    real(rkind) :: dtprob

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use reconnection_data

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

        dx = one/real(nx,rkind)
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
    use constants,        only: zero,third,half,one,two,three,pi,four,eight
    use SolidGrid,        only: u_index,v_index,w_index
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    use SolidMixtureMod,  only: solid_mixture
    
    use reconnection_data

    implicit none
    character(len=*),                                               intent(in)    :: inputfile
    type(decomp_info),                                              intent(in)    :: decomp
    real(rkind),                                                    intent(in)    :: dx,dy,dz
    type(solid_mixture),                                            intent(inout) :: mix
    real(rkind),                                          optional, intent(inout) :: tstop, dt, tviz
    real(rkind), dimension(:,:,:,:),     intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp
    real(rkind) :: t = zero, p_star, rho_star, c_star

    namelist /PROBINPUT/  p_infty, Rgas, gamma, muL, rho_0
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    associate(   u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        p_star = gamma*p_infty + (four/three)*muL
        rho_star = rho_0
        c_star = sqrt(p_star/rho_star)
        ! tstop = tstop * c_star
        ! dt = dt * c_star
        ! tviz = tviz * c_star

        dtprob = dt

        print*, "p_star = ", p_star
        print*, "rho_star = ", rho_star
        print*, "c_star = ", c_star
        print*, "tstop = ", tstop
        print*, "dt = ", dt
        print*, "tviz = ", tviz

        ! Non-dimensionalize problem parameters
        rho_0 = rho_0 / rho_star
        muL = muL / p_star
        p_infty = p_infty / p_star

        lamL = gamma * p_infty - (two/three)*muL
        cL = sqrt( (lamL + two*muL)/ rho_0 )

        ! Set material
        call mix%set_material(1,stiffgas(gamma,Rgas,p_infty),sep1solid(rho_0,muL,1.0D30,1.0D-10))

        tmp = (half*pi)*half*( erf( (x-0.25_rkind)/(dx) ) - erf( (x-0.75_rkind)/(dx) ) )

        u   = one
        v   = zero
        w   = zero

        mix%material(1)%g11 = cos(tmp); mix%material(1)%g12 = -sin(tmp); mix%material(1)%g13 = zero
        mix%material(1)%g21 = sin(tmp); mix%material(1)%g22 =  cos(tmp); mix%material(1)%g23 = zero
        mix%material(1)%g31 = zero;     mix%material(1)%g32 = zero;      mix%material(1)%g33 = one

        mix%material(1)%p = zero

    end associate

end subroutine

subroutine hook_output(decomp,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,eps,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use reconnection_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(solid_mixture),              intent(in) :: mix
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile, str
    integer :: i

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        write(str,'(I4.4,A1,ES7.1E2)') decomp%ysz(1), "_", dtprob
        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/reconnection_"//trim(str)//"_", vizcount, ".dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        write(outputunit,'(6A26)')     'tsim', 'muL', 'lamL', 'cL', 'p_infty'
        write(outputunit,'(6ES26.16)') tsim, muL, lamL, cL, p_infty
        write(outputunit,'(17A26)')    'x', 'rho', 'u', 'e', 'p', &
                                       'mix%material(1)%g11', 'mix%material(1)%g12', 'mix%material(1)%g13', &
                                       'mix%material(1)%g21', 'mix%material(1)%g22', 'mix%material(1)%g23', &
                                       'mix%material(1)%g31', 'mix%material(1)%g32', 'mix%material(1)%g33', &
                                       'mu', 'bulk', 'kap'
        do i=1,decomp%ysz(1)
            write(outputunit,'(17ES26.16)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                           mix%material(1)%g11(i,1,1), mix%material(1)%g12(i,1,1), mix%material(1)%g13(i,1,1), &
                                           mix%material(1)%g21(i,1,1), mix%material(1)%g22(i,1,1), mix%material(1)%g23(i,1,1), &
                                           mix%material(1)%g31(i,1,1), mix%material(1)%g32(i,1,1), mix%material(1)%g33(i,1,1), &
                                           mu(i,1,1), bulk(i,1,1), kap(i,1,1)
        end do
        close(outputunit)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim)
    use kind_parameters,  only: rkind
    use constants,        only: zero,half,one,pi
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use SolidMixtureMod,  only: solid_mixture

    use reconnection_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(solid_mixture),             intent(inout) :: mix

    real(rkind) :: tmp
    integer :: nx

    nx = decomp%ysz(1)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        ! rho( 1,:,:) = rho_0
        ! u  ( 1,:,:) = one
        ! v  ( 1,:,:) = zero
        ! w  ( 1,:,:) = zero

        ! mix%material(1)%p( 1,:,:) = zero

        ! tmp = half*pi
        ! mix%material(1)%g11( 1,:,:) = cos(tmp); mix%material(1)%g12( 1,:,:) = -sin(tmp); mix%material(1)%g13( 1,:,:) = zero
        ! mix%material(1)%g21( 1,:,:) = sin(tmp); mix%material(1)%g22( 1,:,:) =  cos(tmp); mix%material(1)%g23( 1,:,:) = zero
        ! mix%material(1)%g31( 1,:,:) = zero;     mix%material(1)%g32( 1,:,:) = zero;      mix%material(1)%g33( 1,:,:) = one
        ! 
        ! rho(nx,:,:) = rho_0
        ! u  (nx,:,:) = one
        ! v  (nx,:,:) = zero
        ! w  (nx,:,:) = zero

        ! mix%material(1)%p(nx,:,:) = zero

        ! tmp = zero
        ! mix%material(1)%g11(nx,:,:) = cos(tmp); mix%material(1)%g12(nx,:,:) = -sin(tmp); mix%material(1)%g13(nx,:,:) = zero
        ! mix%material(1)%g21(nx,:,:) = sin(tmp); mix%material(1)%g22(nx,:,:) =  cos(tmp); mix%material(1)%g23(nx,:,:) = zero
        ! mix%material(1)%g31(nx,:,:) = zero;     mix%material(1)%g32(nx,:,:) = zero;      mix%material(1)%g33(nx,:,:) = one
    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL
    use SolidMixtureMod,  only: solid_mixture

    use reconnection_data

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

    use reconnection_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    type(solid_mixture),             intent(in)    :: mix

    integer :: i,j,k

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

    use reconnection_data

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

    use reconnection_data

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

    use reconnection_data

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

    use reconnection_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(stiffgas),                  intent(in)    :: hydro
    type(sep1solid),                 intent(in)    :: elastic
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
    real(rkind), dimension(:,:,:),   intent(inout) :: rhs

end subroutine
