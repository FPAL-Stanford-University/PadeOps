module taylorgreen_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none
    
    real(rkind) :: rho_0
    real(rkind) :: tke0, enstrophy0
    real(rkind) :: Re

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one,two,pi
    use decomp_2d,        only: decomp_info

    use taylorgreen_data

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

        dx = two*pi/real(nx,rkind)
        dy = two*pi/real(ny,rkind)
        dz = two*pi/real(nz,rkind)

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
    use constants,        only: zero,third,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,&
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info
    use exits,            only: GracefulExit

    use taylorgreen_data

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
    real(rkind) :: gam

    namelist /PROBINPUT/  Re
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    rho0  = one
    rho_0 = rho0

    if (eostype == 1) then
        gam = eosparams(1)
    else
        call GracefulExit("Can only run this problem with separable EOS",4576)
    end if

    associate( rho => fields(:,:,:,rho_index),   u => fields(:,:,:,  u_index), &
                 v => fields(:,:,:,  v_index),   w => fields(:,:,:,  w_index), &
                 p => fields(:,:,:,  p_index),   T => fields(:,:,:,  T_index), &
                 e => fields(:,:,:,  e_index), g11 => fields(:,:,:,g11_index), &
               g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), & 
               g23 => fields(:,:,:,g23_index), g31 => fields(:,:,:,g31_index), & 
               g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
       
        rho = one
        u   =  sin(x)*cos(y)*cos(z) 
        v   = -cos(x)*sin(y)*cos(z) 
        w   = zero
        p   = 100._rkind/gam + ( (cos(two*z) + two)*(cos(two*x) + cos(two*y)) - two ) / 16._rkind

        g11 = one;  g12 = zero; g13 = zero
        g21 = zero; g22 = one;  g23 = zero
        g31 = zero; g32 = zero; g33 = one

        ! Get rho compatible with det(g) and rho0
        tmp = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        rho = rho0 * tmp

    end associate

end subroutine

subroutine hook_output(decomp,der,fil,dx,dy,dz,outputdir,mesh,fields,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index, &
                                sxx_index,sxy_index,sxz_index,syy_index,syz_index,szz_index
    use decomp_2d,        only: decomp_info, nrank
    use DerivativesMod,   only: derivatives
    use FiltersMod,       only: filters
    use operators,        only: curl
    use reductions,       only: P_MEAN

    use taylorgreen_data

    implicit none
    character(len=clen),             intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(filters),                   intent(in) :: fil
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    real(rkind), dimension(2),       intent(in) :: x_bc, y_bc, z_bc
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile, str
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tke, enstrophy
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),3) :: vorticity

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

        ! Get TKE and Enstrophy to output to file
        tke = half*rho*(u*u + v*v + w*w)
        call curl(decomp, der, u, v, w, vorticity)
        enstrophy = vorticity(:,:,:,1)**2 + vorticity(:,:,:,2)**2 + vorticity(:,:,:,3)**2

        write(str,'(I4.4)') decomp%ysz(2)
        write(outputfile,'(2A)') trim(outputdir),"/taylorgreen_"//trim(str)//".dat"

        if (vizcount == 0) then
            tke0 = P_MEAN( tke )
            enstrophy0 = P_MEAN( vorticity(:,:,:,1)**2 + vorticity(:,:,:,2)**2 + vorticity(:,:,:,3)**2 )
            open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
            write(outputunit,'(3A26)') "Time", "TKE", "Enstrophy"
        else
            open(unit=outputunit, file=trim(outputfile), form='FORMATTED', position='APPEND', status='OLD')
        end if
        write(outputunit,'(3ES26.16)') tsim, P_MEAN(tke)/tke0, P_MEAN(enstrophy)/enstrophy0
        close(outputunit)

        write(str,'(I4.4,A,I4.4,A,I6.6)') decomp%ysz(2), "_", vizcount, "_", nrank
        write(outputfile,'(2A)') trim(outputdir),"/taylorgreen_"//trim(str)//".dat"
        open(unit=outputunit, file=trim(outputfile), form='UNFORMATTED', status='REPLACE')
        write(outputunit) tsim
        write(outputunit) decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)
        write(outputunit) decomp%yst(1), decomp%yst(2), decomp%yst(3)
        write(outputunit) decomp%yen(1), decomp%yen(2), decomp%yen(3)
        write(outputunit) rho
        write(outputunit) u
        write(outputunit) v
        write(outputunit) w
        write(outputunit) vorticity(:,:,:,1)
        write(outputunit) vorticity(:,:,:,2)
        write(outputunit) vorticity(:,:,:,3)
        write(outputunit) p
        close(outputunit)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, one
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info

    use taylorgreen_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind), dimension(2),       intent(in)    :: x_bc, y_bc, z_bc

    integer :: nx, ny

    nx = decomp%ysz(1); ny = decomp%ysz(2)

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

subroutine hook_timestep(decomp,der,mesh,fields,step,tsim,dt,x_bc,y_bc,z_bc,hookcond)
    use kind_parameters,  only: rkind
    use constants,        only: half
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use exits,            only: message
    use reductions,       only: P_MAXVAL, P_MEAN

    use taylorgreen_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim, dt
    real(rkind), dimension(2),       intent(in) :: x_bc,y_bc,z_bc
    logical,            optional, intent(inout) :: hookcond

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message("TKE",P_MEAN(half*rho*(u*u + v*v +w*w)))

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

subroutine hook_finalize

end subroutine
