module PlasticRelaxation_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none

    real(rkind) :: g0 = one
    real(rkind) :: tau_0
    real(rkind) :: sigma_Y
    real(rkind) :: shmod
    real(rkind) :: dtprob

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use PlasticRelaxation_data

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

        dx = one
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

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,rho0,mu,yield,gam,PInf,tau0,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,one,two,three,pi,four,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,&
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    
    use PlasticRelaxation_data

    implicit none
    character(len=*),                                               intent(in)    :: inputfile
    type(decomp_info),                                              intent(in)    :: decomp
    real(rkind),                                                    intent(in)    :: dx,dy,dz
    real(rkind),                                          optional, intent(inout) :: rho0,mu,gam,PInf,tstop,dt,tviz,tau0,yield
    real(rkind), dimension(:,:,:,:),     intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: detg

    integer :: ioUnit

    namelist /PROBINPUT/  g0
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    associate( rho => fields(:,:,:,rho_index),   u => fields(:,:,:,  u_index), &
                 v => fields(:,:,:,  v_index),   w => fields(:,:,:,  w_index), &
                 p => fields(:,:,:,  p_index),   T => fields(:,:,:,  T_index), &
                 e => fields(:,:,:,  e_index), g11 => fields(:,:,:,g11_index), &
               g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), & 
               g23 => fields(:,:,:,g23_index), g31 => fields(:,:,:,g31_index), & 
               g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        tau_0 = tau0
        sigma_Y = yield
        shmod = mu
        dtprob = dt

        u   = zero
        v   = zero
        w   = zero

        g11 = one+g0;  g12 = zero;   g13 = zero
        g21 = zero;    g22 = one;    g23 = zero
        g31 = zero;    g32 = zero;   g33 = one

        p = PInf*(g11**gam - one)

        ! Get rho compatible with det(g) and rho0
        detg = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        rho = rho0 * detg

    end associate

end subroutine

subroutine hook_output(decomp,dx,dy,dz,outputdir,mesh,fields,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,eps,half,twothird,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index, &
                                sxx_index,sxy_index,sxz_index,syy_index,syz_index,szz_index
    use decomp_2d,        only: decomp_info
    use reductions,       only: P_MAXVAL

    use PlasticRelaxation_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile, str
    integer :: i
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: invtaurel

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => fields(:,:,:,g11_index), g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), g23 => fields(:,:,:,g23_index), &
               g31 => fields(:,:,:,g31_index), g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 sxx  => fields(:,:,:, sxx_index), sxy => fields(:,:,:,sxy_index), & 
                 sxz  => fields(:,:,:, sxz_index), syy => fields(:,:,:,syy_index), & 
                 syz  => fields(:,:,:, syz_index), szz => fields(:,:,:,szz_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        ! Get S'S'
        invtaurel = sxx*sxx + two*sxy*sxy + two*sxz*sxz &
                            +     syy*syy + two*syz*syz &
                                          +     szz*szz

        ! 1/tau_rel
        invtaurel = ( invtaurel - (twothird)*sigma_Y**2 ) / shmod**2
        where (invtaurel .LE. zero)
            invtaurel = zero
        end where

        write(str,'(ES8.1E2,A1,ES7.1E2)') g0, "_", dtprob
        write(outputfile,'(2A)') trim(outputdir),"/PlasticRelaxation_"//trim(str)//".dat"

        if (vizcount == 0) then
            open(unit=outputunit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
            write(outputunit,'(12A26)') "Time", "rho", "u", "e", "p", "g11", "g22", "g33", "sxx", "syy", "szz", "yield criterion"
        else
            open(unit=outputunit, file=trim(outputfile), form='FORMATTED', position='APPEND', status='OLD')
        end if
        write(outputunit,'(12ES26.16)') tsim, rho(1,1,1), u(1,1,1), e(1,1,1), p(1,1,1), &
                                       g11(1,1,1), g22(1,1,1), g33(1,1,1), &
                                       sxx(1,1,1), syy(1,1,1), szz(1,1,1), P_MAXVAL(invtaurel)
        close(outputunit)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,tsim)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info

    use PlasticRelaxation_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
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

subroutine hook_timestep(decomp,mesh,fields,step,tsim)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,twothird,two
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,&
                                sxx_index,sxy_index,sxz_index,syy_index,syz_index,szz_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL

    use PlasticRelaxation_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    
    character(len=clen) :: charout
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: invtaurel

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 sxx  => fields(:,:,:, sxx_index), sxy => fields(:,:,:,sxy_index), & 
                 sxz  => fields(:,:,:, sxz_index), syy => fields(:,:,:,syy_index), & 
                 syz  => fields(:,:,:, syz_index), szz => fields(:,:,:,szz_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        ! Get S'S'
        invtaurel = sxx*sxx + two*sxy*sxy + two*sxz*sxz &
                            +     syy*syy + two*syz*syz &
                                          +     szz*szz

        ! 1/tau_rel
        invtaurel = ( invtaurel - (twothird)*sigma_Y**2 ) / shmod**2
        where (invtaurel .LE. zero)
            invtaurel = zero
        end where

        write(charout,'(A,ES26.16,A,ES26.16)') "Time = ", tsim, "    Maximum Yeild criterion = ", P_MAXVAL(invtaurel)
        call message(trim(charout))

    end associate
end subroutine

subroutine hook_source(decomp,mesh,fields,tsim,rhs,rhsg)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info

    use PlasticRelaxation_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhsg

    integer :: i,j,k

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
