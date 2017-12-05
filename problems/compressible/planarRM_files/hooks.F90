module planarRM_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none

    integer, parameter :: ns = 4

    ! Parameters for the 4 materials:               Nitrogen         Oxygen             SF6        Acetone
    !                                        ---------------------------------------------------------------
    real(rkind), dimension(ns) :: gam      = [     1.4_rkind,     1.4_rkind,      1.1_rkind,     1.1_rkind ]
    real(rkind), dimension(ns) :: Pr       = [    0.72_rkind,    0.72_rkind,      0.8_rkind,     0.8_rkind ]
    real(rkind), dimension(ns) :: molwt    = [ 28.0140_rkind, 31.9990_rkind, 146.0570_rkind, 58.0805_rkind ] ! g/mol
    real(rkind), dimension(ns) :: sigma    = [  3.7380_rkind,  3.4800_rkind,   5.1990_rkind,  4.5990_rkind ] ! Angstrom
    real(rkind), dimension(ns) :: eps_by_k = [    82.0_rkind,   102.6_rkind,    212.0_rkind,   458.0_rkind ] ! K
    real(rkind), dimension(ns) :: Rgas
    real(rkind) :: thick = one

    ! Parameters for collision integral for Reid diffusivity
    real(rkind) :: diffA = 1.06036_rkind, diffB =  -0.1561_rkind, diffC = 0.19300_rkind, diffD = -0.47635_rkind, &
                   diffE = 1.03587_rkind, diffF = -1.52996_rkind, diffG = 1.76474_rkind, diffH = -3.89411_rkind
    real(rkind) :: diffConst = real(2.66D-2, rkind)

    ! Parameters for collision integral for Chapman-Enskog viscosity
    real(rkind) :: viscA = 1.16145_rkind, viscB = -0.14874_rkind, viscC =  0.52487_rkind, &
                   viscD = -0.7732_rkind, viscE =  2.16178_rkind, viscF = -2.43787_rkind
    real(rkind) :: viscConst = real(2.6693D-6, rkind)

    real(rkind), parameter :: R_univ = 8.3144598_rkind ! Universal gas constant in SI units

    ! Domain size data
    real(rkind) :: L_x = 0.5_rkind, L_yz = 0.1_rkind
    real(rkind) :: x1 = -0.1_rkind, yz1 = -0.05_rkind

    ! Ambient state
    real(rkind) :: p_amb = real(2.3D4,rkind), T_amb = real(298.0,rkind)
    real(rkind) :: X_N2 = 0.79_rkind, X_O2 = 0.21_rkind
    real(rkind) :: Y_SF6 = 0.8_rkind, Y_Acetone = 0.2_rkind

    ! Initial shock params
    real(rkind) :: x_shock = -0.05_rkind
    real(rkind) :: M_shock = 1.5_rkind

    ! Initial interface params
    real(rkind) :: x_int = 0._rkind
    real(rkind) :: L_int = 0.01_rkind

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use planarRM_data

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

        dx = L_x /real(nx-1,rkind)
        dy = L_yz/real(ny  ,rkind)
        dz = L_yz/real(nz  ,rkind)

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = x1  + real( ix1 - 1 + i - 1, rkind ) * dx
                    y(i,j,k) = yz1 + real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = yz1 + real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz)
    use kind_parameters,             only: rkind
    use constants,                   only: zero,half,one,two,pi,eight
    use CompressibleGrid,            only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,Ys_index
    use decomp_2d,                   only: decomp_info
    use MixtureEOSMod,               only: mixture
    use IdealGasEOS,                 only: idealgas
    use ChapmanEnskogViscosityMod,   only: chapmanEnskogViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use ConstPrandtlConductivityMod, only: constPrandtlConductivity
    use ReidRamshawDiffusivityMod,   only: reidRamshawDiffusivity
    use exits,                       only: GracefulExit
    
    use planarRM_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tstop, dt, tviz

    integer :: i, iounit
    real(rkind) :: Y_N2, Y_O2, rho_air, rho_N2, rho_O2, rho_HG
    real(rkind) :: R_air, Cp_air, gamma_air, Cp_N2, Cp_O2, sos_air, R_HG
    real(rkind) :: rho_shocked, u_shocked, p_shocked
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp

    type(chapmanEnskogViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(constPrandtlConductivity) :: thermcond

    namelist /PROBINPUT/  thick, x_int, L_int, x_shock, M_shock, p_amb, T_amb
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),                    &
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),                    &
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),                    &
                 e => fields(:,:,:,  e_index), Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (mix%ns /= ns) call GracefulExit("Wrong number of species. Check your input file and make ns consistent with the problem file.",4562)

        ! First set the materials
        Rgas = R_univ / (molwt * real(1.D-3,rkind))   ! Set Rgas = R_univ / molwt. Note: converting molwt to SI units

        ! Set each material's transport coefficient object
        do i = 1,mix%ns
            shearvisc = chapmanEnskogViscosity( viscConst, molwt(i), eps_by_k(i), sigma(i), viscA, viscB, viscC, viscD, viscE, viscF )
            bulkvisc  = constRatioBulkViscosity( zero )
            thermcond = constPrandtlConductivity( Pr(i) )
            call mix%set_material( i, idealgas( gam(i), Rgas(i) ),&
                 shearvisc = shearvisc, &
                 bulkvisc  = bulkvisc, &
                 thermcond = thermcond  )
            ! call mix%set_material( i, idealgas( gam(i), Rgas(i) ), &
            !      shearvisc = chapmanEnskogViscosity( viscConst, molwt(i), eps_by_k(i), sigma(i), viscA, viscB, viscC, viscD, viscE, viscF ), &
            !      bulkvisc  = constRatioBulkViscosity( zero ), &
            !      thermcond = constPrandtlConductivity( Pr(i) )  )
        end do

        ! Set mass diffusivity object (Ensure that all units are consistent)
        call mix%set_massdiffusivity( reidRamshawDiffusivity(mix%ns, diffConst, molwt, sigma, eps_by_k, &
             diffA, diffB, diffC, diffD, diffE, diffF, diffG, diffH) )

        rho_N2 = p_amb / (Rgas(1) * T_amb)
        rho_O2 = p_amb / (Rgas(2) * T_amb)
        rho_air = (X_N2 * rho_N2 + X_O2 * rho_O2)
        Y_N2 = rho_N2 * X_N2 / rho_air
        Y_O2 = rho_O2 * X_O2 / rho_air

        ! Y_N2 = 0.767_rkind
        ! Y_O2 = 0.233_rkind
        ! R_air = Y_N2 * Rgas(1) + Y_O2 * Rgas(2)
        ! rho_air = p_amb / (R_air * T_amb)

        print *, "X_N2 = ", X_N2, "X_O2 = ", X_O2
        print *, "Y_N2 = ", Y_N2, "Y_O2 = ", Y_O2

        Cp_N2 = gam(1) * Rgas(1) / (gam(1) - one)
        Cp_O2 = gam(2) * Rgas(2) / (gam(2) - one)
        Cp_air = Y_N2 * Cp_N2   + Y_O2 * Cp_O2
        R_air  = Y_N2 * Rgas(1) + Y_O2 * Rgas(2)
        gamma_air = Cp_air / (Cp_air - R_air)
        print *, "Gamma for air = ", gamma_air

        R_HG = Y_SF6 * Rgas(3) + Y_Acetone * Rgas(4)
        rho_HG = p_amb / (R_HG * T_amb)

        tmp = half * (one + tanh( (x - x_int)/L_int ))

        rho = rho_air * (one-tmp) + rho_HG * tmp

        ! Set the massfractions (must sum to unity)
        Ys(:,:,:,1) = Y_N2 * (one - tmp)
        Ys(:,:,:,2) = Y_O2 * (one - tmp)
        Ys(:,:,:,3) = Y_SF6 * tmp
        Ys(:,:,:,4) = one - sum(Ys(:,:,:,1:mix%ns-1), 4)

        print*, "Max Ys = ", maxval(Ys(:,:,:,1)), ", Min Ys = ", minval(Ys(:,:,:,1))
        print*, "Max Ys = ", maxval(Ys(:,:,:,2)), ", Min Ys = ", minval(Ys(:,:,:,2))
        print*, "Max Ys = ", maxval(Ys(:,:,:,3)), ", Min Ys = ", minval(Ys(:,:,:,3))
        print*, "Max Ys = ", maxval(Ys(:,:,:,4)), ", Min Ys = ", minval(Ys(:,:,:,4))
        print*, "Max sum = ", maxval(sum(Ys, 4)), "Min sum = ", minval(sum(Ys, 4))

        rho_shocked = rho_air * (gamma_air + one) * M_shock**2 / (two + (gamma_air - one) * M_shock**2)

        sos_air = sqrt(gamma_air * p_amb / rho_air)
        u_shocked = M_shock * sos_air * (one - rho_air / rho_shocked)

        p_shocked = p_amb * (one + two * (gamma_air/(gamma_air+one)) * (M_shock**2 - one) )

        print *, "rho_air = ", rho_air
        print *, "rho_shocked = ", rho_shocked
        print *, "rho_HG = ", rho_HG
        print *, "u_shocked = ", u_shocked
        print *, "p_amb = ", p_amb
        print *, "p_shocked = ", p_shocked

        tmp = half * (one + tanh( (x - x_shock)/(thick*dx) ))

        rho = rho_shocked * (one-tmp) + rho * tmp
        u = u_shocked * (one - tmp)
        v   = zero
        w   = zero
        p = p_amb * tmp + p_shocked * (one - tmp)

    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture

    use planarRM_data

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

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/planarRM_", vizcount, ".dat"

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

    use planarRM_data

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

    use planarRM_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(mixture),                   intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields

    real(rkind) :: dx, Ythick, oob
    integer :: ny
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
        
        ny = decomp%ysz(2)
        dx = x(2,1,1) - x(1,1,1)

        write(outputfile,'(A,I4.4,A)') "planarRM_stats_N", ny, ".dat"
        if (step == 1) then
            open(unit=iounit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
            write(iounit,'(3A26)') "Time", "Shock thickness", "MWA"
        else
            open(unit=iounit, file=trim(outputfile), form='FORMATTED', position='APPEND', status='OLD')
        end if
        write(iounit,'(ES26.16)') tsim
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

    use planarRM_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

end subroutine
