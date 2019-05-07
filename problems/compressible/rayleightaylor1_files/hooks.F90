module rayleightaylor1_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, one
    use FiltersMod,       only: filters
    implicit none


    integer, parameter :: ns = 2
    logical :: threeD = .true.

    ! Problem parameters
    real(rkind) :: z_interface = 5.5_rkind    ! Interface location at y=0
    real(rkind) :: p_pre = real(1.0D2, rkind) ! Ambient pressure in Pa
    real(rkind) :: T_pre = 298._rkind         ! Ambient temperature in K
    real(rkind) :: rhoRatio = 3.0_rkind       ! Density ratio
    real(rkind) :: Re = 100.0_rkind           ! Reynolds' number
    real(rkind) :: Sc = 1.0_rkind             ! Schmidt number
    real(rkind) :: gravity = 1.0_rkind        ! gravity
    real(rkind) :: L_int = 0.5_rkind          ! Interface thickness
    real(rkind) :: bulk_Ratio = 1.0_rkind     ! Constant ratio B2/B1 (may later be T-dependent)
    real(rkind) :: mu_Ratio = 1.0_rkind       ! Constant ratio mu2/mu1 (may later be T-dependent)
    real(rkind) :: kap_Ratio = 1.0_rkind      ! Constant ratio k2/k1 (may later be T-dependent)
    ! Parameters for the 2 materials:               Nitrogen            CO2
    !                                        --------------------------------
    real(rkind), dimension(ns) :: gam      = [     1.66667_rkind,    1.66667_rkind ]
    real(rkind), dimension(ns) :: Pr       = [    0.72_rkind,    0.77_rkind ]
    real(rkind), dimension(ns) :: molwt    = [ 28.0130_rkind,  44.010_rkind ] ! g/mol
    real(rkind), dimension(ns) :: sigma    = [  3.6810_rkind,  3.9520_rkind ] ! Angstrom
    real(rkind), dimension(ns) :: eps_by_k = [   91.46_rkind,   200.0_rkind ] ! K
    real(rkind), dimension(ns) :: Rgas
    real(rkind) :: thick = one


    ! Domain size data
    real(rkind) :: L_z = 5.0_rkind
    real(rkind) :: L_x = 1.0_rkind, L_y = 1.0_rkind
    real(rkind) :: x1, y1, z1

    logical :: periodicx = .true., periodicy = .true., periodicz = .false.


    ! Gaussian filter for sponge
    type(filters) :: mygfil


contains


    subroutine get_volumefractions(mix,Ys,Xs)
        use MixtureEOSMod,               only: mixture
        type(mixture),                   intent(in)  :: mix
        real(rkind), dimension(:,:,:,:), intent(in)  :: Ys    ! Species mass fractions
        real(rkind), dimension(size(Ys,1), size(Ys,2), size(Ys,3), size(Ys,4)), intent(out) :: Xs    ! Species volume fractions
        
        integer :: n

        ! Get volume fractons
        do n = 1,mix%ns
            Xs(:,:,:,n) = mix%material(n)%mat%Rgas * Ys(:,:,:,n) / mix%Rgas
        end do
   
    end subroutine

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info

    use rayleightaylor1_data

    implicit none

    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k,ioUnit
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Need to set x, y and z as well as  dx, dy and dz
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )


        dx = L_x/real(nx  ,rkind)
        dy = L_y/real(ny  ,rkind)
        dz = L_z/real(nz-1,rkind)


        x1 = zero
        y1 = zero
        z1 = zero

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = x1 + real( ix1 - 1 + i - 1, rkind ) * dx
                    y(i,j,k) = y1 + real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = z1 + real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
    use kind_parameters,             only: rkind, clen
    use constants,                   only: zero,half,one,two,pi,eight
    use CompressibleGrid,            only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,Ys_index
    use decomp_2d,                   only: decomp_info
    use MixtureEOSMod,               only: mixture
    use IdealGasEOS,                 only: idealgas
    use PowerLawViscosityMod,        only: powerLawViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use ConstPrandtlConductivityMod, only: constPrandtlConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
    use io_hdf5_stuff,               only: io_hdf5
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use mpi
    
    use rayleightaylor1_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    integer :: i, j, k, l, mode, iounit

    type(powerLawViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(constPrandtlConductivity) :: thermcond

    real(rkind) :: rho_ref = 1.0_rkind, T_ref = 1.0_rkind, viscous_exponent = 0.75_rkind, Pr1, Sr, amp, At, A_minus, A_plus
    real(rkind), dimension(ns) :: mu_ref, bulk, Pr_mod

    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: tmp, tmprho, H_minus, H_plus

    integer :: nx, ny, nz

    namelist /PROBINPUT/  threeD, rhoRatio, Re, Sc, Pr1, Sr, gravity, z_interface, L_int, amp, bulk_Ratio, mu_Ratio, kap_Ratio
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),                    &
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),                    &
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),                    &
                 e => fields(:,:,:,  e_index), Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (mix%ns /= ns) call GracefulExit("Wrong number of species. Check your input file and make ns consistent with the problem file.",4562)

        ! First set the material properties and other quantities
        Rgas = [(1.0_rkind-1.0_rkind/gam(1)), (1.0_rkind-1.0_rkind/gam(1))/rhoRatio]
        At = (rhoRatio-1.0_rkind)/(rhoRatio+1.0_rkind)
        A_plus = -Sr/(one+At)
        A_minus = -Sr/(one-At)

        ! Set each material's transport coefficient object
        mu_ref = [one/Re, one*mu_Ratio/Re]
        ! This nondimensionalization requires a modified Pr
        Pr_mod = [Pr1*Rgas(1)*gam(1)/(gam(1)-one), Pr1*Rgas(2)*gam(2)*mu_Ratio/(kap_ratio*(gam(2)-one))]
        do i = 1,mix%ns
            shearvisc = powerLawViscosity( mu_ref(i), T_ref, zero )
            bulkvisc  = constRatioBulkViscosity( zero )
            thermcond = constPrandtlConductivity( Pr_mod(i) )
            call mix%set_material( i, idealgas( gam(i), Rgas(i) ),&
                 shearvisc = shearvisc, &
                 bulkvisc  = bulkvisc, &
                 thermcond = thermcond  )
        end do

        ! Set mass diffusivity object (Ensure that all units are consistent)
        call mix%set_massdiffusivity( constSchmidtDiffusivity(mu_ref(1), rho_ref, Sc) )

        do k = 1,decomp%ysz(3)
            tmp(:,:,k) = z(:,:,k)-z_interface-amp*cos(two*3.1415_rkind/L_x*x(:,:,k))*exp(-two*3.1415_rkind/L_x*abs(z(:,:,k)-z_interface));
            H_plus(:,:,k) = (one+erf(tmp(:,:,k)/L_int))/2
            H_minus(:,:,k) = (one-erf(tmp(:,:,k)/L_int))/2
            rho(:,:,k) = (one+At)*exp(A_minus*tmp(:,:,k))*H_plus(:,:,k)+(one-At)*exp(A_plus*tmp(:,:,k))*H_minus(:,:,k)
        end do

        do k = 1,decomp%ysz(3)
            Ys(:,:,k,2) = (one+At)*exp(A_minus*tmp(:,:,k))*H_plus(:,:,k)/rho(:,:,k)
        end do
         
        Ys(:,:,:,1)  = one - Ys(:,:,:,2)
        call mix%update(Ys)

        u   = zero
        v   = zero
        w   = zero

        T = 414640_rkind
        p = (rho*mix%Rgas)*T

        ! Initialize mygfil
        call mygfil%init(                           decomp, &
                          periodicx,  periodicy, periodicz, &
                         "gaussian", "gaussian", "gaussian" )

    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture

    use rayleightaylor1_data

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

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/rayleightaylor_", vizcount, ".dat"

        ! open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        ! do i=1,decomp%ysz(1)
        !     write(outputunit,'(12ES26.16)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
        !                                    mu(i,1,1), bulk(i,1,1), kap(i,1,1), Ys(i,1,1,1), Ys(i,1,1,2), &
        !                                    diff(i,1,1,1), diff(i,1,1,2)
        ! 
        ! end do
        ! close(outputunit)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use rayleightaylor1_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    integer :: i
    real(rkind) :: dx, filpt, thickT
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
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

    use rayleightaylor1_data

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
        
        ! ny = decomp%ysz(2)
        ! dx = x(2,1,1) - x(1,1,1)

        ! write(outputfile,'(A,I4.4,A)') "rayleightaylor_stats_N", ny, ".dat"
        ! if (step == 1) then
        !     open(unit=iounit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
        !     write(iounit,'(3A26)') "Time", "Shock thickness", "MWA"
        ! else
        !     open(unit=iounit, file=trim(outputfile), form='FORMATTED', position='APPEND', status='OLD')
        ! end if
        ! write(iounit,'(ES26.16)') tsim
        ! close(iounit)

        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message(2,"Maximum diffusivity",P_MAXVAL(diff))

    end associate
end subroutine

subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs,rhsg)
    use kind_parameters, only: rkind
    use decomp_2d,       only: decomp_info
    use MixtureEOSMod,    only: mixture
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index,mom_index, TE_index

    use rayleightaylor1_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    real(rkind), dimension(:,:,:,:), optional, intent(inout) ::rhsg

    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3), 3) :: grav 
    
    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
       
        ! gravity force 
        grav(:,:,:,1) = 0.0_rkind
        grav(:,:,:,2) = 0.0_rkind
        grav(:,:,:,3) = -rho*gravity

        rhs(:,:,:,mom_index:mom_index+2) = rhs(:,:,:,mom_index:mom_index+2) + grav
        rhs(:,:,:,TE_index) = rhs(:,:,:,TE_index) + grav(:,:,:,1)*u + grav(:,:,:,2)*v + grav(:,:,:,3)*w
    end associate
end subroutine
