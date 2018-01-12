module inclinedRM_data
    use kind_parameters,  only: rkind
    use constants,        only: zero, one
    implicit none

    ! Using SI units for this problem
    ! Miranda restart files are in CGS units, so need to make sure
    ! they are converted to SI units properly

    integer, parameter :: ns = 2

    ! Parameters for the 2 materials:               Nitrogen            CO2
    !                                        --------------------------------
    real(rkind), dimension(ns) :: gam      = [     1.4_rkind,    1.28_rkind ]
    real(rkind), dimension(ns) :: Pr       = [    0.72_rkind,    0.77_rkind ]
    real(rkind), dimension(ns) :: molwt    = [ 28.0130_rkind,  44.010_rkind ] ! g/mol
    real(rkind), dimension(ns) :: sigma    = [  3.6810_rkind,  3.9520_rkind ] ! Angstrom
    real(rkind), dimension(ns) :: eps_by_k = [   91.46_rkind,   200.0_rkind ] ! K
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
    real(rkind) :: L_y = 0.1143_rkind
    real(rkind) :: L_x = 16._rkind * L_y, L_z = L_y / 4._rkind
    real(rkind) :: x1 = 2.5146_rkind - Lx, y1 = zero, z1 = -Lz / 2._rkind

    logical :: periodicx = .false., periodicy = .false., periodicz = .true.

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use inclinedRM_data

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
        dy = L_yz/real(ny-1,rkind)
        dz = L_yz/real(nz  ,rkind)

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
    use ChapmanEnskogViscosityMod,   only: chapmanEnskogViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use ConstPrandtlConductivityMod, only: constPrandtlConductivity
    use ReidRamshawDiffusivityMod,   only: reidRamshawDiffusivity
    use miranda_restart_mod,         only: miranda_restart
    use exits,                       only: GracefulExit
    
    use inclinedRM_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    integer :: i, iounit
    real(rkind) :: Y_N2, Y_O2, rho_air, rho_N2, rho_O2, rho_HG
    real(rkind) :: R_air, Cp_air, gamma_air, Cp_N2, Cp_O2, sos_air, R_HG
    real(rkind) :: rho_shocked, u_shocked, p_shocked
    real(rkind), dimension(:,:,:,:), allocatable :: resmesh, resdata
    integer :: resdump = 0
    character(len=clen) :: charout

    type(chapmanEnskogViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(constPrandtlConductivity) :: thermcond

    type(miranda_restart) :: mir

    namelist /PROBINPUT/  resdir, resfile, resdump
    
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

        ! Initialize miranda_restart object
        call mir%init(decomp, resdir, resfile)
        if (mir%ns /= mix%ns) call GracefulExit("Number of species doesn't match that in the restart files",5687)

        ! Allocate restart mesh and data arrays
        allocate( resmesh(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3), 3       ) )
        allocate( resdata(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3), mir%nres) )

        ! Read in the grid
        call mir%read_grid(resmesh)

        ! Read in data
        call mir%read_data(resdump, resdata, tsim, dt)
        
        write(charout,'(A,I0.0,A,A)') "Reading Miranda restart dump ", resdump, " from ", trim(resdir)
        call message(charout)
        call message("Simulation time at restart", tsim)

        rho = resdata(:,:,:,mir%rho_index) * real(1.D3,  rkind) ! g/cm^3 to kg/m^3
        u   = resdata(:,:,:,  mir%u_index) * real(1.D-2, rkind) ! cm/s to m/s
        v   = resdata(:,:,:,  mir%v_index) * real(1.D-2, rkind) ! cm/s to m/s
        w   = resdata(:,:,:,  mir%w_index) * real(1.D-2, rkind) ! cm/s to m/s
        p   = resdata(:,:,:,  mir%p_index) * real(1.D-1, rkind) ! Ba (CGS) to Pa (SI)

        Ys  = resdata(:,:,:,mir%Ys_index:mir%Ys_index+mir%ns-1) ! Non-dimensional

        ! Deallocate temporary arrays and destroy miranda_restart object
        deallocate( resmesh )
        deallocate( resdata )
        call mir%destroy()
    end associate

end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture

    use inclinedRM_data

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

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/inclinedRM_", vizcount, ".dat"

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

subroutine hook_bc(decomp,mesh,fields,mix,tsim)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture

    use inclinedRM_data

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

    use inclinedRM_data

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

        ! write(outputfile,'(A,I4.4,A)') "inclinedRM_stats_N", ny, ".dat"
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

subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters, only: rkind
    use decomp_2d,       only: decomp_info
    use MixtureEOSMod,    only: mixture

    use inclinedRM_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

end subroutine
