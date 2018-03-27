module inclinedRM_data
    use kind_parameters,  only: rkind, clen
    use constants,        only: zero, one
    use FiltersMod,       only: filters
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

    ! real(rkind), parameter :: R_univ = 8.3144598_rkind ! Universal gas constant in SI units
    real(rkind), parameter :: R_univ = 8.314472_rkind ! Universal gas constant in SI units

    ! Domain size data
    real(rkind) :: L_y = 0.1143_rkind
    real(rkind) :: L_x, L_z
    real(rkind) :: x1, y1, z1

    logical :: periodicx = .false., periodicy = .false., periodicz = .true.

    ! Miranda restart parameters
    character(len=clen) :: resdir, resfile
    integer :: resdump = 0

    ! Gaussian filter for sponge
    type(filters) :: mygfil

    ! Exact compatibility with Miranda?
    logical :: miranda_compat = .true.

contains

! ------------------------------------------------------------------------------
! Assign problem-specific viscosity, thermal conductivity,
! species diffusivities and magnetic diffusivity.
! ------------------------------------------------------------------------------ 
    SUBROUTINE prob_properties(rho,p,T,Y,mu,ktc,Diff)
     USE constants, ONLY : epssmall,half,three,two

     real(rkind), dimension(:,:,:),                                             intent(in)  :: rho,p,T
     real(rkind), dimension(:,:,:,:),                                           intent(in)  :: Y
     real(rkind), dimension(size(Y,1), size(Y,2), size(Y,3)),                   intent(out) :: mu,ktc    ! Shear visc and conductivity
     real(rkind), dimension(size(Y,1), size(Y,2), size(Y,3), size(Y,4)),        intent(out) :: Diff      ! Species diffusivity

     real(rkind), DIMENSION(size(Y,1), size(Y,2), size(Y,3)) :: M,CpMix,Rgas,Ts,Omega
     real(rkind), DIMENSION(size(Y,1), size(Y,2), size(Y,3)) :: mu_i,tmp,dum,Dij,Xii,Xjj
     real(rkind) :: Acc,Bcc,Ccc,Dcc,Ecc,Fcc,Gcc,Hcc
     real(rkind) :: Mij,sigij,Teij
     INTEGER :: n,ni,nj
  
    ! This part is from the Miranda problem file to recalculate the physical species diffusivity
    ! Material properties list ------ Property --------  Air ------------ CO2 --------
    real(rkind), DIMENSION(2) :: mySig    = (/ 3.681D0            , 3.952D0   /)
    real(rkind), DIMENSION(2) :: myTeff   = (/ 91.46D0            , 200.0D0   /)
    real(rkind), DIMENSION(2) :: myMolwts = (/ 28.013D0           , 44.010D0  /)
    real(rkind), DIMENSION(2) :: myGammas = (/ 1.4D0              , 1.28D0    /)
    real(rkind), DIMENSION(2) :: myPr     = (/ 0.72D0             , 0.77D0    /)
    real(rkind), DIMENSION(2) :: myMix    = (/ 1.0D0              , 1.0D0     /)
    real(rkind), DIMENSION(2) :: myNuA    = (/ 1.45D-5            , 1.59D-5   /) ! in g/(cm-s-K^0.5)
    real(rkind), DIMENSION(2) :: myNuB    = (/ 1.105D2            , 2.438D2   /) ! in K
    real(rkind), DIMENSION(2) :: myNuC    = (/ 0.0D0              , 0.0D0     /) ! in g/(cm-s)
    real(rkind), DIMENSION(2) :: myNuN    = (/ 1.5D0              , 1.5D0     /)

    real(rkind), PARAMETER :: Runiv  = 8.314472D+7 ! Gas constant              [g*cm^2/(s^2*mol*K)]=[erg/K/mol]

     ! Effective Molecular weight for EOS
     ! M = zero
     ! DO n = 1,ns
     !    M = M + Y(:,:,:,n)/myMolwts(n)
     ! END DO
     ! M = one/M
     ! Rgas = Runiv/M
     Rgas = zero
     DO n = 1,ns
        Rgas = Rgas + Y(:,:,:,n)*(Runiv/myMolwts(n))
     END DO
     M = Runiv/Rgas
   
   
     ! Chapman-Cowling Viscosity
     mu = zero
     ktc = zero
     tmp = zero
     Acc = 1.16145D0
     Bcc = -0.14874D0
     Ccc = 0.52487D0
     Dcc = -0.7732D0
     Ecc = 2.16178D0
     Fcc = -2.43787D0
     DO n=1,ns
         Ts = T/myTeff(n)
         Omega = Acc*Ts**Bcc + Ccc*exp(Dcc*Ts) + Ecc*exp(Fcc*Ts)     ! Edit: Ecc** -> Ecc*
         mu_i = 2.6693D-6*sqrt( myMolwts(n)*T)/(Omega*mySig(n)**2)
         dum = Y(:,:,:,n)/(myMolwts(n)**half)
         mu = mu + mu_i * dum
         tmp = tmp + dum
   
         ! Thermal diff part
         ! CpMix = myGammas(n)*Y(:,:,:,n)/(myMolwts(n)*(myGammas(n)-one))*Runiv
         CpMix = myGammas(n)/(myMolwts(n)*(myGammas(n)-one))*Runiv
         ktc = ktc + CpMix*mu_i*dum/myPr(n)
   
     END DO
     mu = mu/tmp * 10.0D0  ! convert kg/m s to g/cm s
     ktc = ktc/tmp * 10.0D0
   
     ! Diffusivity
     Diff = zero
     tmp = zero
     Acc = 1.06036D0
     Bcc = -0.1561D0
     Ccc = 0.19300D0
     Dcc = -0.47635D0
     Ecc = 1.03587D0
     Fcc = -1.52996D0
     Gcc = 1.76474D0
     Hcc = -3.89411D0
     DO ni=1,ns
         ! Xii = Y(:,:,:,ni)*M/myMolwts(ni)
         Xii = Y(:,:,:,ni)*(Runiv/myMolwts(ni))/Rgas
         tmp = zero
         DO nj=1,ns
            IF (ni .NE. nj) THEN
               ! Xjj = Y(:,:,:,nj)*M/myMolwts(nj)
               Xjj = Y(:,:,:,nj)*(Runiv/myMolwts(nj))/Rgas
               Mij = two/(one/myMolwts(ni)+one/myMolwts(nj))
               sigij = (mySig(ni)+mySig(nj))/two
               Teij = sqrt( myTeff(ni)*myTeff(nj) )
               Ts = T/Teij                                         ! Edit: Add omega calculation
               Omega = Acc*Ts**Bcc + Ccc*exp(Dcc*Ts) + Ecc*exp(Fcc*Ts) + Gcc*exp(Hcc*Ts)
               ! Dij = 0.0266/Omega * T**(three/two)/(p*sqrt(Mij)*sigij**2 )
               Dij = 0.0266D0/Omega * T**(three/two)/(p*sqrt(Mij)*sigij**2 )
               tmp = tmp + Xjj/Dij
            END IF
         END DO
         ! Diff(:,:,:,ni) = (one-Xii)/(tmp + 1.0D-16)
         Diff(:,:,:,ni) = (one-Xii)/(tmp + epssmall)
     END DO 
     ! Diff = Diff * 10.0D4 / 10.0D0 ! Convert m^2/s to cm^2/s  (extra /10 is from pressure)
     Diff = Diff * 1.0D4 * 10.0D0 ! Convert m^2/s to cm^2/s  (extra /10 is from pressure)
   
    END SUBROUTINE prob_properties


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

        L_x = 16._rkind * L_y
        L_z = L_y / 4._rkind

        dx = L_x/real(nx-1,rkind)
        dy = L_y/real(ny-1,rkind)
        dz = L_z/real(nz  ,rkind)

        if (miranda_compat) then
            dx = dy
            dz = dy
            L_z = dz * real(nz-1, rkind)
            L_x = dx * real(nx-1, rkind)
        end if

        x1 = 2.5146_rkind - L_x
        y1 = zero
        z1 = -L_z / 2._rkind

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
    use io_hdf5_stuff,               only: io_hdf5
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use mpi
    
    use inclinedRM_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    integer :: i, j, k, l, iounit
    real(rkind) :: Y_N2, Y_O2, rho_air, rho_N2, rho_O2, rho_HG
    real(rkind) :: R_air, Cp_air, gamma_air, Cp_N2, Cp_O2, sos_air, R_HG
    real(rkind) :: rho_shocked, u_shocked, p_shocked
    real(rkind), dimension(:,:,:,:), allocatable :: resmesh, resdata
    character(len=clen) :: charout

    type(chapmanEnskogViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(constPrandtlConductivity) :: thermcond

    type(miranda_restart) :: mir
    ! type(io_hdf5)         :: viz

    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: mu, bulk, kappa
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3), mix%ns) :: diff

    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: mu_o, bulk_o, kappa_o
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3), mix%ns) :: diff_o

    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: tmp

    integer :: nx, ny, nz

    namelist /PROBINPUT/  resdir, resfile, resdump
    
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

        call message("Error in dx", abs(dx-(resmesh(2,1,1,1)-resmesh(1,1,1,1))*real(1.D-2,rkind)))
        call message("Error in dy", abs(dy-(resmesh(1,2,1,2)-resmesh(1,1,1,2))*real(1.D-2,rkind)))
        call message("Error in dz", abs(dz-(resmesh(1,1,2,3)-resmesh(1,1,1,3))*real(1.D-2,rkind)))

        call message("Max error in x coordinate",P_MAXVAL( abs(x - resmesh(:,:,:,1)*real(1.D-2,rkind))/dx ))
        call message("Min error in x coordinate",P_MINVAL( abs(x - resmesh(:,:,:,1)*real(1.D-2,rkind))/dx ))

        call message("Max error in y coordinate",P_MAXVAL( abs(y - resmesh(:,:,:,2)*real(1.D-2,rkind))/dy ))
        call message("Min error in y coordinate",P_MINVAL( abs(y - resmesh(:,:,:,2)*real(1.D-2,rkind))/dy ))

        call message("Max error in z coordinate",P_MAXVAL( abs(z - resmesh(:,:,:,3)*real(1.D-2,rkind))/dz ))
        call message("Min error in z coordinate",P_MINVAL( abs(z - resmesh(:,:,:,3)*real(1.D-2,rkind))/dz ))

        print *, "x[ 1,1,1] = ", x( 1,1,1), "    resmesh[ 1,1,1,1] = ", resmesh( 1,1,1,1)
        print *, "x[nx,1,1] = ", x(nx,1,1), "    resmesh[nx,1,1,1] = ", resmesh(nx,1,1,1)

        ! Read in data
        call mir%read_data(resdump, resdata, tsim, dt)
        
        write(charout,'(A,I0.0,A,A)') "Reading Miranda restart dump ", resdump, " from ", trim(resdir)
        call message(charout)
        call message("Simulation time at restart", tsim)
        call message("Timestep at restart", dt)

        rho = resdata(:,:,:,mir%rho_index) * real(1.D3,  rkind) ! g/cm^3 to kg/m^3
        u   = resdata(:,:,:,  mir%u_index) * real(1.D-2, rkind) ! cm/s to m/s
        v   = resdata(:,:,:,  mir%v_index) * real(1.D-2, rkind) ! cm/s to m/s
        w   = resdata(:,:,:,  mir%w_index) * real(1.D-2, rkind) ! cm/s to m/s
        p   = resdata(:,:,:,  mir%p_index) * real(1.D-1, rkind) ! Ba (CGS) to Pa (SI)

        Ys  = resdata(:,:,:,mir%Ys_index:mir%Ys_index+mir%ns-1) ! Non-dimensional
        
        if ( nancheck(rho) ) then
            print '(A)', "NaN found in rho"
        end if
        if ( nancheck(u) ) then
            print '(A)', "NaN found in u"
        end if
        if ( nancheck(v) ) then
            print '(A)', "NaN found in v"
        end if
        if ( nancheck(w) ) then
            print '(A)', "NaN found in w"
        end if
        if ( nancheck(p) ) then
            print '(A)', "NaN found in p"
        end if
        if ( nancheck(Ys,i,j,k,l) ) then
            print '(A,4(I0.0,A))', "NaN found at Ys(", i, ",", j, ",", k, ",", l, ")"
        end if
        ! Initialize viz object
        ! call viz%init( mpi_comm_world, decomp, 'y', '.', 'init_test', &
        !                reduce_precision=.true., read_only=.false., jump_to_last=.true.)
        ! call viz%write_coords(mesh)

        ! Initialize mygfil
        call mygfil%init(                           decomp, &
                          periodicx,  periodicy, periodicz, &
                         "gaussian", "gaussian", "gaussian" )

        ! Check for consistency
        call mix%update(Ys)
        call mix%get_e_from_p(rho,p,e)
        call mix%get_T(e, T)

        call message("Max error in e", P_MAXVAL( abs(e - resdata(:,:,:,mir%e_index)*real(1.D-4,rkind)) ))
        call message("Max error in T", P_MAXVAL( abs(T - resdata(:,:,:,mir%T_index)) ))

        call mix%get_transport_properties(p, T, Ys, mu, bulk, kappa, diff)
        call prob_properties( resdata(:,:,:,mir%rho_index), resdata(:,:,:,mir%p_index), &
                              resdata(:,:,:,  mir%T_index), Ys, &
                              mu_o, kappa_o, diff_o)
        bulk_o = zero

        call message("Max error in mu", P_MAXVAL( abs(mu - mu_o*real(1.D-1,rkind)) ))

        call message("Max error in bulk", P_MAXVAL( abs(bulk - bulk_o*real(1.D-1,rkind)) ))
        call message("Max error in kappa", P_MAXVAL( abs(kappa - kappa_o*real(1.D-5,rkind)) ))

        tmp = abs(diff(:,:,:,1) - diff_o(:,:,:,1)*real(1.D-4,rkind))
        where ( Ys(:,:,:,1) > one - real(1.D-5,rkind) )
            tmp = zero
        end where
        call message("Max error in diff(1)", P_MAXVAL(tmp))

        tmp = abs(diff(:,:,:,2) - diff_o(:,:,:,2)*real(1.D-4,rkind))
        where ( Ys(:,:,:,2) > one - real(1.D-5,rkind) )
            tmp = zero
        end where
        call message("Max error in diff(2)", P_MAXVAL(tmp))

        ! call viz%start_viz(tsim)
        ! call viz%write_variable(rho , 'rho' ) 
        ! call viz%write_variable(u   , 'u') 
        ! call viz%write_variable(v   , 'v') 
        ! call viz%write_variable(w   , 'w') 
        ! call viz%write_variable(p   , 'p') 

        ! call viz%write_variable(Ys(:,:,:,1), 'Ys_1') 
        ! call viz%write_variable(Ys(:,:,:,2), 'Ys_2') 

        ! call viz%write_variable(e, 'e') 
        ! call viz%write_variable(resdata(:,:,:,mir%e_index)*real(1.D-4,rkind), 'e_res') 

        ! call viz%write_variable(T, 'T') 
        ! call viz%write_variable(resdata(:,:,:,mir%T_index), 'T_res') 

        ! call viz%write_variable(mu, 'mu') 
        ! call viz%write_variable(mu_o*real(1.D-1,rkind), 'mu_res') 

        ! call viz%write_variable(kappa, 'kappa') 
        ! call viz%write_variable(kappa_o*real(1.D-5,rkind), 'kappa_res') 

        ! call viz%write_variable(diff(:,:,:,1), 'diff_1') 
        ! call viz%write_variable(diff_o(:,:,:,1)*real(1.D-4,rkind), 'diff_res_1') 

        ! call viz%write_variable(diff(:,:,:,2), 'diff_2') 
        ! call viz%write_variable(diff_o(:,:,:,2)*real(1.D-4,rkind), 'diff_res_2') 

        ! call viz%end_viz()

        ! Deallocate temporary arrays and destroy miranda_restart object
        deallocate( resmesh )
        deallocate( resdata )
        call mir%destroy()
        ! call viz%destroy()

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

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use inclinedRM_data

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
        ! Sponge+bulk for exit bc
        ! Gradually apply the exit boundary conditions  
        dx = L_x/real(decomp%xsz(1)-1,rkind)
        filpt = one*real(1.D-2,rkind) / dx ! 1cm width for sponge
        thickT = real(10.D0, rkind)
        do i=1,decomp%ysz(1)
            dumT(i,:,:)=half*(one-tanh( (real( decomp%yst(1) - 1 + i - 1, rkind)-filpt) / thickT ))
        end do
            
        !!!  OUT-FLOW  !!!
        if (decomp%yst(1) == 1) then
            v(1:5,:,:)  = zero      ! Only allow for normal velocities
            w(1:5,:,:)  = zero      ! Only allow for normal velocities
            Ys(1,:,:,1) = one       ! This boundary only has air
            Ys(1,:,:,2) = zero      ! This boundary only has air
        end if   
            
        ! Gussian Filter for last N points in x-direction to act as a sponge (A)
        do i=1,4
            dumF = u
            call filter3D(decomp,mygfil,dumF,1,x_bc,y_bc,z_bc)
            u = u + dumT*(dumF-u) 

            dumF = v
            call filter3D(decomp,mygfil,dumF,1,x_bc,y_bc,z_bc)
            v = v + dumT*(dumF-v)

            dumF = w
            call filter3D(decomp,mygfil,dumF,1,x_bc,y_bc,z_bc)
            w = w + dumT*(dumF-w)

            dumF = p
            call filter3D(decomp,mygfil,dumF,1,x_bc,y_bc,z_bc)
            p = p + dumT*(dumF-p)

            dumF = rho
            call filter3D(decomp,mygfil,dumF,1,x_bc,y_bc,z_bc)
            rho = rho + dumT*(dumF-rho)
        end do

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
