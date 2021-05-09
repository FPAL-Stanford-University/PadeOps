module RTNCfinal_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, one
    use FiltersMod,       only: filters

    implicit none


    integer, parameter :: ns = 2
    logical :: twodim = .true.

    ! Problem parameters
    logical :: velInit = .false.              ! Use nonzero velocity initiation to limit pressure waves?
    integer :: prop_type = 0                  ! Transport property model to use
    real(rkind) :: trans_t = 5.0_rkind        ! Transport property transition time for prop_type=2
    real(rkind) :: trans_width = 0.2_rkind    ! Transport property transition thickness for prop_type=2
    real(rkind) :: T_ref = 75.0_rkind         ! Transport property reference temperature for prop_type=1,2
    real(rkind) :: M_Ratio = 1.5D0            ! Molecular mass ratio M2/M1
    real(rkind) :: At = 2D-1                  ! Atwood number
    real(rkind) :: z_interface = 5.5_rkind    ! Interface location
    real(rkind) :: rhoRatio = 3.0_rkind       ! Density ratio
    real(rkind) :: Re2 = 100.0_rkind          ! Reynolds number
    real(rkind) :: Sc = 1.0_rkind             ! Schmidt number
    real(rkind) :: Pr2 = 1.0_rkind            ! Prandtl number
    real(rkind) :: Sr = 0.04_rkind            ! Stratification Parameter
    real(rkind) :: gravity = 1.0_rkind        ! gravity
    real(rkind) :: L_int = 0.5_rkind          ! Interface thickness
    real(rkind) :: amp = 0.052_rkind          ! Amplitude of perturbation
    real(rkind) :: bulk_Ratio = 1.0_rkind     ! Constant ratio B1/B2
    real(rkind) :: mu_Ratio = 1.0_rkind       ! Constant ratio mu1/mu2
    real(rkind) :: kap_Ratio = 1.0_rkind      ! Constant ratio k1/k2
    real(rkind) :: p_back = zero              ! Background pressure to be added
    ! Parameters for the 2 materials:
    !                                        --------------------------------
    real(rkind), dimension(ns) :: gam      = [     1.66667_rkind,    1.66667_rkind ]
!    real(rkind), dimension(ns) :: molwt    = [ 28.0130_rkind,  44.010_rkind ] ! g/mol
    real(rkind), dimension(ns) :: Rgas

    ! Domain size data
    real(rkind) :: L_z = 5.0_rkind
    real(rkind) :: L_x = 6.0_rkind, L_y = 6.0_rkind
!    real(rkind) :: L_x = 1.0_rkind, L_y = 1.0_rkind
    real(rkind) :: x1, y1, z1

    ! Sponge region parameters
    real(rkind) :: gridLscale
    logical :: useSponge = .false.            ! Use sponge region at boundaries?
    real(rkind) :: spongesize = 0.1_rkind     ! Fraction of domain used as sponge region
    real(rkind) :: spongewidth = 0.1_rkind    ! Transition width for sponge region
    real(rkind) :: forcing_bot = 0.25_rkind   ! Forcing parameter for bottom sponge
    real(rkind) :: forcing_top = 0.25_rkind   ! Forcing parameter for top sponge
    real(rkind) :: forc_scaling = 0.25_rkind  ! Scaling factor for the forcing parameter
    real(rkind) :: diff_scaling = 1.0_rkind   ! Scaling factor for the diffusion parameter

    logical :: periodicx = .true., periodicy = .true., periodicz = .false.

    ! Gaussian filter for sponge
    type(filters) :: mygfil

    ! Perturbation parameters
    character(len=clen) :: x_pertfile, y_pertfile
    integer :: xmodes, ymodes
    real(rkind), dimension(:), allocatable :: kx, amp_x, phi_x
    real(rkind), dimension(:), allocatable :: ky, amp_y, phi_y
    real(rkind), dimension(:), allocatable :: pertz

contains

    subroutine read_perturbation_files()
        use mpi
        use decomp_2d, only: nrank
        use constants, only: zero, two, pi
        integer :: pertunit
        integer :: ierr
        integer :: i

        pertunit = 229

        ! Read in the x perturbation file
        if (nrank == 0) then
          print *, "Reading in the x-perturbation file ", trim(x_pertfile)
          open(unit=pertunit,file=trim(x_pertfile),form='FORMATTED',status='OLD')
          read(pertunit,*) xmodes
        end if

        ! broadcast number of xmodes to all procs
        call mpi_bcast(xmodes, 1, mpi_int, 0, mpi_comm_world, ierr)

        allocate(    kx(xmodes) )
        allocate( amp_x(xmodes) )
        allocate( phi_x(xmodes) )

        if (nrank == 0) then
          do i=1,xmodes
            read(pertunit,*) kx(i),amp_x(i),phi_x(i)
          end do
          close(pertunit)
        end if

        ! broadcast wavenumbers, amplitudes and phases to all procs
        call mpi_bcast(   kx, xmodes, mpirkind, 0, mpi_comm_world, ierr)
        call mpi_bcast(amp_x, xmodes, mpirkind, 0, mpi_comm_world, ierr)
        call mpi_bcast(phi_x, xmodes, mpirkind, 0, mpi_comm_world, ierr)

    end subroutine

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
    use constants,        only: zero, half, one, two
    use decomp_2d,        only: decomp_info

    use RTNCfinal_data

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

        if (twodim) then
            gridLscale = one/((one/dx)**two+(one/dz)**two)
        else
            gridLscale = one/((one/dx)**two+(one/dy)**two+(one/dz)**two)
        end if

    end associate

end subroutine


subroutine initfields(decomp,der,dx,dy,dz,inputfile,mesh,fields,Wbackground,mix,tsim,tstop,dt,tviz)
    use kind_parameters,             only: rkind, clen
    use constants,                   only: zero,half,one,two,pi,eight
    use CompressibleGridNC,          only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index,mom_index,TE_index
    use decomp_2d,                   only: decomp_info,transpose_x_to_y,transpose_y_to_x,transpose_y_to_z,transpose_z_to_y
    use DerivativesMod,              only: derivatives
    use MixtureEOSMod,               only: mixture
    use IdealGasEOS,                 only: idealgas
    use VariableViscosityMod,        only: variableViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use VariableConductivityMod,     only: variableConductivity
    use VariableDiffusivityMod,      only: variableDiffusivity
    use io_hdf5_stuff,               only: io_hdf5
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use mpi
    
    use RTNCfinal_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(derivatives),               intent(in)    :: der
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: Wbackground
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    integer :: i, j, k, l, mode, iounit

    type(variableViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(variableConductivity) :: thermcond

    real(rkind) :: A_minus, A_plus, AA_minus, AA_plus
    real(rkind) :: b_z, b_rho, b_p, H_minus_, H_plus_
    real(rkind), dimension(ns) :: mu_ref, bulk_ref, kap_ref

    real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xbuf1, xbuf2
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: tmp, tmprho, H_minus, H_plus, Rderiv
    real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: zbuf1, zbuf2

    integer :: nx, ny, nz, bc

    namelist /PROBINPUT/ useSponge, forc_scaling, diff_scaling, velInit, prop_type, trans_t, trans_width, T_ref, M_Ratio, At, Re2, Sc, Pr2, Sr, gravity, z_interface, L_int, amp, bulk_Ratio, mu_Ratio, kap_Ratio, p_back, x_pertfile
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,u_index),  &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,w_index),  &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,T_index),  &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:,mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap =>fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),       &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (mix%ns /= ns) call GracefulExit("Wrong number of species. Check your input file and make ns consistent with the problem file.",4562)

        ! First set the material properties and other quantities
        Rgas = [(1.0_rkind-1.0_rkind/gam(2))/M_Ratio, (1.0_rkind-1.0_rkind/gam(2))]
        A_plus = -Sr/(one+At)
        A_minus = -Sr/(one-At)
        AA_plus = A_plus*L_int/2.0_rkind
        AA_minus = A_minus*L_int/2.0_rkind

        ! Set each material's transport coefficient object
        mu_ref = [one*mu_Ratio/Re2, one/Re2]
        bulk_ref = [zero*bulk_Ratio, zero]
        kap_ref = [one*kap_Ratio/(Re2*Pr2), one/(Re2*Pr2)]
        do i = 1,mix%ns
            shearvisc = variableViscosity( mu_ref(i), T_ref, 2.5_rkind, trans_t, trans_width, prop_type )
            bulkvisc  = constRatioBulkViscosity( zero )
            thermcond = variableConductivity( kap_ref(i), T_ref, 2.5_rkind, trans_t, trans_width, prop_type)
            call mix%set_material( i, idealgas( gam(i), Rgas(i) ),&
                 shearvisc = shearvisc, &
                 bulkvisc  = bulkvisc, &
                 thermcond = thermcond  )
        end do

        ! Set mass diffusivity object (Ensure that all units are consistent)
        call mix%set_massdiffusivity( variableDiffusivity( one/(Re2*Sc), T_ref, 2.5_rkind, trans_t, trans_width, prop_type ) )

        call read_perturbation_files()

        tmp = z(:,:,:)-z_interface
        do i = 1,xmodes
          tmp = tmp - amp_x(i)*cos(kx(i)*two*pi/L_x*(x(:,:,:)-phi_x(i)))*exp(-two*pi/one*abs(z(:,:,:)-z_interface))
        end do

        do k = 1,decomp%ysz(3)
!            tmp(:,:,k) = z(:,:,k)-z_interface-amp*cos(two*pi/L_x*x(:,:,k))*exp(-two*pi/one*abs(z(:,:,k)-z_interface));
            H_plus(:,:,k) = (one+erf(tmp(:,:,k)/L_int))/2.0_rkind
            H_minus(:,:,k) = (one-erf(tmp(:,:,k)/L_int))/2.0_rkind
            rho(:,:,k) = (one+At)*exp(A_minus*tmp(:,:,k))*H_plus(:,:,k)+(one-At)*exp(A_plus*tmp(:,:,k))*H_minus(:,:,k)
        end do

        do k = 1,decomp%ysz(3)
            Ys(:,:,k,2) = (one+At)*exp(A_minus*tmp(:,:,k))*H_plus(:,:,k)/rho(:,:,k)
        end do
         
        Ys(:,:,:,1)  = one - Ys(:,:,:,2)
        call mix%update(Ys)

        do k = 1,decomp%ysz(3)
            p(:,:,k) = p_back+(1-At**2)*(exp(A_minus*tmp(:,:,k))*H_plus(:,:,k)+exp(A_plus*tmp(:,:,k))*H_minus(:,:,k) &
              -exp(AA_minus**2.0_rkind)*erf(tmp(:,:,k)/L_int-AA_minus)/2.0_rkind &
              +exp(AA_plus**2.0_rkind)*erf(tmp(:,:,k)/L_int-AA_plus)/2.0_rkind)/Sr
        end do

        T   = p/(rho*mix%Rgas)

        call mix%get_transport_properties(p, T, Ys, zero, mu, bulk, kap, diff)

        if (velInit) then
          tmp = gam(2)/(gam(2)-1)*log(mix%Rgas)
          call transpose_y_to_x(tmp,xbuf1,decomp)
          call der%ddx(xbuf1,xbuf2,0,0)
          call transpose_x_to_y(xbuf2,Rderiv,decomp)
          u   = diff(:,:,:,1)*Rderiv
          call der%ddy(tmp,Rderiv,0,0)
          v   = diff(:,:,:,1)*Rderiv
          call transpose_y_to_z(tmp,zbuf1,decomp)
          call der%ddz(zbuf1,zbuf2,0,0)
          call transpose_z_to_y(zbuf2,Rderiv,decomp)
          w   = diff(:,:,:,1)*Rderiv
        else
          u   = zero
          v   = zero
          w   = zero
        end if

        if (useSponge) then
            do i = 1,mix%ns
                Wbackground(:,:,:,i) = rho*Ys(:,:,:,i)
            end do
            Wbackground(:,:,:,mom_index  ) = rho*u
            Wbackground(:,:,:,mom_index+1) = rho*v
            Wbackground(:,:,:,mom_index+2) = rho*w
            Wbackground(:,:,:, TE_index  ) = rho*( e + half*(u*u + v*v + w*w ) )

            ! Calculates forcing coefficients
            ! bottom fluid speed of sound
            b_z = zero-z_interface
            H_plus_ = (one+erf(b_z/L_int))/two
            H_minus_ = (one-erf(b_z/L_int))/two
            b_rho = (one+At)*exp(A_minus*b_z)*H_plus_+(one-At)*exp(A_plus*b_z)*H_minus_
            b_p = p_back+(1-At**2)*(exp(A_minus*b_z)*H_plus_+exp(A_plus*b_z)*H_minus_ &
              -exp(AA_minus**two)*erf(b_z/L_int-AA_minus)/two &
              +exp(AA_plus**two)*erf(b_z/L_int-AA_plus)/two)/Sr
            forcing_bot = forc_scaling*sqrt(gam(2)*b_p/b_rho)
!            call message(2,"b_z=",b_z) 
!            call message(2,"b_rho=",b_rho)
!            call message(2,"b_p=",b_p)
            ! top fluid speed of sound
            b_z = L_z-z_interface
            H_plus_ = (one+erf(b_z/L_int))/two
            H_minus_ = (one-erf(b_z/L_int))/two
            b_rho = (one+At)*exp(A_minus*b_z)*H_plus_+(one-At)*exp(A_plus*b_z)*H_minus_
            b_p = p_back+(1-At**2)*(exp(A_minus*b_z)*H_plus_+exp(A_plus*b_z)*H_minus_ &
              -exp(AA_minus**two)*erf(b_z/L_int-AA_minus)/two &
              +exp(AA_plus**two)*erf(b_z/L_int-AA_plus)/two)/Sr
            forcing_top = forc_scaling*sqrt(gam(2)*b_p/b_rho)
!            call message(2,"b_z=",b_z)
!            call message(2,"b_rho=",b_rho)
!            call message(2,"b_p=",b_p)
        end if

    end associate

end subroutine


subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGridNC, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture

    use RTNCfinal_data

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

    end associate
end subroutine


subroutine hook_filmask(decomp,mesh,Wcnsrv,mask,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two, pi
    use CompressibleGridNC, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info

    use RTNCfinal_data

    implicit none
    type(decomp_info),               intent(in)      :: decomp
    real(rkind),                     intent(in)      :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)      :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)      :: Wcnsrv
    real(rkind), dimension(:,:,:),   intent(inout)   :: mask
    integer, dimension(2),           intent(in)      :: x_bc, y_bc, z_bc

    integer :: i
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: tmp
    
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        ! makes the mask function
        tmp = one-half*erf((z-L_z*spongesize)/spongewidth)+half*erf((z-L_z+L_z*spongesize)/spongewidth)
        where (tmp < 0.00001_rkind)
            mask = zero
        elsewhere
            mask = tmp
        end where
   end associate
end subroutine


subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two, pi
    use CompressibleGridNC, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use RTNCfinal_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    integer :: i
    real(rkind) :: dx
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (decomp%yst(3) == 1) then
!             Ys(:,:,1,1) = 4.0_rkind/3.0_rkind*Ys(:,:,2,1)-one/3.0_rkind*Ys(:,:,3,1)
!             Ys(:,:,1,2) = 4.0_rkind/3.0_rkind*Ys(:,:,2,2)-one/3.0_rkind*Ys(:,:,3,2)
!             u(:,:,1) = 4.0_rkind/3.0_rkind*u(:,:,2)-one/3.0_rkind*u(:,:,3)
!             T(:,:,1) = 4.0_rkind/3.0_rkind*T(:,:,2)-one/3.0_rkind*T(:,:,3)
             w(:,:,1) = zero
        end if


        ! zsz(3)=nz_global
        if (decomp%yen(3) == decomp%zsz(3)) then
!             Ys(:,:,decomp%yen(3),1) = 4.0_rkind/3.0_rkind*Ys(:,:,decomp%yen(3)-1,1)-one/3.0_rkind*Ys(:,:,decomp%yen(3)-2,1)
!             Ys(:,:,decomp%yen(3),2) = 4.0_rkind/3.0_rkind*Ys(:,:,decomp%yen(3)-1,2)-one/3.0_rkind*Ys(:,:,decomp%yen(3)-2,2)
!             u(:,:,decomp%yen(3)) = 4.0_rkind/3.0_rkind*u(:,:,decomp%yen(3)-1)-one/3.0_rkind*u(:,:,decomp%yen(3)-2)
!             T(:,:,decomp%yen(3)) = 4.0_rkind/3.0_rkind*T(:,:,decomp%yen(3)-1)-one/3.0_rkind*T(:,:,decomp%yen(3)-2)
             w(:,:,decomp%yen(3)) = zero
        end if


    end associate
end subroutine


subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,two,eps
    use CompressibleGridNC, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use exits,            only: message
    use reductions,       only: P_MAXVAL,P_MINVAL

    use RTNCfinal_data

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

        if (useSponge .AND. step==1) then
            call message(2,"Sponge region used with top forcing coefficient of",forcing_top)
            call message(2,"and bottom forcing coefficient of",forcing_bot)
            call message(2,"Diffusion length scale of",gridLscale)
        end if

        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message(2,"Maximum diffusivity",P_MAXVAL(diff))

    end associate
end subroutine


subroutine hook_source(decomp,der,mesh,fields,Wcnsrv,Wbackground,mix,tsim,dt,rhs,rhsg)
    use kind_parameters,    only: rkind
    use constants,          only: zero,half,one,two,four
    use decomp_2d,          only: decomp_info,transpose_x_to_y,transpose_y_to_x,transpose_y_to_z,transpose_z_to_y
    use DerivativesMod,     only: derivatives
    use MixtureEOSMod,      only: mixture
    use CompressibleGridNC, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index,mom_index,TE_index

    use RTNCfinal_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(derivatives),               intent(in)    :: der
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind),                     intent(in)    :: dt
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(in)    :: Wcnsrv
    real(rkind), dimension(:,:,:,:), intent(in)    :: Wbackground
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    real(rkind), dimension(:,:,:,:), optional, intent(inout) ::rhsg

    real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xbuf1, xbuf2
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: Wp, dWdxi2, tmp, maskfunc
    real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: zbuf1, zbuf2
    real(rkind) :: num_diff
    integer :: i

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

        if (useSponge) then
            ! makes the mask function
            tmp = one-half*erf((z-L_z*spongesize)/spongewidth)+half*erf((z-L_z+L_z*spongesize)/spongewidth)
            where (tmp < 0.00001_rkind)
                maskfunc = zero
            elsewhere
                maskfunc = tmp
            end where

            ! now tmp = forcing parameter
            where (z < z_interface)
                tmp = forcing_bot
            elsewhere
                tmp = forcing_top
            end where            

            num_diff = diff_scaling*gridLscale/(two*dt)

            do i = 1,TE_index
                ! add forcing term
                Wp = Wcnsrv(:,:,:,i)-Wbackground(:,:,:,i)
                rhs(:,:,:,i) =  rhs(:,:,:,i) - tmp*maskfunc*Wp
                ! add extra diffusion term
                call transpose_y_to_x(Wp,xbuf1,decomp)
                call der%d2dx2(xbuf1,xbuf2,0,0)
                call transpose_x_to_y(xbuf2,tmp,decomp)
                dWdxi2 = tmp
                call der%d2dy2(Wp,tmp,0,0)
                dWdxi2 = dWdxi2 + tmp
                call transpose_y_to_z(Wp,zbuf1,decomp)
                call der%d2dz2(zbuf1,zbuf2,0,0)
                call transpose_z_to_y(zbuf2,tmp,decomp)
                dWdxi2 = dWdxi2 + tmp
                rhs(:,:,:,i) =  rhs(:,:,:,i) + num_diff*maskfunc*dWdxi2
            end do

        end if

    end associate
end subroutine
