module ShearLayer_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, half, one, two, four, pi, imi
    use FiltersMod,       only: filters
    use decomp_2d,        only: decomp_info, nrank
    use basic_io,         only: read_2d_ascii 
    use reductions,       only: P_MAXVAL,P_MINVAL
    use exits,            only: message
    use mpi
    implicit none
    integer, parameter :: ns = 2

    ! Problem parameters
    real(rkind) :: Mc = 0.7_rkind           ! Convective Mach
    real(rkind) :: Re = 800._rkind          ! Reynolds number
    real(rkind) :: Sc = 1._rkind            ! Schmidt number
    real(rkind) :: p_ref = one              ! reference press
    real(rkind) :: T_ref = one              ! reference temp
    real(rkind) :: rho_ref = one            ! reference density
    real(rkind) :: rho_ratio = one          ! rho2/rho1
    real(rkind) :: dtheta0 = 1._rkind       ! Base profile thickness 
    real(rkind) :: noiseAmp = 1D-6          ! white noise amplitude
    character(len=clen) :: fname_prefix
    logical :: use_lstab = .false. 
    logical :: no_pert = .false. 
    
    ! Parameters for the 2 materials
    real(rkind):: gam=1.4_rkind
    real(rkind):: Pr=0.7_rkind 
    real(rkind), dimension(ns) :: Rgas=[one,one]

    ! Domain size data
    real(rkind) :: Ly = 129._rkind, Lx=172._rkind, Lz=86._rkind
    real(rkind) :: x1, y1, z1
    logical :: periodicx = .true., periodicy = .false., periodicz = .true.

    ! Gaussian filter for sponge
    type(filters) :: mygfil
contains
    subroutine perturb_potential_v2(gp,x,y,z,Lx,Lz,u,v,w)
        use decomp_2d,        only: nrank
        type(decomp_info), intent(in)               :: gp
        real(rkind), dimension(:,:,:), intent(in)   :: x,y,z
        real(rkind), dimension(:,:,:), intent(inout):: u,v,w
        real(rkind), intent(in)                     :: Lx,Lz

        real, dimension(:), allocatable :: A,e,x1d,y1d,z1d
        real :: kx, kz, phx, phz, kmin, Amax, du
        real :: sigma=10,pi = 3.14159265358
        integer :: ni,nj,nk,i, j, k, m, mx, mz, mode_min=4, nmodes_max=4, &
            nxmodes, nzmodes, mpi_ierr, ix1,iz1,ixn,izn
   
        ! Global size
        ni  = gp%xsz(1)
        nj  = gp%ysz(2)
        nk  = gp%zsz(3)

        ! If base decomposition is in Y
        ix1 = gp%yst(1); 
        iz1 = gp%yst(3)
        ixn = gp%yen(1); 
        izn = gp%yen(3)

        ! Store radius, theta, z 
        allocate(x1D(ni),y1D(nj),z1D(nk))
        x1d = x(:,1,1)
        y1D = y(1,:,1) 
        z1d = z(1,1,:)

        ! Transverse mask 
        allocate(e(nj),A(nj))
        e = exp(-sigma*(y1D**2));
        du = P_MAXVAL(u)-P_MINVAL(u)
        Amax = 0.01*(du)
        call mpi_bcast(Amax,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
        A = Amax*e

        call message(0,"Making perturbations")
        call message(2,"Maximum u", P_MAXVAL(abs(u)))
        call message(2,"Maximum v", P_MAXVAL(abs(v)))
        call message(2,"Maximum w", P_MAXVAL(abs(w)))
    
        nxmodes = min(nmodes_max, ni/4)
        nzmodes = min(nmodes_max, nk/4)
        
        do mx=mode_min,mode_min+nxmodes
        do mz=mode_min,mode_min+nzmodes
            kx = 2.D0*pi/Lx * mx
            kz = 2.D0*pi/Lz * mz
            call random_number(phx)
            call random_number(phz)
            call mpi_bcast(phx,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
            call mpi_bcast(phz,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
            phx = phx*2.D0*pi
            phz = phz*2.D0*pi
            do i=ix1,ixn
            do k=iz1,izn
            u(i,:,k) = u(i,:,k) + A*sin(kx*x1D(i)+phx)*sin(kz*z1D(k)+phz)
            w(i,:,k) = w(i,:,k) + A*cos(kx*x1D(i)+phx)*cos(kz*z1D(k)+phz)
            v(i,:,k) = v(i,:,k) + A/kz * sin(kx*x1D(i)+phx) * &
                ( (-2d0*sigma*y1d)*cos(kz*z1D(k)+phz) + kx*sin(kz*z1D(k)+phz) )
            enddo
            enddo
        enddo
        enddo
        deallocate(A,e,x1d,y1d,z1d)
        
        call message(0,"Done making perturbations")
        call message(2,"Maximum u", P_MAXVAL(abs(u)))
        call message(2,"Maximum v", P_MAXVAL(abs(v)))
        call message(2,"Maximum w", P_MAXVAL(abs(w)))
    end subroutine 

end module


subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info, nrank
    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    character(clen) :: inputfile='input.dat'

    namelist /PROBINPUT/ Lx, Ly, Lz,Mc, Re, Pr, Sc,&
                        T_ref, p_ref, rho_ref, rho_ratio,&
                        noiseAmp, fname_prefix, use_lstab, no_pert
    ioUnit = 11
    open(unit=ioUnit, file=inputfile, form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)
    
    ! Global domain size 
    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)

    ! Need to set x, y and z as well as  dx, dy and dz
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        if (nrank == 0) then
            print *, "Domain size: ",Lx,Ly,Lz
        end if
        dx = Lx/real(nx,rkind)
        dy = Ly/real(ny,rkind)
        dz = Lz/real(nz,rkind)
        x1 = 0._rkind! -Lx/2._rkind
        y1 = -Ly/2._rkind
        z1 = 0._rkind!-Lz/2._rkind
        
        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = x1 + real( ix1 + i - 1, rkind ) * dx
                    y(i,j,k) = y1 + real( iy1 + j - 1, rkind ) * dy
                    z(i,j,k) = z1 + real( iz1 + k - 1, rkind ) * dz
                end do
            end do
        end do
    end associate
end subroutine


subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
    use kind_parameters,             only: rkind, clen
    use constants,                   only: zero,half,one,two,four,five,pi,eight
    use CompressibleGrid,            only: rho_index,u_index,v_index,w_index,&
                                           p_index,T_index,e_index,Ys_index
    use decomp_2d,                   only: decomp_info,nrank
    use MixtureEOSMod,               only: mixture
    use IdealGasEOS,                 only: idealgas
    use PowerLawViscosityMod,        only: powerLawViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use ConstPrandtlConductivityMod, only: constPrandtlConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use random,                      only: gaussian_random                  

    use ShearLayer_data
    use mpi

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    type(powerLawViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(constPrandtlConductivity) :: thermcond
    !real(rkind), dimension(:,:), allocatable :: tmp2D
    real(rkind) :: mu_ref, c1,c2,du, Rgas1, Rgas2,lambda,kx, kz, ph 
    integer :: i,iounit, nx, ny, nz
    
    namelist /PROBINPUT/ Lx, Ly, Lz,Mc, Re, Pr, Sc,&
                        T_ref, p_ref, rho_ratio, rho_ref, &
                        noiseAmp, fname_prefix, use_lstab, no_pert 
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    ! Local domain sizes
    nx = decomp%xsz(1) 
    ny = decomp%ysz(2) 
    nz = decomp%zsz(3)

    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),&
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),&
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),&
                 e => fields(:,:,:,  e_index), &
                Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (mix%ns /= ns) call GracefulExit("Wrong number of species. Check your input file and make ns consistent with the problem file.",4562)
        
        ! First set the materials for rho_ratio = rho2/rho1
        Rgas1 = (one+rho_ratio)/two
        Rgas2 = (one+rho_ratio)/(two*rho_ratio)
        Rgas = [Rgas1, Rgas2]
        c1 = sqrt(gam*p_ref/(rho_ref/Rgas1))
        c2 = sqrt(gam*p_ref/(rho_ref/Rgas2))
        du = Mc*(c1+c2)
        mu_ref = one/Re

        ! Set each material's transport coefficient object
        do i = 1,mix%ns
            shearvisc = powerLawViscosity( mu_ref, T_ref, 0._rkind)
            bulkvisc  = constRatioBulkViscosity( zero )
            thermcond = constPrandtlConductivity( Pr )
            call mix%set_material( i, idealgas( gam, Rgas(i) ),&
                 shearvisc = shearvisc, & 
                 bulkvisc  = bulkvisc, &
                 thermcond = thermcond  )
        end do

        ! Set mass diffusivity object (Ensure that all units are consistent)
        lambda = (rho_ratio-1)/(rho_ratio+1)
        call mix%set_massdiffusivity( constSchmidtDiffusivity( mu_ref,rho_ref,Sc))
        Ys(:,:,:,1)  = one - half*(one+lambda*tanh(y/(two*dtheta0)))
        Ys(:,:,:,2)  = one - Ys(:,:,:,1)
        call mix%update(Ys)

        ! Add base flow profiles
        lambda = (rho_ratio-1)/(rho_ratio+1)
        u = u + du*half*tanh(y/(two*dtheta0))
        v = v + zero
        w = w + zero
        p = p + p_ref
        rho = rho + rho_ref*(1 + lambda*tanh(y/(two*dtheta0)))
        T = p/(rho*mix%Rgas) 
		
        ! Perturbations: this must be specific for each problem.
        if (use_lstab) then
            call GracefulExit("LSTAB init is deprecated.",4562)
            !call lstab_pert(decomp,x,z,fname_prefix,u,2)
            !call lstab_pert(decomp,x,z,fname_prefix,v,3)
            !call lstab_pert(decomp,x,z,fname_prefix,w,4)
            !call lstab_pert(decomp,x,z,fname_prefix,rho,1)
        else if (no_pert) then
            call message(0,"No perturbations")
        else
            call perturb_potential_v2(decomp,x,y,z,Lx,Lz,u,v,w)
        endif

        ! Initialize gaussian filter mygfil
        call mygfil%init(decomp, periodicx, periodicy, periodicz, &
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
    use reductions,       only: P_MEAN
    use ShearLayer_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(mixture),                   intent(in) :: mix
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) ::&
    tke,ubase
    real(rkind) :: tke_mean,tke0
    character(len=clen) :: outputfile,str
    integer :: i,outputunit=229

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ShearLayer_", &
            vizcount, ".dat"

        ! Get TKE
        !ubase = 0.5*(1+tanh(y))
        !tke = half*rho*((u-ubase)**2 + v*v + w*w)

        !write(str,'(I4.4)') decomp%ysz(2)
        !write(outputfile,'(2A)') trim(outputdir),"/ShearLayer_"//trim(str)//".dat"

        !if (vizcount == 0) then
        !    ! On the first step, get initial disturbance energy and
        !    ! write the header for the output file.
        !    tke0 = P_MEAN( tke )
        !    if (nrank == 0) then
        !        open(unit=outputunit, file=trim(outputfile), &
        !            form='FORMATTED', status='REPLACE')
        !        write(outputunit,'(3A26)') "Time", "TKE"
        !    end if
        !else
        !    ! Open the previously started file
        !    if (nrank == 0) then
        !        open(unit=outputunit, file=trim(outputfile), &
        !            form='FORMATTED', position='APPEND', status='OLD')
        !    end if
        !end if

        !tke_mean = P_MEAN(tke)
        !if (nrank == 0) then
        !    write(outputunit,'(3ES26.16)') tsim, tke_mean/tke0
        !    close(outputunit)
        !end if
    end associate
end subroutine


subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    integer :: i
    real(rkind) :: dy, filpt, thickT
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
        dy = Ly/real(decomp%ysz(2)-1,rkind)
        filpt = 5.0_rkind/dy 
        thickT = real(5.D0, rkind)
        
        ! Gussian Filter for 
        do i=1,decomp%ysz(2)
            dumT(:,i,:)=half*(one-tanh( (real( decomp%yst(2) - 1 + i - 1, rkind)-filpt) / thickT ))
        end do
            
        dumF = u
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        u = u + dumT*(dumF-u) 

        dumF = v
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        v = v + dumT*(dumF-v)

        dumF = w
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        w = w + dumT*(dumF-w)

        dumF = p
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        p = p + dumT*(dumF-p)

        dumF = rho
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        rho = rho + dumT*(dumF-rho)

        ! Gussian Filter for the top
        do i=1,decomp%ysz(2)
            dumT(:,i,:)=half*(one-tanh( (real(decomp%ysz(2)- (decomp%yst(2) - 1 + i - 1), rkind)-filpt) / thickT ))
        end do

        dumF = u
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        u = u + dumT*(dumF-u) 

        dumF = v
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        v = v + dumT*(dumF-v)

        dumF = w
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        w = w + dumT*(dumF-w)

        dumF = p
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        p = p + dumT*(dumF-p)

        dumF = rho
        call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
        rho = rho + dumT*(dumF-rho)

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

    use ShearLayer_data

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
    integer :: ETDNS=0

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        call message(2,"Maximum v-velocity",P_MAXVAL(v))
        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message(2,"Maximum diffusivity",P_MAXVAL(diff))

    end associate

end subroutine


subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs,Wcnsrv,tkeb,tsim_0,dtheta_0)
    use CompressibleGrid,   only: rho_index,u_index,v_index,w_index,&
                                  p_index,T_index,e_index,Ys_index
    use kind_parameters,    only: rkind
    use decomp_2d,          only: decomp_info
    use MixtureEOSMod,      only: mixture
    use TKEBudgetMod,       only: tkeBudget
    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

    real(rkind), dimension(:,:,:,:), optional,intent(in)    :: Wcnsrv
    type(tkeBudget), optional,       intent(inout) :: tkeb
    real(rkind), optional,           intent(in)    :: tsim_0 ! the previous time 
    real(rkind), optional,           intent(in)    :: dtheta_0 ! the previous L99 
    real(rkind) :: factor=0.d0, dtheta=0.d0, rate=0.d0
    integer :: i,mpi_ierr
    integer :: mass_index, mom_index, TE_index
    logical :: forcing
    
    if (present(tkeb)) then
        forcing=.true.
    else 
        forcing = .false.
    endif
    if (forcing) then
    
    call tkeb%get_dtheta(decomp,mesh(:,:,:,2),fields(:,:,:,rho_index),&
        fields(:,:,:,u_index),dtheta,rate,fields(:,:,:,v_index))
    factor  = rate/dtheta
    call mpi_bcast(factor,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
    call message(3,'Factor',factor)
    
    ! Set mass, momentum and energy indices in Wcnsrv
    mass_index = 1
    mom_index  = mass_index + ns
    TE_index   = mom_index + 3
    do i = 1,ns
        rhs(:,:,:,i) = rhs(:,:,:,i) - factor*Wcnsrv(:,:,:,Ys_index+i) 
    enddo 
    rhs(:,:,:,mom_index  ) = rhs(:,:,:,mom_index  ) - factor*Wcnsrv(:,:,:,mom_index  )
    rhs(:,:,:,mom_index+1) = rhs(:,:,:,mom_index+1) - factor*Wcnsrv(:,:,:,mom_index+1)
    rhs(:,:,:,mom_index+2) = rhs(:,:,:,mom_index+2) - factor*Wcnsrv(:,:,:,mom_index+2)
    rhs(:,:,:,TE_index   ) = rhs(:,:,:,TE_index   ) - factor*Wcnsrv(:,:,:,TE_index   )

    endif
end subroutine
