module ShearLayer_sCO2_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, half, one, two, four, pi, imi
    use FiltersMod,       only: filters
    use decomp_2d,        only: decomp_info, nrank
    use basic_io,         only: read_2d_ascii 
    use exits,            only: message
    use mpi
    implicit none
    integer, parameter :: ns = 1

    ! Problem parameters
    real(rkind) :: dtheta0 = 1._rkind       ! Base profile thickness 
    real(rkind) :: MW, Tc, Pc, omega, T0, rho0, noiseAmp
    real(rkind) :: mu_ref, T_ref
    real(rkind) :: U1, T1, P1, rho1 ! rho1 and rho2 are guesses
    real(rkind) :: U2, T2, P2, rho2 ! to be used for Newton Raphson iteration
    real(rkind) :: Pr, Sc, Re
    real(rkind) :: Ly, Lx, Lz
    real(rkind) :: x1, y1, z1
    logical :: periodicx = .true., periodicy = .false., periodicz = .true.
    integer(rkind) :: mpi_ierr

    ! Gaussian filter for sponge
    type(filters) :: mygfil
contains
    subroutine perturb_potential(x,y,z,Lx,Lz,u,v,w)
        !type(decomp_info), intent(in)               :: gp
        real(rkind), dimension(:,:,:), intent(in)   :: x,y,z
        real(rkind), dimension(:,:,:), intent(inout):: u,v,w 
        real(rkind), intent(in)                     :: Lx,Lz

        real(rkind) :: kx, kz, phx=0, phz=0, A, eps=0.15, sigma=5.D0
        integer(rkind) :: mx, mz, nmodes=10

        do mx = 4, 4+nmodes 
        do mz = 4, 4+nmodes

            kx = two*pi/Lx * mx 
            kz = two*pi/Lz * mz
            !call random_number(phx)
            !call random_number(phz)
            !call mpi_bcast(phx,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
            !call mpi_bcast(phz,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)

            A = eps/(four*pi**2*kx*kz)
            u = u + A * kx*sin(kx*x+phx*two*pi) * cos(kz*z+phz*two*pi) *&
                exp(-(sigma*y**2))*( sin(y) + 2*sigma*y*cos(y) )
            w = w + A * cos(kx*x+phx*two*pi) * kz*sin(kz*z+phz*two*pi) *&
                exp(-(sigma*y**2))*( sin(y) + 2*sigma*y*cos(y) )
            v = v + A * cos(kx*x+phx*two*pi) * kz*sin(kz*z+phz*two*pi) *&
                ( -two*sigma*y*exp(-(sigma*y**2))*( sin(y) + 2*sigma*y*cos(y) ) + &
                exp(-(sigma*y**2))*( cos(y) - 2*sigma*y*sin(y) ) )

        enddo
        enddo
    end subroutine 
end module


subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info, nrank
    use ShearLayer_sCO2_data

    implicit none
    character(clen)                                :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    namelist /PROBINPUT/ & 
        dtheta0,&
        MW, Tc, Pc, omega, mu_ref, T_ref, noiseAmp,&
        U1, T1, P1, rho1,&
        U2, T2, P2, rho2,&
        Pr, Sc, Re, &
        Ly, Lx, Lz
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
    use MixtureEOSMod_Real,          only: mixture
    use RealGasEOS,                  only: realgas
    use PowerLawViscosityMod,        only: powerLawViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use ConstPrandtlConductivityMod, only: constPrandtlConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use random,                      only: gaussian_random                  

    use ShearLayer_sCO2_data
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
    
    real(rkind), dimension(:), allocatable :: tmp1,tmp2,p1D,T1D,rho1D,alpha,f,df
    real(rkind) :: A,B,c,R,Ti,rhoi
    integer :: i,j, iounit, nx, ny, nz
    
    namelist /PROBINPUT/ & 
        dtheta0,&
        MW, Tc, Pc, omega, mu_ref, T_ref, noiseAmp,&
        U1, T1, P1, rho1,&
        U2, T2, P2, rho2,&
        Pr, Sc, Re, &
        Ly, Lx, Lz
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


    ! Set each material's transport coefficient object.
    ! Also initialize guesses for Ti,rhoi used in get_T_from_e
    if (mix%ns /= ns) call GracefulExit("Wrong number of species.",4562)
    Ti = 0.5*(T1+T2)
    rhoi = 0.5*(rho1+rho2)
    do i = 1,mix%ns
        shearvisc = powerLawViscosity( mu_ref, T_ref, 0.7)
        bulkvisc  = constRatioBulkViscosity( zero )
        thermcond = constPrandtlConductivity( Pr ) ! TODO: check if this is ok
        call mix%set_material( i, realgas(decomp,MW,Tc,Pc,omega,Ti,rhoi),&
             shearvisc = shearvisc, & 
             bulkvisc  = bulkvisc, &
             thermcond = thermcond  )
    end do

    ! Add base flow profiles.
    ! Set P and T, then iterate for consistent rho 
    if (allocated(tmp1)) deallocate(tmp1); allocate(tmp1(ny))
    if (allocated(tmp2)) deallocate(tmp2); allocate(tmp2(ny))
    if (allocated(p1D))  deallocate(p1D);  allocate(p1D(ny))
    if (allocated(T1D))  deallocate(T1D);  allocate(T1D(ny))
    if (allocated(rho1D)) deallocate(rho1D); allocate(rho1D(ny))
    if (allocated(alpha))   deallocate(alpha);  allocate(alpha(ny))
    if (allocated(f))       deallocate(f);      allocate(f(ny))
    if (allocated(df))      deallocate(df);     allocate(df(ny))
    
    tmp1 = tanh(y(0,:,0)/dtheta0)
    tmp2 = 1.d0-tmp1 
    do j=1,ny 
        u(:,j,:) = U1*tmp1(j) + U2*tmp2(j)
    enddo
    v = zero
    w = zero
    p1D = P1*tmp1 + P2*tmp2 
    T1D = T1*tmp1 + T2*tmp2 
    rho1D = rho1D + rho1*tmp1 + rho2*tmp2

    ! Iteratre to get density for corresponding P1,P2,T1,T2
    ! tmp1,tmp2 are PR denominators
    ! tmp1 = D1 = vm-B
    ! tmp2 = D2 = vm^2 - 2*vm*B - B^2
    A = mix%material(1)%mat%PR_B
    B = mix%material(1)%mat%PR_B
    c = mix%material(1)%mat%PR_kappa
    R = mix%material(1)%mat%Rgas
    do i=1,100
        tmp1 = 1d0/rho1D-B
        tmp2 = rho1D**-2 + 2*B/rho1D - B**2
        alpha = (1 + c*(1-(T1D/Tc)**0.5))**2 
        f = R*T1D/tmp1 - A*alpha/tmp2 - p1d
        df = R*T1D/tmp1*(rho1D**-2) + A*a/tmp2*(-2)*(rho1D**-3-rho1D**-2*B)
        rho1D = rho1D - f/df
        if (maxval(abs(f/df)).lt.1e-6) then
            exit
        endif
    enddo

    ! Only set temperature and density
    ! Pressure and e will be set in cgrid initialization
    do i = 1,ny
        rho(:,j,:) = rho1D(j)
        T(:,j,:) = T1D(j)
    enddo

    ! Perturbations: this must be specific for each problem.
    if (noiseAmp.gt.0d0) then
        call perturb_potential(x,y,z,Lx,Lz,u,v,w)
    endif

    deallocate(tmp1,tmp2,p1D,T1D,rho1D,alpha,f,df)
    

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
    use ShearLayer_sCO2_data

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

    use ShearLayer_sCO2_data

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

    use ShearLayer_sCO2_data

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
        call message(2,"Maximum v-velocity",P_MAXVAL(v))
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

    use ShearLayer_sCO2_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
end subroutine
