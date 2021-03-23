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
    subroutine lstab_pert(gp,x,z,fname_prefix,q,qi)
        type(decomp_info), intent(in)   :: gp
        character(len=*), intent(in)    :: fname_prefix 
        real(rkind), dimension(:,:,:), intent(in) :: x, z
        real(rkind), dimension(:,:,:), intent(inout) :: q!the primitive var
        integer, intent(in)       :: qi ! the primitive var col idx in input
    
        complex(rkind), dimension(:,:), allocatable :: qhat
        real(rkind), dimension(:,:), allocatable :: &
            qhat_real, qhat_imag, data2read
        real(rkind), dimension(:), allocatable :: kx, kz
        real(rkind) :: ph 
        complex(rkind) :: e
        character(len=clen) :: fname 
        integer :: nmodes=1, ny, modeID, i, j, k, nx, nz
       
        ! Read mode info in x and y
        fname = trim(fname_prefix)//"_mode_info.dat"
        call read_2d_ascii(data2read,fname)
        nmodes = size(data2read,1)
        allocate(kx(nmodes), kz(nmodes))
        kx = data2read(:,1)!alpha
        kz = data2read(:,2)!beta
        deallocate(data2read)
       
        ! Local sizes of the chunk of domain on this proc:
        ! Assuming y-decomp
        ny = gp%ysz(2)
        nx = size(x,1)
        nz = size(x,3)
        
        ! Allocate based on global domain height
        allocate(qhat_real(ny,nmodes),qhat_imag(ny,nmodes),qhat(ny,nmodes))
    
        ! Read imaginary mode datai
        fname = trim(fname_prefix)//"_init_info_imag.dat"
        call read_2d_ascii(data2read,fname)
        qhat_imag = reshape(data2read(:,qi),[ny,nmodes])
        deallocate(data2read)
    
        ! Read real mode data 
        fname = trim(fname_prefix)//"_init_info_real.dat"
        call read_2d_ascii(data2read,fname)
        qhat_real = reshape(data2read(:,qi),[ny,nmodes])
        deallocate(data2read)
        
        ! Init perturbation fields
        qhat = qhat_real + imi*qhat_imag

        do modeID = 1,nmodes
            do j = 1,ny
                do k = 1,nz
                    do i = 1,nx
                        e = exp(imi*(kx(modeID)*x(i,1,1) + kz(modeID)*z(1,1,k)))
                        q(i,j,k) = q(i,j,k) + real(qhat(j,modeID)*e,rkind) 
                    end do !x 
                end do !z 
            end do !y
        end do !modes
    
        deallocate(kx, kz)
        deallocate(qhat_real,qhat_imag,qhat)
    end subroutine 
    
    subroutine make_pert(gp,x,y,z,Lx,Lz,q,maskWidth,maxAmp,nmodes,Mc)
        type(decomp_info), intent(in)               :: gp
        real(rkind), dimension(:,:,:), intent(in)   :: x,y,z ! a 3D grid
        real(rkind), intent(in)                     :: Lx,Lz  ! grid length
        real(rkind), intent(in)                     :: maskWidth,maxAmp 
        integer, intent(in)                         :: nmodes 
        real(rkind), intent(in)                     :: Mc ! should we use supersonic modes? 
        real(rkind), dimension(:,:,:), intent(inout) :: q !the primitive var

        real(rkind), dimension(gp%ysz(1),gp%ysz(2),gp%ysz(3)) :: tmp 
        complex(rkind), dimension(gp%ysz(1),gp%ysz(2),gp%ysz(3)) :: e 
        real(rkind) :: kx, kz, k, ph, qmax_local, qmax
        integer(rkind) :: j, m, mx, mz
        integer :: mpi_ierr

        ! Add modes to field q
        do m = 3,3+nmodes 
            kx = two*pi/Lx * m
            kz = two*pi/Lz * m
            k  = (kx**two+kz**two)**half
            call random_number(ph)
            call mpi_bcast(ph,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
            e = exp(imi*(kx*x + kz*z)+ph*2*pi)
            q = q + real(k**(-5._rkind/3._rkind)*e,rkind) 
        enddo

        ! Get maxval and scale
        qmax_local = maxval(q)
        call mpi_allreduce(qmax_local, qmax, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD,mpi_ierr)
        tmp = maxAmp*exp(-abs(y)/maskWidth)
        q = q/qmax * tmp

        ! We also need more oscillatory modes at higher Mc
        if (Mc > 0.8) then
            tmp = 0.D0
            do mx = 3,6
            do mz = 3,6
                kx = two*pi/Lx * mx
                kz = two*pi/Lz * mz
                call random_number(ph)
                call mpi_bcast(ph,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
                e = exp(imi*(kx*x + kz*z)+ph*2*pi)
                tmp = tmp + real(e,rkind) 
            enddo
            enddo

            ! Get maxval and scale. Make oscillatory mask
            qmax_local = maxval(tmp)
            call mpi_allreduce(qmax_local, qmax, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD,mpi_ierr)
            q = q + tmp/qmax * 0.5*maxAmp*exp(-abs(y)/maxval(y))*sin(y*2.D0*pi/maskWidth/2)
        endif

    end subroutine 
    
    subroutine perturb_potential(gp,x,y,z,Lx,Lz,u,v,w)
        type(decomp_info), intent(in)               :: gp
        real(rkind), dimension(:,:,:), intent(in)   :: x,y,z
        real(rkind), dimension(:,:,:), intent(inout):: u,v,w 
        real(rkind), intent(in)                     :: Lx,Lz

        real(rkind) :: kx, kz, phx, phz, A, eps=0.15, sigma=5.D0
        integer(rkind) :: i,j, m, mx, mz, nmodes=10
        integer :: mpi_ierr

        do i = 1, nmodes 
        do j = 1, nmodes

            kx = two*pi/Lx * i 
            kz = two*pi/Lx * j
            call random_number(phx)
            call random_number(phz)
            call mpi_bcast(phx,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
            call mpi_bcast(phz,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)

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
    
    subroutine perturb_potential_v2(gp,x,y,z,Lx,Lz,u,v,w)
        type(decomp_info), intent(in)               :: gp
        real(rkind), dimension(:,:,:), intent(in)   :: x,y,z
        real(rkind), dimension(:,:,:), intent(inout):: u,v,w
        real(rkind), intent(in)                     :: Lx,Lz

        real(rkind), dimension(gp%ysz(1),gp%ysz(2),gp%ysz(3)) :: &
            siny,cosy,tmpx,tmpz,tmpy,tmpy2 
        real(rkind) :: kx, kz, phx, phz, A, eps=0.15, sigma=5.D0
        integer(rkind) :: i,j, m, mx, mz, nmodes=50
        integer :: mpi_ierr
        character(len=8) :: c1,c2

        tmpy = exp(-sigma*y*y)
        siny = sin(y)
        cosy = cos(y)
        tmpy2 = -two*sigma*y*(siny+two*sigma*y*cosy ) + (cosy+two*sigma*(cosy-y*siny))
        tmpy2 = tmpy2*tmpy

        call message(0,"Making perturbations")
        do i = 1, nmodes 
        do j = 1, nmodes
            call message(0,"")

            kx = two*pi/Lx * (5+i) 
            kz = two*pi/Lz * (5+j)
            call random_number(phx)
            call random_number(phz)
            call mpi_bcast(phx,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
            call mpi_bcast(phz,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)

            A = eps/(four*pi**two*kx*kz)
            tmpx = kx*x+phx*two*pi
            tmpz = kz*z+phz*two*pi
            u = u + A*kx*sin(tmpx)*cos(tmpz)*tmpy*( siny + two*sigma*y*cosy )
            w = w + A*cos(tmpx)*kz*sin(tmpz)*tmpy*( siny + two*sigma*y*cosy )
            v = v + A*cos(tmpx)*cos(tmpz)*tmpy2
            call message(2,"Maximum tmpx", P_MAXVAL(abs(tmpz)))
            call message(2,"Maximum tmpz", P_MAXVAL(abs(tmpx)))
            call message(2,"Maximum u", P_MAXVAL(abs(u)))
            call message(2,"Maximum v", P_MAXVAL(abs(v)))
            call message(2,"Maximum w", P_MAXVAL(abs(w)))
        enddo
        enddo

        call message(0,"Done making perturbations")
    end subroutine 
end module


subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info, nrank
    use ShearLayer_data

    implicit none
    character(clen)                                :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    namelist /PROBINPUT/ Lx, Ly, Lz,Mc, Re, Pr, Sc,&
                        T_ref, p_ref, rho_ref, rho_ratio,&
                        noiseAmp, fname_prefix, use_lstab, no_pert
    ioUnit = 11
    open(unit=ioUnit, file='input.dat', form='FORMATTED')
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
    real(rkind), dimension(:,:), allocatable :: tmp2D
    real(rkind) :: mu_ref, c1,c2,du, Rgas1, Rgas2,lambda,maskWidth,maxAmp
    integer :: i,j, iounit, nx, ny, nz, nModes
        
    real(rkind) :: kx, kz, ph 
    integer(rkind) :: m, mpi_ierr
    
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
		
        ! Perturbations: this must be specific for each problem.
        if (use_lstab) then
            call lstab_pert(decomp,x,z,fname_prefix,u,2)
            call lstab_pert(decomp,x,z,fname_prefix,v,3)
            call lstab_pert(decomp,x,z,fname_prefix,w,4)
            call lstab_pert(decomp,x,z,fname_prefix,rho,1)
        else if (no_pert) then
            call message(0,"No perturbations")
        else
            call perturb_potential_v2(decomp,x,y,z,Lx,Lz,u,v,w)
            !maskWidth = two*dtheta0
            !maxAmp = 1.D-3*du
            !nModes = 10 
            !call make_pert(decomp,x,y,z,Lx,Lz,u,maskWidth,maxAmp,nModes,Mc)
            !call make_pert(decomp,x,y,z,Lx,Lz,v,maskWidth,maxAmp,nModes,Mc)
            !call make_pert(decomp,x,y,z,Lx,Lz,w,maskWidth,maxAmp,nModes,Mc)
            !call make_pert(decomp,x,y,z,Lx,Lz,p,maskWidth,maxAmp,nModes,Mc)
            !call make_pert(decomp,x,y,z,Lx,Lz,rho,maskWidth,maxAmp,nModes,Mc)
        endif

        ! Add base flow profiles
        lambda = (rho_ratio-1)/(rho_ratio+1)
        u = u + du*half*tanh(y/(two*dtheta0))
        v = v + zero
        w = w + zero
        p = p + p_ref
        rho = rho + rho_ref*(1 + lambda*tanh(y/(two*dtheta0)))
        T = p/(rho*mix%Rgas) 

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

    endif
end subroutine


subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs,tkeb,tsim0,dtheta_0)
    use kind_parameters,    only: rkind
    use decomp_2d,          only: decomp_info
    use MixtureEOSMod,      only: mixture
    use AveragingMod,       only: averaging
    use TKEBudgetMod,       only: tkeBudget
    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    
    !type(averaging), optional,       intent(in)    :: avg
    type(tkeBudget), optional,       intent(inout) :: tkeb
    real(rkind), optional,           intent(inout) :: tsim0 ! the previous time 
    real(rkind), optional,           intent(inout) :: dtheta_0 ! the previous L99 
    
    !real(rkind), dimension(1,size(fields,2),1) :: utilde_p ! on each proc
    real(rkind), dimension(1,decomp%ysz(2),1) :: rbar,utilde,buf ! full size, assuming ystencil  
    real(rkind) :: factor=0.d0, U1,U2,rho0, dtheta=0.d0
    integer :: j,dy,ny,iavg=5
   
    ny = decomp%ysz(2)
    dy = abs(mesh(1,1,1,2)-mesh(1,2,1,2))
    print *, 'dy:',dy 

    ! Compute the momentum thickness:
    ! dtheta = 1/(r0*dU^2) \int{ rbar(U1-utilde)(U2-utilde) }dy
    call tkeb%reynolds_avg(fields(:,:,:,1),rbar)
    call tkeb%favre_avg(fields(:,:,:,1),fields(:,:,:,2)/fields(:,:,:,1),utilde)
    print *, 'Shape of utilde and buf:', size(utilde), size(buf)
    rho0    = rbar  (1,1 ,1)
    U1      = Utilde(1,1 ,1)
    U2      = Utilde(1,ny,1)
    print *, 'rho0:',rho0 
    print *, 'U1:',U2 
    print *, 'U2:',U2 
    buf     = rbar*(U1-utilde)*(U2-utilde) 
    dtheta  = sum(buf)*dy 
    factor  = (dtheta-dtheta0)/(tsim-tsim0) / dtheta
    print *, 'dtheta0:',dtheta_0
    print *, 'dtheta:',dtheta 
    print *, 'factor:',factor
    pause

    rhs(:,:,:,1) = rhs(:,:,:,1) - factor*fields(:,:,:,1)
    rhs(:,:,:,2) = rhs(:,:,:,2) - factor*fields(:,:,:,2)
    rhs(:,:,:,3) = rhs(:,:,:,3) - factor*fields(:,:,:,3)
    rhs(:,:,:,4) = rhs(:,:,:,4) - factor*fields(:,:,:,4)
    rhs(:,:,:,5) = rhs(:,:,:,5) - factor*fields(:,:,:,5)
end subroutine
