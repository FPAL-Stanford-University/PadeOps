module ShearLayer_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, one, two, pi, imi
    use FiltersMod,       only: filters
    use decomp_2d,        only: decomp_info, nrank
    use basic_io,         only: read_2d_ascii 
    use exits,            only: message
    use mpi
    implicit none
    integer, parameter :: ns = 2

    ! Problem parameters
    real(rkind) :: Mc = 0.7_rkind	          ! Convective Mach
    real(rkind) :: Re = 800._rkind	          ! Reynolds number
    real(rkind) :: Sc = 1._rkind	          ! Schmidt number
    real(rkind) :: p_ref = one                ! reference press
    real(rkind) :: T_ref = one                ! reference temp
    real(rkind) :: rho_ref = one              ! reference density
    real(rkind) :: rho_ratio = one            ! rho2/rho1
    real(rkind) :: dtheta0 = 1._rkind          ! Base profile thickness 
    real(rkind) :: noiseAmp = 1D-6              ! white noise amplitude
    character(len=clen) :: InitFileTag, InitFileDirectory
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
    subroutine get_pert(gp,x,z,InitFileTag,InitFileDirectory,q,qi)
        type(decomp_info), intent(in)   :: gp
        character(len=*), intent(in)    :: InitFileTag, InitFileDirectory
        real(rkind), dimension(:,:,:), intent(in) :: x, z
        real(rkind), dimension(:,:,:), intent(inout) :: q!the primitive var
        integer, intent(in)       :: qi ! the primitive var col idx in input
    
        real(rkind), dimension(:,:), allocatable :: &
            qhat_real, qhat_imag, qhat, data2read
        real(rkind), dimension(:), allocatable :: kx, kz
        real(rkind) :: arg1
        complex(rkind) :: e
        character(len=clen) :: fname 
        integer :: nmodes=1, ny, modeID, i, j, k, nx, nz
       
        ! Read mode info in x and y
        fname = trim(InitFileDirectory)//"/"//trim(InitFileTag)//"_mode_info.dat"
        call read_2d_ascii(data2read,fname)
        !nmodes = size(data2read,1)
        allocate(kx(nmodes), kz(nmodes))
        kx = data2read(:,1)!alpha
        kz = data2read(:,2)!beta
        deallocate(data2read)
        call message(0, "Number of normal modes being used:",nmodes)
       
        ! Local sizes of the chunk of domain on this proc:
        ! Assuming y-decomp
        ny = gp%ysz(2)
        nx = size(x,1)
        nz = size(x,3)
        
        ! Allocate based on global domain height
        allocate(qhat_real(ny,nmodes),qhat_imag(ny,nmodes),qhat(ny,nmodes))
    
        ! Read imaginary mode data 
        fname = trim(InitFileDirectory)//"/"//&
                trim(InitFileTag)//"_init_info_imag.dat"
        call read_2d_ascii(data2read,fname)
        qhat_imag = reshape(data2read(:,qi),[ny,nmodes])
        deallocate(data2read)
    
        ! Read real mode data 
        fname = trim(InitFileDirectory)//"/"//&
                trim(InitFileTag)//"_init_info_real.dat"
        call read_2d_ascii(data2read,fname)
        qhat_real = reshape(data2read(:,qi),[ny,nmodes])
        deallocate(data2read)
        
        ! Init perturbation fields
        qhat = qhat_real + imi*qhat_imag
        q = 0.d0

        do modeID = 1,nmodes
            do j = 1,ny
                do k = 1,nz
                    do i = 1,nx
                        e = exp(imi*(kx(modeID)*x(i,1,1) + kz(modeID)*z(1,1,k) ))
                        q(i,j,k) = q(i,j,k) + real(qhat(j,modeID)*e,rkind) 
                    end do !x 
                end do !z 
            end do !y
        end do !modes
    
        deallocate(kx, kz)
        deallocate(qhat_real,qhat_imag,qhat)
    end subroutine 
    
    subroutine make_pert(gp,x,z,Lx,Lz,q)
        type(decomp_info), intent(in)           :: gp
        real(rkind), dimension(:,:), intent(in) :: x,z
        real(rkind), intent(in)                 :: Lx, Lz
        real(rkind), dimension(:,:,:), intent(inout) :: q!the primitive var

        real(rkind) :: kx, kz
        complex(rkind) :: coeff
        complex(rkind), dimension(gp%ysz(1),gp%ysz(3)) :: e
        integer :: ny, i, j, k, nx, nz, gnx, gnz, wx, wz
      
        ! some global properties...
        gNx = gp%xsz(1) 
        gNz = gp%zsz(3)
       
        ! Local sizes of the chunk of domain on this block: Assuming y-decomp
        ny = gp%ysz(2)
        nx = gp%ysz(1)
        nz = gp%ysz(3)
       
        ! Init perturbation fields
        q = 0.d0
        do j = 1,ny
            do wx = -gNx/2,gNx/2-1
            do wz = -gNz/2,gNz/2-1
                kx = wx*2*pi/Lx
                kz = wx*2*pi/Lz
                e  = exp(imi*(kx*x + kz*z))
                coeff = (kx**(-5/3)+kz**(-5/3))!should be complex for phase speed
                q(:,j,:) = q(:,j,:) + real( coeff*e,rkind )
            end do
            end do
        end do
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
                        noiseAmp, InitFileTag, InitFileDirectory
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
    use constants,                   only: zero,half,one,two,five,pi,eight
    use CompressibleGrid,            only: rho_index,u_index,v_index,w_index,&
                                           p_index,T_index,e_index,Ys_index
    use decomp_2d,                   only: decomp_info,nrank
    use MixtureEOSMod,               only: mixture
    use IdealGasEOS,                 only: idealgas
    use PowerLawViscosityMod,		 only: powerLawViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use ConstPrandtlConductivityMod, only: constPrandtlConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use random,                      only: gaussian_random                  
    use mpi
    use ShearLayer_data

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    integer :: i, iounit, nx, ny, nz
    integer :: seedu=321341, seedv=423424, seedw=131344, &
            seedr=452123,seedp=456321, seedT=321644
    type(powerLawViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(constPrandtlConductivity) :: thermcond
    real(rkind), dimension(:,:,:), allocatable :: tmp, pert
    real(rkind) :: mu_ref, c1,c2,du, Rgas1, Rgas2,lambda
    
    namelist /PROBINPUT/ Lx, Ly, Lz,Mc, Re, Pr, Sc,&
                        T_ref, p_ref, rho_ratio, rho_ref, &
                        noiseAmp, InitFileTag, InitFileDirectory
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    nx = decomp%xsz(1) 
    ny = decomp%ysz(2) 
    nz = decomp%zsz(3)
    allocate(tmp(nx,ny,nz))
    allocate(pert(nx,ny,nz))

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

        ! Also add two passive tracers for each species
        do i = mix%ns+1,mix%ns*2
            shearvisc = powerLawViscosity( mu_ref, T_ref, 0._rkind)
            bulkvisc  = constRatioBulkViscosity( zero )
            thermcond = constPrandtlConductivity( Pr )
            call mix%set_material( i, idealgas( gam, Rgas(i-mix%ns) ),&
                 shearvisc = shearvisc, & 
                 bulkvisc  = bulkvisc, &
                 thermcond = thermcond  )
        end do

        ! Set mass diffusivity object (Ensure that all units are consistent)
        lambda = (rho_ratio-1)/(rho_ratio+1)
        call mix%set_massdiffusivity( constSchmidtDiffusivity( mu_ref,rho_ref,Sc))
        tmp = half*(one+lambda*tanh(y/(two*dtheta0)))
        Ys(:,:,:,1)  = one - tmp
        Ys(:,:,:,2)  = one - Ys(:,:,:,1)
        Ys(:,:,:,3)  = zero 
        Ys(:,:,:,4)  = one 
        call mix%update(Ys)
		
        ! Base flow profiles
        u = du*half*tanh(y/(two*dtheta0))
        v = zero
        w = zero
        rho = rho_ref*(1 + lambda*tanh(y/(two*dtheta0)))
        p = p_ref
        T = p/(rho*mix%Rgas) 
        
        ! Modal perturbations: this must be specific for each problem.
        call make_pert(decomp,x(:,1,:),z(:,1,:),Lx,Lz,pert)
        rho = rho+pert
        !call get_pert(decomp, x, z, InitFileTag, InitFileDirectory, pert, 1)
        !rho = rho + pert
        !call get_pert(decomp, x, z, InitFileTag, InitFileDirectory, pert, 2)
        !u = u + pert
        !call get_pert(decomp, x, z, InitFileTag, InitFileDirectory, pert, 3)
        !v = v + pert
        !call get_pert(decomp, x, z, InitFileTag, InitFileDirectory, pert, 4)
        !w = w + pert
        !call get_pert(decomp, x, z, InitFileTag, InitFileDirectory, pert, 5)
        !T = T + pert
        !call get_pert(decomp, x, z, InitFileTag, InitFileDirectory, pert, 6)
        !p = p + pert
        
        ! Gaussian noise decaying exponentially
        tmp = exp(-abs(y/(two*dtheta0)))
        call gaussian_random(pert,zero,one,seedu+100*nrank)
        u = u + noiseAmp*pert*tmp
        call gaussian_random(pert,zero,one,seedv+100*nrank)
        v = v + noiseAmp*pert*tmp
        call gaussian_random(pert,zero,one,seedw+100*nrank)
        w = w + noiseAmp*pert*tmp
        call gaussian_random(pert,zero,one,seedr+100*nrank)
        rho = rho + noiseAmp*pert*tmp
        call gaussian_random(pert,zero,one,seedp+100*nrank)
        p = p + noiseAmp*pert*tmp
        call gaussian_random(pert,zero,one,seedT+100*nrank)
        T = T + noiseAmp*pert*tmp

        ! Initialize gaussian filter mygfil
        call mygfil%init(decomp, periodicx, periodicy, periodicz, &
                        "gaussian", "gaussian", "gaussian" )
        
        deallocate(tmp,pert)
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


subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs,rhsg)
    use kind_parameters, only: rkind
    use decomp_2d,       only: decomp_info
    use MixtureEOSMod,    only: mixture

    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    real(rkind), dimension(:,:,:,:), optional, intent(inout) ::rhsg
end subroutine
