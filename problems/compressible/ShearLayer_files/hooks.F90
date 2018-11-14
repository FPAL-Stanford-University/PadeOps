module ShearLayer_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, one
    use FiltersMod,       only: filters
    use decomp_2d,        only: decomp_info
    use basic_io,         only: read_2d_ascii 
    use exits,            only: message
    implicit none
    integer, parameter :: ns = 2

    ! Problem parameters
    real(rkind) :: Mach = 0.7_rkind	          ! Convective Mach
    real(rkind) :: Re = 800._rkind	          ! Reynolds number
    real(rkind) :: Sc = 1._rkind	          ! Schmidt number
    real(rkind) :: rho_ratio = one            ! rho2/rho1
    real(rkind) :: dtheta = 1._rkind          ! Base profile thickness 

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

    ! Perturbation parameters
    character(len=clen) :: pertfile='pertfile.dat'
    real(rkind) :: kx, amp_x
    real(rkind) :: kz, amp_z
    real(rkind), dimension(:), allocatable :: perty
    real(rkind), dimension(:,:,:), allocatable :: perty_field
contains

    subroutine read_perturbation_files(nx,ny,nz)
        use mpi 
        use decomp_2d, only: nrank
        use constants, only: zero, two, pi
        integer, intent(in) :: nx,ny,nz
        integer :: pertunit
        integer :: ierr
        integer :: i, L_by_lambda_z

        pertunit = 229
		allocate(perty(ny))
		allocate(perty_field(nx,ny,nz))
        
        ! Read in the y perturbation file
        if (nrank == 0) then
            print *,  "Reading in the perturbation file ", trim(pertfile)
            open(unit=pertunit,file=trim(pertfile),form='FORMATTED',status='OLD')
            read(pertunit,*) Lx, Ly, Lz
        end if

        ! broadcast perturbations to all procs
        call mpi_bcast(kx, 1, mpirkind, 0, mpi_comm_world, ierr)
        call mpi_bcast(kz, 1, mpirkind, 0, mpi_comm_world, ierr)
        call mpi_bcast(amp_x, 1, mpirkind, 0, mpi_comm_world, ierr)
        call mpi_bcast(amp_z, 1, mpirkind, 0, mpi_comm_world, ierr)

        if (nrank == 0) then
            do i=1,ny
                read(pertunit,*) perty(i)
                perty_field(:,i,:)=perty(i)
            end do
            close(pertunit)
        end if

        call mpi_bcast(perty, ny, mpirkind, 0, mpi_comm_world, ierr)
    end subroutine

    subroutine read_domain_info(fname,Lx,Ly,Lz)
        use mpi 
        use decomp_2d, only: nrank
        character(len=*),   intent(in) :: fname 
        real(rkind),        intent(out) :: Lx, Ly, Lz
        integer :: funit=229
        
        ! Read in the y perturbation file
        if (nrank == 0) then
            print *,  "Reading domain info", trim(fname)
            open(unit=funit,file=trim(fname),form='FORMATTED',status='OLD')
            read(funit,*) Lx
            read(funit,*) Ly
            read(funit,*) Lz
            print *, "Domain size: ",Lx,Ly,Lz
        end if
    end subroutine


    subroutine get_perturbations(gp,x,z,InitFileTag,InitFileDirectory,u,v,w)
        use constants, only: imi, pi    
        type(decomp_info),  intent(in) :: gp
        character(len=*),   intent(in) :: InitFileTag, InitFileDirectory
        real(rkind), dimension(:,:,:),    intent(in) :: x, z
        real(rkind), dimension(:,:,:),    intent(out) :: u, v, w
    
        real(rkind), dimension(:,:), allocatable :: & 
                      data2read, uhat_real, uhat_imag, vhat_real, vhat_imag, &
                      what_real, what_imag 
        complex(rkind), dimension(:,:), allocatable :: uhat, vhat, what
        real(rkind), dimension(:), allocatable :: kx, kz, kmode, beta
        real(rkind) :: arg1
        complex(rkind) :: expfact, uhatn, vhatn, whatn
        character(len=clen) :: fname 
        integer :: nmodes, ny, modeID, i, j, k
       
        ! Read mode info in x and y
        fname = trim(InitFileDirectory)//"/"//trim(InitFileTag)//"_mode_info.dat"
        call read_2d_ascii(data2read,fname)
        nmodes = size(data2read,1)
        allocate(kx(nmodes), kz(nmodes))
        allocate(beta(nmodes), kmode(nmodes))
        beta  = data2read(:,1)
        kmode = data2read(:,2)
        kx = -kmode*sin(beta*pi/180.d0)
        kz =  kmode*cos(beta*pi/180.d0)
        deallocate(data2read)
        call message(0, "Number of normal modes being used:",nmodes)
        
        ! Allocate based on global domain height
        ny = gp%ysz(2) 
        allocate(uhat_real(ny,nmodes),vhat_real(ny,nmodes),what_real(ny,nmodes))
        allocate(uhat_imag(ny,nmodes),vhat_imag(ny,nmodes),what_imag(ny,nmodes))
        allocate(uhat(ny,nmodes),vhat(ny,nmodes),what(ny,nmodes))
    
        ! Read imaginary mode data 
        fname = trim(InitFileDirectory)//"/"//trim(InitFileTag)//"_init_info_imag.dat"
        call read_2d_ascii(data2read,fname)
        uhat_imag = reshape(data2read(:,1),[ny,nmodes])
        vhat_imag = reshape(data2read(:,2),[ny,nmodes])
        what_imag = reshape(data2read(:,3),[ny,nmodes])
        deallocate(data2read)
    
        ! Read real mode data 
        fname = trim(InitFileDirectory)//"/"//trim(InitFileTag)//"_init_info_real.dat"
        call read_2d_ascii(data2read,fname)
        uhat_real = reshape(data2read(:,1),[ny,nmodes])
        vhat_real = reshape(data2read(:,2),[ny,nmodes])
        what_real = reshape(data2read(:,3),[ny,nmodes])
        deallocate(data2read)
       
        ! Init perturbation fields
        uhat = uhat_real + imi*uhat_imag
        vhat = vhat_real + imi*vhat_imag
        what = what_real + imi*what_imag
        u = 0.d0
        v = 0.d0
        w = 0.d0
    
        do modeID = 1,nmodes
            arg1  = 0!beta(modeID)*pi/180.d0
            do j = gp%yst(2),gp%yen(2)
                do k = 1,gp%ysz(3)
                    do i = 1,gp%ysz(1)
                        expfact = exp(imi*kz(modeID)*z(1,1,k) + imi*kx(modeID)*x(i,1,1) )
                        uhatn = uhat(j,modeID)*cos(arg1) - vhat(j,modeID)*sin(arg1)
                        vhatn = uhat(j,modeID)*sin(arg1) + vhat(j,modeID)*cos(arg1)
                        whatn = what(j,modeID) 
                        u(i,j-gp%yst(2)+1,k) = u(i,j-gp%yst(2)+1,k) + real(uhatn*expfact,rkind) 
                        v(i,j-gp%yst(2)+1,k) = v(i,j-gp%yst(2)+1,k) + real(vhatn*expfact,rkind) 
                        w(i,j-gp%yst(2)+1,k) = w(i,j-gp%yst(2)+1,k) + real(whatn*expfact,rkind)
                    end do !x 
                end do !z 
            end do !y
        end do !modes
    
        deallocate(uhat_real,vhat_real,what_real,uhat_imag,vhat_imag,what_imag)
        deallocate(beta,kmode)
        deallocate(uhat, vhat, what)
        deallocate(kx, kz)
    end subroutine 

end module


subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info
    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    !character(len=*),                intent(in)    :: inputfile
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    !real(rkind) :: p_ref=one, mu_ref=one, T_ref=one, rho_ref=one, &
    !               c1,c2,du, noiseAmp
    !character(len=clen) :: InitFileTag, InitFileDirectory
    !
    !namelist /PROBINPUT/ Mach, Re, Pr, Sc, rho_ref, p_ref, T_ref, rho_ratio,&
    !                    gam, dtheta, InitFileTag, InitFileDirectory, noiseAmp,&
    !                    Lx, Ly, Lz
    !ioUnit = 11
    !open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    !read(unit=ioUnit, NML=PROBINPUT)
    !close(ioUnit)

    ! domain size for this proc
    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Need to set x, y and z as well as  dx, dy and dz
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        call read_domain_info("/home/kmatsuno/ShearLayerInput/PadeInput/ShearLayer_domain_info.dat",Lx,Ly,Lz)
        dx = Lx/real(nx,rkind)
        dy = Ly/real(ny-1,rkind)
        dz = Lz/real(nz,rkind)
        x1 = -Lx/2._rkind
        y1 = -Ly/2._rkind
        z1 = -Lz/2._rkind

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
    use kind_parameters ,             only: rkind, clen
    use constants,                   only: zero,half,one,two,pi,eight
    use CompressibleGrid,            only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,Ys_index
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
    integer :: seedu=321341, seedv=423424, seedw=131344
    type(powerLawViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(constPrandtlConductivity) :: thermcond

    real(rkind), dimension(:,:,:), allocatable :: tmp, upert, vpert, wpert 
    real(rkind) :: p_ref=one, mu_ref=one, T_ref=one, rho_ref=one, &
                   c1,c2,du, Rgas1, Rgas2, noiseAmp
    character(len=clen) :: InitFileTag, InitFileDirectory
    
    namelist /PROBINPUT/ Mach, Re, Pr, Sc, rho_ref, p_ref, T_ref, rho_ratio,&
                        gam, dtheta, InitFileTag, InitFileDirectory, noiseAmp,&
                        Lx, Ly, Lz
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    nx = decomp%xsz(1) 
    ny = decomp%ysz(2) 
    nz = decomp%zsz(3)
    allocate(tmp(nx,ny,nz))
    allocate(upert(nx,ny,nz))
    allocate(vpert(nx,ny,nz))
    allocate(wpert(nx,ny,nz))

    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),                    &
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),                    &
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),                    &
                 e => fields(:,:,:,  e_index), Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (mix%ns /= ns) call GracefulExit("Wrong number of species. Check your input file and make ns consistent with the problem file.",4562)

        ! First set the materials for rho_ratio = rho2/rho1
        Rgas1 = (1+rho_ratio)/two			
        Rgas2 = (1+rho_ratio)/(two*rho_ratio) 
		Rgas = [Rgas1, Rgas2]
		mu_ref = 1/Re

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
        call mix%set_massdiffusivity( constSchmidtDiffusivity( mu_ref,rho_ref,Sc))
        tmp = half*(one-tanh(y/(2*dtheta)))
        Ys(:,:,:,1)  = one - tmp
        Ys(:,:,:,2)  = one - Ys(:,:,:,1)
        call mix%update(Ys)
		
        ! Base flow profiles 
        c1 = sqrt(gam*p_ref/(rho_ref/Rgas1))
        c2 = sqrt(gam*p_ref/(rho_ref/Rgas2))
		du = Mach*(c1+c2)
        u   = zero;!du*(tmp-half)
        v   = zero
        w   = zero
        p   = p_ref
        T   = T_ref
        rho = p / (mix%Rgas * T)

        ! Modal perturbations
        call get_perturbations(decomp, x, z, InitFileTag, InitFileDirectory, &
                 upert, vpert, wpert)
        u = u + upert
        v = v + vpert
        w = w + wpert

        ! Gaussian noise
        !call gaussian_random(upert,zero,one,seedu+100*nrank)
        !call gaussian_random(vpert,zero,one,seedu+100*nrank)
        !call gaussian_random(wpert,zero,one,seedu+100*nrank)
        !u = u + noiseAmp*upert
        !v = v + noiseAmp*vpert
        !w = w + noiseAmp*wpert
        deallocate(upert, vpert, wpert)

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

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ShearLayer_", vizcount, ".dat"

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
        filpt = 10.0_rkind/dy 
        thickT = real(10.D0, rkind)
        do i=1,decomp%ysz(2)
            dumT(:,i,:)=half*(one-tanh( (real( decomp%yst(2) - 1 + i - 1, rkind)-filpt) / thickT ))
        end do
            
            
        ! Gussian Filter for last N points in x-direction to act as a sponge (A)
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

        do i=1,decomp%ysz(2)
            dumT(:,i,:)=half*(one-tanh( (real(decomp%ysz(2)- (decomp%yst(2) - 1 + i - 1), rkind)-filpt) / thickT ))
        end do
            
            
        ! Gussian Filter for last N points in x-direction to act as a sponge (A)
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
        !call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        !call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        !call message(2,"Maximum conductivity",P_MAXVAL(kap))
        !call message(2,"Maximum diffusivity",P_MAXVAL(diff))

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
