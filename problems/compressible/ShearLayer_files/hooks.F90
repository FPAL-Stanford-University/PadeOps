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
    real(rkind) :: rho_ref = one              ! rho1 = 1/R?
    real(rkind) :: rho_ratio = one            ! rho2/rho1
    real(rkind) :: vel_ratio = one            ! U2/U1
    real(rkind) :: dtheta0 = 1._rkind          ! Base profile thickness 
    real(rkind) :: noiseAmp = 1D-3              ! white noise amplitude
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

    subroutine read_domain_info(fname,Lx,Ly,Lz)
        character(len=*),   intent(in) :: fname 
        real(rkind),        intent(out) :: Lx, Ly, Lz
        integer :: funit=229
        
        ! Read in the y perturbation file
        if (nrank == 0) then
            open(unit=funit,file=trim(fname),form='FORMATTED',status='OLD')
            read(funit,*) Lx
            read(funit,*) Ly
            read(funit,*) Lz
            print *, "Domain size: ",Lx,Ly,Lz
        end if
    end subroutine


    subroutine get_perturbations(gp,x,z,InitFileTag,InitFileDirectory,u,v,w,T,p)
        type(decomp_info),  intent(in) :: gp
        character(len=*),   intent(in) :: InitFileTag, InitFileDirectory
        real(rkind), dimension(:,:,:),    intent(in)    :: x, z
        real(rkind), dimension(:,:,:),    intent(inout) :: u, v, w, T, p
    
        real(rkind), dimension(:,:), allocatable :: & 
                      data2read, uhat_real, uhat_imag, vhat_real, vhat_imag, &
                      what_real, what_imag, That_real, That_imag, phat_real,phat_imag 
        complex(rkind), dimension(:,:), allocatable :: uhat, vhat, what, That,phat
        real(rkind), dimension(:), allocatable :: kx, kz
        real(rkind) :: arg1
        complex(rkind) :: e
        character(len=clen) :: fname 
        integer :: nmodes, ny, modeID, i, j, k, nx, nz
       
        ! Read mode info in x and y
        fname = trim(InitFileDirectory)//"/"//trim(InitFileTag)//"_mode_info.dat"
        call read_2d_ascii(data2read,fname)
        nmodes = size(data2read,1)
        allocate(kx(nmodes), kz(nmodes))
        kx = data2read(:,1)!alpha
        kz = data2read(:,2)!beta
        deallocate(data2read)
        call message(0, "Number of normal modes being used:",nmodes)
       
        ! Local sizes of the chunk of domain on this proc:
        ny = gp%ysz(2)
        nx = size(x,1)
        nz = size(x,3)

        ! Allocate based on global domain height
        allocate(uhat_real(ny,nmodes),uhat_imag(ny,nmodes),uhat(ny,nmodes))
        allocate(vhat_real(ny,nmodes),vhat_imag(ny,nmodes),vhat(ny,nmodes))
        allocate(what_real(ny,nmodes),what_imag(ny,nmodes),what(ny,nmodes))
        allocate(That_real(ny,nmodes),That_imag(ny,nmodes),That(ny,nmodes))
        allocate(phat_real(ny,nmodes),phat_imag(ny,nmodes),phat(ny,nmodes))
    
        ! Read imaginary mode data 
        fname = trim(InitFileDirectory)//"/"//&
                trim(InitFileTag)//"_init_info_imag.dat"
        call read_2d_ascii(data2read,fname)
        uhat_imag = reshape(data2read(:,2),[ny,nmodes])
        vhat_imag = reshape(data2read(:,3),[ny,nmodes])
        what_imag = reshape(data2read(:,4),[ny,nmodes])
        That_imag = reshape(data2read(:,5),[ny,nmodes])
        phat_imag = reshape(data2read(:,6),[ny,nmodes])
        deallocate(data2read)
    
        ! Read real mode data 
        fname = trim(InitFileDirectory)//"/"//&
                trim(InitFileTag)//"_init_info_real.dat"
        call read_2d_ascii(data2read,fname)
        uhat_real = reshape(data2read(:,2),[ny,nmodes])
        vhat_real = reshape(data2read(:,3),[ny,nmodes])
        what_real = reshape(data2read(:,4),[ny,nmodes])
        That_real = reshape(data2read(:,5),[ny,nmodes])
        phat_real = reshape(data2read(:,6),[ny,nmodes])
        deallocate(data2read)
       
        ! Init perturbation fields
        uhat = uhat_real + imi*uhat_imag
        vhat = vhat_real + imi*vhat_imag
        what = what_real + imi*what_imag
        That = That_real + imi*That_imag
        phat = phat_real + imi*phat_imag
        u = 0.d0
        v = 0.d0
        w = 0.d0
        T = 0.d0
        p = 0.d0
      
        do modeID = 1,nmodes
            do j = 1,ny
                do k = 1,nz
                    do i = 1,nx
                        e = exp(imi*(kx(modeID)*x(i,1,1) + kz(modeID)*z(1,1,k) ))
                        u(i,j,k) = u(i,j,k) + real(uhat(j,modeID)*e,rkind) 
                        v(i,j,k) = v(i,j,k) + real(vhat(j,modeID)*e,rkind) 
                        w(i,j,k) = w(i,j,k) + real(what(j,modeID)*e,rkind)
                        T(i,j,k) = T(i,j,k) + real(That(j,modeID)*e,rkind) 
                        p(i,j,k) = p(i,j,k) + real(phat(j,modeID)*e,rkind)
                    end do !x 
                end do !z 
            end do !y
        end do !modes
    
        deallocate(kx, kz)
        deallocate(uhat_real,uhat_imag,uhat)
        deallocate(vhat_real,vhat_imag,vhat)
        deallocate(what_real,what_imag,what)
        deallocate(That_real,That_imag,That)
        deallocate(phat_real,phat_imag,phat)
    end subroutine 

end module


subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info, nrank
    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    character(clen)                   :: inputfile
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    inputfile = './input.dat'
    namelist /PROBINPUT/ Lx, Ly, Lz,Mc, Re, Pr, Sc,&
                        rho_ref, rho_ratio, vel_ratio,&
                        noiseAmp, InitFileTag, InitFileDirectory
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
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

    real(rkind), dimension(:,:,:), allocatable :: &
        tmp, upert, vpert, wpert, Tpert, ppert
    real(rkind) :: p_ref, mu_ref, T_ref, c1,c2,du, Rgas1, Rgas2,lambda
    
    namelist /PROBINPUT/ Lx, Ly, Lz,Mc, Re, Pr, Sc,&
                        rho_ref, rho_ratio, vel_ratio,&
                        noiseAmp, InitFileTag, InitFileDirectory
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
    allocate(Tpert(nx,ny,nz))
    allocate(ppert(nx,ny,nz))

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
        
        ! Constant for the problem
        c1 = one
        c2 = c1/sqrt(rho_ratio)
		du = Mc*(c1+c2)
        lambda = (1-rho_ratio)/(1+rho_ratio)
        T_ref = one/gam
        mu_ref = rho_ref*dU*5*dtheta0/Re

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
        tmp = half*(one-tanh(y/(2*dtheta0)))
        Ys(:,:,:,1)  = one - tmp
        Ys(:,:,:,2)  = one - Ys(:,:,:,1)
        call mix%update(Ys)
		
        ! Base flow profiles
        rho = rho_ref*( 1+lambda*tanh(-y/(2*dtheta0)) ) 
        u = dU*half*tanh(-y/(2*dtheta0))!du*(tmp-half)
        v = zero
        w = zero
        p = rho*c1**2/gam
        T = T_ref

        ! Modal perturbations: this must be specific for each problem.
        call get_perturbations(decomp, x, z, InitFileTag, InitFileDirectory, &
                 upert, vpert, wpert, Tpert, ppert)
        u = u + noiseAmp*upert
        v = v + noiseAmp*vpert
        w = w + noiseAmp*wpert
        T = T + noiseAmp*Tpert
        p = p + noiseAmp*ppert
        
        ! Gaussian noise
        call gaussian_random(upert,zero,one,seedu+100*nrank)
        call gaussian_random(vpert,zero,one,seedv+100*nrank)
        call gaussian_random(wpert,zero,one,seedw+100*nrank)
        u = u + noiseAmp*10*upert
        v = v + noiseAmp*10*vpert
        w = w + noiseAmp*10*wpert
        deallocate(upert, vpert, wpert, Tpert, ppert)

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
