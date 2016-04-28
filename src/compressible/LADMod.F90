module LADMod

    use kind_parameters, only: rkind
    use constants,       only: zero,half,one
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters

    implicit none

    type :: ladobject

        real(rkind) :: Cbeta
        real(rkind) :: Cmu
        real(rkind) :: Ckap
        real(rkind) :: Cdiff
        real(rkind) :: CY

        type(decomp_info), pointer :: decomp
        type(derivatives), pointer :: der
        type(filters),     pointer :: fil

        integer :: nfils
        real(rkind) :: dx, dy, dz

    contains

        procedure          :: init
        procedure          :: get_viscosities
        procedure          :: get_conductivity
        procedure          :: get_diffusivity
        procedure, private :: filter
        procedure          :: destroy

    end type

contains

    subroutine init(this,decomp,der,fil,nfils,dx,dy,dz,Cbeta,Cmu,Ckap,Cdiff,CY)
        class(ladobject),        intent(in) :: this
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        type(filters),     intent(in) :: fil
        integer,           intent(in) :: nfils
        real(rkind),       intent(in) :: Cbeta,Cmu,Ckap,Cdiff,CY

        ! Set all coefficients
        this%Cbeta = Cbeta
        this%Cmu = Cmu
        this%Ckap = Ckap
        this%Cdiff = Cdiff
        this%CY = CY

        ! Point type pointers to external types
        this%decomp => decomp
        this%der => der
        this%fil => fil

        ! Set number of times to filter
        this%nfils = nfils

        ! Set grid spacing
        this%dx = dx
        this%dy = dy
        this%dz = dz
    end subroutine

    subroutine get_viscosities(this,rho,duidxj,mu,bulk,x_bc,y_bc,z_bc)
        class(ladobject),        intent(in) :: this
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)),           intent(in)  :: rho
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),9), target, intent(in)  :: duidxj
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)),           intent(inout) :: mu, bulk
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), pointer, intent(in) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: mustar, bulkstar
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xtmp1,xtmp2
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5,func
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: ztmp1,ztmp2

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
       
        ! -------- Artificial Shear Viscosity --------

        ! Magnitude of strain rate
        func = sqrt(dudx**2 + half*(dvdx+dudy)**2 + half*(dwdx+dudz)**2 &
                            +             dvdy**2 + half*(dwdy+dvdz)**2 &
                                                  +             dwdz**2 )
        
        ! Get 4th derivative in X
        call transpose_y_to_x(func,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
        xtmp2 = xtmp1*this%dx**6
        call transpose_x_to_y(xtmp2,mustar,this%decomp)
        
        ! Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
        ztmp2 = ztmp1*this%dz**6
        call transpose_z_to_y(ztmp2,ytmp1,this%decomp)
        mustar = mustar + ytmp1
        
        ! Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp1,y_bc(1),y_bc(2))
        call this%der%d2dy2(ytmp1,ytmp2,y_bc(1),y_bc(2))
        ytmp1 = ytmp2*this%dy**6
        mustar = mustar + ytmp1

        mustar = this%Cmu*rho*abs(mustar)
        
        ! Filter mustar
        call this%filter(mustar, x_bc, y_bc, z_bc)
        
        mu = mu + mustar

        ! -------- Artificial Bulk Viscosity --------
        
        func = dudx + dvdy + dwdz      ! dilatation
        
        ! Step 1: Get components of grad(rho) squared individually
        call this%gradient(rho,ytmp1,ytmp2,ytmp3,x_bc,y_bc,z_bc) ! Does not use any Y buffers
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(func,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
        xtmp2 = xtmp1*this%dx**4
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        bulkstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
        ztmp2 = ztmp1*this%dz**4
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        bulkstar = bulkstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp4,y_bc(1),y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
        ytmp4 = ytmp5*this%dy**4
        bulkstar = bulkstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2

        ! Now, all ytmps are free to use
        ytmp1 = dwdy-dvdz; ytmp2 = dudz-dwdx; ytmp3 = dvdx-dudy
        ytmp4 = ytmp1*ytmp1 + ytmp2*ytmp2 + ytmp3*ytmp3 ! |curl(u)|^2
        ytmp2 = func*func ! dilatation^2

        ! Calculate the switching function
        ytmp1 = ytmp2 / (ytmp2 + ytmp4 + real(1.0D-32,rkind)) ! Switching function f_sw
        where (func .GE. zero)
            ytmp1 = zero
        end where

        bulkstar = this%Cbeta*rho*ytmp1*abs(bulkstar)

        ! Filter bulkstar
        call this%filter(bulkstar, x_bc, y_bc, z_bc)

        bulk = bulk + bulkstar
    end subroutine

    subroutine get_conductivity(this,rho,e,T,sos,kap,x_bc,y_bc,z_bc)
        class(ladobject),  intent(in) :: this
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(in)  :: rho,e,T,sos
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(inout) :: kap
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: kapstar
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xtmp1,xtmp2
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: ztmp1,ztmp2

        ! -------- Artificial Conductivity --------

        ! Step 1: Get components of grad(e) squared individually
        call this%gradient(e,ytmp1,ytmp2,ytmp3,x_bc,y_bc,z_bc) ! Does not use any Y buffers
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(e,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
        xtmp2 = xtmp1*this%dx**4
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        kapstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(e,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
        ztmp2 = ztmp1*this%dz**4
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        kapstar = kapstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(e,ytmp4,y_bc(1),y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
        ytmp4 = ytmp5*this%dy**4
        kapstar = kapstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        kapstar = this%Ckap*rho*sos*abs(kapstar)/T

        ! Filter kapstar
        call this%filter(kapstar, x_bc, y_bc, z_bc)

        kap = kap + kapstar
    end subroutine

    subroutine get_diffusivity(this,Ys,dYsdx,dYsdy,dYsdz,sos,diff)
        class(ladobject),  intent(in) :: this
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(in)  :: rho,Ys,dYsdx,dYsdy,dYsdz,sos
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(inout) :: diff
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: diffstar
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xtmp1,xtmp2
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: ztmp1,ztmp2

        ! -------- Artificial Diffusivity ---------

        ! Step 1: Get components of grad(Ys) squared individually
        ytmp1 = dYsdx*dYsdx
        ytmp2 = dYsdy*dYsdy
        ytmp3 = dYsdz*dYsdz

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(Ys,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
        xtmp2 = xtmp1*this%dx**4
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        diffstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(Ys,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
        ztmp2 = ztmp1*this%dz**4
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        diffstar = diffstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(Ys,ytmp4,y_bc(1),y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
        ytmp4 = ytmp5*this%dy**4
        diffstar = diffstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        diffstar = this%Cdiff*sos*abs(diffstar) ! CD part of diff

        ytmp4 = (sqrt(ytmp1)*this%dx + sqrt(ytmp2)*this%dy + sqrt(ytmp3)*this%dz) &
                        / ( sqrt(ytmp1+ytmp2+ytmp3) + real(1.0D-32,rkind) ) ! grid scale
        ytmp5 = this%CY*sos*( half*(abs(Ys)-one + abs(Ys-one)) )*ytmp4 ! CY part of diff

        ! Filter each part
        call this%filter(diffstar, x_bc, y_bc, z_bc)
        call this%filter(ytmp5, x_bc, y_bc, z_bc)

        diffstar = max(diffstar, ytmp5) ! Take max of both terms instead of add to minimize the dissipation

        diff = diff + diffstar
    end subroutine

    subroutine filter(this,arr,x_bc,y_bc,z_bc)
        class(ladobject), intent(in) :: this
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(inout) :: arr
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: tmp_in_y
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: tmp1_in_x, tmp2_in_x
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: tmp1_in_z, tmp2_in_z
        
        ! First filter in y
        call this%fil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
        ! Subsequent refilters 
        do idx = 1,this%nfils-1
            arr = tmp_in_y
            call this%fil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
        end do
        
        ! Then transpose to x
        call transpose_y_to_x(tmp_in_y,tmp1_in_x,this%decomp)

        ! First filter in x
        call this%fil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
        ! Subsequent refilters
        do idx = 1,this%nfils-1
            tmp1_in_x = tmp2_in_x
            call this%fil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
        end do 

        ! Now transpose back to y
        call transpose_x_to_y(tmp2_in_x,tmp_in_y,this%decomp)

        ! Now transpose to z
        call transpose_y_to_z(tmp_in_y,tmp1_in_z,this%decomp)

        !First filter in z
        call this%fil%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
        ! Subsequent refilters
        do idx = 1,this%nfils-1
            tmp1_in_z = tmp2_in_z
            call this%fil%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
        end do 

        ! Now transpose back to y
        call transpose_z_to_y(tmp2_in_z,arr,this%decomp)

        ! Finished
    end subroutine
   
end module
