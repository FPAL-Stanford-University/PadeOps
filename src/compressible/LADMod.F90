module LADMod

    use kind_parameters, only: rkind
    use constants,       only: zero,half,one
    use decomp_2d,       only: decomp_info, transpose_y_to_x, transpose_x_to_y, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use operators,       only: gradient

    implicit none

    type :: ladobject

        real(rkind) :: Cbeta
        real(rkind) :: CbetaP
        real(rkind) :: Cmu
        real(rkind) :: Ckap
        real(rkind) :: CkapP
        real(rkind) :: Cdiff
        real(rkind) :: CY
        real(rkind) :: Cdiff_g
        real(rkind) :: Cdiff_gt
        real(rkind) :: Cdiff_gp
        real(rkind) :: Cdiff_pe
        real(rkind) :: Cdiff_pe_2

        type(decomp_info), pointer :: decomp
        type(derivatives), pointer :: der
        type(filters),     pointer :: fil

        integer :: nfils
        real(rkind) :: dx, dy, dz

        logical :: use_gTg

    contains

        procedure          :: init
        procedure          :: get_viscosities
        procedure          :: get_conductivity
        procedure          :: get_diffusivity
        procedure          :: get_diff_g
        procedure          :: get_diff_pe
        procedure, private :: filter
        final              :: destroy

    end type

contains

    subroutine init(this,decomp,der,fil,nfils,dx,dy,dz,Cbeta,CbetaP,Cmu,Ckap,CkapP,Cdiff,CY,Cdiff_g,Cdiff_gt,Cdiff_gp,Cdiff_pe,Cdiff_pe_2)
        class(ladobject),        intent(inout) :: this
        type(decomp_info), target, intent(in) :: decomp
        type(derivatives), target, intent(in) :: der
        type(filters),     target, intent(in) :: fil
        integer,           intent(in) :: nfils
        real(rkind),       intent(in) :: Cbeta,CbetaP,Cmu,Ckap,CkapP,Cdiff,CY,Cdiff_g,Cdiff_gt,Cdiff_gp,Cdiff_pe,Cdiff_pe_2,dx,dy,dz    

        ! Set all coefficients
        this%Cbeta = Cbeta
        this%CbetaP = CbetaP
        this%Cmu = Cmu
        this%Ckap = Ckap
        this%CkapP = CkapP
        this%Cdiff = Cdiff
        this%CY = CY
        this%Cdiff_g = Cdiff_g
        this%Cdiff_gt = Cdiff_gt
        this%Cdiff_gp = Cdiff_gp
        this%Cdiff_pe = Cdiff_pe
        this%Cdiff_pe_2 = Cdiff_pe_2

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


    pure elemental subroutine destroy(this)
        type(ladobject),        intent(inout) :: this

        nullify(this%fil)
        nullify(this%der)
        nullify(this%decomp)

    end subroutine

    subroutine get_viscosities(this,rho,p,sos,duidxj,mu,bulk,x_bc,y_bc,z_bc,dt,pfloor)
        class(ladobject),        intent(in) :: this
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)),           intent(in)  :: rho,p,sos
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),9), target, intent(in)  :: duidxj
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)),           intent(inout) :: mu, bulk
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), intent(in) :: pfloor

        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: mustar, bulkstar
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xtmp1,xtmp2
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5,func,bulkP
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: ztmp1,ztmp2

        real(rkind), intent(in) :: dt

        integer :: i,j,k

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
        call gradient(this%decomp,this%der,rho,ytmp1,ytmp2,ytmp3,x_bc,y_bc,z_bc) ! Does not use any Y buffers
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

        ! original
        ! Calculate the switching function
        ytmp1 = ytmp2 / (ytmp2 + ytmp4 + real(1.0D-32,rkind)) ! Switching function f_sw
        where (func .GE. zero)
            ytmp1 = zero
        end where

        ! ! new
        ! ! Calculate the switching function
        ! ytmp1 = one
        ! !ytmp1 = ytmp2 / (ytmp2 + ytmp4 + real(1.0D-32,rkind)) ! Switching function f_sw
        ! !where ((func .GE. zero).AND.(p.gt.1.D1*pfloor))
        ! !   ytmp1 = zero
        ! !end where

        bulkstar = this%Cbeta*rho*ytmp1*abs(bulkstar)



        ! !####################################################################

        ! if(this%CbetaP.gt.zero) then
        !    do k=1,this%decomp%ysz(3)
        !       do j=1,this%decomp%ysz(2)
        !          do i=1,this%decomp%ysz(1)
        !             !bulkstar(i,j,k) = max(bulkstar(i,j,k) , this%CbetaP * max((1.0/max(p(i,j,k),1.0D-32)-pfloor) - (1.0/pfloor-pfloor),zero))
        !             bulkstar(i,j,k) = max(bulkstar(i,j,k) , this%CbetaP * max((1.0-p(i,j,k)/sqrt(max(pfloor,1.0D-32))),0.0))
        !          enddo
        !       enddo
        !    enddo
        ! endif
        


        ! !test modifications -- active only for CbetaP > 1.0D-32 : not default on

        ! !un-normalized bulkstar by local density
        ! if(this%CbetaP.gt.1.0D-32) then
        !    bulkstar = abs(bulkstar) / rho
        ! endif


        ! !New for pressure based bulk viscosity

        ! func = p

        ! ! Step 2: Get 4th derivative in X
        ! call transpose_y_to_x(func,xtmp1,this%decomp)
        ! call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
        ! call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
        ! xtmp2 = xtmp1*this%dx**4
        ! call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        ! bulkP = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )!**2

        ! ! Step 3: Get 4th derivative in Z
        ! call transpose_y_to_z(func,ztmp1,this%decomp)
        ! call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
        ! call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
        ! ztmp2 = ztmp1*this%dz**4
        ! call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        ! bulkP = bulkP + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )!**2

        ! ! Step 4: Get 4th derivative in Y
        ! call this%der%d2dy2(func,ytmp4,y_bc(1),y_bc(2))
        ! call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
        ! ytmp4 = ytmp5*this%dy**4
        ! bulkP = bulkP + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )!**2


        ! !bulkstar = bulkstar + this%CbetaP*rho*ytmp1*abs(bulkP)/(abs(p)*dt)
        ! do k=1,this%decomp%ysz(3)
        !    do j=1,this%decomp%ysz(2)
        !       do i=1,this%decomp%ysz(1)
        !         ! bulkstar(i,j,k) = max(bulkstar(i,j,k), this%CbetaP*rho(i,j,k)*abs(bulkP(i,j,k)) *sos(i,j,k)/max(abs(func(i,j,k))**4,1.0D-32))
        !          bulkstar(i,j,k) = max(bulkstar(i,j,k), this%CbetaP*abs(bulkP(i,j,k)) *sos(i,j,k)/max(abs(func(i,j,k))**one,1.0D-32))
        !       enddo
        !    enddo
        ! enddo

        ! !####################################################################

        ! Filter bulkstar
        call this%filter(bulkstar, x_bc, y_bc, z_bc)

        bulk = bulk + bulkstar
    end subroutine

    subroutine get_conductivity(this,rho,p,e,T,sos,kap,x_bc,y_bc,z_bc,tfloor)
        class(ladobject),  intent(in) :: this
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(in)  :: rho,e,T,sos,p
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(inout) :: kap
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: kapstar
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xtmp1,xtmp2
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5,kapP,func
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: ztmp1,ztmp2
        real(rkind), intent(in) :: tfloor

        integer :: i,j,k

        ! -------- Artificial Conductivity --------

        ! Step 1: Get components of grad(e) squared individually
        call gradient(this%decomp,this%der,e,ytmp1,ytmp2,ytmp3,x_bc,y_bc,z_bc) ! Does not use any Y buffers
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

        ! ! !####################################################
        ! ! !New for pressure based thermal diffusivity

        ! ! if(this%CkapP.gt.zero) then
        ! !    do k=1,this%decomp%ysz(3)
        ! !       do j=1,this%decomp%ysz(2)
        ! !          do i=1,this%decomp%ysz(1)
        ! !             kapstar(i,j,k) = max(kapstar(i,j,k) , this%CkapP * max((1.0-T(i,j,k)/sqrt(max(tfloor,1.0D-32))),0.0))
        ! !          enddo
        ! !       enddo
        ! !    enddo
        ! ! endif


        ! !func = T
        ! do k=1,this%decomp%ysz(3)
        !    do j=1,this%decomp%ysz(2)
        !       do i=1,this%decomp%ysz(1)
        !          func(i,j,k) = one/max(T(i,j,k),tfloor)
        !       enddo
        !    enddo
        ! enddo

        ! ! Step 2: Get 4th derivative in X
        ! call transpose_y_to_x(func,xtmp1,this%decomp)
        ! call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
        ! call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
        ! xtmp2 = xtmp1*this%dx**4
        ! call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        ! kapP = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! ! Step 3: Get 4th derivative in Z
        ! call transpose_y_to_z(func,ztmp1,this%decomp)
        ! call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
        ! call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
        ! ztmp2 = ztmp1*this%dz**4
        ! call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        ! kapP = kapP + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! ! Step 4: Get 4th derivative in Y
        ! call this%der%d2dy2(func,ytmp4,y_bc(1),y_bc(2))
        ! call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
        ! ytmp4 = ytmp5*this%dy**4
        ! kapP = kapP + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero


        ! do k=1,this%decomp%ysz(3)
        !    do j=1,this%decomp%ysz(2)
        !       do i=1,this%decomp%ysz(1)
        !          !kapstar(i,j,k) = max(kapstar(i,j,k), this%CkapP*rho(i,j,k)*e(i,j,k)*sos(i,j,k)*abs(kapP(i,j,k))/max(func(i,j,k),1.0D-32)**2)
        !          kapstar(i,j,k) = max(kapstar(i,j,k), this%CkapP*rho(i,j,k)*e(i,j,k)*sos(i,j,k)*abs(kapP(i,j,k)))
        !       enddo
        !    enddo
        ! enddo

        ! ! !####################################################

        ! Filter kapstar
        call this%filter(kapstar, x_bc, y_bc, z_bc)

        kap = kap + kapstar
    end subroutine

    subroutine get_diffusivity(this,Ys,dYsdx,dYsdy,dYsdz,sos,diff,x_bc,y_bc,z_bc)
        class(ladobject),  intent(in) :: this
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(in)  :: Ys,dYsdx,dYsdy,dYsdz,sos
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


    subroutine get_diff_pe(this,pe,sos,diff_pe,x_bc,y_bc,z_bc)
        class(ladobject),  intent(in) :: this
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(in)  :: pe,sos
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(inout) :: diff_pe
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: diffstar
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xtmp1,xtmp2
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: ztmp1,ztmp2

        ! -------- Artificial plastic entropy (pe) Diffusivity ---------

        ! Step 1: Get components of grad(e_g) squared individually
        call gradient(this%decomp,this%der,pe,ytmp1,ytmp2,ytmp3,x_bc,y_bc,z_bc) ! Does not use any Y buffers
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(pe,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
        xtmp2 = xtmp1*this%dx**4
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        diffstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(pe,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
        ztmp2 = ztmp1*this%dz**4
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        diffstar = diffstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(pe,ytmp4,y_bc(1),y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
        ytmp4 = ytmp5*this%dy**4
        diffstar = diffstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        diffstar = this%Cdiff_pe*sos*abs(diffstar) ! CD part of diff

        ytmp4 = (sqrt(ytmp1)*this%dx + sqrt(ytmp2)*this%dy + sqrt(ytmp3)*this%dz) &
                        / ( sqrt(ytmp1+ytmp2+ytmp3) + real(1.0D-32,rkind) ) ! grid scale
        ytmp5 = this%Cdiff_pe_2*sos*( abs(min(pe,0.0d0)) )*ytmp4 ! CY part of diff

        ! Filter each part
        call this%filter(diffstar, x_bc, y_bc, z_bc)
        call this%filter(ytmp5, x_bc, y_bc, z_bc)

        diffstar = max(diffstar, ytmp5) ! Take max of both terms instead of add to minimize the dissipation

        diff_pe = diff_pe + diffstar
    end subroutine


    !full-component version
    ! subroutine get_diff_g(this,g,g_t,diff_g,diff_gt,x_bc,y_bc,z_bc)
    !     class(ladobject),  intent(in) :: this
    !     !real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(in)  :: rho,e,T,sos
    !     real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),9), intent(in) :: g,g_t
    !     real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),9), intent(inout) :: diff_g,diff_gt
    !     integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
    !     !logical, intent(in) :: strainHard

    !     real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: diffstar_g
    !     real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xtmp1,xtmp2
    !     real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5
    !     real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: ztmp1,ztmp2

    !     integer :: n


    !     ! -------- g_e diffusivity --------
    !     do n=1,9

    !        ! Step 1: Get components of grad(g(:,:,:,n)) squared individually
    !        call gradient(this%decomp,this%der,g(:,:,:,n),ytmp1,ytmp2,ytmp3,x_bc,y_bc,z_bc) ! Does not use any Y buffers
    !        ytmp1 = ytmp1*ytmp1
    !        ytmp2 = ytmp2*ytmp2
    !        ytmp3 = ytmp3*ytmp3
           
    !        ! Step 2: Get 4th derivative in X
    !        call transpose_y_to_x(g(:,:,:,n),xtmp1,this%decomp)
    !        call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
    !        call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
    !        xtmp2 = xtmp1*this%dx**4
    !        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
    !        diffstar_g = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
           
    !        ! Step 3: Get 4th derivative in Z
    !        call transpose_y_to_z(g(:,:,:,n),ztmp1,this%decomp)
    !        call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
    !        call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
    !        ztmp2 = ztmp1*this%dz**4
    !        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
    !        diffstar_g = diffstar_g + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
           
    !        ! Step 4: Get 4th derivative in Y
    !        call this%der%d2dy2(g(:,:,:,n),ytmp4,y_bc(1),y_bc(2))
    !        call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
    !        ytmp4 = ytmp5*this%dy**4
    !        diffstar_g = diffstar_g + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero
           
    !        diffstar_g = this%Cdiff_g*abs(diffstar_g)
           
    !        ! Filter diffstar_g
    !        call this%filter(diffstar_g, x_bc, y_bc, z_bc)
           
    !        diff_g(:,:,:,n) = diff_g(:,:,:,n) + diffstar_g

    !     enddo


    !     ! -------- g_t diffusivity --------
    !     if(this%strainHard) then
    !        do n=1,9

    !           ! Step 1: Get components of grad(g(:,:,:,n)) squared individually
    !           call gradient(this%decomp,this%der,g_t(:,:,:,n),ytmp1,ytmp2,ytmp3,x_bc,y_bc,z_bc) ! Does not use any Y buffers
    !           ytmp1 = ytmp1*ytmp1
    !           ytmp2 = ytmp2*ytmp2
    !           ytmp3 = ytmp3*ytmp3

    !           ! Step 2: Get 4th derivative in X
    !           call transpose_y_to_x(g_t(:,:,:,n),xtmp1,this%decomp)
    !           call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
    !           call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
    !           xtmp2 = xtmp1*this%dx**4
    !           call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
    !           diffstar_g = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

    !           ! Step 3: Get 4th derivative in Z
    !           call transpose_y_to_z(g_t(:,:,:,n),ztmp1,this%decomp)
    !           call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
    !           call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
    !           ztmp2 = ztmp1*this%dz**4
    !           call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
    !           diffstar_g = diffstar_g + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

    !           ! Step 4: Get 4th derivative in Y
    !           call this%der%d2dy2(g_t(:,:,:,n),ytmp4,y_bc(1),y_bc(2))
    !           call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
    !           ytmp4 = ytmp5*this%dy**4
    !           diffstar_g = diffstar_g + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

    !           diffstar_g = this%Cdiff_gt*abs(diffstar_g)

    !           ! Filter diffstar_g
    !           call this%filter(diffstar_g, x_bc, y_bc, z_bc)

    !           diff_gt(:,:,:,n) = diff_gt(:,:,:,n) + diffstar_g

    !        enddo
    !     endif


    ! end subroutine



    ! Eulerian-Almansi strain version
    subroutine get_diff_g(this,g,g_t,g_p,sos,diff_g,diff_gt,diff_gp,use_gTg,strainHard, x_bc,y_bc,z_bc)
        class(ladobject),  intent(in) :: this
        !real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(in)   :: rho,e,T,sos
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),9), intent(in)  :: g,g_t,g_p
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(in)    :: sos
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(inout) :: diff_g,diff_gt,diff_gp
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: diffstar_g,e_g,e_gt,e_gp
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xtmp1,xtmp2
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ytmp1,ytmp2,ytmp3,ytmp4,ytmp5
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: ztmp1,ztmp2

        integer :: n

        logical, intent(in) :: use_gTg,strainHard

        diff_g  = zero
        diff_gt = zero
        diff_gp = zero

        ! -------- g_e diffusivity --------
        ! !Eulerian-Almansi strain tensor: ea = (I-(g_e)^T.g_e)/2 --- not explicitly calculated
        ! !strain norm: e_g = sqrt(2/3*ea_ij*ea_ij)
        !if (use_gTg) then
        if ((use_gTg).and.(.not. strainHard)) then
           e_g = sqrt(( (1.0d0 - g(:,:,:,1))**2 + g(:,:,:,2)**2 + g(:,:,:,3)**2 + g(:,:,:,4)**2 + (1.0d0 - g(:,:,:,5))**2 + g(:,:,:,6)**2 + g(:,:,:,7)**2 + g(:,:,:,8)**2 + (1.0d0 - g(:,:,:,9))**2 )/6.0d0)
        else
           e_g = sqrt( ( (1.0d0 - g(:,:,:,1)**2 - g(:,:,:,4)**2 - g(:,:,:,7)**2)**2 + (1.0d0 - g(:,:,:,2)**2 - g(:,:,:,5)**2 - g(:,:,:,8)**2)**2 + (1.0d0 - g(:,:,:,3)**2 - g(:,:,:,6)**2 - g(:,:,:,9)**2)**2 )/6.0d0 + ( (-g(:,:,:,1)*g(:,:,:,2) - g(:,:,:,4)*g(:,:,:,5) - g(:,:,:,7)*g(:,:,:,8))**2 + (-g(:,:,:,1)*g(:,:,:,3) - g(:,:,:,4)*g(:,:,:,6) - g(:,:,:,7)*g(:,:,:,9))**2 + (-g(:,:,:,2)*g(:,:,:,3) - g(:,:,:,5)*g(:,:,:,6) - g(:,:,:,8)*g(:,:,:,9))**2 )/3.0d0 )
        endif

        
        ! Step 1: Get components of grad(e_g) squared individually
        call gradient(this%decomp,this%der,e_g,ytmp1,ytmp2,ytmp3,x_bc,y_bc,z_bc) ! Does not use any Y buffers
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(e_g,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
        xtmp2 = xtmp1*this%dx**4
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        diffstar_g = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(e_g,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
        ztmp2 = ztmp1*this%dz**4
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        diffstar_g = diffstar_g + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(e_g,ytmp4,y_bc(1),y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
        ytmp4 = ytmp5*this%dy**4
        diffstar_g = diffstar_g + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        diffstar_g = this%Cdiff_g*sos*abs(diffstar_g)

        ! Filter diffstar_g
        call this%filter(diffstar_g, x_bc, y_bc, z_bc)

        diff_g = diff_g + diffstar_g

        
        ! -------- g_t diffusivity --------
        ! !Eulerian-Almansi strain tensor: ea = (I-(g_t)^T.g_t)/2 --- not explicitly calculated
        ! !strain norm: e_gt = sqrt(2/3*ea_ij*ea_ij)
        if(strainHard) then
           if (use_gTg) then
              e_gt = sqrt(( (1.0d0 - g_t(:,:,:,1))**2 + g_t(:,:,:,2)**2 + g_t(:,:,:,3)**2 + g_t(:,:,:,4)**2 + (1.0d0 - g_t(:,:,:,5))**2 + g_t(:,:,:,6)**2 + g_t(:,:,:,7)**2 + g_t(:,:,:,8)**2 + (1.0d0 - g_t(:,:,:,9))**2 )/6.0d0)
           else
              e_gt = sqrt( ( (1.0d0 - g_t(:,:,:,1)**2 - g_t(:,:,:,4)**2 - g_t(:,:,:,7)**2)**2 + (1.0d0 - g_t(:,:,:,2)**2 - g_t(:,:,:,5)**2 - g_t(:,:,:,8)**2)**2 + (1.0d0 - g_t(:,:,:,3)**2 - g_t(:,:,:,6)**2 - g_t(:,:,:,9)**2)**2 )/6.0d0 + ( (-g_t(:,:,:,1)*g_t(:,:,:,2) - g_t(:,:,:,4)*g_t(:,:,:,5) - g_t(:,:,:,7)*g_t(:,:,:,8))**2 + (-g_t(:,:,:,1)*g_t(:,:,:,3) - g_t(:,:,:,4)*g_t(:,:,:,6) - g_t(:,:,:,7)*g_t(:,:,:,9))**2 + (-g_t(:,:,:,2)*g_t(:,:,:,3) - g_t(:,:,:,5)*g_t(:,:,:,6) - g_t(:,:,:,8)*g_t(:,:,:,9))**2 )/3.0d0 )
           endif


           ! Step 1: Get components of grad(e_gt) squared individually
           call gradient(this%decomp,this%der,e_gt,ytmp1,ytmp2,ytmp3,x_bc,y_bc,z_bc) ! Does not use any Y buffers
           ytmp1 = ytmp1*ytmp1
           ytmp2 = ytmp2*ytmp2
           ytmp3 = ytmp3*ytmp3

           ! Step 2: Get 4th derivative in X
           call transpose_y_to_x(e_gt,xtmp1,this%decomp)
           call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
           call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
           xtmp2 = xtmp1*this%dx**4
           call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
           diffstar_g = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

           ! Step 3: Get 4th derivative in Z
           call transpose_y_to_z(e_gt,ztmp1,this%decomp)
           call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
           call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
           ztmp2 = ztmp1*this%dz**4
           call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
           diffstar_g = diffstar_g + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

           ! Step 4: Get 4th derivative in Y
           call this%der%d2dy2(e_gt,ytmp4,y_bc(1),y_bc(2))
           call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
           ytmp4 = ytmp5*this%dy**4
           diffstar_g = diffstar_g + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

           diffstar_g = this%Cdiff_gt*sos*abs(diffstar_g)

           ! Filter diffstar_g
           call this%filter(diffstar_g, x_bc, y_bc, z_bc)

           diff_gt = diff_gt + diffstar_g


           !!print*,this%Cdiff_g,this%Cdiff_gt,maxval(sos)
           
           ! print '(6(E14.7E3,3X))',maxval(diff_g),maxval(diff_gt),maxval(g),maxval(g_t),maxval(e_g),maxval(e_gt)




        ! -------- g_p diffusivity --------
        ! !Eulerian-Almansi strain tensor: ea = (I-(g_p)^T.g_p)/2 --- not explicitly calculated
        ! !strain norm: e_gp = sqrt(2/3*ea_ij*ea_ij)

           if (use_gTg) then
              e_gp = sqrt(( (1.0d0 - g_p(:,:,:,1))**2 + g_p(:,:,:,2)**2 + g_p(:,:,:,3)**2 + g_p(:,:,:,4)**2 + (1.0d0 - g_p(:,:,:,5))**2 + g_p(:,:,:,6)**2 + g_p(:,:,:,7)**2 + g_p(:,:,:,8)**2 + (1.0d0 - g_p(:,:,:,9))**2 )/6.0d0)
           else
              e_gp = sqrt( ( (1.0d0 - g_p(:,:,:,1)**2 - g_p(:,:,:,4)**2 - g_p(:,:,:,7)**2)**2 + (1.0d0 - g_p(:,:,:,2)**2 - g_p(:,:,:,5)**2 - g_p(:,:,:,8)**2)**2 + (1.0d0 - g_p(:,:,:,3)**2 - g_p(:,:,:,6)**2 - g_p(:,:,:,9)**2)**2 )/6.0d0 + ( (-g_p(:,:,:,1)*g_p(:,:,:,2) - g_p(:,:,:,4)*g_p(:,:,:,5) - g_p(:,:,:,7)*g_p(:,:,:,8))**2 + (-g_p(:,:,:,1)*g_p(:,:,:,3) - g_p(:,:,:,4)*g_p(:,:,:,6) - g_p(:,:,:,7)*g_p(:,:,:,9))**2 + (-g_p(:,:,:,2)*g_p(:,:,:,3) - g_p(:,:,:,5)*g_p(:,:,:,6) - g_p(:,:,:,8)*g_p(:,:,:,9))**2 )/3.0d0 )
           endif


           ! Step 1: Get components of grad(e_gp) squared individually
           call gradient(this%decomp,this%der,e_gp,ytmp1,ytmp2,ytmp3,x_bc,y_bc,z_bc) ! Does not use any Y buffers
           ytmp1 = ytmp1*ytmp1
           ytmp2 = ytmp2*ytmp2
           ytmp3 = ytmp3*ytmp3

           ! Step 2: Get 4th derivative in X
           call transpose_y_to_x(e_gp,xtmp1,this%decomp)
           call this%der%d2dx2(xtmp1,xtmp2,x_bc(1),x_bc(2))
           call this%der%d2dx2(xtmp2,xtmp1,x_bc(1),x_bc(2))
           xtmp2 = xtmp1*this%dx**4
           call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
           diffstar_g = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

           ! Step 3: Get 4th derivative in Z
           call transpose_y_to_z(e_gp,ztmp1,this%decomp)
           call this%der%d2dz2(ztmp1,ztmp2,z_bc(1),z_bc(2))
           call this%der%d2dz2(ztmp2,ztmp1,z_bc(1),z_bc(2))
           ztmp2 = ztmp1*this%dz**4
           call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
           diffstar_g = diffstar_g + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

           ! Step 4: Get 4th derivative in Y
           call this%der%d2dy2(e_gp,ytmp4,y_bc(1),y_bc(2))
           call this%der%d2dy2(ytmp4,ytmp5,y_bc(1),y_bc(2))
           ytmp4 = ytmp5*this%dy**4
           diffstar_g = diffstar_g + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

           diffstar_g = this%Cdiff_gp*sos*abs(diffstar_g)

           ! Filter diffstar_g
           call this%filter(diffstar_g, x_bc, y_bc, z_bc)

           diff_gp = diff_gp + diffstar_g
        endif

           !!print*,this%Cdiff_g,this%Cdiff_gt,maxval(sos)
           
           ! print '(6(E14.7E3,3X))',maxval(diff_g),maxval(diff_gt),maxval(g),maxval(g_t),maxval(e_g),maxval(e_gt)


    end subroutine



    subroutine filter(this,arr,x_bc,y_bc,z_bc)
        class(ladobject), intent(in) :: this
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)), intent(inout) :: arr
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: tmp_in_y
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: tmp1_in_x, tmp2_in_x
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: tmp1_in_z, tmp2_in_z
       
        integer :: idx

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
