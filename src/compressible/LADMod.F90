module LADMod

    use kind_parameters, only: rkind
    use constants,       only: half,one
    use DerivativesMod,  only: derivatives

    implicit none

    type :: lad
        real(rkind) :: Cbeta
        real(rkind) :: Cmu
        real(rkind) :: Ckap
        real(rkind) :: Cdiff
        real(rkind) :: CY
    contains
        procedure :: init
        procedure :: get_viscosities
        procedure :: get_conductivity
        procedure :: get_diffusivity
        procedure :: destroy

contains

    subroutine get_viscosities(this,duidxj,mustar,bulkstar)    
        ! -------- Artificial Shear Viscosity --------

        ! Magnitude of strain rate
        func = sqrt(dudx**2 + half*(dvdx+dudy)**2 + half*(dwdx+dudz)**2 &
                            +             dvdy**2 + half*(dwdy+dvdz)**2 &
                                                  +             dwdz**2 )
        
        ! Get 4th derivative in X
        call transpose_y_to_x(func,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,this%x_bc(1),this%x_bc(2))
        xtmp2 = xtmp1*this%dx**6
        call transpose_x_to_y(xtmp2,mustar,this%decomp)
        
        ! Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
        ztmp2 = ztmp1*this%dz**6
        call transpose_z_to_y(ztmp2,ytmp1,this%decomp)
        mustar = mustar + ytmp1
        
        ! Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp1,this%y_bc(1),this%y_bc(2))
        call this%der%d2dy2(ytmp1,ytmp2,this%y_bc(1),this%y_bc(2))
        ytmp1 = ytmp2*this%dy**6
        mustar = mustar + ytmp1

        mustar = this%Cmu*this%rho*abs(mustar)
        
        ! Filter mustar
        call this%filter(mustar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)
        
        ! -------- Artificial Bulk Viscosity --------
        
        func = dudx + dvdy + dwdz      ! dilatation
        
        ! Step 1: Get components of grad(rho) squared individually
        call this%gradient(this%rho,ytmp1,ytmp2,ytmp3,this%x_bc,this%y_bc,this%z_bc) ! Does not use any Y buffers
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(func,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,this%x_bc(1),this%x_bc(2))
        xtmp2 = xtmp1*this%dx**4
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        bulkstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(func,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
        ztmp2 = ztmp1*this%dz**4
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        bulkstar = bulkstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) )**2

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(func,ytmp4,this%y_bc(1),this%y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,this%y_bc(1),this%y_bc(2))
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

        bulkstar = this%Cbeta*this%rho*ytmp1*abs(bulkstar)

        ! Filter bulkstar
        call this%filter(bulkstar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)

    end subroutine

    subroutine get_conductivity(this,e)
        ! -------- Artificial Conductivity --------

        ! Step 1: Get components of grad(e) squared individually
        call this%gradient(this%e,ytmp1,ytmp2,ytmp3,this%x_bc,this%y_bc,this%z_bc) ! Does not use any Y buffers
        ytmp1 = ytmp1*ytmp1
        ytmp2 = ytmp2*ytmp2
        ytmp3 = ytmp3*ytmp3

        ! Step 2: Get 4th derivative in X
        call transpose_y_to_x(this%e,xtmp1,this%decomp)
        call this%der%d2dx2(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
        call this%der%d2dx2(xtmp2,xtmp1,this%x_bc(1),this%x_bc(2))
        xtmp2 = xtmp1*this%dx**4
        call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
        kapstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 3: Get 4th derivative in Z
        call transpose_y_to_z(this%e,ztmp1,this%decomp)
        call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
        call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
        ztmp2 = ztmp1*this%dz**4
        call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
        kapstar = kapstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Step 4: Get 4th derivative in Y
        call this%der%d2dy2(this%e,ytmp4,this%y_bc(1),this%y_bc(2))
        call this%der%d2dy2(ytmp4,ytmp5,this%y_bc(1),this%y_bc(2))
        ytmp4 = ytmp5*this%dy**4
        kapstar = kapstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

        ! Now, all ytmps are free to use
        call this%mix%get_sos(this%rho,this%p,ytmp1)  ! Speed of sound

        kapstar = this%Ckap*this%rho*ytmp1*abs(kapstar)/this%T

        ! Filter kapstar
        call this%filter(kapstar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)

    end subroutine

    subroutine get_diffusivity(this,Ys,dYsdx,dYsdy,dYsdz)
        ! -------- Artificial Diffusivity ---------
        if (this%mix%ns .GT. 1) then
            do i = 1,this%mix%ns
                ! Step 1: Get components of grad(Ys) squared individually
                ytmp1 = dYsdx(:,:,:,i)*dYsdx(:,:,:,i)
                ytmp2 = dYsdy(:,:,:,i)*dYsdy(:,:,:,i)
                ytmp3 = dYsdz(:,:,:,i)*dYsdz(:,:,:,i)

                call this%mix%get_sos(this%rho,this%p,ytmp6)  ! Speed of sound

                ! Step 2: Get 4th derivative in X
                call transpose_y_to_x(this%Ys(:,:,:,i),xtmp1,this%decomp)
                call this%der%d2dx2(xtmp1,xtmp2,this%x_bc(1),this%x_bc(2))
                call this%der%d2dx2(xtmp2,xtmp1,this%x_bc(1),this%x_bc(2))
                xtmp2 = xtmp1*this%dx**4
                call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
                diffstar = ytmp4 * ( this%dx * ytmp1 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

                ! Step 3: Get 4th derivative in Z
                call transpose_y_to_z(this%Ys(:,:,:,i),ztmp1,this%decomp)
                call this%der%d2dz2(ztmp1,ztmp2,this%z_bc(1),this%z_bc(2))
                call this%der%d2dz2(ztmp2,ztmp1,this%z_bc(1),this%z_bc(2))
                ztmp2 = ztmp1*this%dz**4
                call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
                diffstar = diffstar + ytmp4 * ( this%dz * ytmp3 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero

                ! Step 4: Get 4th derivative in Y
                call this%der%d2dy2(this%Ys(:,:,:,i),ytmp4,this%y_bc(1),this%y_bc(2))
                call this%der%d2dy2(ytmp4,ytmp5,this%y_bc(1),this%y_bc(2))
                ytmp4 = ytmp5*this%dy**4
                diffstar = diffstar + ytmp4 * ( this%dy * ytmp2 / (ytmp1 + ytmp2 + ytmp3 + real(1.0D-32,rkind)) ) ! Add eps in case denominator is zero


                diffstar = this%Cdiff*ytmp6*abs(diffstar) ! CD part of diff

                ytmp4 = (sqrt(ytmp1)*this%dx + sqrt(ytmp2)*this%dy + sqrt(ytmp3)*this%dz) &
                                / ( sqrt(ytmp1+ytmp2+ytmp3) + real(1.0D-32,rkind) ) ! grid scale
                ytmp5 = this%CY*ytmp6*( half*(abs(this%Ys(:,:,:,i))-one + abs(this%Ys(:,:,:,i)-one)) )*ytmp4 ! CY part of diff

                diffstar = max(diffstar, ytmp5) ! Take max of both terms instead of add to minimize the dissipation

                ! Filter diffstar
                call this%filter(diffstar, this%gfil, 2, this%x_bc, this%y_bc, this%z_bc)

                ! Add to physical diffusivity
                this%diff(:,:,:,i) = this%diff(:,:,:,i) + diffstar
            end do
        end if
    end subroutine

        ! Now, add to physical fluid properties
        this%mu   = this%mu   + mustar
        this%bulk = this%bulk + bulkstar
        this%kap  = this%kap  + kapstar
end module
