module GaborModeMod
    use constants
    use kind_parameters, only: rkind, clen
    use fortran_assert, only: assert
    
    implicit none 
    logical :: useStrain = .true.
    real(rkind) :: AnuRNG = 4.0d-4 ! spectral eddy viscosity coefficient
    real(rkind) :: nu = 0.d0       ! kinematic viscosity of working fluid
    type :: gaborMode 
        !private   
        real(rkind) :: uhatR, uhatI, vhatR, vhatI, whatR, whatI, kx, ky, kz, x, y, z
        real(rkind) :: wSupport
         
        contains 
            procedure          :: init 
            procedure          :: evolve
            procedure, private :: RK4
            procedure, private :: periodicBClocation
            procedure          :: render 
    end type 

contains 

#include "GaborMode_files/timeSteppingStuff.F90"

    subroutine init(this, uhat, vhat, what, x, y, z, kx, ky, kz, wSupport)
        class(gaborMode), intent(inout) :: this
        complex(kind=8), intent(in) :: uhat, vhat, what
        real(rkind), intent(in) :: x, y, z, kx, ky, kz
        real(rkind), intent(in), optional :: wSupport 

        this%uhatR = real(uhat, kind=rkind)
        this%uhatI = dimag(uhat)
        
        this%vhatR = real(vhat, kind=rkind)
        this%vhatI = dimag(vhat)
        
        this%whatR = real(what, kind=rkind)
        this%whatI = dimag(what)

        this%x = x
        this%y = y
        this%z = z

        this%kx = kx
        this%ky = ky
        this%kz = kz

        if (present(wSupport)) then
            this%wSupport = wSupport 
        else
            this%wSupport = 1.d0 ! temporary  
        end if 
    end subroutine
    
    subroutine evolve(this, ules, vles, wles, dudx, dudy, dudz, dvdx, dvdy, &
                      & dvdz, dwdx, dwdy, dwdz, epsKE, dt, kmin)
        class(gaborMode), intent(inout) :: this 
        real(rkind), intent(in) :: ules, vles, wles, dudx, dudy, dudz, dvdx
        real(rkind), intent(in) :: dvdy, dvdz, dwdx, dwdy, dwdz, dt, epsKE, kmin
        real(rkind) :: uhatR, uhatI, vhatR, vhatI, whatR, whatI
        real(rkind) :: ksq, onebyksq
        
        ! Inputs:
        !   ules, vles, wles --> large-scale velocity components at mode location
        !   dudx, ..., dwdz  --> large-scale velocity-gradient components at mode location
        !   epsKE            --> Kinetic energy (KE) dissipation rate
        !   dt               --> time-step size
        !   kmin             --> The minimum wavenumber for enrichment (i.e.
        !                        Nyquist wavenumber of LES to be enriched)

        ! TODO : Finish this 
        this%x = this%x + dt*ules 
        this%y = this%y + dt*vles 
        this%z = this%z + dt*wles
        call this%periodicBClocation()
      
        if(useStrain) then
            ! Step 2: Straining + relaxation
            call this%RK4(AnuRNG, nu, epsKE, dudx, dudy, dudz, dvdx, dvdy, &
              & dvdz, dwdx, dwdy, dwdz, kmin, dt)
            
            ksq = (this%kx*this%kx + this%ky*this%ky + this%kz*this%kz)
            onebyksq = 1.0d0/ksq

            ! Step 4: Projection
            this%uhatR = this%uhatR - onebyksq*(this%kx*this%kx*this%uhatR + &
                          this%ky*this%kx*this%vhatR + this%kz*this%kx*this%whatR)
            this%uhatI = this%uhatI - onebyksq*(this%kx*this%kx*this%uhatI + &
                          this%ky*this%kx*this%vhatI + this%kz*this%kx*this%whatI)

            this%vhatR = this%vhatR - onebyksq*(this%kx*this%ky*this%uhatR + &
                          this%ky*this%ky*this%vhatR + this%kz*this%ky*this%whatR)
            this%vhatI = this%vhatI - onebyksq*(this%kx*this%ky*this%uhatI + &
                          this%ky*this%ky*this%vhatI + this%kz*this%ky*this%whatI)
            
            this%whatR = this%whatR - onebyksq*(this%kx*this%kz*this%uhatR + &
                          this%ky*this%kz*this%vhatR + this%kz*this%kz*this%whatR)
            this%whatI = this%whatI - onebyksq*(this%kx*this%kz*this%uhatI + &
                            this%ky*this%kz*this%vhatI + this%kz*this%kz*this%whatI)
        end if ! if (useStrain)
    end subroutine 

    subroutine render(this, u, v, w, xRange, yRange, zRange, delta)
        class(gaborMode), intent(in) :: this
        real(rkind), dimension(2), intent(in) :: xRange, yRange, zRange
        real(rkind), dimension(3), intent(in) :: delta  
        real(kind=4), dimension(:,:,:), intent(inout) :: u, v, w
        integer :: xst, xen, yst, yen, zst, zen, i, j, k 
        real(rkind) :: xmin, xmax, ymin, ymax, zmin, zmax 
        real(kind=4) :: Lw, uhatR, vhatR, whatR, uhatI, vhatI, whatI
        real(kind=4) :: kx_sp, ky_sp, kz_sp, kdotx, cs, weight, weight_yz
        real(kind=4) :: xshift, yshift, zshift, ss 
      
        xmin = max(this%x - half*this%wSupport, xRange(1))
        ymin = max(this%y - half*this%wSupport, yRange(1))
        zmin = max(this%z - half*this%wSupport, zRange(1))
        
        xmax = min(this%x + half*this%wSupport, xRange(2))
        ymax = min(this%y + half*this%wSupport, yRange(2))
        zmax = min(this%z + half*this%wSupport, zRange(2))
       
        xst = floor(((xmin - xRange(1)))/delta(1)) + 1 
        yst = floor(((ymin - yRange(1)))/delta(2)) + 1 
        zst = floor(((zmin - zRange(1)))/delta(3)) + 1 
        
        xen = floor(((xmax - xRange(1)))/delta(1)) + 1 
        yen = floor(((ymax - yRange(1)))/delta(2)) + 1 
        zen = floor(((zmax - zRange(1)))/delta(3)) + 1 

        Lw = real(pi/this%wSupport ,kind=4)

        uhatR = real(this%uhatR, kind=4)
        vhatR = real(this%vhatR, kind=4)
        whatR = real(this%whatR, kind=4)
        
        uhatI = real(this%uhatI, kind=4)
        vhatI = real(this%vhatI, kind=4)
        whatI = real(this%whatI, kind=4)

        kx_sp = real(this%kx, kind=4)
        ky_sp = real(this%ky, kind=4)
        kz_sp = real(this%kz, kind=4)

        ! Induce velocity 
        do k = zst,zen
            zshift = real( ((k-1)*delta(3) + zRange(1)) - this%z, kind=4 )
            
            do j = yst,yen
                yshift = real( ((j-1)*delta(2) + yRange(1)) - this%y, kind=4 )
                weight_yz = cos(yshift*Lw)*cos(zshift*Lw) 
                
                do i = xst,xen
                    xshift = real( ((i-1)*delta(1) + xRange(1)) - this%x, kind=4 )
                    weight = cos(xshift*Lw)*weight_yz 

                    kdotx = kx_sp*xshift + ky_sp*yshift + kz_sp*zshift 
                    
                    cs = cos(kdotx)
                    ss = sin(kdotx)

                    u(i,j,k) = u(i,j,k) + weight*(uhatR*cs - uhatI*ss)  
                    v(i,j,k) = v(i,j,k) + weight*(vhatR*cs - vhatI*ss) 
                    w(i,j,k) = w(i,j,k) + weight*(whatR*cs - whatI*ss) 
                end do 
            end do 
        end do 

    end subroutine 


end module 
