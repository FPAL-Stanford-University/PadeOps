module GaborModeMod
    use constants
    use kind_parameters, only: rkind, clen
    use fortran_assert, only: assert
    use decomp_2d
    
    implicit none 
    logical :: useStrain = .true.
    real(rkind) :: AnuRNG = 4.0d-4 ! spectral eddy viscosity coefficient
    real(rkind) :: nu = 0.d0       ! kinematic viscosity of working fluid
    logical :: MMS = .false.       ! Do method of manufactured solution test
    real(rkind) :: nutMMS = 0.1d0  ! Used if MMS = .true.
    type :: gaborMode 
        !private   
        real(rkind) :: uhatR, uhatI, vhatR, vhatI, whatR, whatI, kx, ky, kz, x, y, z
        real(rkind) :: Wx, Wy, Wz ! Window function support widths
         
        contains 
            procedure          :: init 
            procedure          :: evolve
            procedure, private :: RK4
            procedure, private :: periodicBClocation
            procedure, private :: projectDivergenceFree
            procedure          :: render 
    end type 

contains 

#include "GaborMode_files/timeSteppingStuff.F90"

    subroutine init(this, uhat, vhat, what, x, y, z, kx, ky, kz, Wx, Wy, Wz)
        class(gaborMode), intent(inout) :: this
        complex(kind=8), intent(in) :: uhat, vhat, what
        real(rkind), intent(in) :: x, y, z, kx, ky, kz
        real(rkind), intent(in) :: Wx, Wy, Wz

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

        this%Wx = Wx 
        this%Wy = Wy 
        this%Wz = Wz 
    end subroutine
    
    subroutine evolve(this, ules, vles, wles, dudx, dudy, dudz, dvdx, dvdy, &
                      & dvdz, dwdx, dwdy, dwdz, epsKE, dt, kmin)
        class(gaborMode), intent(inout) :: this 
        real(rkind), intent(in) :: ules, vles, wles, dudx, dudy, dudz, dvdx
        real(rkind), intent(in) :: dvdy, dvdz, dwdx, dwdy, dwdz, dt, epsKE, kmin
        
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
            call this%RK4(AnuRNG, nu, epsKE, dudx, dudy, dudz, dvdx, dvdy, &
              & dvdz, dwdx, dwdy, dwdz, kmin, dt)
            if (.not. MMS) call this%projectDivergenceFree()            
        end if
    end subroutine 

    subroutine render(this, u, v, w, xRange, yRange, zRange, delta)
        ! Compute the real-space velocity contribution from a Gabor Mode
        ! Inputs:
        !   u, v, w, --> velocity field arrays. These use global indices
        !   xRange, yRange, zRange --> The domain bounds for the current MPI rank
        !   delta --> Grid spacing of the Eulerian mesh

        class(gaborMode), intent(in) :: this
        real(rkind), dimension(2), intent(in) :: xRange, yRange, zRange
        real(rkind), dimension(3), intent(in) :: delta  
        real(kind=4), dimension(:,:,:), intent(inout) :: u, v, w
        integer :: xst, xen, yst, yen, zst, zen, i, j, k 
        real(rkind) :: xmin, xmax, ymin, ymax, zmin, zmax 
        real(kind=4) :: Lwx, Lwy, Lwz
        real(kind=4) :: uhatR, vhatR, whatR, uhatI, vhatI, whatI
        real(kind=4) :: kx_sp, ky_sp, kz_sp, kdotx, cs, weight, weight_z, weight_yz
        real(kind=4) :: xshift, yshift, zshift, ss 
      
        xmin = max(this%x - half*this%Wx, xRange(1))
        ymin = max(this%y - half*this%Wy, yRange(1))
        zmin = max(this%z - half*this%Wz, zRange(1))
        
        xmax = min(this%x + half*this%Wx, xRange(2))
        ymax = min(this%y + half*this%Wy, yRange(2))
        zmax = min(this%z + half*this%Wz, zRange(2))
       
        xst = ceiling(((xmin - xRange(1)))/delta(1)) + 1 
        yst = ceiling(((ymin - yRange(1)))/delta(2)) + 1 
        zst = ceiling(((zmin - zRange(1)))/delta(3)) + 1 
        
        xen = floor(((xmax - xRange(1)))/delta(1)) + 1 
        yen = floor(((ymax - yRange(1)))/delta(2)) + 1 
        zen = floor(((zmax - zRange(1)))/delta(3)) + 1 

        Lwx = real(pi/this%Wx ,kind=4)
        Lwy = real(pi/this%Wy ,kind=4)
        Lwz = real(pi/this%Wz ,kind=4)

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
            weight_z = cos(zshift*Lwz)

            do j = yst,yen
                yshift = real( ((j-1)*delta(2) + yRange(1)) - this%y, kind=4 )
                weight_yz = cos(yshift*Lwy)*weight_z
                
                do i = xst,xen
                    xshift = real( ((i-1)*delta(1) + xRange(1)) - this%x, kind=4 )
                    weight = cos(xshift*Lwx)*weight_yz 

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
