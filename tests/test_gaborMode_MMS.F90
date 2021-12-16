program test_gaborMode_MMS
    use domainSetup, only: xDom, yDom, zDom
    use gaborModeMod, only : gaborMode, MMS
    use kind_parameters, only: rkind, clen 
    use constants 
    use timer, only: tic, toc 
    use exits, only: message
    use MMS_mod, only: tMMS
    use fortran_assert, only: assert

    real(rkind), dimension(2), parameter :: DomX = [0.d0, two*pi], DomY = [0.d0, two*pi], DomZ = [0.d0, two*pi]
    integer, parameter :: levels = 4
    type(gaborMode) :: gm
    complex(kind=8), parameter :: uhat = (1.d0,1.d0), vhat = (0.d0,0.d0), what = (1.d0,1.d0)
    integer :: n, ierr
    real(rkind) :: small = 1.1d-2
    real(rkind) :: a1, a2, a3, k1, k2, k3
    
    ! Test with these inputs  
    real(rkind), parameter :: wSupport = 1.5d0  
    real(rkind), parameter :: kx = 2.d0, ky = 0.d0, kz = 0.d0 
    real(rkind), parameter :: xMode = 1.d0, yMode = 1.d0, zMode = 1.d0   
    integer, parameter :: modeAtLevel = 4
    real(rkind) :: ules = 1.d0, vles = 2.d0, wles = 3.d0
    real(rkind) :: dudx = 1.1d0, dudy = 1.2d0, dudz = 1.3d0
    real(rkind) :: dvdx = 2.1d0, dvdy = 2.2d0, dvdz = 2.3d0
    real(rkind) :: dwdx = 3.1d0, dwdy = 3.2d0, dwdz = 3.3d0
    real(rkind) :: dt = 0.001, epsKE = 0.d0, kmin = 0.d0
    integer :: nsteps = 1000
    character(len=clen) :: mssg

    call MPI_Init(ierr)

    xDom = DomX; yDom = DomY; zDom = DomZ

    ! Do gabor mode stuff     
    call gm%init(uhat, vhat, what, xMode, yMode, zMode, kx, ky, kz, wSupport)
    
    call tic

    MMS = .true.
    do n = 1,nsteps
        tMMS = real(n,rkind)*dt
        call gm%evolve(ules,vles,wles,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,epsKE,dt,kmin)
        
        ! Analytical solution
        a1 = cos(2.d0*pi*tMMS)
        a2 = sin(2.d0*pi*tMMS)
        a3 = cos(4.d0*pi*tMMS)
        k1 = 2.d0*cos(2.d0*pi*tMMS)
        k2 = 0.5d0*sin(4.d0*pi*tMMS)
        k3 = 1.1d0*sin(2.d0*pi*tMMS)

        ! Compare analytical and numerical solutions
        write(mssg,'(A,F20.16,A,F20.16)')'abs(a1 - gm%uhatR) < small | a1 = ',a1,' & gm%uhatR = ', gm%uhatR
        call assert(abs(a1 - gm%uhatR) < small, trim(mssg))
        call assert(abs(a1 - gm%uhatI) < small, 'abs(a1 - gm%uhatI) < small')
        call assert(abs(a2 - gm%vhatR) < small, 'abs(a2 - gm%vhatR) < small')
        call assert(abs(a2 - gm%vhatI) < small, 'abs(a2 - gm%vhatI) < small')
        call assert(abs(a3 - gm%whatR) < small, 'abs(a3 - gm%whatR) < small')
        call assert(abs(a3 - gm%whatI) < small, 'abs(a3 - gm%whatI) < small')
        call assert(abs(k1 - gm%kx) < small, 'abs(k1 - gm%kx) < small')
        call assert(abs(k2 - gm%ky) < small, 'abs(k2 - gm%ky) < small')
        call assert(abs(k3 - gm%kz) < small, 'abs(k3 - gm%kz) < small')

        if (mod(n,10) == 0) then
          write(mssg,'(A,I4,A,F10.5)')"Time step ", n, '; gm%uhatR: ', gm%uhatR
          call message(trim(mssg))
        end if
    end do
    call MPI_Finalize(ierr)    
end program 
