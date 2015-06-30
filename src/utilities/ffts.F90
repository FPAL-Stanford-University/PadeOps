module ffts

    use kind_parameters, only: rkind
    implicit none

    real(rkind), allocatable, dimension(:)   :: k1, k2, k3        ! Wavenumbers

contains

    function InitFFT(nx, ny, nz) result(ierr)
        integer, intent(in) :: nx,ny,nz
        integer, intent(out) :: ierr

        ! Allocate wavenumbers
        if(allocated( k1 )) deallocate( k1 ); allocate( k1(nx) )
        if(allocated( k2 )) deallocate( k2 ); allocate( k2(ny) )
        if(allocated( k3 )) deallocate( k3 ); allocate( k3(nz) )

        ierr = 0
    end function

end module
