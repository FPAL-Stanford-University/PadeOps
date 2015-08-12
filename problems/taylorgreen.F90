module hooks

    use kind_parameters, only: rkind
    use constants,       only: two,pi
    implicit none

contains

    subroutine meshgen(cgp)
        class(cgrid), intent(inout) :: cgp
        integer :: i,j,k

        associate( nx => cgp%nx, ny => cgp%ny, nz => cgp%nz, &
                   dx => cgp%dx, dy => cgp%dy, dz => cgp%dz, &
                   x => cgp%mesh(:,:,:,1), y => cgp%mesh(:,:,:,2), z => cgp%mesh(:,:,:,3), &
                   proc_i1 => cgp%decomp%yst(1), proc_in => cgp%decomp%yen(1), proc_ni => cgp%decomp%ysz(1), & 
                   proc_j1 => cgp%decomp%yst(2), proc_jn => cgp%decomp%yen(2), proc_nj => cgp%decomp%ysz(2), & 
                   proc_k1 => cgp%decomp%yst(3), proc_kn => cgp%decomp%yen(3), proc_nk => cgp%decomp%ysz(3) )

            ! Create mesh from [0,2*pi)x[0,2*pi)x[0,2*pi) using nx, ny, nz points in x, y and z respectively
            ! Need to set x, y and z as well as  dx, dy and dz

            dx = two*pi/real(nx,rkind)
            dy = two*pi/real(ny,rkind)
            dz = two*pi/real(nz,rkind)

            do k=1,proc_nk
                do j=1,proc_nj
                    do i=1,proc_ni
                        x(i,j,k) = real( proc_i1 + i - 1, rkind ) * dx
                        y(i,j,k) = real( proc_j1 + j - 1, rkind ) * dy
                        z(i,j,k) = real( proc_k1 + k - 1, rkind ) * dz
                    end do
                end do
            end do

        end associate

    end subroutine

    subroutine initfields(cgp)
        class(cgrid), intent(inout) :: cgp
        integer :: i,j,k

    end subroutine

end module

program taylorgreen

    use kind_parameters,  only: clen,stdout
    use CompressibleGrid, only: cgrid
    implicit none

    type(cgrid) :: cgp
    character(len=clen) :: inputfile
    integer :: ierr

    ! Start MPI
    call MPI_Init(ierr)

    call GETARG(1,inputfile)

    ! Initialize the grid object
    call cgp%init(inputfile)

    associate( x => cgp%mesh(:,:,:,1) )
        write(stdout,*) x(:,1,1)
    end associate

    ! Destroy everythin before ending
    call cgp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
