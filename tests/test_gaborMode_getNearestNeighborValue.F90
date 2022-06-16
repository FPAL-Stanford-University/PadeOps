program test_gaborMode_getNearestNeighborValue
    use GaborModeRoutines, only: getNearestNeighborValue 
    use kind_parameters,   only: rkind
    use decomp_2d
    use mpi  
    implicit none 

    real(rkind) :: dx = 1.d0, dy = 1.d0, dz = 1.d0 
    real(rkind) :: x, y, z 
    type(decomp_info) :: gp
    real(rkind), dimension(2,2,2) :: datIn
    real(rkind) :: datOut 
    integer :: ierr, i, j, k, idx
 
    call MPI_Init(ierr)
    call decomp_2d_init(2,2,2,0,0)
    call get_decomp_info(gp)
    x = 0.3d0
    y = 0.3d0
    z = 0.9d0

    idx = 1
    do k = 1,2 
        do j = 1,2
            do i = 1,2
                datIn(i,j,k) = real(idx,rkind)
                idx = idx + 1
            end do 
        end do 
    end do

    call getNearestNeighborValue(datIn,datOut,gp,dx,dy,dz,x,y,z)
    print*, "datOut = ", datOut

    call decomp_info_finalize(gp)
    call decomp_2d_finalize()
    call MPI_Finalize(ierr) 
end program