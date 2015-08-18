subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: two,pi
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,2*pi)x[0,2*pi)x[0,2*pi) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = two*pi/real(nx,rkind)
        dy = two*pi/real(ny,rkind)
        dz = two*pi/real(nz,rkind)

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( proc_st(1) + i - 1, rkind ) * dx
                    y(i,j,k) = real( proc_st(2) + j - 1, rkind ) * dy
                    z(i,j,k) = real( proc_st(3) + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inpDirectory,mesh,fields)
    use kind_parameters,  only: rkind
    use constants,        only: zero,one,two,pi
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index
    use decomp_2d,          only: decomp_info
    
    implicit none
    character(len=*),                                               intent(in)    :: inpDirectory
    type(decomp_info),                                              intent(in)    :: decomp
    real(rkind),                                                    intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:),     intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    associate( rho => fields(:,:,:,rho_index), u => fields(:,:,:,u_index), &
                 v => fields(:,:,:,  v_index), w => fields(:,:,:,w_index), &
                 p => fields(:,:,:,  p_index), T => fields(:,:,:,T_index), &
                 e => fields(:,:,:,  e_index),                             &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        rho = one
        u   =  sin(x)*cos(y)*cos(z)
        v   = -cos(x)*sin(y)*cos(z)
        w   = zero
        p   = 100._rkind + ( (cos(two*z) + two)*(cos(two*x) + cos(two*y)) - two ) / 16._rkind
           
    end associate

end subroutine


program taylorgreen

    use kind_parameters,  only: rkind,clen,stdout
    use CompressibleGrid, only: cgrid,u_index
    use GridMod,          only: alloc_buffs, destroy_buffs
    use reductions,       only: P_MAXVAL
    use exits,            only: message
    use timer,            only: tic, toc
    implicit none

    type(cgrid) :: cgp
    character(len=clen) :: inputfile, iters
    integer :: ierr, niters, idx 

    real(rkind), dimension(:,:,:,:), allocatable :: grad_u

    ! Start MPI
    call MPI_Init(ierr)

    ! Get file location 
    call GETARG(1,inputfile)
    
    ! Initialize the grid object
    call cgp%init(inputfile)

    
    ! Destroy everythin before ending
    call cgp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
