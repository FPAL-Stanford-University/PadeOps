
module kdtree_wrapper
  use kdtree2_module

  implicit none

  integer, private :: nx, ny, nz, nquery
  integer, private :: i, j, k, l, index
  integer, private, parameter :: Dim = 3, np = 6 
  integer, private :: npoints

  real(kdkind), dimension(:, :), allocatable, private :: grid
  real(kdkind), dimension(:), allocatable, private :: qvec 

  type(kdtree2), pointer, private :: tree, ibtree 
  type(kdtree2_result), dimension(:), allocatable, private :: results

  type lmesh_point_t
     real(kdkind), dimension(0:np) :: val
  end type lmesh_point_t

  type lmesh_inear_t
     integer :: inear, jnear, knear
  end type lmesh_inear_t

  type lmesh_rnear_t
     real(kdkind) :: xnear, ynear, znear
  end type lmesh_rnear_t

  public :: initialize_kdtree, &
       &    create_kdtree, &
       &    probe_nearest_points_kdtree, &
       &    finalize_kdtree 


contains 

  subroutine initialize_kdtree_ib(npts, x, y, z)
    implicit none
    integer, intent(in) :: npts
    real(kdkind), dimension(1:npts), intent(in) :: x, y, z
    integer :: i

    npoints = npts 

    allocate (grid(Dim, 1:npoints))
    allocate (qvec(Dim))
    allocate (results(1))

    do i = 1, npoints 
       grid(1, i) = x(i)
       grid(2, i) = y(i)
       grid(3, i) = z(i)
    end do

  end subroutine initialize_kdtree_ib

  subroutine initialize_kdtree(nptsx, nptsy, nptsz, nprobe, x, y, z)
    implicit none
    integer, intent(in) :: nptsx, nptsy, nptsz, nprobe
    real(kdkind), dimension(0:nptsx), intent(in) :: x
    real(kdkind), dimension(0:nptsy), intent(in) :: y
    real(kdkind), dimension(0:nptsz), intent(in) :: z 

    nx = nptsx
    ny = nptsy
    nz = nptsz

    nquery = nprobe

    npoints = (nx+1)*(ny+1)*(nz+1)
    ! Dim = 3

    allocate (grid(Dim, 1:npoints))
    allocate (qvec(Dim))
    allocate (results(1))

    do k = 0, nz
       do j = 0, ny
          do i = 0, nx 
             index = i + j*(nx+1) + k*(nx+1)*(ny+1) + 1
             grid(1, index) = x(i)
             grid(2, index) = y(j)
             grid(3, index) = z(k)
          end do
       end do
    end do

    write(*, *), 'Initialized KDTREE..'

  end subroutine initialize_kdtree

  subroutine create_kdtree_ib()
    implicit none

    ibtree =>kdtree2_create(grid, sort = .false., rearrange = .false.)

    write(*, *), 'Created kdtree for the immersed boundary'

  end subroutine create_kdtree_ib

  subroutine create_kdtree()
    implicit none
    
    tree => kdtree2_create(grid, sort = .false., rearrange = .false.)

    write(*, *), 'Created the kdtree...'

  end subroutine create_kdtree

  subroutine probe_nearest_points_kdtree_ib (imin, imax, jmin, jmax, kmin, kmax, &
       &                                     x, y, z, lnearlist)
    implicit none
    integer, intent(in) :: imin, imax, jmin, jmax, kmin, kmax 
    real(kdkind), dimension(imin:imax), intent(in) :: x 
    real(kdkind), dimension(jmin:jmax), intent(in) :: y
    real(kdkind), dimension(kmin:kmax), intent(in) :: z
    integer, dimension(imin:imax, jmin:jmax, kmin:kmax), intent(out) :: lnearlist
    integer :: i, j, k, index

    do k = kmin, kmax 
       do j = jmin, jmax 
          do i = imin, imax 

             qvec(1) = x(i)
             qvec(2) = y(j)
             qvec(3) = z(k)

             call kdtree2_n_nearest (tp = ibtree, qv = qvec, nn = 1, results = results)

             index = results(1)%idx 

             lnearlist(i, j, k) = index

          end do
       end do
    end do

  end subroutine probe_nearest_points_kdtree_ib

  subroutine probe_nearest_points_kdtree (nptsb, xpoints, ypoints, zpoints,  &
       &                                  inearlist, jnearlist, knearlist, &
       &                                  xnearlist, ynearlist, znearlist)
    implicit none
    integer, intent(in) :: nptsb
    real(kdkind), dimension(1:nptsb), intent(in) :: xpoints, ypoints, zpoints 
    real(kdkind), dimension(1:nptsb), intent(out) :: xnearlist, ynearlist, znearlist
    integer, dimension(1:nptsb), intent(out) :: inearlist, jnearlist, knearlist
    
    do l = 1, nptsb 

       qvec(1) = xpoints(l)
       qvec(2) = ypoints(l)
       qvec(3) = zpoints(l)

       call kdtree2_n_nearest (tp = tree, qv = qvec, nn = 1, results = results)

       index = results(1)%idx -1
       k = index/((nx+1)*(ny+1)) 
       j = (index - k*(nx+1)*(ny+1))/(nx+1)
       i = (index - k*(nx+1)*(ny+1) - j*(nx+1))

       inearlist(l) = i
       jnearlist(l) = j 
       knearlist(l) = k 

       xnearlist(l) = grid(1, index)
       ynearlist(l) = grid(2, index)
       znearlist(l) = grid(3, index)

    end do

  end subroutine probe_nearest_points_kdtree

  subroutine finalize_kdtree()
    implicit none

    call kdtree2_destroy(tree)

    deallocate (grid)
    deallocate (qvec)
    deallocate (results)

  end subroutine finalize_kdtree

  subroutine finalize_kdtree_ib()
    implicit none

    call kdtree2_destroy(ibtree)

    deallocate (grid)
    deallocate (qvec)
    deallocate (results)

  end subroutine finalize_kdtree_ib

end module kdtree_wrapper
