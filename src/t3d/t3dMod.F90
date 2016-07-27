module t3dMod
    use kind_parameters, only: rkind, mpirkind
    use exits,           only: GracefulExit
    
    use mpi

    implicit none
    private
    public :: t3d, square_factor, roundrobin_split, xnumbering
        
    logical :: xnumbering = .true.
    integer, dimension(0:2) :: perm = [0, 1, 2]

    type :: t3d
        private
        integer, public :: px, py, pz
        integer, public :: comm3d, commX, commY, commZ
        integer, public :: commXY, commYZ, commXZ
        integer, public :: rank3d, rankX, rankY, rankZ
        integer, public :: nprocs
        integer, dimension(0:2) :: coords3d
        integer, dimension(0:1) :: coordsX, coordsY, coordsZ
        integer, dimension(:,:), allocatable :: coordsXall, coordsYall, coordsZall
        integer :: px_y, px_z
        integer :: py_x, py_z
        integer :: pz_y, pz_x
        integer, dimension(3), public :: sz3d, st3d, en3d 
        integer, dimension(3), public :: szX , stX , enX  
        integer, dimension(3), public :: szY , stY , enY  
        integer, dimension(3), public :: szZ , stZ , enZ  
        integer, dimension(:,:), allocatable :: sz3DX, sz3DY, sz3DZ 
        integer, dimension(:,:), allocatable :: st3DX, st3DY, st3DZ 
        integer, dimension(:,:), allocatable :: en3DX, en3DY, en3DZ 
        integer, dimension(:,:), allocatable :: szXall , stXall , enXall
        integer, dimension(:,:), allocatable :: szYall , stYall , enYall
        integer, dimension(:,:), allocatable :: szZall , stZall , enZall
        integer, dimension(:), allocatable :: disp3DX, count3DX, dispX, countX
        integer, dimension(:), allocatable :: disp3DY, count3DY, dispY, countY
        integer, dimension(:), allocatable :: disp3DZ, count3DZ, dispZ, countZ
        logical :: unequalX = .true., unequalY = .true., unequalZ = .true.
        integer :: xleft, xright 
        integer :: yleft, yright 
        integer :: zleft, zright
        
        integer, dimension(:), allocatable :: splitx_y, splitx_z
        integer, dimension(:), allocatable :: splity_x, splity_z
        integer, dimension(:), allocatable :: splitz_y, splitz_x
        
    contains
        procedure :: transpose_3d_to_x
        procedure :: transpose_x_to_3d
        procedure :: transpose_3d_to_y
        procedure :: transpose_y_to_3d
        procedure :: transpose_3d_to_z
        procedure :: transpose_z_to_3d
        procedure :: print_summary
        procedure, private :: timed_transpose
        procedure :: time
        final :: destroy
    end type 
    
    interface t3d
        module procedure init, optimize_decomposition
    end interface

contains 

    function init(comm3d, nx, ny, nz, px, py, pz, periodic_, reorder, fail, createCrossCommunicators) result(this)
        use reductions, only: P_OR, P_MAXVAL
        type(t3d) :: this
        integer, intent(in) :: comm3d, nx, ny, nz, px, py, pz
        logical, dimension(3), intent(in) :: periodic_
        logical, intent(in) :: reorder
        logical, intent(out) :: fail
        logical, optional, intent(in) :: createCrossCommunicators
        logical, dimension(0:2) :: remain
        logical, dimension(0:2) :: periodic
        integer, dimension(0:2) :: procgrid
        integer :: ierr, dummycomm, proc
        integer, dimension(0:px-1) :: splitx
        integer, dimension(0:py-1) :: splity
        integer, dimension(0:pz-1) :: splitz
        ! real(rkind), dimension(20) :: times
        ! integer :: tind, i
        logical :: createCrossCommunicators_

        createCrossCommunicators_ = .true.
        if (present(createCrossCommunicators)) createCrossCommunicators_ = createCrossCommunicators

        ! if ((px == 0) .and. (py == 0) .and. (pz == 0)) then
        !     call optimize_decomposition(comm3d, nx, ny, nz, px, py, pz)
        ! end if 

        this%comm3d = comm3d
        ! tind = 0

        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        fail = .false.

        ! Safegaurd
        if ((px > nx) .or. (py > ny) .or. (pz > nz)) then
            fail = .true.
            return
        end if 
        if ((px == 0) .or. (py == 0) .or. (pz == 0)) then
            call GracefulExit("Cannot have zero processors in any direction", 0)
        end if 

        call mpi_comm_size(comm3d, this%nprocs, ierr)
        call mpi_comm_rank(comm3d, this%rank3d, ierr)
         
        if (px*py*pz .ne. this%nprocs) then
            call GracefulExit("(PX x PY x PZ) must equal the number of processes in &
            & comm3d sent in.", 1)
        end if

        this%szX(1) = nx; this%szY(2) = ny; this%szZ(3) = nz
        this%stX(1) =  1; this%stY(2) =  1; this%stZ(3) =  1
        this%enX(1) = nx; this%enY(2) = ny; this%enZ(3) = nz

        if ( xnumbering ) then
            perm = [2, 1, 0]
        else
            perm = [0, 1, 2]
        end if

        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        procgrid(perm) = [px, py, pz]
        periodic(perm) = periodic_(:)
        call mpi_cart_create(comm3d, 3, procgrid, periodic, reorder, this%comm3d, ierr)
        call mpi_cart_coords(this%comm3d, this%rank3d, 3, this%coords3d, ierr)

        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        ! Create X communicator and get properties
        remain(perm) = [.true., .false., .false.]
        call mpi_cart_sub(this%comm3d, remain, this%commX, ierr)
        call mpi_comm_rank(this%commX, this%rankX, ierr)
        call mpi_comm_size(this%commX, this%px, ierr)
        call mpi_cart_shift(this%commX, 0, 1, this%xleft, this%xright, ierr)

        ! Get points in X
        call roundrobin_split( nx, this%px, splitx )
        this%sz3d(1) = splitx( this%coords3D(perm(0)) )
        this%st3d(1) = sum( splitx(0:this%coords3D(perm(0))-1) ) + 1
        this%en3d(1) = sum( splitx(0:this%coords3D(perm(0))  ) )

        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        ! Create Y communicator and get properties
        remain(perm) = [.false., .true., .false.]
        call mpi_cart_sub(this%comm3d, remain, this%commY, ierr)
        call mpi_comm_rank(this%commY, this%rankY, ierr)
        call mpi_comm_size(this%commY, this%py, ierr)
        call mpi_cart_shift(this%commY, 0, 1, this%yleft, this%yright, ierr)

        ! Get points in Y
        call roundrobin_split( ny, this%py, splity )
        this%sz3d(2) = splity( this%coords3D(perm(1)) )
        this%st3d(2) = sum( splity(0:this%coords3D(perm(1))-1) ) + 1
        this%en3d(2) = sum( splity(0:this%coords3D(perm(1))  ) )

        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        ! Create Z communicator and get properties
        remain(perm) = [.false., .false., .true.]
        call mpi_cart_sub(this%comm3d, remain, this%commZ, ierr)
        call mpi_comm_rank(this%commZ, this%rankZ, ierr)
        call mpi_comm_size(this%commZ, this%pz, ierr)
        call mpi_cart_shift(this%commZ, 0, 1, this%zleft, this%zright, ierr)

        ! Get points in Z
        call roundrobin_split( nz, this%pz, splitz )
        this%sz3d(3) = splitz( this%coords3D(perm(2)) )
        this%st3d(3) = sum( splitz(0:this%coords3D(perm(2))-1) ) + 1
        this%en3d(3) = sum( splitz(0:this%coords3D(perm(2))  ) )

        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        if ( allocated(this%sz3DX) ) deallocate(this%sz3DX); allocate( this%sz3DX(3, 0:this%px-1) )
        if ( allocated(this%st3DX) ) deallocate(this%st3DX); allocate( this%st3DX(3, 0:this%px-1) )
        if ( allocated(this%en3DX) ) deallocate(this%en3DX); allocate( this%en3DX(3, 0:this%px-1) )
        call mpi_allgather(this%sz3D, 3, MPI_INTEGER, this%sz3DX, 3, MPI_INTEGER, this%commX, ierr)
        call mpi_allgather(this%st3D, 3, MPI_INTEGER, this%st3DX, 3, MPI_INTEGER, this%commX, ierr)
        call mpi_allgather(this%en3D, 3, MPI_INTEGER, this%en3DX, 3, MPI_INTEGER, this%commX, ierr)
        
        if ( allocated(this%sz3DY) ) deallocate(this%sz3DY); allocate( this%sz3DY(3, 0:this%py-1) )
        if ( allocated(this%st3DY) ) deallocate(this%st3DY); allocate( this%st3DY(3, 0:this%py-1) )
        if ( allocated(this%en3DY) ) deallocate(this%en3DY); allocate( this%en3DY(3, 0:this%py-1) )
        call mpi_allgather(this%sz3D, 3, MPI_INTEGER, this%sz3DY, 3, MPI_INTEGER, this%commY, ierr)
        call mpi_allgather(this%st3D, 3, MPI_INTEGER, this%st3DY, 3, MPI_INTEGER, this%commY, ierr)
        call mpi_allgather(this%en3D, 3, MPI_INTEGER, this%en3DY, 3, MPI_INTEGER, this%commY, ierr)
        
        if ( allocated(this%sz3DZ) ) deallocate(this%sz3DZ); allocate( this%sz3DZ(3, 0:this%pz-1) )
        if ( allocated(this%st3DZ) ) deallocate(this%st3DZ); allocate( this%st3DZ(3, 0:this%pz-1) )
        if ( allocated(this%en3DZ) ) deallocate(this%en3DZ); allocate( this%en3DZ(3, 0:this%pz-1) )
        call mpi_allgather(this%sz3D, 3, MPI_INTEGER, this%sz3DZ, 3, MPI_INTEGER, this%commZ, ierr)
        call mpi_allgather(this%st3D, 3, MPI_INTEGER, this%st3DZ, 3, MPI_INTEGER, this%commZ, ierr)
        call mpi_allgather(this%en3D, 3, MPI_INTEGER, this%en3DZ, 3, MPI_INTEGER, this%commZ, ierr)
        
        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        if (createCrossCommunicators_) then
            ! Create XY communicator
            remain(perm) = [.true., .true., .false.]
            call mpi_cart_sub(this%comm3d, remain, this%commXY, ierr)

            ! Create YZ communicator
            remain(perm) = [.false., .true., .true.]
            call mpi_cart_sub(this%comm3d, remain, this%commYZ, ierr)

            ! Create XZ communicator
            remain(perm) = [.true., .false., .true.]
            call mpi_cart_sub(this%comm3d, remain, this%commXZ, ierr)
        end if

        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        ! Now get local transpose pencil grid dimensions
        ! Optimized for vectorization, not thread performance

        ! X split
        fail = square_factor( this%px, this%sz3d(2), this%sz3d(3), this%px_y, this%px_z)
        if ( P_OR(fail, this%comm3d) ) then
            fail = .true.
            return
        end if
        call mpi_cart_create(this%commX, 2, [this%px_y, this%px_z], [.FALSE., .FALSE.], .FALSE., dummycomm, ierr)
        call mpi_cart_coords(dummycomm, this%rankX, 2, this%coordsX, ierr)

        if ( allocated(this%coordsXall) ) deallocate(this%coordsXall); allocate( this%coordsXall(0:1, 0:this%px-1) )
        call mpi_allgather(this%coordsX, 2, MPI_INTEGER, this%coordsXall, 2, MPI_INTEGER, this%commX, ierr)
        
        if ( allocated(this%szXall) ) deallocate(this%szXall); allocate( this%szXall(3, 0:this%px-1) )
        if ( allocated(this%stXall) ) deallocate(this%stXall); allocate( this%stXall(3, 0:this%px-1) )
        if ( allocated(this%enXall) ) deallocate(this%enXall); allocate( this%enXall(3, 0:this%px-1) )

        this%szXall(1,:) = nx; this%stXall(1,:) = 1; this%enXall(1,:) = nx

        if ( allocated(this%splitx_y) ) deallocate(this%splitx_y); allocate( this%splitx_y(0:this%px_y-1) )
        call roundrobin_split( this%sz3d(2), this%px_y, this%splitx_y)
        do proc = 0,this%px-1
            this%szXall(2,proc) = this%splitx_y( this%coordsXall(0,proc) )
            this%stXall(2,proc) = sum( this%splitx_y(0:this%coordsXall(0,proc)-1) ) + 1
            this%enXall(2,proc) = sum( this%splitx_y(0:this%coordsXall(0,proc)  ) )
        end do
        this%szX(2) = this%szXall( 2, this%rankX )
        this%stX(2) = this%stXall( 2, this%rankX )
        this%enX(2) = this%enXall( 2, this%rankX )

        if ( allocated(this%splitx_z) ) deallocate(this%splitx_z); allocate( this%splitx_z(0:this%px_z-1) )
        call roundrobin_split( this%sz3d(3), this%px_z, this%splitx_z)
        do proc = 0,this%px-1
            this%szXall(3,proc) = this%splitx_z( this%coordsXall(1,proc) )
            this%stXall(3,proc) = sum( this%splitx_z(0:this%coordsXall(1,proc)-1) ) + 1
            this%enXall(3,proc) = sum( this%splitx_z(0:this%coordsXall(1,proc)  ) )
        end do
        this%szX(3) = this%szXall( 3, this%rankX )
        this%stX(3) = this%stXall( 3, this%rankX )
        this%enX(3) = this%enXall( 3, this%rankX )

        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        ! Y split
        fail = square_factor(this%py, this%sz3d(1), this%sz3d(3),this%py_x, this%py_z)
        if ( P_OR(fail, this%comm3d) ) then
            fail = .true.
            return
        end if
        call mpi_cart_create(this%commY, 2, [this%py_x, this%py_z], [.FALSE., .FALSE.], .FALSE., dummycomm, ierr)
        call mpi_cart_coords(dummycomm, this%rankY, 2, this%coordsY, ierr)

        if ( allocated(this%coordsYall) ) deallocate(this%coordsYall); allocate( this%coordsYall(0:1, 0:this%py-1) )
        call mpi_allgather(this%coordsY, 2, MPI_INTEGER, this%coordsYall, 2, MPI_INTEGER, this%commY, ierr)
        
        if ( allocated(this%szYall) ) deallocate(this%szYall); allocate( this%szYall(3, 0:this%py-1) )
        if ( allocated(this%stYall) ) deallocate(this%stYall); allocate( this%stYall(3, 0:this%py-1) )
        if ( allocated(this%enYall) ) deallocate(this%enYall); allocate( this%enYall(3, 0:this%py-1) )

        this%szYall(2,:) = ny; this%stYall(2,:) = 1; this%enYall(2,:) = ny

        if ( allocated(this%splity_x) ) deallocate(this%splity_x); allocate( this%splity_x(0:this%py_x-1) )
        call roundrobin_split( this%sz3d(1), this%py_x, this%splity_x)
        do proc = 0,this%py-1
            this%szYall(1,proc) = this%splity_x( this%coordsYall(0,proc) )
            this%stYall(1,proc) = sum( this%splity_x(0:this%coordsYall(0,proc)-1) ) + 1
            this%enYall(1,proc) = sum( this%splity_x(0:this%coordsYall(0,proc)  ) )
        end do
        this%szY(1) = this%szYall( 1, this%rankY )
        this%stY(1) = this%stYall( 1, this%rankY )
        this%enY(1) = this%enYall( 1, this%rankY )

        if ( allocated(this%splity_z) ) deallocate(this%splity_z); allocate( this%splity_z(0:this%py_z-1) )
        call roundrobin_split( this%sz3d(3), this%py_z, this%splity_z)
        do proc = 0,this%py-1
            this%szYall(3,proc) = this%splity_z( this%coordsYall(1,proc) )
            this%stYall(3,proc) = sum( this%splity_z(0:this%coordsYall(1,proc)-1) ) + 1
            this%enYall(3,proc) = sum( this%splity_z(0:this%coordsYall(1,proc)  ) )
        end do
        this%szY(3) = this%szYall( 3, this%rankY )
        this%stY(3) = this%stYall( 3, this%rankY )
        this%enY(3) = this%enYall( 3, this%rankY )

        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        ! Z split
        fail = square_factor(this%pz, this%sz3d(1), this%sz3d(2),this%pz_x, this%pz_y) 
        if ( P_OR(fail, this%comm3d) ) then
            fail = .true.
            return
        end if
        call mpi_cart_create(this%commZ, 2, [this%pz_x, this%pz_y], [.FALSE., .FALSE.], .FALSE., dummycomm, ierr)
        call mpi_cart_coords(dummycomm, this%rankZ, 2, this%coordsZ, ierr)

        if ( allocated(this%coordsZall) ) deallocate(this%coordsZall); allocate( this%coordsZall(0:1, 0:this%pz-1) )
        call mpi_allgather(this%coordsZ, 2, MPI_INTEGER, this%coordsZall, 2, MPI_INTEGER, this%commZ, ierr)
        
        if ( allocated(this%szZall) ) deallocate(this%szZall); allocate( this%szZall(3, 0:this%pz-1) )
        if ( allocated(this%stZall) ) deallocate(this%stZall); allocate( this%stZall(3, 0:this%pz-1) )
        if ( allocated(this%enZall) ) deallocate(this%enZall); allocate( this%enZall(3, 0:this%pz-1) )

        this%szZall(3,:) = nz; this%stZall(3,:) = 1; this%enZall(3,:) = nz

        if ( allocated(this%splitz_x) ) deallocate(this%splitz_x); allocate( this%splitz_x(0:this%pz_x-1) )
        call roundrobin_split( this%sz3d(1), this%pz_x, this%splitz_x)
        do proc = 0,this%pz-1
            this%szZall(1,proc) = this%splitz_x( this%coordsZall(0,proc) )
            this%stZall(1,proc) = sum( this%splitz_x(0:this%coordsZall(0,proc)-1) ) + 1
            this%enZall(1,proc) = sum( this%splitz_x(0:this%coordsZall(0,proc)  ) )
        end do
        this%szZ(1) = this%szZall( 1, this%rankZ )
        this%stZ(1) = this%stZall( 1, this%rankZ )
        this%enZ(1) = this%enZall( 1, this%rankZ )

        if ( allocated(this%splitz_y) ) deallocate(this%splitz_y); allocate( this%splitz_y(0:this%pz_y-1) )
        call roundrobin_split( this%sz3d(2), this%pz_y, this%splitz_y)
        do proc = 0,this%pz-1
            this%szZall(2,proc) = this%splitz_y( this%coordsZall(1,proc) )
            this%stZall(2,proc) = sum( this%splitz_y(0:this%coordsZall(1,proc)-1) ) + 1
            this%enZall(2,proc) = sum( this%splitz_y(0:this%coordsZall(1,proc)  ) )
        end do
        this%szZ(2) = this%szZall( 2, this%rankZ )
        this%stZ(2) = this%stZall( 2, this%rankZ )
        this%enZ(2) = this%enZall( 2, this%rankZ )

        ! tind = tind + 1; times(tind) = this%time(barrier=.true.)

        ! do i = 2,tind
        !     times(i-1) = P_MAXVAL(times(i) - times(i-1))
        !     if (this%rank3d == 0) print*, "Time ", i-1, times(i-1)
        ! end do
        ! if (this%rank3d == 0) print*, "Total time ", sum(times(1:tind-1))
        
        ! Allocate X displacements and counts
        if ( allocated(this%disp3DX) ) deallocate(this%disp3DX); allocate( this%disp3DX(0:this%px-1) )
        if ( allocated(this%count3DX) ) deallocate(this%count3DX); allocate( this%count3DX(0:this%px-1) )

        if ( allocated(this%dispX) ) deallocate(this%dispX); allocate( this%dispX(0:this%px-1) )
        if ( allocated(this%countX) ) deallocate(this%countX); allocate( this%countX(0:this%px-1) )

        do proc = 0,this%px-1
            this%disp3DX(proc) = this%sz3D(1)*sum( this%szXall(2,0:proc-1) * this%szXall(3,0:proc-1) )
            this%count3DX(proc) = this%sz3D(1) * this%szXall(2,proc) * this%szXall(3,proc)
            
            this%dispX(proc) = sum( this%sz3DX(1,0:proc-1) * this%szX(2) * this%szX(3) )
            this%countX(proc) = this%sz3DX(1,proc) * this%szX(2) * this%szX(3)
        end do

        if ( (maxval(this%count3DX) == minval(this%count3DX)) .and. &
             (maxval(this%countX  ) == minval(this%countX  )) .and. &
             (maxval(this%count3DX) == minval(this%countX  )) ) this%unequalX = .false.

        ! Allocate Y displacements and counts
        if ( allocated(this%disp3DY) ) deallocate(this%disp3DY); allocate( this%disp3DY(0:this%py-1) )
        if ( allocated(this%count3DY) ) deallocate(this%count3DY); allocate( this%count3DY(0:this%py-1) )

        if ( allocated(this%dispY) ) deallocate(this%dispY); allocate( this%dispY(0:this%py-1) )
        if ( allocated(this%countY) ) deallocate(this%countY); allocate( this%countY(0:this%py-1) )

        do proc = 0,this%py-1
            this%disp3DY(proc) = this%sz3D(2)*sum( this%szYall(1,0:proc-1) * this%szYall(3,0:proc-1) )
            this%count3DY(proc) = this%sz3D(2) * this%szYall(1,proc) * this%szYall(3,proc)
            
            this%dispY(proc) = sum( this%sz3DY(2,0:proc-1) * this%szY(1) * this%szY(3) )
            this%countY(proc) = this%sz3DY(2,proc) * this%szY(1) * this%szY(3)
        end do

        if ( (maxval(this%count3DY) == minval(this%count3DY)) .and. &
             (maxval(this%countY  ) == minval(this%countY  )) .and. &
             (maxval(this%count3DY) == minval(this%countY  )) ) this%unequalY = .false.

        ! Allocate Z displacements and counts
        if ( allocated(this%disp3DZ) ) deallocate(this%disp3DZ); allocate( this%disp3DZ(0:this%pz-1) )
        if ( allocated(this%count3DZ) ) deallocate(this%count3DZ); allocate( this%count3DZ(0:this%pz-1) )

        if ( allocated(this%dispZ) ) deallocate(this%dispZ); allocate( this%dispZ(0:this%pz-1) )
        if ( allocated(this%countZ) ) deallocate(this%countZ); allocate( this%countZ(0:this%pz-1) )

        do proc = 0,this%pz-1
            this%disp3DZ(proc) = this%sz3D(3)*sum( this%szZall(1,0:proc-1) * this%szZall(2,0:proc-1) )
            this%count3DZ(proc) = this%sz3D(3) * this%szZall(1,proc) * this%szZall(2,proc)
            
            this%dispZ(proc) = sum( this%sz3DZ(3,0:proc-1) * this%szZ(1) * this%szZ(2) )
            this%countZ(proc) = this%sz3DZ(3,proc) * this%szZ(1) * this%szZ(2)
        end do

        if ( (maxval(this%count3DZ) == minval(this%count3DZ)) .and. &
             (maxval(this%countZ  ) == minval(this%countZ  )) .and. &
             (maxval(this%count3DZ) == minval(this%countZ  )) ) this%unequalZ = .false.

    end function
   

    pure elemental subroutine destroy(this)
        type(t3d), intent(inout) :: this

        if ( allocated(this%splitx_y) ) deallocate( this%splitx_y )
        if ( allocated(this%splitx_z) ) deallocate( this%splitx_z )

        if ( allocated(this%splity_x) ) deallocate( this%splity_x )
        if ( allocated(this%splity_z) ) deallocate( this%splity_z )

        if ( allocated(this%splitz_y) ) deallocate( this%splitz_y )
        if ( allocated(this%splitz_x) ) deallocate( this%splitz_x )
    end subroutine

    subroutine transpose_3d_to_x(this, input, output)
        class(t3d), intent(in) :: this
        real(rkind), dimension(this%sz3d(1),this%sz3d(2),this%sz3d(3)), intent(in)  :: input
        real(rkind), dimension(this%szX (1),this%szX (2),this%szX (3)), intent(out) :: output
        real(rkind), dimension(this%sz3d(1)*this%sz3d(2)*this%sz3d(3))              :: buffer3D
        real(rkind), dimension(this%szX (1)*this%szX (2)*this%szX (3))              :: bufferX
        integer :: proc, i, j, k, pos, ierr
        ! real(rkind) :: start, endt

        ! start = this%time(barrier=.false.)
        do proc = 0,this%px-1
            do k = this%stXall(3,proc),this%enXall(3,proc)
                do j = this%stXall(2,proc),this%enXall(2,proc)
                    do i = 1,this%sz3D(1)
                        pos = ( 1 + (i-1) + this%sz3D(1)*(j-this%stXall(2,proc)) + &
                              this%sz3D(1)*this%szXall(2,proc)*(k-this%stXall(3,proc)) ) + this%disp3DX(proc)
                        buffer3D(pos) = input(i,j,k)
                    end do
                end do
            end do
        end do
        ! endt = this%time(start,reduce=.false.)
        ! if (this%rank3d == 0) print*, "Do 1", endt

        ! start = this%time(barrier=.false.)
        select case(this%unequalX)
        case (.true.)
            call mpi_alltoallv(buffer3D,this%count3DX,this%disp3DX,mpirkind, &
                               bufferX, this%countX,  this%dispX,  mpirkind, this%commX, ierr)
        case (.false.)
            call mpi_alltoall (buffer3D,this%count3DX(0), mpirkind, &
                               bufferX, this%countX  (0), mpirkind, this%commX, ierr)
        end select
        ! endt = this%time(start,reduce=.false.)
        ! if (this%rank3d == 0) print*, "2", endt

        ! start = this%time(barrier=.false.)
        do proc = 0,this%px-1
            do k = this%stX(3),this%enX(3)
                do j = this%stX(2),this%enX(2)
                    do i = this%st3DX(1,proc),this%en3DX(1,proc)
                        pos = ( 1 + (i-this%st3DX(1,proc)) + &
                               this%sz3DX(1,proc)*(j-this%stX(2)) + & 
                               this%sz3DX(1,proc)*this%szX(2)*(k-this%stX(3)) ) + this%dispX(proc)
                        output(i,j-this%stX(2)+1,k-this%stX(3)+1) = bufferX(pos)
                    end do
                end do
            end do
        end do
        ! endt = this%time(start,reduce=.false.)
        ! if (this%rank3d == 0) print*, "Do 3", endt

    end subroutine 

    subroutine transpose_x_to_3d(this, input, output)
        class(t3d), intent(in) :: this
        real(rkind), dimension(this%szX (1),this%szX (2),this%szX (3)), intent(in)  :: input
        real(rkind), dimension(this%sz3d(1),this%sz3d(2),this%sz3d(3)), intent(out) :: output
        real(rkind), dimension(this%sz3d(1)*this%sz3d(2)*this%sz3d(3))              :: buffer3D
        real(rkind), dimension(this%szX (1)*this%szX (2)*this%szX (3))              :: bufferX
        integer :: proc, i, j, k, pos, ierr

        do proc = 0,this%px-1
            do k = this%stX(3),this%enX(3)
                do j = this%stX(2),this%enX(2)
                    do i = this%st3DX(1,proc),this%en3DX(1,proc)
                        pos = ( 1 + (i-this%st3DX(1,proc)) + this%sz3DX(1,proc)*(j-this%stX(2)) + &
                              this%sz3DX(1,proc)*this%szX(2)*(k-this%stX(3)) ) + this%dispX(proc)
                        bufferX(pos) = input(i,j-this%stX(2)+1,k-this%stX(3)+1)
                    end do
                end do
            end do
        end do

        select case(this%unequalX)
        case (.true.)
            call mpi_alltoallv(bufferX, this%countX,  this%dispX,  mpirkind, &
                               buffer3D,this%count3DX,this%disp3DX,mpirkind, this%commX, ierr)
        case (.false.)
            call mpi_alltoall (bufferX ,this%countX  (0), mpirkind, &
                               buffer3D,this%count3DX(0), mpirkind, this%commX, ierr)
        end select

        do proc = 0,this%px-1
            do k = this%stXall(3,proc),this%enXall(3,proc)
                do j = this%stXall(2,proc),this%enXall(2,proc)
                    do i = 1,this%sz3D(1)
                        pos = ( 1 + (i-1) + this%sz3D(1)*(j-this%stXall(2,proc)) + &
                              this%sz3D(1)*this%szXall(2,proc)*(k-this%stXall(3,proc)) ) + this%disp3DX(proc)
                        output(i,j,k) = buffer3D(pos)
                    end do
                end do
            end do
        end do

    end subroutine 

    subroutine transpose_3d_to_y(this, input, output)
        class(t3d), intent(in) :: this
        real(rkind), dimension(this%sz3d(1),this%sz3d(2),this%sz3d(3)), intent(in)  :: input
        real(rkind), dimension(this%szY (1),this%szY (2),this%szY (3)), intent(out) :: output
        real(rkind), dimension(this%sz3d(1)*this%sz3d(2)*this%sz3d(3))              :: buffer3D
        real(rkind), dimension(this%szY (1)*this%szY (2)*this%szY (3))              :: bufferY
        integer :: proc, i, j, k, pos, ierr
        ! real(rkind) :: start, endt

        ! start = this%time(barrier=.false.)
        do proc = 0,this%py-1
            do k = this%stYall(3,proc),this%enYall(3,proc)
                do j = 1,this%sz3D(2)
                    do i = this%stYall(1,proc),this%enYall(1,proc)
                        pos = ( 1 + (i-this%stYall(1,proc)) + this%szYall(1,proc)*(j-1) + &
                              this%szYall(1,proc)*this%sz3D(2)*(k-this%stYall(3,proc)) ) + this%disp3DY(proc)
                        buffer3D(pos) = input(i,j,k)
                    end do
                end do
            end do
        end do
        ! endt = this%time(start,reduce=.false.)
        ! if (this%rank3d == 0) print*, "Do 1", endt

        ! start = this%time(barrier=.false.)
        select case(this%unequalY)
        case (.true.)
            call mpi_alltoallv(buffer3D,this%count3DY,this%disp3DY,mpirkind, &
                               bufferY, this%countY,  this%dispY,  mpirkind, this%commY, ierr)
        case (.false.)
            call mpi_alltoall (buffer3D,this%count3DY(0), mpirkind, &
                               bufferY, this%countY  (0), mpirkind, this%commY, ierr)
        end select
        ! endt = this%time(start,reduce=.false.)
        ! if (this%rank3d == 0) print*, "2", endt

        ! start = this%time(barrier=.false.)
        do proc = 0,this%py-1
            do k = this%stY(3),this%enY(3)
                do j = this%st3DY(2,proc),this%en3DY(2,proc)
                    do i = this%stY(1),this%enY(1)
                        pos = ( 1 + (i-this%stY(1)) + &
                               this%szY(1)*(j-this%st3DY(2,proc)) + & 
                               this%szY(1)*this%sz3DY(2,proc)*(k-this%stY(3)) ) + this%dispY(proc)
                        output(i-this%stY(1)+1,j,k-this%stY(3)+1) = bufferY(pos)
                    end do
                end do
            end do
        end do
        ! endt = this%time(start,reduce=.false.)
        ! if (this%rank3d == 0) print*, "Do 3", endt

    end subroutine 

    subroutine transpose_y_to_3d(this, input, output)
        class(t3d), intent(in) :: this
        real(rkind), dimension(this%szY (1),this%szY (2),this%szY (3)), intent(in)  :: input
        real(rkind), dimension(this%sz3d(1),this%sz3d(2),this%sz3d(3)), intent(out) :: output
        real(rkind), dimension(this%sz3d(1)*this%sz3d(2)*this%sz3d(3))              :: buffer3D
        real(rkind), dimension(this%szY (1)*this%szY (2)*this%szY (3))              :: bufferY
        integer :: proc, i, j, k, pos, ierr

        do proc = 0,this%py-1
            do k = this%stY(3),this%enY(3)
                do j = this%st3DY(2,proc),this%en3DY(2,proc)
                    do i = this%stY(1),this%enY(1)
                        pos = ( 1 + (i-this%stY(1)) + &
                               this%szY(1)*(j-this%st3DY(2,proc)) + & 
                               this%szY(1)*this%sz3DY(2,proc)*(k-this%stY(3)) ) + this%dispY(proc)
                        bufferY(pos) = input(i-this%stY(1)+1,j,k-this%stY(3)+1)
                    end do
                end do
            end do
        end do

        select case(this%unequalY)
        case (.true.)
            call mpi_alltoallv(bufferY, this%countY,  this%dispY,  mpirkind, &
                               buffer3D,this%count3DY,this%disp3DY,mpirkind, this%commY, ierr)
        case (.false.)
            call mpi_alltoall (bufferY, this%countY  (0), mpirkind, &
                               buffer3D,this%count3DY(0), mpirkind, this%commY, ierr)
        end select

        do proc = 0,this%py-1
            do k = this%stYall(3,proc),this%enYall(3,proc)
                do j = 1,this%sz3D(2)
                    do i = this%stYall(1,proc),this%enYall(1,proc)
                        pos = ( 1 + (i-this%stYall(1,proc)) + this%szYall(1,proc)*(j-1) + &
                              this%szYall(1,proc)*this%sz3D(2)*(k-this%stYall(3,proc)) ) + this%disp3DY(proc)
                        output(i,j,k) = buffer3D(pos) 
                    end do
                end do
            end do
        end do

    end subroutine 

    subroutine transpose_3d_to_z(this, input, output)
        class(t3d), intent(in) :: this
        real(rkind), dimension(this%sz3d(1),this%sz3d(2),this%sz3d(3)), intent(in)  :: input
        real(rkind), dimension(this%szZ (1),this%szZ (2),this%szZ (3)), intent(out) :: output
        real(rkind), dimension(this%sz3d(1)*this%sz3d(2)*this%sz3d(3))              :: buffer3D
        real(rkind), dimension(this%szZ (1)*this%szZ (2)*this%szZ (3))              :: bufferZ
        integer :: proc, i, j, k, pos, ierr
        ! real(rkind) :: start, endt

        ! start = this%time(barrier=.false.)
        do proc = 0,this%pz-1
            do k = 1,this%sz3D(3)
                do j = this%stZall(2,proc),this%enZall(2,proc)
                    do i = this%stZall(1,proc),this%enZall(1,proc)
                        pos = ( 1 + (i-this%stZall(1,proc)) + this%szZall(1,proc)*(j-this%stZall(2,proc)) + &
                              this%szZall(1,proc)*this%szZall(2,proc)*(k-1) ) + this%disp3DZ(proc)
                        buffer3D(pos) = input(i,j,k)
                    end do
                end do
            end do
        end do
        ! endt = this%time(start,reduce=.false.)
        ! if (this%rank3d == 0) print*, "Do 1", endt

        ! start = this%time(barrier=.false.)
        select case(this%unequalZ)
        case (.true.)
            call mpi_alltoallv(buffer3D,this%count3DZ,this%disp3DZ,mpirkind, &
                               bufferZ, this%countZ,  this%dispZ,  mpirkind, this%commZ, ierr)
        case (.false.)
            call mpi_alltoall (buffer3D,this%count3DZ(0), mpirkind, &
                               bufferZ, this%countZ  (0), mpirkind, this%commZ, ierr)
        end select
        ! endt = this%time(start,reduce=.false.)
        ! if (this%rank3d == 0) print*, "2", endt

        ! start = this%time(barrier=.false.)
        do proc = 0,this%pz-1
            do k = this%st3DZ(3,proc),this%en3DZ(3,proc)
                do j = this%stZ(2),this%enZ(2)
                    do i = this%stZ(1),this%enZ(1)
                        pos = ( 1 + (i-this%stZ(1)) + &
                               this%szZ(1)*(j-this%stZ(2)) + & 
                               this%szZ(1)*this%szZ(2)*(k-this%st3DZ(3,proc)) ) + this%dispZ(proc)
                        output(i-this%stZ(1)+1,j-this%stZ(2)+1,k) = bufferZ(pos)
                    end do
                end do
            end do
        end do
        ! endt = this%time(start,reduce=.false.)
        ! if (this%rank3d == 0) print*, "Do 3", endt

    end subroutine 

    subroutine transpose_z_to_3d(this, input, output)
        class(t3d), intent(in) :: this
        real(rkind), dimension(this%szZ (1),this%szZ (2),this%szZ (3)), intent(in)  :: input
        real(rkind), dimension(this%sz3d(1),this%sz3d(2),this%sz3d(3)), intent(out) :: output
        real(rkind), dimension(this%sz3d(1)*this%sz3d(2)*this%sz3d(3))              :: buffer3D
        real(rkind), dimension(this%szZ (1)*this%szZ (2)*this%szZ (3))              :: bufferZ
        integer :: proc, i, j, k, pos, ierr

        do proc = 0,this%pz-1
            do k = this%st3DZ(3,proc),this%en3DZ(3,proc)
                do j = this%stZ(2),this%enZ(2)
                    do i = this%stZ(1),this%enZ(1)
                        pos = ( 1 + (i-this%stZ(1)) + &
                               this%szZ(1)*(j-this%stZ(2)) + & 
                               this%szZ(1)*this%szZ(2)*(k-this%st3DZ(3,proc)) ) + this%dispZ(proc)
                        bufferZ(pos) = input(i-this%stZ(1)+1,j-this%stZ(2)+1,k)
                    end do
                end do
            end do
        end do

        select case(this%unequalZ)
        case (.true.)
            call mpi_alltoallv(bufferZ, this%countZ,  this%dispZ,  mpirkind, &
                               buffer3D,this%count3DZ,this%disp3DZ,mpirkind, this%commZ, ierr)
        case (.false.)
            call mpi_alltoall (bufferZ, this%countZ  (0), mpirkind, &
                               buffer3D,this%count3DZ(0), mpirkind, this%commZ, ierr)
        end select

        do proc = 0,this%pz-1
            do k = 1,this%sz3D(3)
                do j = this%stZall(2,proc),this%enZall(2,proc)
                    do i = this%stZall(1,proc),this%enZall(1,proc)
                        pos = ( 1 + (i-this%stZall(1,proc)) + this%szZall(1,proc)*(j-this%stZall(2,proc)) + &
                              this%szZall(1,proc)*this%szZall(2,proc)*(k-1) ) + this%disp3DZ(proc)
                        output(i,j,k) = buffer3D(pos)
                    end do
                end do
            end do
        end do

    end subroutine 

    logical function square_factor(nprocs,nrow,ncol,prow,pcol) result(fail)
        use constants, only: eps
        integer, intent(in) :: nprocs, nrow, ncol
        integer, intent(out) :: prow, pcol
        integer :: incr

        fail = .false.
        if ( nrow < ncol ) then
            incr = -1
            prow = int(sqrt( real(nprocs,rkind) ) + 1000*eps)
            pcol = nprocs / prow
        else
            incr = 1
            pcol = int(sqrt( real(nprocs,rkind) ) + 1000*eps)
            prow = nprocs / pcol
        end if

        do while ( (prow > 1) .and. (prow < nprocs) .and. &
                 ( (mod(nprocs,prow) /= 0) .or. (prow > nrow) .or. (pcol > ncol) ) )
            prow = prow + incr
            pcol = nprocs / prow
        end do
        
        if ( (mod(nprocs,prow) /= 0) .or. (prow > nrow) .or. (pcol > ncol) ) then
            fail = .true.
        end if

    end function

    pure subroutine roundrobin_split(na, nb, split)
        integer, intent(in) :: na, nb
        integer, dimension(0:nb-1), intent(out) :: split
        integer :: i, idx

        split = 0
        do i = 0,na-1
            idx = mod(i,nb)
            split(idx) = split(idx) + 1
        end do
    end subroutine

    function optimize_decomposition(comm3d, nx, ny, nz, periodic) result(this)
        use constants,       only: rhuge
        use kind_parameters, only: stdout
        type(t3d) :: this
        integer, intent(in)  :: comm3d, nx, ny, nz
        logical, dimension(3), intent(in) :: periodic
        integer :: px, py, pz, pxopt, pyopt, pzopt
        real(rkind) :: t, topt
        integer :: pypz, nprocs, ierr, rank, niters, siters
        logical :: fail
        logical, parameter :: reorder = .false.

        call mpi_comm_size(comm3d,nprocs,ierr)
        call mpi_comm_rank(comm3d,rank  ,ierr)

        if (rank == 0) then
            write(stdout,'(A)') " "
            write(stdout,'(A)') " ================ Optimizing t3d ================ "
        end if

        niters = 0; siters = 0
        topt = rhuge
        do px = 1,nprocs
            if (mod(nprocs,px) /= 0) cycle
            pypz = nprocs / px
            do py = 1,pypz
                if (mod(pypz,py) /= 0) cycle
                pz = pypz / py
               
                niters = niters + 1
                call mpi_barrier(comm3d,ierr)
                this = t3d(comm3d, nx, ny, nz, px, py, pz, periodic, reorder, fail, createCrossCommunicators=.false.)
                if (.not. fail) then
                    siters = siters + 1
                    t = this%timed_transpose()
                    if (rank == 0) then
                        write(stdout,'(A,3(I5,A),ES12.3E3,A)') "Processor decomposition: ", &
                                        px, " x", py, " x", pz, ". Time = ", t, " seconds"
                    end if
                    if (t < topt) then
                        topt = t; pxopt = px; pyopt = py; pzopt = pz;
                        if (rank == 0) then
                            write(stdout,'(A,3(I0,A),ES12.3E3,A)') " >>>> Found a better processor decomposition ", &
                                            pxopt, " x ", pyopt, " x ", pzopt, " with time ", topt, " seconds"
                        end if
                    end if
                else
                    if (rank == 0) then
                        write(stdout,'(A,3(I3,A))') "Processor decomposition ", &
                                        px, " x ", py, " x ", pz, " infeasible."
                    end if
                end if
            end do
        end do

        this = t3d(comm3d, nx, ny, nz, pxopt, pyopt, pzopt, periodic, reorder, fail, createCrossCommunicators=.true.)
        if (fail) then
            print*, pxopt, pyopt, pzopt
            call GracefulExit("Couldn't find a working decomposition for t3d.",457)
        end if

        if (rank == 0) then
            print '(3(A,I0))', ">>>> Using 3D processor decomposition ", pxopt, 'x', pyopt, 'x', pzopt
            print '(1(A,I0))', ">>>> Total decompositions tried = ", niters
            print '(1(A,I0))', ">>>> Total feasible decompositions = ", siters
        end if

        if (rank == 0) then
            write(stdout,'(A)') " ================================================ "
            write(stdout,'(A)') " "
        end if

    end function

    function timed_transpose(gp) result(time)
        real(rkind) :: time
        class(t3d), intent(in) :: gp
        real(rkind), dimension(gp%sz3d(1),gp%sz3d(2),gp%sz3d(3)) :: array3D
        real(rkind), dimension(gp%szX (1),gp%szX (2),gp%szX (3)) :: arrayX 
        real(rkind), dimension(gp%szY (1),gp%szY (2),gp%szY (3)) :: arrayY 
        real(rkind), dimension(gp%szZ (1),gp%szZ (2),gp%szZ (3)) :: arrayZ 
        real(rkind) :: start
        
        array3D = real(gp%sz3d(1),rkind)

        start = gp%time(barrier=.true.,reduce=.false.)

        call gp%transpose_3D_to_x(array3D, arrayX)
        call gp%transpose_x_to_3d(arrayX, array3D)
        
        call gp%transpose_3D_to_y(array3D, arrayY)
        call gp%transpose_y_to_3d(arrayY, array3D)
        
        call gp%transpose_3D_to_z(array3D, arrayZ)
        call gp%transpose_z_to_3d(arrayZ, array3D)
        
        time = gp%time(start=start, barrier=.true., reduce=.true.)

    end function

    real(rkind) function time(this,start,barrier,reduce)
        class(t3d), intent(in) :: this
        real(rkind), optional, intent(in) :: start
        logical, optional, intent(in) :: barrier
        logical, optional, intent(in) :: reduce 
        logical :: barrier_, reduce_
        real(rkind) :: ptime
        integer :: ierr

        barrier_ = .false.
        if ( present(barrier) ) barrier_ = barrier

        reduce_ = .false.
        if ( present(reduce) ) reduce_ = reduce

        if (barrier_) call mpi_barrier(this%comm3d,ierr)

        ptime = mpi_wtime()
        if ( present(start) ) ptime = ptime - start

        select case (reduce_)
        case(.true.)
            call mpi_allreduce(ptime, time, 1, mpirkind, MPI_MAX, this%comm3d, ierr)
        case(.false.)
            time = ptime
        end select
    end function

    subroutine print_summary(this)
        use kind_parameters, only: stdout
        class(t3d), intent(in) :: this
        integer :: ierr

        if (this%rank3D == 0) then
            write(stdout,'(A)') "========== SUMMARY =========="
            write(stdout,'(3(A,I0))') "Grid size: ", this%szX(1), ' x ', this%szY(2), ' x ', this%szZ(3)

            ! Processors
            write(stdout,'(3(A,I0))') "Processor decomposition: ", this%px, ' x ', this%py, ' x ', this%pz
            write(stdout,'(A)') " "
        end if

        call mpi_barrier(this%comm3D, ierr)
        call sleep(this%rank3d)

        write(stdout,'(4(A,I0))') "Process ", this%rank3d, " Process coordinate: ", this%coords3D(perm(0)), ' x ', this%coords3D(perm(1)), ' x ', this%coords3D(perm(2))
        write(stdout,'(4(A,I0))') "Process ", this%rank3d, " Grid size: ", this%sz3d(1), ' x ', this%sz3d(2), ' x ', this%sz3d(3)
        write(stdout,'(4(A,I0))') "Process ", this%rank3d, " Grid start index: ", this%st3d(1), ' x ', this%st3d(2), ' x ', this%st3d(3)
        write(stdout,'(4(A,I0))') "Process ", this%rank3d, " Grid last index: ", this%en3d(1), ' x ', this%en3d(2), ' x ', this%en3d(3)
        write(stdout,'(A)') " "
        
        write(stdout,'(2(A,I0))') "X pencil processor grid: ", this%px_y, ' x ', this%px_z
        write(stdout,'(4(A,I0))') "X pencil size: ", this%szX(1), ' x ', this%szX(2), ' x ', this%szX(3)
        write(stdout,'(6(A,I0))') "X pencil start: ", this%stX(1), ':',this%enX(1), ' x ',  this%stX(2), ':', this%enX(2), ' x ', this%stX(3) , ':', this%enX(3)
        
        write(stdout,'(A)') " "
        write(stdout,'(2(A,I0))') "Y pencil processor grid: ", this%py_x, ' x ', this%py_z
        write(stdout,'(4(A,I0))') "Y pencil size: ", this%szY(1), ' x ', this%szY(2), ' x ', this%szY(3)
        write(stdout,'(6(A,I0))') "Y pencil start: ", this%stY(1), ':',this%enY(1), ' x ',  this%stY(2), ':', this%enY(2), ' x ', this%stY(3) , ':', this%enY(3)
        
        write(stdout,'(A)') " "
        write(stdout,'(2(A,I0))') "Z pencil processor grid: ", this%pz_x, ' x ', this%pz_y
        write(stdout,'(4(A,I0))') "Z pencil size: ", this%szZ(1), ' x ', this%szZ(2), ' x ', this%szZ(3)
        write(stdout,'(6(A,I0))') "Z pencil start: ", this%stZ(1), ':',this%enZ(1), ' x ',  this%stZ(2), ':', this%enZ(2), ' x ', this%stZ(3) , ':', this%enZ(3)
        write(stdout,'(A)') " -------------------------"
        write(stdout,'(A)') " "

        call sleep(1)
        call mpi_barrier(this%comm3D, ierr)
       
    end subroutine

end module 
