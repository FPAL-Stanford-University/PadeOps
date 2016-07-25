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
        integer :: xleft, xright 
        integer :: yleft, yright 
        integer :: zleft, zright
        
        integer, dimension(:), allocatable :: splitx_y, splitx_z
        integer, dimension(:), allocatable :: splity_x, splity_z
        integer, dimension(:), allocatable :: splitz_y, splitz_x
        
    contains
        procedure :: transpose_3d_to_x
        procedure :: transpose_x_to_3d
        procedure :: print_summary
        final :: destroy
    end type 
    
    interface t3d
        module procedure init
    end interface

contains 

    function init(comm3d, nx, ny, nz, px, py, pz, periodic_, reorder, fail) result(this)
        type(t3d) :: this
        integer, intent(in) :: comm3d, nx, ny, nz, px, py, pz
        logical, dimension(3), intent(in) :: periodic_
        logical, intent(in) :: reorder
        logical, intent(out) :: fail
        logical, dimension(0:2) :: remain
        logical, dimension(0:2) :: periodic
        integer, dimension(0:2) :: procgrid
        integer :: ierr, dummycomm, proc
        integer, dimension(0:px-1) :: splitx
        integer, dimension(0:py-1) :: splity
        integer, dimension(0:pz-1) :: splitz

        ! if ((px == 0) .and. (py == 0) .and. (pz == 0)) then
        !     call optimize_decomposition(comm3d, nx, ny, nz, px, py, pz)
        ! end if 

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

        procgrid(perm) = [px, py, pz]
        periodic(perm) = periodic_(:)
        call mpi_cart_create(comm3d, 3, procgrid, periodic, reorder, this%comm3d, ierr)
        call mpi_cart_coords(this%comm3d, this%rank3d, 3, this%coords3d, ierr)

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
        
        ! Create XY communicator
        remain(perm) = [.true., .true., .false.]
        call mpi_cart_sub(this%comm3d, remain, this%commXY, ierr)

        ! Create YZ communicator
        remain(perm) = [.false., .true., .true.]
        call mpi_cart_sub(this%comm3d, remain, this%commYZ, ierr)

        ! Create XZ communicator
        remain(perm) = [.true., .false., .true.]
        call mpi_cart_sub(this%comm3d, remain, this%commXZ, ierr)

        ! Now get local transpose pencil grid dimensions
        ! Optimized for vectorization, not thread performance

        ! X split
        if ( square_factor( this%px, this%sz3d(2), this%sz3d(3), this%px_y, this%px_z) ) then
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

        ! Y split
        if ( square_factor(this%py, this%sz3d(1), this%sz3d(3),this%py_x, this%py_z) ) then
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

        ! Z split
        if ( square_factor(this%pz, this%sz3d(1), this%sz3d(2),this%pz_x, this%pz_y) ) then
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

    subroutine transpose_3d_to_x(this, input, output)
        class(t3d), intent(in) :: this
        real(rkind), dimension(this%sz3d(1),this%sz3d(2),this%sz3d(3)), intent(inout)  :: input
        real(rkind), dimension(this%szX (1),this%szX (2),this%szX (3)), intent(out) :: output
        real(rkind), dimension(this%sz3d(1)*this%sz3d(2)*this%sz3d(3))              :: buffer3D
        real(rkind), dimension(this%szX (1)*this%szX (2)*this%szX (3))              :: bufferX
        integer :: proc, i, j, k, pos, ierr
        integer, dimension(0:this%px-1) :: disp3D, count3D
        integer, dimension(0:this%px-1) :: dispX, countX

        do proc = 0,this%px-1
            input(:,this%stXall(2,proc):this%enXall(2,proc),this%stXall(3,proc):this%enXall(3,proc)) = real(this%stXall(2,proc)*this%stXall(3,proc)*this%rank3D, rkind)
            
            disp3D(proc) = this%sz3D(1)*sum( this%szXall(2,0:proc-1) * this%szXall(3,0:proc-1) )
            count3D(proc) = this%sz3D(1) * this%szXall(2,proc) * this%szXall(3,proc)
            
            dispX(proc) = sum( this%sz3DX(1,0:proc-1) * this%szX(2) * this%szX(3) )
            countX(proc) = this%sz3DX(1,proc) * this%szX(2) * this%szX(3)
        end do

        do proc = 0,this%px-1
            do k = this%stXall(3,proc),this%enXall(3,proc)
                do j = this%stXall(2,proc),this%enXall(2,proc)
                    do i = 1,this%sz3D(1)
                        pos = ( 1 + (i-1) + this%sz3D(1)*(j-this%stXall(2,proc)) + this%sz3D(1)*this%szXall(2,proc)*(k-this%stXall(3,proc)) ) + disp3D(proc)
                        buffer3D(pos) = input(i,j,k)
                    end do
                end do
            end do
        end do

        call mpi_alltoallv(buffer3D,count3D,disp3D,mpirkind,bufferX,countX,dispX,mpirkind,this%commX,ierr)

        do proc = 0,this%px-1
            do k = this%stX(3),this%enX(3)
                do j = this%stX(2),this%enX(2)
                    do i = this%st3DX(1,proc),this%en3DX(1,proc)
                        pos = ( 1 + (i-this%st3DX(1,proc)) + this%sz3DX(1,proc)*(j-this%stX(2)) + this%sz3DX(1,proc)*this%szX(2)*(k-this%stX(3)) ) + dispX(proc)
                        output(i,j-this%stX(2)+1,k-this%stX(3)+1) = bufferX(pos)
                    end do
                end do
            end do
        end do

        ! call sleep(this%rank3D)
        ! print*, "---------------------------------------"
        ! print*, "rank = ", this%rank3D, this%rankX
        ! print*, count3D, sum(count3D)
        ! print*, disp3D
        ! print*, countX, sum(countX)
        ! print*, dispX
        ! print*, output
        ! print*, "---------------------------------------"

    end subroutine 

    subroutine transpose_x_to_3d(this, input, output)
        class(t3d), intent(in) :: this
        real(rkind), dimension(this%szX (1),this%szX (2),this%szX (3)), intent(in)  :: input
        real(rkind), dimension(this%sz3d(1),this%sz3d(2),this%sz3d(3)), intent(out) :: output
        real(rkind), dimension(this%sz3d(1)*this%sz3d(2)*this%sz3d(3))              :: buffer3D
        real(rkind), dimension(this%szX (1)*this%szX (2)*this%szX (3))              :: bufferX
        integer :: proc, i, j, k, pos, ierr
        integer, dimension(0:this%px-1) :: disp3D, count3D
        integer, dimension(0:this%px-1) :: dispX, countX

        do proc = 0,this%px-1
            disp3D(proc) = this%sz3D(1)*sum( this%szXall(2,0:proc-1) * this%szXall(3,0:proc-1) )
            count3D(proc) = this%sz3D(1) * this%szXall(2,proc) * this%szXall(3,proc)
            
            dispX(proc) = sum( this%sz3DX(1,0:proc-1) * this%szX(2) * this%szX(3) )
            countX(proc) = this%sz3DX(1,proc) * this%szX(2) * this%szX(3)
        end do

        do proc = 0,this%px-1
            do k = this%stX(3),this%enX(3)
                do j = this%stX(2),this%enX(2)
                    do i = this%st3DX(1,proc),this%en3DX(1,proc)
                        pos = ( 1 + (i-this%st3DX(1,proc)) + this%sz3DX(1,proc)*(j-this%stX(2)) + this%sz3DX(1,proc)*this%szX(2)*(k-this%stX(3)) ) + dispX(proc)
                        bufferX(pos) = input(i,j-this%stX(2)+1,k-this%stX(3)+1)
                    end do
                end do
            end do
        end do

        call mpi_alltoallv(bufferX,countX,dispX,mpirkind,buffer3D,count3D,disp3D,mpirkind,this%commX,ierr)

        do proc = 0,this%px-1
            do k = this%stXall(3,proc),this%enXall(3,proc)
                do j = this%stXall(2,proc),this%enXall(2,proc)
                    do i = 1,this%sz3D(1)
                        pos = ( 1 + (i-1) + this%sz3D(1)*(j-this%stXall(2,proc)) + this%sz3D(1)*this%szXall(2,proc)*(k-this%stXall(3,proc)) ) + disp3D(proc)
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
        else
            incr = 1
        end if

        prow = int(sqrt( real(nprocs,rkind) ) + 1000*eps)
        pcol = nprocs / prow
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

    subroutine optimize_decomposition(comm3d, nx, ny, nz, px, py, pz)
        integer, intent(in)  :: comm3d, nx, ny, nz
        integer, intent(out) :: px, py, pz
        integer :: pypz, nprocs, ierr, rank

        call mpi_comm_size(comm3d,nprocs,ierr)
        call mpi_comm_rank(comm3d,rank  ,ierr)

        px = int( ( real(nprocs,rkind) )**(1._rkind/3._rkind) )
        do while (MOD(nprocs,px) /= 0)
            px = px - 1
        end do
        pypz = nprocs / px

        py = int( sqrt( real(pypz,rkind) ) )
        do while (MOD(pypz,py) /= 0)
            py = py - 1
        end do
        pz = pypz / py

        if (rank == 0) then
            print '(3(A,I0))', "Using 3D processor decomposition ", px, 'x', py, 'x', pz
        end if

    end subroutine

end module 
