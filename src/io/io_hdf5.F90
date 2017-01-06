module IO_HDF5_stuff
    use hdf5
    use mpi
    use kind_parameters, only: rkind, clen
    use decomp_2d,       only: decomp_info, nrank, nproc
    use exits,           only: GracefulExit
    implicit none

    type :: io_hdf5

        integer :: vizcount
        character(len=clen) :: filename
        character(len=clen) :: basename
        character(len=clen) :: vizdir
        character(len=clen) :: xdmf_filename

        integer :: comm                  ! Communicator for parallel I/O
        character(len=1) :: pencil       ! Which pencil to use for I/O

        integer(hid_t) :: file_id        ! File identifier
        integer        :: xdmf_file_id   ! XDMF file identifier
        integer(hid_t) :: dset_id        ! Dataset identifier
        integer(hid_t) :: filespace      ! Dataspace identifier in file
        integer(hid_t) :: memspace       ! Dataspace identifier in memory
        integer(hid_t) :: plist_id       ! Property list identifier

        integer :: rank = 3              ! Dataset rank

        integer(hsize_t), dimension(3) :: dimsf      ! Dataset dimensions in file
        integer(hsize_t), dimension(3) :: chunk_dims ! Chunk dimensions

        integer(hsize_t),  dimension(3) :: count
        integer(hssize_t), dimension(3) :: offset
        integer(hsize_t),  dimension(3) :: stride
        integer(hsize_t),  dimension(3) :: block

        logical :: master

    contains

        procedure :: init
        procedure :: write_dataset
        procedure :: write_coords
        procedure :: start_viz
        procedure :: end_viz
        procedure :: write_variable
        final     :: destroy

    end type

contains

    subroutine init(this, comm_, gp, pencil_, vizdir_, filename_)
        class(io_hdf5),     intent(inout) :: this
        integer,            intent(in)    :: comm_
        class(decomp_info), intent(in)    :: gp
        character(len=1),   intent(in)    :: pencil_
        character(len=*),   intent(in)    :: vizdir_
        character(len=*),   intent(in)    :: filename_

        integer :: nrank
        integer :: info
        integer :: error

        this%comm = comm_
        this%pencil = pencil_

        this%vizcount = 0
        this%vizdir = adjustl(trim(vizdir_))

        if (trim(filename_) == '') call GracefulExit("Cannot have empty string for HDF5 output filename. &
                                                     &Re-run with a non-empty string", 7356)
        this%filename = adjustl(trim(vizdir_)) // "/" // adjustl(trim(filename_)) // '.h5'
        this%basename = adjustl(trim(filename_)) // '.h5'

        call mpi_comm_rank(this%comm, nrank, error)
        if (nrank == 0) then
            this%master = .true.
        else
            this%master = .false.
        end if

        ! Create vizdir if it does not exist
        call system('mkdir -p ' // adjustl(trim(this%vizdir)))

        ! Initialize the HDF5 library and Fortran interfaces
        call h5open_f(error)
        if (error /= 0) call GracefulExit("Could not initialize HDF5 library and Fortran interfaces.",7356)

        ! Setup file access property list with parallel I/O access.
        call h5pcreate_f(H5P_FILE_ACCESS_F, this%plist_id, error)
        call h5pset_fapl_mpio_f(this%plist_id, this%comm, info, error)

        ! Create the file collectively
        call h5fcreate_f(this%filename, H5F_ACC_TRUNC_F, this%file_id, error, access_prp = this%plist_id)
        call h5pclose_f(this%plist_id, error)

        ! Set the dataset dimensions
        this%dimsf = [gp%xsz(1), gp%ysz(2), gp%zsz(3)]

        ! Set dimensions of each chunk
        select case(this%pencil)
        case('x')
            this%chunk_dims = gp%xsz
            this%offset     = gp%xst - 1
        case('y')
            this%chunk_dims = gp%ysz
            this%offset     = gp%yst - 1
        case('z')
            this%chunk_dims = gp%zsz
            this%offset     = gp%zst - 1
        case default
            call GracefulExit("Pencil variable for HDF5 I/O should be 'x', 'y' or 'z'",7356)
        end select

        ! Set stride, count, block and offset variables (assuming contiguous chunk of memory. Have the flexibility to change this
        ! later)
        this%stride = 1
        this%count  = 1
        this%block  = this%chunk_dims
    end subroutine

    impure elemental subroutine destroy(this)
        type(io_hdf5), intent(inout) :: this

        integer :: error

        this%vizcount = 0
        this%filename = ''
        this%vizdir = ''

        ! Close the file
        call h5fclose_f(this%file_id, error)

        ! Close Fortran interfaces and HDF5 library
        call h5close_f(error)
    end subroutine

    subroutine write_dataset(this,field,dsetname)
        class(io_hdf5), intent(inout) :: this
        real(rkind), dimension(this%chunk_dims(1),this%chunk_dims(2),this%chunk_dims(3)), intent(in) :: field
        character(len=*), intent(in) :: dsetname

        integer :: error

        ! Create the dataspace for the dataset
        call h5screate_simple_f(this%rank, this%dimsf, this%filespace, error)
        call h5screate_simple_f(this%rank, this%chunk_dims, this%memspace, error)

        ! Create chunked dataset
        call h5pcreate_f(H5P_DATASET_CREATE_F, this%plist_id, error)
        call h5pset_chunk_f(this%plist_id, this%rank, this%chunk_dims, error)
        call h5dcreate_f(this%file_id, adjustl(trim(dsetname)), H5T_NATIVE_DOUBLE, this%filespace, &
                         this%dset_id, error, this%plist_id)
        call h5sclose_f(this%filespace, error)

        ! Select hyperslab in file
        call h5dget_space_f(this%dset_id, this%filespace, error)
        call h5sselect_hyperslab_f(this%filespace, H5S_SELECT_SET_F, this%offset, this%count, error, &
                                   this%stride, this%block)

        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, this%plist_id, error)
        call h5pset_dxpl_mpio_f(this%plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! Write dataset collectively
        call h5dwrite_f(this%dset_id, H5T_NATIVE_DOUBLE, field, this%dimsf, error, &
                        file_space_id = this%filespace, mem_space_id = this%memspace, xfer_prp = this%plist_id)

        ! Close dataspaces
        call h5sclose_f(this%filespace, error)
        call h5sclose_f(this%memspace, error)

        ! Close the dataset
        call h5dclose_f(this%dset_id, error)

        ! Close the property list
        call h5pclose_f(this%plist_id, error)

    end subroutine

    subroutine write_coords(this, coords)
        class(io_hdf5), intent(inout) :: this
        real(rkind), dimension(this%chunk_dims(1),this%chunk_dims(2),this%chunk_dims(3),3), intent(in) :: coords

        call this%write_dataset(coords(:,:,:,1), 'coords/X')
        call this%write_dataset(coords(:,:,:,2), 'coords/Y')
        call this%write_dataset(coords(:,:,:,3), 'coords/Z')

    end subroutine

    subroutine start_viz(this, time)
        class(io_hdf5), intent(inout) :: this
        real(rkind),    intent(in)    :: time

        this%xdmf_filename = ''
        write(this%xdmf_filename,'(A,A,I4.4,A)') adjustl(trim(this%vizdir)), "/", this%vizcount, '.xmf'

        if (this%master) then
            open(unit=this%xdmf_file_id, file=adjustl(trim(this%xdmf_filename)), form='FORMATTED', status='REPLACE')

            write(this%xdmf_file_id,'(A)')           '<?xml version="1.0" ?>'
            write(this%xdmf_file_id,'(A)')           '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
            write(this%xdmf_file_id,'(A)')           '<Xdmf Version="2.0">'
            write(this%xdmf_file_id,'(A,ES26.16,A)') '<Time Value="', time, '" />'
            write(this%xdmf_file_id,'(A)')           ' <Domain>'
            write(this%xdmf_file_id,'(A)')           '   <Grid Name="mesh" GridType="Uniform">'
            write(this%xdmf_file_id,'(A,3(I0,A))')   '     <Topology TopologyType="3DSMesh" NumberOfElements="', this%dimsf(1), ' ', this%dimsf(2),' ', this%dimsf(3), '"/>'
            write(this%xdmf_file_id,'(A)')           '     <Geometry GeometryType="X_Y_Z">'
            write(this%xdmf_file_id,'(A,3(I0,A))')   '       <DataItem Dimensions="', this%dimsf(1), ' ', this%dimsf(2),' ', this%dimsf(3), '" NumberType="Float" Precision="8" Format="HDF">'
            write(this%xdmf_file_id,'(3A)')          '        ', adjustl(trim(this%basename)), ':/coords/X'
            write(this%xdmf_file_id,'(A)')           '       </DataItem>'
            write(this%xdmf_file_id,'(A,3(I0,A))')   '       <DataItem Dimensions="', this%dimsf(1), ' ', this%dimsf(2),' ', this%dimsf(3), '" NumberType="Float" Precision="8" Format="HDF">'
            write(this%xdmf_file_id,'(3A)')          '        ', adjustl(trim(this%basename)), ':/coords/Y'
            write(this%xdmf_file_id,'(A)')           '       </DataItem>'
            write(this%xdmf_file_id,'(A,3(I0,A))')   '       <DataItem Dimensions="', this%dimsf(1), ' ', this%dimsf(2),' ', this%dimsf(3), '" NumberType="Float" Precision="8" Format="HDF">'
            write(this%xdmf_file_id,'(3A)')          '        ', adjustl(trim(this%basename)), ':/coords/Z'
            write(this%xdmf_file_id,'(A)')           '       </DataItem>'
            write(this%xdmf_file_id,'(A)')           '     </Geometry>'
        end if

    end subroutine

    subroutine write_variable(this, field, varname)
        class(io_hdf5),                                                                   intent(inout) :: this
        real(rkind), dimension(this%chunk_dims(1),this%chunk_dims(2),this%chunk_dims(3)), intent(in)    :: field
        character(len=*),                                                                 intent(in)    :: varname

        character(len=clen) :: dsetname

        write(dsetname,'(I4.4,A,A)') this%vizcount, '/', adjustl(trim(varname))
        call this%write_dataset(field, dsetname)

        if (this%master) then
            write(this%xdmf_file_id,'(A)')           '     <Attribute Name="', adjustl(trim(varname)), '" AttributeType="Scalar" Center="Node">'
            write(this%xdmf_file_id,'(A,3(I0,A))')   '       <DataItem Dimensions="', this%dimsf(1), ' ', this%dimsf(2),' ', this%dimsf(3), '" NumberType="Float" Precision="8" Format="HDF">'
            write(this%xdmf_file_id,'(4A)')          '        ', adjustl(trim(this%basename)), ':/', adjustl(trim(dsetname))
            write(this%xdmf_file_id,'(A)')           '       </DataItem>'
            write(this%xdmf_file_id,'(A)')           '     </Attribute>'
        end if

    end subroutine

    subroutine end_viz(this)
        class(io_hdf5), intent(inout) :: this

        if (this%master) then
            write(this%xdmf_file_id,'(A)')           '   </Grid>'
            write(this%xdmf_file_id,'(A)')           ' </Domain>'
            write(this%xdmf_file_id,'(A)')           '</Xdmf>'

            close(this%xdmf_file_id)
        end if

        this%vizcount = this%vizcount + 1
    end subroutine

end module
