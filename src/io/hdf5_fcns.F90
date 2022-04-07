module hdf5_fcns
  !-----------------------------------------------------------------------------!
  ! The below code was modifeid from the example at the HDF Group found at
  ! https://raw.githubusercontent.com/HDFGroup/hdf5/develop/fortran/examples/
  ! ph5example.f90
  !-----------------------------------------------------------------------------!
  use kind_parameters, only: clen, rkind
  use mpi
  use hdf5
  use fortran_assert, only: assert
  implicit none
  
  interface read_h5_chunk_data
    module procedure read_h5_chunk_data2, read_h5_chunk_data3, &
                     read_h5_chunk_data4, read_h5_chunk_data5
  end interface
  interface write_h5_chunk_data
    module procedure write_h5_chunk_data2, write_h5_chunk_data3, &
                     write_h5_chunk_data4
  end interface
contains
  
  subroutine read_h5_chunk_data2(block, chunk_dims, count, dimsf, dsetname, &
                                dset_rank, fname, memtype, offset, stride, &
                                d_read, communicator)
    !-------------------------------------------------------------------------!
    ! Read HDF5 data in parallel. This will read a 2 dimensional dataset      !
    !-------------------------------------------------------------------------!
    character(len=*), intent(in) :: fname
    real(rkind), dimension(:,:), intent(inout) :: d_read
    integer(HSIZE_T), dimension(:), intent(in) :: dimsf, chunk_dims, count, &
                                                  offset, stride, block
    integer(HID_T), intent(in) :: memtype
    integer, intent(in) :: dset_rank
    character(len=*), intent(in) :: dsetname
    integer, intent(in), optional :: communicator
    integer :: cart_comm
    integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace 
    integer :: error

    cart_comm = MPI_COMM_WORLD
    if (present(communicator)) cart_comm = communicator

    ! Setup file access property list with parllel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error) ! 1
    call h5pset_fapl_mpio_f(plist_id, cart_comm, MPI_INFO_NULL, error) ! 2

    ! Create the file collectively
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id) ! 3
    call h5pclose_f(plist_id, error) ! 4

    ! Create the data space for the dataset
    call h5screate_simple_f(dset_rank, dimsf, dataspace, error) ! 5
    call h5screate_simple_f(dset_rank, dimsf, filespace, error) ! 7
    call h5screate_simple_f(dset_rank, chunk_dims, memspace, error) ! 6

    ! Select hyperslab in the file
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, &
                               error, stride, block) ! 8

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    ! Read the dataset collectively
    call h5dopen_f(file_id, trim(dsetname), dset_id, error, H5P_DEFAULT_F) ! 9
    call h5dread_f(dset_id, memtype, d_read, dimsf, error, &
                   mem_space_id = memspace, file_space_id = filespace,&
                   xfer_prp = plist_id) ! 10

    ! Close resources
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dataspace, error)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)    
  end subroutine
  
  subroutine read_h5_chunk_data3(block, chunk_dims, count, dimsf, dsetname, &
                                dset_rank, fname, memtype, offset, stride, &
                                d_read, communicator)
    !-------------------------------------------------------------------------!
    ! Read HDF5 data in parallel. This will read a 3 dimensional dataset      !
    !-------------------------------------------------------------------------!
    character(len=*), intent(in) :: fname
    real(rkind), dimension(:,:,:), intent(inout) :: d_read
    integer(HSIZE_T), dimension(:), intent(in) :: dimsf, chunk_dims, count, &
                                                  offset, stride, block
    integer(HID_T), intent(in) :: memtype
    integer, intent(in) :: dset_rank
    character(len=*), intent(in) :: dsetname
    integer, intent(in), optional :: communicator
    integer :: cart_comm
    integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace 
    integer :: error

    cart_comm = MPI_COMM_WORLD
    if (present(communicator)) cart_comm = communicator

    ! Setup file access property list with parllel I/O access
    call h5pcreate_f(H5P_FILE_CREATE_F, plist_id, error) ! 1
    call assert(error == 0,'Error with h5pcreate_f() in hdf5_fcns.F90')
    call h5pset_fapl_mpio_f(plist_id, cart_comm, MPI_INFO_NULL, error) ! 2

    ! Create the file collectively
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id) ! 3
    call h5pclose_f(plist_id, error) ! 4

    ! Create the data space for the dataset
    call h5screate_simple_f(dset_rank, dimsf, dataspace, error) ! 5
    call h5screate_simple_f(dset_rank, dimsf, filespace, error) ! 7
    call h5screate_simple_f(dset_rank, chunk_dims, memspace, error) ! 6

    ! Select hyperslab in the file
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, &
                               error, stride, block) ! 8

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    ! Read the dataset collectively
    call h5dopen_f(file_id, trim(dsetname), dset_id, error, H5P_DEFAULT_F) ! 9
    call h5dread_f(dset_id, memtype, d_read, dimsf, error, &
                   mem_space_id = memspace, file_space_id = filespace,&
                   xfer_prp = plist_id) ! 10

    ! Close resources
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dataspace, error)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)    
  end subroutine

  subroutine read_h5_chunk_data4(block, chunk_dims, count, dimsf, dsetname, &
                                dset_rank, fname, memtype, offset, stride, &
                                d_read, communicator)
    !-------------------------------------------------------------------------!
    ! Read HDF5 data in parallel. This will read a 4 dimensional dataset      !
    !-------------------------------------------------------------------------!
    character(len=*), intent(in) :: fname
    real(rkind), dimension(:,:,:,:), intent(inout) :: d_read
    integer(HSIZE_T), dimension(:), intent(in) :: dimsf, chunk_dims, count, &
                                                  offset, stride, block
    integer(HID_T), intent(in) :: memtype
    integer, intent(in) :: dset_rank
    character(len=*), intent(in) :: dsetname
    integer, intent(in), optional :: communicator
    integer :: cart_comm
    integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace 
    integer :: error

    cart_comm = MPI_COMM_WORLD
    if (present(communicator)) cart_comm = communicator

    ! Setup file access property list with parllel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error) ! 1
    call h5pset_fapl_mpio_f(plist_id, cart_comm, MPI_INFO_NULL, error) ! 2

    ! Create the file collectively
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id) ! 3
    call h5pclose_f(plist_id, error) ! 4
    
    ! Create the data space for the dataset
    call h5screate_simple_f(dset_rank, dimsf, dataspace, error) ! 5
    call h5screate_simple_f(dset_rank, dimsf, filespace, error) ! 7
    call h5screate_simple_f(dset_rank, chunk_dims, memspace, error) ! 6

    ! Select hyperslab in the file
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, &
                               error, stride, block) ! 8

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    ! Read the dataset collectively
    call h5dopen_f(file_id, trim(dsetname), dset_id, error, H5P_DEFAULT_F) ! 9
    call h5dread_f(dset_id, memtype, d_read, dimsf, error, &
                   mem_space_id = memspace, file_space_id = filespace,&
                   xfer_prp = plist_id) ! 10

    ! Close resources
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dataspace, error)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)    
  end subroutine

  subroutine read_h5_chunk_data5(block, chunk_dims, count, dimsf, dsetname, &
                                dset_rank, fname, memtype, offset, stride, &
                                d_read, communicator)
    !-------------------------------------------------------------------------!
    ! Read HDF5 data in parallel. This will read a 5 dimensional dataset      !
    !-------------------------------------------------------------------------!
    character(len=*), intent(in) :: fname
    real(rkind), dimension(:,:,:,:,:), intent(inout) :: d_read
    integer(HSIZE_T), dimension(:), intent(in) :: dimsf, chunk_dims, count, &
                                                  offset, stride, block
    integer(HID_T), intent(in) :: memtype
    integer, intent(in) :: dset_rank
    character(len=*), intent(in) :: dsetname
    integer, intent(in), optional :: communicator
    integer :: cart_comm
    integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace 
    integer :: error

    cart_comm = MPI_COMM_WORLD
    if (present(communicator)) cart_comm = communicator

    ! Setup file access property list with parllel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error) ! 1
    call h5pset_fapl_mpio_f(plist_id, cart_comm, MPI_INFO_NULL, error) ! 2

    ! Create the file collectively
    call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id) ! 3
    call h5pclose_f(plist_id, error) ! 4

    ! Create the data space for the dataset
    call h5screate_simple_f(dset_rank, dimsf, dataspace, error) ! 5
    call h5screate_simple_f(dset_rank, dimsf, filespace, error) ! 7
    call h5screate_simple_f(dset_rank, chunk_dims, memspace, error) ! 6

    ! Select hyperslab in the file
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, &
                               error, stride, block) ! 8

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    ! Read the dataset collectively
    call h5dopen_f(file_id, trim(dsetname), dset_id, error, H5P_DEFAULT_F) ! 9
    call h5dread_f(dset_id, memtype, d_read, dimsf, error, &
                   mem_space_id = memspace, file_space_id = filespace,&
                   xfer_prp = plist_id) ! 10

    ! Close resources
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dataspace, error)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)    
  end subroutine

  subroutine write_h5_chunk_data2(block, chunk_dims, count, dimsf, dsetname, &
                                 dset_rank, fname, memtype, offset, stride, &
                                 d_write, communicator)
    character(len=*), intent(in) :: fname
    real(rkind), dimension(:,:), intent(in) :: d_write
    integer(HSIZE_T), dimension(:), intent(in) :: dimsf, chunk_dims, count, &
                                                  offset, stride, block
    integer(HID_T), intent(in) :: memtype
    character(len=*), intent(in) :: dsetname
    integer, intent(in), optional :: communicator
    integer :: cart_comm
    integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace
    integer, intent(in) :: dset_rank
    integer :: error
    
    cart_comm = MPI_COMM_WORLD
    if (present(communicator)) cart_comm = communicator

    ! Setup file access property list with parllel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error) ! 1
    call h5pset_fapl_mpio_f(plist_id, cart_comm, MPI_INFO_NULL, error) ! 2

    ! Create the file collectively
    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id) ! 3
    call h5pclose_f(plist_id, error) ! 4

    ! Create the data space for the dataset
    call h5screate_simple_f(dset_rank, dimsf, dataspace, error) ! 5 This is in
                                                          ! Hang's code, but not mine.
    call h5screate_simple_f(dset_rank, dimsf, filespace, error) ! 6
    call h5screate_simple_f(dset_rank, chunk_dims, memspace, error) ! 7

    ! Select hyperslab in the file
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, &
                               error, stride, block) ! 8

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) ! 9
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error) ! 10
    
    ! Write the dataset collectively
    call h5dcreate_f(file_id, trim(dsetname), memtype, dataspace, dset_id,&
                     error) ! 11 This is in Hang's code, but not mine.
    call h5dwrite_f(dset_id, memtype, d_write, dimsf, error, &
                    file_space_id = filespace, mem_space_id = memspace, &
                    xfer_prp = plist_id) ! 12

    ! Close resources
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)    

  end subroutine

  subroutine write_h5_chunk_data3(block, chunk_dims, count, dimsf, dsetname, &
                                 dset_rank, fname, memtype, offset, stride, &
                                 d_write, communicator)
    character(len=*), intent(in) :: fname
    real(rkind), dimension(:,:,:), intent(in) :: d_write
    integer(HSIZE_T), dimension(:), intent(in) :: dimsf, chunk_dims, count, &
                                                  offset, stride, block
    integer(HID_T), intent(in) :: memtype
    character(len=*), intent(in) :: dsetname
    integer, intent(in), optional :: communicator
    integer :: cart_comm
    integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace
    integer, intent(in) :: dset_rank
    integer :: error
    
    cart_comm = MPI_COMM_WORLD
    if (present(communicator)) cart_comm = communicator

    ! Setup file access property list with parllel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error) ! 1
    call h5pset_fapl_mpio_f(plist_id, cart_comm, MPI_INFO_NULL, error) ! 2

    ! Create the file collectively
    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id) ! 3
    call h5pclose_f(plist_id, error) ! 4

    ! Create the data space for the dataset
    call h5screate_simple_f(dset_rank, dimsf, dataspace, error) ! 5 This is in
                                                          ! Hang's code, but not mine.
    call h5screate_simple_f(dset_rank, dimsf, filespace, error) ! 6
    call h5screate_simple_f(dset_rank, chunk_dims, memspace, error) ! 7

    ! Select hyperslab in the file
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, &
                               error, stride, block) ! 8

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) ! 9
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error) ! 10
    
    ! Write the dataset collectively
    call h5dcreate_f(file_id, trim(dsetname), memtype, dataspace, dset_id,&
                     error) ! 11 This is in Hang's code, but not mine.
    call h5dwrite_f(dset_id, memtype, d_write, dimsf, error, &
                    file_space_id = filespace, mem_space_id = memspace, &
                    xfer_prp = plist_id) ! 12

    ! Close resources
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)    

  end subroutine
  
  subroutine write_h5_chunk_data4(block, chunk_dims, count, dimsf, dsetname, &
                                 dset_rank, fname, memtype, offset, stride, &
                                 d_write, communicator)
    character(len=*), intent(in) :: fname
    real(rkind), dimension(:,:,:,:), intent(in) :: d_write
    integer(HSIZE_T), dimension(:), intent(in) :: dimsf, chunk_dims, count, &
                                                  offset, stride, block
    integer(HID_T), intent(in) :: memtype
    character(len=*), intent(in) :: dsetname
    integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace
    integer, intent(in) :: dset_rank
    integer, intent(in), optional :: communicator
    integer :: cart_comm
    integer :: error
    
    cart_comm = MPI_COMM_WORLD
    if (present(communicator)) cart_comm = communicator

    ! Setup file access property list with parllel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error) ! 1
    call h5pset_fapl_mpio_f(plist_id, cart_comm, MPI_INFO_NULL, error) ! 2

    ! Create the file collectively
    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id) ! 3
    call h5pclose_f(plist_id, error) ! 4

    ! Create the data space for the dataset
    call h5screate_simple_f(dset_rank, dimsf, dataspace, error) ! 5 This is in
                                                          ! Hang's code, but not mine.
    call h5screate_simple_f(dset_rank, dimsf, filespace, error) ! 6
    call h5screate_simple_f(dset_rank, chunk_dims, memspace, error) ! 7

    ! Select hyperslab in the file
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, &
                               error, stride, block) ! 8

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) ! 9
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error) ! 10
    
    ! Write the dataset collectively
    call h5dcreate_f(file_id, trim(dsetname), memtype, dataspace, dset_id,&
                     error) ! 11 This is in Hang's code, but not mine.
    call h5dwrite_f(dset_id, memtype, d_write, dimsf, error, &
                    file_space_id = filespace, mem_space_id = memspace, &
                    xfer_prp = plist_id) ! 12

    ! Close resources
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)    

  end subroutine

  subroutine get_dataset_dimensions(fname,dsetname,ndims,dims)
      integer(HSIZE_T), dimension(ndims), intent(out) :: dims
      character(len=clen), intent(in) :: fname, dsetname
      integer(HSIZE_T), dimension(ndims) :: maxdims
      integer(HID_T) :: dspace_id, dset_id, file_id
      integer :: ndims, error

      call H5Fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
      call H5Dopen_f(file_id, trim(dsetname), dset_id, error)
      call H5Dget_space_f(dset_id, dspace_id, error) 
      call H5Sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)

      call H5Fclose_f(file_id,error)
      call H5Dclose_f(dset_id,error)
      call H5Sclose_f(dspace_id,error)

  end subroutine
end module hdf5_fcns
