module gaborIO_mod
  ! This module is miscellaneous IO routines specific to the enrichment
  ! algorithm. This includes reading large-scale velocity fields. Perhaps it
  ! makes more sense for that to reside in the largeScales.F90 file, but for
  ! now, it is here.
  use kind_parameters, only: rkind
  use decomp_2d
  use hdf5
  use hdf5_fcns, only: read_h5_chunk_data, write_h5_chunk_data, &
    createAndOpenFile, closeFileResources
  use mpi, only: MPI_COMM_WORLD
  use domainSetup, only: getStartAndEndIndices
  use fortran_assert, only: assert

  implicit none
  contains

    subroutine writeModes(fname)
      use GaborModeRoutines, only: nmodes, nmodesALL, &
        uhatR, uhatI, vhatR, vhatI, whatR, whatI, &
        kx, ky, kz, gmxloc, gmyloc, gmzloc
      use decomp_2D, only: nrank
      ! Write Gabor mode info to disk
      character(len=*), intent(in) :: fname
  
      ! HDF5 variables
      integer(HID_T) :: memtype
      integer(HSIZE_T), dimension(1) :: dimsf, chunk_dims, stride, count, block, offset
      integer :: dset_rank, error
      integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace

      memtype = H5T_NATIVE_DOUBLE
      dset_rank = 1
      dimsf = [sum(nmodesALL)]

      ! Get the partition dimensions. Note, we call getStartAndEndIndices for
      ! generality since we do not assume a particular orientation for the
      ! domain decomposition a priori
      chunk_dims = [nmodes]

      stride = 1
      count = 1
      block = chunk_dims
      offset = sum(nmodesALL(1:(nrank+1)))-nmodes
 
      ! Create the file and other HDF5 stuff without creating the dataspace or
      ! writing the data
      call createAndOpenFile(block, chunk_dims, count, dimsf, dset_rank, fname, &
        offset, stride, MPI_COMM_WORLD, file_id, plist_id, filespace, &
        memspace, dataspace)

      call write_h5_chunk_data(file_id, 'uhatR', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, uhatR)
      call write_h5_chunk_data(file_id, 'uhatI', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, uhatI)
  
      call write_h5_chunk_data(file_id, 'vhatR', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, vhatR)
      call write_h5_chunk_data(file_id, 'vhatI', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, vhatI)
  
      call write_h5_chunk_data(file_id, 'whatR', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, whatR)
      call write_h5_chunk_data(file_id, 'whatI', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, whatI)
  
      call write_h5_chunk_data(file_id, 'kx', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, kx)
      call write_h5_chunk_data(file_id, 'ky', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, ky)
      call write_h5_chunk_data(file_id, 'kz', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, kz)
  
      call write_h5_chunk_data(file_id, 'gmxloc', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, gmxloc)
      call write_h5_chunk_data(file_id, 'gmyloc', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, gmyloc)
      call write_h5_chunk_data(file_id, 'gmzloc', memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, gmzloc)
      
      ! Close resources
      call closeFileResources(filespace,memspace,plist_id)
      call h5fclose_f(file_id, error)    
    end subroutine

    subroutine readModes(fname,uR,uI,vR,vI,wR,wI,k1,k2,k3,x,y,z)
      use GaborModeRoutines, only: nmodes, nmodesALL, &
        uhatR, uhatI, vhatR, vhatI, whatR, whatI, &
        kx, ky, kz, gmxloc, gmyloc, gmzloc
      use decomp_2D, only: nrank
      ! Read Gabor info from disk
      character(len=*), intent(in) :: fname
      real(rkind), dimension(:), intent(inout), optional :: uR, uI, vR, vI, wR,&
        wI, k1, k2, k3, x, y, z
  
      ! HDF5 variables
      integer(HID_T) :: memtype
      integer(HSIZE_T), dimension(1) :: dimsf, chunk_dims, stride, count, block, offset
      integer :: dset_rank

      memtype = H5T_NATIVE_DOUBLE
      dset_rank = 1
      dimsf = [sum(nmodesALL)]

      ! Get the partition dimension
      chunk_dims = [nmodes]

      stride = 1
      count = 1
      block = chunk_dims
      offset = sum(nmodesALL(1:(nrank+1)))-nmodes

      if (present(uR)) then
        call assert(present(uI),'Must provide inout array for uI data') 
        call assert(present(vR),'Must provide inout array for vR data') 
        call assert(present(vI),'Must provide inout array for vI data') 
        call assert(present(wR),'Must provide inout array for wR data') 
        call assert(present(wI),'Must provide inout array for wI data') 
        call assert(present(k1),'Must provide inout array for k1 data') 
        call assert(present(k2),'Must provide inout array for k2 data') 
        call assert(present(k3),'Must provide inout array for k3 data') 
        call assert(present(x),'Must provide inout array for x data') 
        call assert(present(y),'Must provide inout array for y data') 
        call assert(present(z),'Must provide inout array for z data')

        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/uhatR', &
          dset_rank, trim(fname), memtype, offset, stride, uR, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/uhatI', &
          dset_rank, trim(fname), memtype, offset, stride, uI, MPI_COMM_WORLD)

        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/vhatR', &
          dset_rank, trim(fname), memtype, offset, stride, vR, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/vhatI', &
          dset_rank, trim(fname), memtype, offset, stride, vI, MPI_COMM_WORLD)

        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/whatR', &
          dset_rank, trim(fname), memtype, offset, stride, wR, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/whatI', &
          dset_rank, trim(fname), memtype, offset, stride, wI, MPI_COMM_WORLD)

        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/kx', &
          dset_rank, trim(fname), memtype, offset, stride, k1, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/ky', &
          dset_rank, trim(fname), memtype, offset, stride, k2, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/kz', &
          dset_rank, trim(fname), memtype, offset, stride, k3, MPI_COMM_WORLD)

        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/gmxloc', &
          dset_rank, trim(fname), memtype, offset, stride, x, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/gmyloc', &
          dset_rank, trim(fname), memtype, offset, stride, y, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/gmzloc', &
          dset_rank, trim(fname), memtype, offset, stride, z, MPI_COMM_WORLD)
      else
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/uhatR', &
          dset_rank, trim(fname), memtype, offset, stride, uhatR, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/uhatI', &
          dset_rank, trim(fname), memtype, offset, stride, uhatI, MPI_COMM_WORLD)

        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/vhatR', &
          dset_rank, trim(fname), memtype, offset, stride, vhatR, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/vhatI', &
          dset_rank, trim(fname), memtype, offset, stride, vhatI, MPI_COMM_WORLD)

        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/whatR', &
          dset_rank, trim(fname), memtype, offset, stride, whatR, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/whatI', &
          dset_rank, trim(fname), memtype, offset, stride, whatI, MPI_COMM_WORLD)

        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/kx', &
          dset_rank, trim(fname), memtype, offset, stride, kx, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/ky', &
          dset_rank, trim(fname), memtype, offset, stride, ky, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/kz', &
          dset_rank, trim(fname), memtype, offset, stride, kz, MPI_COMM_WORLD)

        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/gmxloc', &
          dset_rank, trim(fname), memtype, offset, stride, gmxloc, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/gmyloc', &
          dset_rank, trim(fname), memtype, offset, stride, gmyloc, MPI_COMM_WORLD)
        call read_h5_chunk_data(block, chunk_dims, count, dimsf, '/gmzloc', &
          dset_rank, trim(fname), memtype, offset, stride, gmzloc, MPI_COMM_WORLD)
      end if
    end subroutine
end module
