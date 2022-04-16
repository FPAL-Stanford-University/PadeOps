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
  use domainSetup, only: getStartAndEndIndices, decomp2Dpencil
  use exits, only: message
  use fortran_assert, only: assert

  implicit none
  interface readFields
    module procedure read3Fields3D, read2Fields3D, read1Field5D
  end interface
  interface writeFields
    module procedure write3Fields3D, write1Field3D, write1Field5D
  end interface
  contains
    subroutine read1Field5D(fname,f1,dsetnm1,gp)
      ! Use this to read in large scale velocity data

      character(len=*), intent(in) :: fname
      real(rkind), dimension(:,:,:,:,:), intent(inout) :: f1
      character(len=*), intent(in) :: dsetnm1
      type(decomp_info), intent(in) :: gp
      integer :: nx, ny, nz
      integer :: ist, ien, jst, jen, kst, ken
      integer :: isz, jsz, ksz
  
      ! HDF5 variables
      integer(HID_T) :: memtype
      integer(HSIZE_T), dimension(5) :: dimsf, chunk_dims, stride, count, block, offset
      integer :: dset_rank

      call assert(decomp2Dpencil == 'x','IO routines only support x-pencil'//&
        ' decomposition')
      nx = gp%xsz(1)
      ny = gp%ysz(2)
      nz = gp%zsz(3)

      memtype = H5T_NATIVE_DOUBLE
      dset_rank = 5
      dimsf = [3,3,nx,ny,nz]
      
      ! Get the partition dimensions. Note, we call getStartAndEndIndices for
      ! generality since we do not assume a particular orientation for the
      ! domain decomposition a priori
      call getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken,isz,jsz,ksz) 
      chunk_dims = [3,3,isz,jsz,ksz]

      stride = 1
      count = 1
      block = chunk_dims
      offset(1:3) = 0
      offset(4) = gp%xst(2)-1
      offset(5) = gp%xst(3)-1
      
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetnm1, &
        dset_rank, trim(fname), memtype, offset, stride, f1, MPI_COMM_WORLD)
    end subroutine
    
    subroutine read2Fields3D(fname,f1,f2,dsetnm1,dsetnm2,gp)
      ! Use this to read in large scale velocity data

      character(len=*), intent(in) :: fname
      real(rkind), dimension(:,:,:), intent(inout) :: f1, f2
      character(len=*), intent(in) :: dsetnm1, dsetnm2
      type(decomp_info), intent(in) :: gp
      integer :: nx, ny, nz
      integer :: ist, ien, jst, jen, kst, ken
      integer :: isz, jsz, ksz
  
      ! HDF5 variables
      integer(HID_T) :: memtype
      integer(HSIZE_T), dimension(3) :: dimsf, chunk_dims, stride, count, block, offset
      integer :: dset_rank

      call assert(decomp2Dpencil == 'x','IO routines only support x-pencil'//&
        ' decomposition')
      nx = gp%xsz(1)
      ny = gp%ysz(2)
      nz = gp%zsz(3)

      memtype = H5T_NATIVE_DOUBLE
      dset_rank = 3
      dimsf = [nx,ny,nz]
      
      ! Get the partition dimensions. Note, we call getStartAndEndIndices for
      ! generality since we do not assume a particular orientation for the
      ! domain decomposition a priori
      call getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken,isz,jsz,ksz) 
      chunk_dims = [isz,jsz,ksz]

      stride = 1
      count = 1
      block = chunk_dims
      offset(1) = 0
      offset(2) = gp%xst(2)-1
      offset(3) = gp%xst(3)-1
      
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetnm1, &
        dset_rank, trim(fname), memtype, offset, stride, f1, MPI_COMM_WORLD)
    
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetnm2, &
        dset_rank, trim(fname), memtype, offset, stride, f2, MPI_COMM_WORLD)
    end subroutine
    
    subroutine read3Fields3D(fname,f1,f2,f3,dsetnm1,dsetnm2,dsetnm3,gp)
      ! Use this to read in large scale velocity data

      character(len=*), intent(in) :: fname
      real(rkind), dimension(:,:,:), intent(inout) :: f1, f2, f3
      character(len=*), intent(in) :: dsetnm1, dsetnm2, dsetnm3
      type(decomp_info), intent(in) :: gp
      integer :: nx, ny, nz
      integer :: ist, ien, jst, jen, kst, ken
      integer :: isz, jsz, ksz
  
      ! HDF5 variables
      integer(HID_T) :: memtype
      integer(HSIZE_T), dimension(3) :: dimsf, chunk_dims, stride, count, block, offset
      integer :: dset_rank

      call assert(decomp2Dpencil == 'x','IO routines only support x-pencil'//&
        ' decomposition')
      nx = gp%xsz(1)
      ny = gp%ysz(2)
      nz = gp%zsz(3)

      memtype = H5T_NATIVE_DOUBLE
      dset_rank = 3
      dimsf = [nx,ny,nz]
      
      ! Get the partition dimensions. Note, we call getStartAndEndIndices for
      ! generality since we do not assume a particular orientation for the
      ! domain decomposition a priori
      call getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken,isz,jsz,ksz) 
      chunk_dims = [isz,jsz,ksz]

      stride = 1
      count = 1
      block = chunk_dims
      offset(1) = 0
      offset(2) = gp%xst(2)-1
      offset(3) = gp%xst(3)-1
      
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetnm1, &
        dset_rank, trim(fname), memtype, offset, stride, f1, MPI_COMM_WORLD)
    
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetnm2, &
        dset_rank, trim(fname), memtype, offset, stride, f2, MPI_COMM_WORLD)

      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetnm3, &
        dset_rank, trim(fname), memtype, offset, stride, f3, MPI_COMM_WORLD)
    end subroutine
    
    subroutine write1field3D(fname,f1,dsetnm1,gp)
      ! See read3Fields3D for variable descriptions

      character(len=*), intent(in) :: fname
      real(rkind), dimension(:,:,:), intent(in) :: f1
      type(decomp_info), intent(in) :: gp
      character(len=*), intent(in) :: dsetnm1
      integer :: ist, ien, jst, jen, kst, ken
      integer :: isz, jsz, ksz
      integer :: nx, ny, nz ! <-- global mesh size
  
      ! HDF5 variables
      integer(HID_T) :: memtype
      integer(HSIZE_T), dimension(3) :: dimsf, chunk_dims, stride, count, block, offset
      character(len=2) :: dsetname
      integer :: dset_rank, error
      integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace

      call assert(decomp2Dpencil == 'x','IO routines only support x-pencil'//&
        ' decomposition')
      nx = gp%xsz(1)
      ny = gp%ysz(2)
      nz = gp%zsz(3)

      memtype = H5T_NATIVE_DOUBLE
      dset_rank = 3
      dimsf = [nx,ny,nz]

      ! Get the partition dimensions. Note, we call getStartAndEndIndices for
      ! generality since we do not assume a particular orientation for the
      ! domain decomposition a priori
      call getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      chunk_dims = [isz,jsz,ksz]

      stride = 1
      count = 1
      block = chunk_dims
      offset(1) = 0
      offset(2) = gp%xst(2)-1
      offset(3) = gp%xst(3)-1
 
      ! Create the file and other HDF5 stuff without creating the dataspace or
      ! writing the data
      call createAndOpenFile(block, chunk_dims, count, dimsf, dset_rank, fname, &
        offset, stride, MPI_COMM_WORLD, file_id, plist_id, filespace, &
        memspace, dataspace)

      ! Now write all velocity components to a single file
      dsetname = dsetnm1 
      call write_h5_chunk_data(file_id, dsetname, memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, f1)
  
      ! Close resources
      call closeFileResources(filespace,memspace,dset_id,plist_id)
      call h5fclose_f(file_id, error)    
    end subroutine

    subroutine write3fields3D(fname,f1,f2,f3,dsetnm1,dsetnm2,dsetnm3,gp)
      ! See read3Fields3D for variable descriptions

      character(len=*), intent(in) :: fname
      real(rkind), dimension(:,:,:), intent(in) :: f1, f2, f3
      type(decomp_info), intent(in) :: gp
      character(len=*), intent(in) :: dsetnm1, dsetnm2, dsetnm3
      integer :: ist, ien, jst, jen, kst, ken
      integer :: isz, jsz, ksz
      integer :: nx, ny, nz ! <-- global mesh size
  
      ! HDF5 variables
      integer(HID_T) :: memtype
      integer(HSIZE_T), dimension(3) :: dimsf, chunk_dims, stride, count, block, offset
      character(len=2) :: dsetname
      integer :: dset_rank, error
      integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace

      call assert(decomp2Dpencil == 'x','IO routines only support x-pencil'//&
        ' decomposition')
      nx = gp%xsz(1)
      ny = gp%ysz(2)
      nz = gp%zsz(3)

      memtype = H5T_NATIVE_DOUBLE
      dset_rank = 3
      dimsf = [nx,ny,nz]

      ! Get the partition dimensions. Note, we call getStartAndEndIndices for
      ! generality since we do not assume a particular orientation for the
      ! domain decomposition a priori
      call getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      chunk_dims = [isz,jsz,ksz]

      stride = 1
      count = 1
      block = chunk_dims
      offset(1) = 0
      offset(2) = gp%xst(2)-1
      offset(3) = gp%xst(3)-1
 
      ! Create the file and other HDF5 stuff without creating the dataspace or
      ! writing the data
      call createAndOpenFile(block, chunk_dims, count, dimsf, dset_rank, fname, &
        offset, stride, MPI_COMM_WORLD, file_id, plist_id, filespace, &
        memspace, dataspace)

      ! Now write all velocity components to a single file
      dsetname = dsetnm1 
      call write_h5_chunk_data(file_id, dsetname, memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, f1)
  
      dsetname = dsetnm2
      call write_h5_chunk_data(file_id, dsetname, memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, f2)
  
      dsetname = dsetnm3
      call write_h5_chunk_data(file_id, dsetname, memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, f3)
    
      ! Close resources
      call closeFileResources(filespace,memspace,dset_id,plist_id)
      call h5fclose_f(file_id, error)    
    end subroutine

    subroutine write1Field5D(fname,f,dsetname,gp)
      ! See read3Fields3D for variable descriptions

      character(len=*), intent(in) :: fname
      real(rkind), dimension(:,:,:,:,:), intent(in) :: f
      type(decomp_info), intent(in) :: gp
      integer :: ist, ien, jst, jen, kst, ken
      integer :: isz, jsz, ksz
      integer :: nx, ny, nz ! <-- global mesh size
  
      ! HDF5 variables
      integer(HID_T) :: memtype
      integer(HSIZE_T), dimension(5) :: dimsf, chunk_dims, stride, count, block, offset
      character(len=*), intent(in) :: dsetname
      integer :: dset_rank, error
      integer(HID_T) :: file_id, dset_id, plist_id, filespace, memspace, dataspace

      call assert(decomp2Dpencil == 'x','IO routines only support x-pencil'//&
        ' decomposition')
      nx = gp%xsz(1)
      ny = gp%ysz(2)
      nz = gp%zsz(3)

      memtype = H5T_NATIVE_DOUBLE
      dset_rank = 5
      dimsf = [3,3,nx,ny,nz]

      ! Get the partition dimensions. Note, we call getStartAndEndIndices for
      ! generality since we do not assume a particular orientation for the
      ! domain decomposition a priori
      call getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken,isz,jsz,ksz) 
      chunk_dims = [3,3,isz,jsz,ksz]

      stride = 1
      count = 1
      block = chunk_dims
      offset(1:3) = 0
      offset(4) = gp%xst(2)-1
      offset(5) = gp%xst(3)-1
 
      ! Create the file and other HDF5 stuff without creating the dataspace or
      ! writing the data
      call createAndOpenFile(block, chunk_dims, count, dimsf, dset_rank, fname, &
        offset, stride, MPI_COMM_WORLD, file_id, plist_id, filespace, &
        memspace, dataspace)

      ! Now write all velocity components to a single file
      call write_h5_chunk_data(file_id, dsetname, memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, f)
  
      ! Close resources
      call closeFileResources(filespace,memspace,dset_id,plist_id)
      call h5fclose_f(file_id, error)    
    end subroutine
end module
