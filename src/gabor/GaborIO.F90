module gaborIO_mod
  use kind_parameters, only: rkind
  use decomp_2d
  use hdf5
  use hdf5_fcns, only: read_h5_chunk_data, write_h5_chunk_data, &
    createAndOpenFile, closeFileResources
  use mpi, only: MPI_COMM_WORLD
  use domainSetup, only: getStartAndEndIndices
  use exits, only: message

  implicit none
  interface writeFields
    module procedure write3Fields3D, write1Field3D, write1Field5D
  end interface
  contains
    subroutine readVelocityField(fname,u,v,w,gp)
      ! Use this to read in large scale velocity data

      character(len=*), intent(in) :: fname
      real(rkind), dimension(:,:,:), intent(inout) :: u, v, w
      type(decomp_info), intent(in) :: gp
      integer :: nx, ny, nz
      integer :: ist, ien, jst, jen, kst, ken
      integer :: isz, jsz, ksz
  
      ! HDF5 variables
      integer(HID_T) :: memtype
      integer(HSIZE_T), dimension(3) :: dimsf, chunk_dims, stride, count, block, offset
      character(len=2) :: dsetname
      character(len=6) :: dsetname2
      integer :: dset_rank

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

      chunk_dims = [gp%xsz(1),gp%xsz(2),gp%xsz(3)]
      stride = 1
      count = 1
      block = chunk_dims
      call message("NOTE: 'x'-pencil decomposition assumed in HDF5 routines")
      offset(1) = 0
      offset(2) = gp%xst(2)-1
      offset(3) = gp%xst(3)-1
      
      dsetname = '/u'
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetname, &
        dset_rank, trim(fname), memtype, offset, stride, u, MPI_COMM_WORLD)
    
      dsetname = '/v'
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetname, &
        dset_rank, trim(fname), memtype, offset, stride, v, MPI_COMM_WORLD)

      dsetname = '/w'
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetname, &
        dset_rank, trim(fname), memtype, offset, stride, w, MPI_COMM_WORLD)
    end subroutine
    
    subroutine write1field3D(fname,f1,dsetnm1,gp)
      ! See readVelocityField for variable describptions

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
      call message("NOTE: 'x'-pencil decomposition assumed in HDF5 routines")
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
      ! See readVelocityField for variable describptions

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
      call message("NOTE: 'x'-pencil decomposition assumed in HDF5 routines")
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
      ! See readVelocityField for variable describptions

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
      call message("NOTE: 'x'-pencil decomposition assumed in HDF5 routines")
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
