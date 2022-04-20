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
      offset(1:2) = 0
      offset(3) = ist - 1
      offset(4) = jst - 1
      offset(5) = kst - 1
      
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
      offset(1) = ist - 1
      offset(2) = jst - 1
      offset(3) = kst - 1
      
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
      offset(1) = ist - 1
      offset(2) = jst - 1
      offset(3) = kst - 1
      
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
      offset(1) = ist - 1
      offset(2) = jst - 1
      offset(3) = kst - 1
 
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
      call closeFileResources(filespace,memspace,plist_id)
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
      offset(1) = ist - 1
      offset(2) = jst - 1
      offset(3) = kst - 1
 
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
      call closeFileResources(filespace,memspace,plist_id)
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
      offset(1:2) = 0
      offset(3) = ist - 1
      offset(4) = jst - 1
      offset(5) = kst - 1
 
      ! Create the file and other HDF5 stuff without creating the dataspace or
      ! writing the data
      call createAndOpenFile(block, chunk_dims, count, dimsf, dset_rank, fname, &
        offset, stride, MPI_COMM_WORLD, file_id, plist_id, filespace, &
        memspace, dataspace)

      ! Now write all velocity components to a single file
      call write_h5_chunk_data(file_id, dsetname, memtype, dataspace, &
      dset_id, dimsf, filespace, memspace, plist_id, f)
  
      ! Close resources
      call closeFileResources(filespace,memspace,plist_id)
      call h5fclose_f(file_id, error)    
    end subroutine

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

      ! Get the partition dimensions. Note, we call getStartAndEndIndices for
      ! generality since we do not assume a particular orientation for the
      ! domain decomposition a priori
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
