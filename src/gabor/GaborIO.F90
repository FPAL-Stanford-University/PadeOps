module gaborIO_mod
  use kind_parameters, only: rkind
  use decomp_2d
  use hdf5
  use hdf5_fcns, only: read_h5_chunk_data
  use mpi, only: MPI_COMM_WORLD

  implicit none
  contains
    subroutine readLargeScaleData(fname,U,V,W,nxG,nyG,nzG,gp)
      ! Use this to read in large scale velocity data
      ! nxG, nyG, nzG --> total number of grid points for "G"lobal domain

      character(len=*), intent(in) :: fname
      real(rkind), dimension(:,:,:), intent(inout) :: U, V, W
      integer, intent(in) :: nxG, nyG, nzG
      type(decomp_info), intent(in) :: gp
  
      ! HDF5 variables
      integer(HID_T) :: memtype
      integer(HSIZE_T), dimension(3) :: dimsf, chunk_dims, stride, count, block, offset
      integer(HSIZE_T), dimension(5) :: dimsf2, chunk_dims2, stride2, count2, block2, offset2
      character(len=2) :: dsetname
      character(len=6) :: dsetname2
      integer :: dset_rank

      memtype = H5T_NATIVE_DOUBLE
      dset_rank = 3
      dimsf = [nxG,nyG,nzG]
      chunk_dims = [gp%xsz(1),gp%xsz(2),gp%xsz(3)]
      stride = 1
      count = 1
      block = chunk_dims
      offset(1) = 0
      offset(2) = gp%xst(2)-1
      offset(3) = gp%xst(3)-1
      
      dsetname = '/U'
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetname, &
        dset_rank, trim(fname), memtype, offset, stride, U, MPI_COMM_WORLD)
    
      dsetname = '/V'
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetname, &
        dset_rank, trim(fname), memtype, offset, stride, V, MPI_COMM_WORLD)

      dsetname = '/W'
      call read_h5_chunk_data(block, chunk_dims, count, dimsf, dsetname, &
        dset_rank, trim(fname), memtype, offset, stride, W, MPI_COMM_WORLD)
    end subroutine
end module
