module StratifiedShearLayer_IO

    use kind_parameters, only: rkind, clen
    use decomp_2d,       only: decomp_info,nrank,nproc
    use basic_io,        only: read_2d_ascii 
    use exits,           only: message
    implicit none

contains

subroutine read_Domain_info(Lx,Ly,Lz,fname)
   real(rkind), intent(out) :: Lx, Ly, Lz
   character(len=*), intent(in) :: fname
   
   real(rkind), dimension(:,:), allocatable :: data2read

   call read_2d_ascii(data2read,fname) 

   Lx = data2read(1,1)
   Ly = data2read(2,1)
   Lz = data2read(3,1)

   deallocate(data2read)
end subroutine


end module 
