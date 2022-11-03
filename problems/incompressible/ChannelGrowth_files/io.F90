module ChannelGrowth_IO

    use kind_parameters, only: rkind, clen
    use decomp_2d,       only: decomp_info,nrank,nproc
    use basic_io,        only: read_2d_ascii 
    use exits,           only: message, gracefulExit
    implicit none

contains

subroutine get_perturbations(gp, x, y, kx, ky, InitFileName, u, v, w)
   use constants, only: imi, pi   
   use decomp_2d
   type(decomp_info), intent(in) :: gp
   character(len=*), intent(in) :: InitFileName
   real(Rkind), intent(in) :: kx, ky 
   real(rkind), dimension(:,:,:), intent(in) :: x, y
   real(rkind), dimension(:,:,:), intent(out) :: u, v, w
   real(rkind), dimension(:,:), allocatable :: data2read
   integer :: k
   complex(rkind)  :: uhat, what
   real(rkind), dimension(:,:), allocatable :: xinZ, yinZ
   real(rkind), dimension(:,:,:), allocatable :: buffy, buffz1, buffz2


   call read_2d_ascii(data2read,InitFileName)
    
   if (size(data2read,1) .ne. gp%zsz(3)) then
    call gracefulExit("Incorrect nz for the mode file",23)
   end if 

   allocate(buffy(gp%ysz(1),gp%ysz(2),gp%ysz(3)))
   allocate(buffz1(gp%zsz(1),gp%zsz(2),gp%zsz(3)))
   allocate(buffz2(gp%zsz(1),gp%zsz(2),gp%zsz(3)))
   allocate(xinZ(gp%zsz(1),gp%zsz(2)))
   allocate(yinZ(gp%zsz(1),gp%zsz(2)))

   call transpose_x_to_y(x,buffy,gp)
   call transpose_y_to_z(buffy,buffz1,gp)
   xinZ = buffz1(:,:,1) 
   
   call transpose_x_to_y(y,buffy,gp)
   call transpose_y_to_z(buffy,buffz1,gp)
   yinZ = buffz1(:,:,1) 

   v = 0.d0
   do k =  1,gp%zsz(3)
       uhat = data2read(k,2) + imi*data2read(k,3)
       what = data2read(k,4) + imi*data2read(k,5)
       
       buffz1(:,:,k) = real(uhat*exp(imi*kx*xinZ + imi*ky*yinZ),rkind)
       buffz2(:,:,k) = real(what*exp(imi*kx*xinZ + imi*ky*yinZ),rkind)

   end do 

   call transpose_z_to_y(buffz1,buffy,gp)
   call transpose_y_to_x(buffy,u)
   
   call transpose_z_to_y(buffz2,buffy,gp)
   call transpose_y_to_x(buffy,w)

   deallocate(buffy, buffz1, buffz2)
   deallocate(xinZ, yinZ)

end subroutine 

end module 
