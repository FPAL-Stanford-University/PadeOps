module initprocedures
   implicit none
contains
   subroutine meshgen(decomp, dx, dy, dz, mesh, inputfile)
       use kind_parameters,  only: rkind
       use constants,        only: one, two
       use decomp_2d,        only: decomp_info
       implicit none
   
       type(decomp_info),                                          intent(in)    :: decomp
       real(rkind),                                                intent(inout) :: dx,dy,dz
       real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
       real(rkind) :: z0init
       integer :: i,j,k, ioUnit
       character(len=*),                intent(in)    :: inputfile
       integer :: ix1, ixn, iy1, iyn, iz1, izn, nxg, nyg, nzg
       real(rkind)  :: Lx = one, Ly = one, Lz = one
       namelist /PBLINPUT/ Lx, Ly, Lz, z0init 
   
       ioUnit = 11
       open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
       read(unit=ioUnit, NML=PBLINPUT)
       close(ioUnit)    
   
       !Lx = two*pi; Ly = two*pi; Lz = one
   
       nxg = decomp%xsz(1); nyg = decomp%ysz(2); nzg = decomp%zsz(3)
   
       ! If base decomposition is in Y
       ix1 = decomp%xst(1); iy1 = decomp%xst(2); iz1 = decomp%xst(3)
       ixn = decomp%xen(1); iyn = decomp%xen(2); izn = decomp%xen(3)
       
       associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
   
           dx = Lx/real(nxg,rkind)
           dy = Ly/real(nyg,rkind)
           dz = Lz/real(nzg,rkind)
   
           do k=1,size(mesh,3)
               do j=1,size(mesh,2)
                   do i=1,size(mesh,1)
                       x(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                       y(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                       z(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
                   end do
               end do
           end do
   
           ! Shift everything to the origin 
           x = x - dx
           y = y - dy
           z = z - dz 
   
       end associate
   
   end subroutine


   subroutine writeVisualizationFile(tid, rid, f, gpC, label, InputDir)
       use decomp_2d_io
       use decomp_2d
       use kind_parameters, only: rkind, clen
       use mpi
       use exits, only: message
       type(decomp_info), intent(in) :: gpC
       integer, intent(in) :: tid, rid
       real(rkind), dimension(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)), intent(in) :: f
       character(len=clen) :: tempname, fname
       character(len=4), intent(in) :: label
       character(len=*), intent(in) :: InputDir

       write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_",label,"_t",tid,".out"
       fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(1,f,fname, gpC)
       
   end subroutine


   subroutine readVisualizationFile(tid, rid, u, v, wC, gpC, InputDir)
       use decomp_2d_io
       use decomp_2d
       use kind_parameters, only: rkind, clen
       use mpi
       use exits, only: message
       type(decomp_info), intent(in) :: gpC
       integer, intent(in) :: tid, rid
       real(rkind), dimension(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)), intent(out) :: u, v, wC
       character(len=clen) :: tempname, fname
       character(len=4) :: label
       character(len=*), intent(in) :: InputDir

       label = "uVel"
       write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_",label,"_t",tid,".out"
       fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
       call decomp_2d_read_one(1,u,fname, gpC)
       
       label = "vVel"
       write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_",label,"_t",tid,".out"
       fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
       call decomp_2d_read_one(1,v,fname, gpC)

       label = "wVel"
       write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_",label,"_t",tid,".out"
       fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
       call decomp_2d_read_one(1,wC,fname, gpC)
   
       call message("================= VISUALIZATION FILE READ  ======================")

   end subroutine

end module 
