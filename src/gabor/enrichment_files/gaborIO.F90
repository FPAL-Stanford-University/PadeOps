subroutine dumpDataAll(this,x,y,z,kx,ky,kz,uR,uI,vR,vI,wR,wI,KE,L)
  use decomp_2d, only: nrank
  use basic_io,  only: write_1D_ascii
  class(enrichmentOperator), intent(inout) :: this
  real(rkind), dimension(:), intent(in) :: x, y, z, kx, ky, kz, uR, uI, vR, &
    vI, wR, wI, KE, L
  character(len=clen) :: fname

  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_uhatR.txt'
  call write_1D_ascii(uR,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_uhatI.txt'
  call write_1D_ascii(uI,fname)
  
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_vhatR.txt'
  call write_1D_ascii(vR,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_vhatI.txt'
  call write_1D_ascii(vI,fname)
  
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_whatR.txt'
  call write_1D_ascii(wR,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_whatI.txt'
  call write_1D_ascii(wI,fname)

  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_kx.txt'
  call write_1D_ascii(kx,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_ky.txt'
  call write_1D_ascii(ky,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_kz.txt'
  call write_1D_ascii(kz,fname)

  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_x.txt'
  call write_1D_ascii(x,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_y.txt'
  call write_1D_ascii(y,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_z.txt'
  call write_1D_ascii(z,fname)
  
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_L.txt'
  call write_1D_ascii(L,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_KE.txt'
  call write_1D_ascii(KE,fname)
end subroutine

subroutine dumpDataLoc(this,x,y,z,fdesc)
  use decomp_2d, only: nrank
  use basic_io,  only: write_1D_ascii
  class(enrichmentOperator), intent(inout) :: this
  real(rkind), dimension(:), intent(in) :: x, y, z
  character(len=*), intent(in), optional :: fdesc
  character(len=clen) :: fname

  if (present(fdesc)) then
    write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_'//trim(fdesc)//'_x.txt'
    call write_1D_ascii(x,fname)
    write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_'//trim(fdesc)//'_y.txt'
    call write_1D_ascii(y,fname)
    write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_'//trim(fdesc)//'_z.txt'
    call write_1D_ascii(z,fname)
  else
    write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_x.txt'
    call write_1D_ascii(x,fname)
    write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_y.txt'
    call write_1D_ascii(y,fname)
    write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_z.txt'
    call write_1D_ascii(z,fname)
  end if
end subroutine

subroutine dumpMeshDetails(QHxE, QHyE, QHzE, QHgID, largeScalesMesh, smallScalesMesh, outputdir)
  use basic_io,  only: write_1D_ascii
  real(rkind), dimension(:), intent(in) :: QHxE, QHyE, QHzE
  integer, dimension(:,:,:), intent(in) :: QHgID
  real(rkind), dimension(:,:,:,:), intent(in), target :: largeScalesMesh, smallScalesMesh
  character(len=*), intent(in) :: outputdir
  real(rkind), dimension(:), pointer :: xL, yL, zL, xS, yS, zS
  real(rkind), dimension(:), allocatable :: QHgIDvec
  character(len=clen) :: fname

  xL => largeScalesMesh(:,1,1,1)
  yL => largeScalesMesh(1,:,1,2)
  zL => largeScalesMesh(1,1,:,3)

  xS => smallScalesMesh(:,1,1,1)
  yS => smallScalesMesh(1,:,1,2)
  zS => smallScalesMesh(1,1,:,3)

  QHgIDvec = real(reshape(QHgID,[size(QHgID)]),rkind)
  write(fname,'(A,I2.2,A)')trim(outputdir)//'/rank',nrank,'_QHxE.txt'
  call write_1D_ascii(QHxE,fname)
  write(fname,'(A,I2.2,A)')trim(outputdir)//'/rank',nrank,'_QHyE.txt'
  call write_1D_ascii(QHyE,fname)
  write(fname,'(A,I2.2,A)')trim(outputdir)//'/rank',nrank,'_QHzE.txt'
  call write_1D_ascii(QHzE,fname)
  write(fname,'(A,I2.2,A)')trim(outputdir)//'/rank',nrank,'_QHgID.txt'
  call write_1D_ascii(QHgIDvec,fname)

  write(fname,'(A,I2.2,A)')trim(outputdir)//'/rank',nrank,'_xL.txt'
  call write_1D_ascii(xL,fname)
  write(fname,'(A,I2.2,A)')trim(outputdir)//'/rank',nrank,'_yL.txt'
  call write_1D_ascii(yL,fname)
  write(fname,'(A,I2.2,A)')trim(outputdir)//'/rank',nrank,'_zL.txt'
  call write_1D_ascii(zL,fname)

  write(fname,'(A,I2.2,A)')trim(outputdir)//'/rank',nrank,'_xS.txt'
  call write_1D_ascii(xS,fname)
  write(fname,'(A,I2.2,A)')trim(outputdir)//'/rank',nrank,'_yS.txt'
  call write_1D_ascii(yS,fname)
  write(fname,'(A,I2.2,A)')trim(outputdir)//'/rank',nrank,'_zS.txt'
  call write_1D_ascii(zS,fname)

  nullify(xL, yL, zL, xS, yS, zS)
  deallocate(QHgIDvec)
end subroutine
