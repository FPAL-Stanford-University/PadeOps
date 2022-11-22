subroutine dumpData(this,x,y,z,kx,ky,kz,uR,uI,vR,vI,wR,wI)
  use decomp_2d, only: nrank
  use basic_io,  only: write_1D_ascii
  class(enrichmentOperator), intent(inout) :: this
  real(rkind), dimension(:), intent(in) :: x, y, z, kx, ky, kz, uR, uI, vR, &
    vI, wR, wI
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
end subroutine
