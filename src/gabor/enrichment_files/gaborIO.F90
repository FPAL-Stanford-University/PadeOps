subroutine dumpData(this)
  use decomp_2d, only: nrank
  use basic_io,  only: write_1D_ascii
  class(enrichmentOperator), intent(inout) :: this
  character(len=clen) :: fname

  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_uhatR.txt'
  call write_1D_ascii(this%uhatR,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_uhatI.txt'
  call write_1D_ascii(this%uhatI,fname)
  
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_vhatR.txt'
  call write_1D_ascii(this%vhatR,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_vhatI.txt'
  call write_1D_ascii(this%vhatI,fname)
  
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_whatR.txt'
  call write_1D_ascii(this%whatR,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_whatI.txt'
  call write_1D_ascii(this%whatI,fname)

  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_kx.txt'
  call write_1D_ascii(this%kx,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_ky.txt'
  call write_1D_ascii(this%ky,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_kz.txt'
  call write_1D_ascii(this%kz,fname)

  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_x.txt'
  call write_1D_ascii(this%x,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_y.txt'
  call write_1D_ascii(this%y,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_z.txt'
  call write_1D_ascii(this%z,fname)
end subroutine
