subroutine dumpData(this)
  use decomp_2d, only: nrank
  use basic_io,  only: write_1D_ascii
  class(enrichmentOperator), intent(inout) :: this
  character(len=clen) :: fname

  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_uhatR.out'
  call write_1D_ascii(this%uhatR,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_uhatI.out'
  call write_1D_ascii(this%uhatI,fname)
  
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_vhatR.out'
  call write_1D_ascii(this%vhatR,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_vhatI.out'
  call write_1D_ascii(this%vhatI,fname)
  
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_whatR.out'
  call write_1D_ascii(this%whatR,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_whatI.out'
  call write_1D_ascii(this%whatI,fname)

  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_kx.out'
  call write_1D_ascii(this%kx,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_ky.out'
  call write_1D_ascii(this%ky,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_kz.out'
  call write_1D_ascii(this%kz,fname)

  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_x.out'
  call write_1D_ascii(this%x,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_y.out'
  call write_1D_ascii(this%y,fname)
  write(fname,'(A,I2.2,A)')trim(this%outputdir)//'/rank',nrank,'_z.out'
  call write_1D_ascii(this%z,fname)

end subroutine
