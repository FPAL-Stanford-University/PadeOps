program JHU_channel_dumpPlanes
   use kind_parameters, only: rkind, clen
   use mpi
   use decomp_2d
   use timer, only: tic, toc
   use exits, only: message, gracefulexit
   use decomp_2d_io

   implicit none
   character(len=clen) :: inputfile, inputdir, outputdir, tempname, fname
   integer :: xyPlaneID, yzPlaneID, xzPlaneID, tid
   type(decomp_info) :: gp
   integer, parameter :: nx = 2048, ny = 1536, nz = 512
   real(rkind), dimension(:,:,:), allocatable :: field   
   integer :: ierr, ioUnit
   logical :: alsoDumpFiltered = .false., dump_xyPlane, dump_yzPlane, dump_xzPlane 

   namelist /INPUT/ inputdir, outputdir, tid, dump_xyPlane, dump_xzPlane, dump_yzPlane, xyPlaneID, &
                   & xzPlaneID, yzPlaneID, alsoDumpFiltered

   
   call MPI_init(ierr)

   call GETARG(1,inputfile)
   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=ioUnit, NML=INPUT)
   close(ioUnit)
   
   
   call decomp_2d_init(nx,ny,nz,0,0,[.TRUE.,.TRUE.,.FALSE.])
   call get_decomp_info(gp)
   allocate(field(gp%xsz(1),gp%xsz(2),gp%xsz(3)))

   call message(0,"Processing U velocity")
   call doStuff("u")

   call message(0,"Processing V velocity")
   call doStuff("v")

   call message(0,"Processing W velocity")
   call doStuff("w")
   
   call MPI_Barrier(mpi_comm_world, ierr)
   call MPI_Finalize(ierr)   

contains

   subroutine doStuff(fieldlabel)
      character(len=1), intent(in) :: fieldlabel

      call message(1,"Reading Baseline Field.")
      write(tempname,"(A1,A2,I3.3,A15)") fieldLabel,"_t",tid,"_JHUchannel.bin"
      fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
      call decomp_2d_read_one(1,field,fname, gp)

      if (dump_xyPlane) then
         call message(1,"Dumping xy plane (Baseline).")
         write(tempname,"(A1,A2,I3.3,A8,I4.4,A15)") fieldLabel,"_t",tid,"_xyPlane",xyPlaneID,"_JHUchannel.pln"
         fname = OutputDir(:len_trim(InputDir))//"/"//trim(tempname)
         call decomp_2d_write_plane(1,field,3, xyPlaneID, fname)
      end if

      if (dump_xzPlane) then
         call message(1,"Dumping xz plane (Baseline).")
         write(tempname,"(A1,A2,I3.3,A8,I4.4,A15)") fieldLabel,"_t",tid,"_xzPlane",xzPlaneID,"_JHUchannel.pln"
         fname = OutputDir(:len_trim(InputDir))//"/"//trim(tempname)
         call decomp_2d_write_plane(1,field,2, xzPlaneID, fname)
      end if

      if (dump_yzPlane) then
         call message(1,"Dumping yz plane (Baseline).")
         write(tempname,"(A1,A2,I3.3,A8,I4.4,A15)") fieldLabel,"_t",tid,"_yzPlane",xyPlaneID,"_JHUchannel.pln"
         fname = OutputDir(:len_trim(InputDir))//"/"//trim(tempname)
         call decomp_2d_write_plane(1,field,1, yzPlaneID, fname)
      end if

      if (alsoDumpFiltered) then
         call message(1,"Reading Filtered Field.")
         write(tempname,"(A1,A1,A2,I3.3,A15)") fieldLabel,"F","_t",tid,"_JHUchannel.bin"
         fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
         call decomp_2d_read_one(1,field,fname, gp)

         if (dump_xyPlane) then
            call message(1,"Dumping xy plane (Filtered).")
            write(tempname,"(A1,A1,A2,I3.3,A8,I4.4,A15)") fieldLabel,"F","_t",tid,"_xyPlane",xyPlaneID,"_JHUchannel.pln"
            fname = OutputDir(:len_trim(InputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,field,3, xyPlaneID, fname)
         end if

         if (dump_xzPlane) then
            call message(1,"Dumping xz plane (Filtered).")
            write(tempname,"(A1,A1,A2,I3.3,A8,I4.4,A15)") fieldLabel,"F","_t",tid,"_xzPlane",xzPlaneID,"_JHUchannel.pln"
            fname = OutputDir(:len_trim(InputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,field,2, xzPlaneID, fname)
         end if

         if (dump_yzPlane) then
            call message(1,"Dumping yz plane (Filtered).")
            write(tempname,"(A1,A1,A2,I3.3,A8,I4.4,A15)") fieldLabel,"F","_t",tid,"_yzPlane",xyPlaneID,"_JHUchannel.pln"
            fname = OutputDir(:len_trim(InputDir))//"/"//trim(tempname)
            call decomp_2d_write_plane(1,field,1, yzPlaneID, fname)
         end if
      end if
   end subroutine
end program
