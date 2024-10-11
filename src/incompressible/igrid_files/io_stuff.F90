   subroutine dumpRestartFile(this)
       use decomp_2d_io
       use mpi
       use exits, only: message
       use basic_io, only: write_1d_ascii
       class(igrid), intent(in) :: this
       character(len=clen) :: tempname, fname
       integer :: ierr, idx
       real(rkind), dimension(:), allocatable :: tmp 

       write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_u.",this%step
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(1,this%u,fname, this%gpC)
       ! Link "LATEST" restart file to the recently dumped file
       if (nrank == 0) then
           write(tempname,"(A7,A4,I2.2,A9)") "RESTART", "_Run",this%runID, "_u.LATEST"
           call execute_command_line('ln -sf '//trim(fname)//' '//trim(this%outputdir)//'/'//trim(tempname))
       end if

       write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_v.",this%step
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(1,this%v,fname, this%gpC)
       if (nrank == 0) then
           ! Link "LATEST" restart file to the recently dumped file
           write(tempname,"(A7,A4,I2.2,A9)") "RESTART", "_Run",this%runID, "_v.LATEST"
           call execute_command_line('ln -sf '//trim(fname)//' '//trim(this%outputdir)//'/'//trim(tempname))
       end if

       write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_w.",this%step
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(1,this%w,fname, this%gpE)
       if (nrank == 0) then
           ! Link "LATEST" restart file to the recently dumped file
           write(tempname,"(A7,A4,I2.2,A9)") "RESTART", "_Run",this%runID, "_w.LATEST"
           call execute_command_line('ln -sf '//trim(fname)//' '//trim(this%outputdir)//'/'//trim(tempname))
       end if

       if (this%isStratified .or. this%initspinup) then
           write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_T.",this%step
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(1,this%T,fname, this%gpC)
           if (nrank == 0) then
               ! Link "LATEST" restart file to the recently dumped file
               write(tempname,"(A7,A4,I2.2,A9)") "RESTART", "_Run",this%runID, "_T.LATEST"
               call execute_command_line('ln -sf '//trim(fname)//' '//trim(this%outputdir)//'/'//trim(tempname))
           end if
       end if 

       if (this%usescalars) then
           do idx = 1,this%n_scalars
              call this%scalars(idx)%dumprestart(this%step)
           end do
       end if

       if (this%localizedForceLayer > 0) then
           if (this%localizedForceLayer == 1) then
               if (nrank == 0) then
                   allocate(tmp(2))
                   tmp(1) = this%forceLayer%seedFact
                   tmp(2) = this%forceLayer%ampFact
                   write(tempname,"(A7,A4,I2.2,A7,I6.6)") "RESTART", "_Run",this%runID, "_finfo.",this%step
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call write_1d_ascii(tmp,fname)
                   deallocate(tmp)
               end if
               call MPI_Barrier(MPI_COMM_WORLD,ierr)
           
               write(tempname,"(A7,A4,I2.2,A4,I6.6)") "RESTART", "_Run",this%runID, "_fx.",this%step
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_one(1,this%forceLayer%fx,trim(fname), this%gpC)

               write(tempname,"(A7,A4,I2.2,A4,I6.6)") "RESTART", "_Run",this%runID, "_fy.",this%step
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_one(1,this%forceLayer%fy,fname, this%gpC)

               write(tempname,"(A7,A4,I2.2,A4,I6.6)") "RESTART", "_Run",this%runID, "_fz.",this%step
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_one(1,this%forceLayer%fz,fname, this%gpE)

           end if
       end if

       if (nrank == 0) then
           write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",this%runID, "_info.",this%step
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           OPEN(UNIT=10, FILE=trim(fname))
           write(10,"(100g15.5)") this%tsim
           write(10,"(100g15.5)") this%step
           if (this%useControl) then
              write(10,"(100g15.5)") this%G_alpha
           end if
           close(10)
           ! Link "LATEST" restart file to the recently dumped file
           write(tempname,"(A7,A4,I2.2,A12)") "RESTART", "_Run",this%runID, "_info.LATEST"
           call execute_command_line('ln -sf '//trim(fname)//' '//trim(this%outputdir)//'/'//trim(tempname))
       end if 

       call mpi_barrier(mpi_comm_world, ierr)
       call message(1, "Just Dumped a RESTART file")

   end subroutine

   ! NOTE: If you want to dump an edge field, you need to call in dumpFullField
   ! routine with this%gpE passed in as the 3rd argument. If it's a cell field,
   ! then you don't need to pass in any gp since the default gp is this%gpC
   subroutine dumpFullField(this,arr,label,gp2use,step,pencil)
       use decomp_2d_io
       use mpi
       use exits, only: message
       class(igrid), intent(in) :: this
       character(len=clen) :: tempname, fname
       real(rkind), dimension(:,:,:), intent(in) :: arr
       character(len=4), intent(in) :: label
       type(decomp_info), intent(in), optional :: gp2use
       integer, intent(in), optional :: step
       character(len=1), intent(in), optional :: pencil
       integer :: decomp_dir

       decomp_dir = 1
       if (present(pencil)) then
           select case (pencil)
           case ('x')
               decomp_dir = 1
           case ('y')
               decomp_dir = 2
           case ('z')
               decomp_dir = 3
           end select
       end if

       if (present(step)) then
           write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",     step,".out"
       else
           write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",this%step,".out"
       end if
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       if (present(gp2use)) then
          call decomp_2d_write_one(decomp_dir,arr,fname,gp2use)
       else
          call decomp_2d_write_one(decomp_dir,arr,fname,this%gpC)
       end if

   end subroutine
   subroutine dumpFullField_cmplx(this,arr,label,gp2use,step,pencil)
       use decomp_2d_io
       use mpi
       use exits, only: message
       class(igrid), intent(in) :: this
       character(len=clen) :: tempname, fname
       complex(rkind), dimension(:,:,:), intent(in) :: arr
       character(len=4), intent(in) :: label
       type(decomp_info), intent(in), optional :: gp2use
       integer, intent(in), optional :: step
       character(len=1), intent(in), optional :: pencil
       integer :: decomp_dir

       decomp_dir = 2
       if (present(pencil)) then
           select case (pencil)
           case ('x')
               decomp_dir = 1
           case ('y')
               decomp_dir = 2
           case ('z')
               decomp_dir = 3
           end select
       end if

        if (present(step)) then
            write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",     step,".out"
        else
            write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",this%step,".out"
        end if
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        if (present(gp2use)) then
           call decomp_2d_write_one(decomp_dir,arr,fname,gp2use)
        else
           call decomp_2d_write_one(decomp_dir,arr,fname,this%sp_gpC)
        end if

   end subroutine

   ! NOTE: If you want to dump an edge field, you need to call in dumpFullField
   ! routine with this%gpE passed in as the 3rd argument. If it's a cell field,
   ! then you don't need to pass in any gp since the default gp is this%gpC
   subroutine dumpSpectralField(this,arr,label,gp2use)
       use decomp_2d_io
       use mpi
       use exits, only: message
       class(igrid), intent(in) :: this
       character(len=clen) :: tempname, fname
       real(rkind), dimension(:,:,:), intent(in) :: arr
       character(len=4), intent(in) :: label
       type(decomp_info), intent(in), optional :: gp2use

        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",this%step,".out"
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        if (present(gp2use)) then
           call decomp_2d_write_one(2,arr,fname,gp2use)
        else
           call decomp_2d_write_one(2,arr,fname,this%gpC)
        end if

   end subroutine

   subroutine dumpVisualizationInfo(this)
       class(igrid), intent(in) :: this
       character(len=clen) :: tempname, fname

       
     if (nrank == 0) then
           write(tempname,"(a3,i2.2,a1,a4,a2,i6.6,a4)") "Run",this%runid, "_","info","_t",this%step,".out"
           fname = this%outputdir(:len_trim(this%outputdir))//"/"//trim(tempname)
           open(unit=10, file=trim(fname))
           write(10,"(100g17.9)") this%tsim
           write(10,"(100g17.9)") real(this%nx,rkind)
           write(10,"(100g17.9)") real(this%ny,rkind)
           write(10,"(100g17.9)") real(this%nz,rkind)
           close(10)

           if (this%localizedForceLayer == 2) then
               write(tempname,"(A3,I2.2,A1,A10,A2,I6.6,A4)") "Run",this%runID, "_","force_info","_t",this%step,".out"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               write(10,"(100g17.9)") this%spectForceLayer%ampFact_x
               write(10,"(100g17.9)") this%spectForceLayer%ampFact_y
               write(10,"(100g17.9)") this%spectForceLayer%ampFact_z
               close(10)
           end if
       end if 
   end subroutine

   subroutine dump_pointProbes(this)
       !use kind_parameters, only: mpirkind
       !class(igrid), intent(inout) :: this
       class(igrid), intent(in) :: this
       !character(len=clen) :: fname
       !character(len=clen) :: tempname
       !integer :: ierr

       ! ADITYA -> NIRANJAN: Why isn't this quantity inside turbarray? It
       ! segfaults for certain specific casses. 
       if(this%useWindTurbines) then
        !   this%runningSum_turb = zero
           !call MPI_reduce(this%inst_horz_avg_turb, this%runningSum_turb, 8*this%WindTurbineArr%nTurbines, mpirkind, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
           !if(nrank == 0) then
           !    write(tempname,"(A3,I2.2,A15)") "Run", this%RunID,"_timeseries.prb"
           !    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           !    open(unit=10,file=fname,status='old',action='write',position='append',iostat=ierr)
           !    if(ierr .ne. 0) open(unit=10,file=fname,status='replace')
           !    write(10,'(1000(e19.12,1x))') this%tsim, this%inst_horz_avg, this%runningSum_turb
           !    close(10)
           !end if
       endif

   end subroutine 

   subroutine dump_planes(this)
       use decomp_2d_io
       class(igrid), intent(inout) :: this
       integer :: nxplanes, nyplanes, nzplanes
       integer :: idx, pid, dirid, tid, sid
       character(len=clen) :: fname
       character(len=clen) :: tempname



       tid = this%step 
       if (allocated(this%xplanes)) then
           nxplanes = size(this%xplanes)
           dirid = 1
           do idx = 1,nxplanes
               pid = this%xplanes(idx)
           
               write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plu"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(1,this%u,dirid, pid, fname, this%gpC)

               write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plv"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(1,this%v,dirid, pid, fname, this%gpC)

               write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plw"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(1,this%wC,dirid, pid, fname, this%gpC)
               
               if (this%isStratified) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plT"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%T,dirid, pid, fname, this%gpC)
               end if

               if (this%fastCalcPressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plP"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure,dirid, pid, fname, this%gpC)
               end if 

               if (this%computeDNSPressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plD"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure_dns,dirid, pid, fname, this%gpC)

                   if (this%computeRapidSlowPressure) then
                       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".pRp"
                       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                       call decomp_2d_write_plane(1,this%prapid,dirid, pid, fname, this%gpC)
                       
                       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".pSl"
                       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                       call decomp_2d_write_plane(1,this%pslow,dirid, pid, fname, this%gpC)
                   end if 
               end if 

               if (this%computeturbinePressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".ptb"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure_turbine,dirid, pid, fname, this%gpC)
               end if 
               
               if (this%computeFringePressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".plF"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure_fringe,dirid, pid, fname, this%gpC)
               end if 


               if (this%computevorticity) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".pox"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%ox,dirid, pid, fname, this%gpC)
                   
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".poy"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%oy,dirid, pid, fname, this%gpC)
                   
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".poz"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%oz,dirid, pid, fname, this%gpC)
                   if (this%isStratified) then ! Compute potential vorticity
                       call this%compute_potential_vorticity(this%rbuffxC(:,:,:,1))
                       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A8)") "Run", this%RunID,"_t",&
                         tid,"_x",pid,".potvort"
                       call decomp_2d_write_plane(1,this%rbuffxC(:,:,:,1),dirid, pid, fname, this%gpC)
                   end if
               end if 

               if (this%WriteTurbineForce) then 
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".pTx"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%WindTurbineArr%fx,dirid, pid, fname, this%gpC)
                    
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".pTy"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%WindTurbineArr%fy,dirid, pid, fname, this%gpC)

                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".pTz"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%WindTurbineArr%fz,dirid, pid, fname, this%gpC)
               end if 

               ! planes for KS preprocess
               if (this%PreProcessForKS) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".ksu"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%uFil4KS,dirid, pid, fname, this%gpC)

                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".ksv"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%vFil4KS,dirid, pid, fname, this%gpC)

                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".ksw"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%wFil4KS,dirid, pid, fname, this%gpC)
               end if
              
               if (this%usescalars) then
                 do sid = 1,this%n_scalars
                    call this%scalars(sid)%dump_planes(tid, pid, dirid, "_x")
                 end do 
               end if 


           end do 
       end if 
           
           
       if (allocated(this%yplanes)) then
           nyplanes = size(this%yplanes)
           dirid = 2
           do idx = 1,nyplanes
               pid = this%yplanes(idx)

               write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plu"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(1,this%u,dirid, pid, fname, this%gpC)

               write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plv"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(1,this%v,dirid, pid, fname, this%gpC)

               write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plw"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(1,this%wC,dirid, pid, fname, this%gpC)
               
               if (this%isStratified) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plT"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%T,dirid, pid, fname, this%gpC)
               end if
               
               if (this%fastCalcPressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plP"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure,dirid, pid, fname, this%gpC)
               end if 

               if (this%computeDNSPressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plD"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure_dns,dirid, pid, fname, this%gpC)
                   
                   if (this%computeRapidSlowPressure) then
                       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".pRp"
                       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                       call decomp_2d_write_plane(1,this%prapid,dirid, pid, fname, this%gpC)
                       
                       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".pSl"
                       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                       call decomp_2d_write_plane(1,this%pslow,dirid, pid, fname, this%gpC)
                   end if 
               end if 
               
               if (this%computeturbinePressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".ptb"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure_turbine,dirid, pid, fname, this%gpC)
               end if 

               if (this%computeFringePressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".plF"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure_fringe,dirid, pid, fname, this%gpC)
               end if 
               
               if (this%computevorticity) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".pox"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%ox,dirid, pid, fname, this%gpC)
                   
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".poy"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%oy,dirid, pid, fname, this%gpC)
                   
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".poz"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%oz,dirid, pid, fname, this%gpC)
                   if (this%isStratified) then ! Compute potential vorticity
                       call this%compute_potential_vorticity(this%rbuffxC(:,:,:,1))
                       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A8)") "Run", this%RunID,"_t",&
                         tid,"_y",pid,".potvort"
                       call decomp_2d_write_plane(1,this%rbuffxC(:,:,:,1),dirid, pid, fname, this%gpC)
                   end if
               end if 
               
               if (this%WriteTurbineForce) then 
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".pTx"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%WindTurbineArr%fx,dirid, pid, fname, this%gpC)
                    
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".pTy"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%WindTurbineArr%fy,dirid, pid, fname, this%gpC)

                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".pTz"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%WindTurbineArr%fz,dirid, pid, fname, this%gpC)
               end if 
               
               ! planes for KS preprocess
               if (this%PreProcessForKS) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".ksu"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%uFil4KS,dirid, pid, fname, this%gpC)

                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".ksv"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%vFil4KS,dirid, pid, fname, this%gpC)

                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_y",pid,".ksw"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%wFil4KS,dirid, pid, fname, this%gpC)

               end if

               if (this%usescalars) then
                 do sid = 1,this%n_scalars
                    call this%scalars(sid)%dump_planes(tid, pid, dirid, "_y")
                 end do 
               end if 
           end do 
       end if 
       
       
       if (allocated(this%zplanes)) then
           nzplanes = size(this%zplanes)
           dirid = 3
           do idx = 1,nzplanes
               pid = this%zplanes(idx)

               write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plu"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(1,this%u,dirid, pid, fname, this%gpC)

               write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plv"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(1,this%v,dirid, pid, fname, this%gpC)

               write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plw"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(1,this%wC,dirid, pid, fname, this%gpC)
               
               if (this%isStratified) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plT"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%T,dirid, pid, fname, this%gpC)
               end if
               
               if (this%fastCalcPressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plP"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure,dirid, pid, fname, this%gpC)
               end if 

               if (this%computeDNSPressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plD"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure_dns,dirid, pid, fname, this%gpC)
                   
                   if (this%computeRapidSlowPressure) then
                       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".pRp"
                       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                       call decomp_2d_write_plane(1,this%prapid,dirid, pid, fname, this%gpC)
                       
                       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".pSl"
                       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                       call decomp_2d_write_plane(1,this%pslow,dirid, pid, fname, this%gpC)
                   end if 
               end if 
               
               if (this%computeturbinePressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_xz",pid,".ptb"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure_turbine,dirid, pid, fname, this%gpC)
               end if 

               if (this%computeFringePressure) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".plF"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%Pressure_fringe,dirid, pid, fname, this%gpC)
               end if 
               
               if (this%computevorticity) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".pox"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%ox,dirid, pid, fname, this%gpC)
                   
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".poy"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%oy,dirid, pid, fname, this%gpC)
                   
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_x",pid,".poz"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%oz,dirid, pid, fname, this%gpC)
                   if (this%isStratified) then ! Compute potential vorticity
                       call this%compute_potential_vorticity(this%rbuffxC(:,:,:,1))
                       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A8)") "Run", this%RunID,"_t",&
                         tid,"_z",pid,".potvort"
                       call decomp_2d_write_plane(1,this%rbuffxC(:,:,:,1),dirid, pid, fname, this%gpC)
                   end if
               end if 

               if (this%WriteTurbineForce) then 
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".pTx"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%WindTurbineArr%fx,dirid, pid, fname, this%gpC)
                    
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".pTy"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%WindTurbineArr%fy,dirid, pid, fname, this%gpC)

                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".pTz"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%WindTurbineArr%fz,dirid, pid, fname, this%gpC)
               end if 
               
               ! planes for KS preprocess
               if (this%PreProcessForKS) then
                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".ksu"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%uFil4KS,dirid, pid, fname, this%gpC)

                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".ksv"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%vFil4KS,dirid, pid, fname, this%gpC)

                   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A4)") "Run", this%RunID,"_t",tid,"_z",pid,".ksw"
                   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                   call decomp_2d_write_plane(1,this%wFil4KS,dirid, pid, fname, this%gpC)
               end if
               
               if (this%usescalars) then
                 do sid = 1,this%n_scalars
                    call this%scalars(sid)%dump_planes(tid, pid, dirid, "_z")
                 end do 
               end if 
           end do 
       end if 
       call message(1, "Dumped Planes.")        
   end subroutine 

   
   subroutine finalize_io(this)
       class(igrid), intent(in) :: this

       if (nrank == 0) then
           write(this%headerfid,*) "--------------------------------------------------------------"
           write(this%headerfid,*) "------------------ END OF HEADER FILE ------------------------"
           close(this%headerfid)
       end if 
   end subroutine

   subroutine readRestartFile(this, tid, rid, u, v, w, T, gpC, gpE, tid_out)
       use mpi
       use exits, only: message
       use kind_parameters, only: mpirkind
       class(igrid), intent(inout) :: this
       integer, intent(in) :: tid, rid
       real(rkind), dimension(:,:,:), intent(inout) :: u, v, w, T
       class(decomp_info), intent(in) :: gpC, gpE
       integer, intent(out) :: tid_out
       character(len=clen) :: tempname, fname
       integer :: ierr, fid, idx, ios 
       real(rkind) :: tmp

       if (tid < 0) then
           write(tempname,"(A7,A4,I2.2,A9)") "RESTART", "_Run",rid, "_u.LATEST"
       else
           write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_u.",tid
       end if
       fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
       call check_file_existence(trim(fname))
       call decomp_2d_read_one(1,u,fname, gpC)

       if (tid < 0) then
           write(tempname,"(A7,A4,I2.2,A9)") "RESTART", "_Run",rid, "_v.LATEST"
       else
           write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_v.",tid
       end if
       fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
       call check_file_existence(trim(fname))
       call decomp_2d_read_one(1,v,fname, gpC)

       if (tid < 0) then
           write(tempname,"(A7,A4,I2.2,A9)") "RESTART", "_Run",rid, "_w.LATEST"
       else
           write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_w.",tid
       end if
       fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
       call check_file_existence(trim(fname))
       call decomp_2d_read_one(1,w,fname, gpE)

       if (this%isStratified) then
           if (tid < 0) then
               write(tempname,"(A7,A4,I2.2,A9)") "RESTART", "_Run",rid, "_T.LATEST"
           else
               write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_T.",tid
           end if
           fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
           call check_file_existence(trim(fname))
           call decomp_2d_read_one(1,T,fname, gpC)
       end if

       if (nrank == 0) then
           if (tid < 0) then
               write(tempname,"(A7,A4,I2.2,A12)") "RESTART", "_Run",rid, "_info.LATEST"
           else
               write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",rid, "_info.",tid
           end if
           fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
           call check_file_existence(trim(fname))
           fid = 10
           open(unit=fid,file=trim(fname),status="old",action="read")
           read (fid, *, iostat=ios)  this%tsim
           if (ios /= 0) then
               print*, 'Error reading simulation time: ',ios
               close(fid)
               stop
           end if
           read (fid, *, iostat=ios)  tid_out
           if (ios /= 0) then
               print*, 'Error reading simulation time ID: ',ios
               close(fid)
               stop
           end if
           if (this%useControl) then
              read(fid,"(100g15.5)") this%restartPhi
           end if
           close(fid)
       end if

       call MPI_BCast(tid_out,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       call mpi_barrier(mpi_comm_world, ierr)
       call mpi_bcast(this%tsim,1,mpirkind,0,mpi_comm_world,ierr)
       call mpi_barrier(mpi_comm_world, ierr)
       call message("================= RESTART FILE USED ======================")
       call message(0, "Simulation Time at restart:", this%tsim)
       call message_min_max(0,"u bounds:",p_minval(minval(u)),p_maxval(maxval(u)))
       call message_min_max(0,"v bounds:",p_minval(minval(v)),p_maxval(maxval(v)))
       call message_min_max(0,"w bounds:",p_minval(minval(w)),p_maxval(maxval(w)))
       call message("=================================== ======================")

   end subroutine

   subroutine readVizForRestart(this, tid, rid, u, v, w, wC, T, gpC)
       use decomp_2d_io
       use mpi
       use exits, only: message
       use kind_parameters, only: mpirkind
       class(igrid), intent(inout) :: this
       integer, intent(in) :: tid, rid
       real(rkind), dimension(:,:,:), intent(inout) :: u, v, w, wC, T
       class(decomp_info), intent(in) :: gpC
       character(len=clen) :: tempname, fname
       integer :: ierr, fid 

       write(tempname,"(A3,I2.2,A7,I6.6,A4)") "Run",rid, "_uVel_t",tid,".out"
       fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
       call check_file_existence(trim(fname))
       call decomp_2d_read_one(1,u,fname, gpC)

       write(tempname,"(A3,I2.2,A7,I6.6,A4)") "Run",rid, "_vVel_t",tid,".out"
       fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
       call check_file_existence(trim(fname))
       call decomp_2d_read_one(1,v,fname, gpC)

       write(tempname,"(A3,I2.2,A7,I6.6,A4)") "Run",rid, "_wVel_t",tid,".out"
       fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
       call check_file_existence(trim(fname))
       call decomp_2d_read_one(1,wC,fname, gpC)
       call transpose_x_to_y(this%wC, this%rbuffyC(:,:,:,1), gpC)
       call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), gpC)
       call this%Pade6opZ%interpz_C2E(this%rbuffzC(:,:,:,1), this%wC, wBC_bottom, wBC_top)

       if (this%isStratified) then
           write(tempname,"(A3,I2.2,A7,I6.6,A4)") "Run",rid,"_potT_t",tid,".out"
           fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
           call check_file_existence(trim(fname))
           call decomp_2d_read_one(1,T,fname, gpC)
       end if
       
       if (nrank == 0) then
         write(tempname,"(A3,I2.2,A7,I6.6,A4)") "Run",rid, "_info_t",tid,".out"
         fname = this%inputdir(:len_trim(this%inputdir))//"/"//trim(tempname)
         call check_file_existence(trim(fname))
         open(unit=10,file=fname,access='sequential',form='formatted')
         read (10, *)  this%tsim
         close(10)
         call message(0, "Reading visualization data dumped at tSIM=", this%tsim)
       end if 

       call mpi_barrier(mpi_comm_world, ierr)
       call mpi_bcast(this%tsim,1,mpirkind,0,mpi_comm_world,ierr)
       call mpi_barrier(mpi_comm_world, ierr)
       call message("============== RESTARTING FROM VISUALIZATION  ============")
       call message(0, "Simulation Time at restart:", this%tsim)
       call message_min_max(0,"u bounds:",p_minval(minval(u)),p_maxval(maxval(u)))
       call message_min_max(0,"v bounds:",p_minval(minval(v)),p_maxval(maxval(v)))
       call message_min_max(0,"w bounds:",p_minval(minval(w)),p_maxval(maxval(w)))
       call message("=================================== ======================")

   end subroutine

   !subroutine initialize_hdf5_io(this)
   !    class(igrid), intent(inout) :: this
   !    character(len=5) :: filename_prefix

   !    write(filename_prefix,"(A3,I2.2)") "Run", this%runID
   !    if (this%ioType == 1) then
   !        call this%viz_hdf5%init(MPI_COMM_WORLD,this%gpC, "x",this%outputdir, filename_prefix, &
   !           reduce_precision=.false.,write_xdmf=.true., read_only=.false., wider_time_format=.true.) 
   !    elseif (this%ioType == 2) then
   !        call this%viz_hdf5%init(MPI_COMM_WORLD,this%gpC, "x",this%outputdir, filename_prefix, &
   !           reduce_precision=.true.,write_xdmf=.true., read_only=.false., wider_time_format=.true.)  
   !    else
   !        call GracefulExit("Invalid choice for IOTYPE", 312)
   !    end if 

   !    call this%viz_hdf5%write_coords(this%mesh)
   !end subroutine 

   !subroutine destroy_hdf5_io(this)
   !    class(igrid), intent(inout) :: this

   !    call this%viz_hdf5%destroy()
   !end subroutine 

   subroutine append_visualization_info(this)
       class(igrid), intent(in) :: this 
       character(len=clen) :: tempname, fname 
       logical :: exists 

       write(tempname,"(A3,I2.2,A12,A4)") "Run",this%runID, "_vis_summary",".smm"
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)

       if (nrank == 0) then
           inquire(file=fname, exist=exists)
           if (exists) then
               open(12, file=fname, status="old", position="append", action="write")
           else
               open(12, file=fname, status="new", action="write")
           end if
           write(12,*) this%tsim, this%step
           close(12)
       end if 
   end subroutine 

   subroutine dump_scalar_fields(this)
     class(igrid), intent(inout) :: this
     integer :: idx

     if ((this%usescalars) .and. allocated(this%scalars)) then
        do idx = 1,this%n_scalars
           if (this%ioType == 0) then
              call this%scalars(idx)%dumpScalarField(this%step,this%dump_KAPPA_SGS)
           else
              !call this%scalars(idx)%dumpScalarField(this%step, this%viz_hdf5)
           end if 
        end do 
     end if 
   end subroutine 

   subroutine dump_visualization_files(this,ddt,step)
       use exits, only: message_min_max, message
       class(igrid), intent(inout) :: this
       logical, intent(in), optional :: ddt
       integer, intent(in), optional :: step
       integer :: tstep
       logical :: dump_ddt_terms

       dump_ddt_terms = .false.
       tstep          = this%step
       if (present(ddt)) dump_ddt_terms = ddt
       if (present(step)) tstep = step

       select case (this%ioType) 
       case(0)
           if (dump_ddt_terms) then
               call this%dumpFullField(this%dudt,'dudt',step=tstep)
               call this%dumpFullField(this%dvdt,'dvdt',step=tstep)
               call this%dumpFullField(this%dwdt,'dwdt',step=tstep)
               if (this%isStratified) then
                   call this%dumpFullField(this%dTdt,'dTdt')
               end if
           else
               call this%dumpFullField(this%u,'uVel')
               call this%dumpFullField(this%v,'vVel')
               call this%dumpFullField(this%wC,'wVel')
               call this%dump_scalar_fields()
               call this%dumpVisualizationInfo()
               if (this%isStratified .or. this%initspinup) call this%dumpFullField(this%T,'potT')
               if (this%fastCalcPressure) call this%dumpFullField(this%pressure,'prss')
               if (this%computeDNSpressure) call this%dumpFullField(this%pressure_dns,'pdns')
               if (this%computeturbinepressure) call this%dumpFullField(this%pressure_turbine,'ptrb')
               if (this%computefringepressure) call this%dumpFullField(this%pressure_fringe,'pfrn')
               if ((this%useSGS) .and. (this%dump_NU_SGS)) call this%dumpFullField(this%nu_SGS,'nSGS')
               if ((this%useSGS) .and. (this%dump_KAPPA_SGS) .and. (this%isStratified)) call this%dumpFullField(this%kappaSGS,'kSGS')
               if ((this%useSGS) .and. (this%dump_KAPPA_SGS) .and. (this%isStratified) .and. associated(this%kappa_bounding)) then
                  if (this%augment_SGS_with_scalar_bounding .or. this%use_scalar_bounding_as_SGS) then
                      call this%dumpFullField(this%kappa_bounding,'kBND')
                  else
                      this%rbuffxC(:,:,:,1) = this%kappa_bounding*this%kappa_bounding_mask
                      call this%dumpFullField(this%rbuffxC(:,:,:,1),'kBND')
                  end if
               end if 
               if (this%computeRapidSlowPressure) then
                   call this%dumpFullField(this%prapid,'prap')
                   call this%dumpFullField(this%pslow,'pslo')
               end if 
               if (this%computevorticity) then
                   call this%dumpFullField(this%ox,'omgX')
                   call this%dumpFullField(this%oy,'omgY')
                   call this%dumpFullField(this%oz,'omgZ')
                   if (this%isStratified) then ! dump potential vorticity
                       call this%compute_potential_vorticity(this%rbuffxC(:,:,:,1))
                       call this%dumpFullField(this%rbuffxC(:,:,:,1),'pVrt')
                   end if
               end if
               if (this%WriteTurbineForce) then 
                    call this%dumpFullField(this%WindTurbineArr%fx, "TrbX")
                    call this%dumpFullField(this%WindTurbineArr%fy, "TrbY")
                    call this%dumpFullField(this%WindTurbineArr%fz, "TrbZ")
               end if

               if (this%localizedForceLayer == 1) then
                   if (this%forceLayer%dumpForce) then
                       this%rbuffxC(:,:,:,1) = this%forceLayer%ampfact*this%forceLayer%fx
                       this%rbuffxC(:,:,:,2) = this%forceLayer%ampfact*this%forceLayer%fy
                       this%rbuffxE(:,:,:,1) = this%forceLayer%ampfact*this%forceLayer%fz
                       call this%dumpFullField(this%rbuffxC(:,:,:,1), "frcx")
                       call this%dumpFullField(this%rbuffxC(:,:,:,2), "frcy")
                       call this%forceLayer%interpE2C(this%rbuffxE(:,:,:,1),&
                         this%rbuffxC(:,:,:,3), this%rbuffyC(:,:,:,1), &
                         this%rbuffzC(:,:,:,1), this%rbuffyE(:,:,:,1), &
                         this%rbuffzE(:,:,:,1))
                       call this%dumpFullField(this%rbuffxC(:,:,:,3), "frcz")
                       if (this%forceLayer%dumpSplines) then
                           call this%dumpFullField(this%forceLayer%phixC, "phxC", this%gpC)
                           call this%dumpFullField(this%forceLayer%phiyC, "phyC", this%gpC)
                           call this%dumpFullField(this%forceLayer%phizC, "phzC", this%gpC)
                           
                           call this%dumpFullField(this%forceLayer%dphixC, "dpxC", this%gpC)
                           call this%dumpFullField(this%forceLayer%dphiyC, "dpyC", this%gpC)
                           call this%dumpFullField(this%forceLayer%dphizC, "dpzC", this%gpC)
                           
                           call this%dumpFullField(this%forceLayer%phixE, "phxE", this%gpE)
                           call this%dumpFullField(this%forceLayer%phiyE, "phyE", this%gpE)
                           call this%dumpFullField(this%forceLayer%phizE, "phzE", this%gpE)
                           
                           call this%dumpFullField(this%forceLayer%dphixE, "dpxE", this%gpE)
                           call this%dumpFullField(this%forceLayer%dphiyE, "dpyE", this%gpE)
                           call this%dumpFullField(this%forceLayer%dphizE, "dpzE", this%gpE)
                       end if
                   end if
               elseif (this%localizedForceLayer == 2) then
                   if (this%spectForceLayer%dumpForce) then
                       call this%dumpFullField(this%spectForceLayer%ampFact_x*this%spectForceLayer%fx, "frcx")
                       call this%dumpFullField(this%spectForceLayer%ampFact_y*this%spectForceLayer%fy, "frcy")
                       call this%spectForceLayer%interpE2C(this%spectForceLayer%fz,&
                         this%rbuffxC(:,:,:,3), this%rbuffyC(:,:,:,1), &
                         this%rbuffzC(:,:,:,1), this%rbuffyE(:,:,:,1), &
                         this%rbuffzE(:,:,:,1))
                       call this%dumpFullField(this%spectForceLayer%ampFact_z*this%rbuffxC(:,:,:,3), "frcz")
                       if (this%isStratified) then
                           call this%dumpFullField(this%spectForceLayer%fT, "frcT")
                       end if

                   end if
               end if
           end if
       case default
           !call this%viz_hdf5%update_vizcount(this%step)
           !call this%viz_hdf5%start_viz(this%tsim)
           !call this%viz_hdf5%write_variable(this%u, "uVel")
           !call this%viz_hdf5%write_variable(this%v, "vVel")
           !call this%viz_hdf5%write_variable(this%wC, "wVel")
           !call this%dump_scalar_fields()           
           !if (this%isStratified .or. this%initspinup) call this%viz_hdf5%write_variable(this%T,'potT')
           !if (this%fastCalcPressure) call this%viz_hdf5%write_variable(this%pressure,'prss')
           !if (this%computeDNSpressure) call this%viz_hdf5%write_variable(this%pressure_dns,'pdns')
           !if (this%computeturbinepressure) call this%viz_hdf5%write_variable(this%pressure_turbine,'ptrb')
           !if (this%computefringepressure) call this%viz_hdf5%write_variable(this%pressure_fringe,'pfrn')
           !if ((this%useSGS) .and. (this%dump_NU_SGS)) call this%viz_hdf5%write_variable(this%nu_SGS,'nSGS')
           !if ((this%useSGS) .and. (this%dump_KAPPA_SGS) .and. (this%isStratified)) call this%viz_hdf5%write_variable(this%kappaSGS,'kSGS')
           !if ((this%useSGS) .and. (this%dump_KAPPA_SGS) .and. (this%isStratified) .and. associated(this%kappa_bounding)) then
           !   call this%viz_hdf5%write_variable(this%kappa_bounding,'kBND')
           !end if 
           !if (this%computeRapidSlowPressure) then
           !    call this%viz_hdf5%write_variable(this%prapid,'prap')
           !    call this%viz_hdf5%write_variable(this%pslow,'pslo')
           !end if 
           !if (this%computevorticity) then
           !    call this%viz_hdf5%write_variable(this%ox,'omgX')
           !    call this%viz_hdf5%write_variable(this%oy,'omgY')
           !    call this%viz_hdf5%write_variable(this%oz,'omgZ')
           !end if 
           !if (this%WriteTurbineForce) then 
           !     call this%viz_hdf5%write_variable(this%WindTurbineArr%fx, "TrbX")
           !     call this%viz_hdf5%write_variable(this%WindTurbineArr%fy, "TrbY")
           !     call this%viz_hdf5%write_variable(this%WindTurbineArr%fz, "TrbZ")
           !end if

           !call this%viz_hdf5%end_viz()
       end select 
       call this%append_visualization_info()
       
       if (this%PreProcessForKS) then
            if (.not. this%KSupdated) then
                call this%LES2KS%applyFilterForKS(this%u, this%v, this%wC)
                this%KSupdated = .true. 
            end if
            call this%dumpFullField(this%uFil4KS,'uFks')
            call this%dumpFullField(this%vFil4KS,'vFks')
            call this%dumpFullField(this%wFil4KS,'wFks')
       end if
       if (this%useWindTurbines) then
           call this%WindTurbineArr%write_turbine_power(this%step, this%outputdir, this%runID)
       end if
   end subroutine 
   
   subroutine start_io(this, dumpInitField)
       class(igrid), target, intent(inout) :: this 
       character(len=clen) :: fname
       character(len=clen) :: tempname
       !character(len=clen) :: command
       character(len=clen) :: OutputDir
       !integer :: system 
       integer :: runIDX 
       logical :: isThere
       integer :: tag, idx, status(MPI_STATUS_SIZE), ierr
       integer, dimension(:,:), allocatable        :: xst,xen,xsz
       logical, optional, intent(in) :: dumpInitField

       ! Create data sharing info
       !if (nrank == 0) then
           allocate(xst(0:nproc-1,3),xen(0:nproc-1,3),xsz(0:nproc-1,3))
           xst = 0; xen = 0; xsz = 0;
       !end if


       ! communicate local processor grid info (Assume x-decomposition)
       if (nrank == 0) then
           xst(0,:) = this%gpC%xst
           xen(0,:) = this%gpC%xen
           
           tag = 0
           do idx = 1,nproc-1
               call MPI_RECV(xst(idx,:), 3, MPI_INTEGER, idx, tag,&
                   MPI_COMM_WORLD, status, ierr)
           end do 
           tag = 1
           do idx = 1,nproc-1
               call MPI_RECV(xen(idx,:), 3, MPI_INTEGER, idx, tag,&
                   MPI_COMM_WORLD, status, ierr)
           end do
          tag = 2
           do idx = 1,nproc-1
               call MPI_RECV(xsz(idx,:), 3, MPI_INTEGER, idx, tag,&
                   MPI_COMM_WORLD, status, ierr)
           end do

       else
           tag = 0
           call MPI_SEND(this%gpC%xst, 3, MPI_INTEGER, 0, tag, &
                &      MPI_COMM_WORLD, ierr)
           tag = 1
           call MPI_SEND(this%gpC%xen, 3, MPI_INTEGER, 0, tag, &
                &      MPI_COMM_WORLD, ierr)
           tag = 2
           call MPI_SEND(this%gpC%xsz, 3, MPI_INTEGER, 0, tag, &
                &      MPI_COMM_WORLD, ierr)

       end if 

       OutputDir = this%outputdir
       runIDX = this%runID
       
       inquire(FILE=trim(OutputDir), exist=isThere)
       if (nrank == 0) then
           write(tempname,"(A3,I2.2,A6,A4)") "Run", runIDX, "HEADER",".txt"
           fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
           open (this%headerfid, file=trim(fname), FORM='formatted', STATUS='replace',ACTION='write')
           write(this%headerfid,*)"========================================================================="
           write(this%headerfid,*)"---------------------  Header file for MATLAB ---------------------------"
           write(this%headerfid,"(A9,A10,A10,A10,A10,A10,A10)") "PROC", "xst", "xen", "yst", "yen","zst","zen"
           write(this%headerfid,*)"-------------------------------------------------------------------------"
           do idx = 0,nproc-1
               write(this%headerfid,"(I8,6I10)") idx, xst(idx,1), xen(idx,1), xst(idx,2), xen(idx,2), xst(idx,3), xen(idx,3)
           end do 
           write(this%headerfid,*)"-------------------------------------------------------------------------"
           write(this%headerfid,*)"Dumps made at:"
       end if
       call mpi_barrier(mpi_comm_world,ierr)
       
       !if (nrank == 0) then
           deallocate(xst, xen, xsz)
       !end if 

       if (present(dumpInitField)) then
           if (dumpInitField) then
               call message(0,"Performing initialization data dump.")
               call this%dump_visualization_files()
               call message(0,"Done with the initialization data dump.")
           end if
       end if
       call message(0,"Done with io_startup stuff")
   end subroutine
   
   subroutine readField3D(RunID, TIDX, inputDir, label, field, gpC)
       use exits, only: GracefulExit
       use decomp_2d_io
       
       integer, intent(in) :: RunID, TIDX
       character(len=4), intent(in) :: label
       character(len=*), intent(in) :: inputDir
       type(decomp_info), intent(in) :: gpC
       real(rkind), dimension(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)), intent(out) :: field
       character(len=clen) :: tempname, fname

       write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",RunID, "_",label,"_t",TIDX,".out"
       fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
       
       open(777,file=trim(fname),status='old',iostat=ierr)
       if(ierr .ne. 0) then
        print*, "Rank:", nrank, ". File:", fname, " not found"
        call mpi_barrier(mpi_comm_world, ierr)
        call gracefulExit("File I/O issue.",44)
       end if 
       close(777)

       call decomp_2d_read_one(1,field,fname,gpC)

   end subroutine 


   subroutine dumpProbes(this)
       use basic_io, only: write_2d_ascii, append_2D_ascii
       class(igrid), intent(in) :: this
       character(len=clen) :: tempname, fname
       integer :: pid, idx
       logical :: exists

       do idx = 1,this%nprobes
           pid = this%probes(4,idx)
           !write(tempname,"(A3,I2.2,A6,I3.3,A4,I6.6,A4,I6.6,A4)") "Run",this%runID, "_PROBE",pid,"_tst",this%probeStartStep,"_ten",this%step,".out"
           !fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           !call write_2d_ascii(transpose(this%probe_data(:,idx,this%probeStartStep:this%step)), fname)

           write(tempname,"(A3,I2.2,A6,I3.3,A4,I6.6,A4,I6.6,A4)") "Run",this%runID, "_PROBE",pid,"_tst",this%probeStartStep,".out"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           inquire(file=trim(fname),exist=exists) 
           if (exists) then
             call append_2D_ascii(transpose(this%probe_data(:,:,idx)),trim(fname))
           else
             call write_2D_ascii(transpose(this%probe_data(:,:,idx)),trim(fname))
           end if

       end do 

       ! KS - preprocess
       if (this%PreprocessForKS) then
           do idx = 1,this%nprobes
               pid = this%probes(4,idx)
               write(tempname,"(A3,I2.2,A9,I3.3,A4,I6.6,A4,I6.6,A4)") "Run",this%runID, "_PROBE_KS",pid,"_tst",this%probeStartStep,"_ten",this%step,".out"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call write_2d_ascii(transpose(this%KS_probe_data(:,idx,this%probeStartStep:this%step)), fname)
           end do 
       end if
       
   end subroutine

   function vizFileExists(this,varname,tid) result(res)
       class(igrid), intent(in) :: this
       character(len=4), intent(in) :: varname
       integer, intent(in) :: tid
       character(len=clen) :: fname
       logical :: res

       write(fname,"(A,I2.2,A1,A4,A2,I6.6,A4)")&
         trim(this%inputdir)//'/Run',this%runID,'_',varname,'_t',tid,'.out'
       inquire( file=trim(fname), exist=res )
   end function

   subroutine readVizFile(this,varname,tid,field)
       class(igrid), intent(in) :: this
       character(len=4), intent(in) :: varname
       integer, intent(in) :: tid
       real(rkind), dimension(:,:,:), intent(inout) :: field
       character(len=clen) :: fname

       write(fname,"(A,I2.2,A1,A4,A2,I6.6,A4)")&
         trim(this%inputdir)//'/Run',this%runID,'_',varname,'_t',tid,'.out'
       call decomp_2d_read_one(1,field,trim(fname),this%gpC)

   end subroutine
