   subroutine dumpRestartFile(this)
       use decomp_2d_io
       use mpi
       use exits, only: message
       class(igrid), intent(in) :: this
       character(len=clen) :: tempname, fname
       integer :: ierr, idx 

       write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_u.",this%step
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(1,this%u,fname, this%gpC)

       write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_v.",this%step
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(1,this%v,fname, this%gpC)

       write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_w.",this%step
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(1,this%w,fname, this%gpE)

       if (this%isStratified .or. this%initspinup) then
           write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",this%runID, "_T.",this%step
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(1,this%T,fname, this%gpC)
       end if 

       if (this%usescalars) then
           do idx = 1,this%n_scalars
              call this%scalars(idx)%dumprestart(this%step)
           end do
       end if

       if (nrank == 0) then
           write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",this%runID, "_info.",this%step
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           OPEN(UNIT=10, FILE=trim(fname))
           write(10,"(100g15.5)") this%tsim
           if (this%useControl) then
              write(10,"(100g15.5)") this%G_alpha
           end if
           close(10)
       end if 

       call mpi_barrier(mpi_comm_world, ierr)
       call message(1, "Just Dumped a RESTART file")

   end subroutine

   ! NOTE: If you want to dump an edge field, you need to call in dumpFullField
   ! routine with this%gpE passed in as the 3rd argument. If it's a cell field,
   ! then you don't need to pass in any gp since the default gp is this%gpC
   subroutine dumpFullField(this,arr,label,gp2use)
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
           call decomp_2d_write_one(1,arr,fname,gp2use)
        else
           call decomp_2d_write_one(1,arr,fname,this%gpC)
        end if

   end subroutine


   subroutine dumpVisualizationInfo(this)
       class(igrid), intent(in) :: this
       character(len=clen) :: tempname, fname

       
     if (nrank == 0) then
           write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_","info","_t",this%step,".out"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           OPEN(UNIT=10, FILE=trim(fname))
           write(10,"(100g17.9)") this%tsim
           write(10,"(100g17.9)") real(this%nx,rkind)
           write(10,"(100g17.9)") real(this%ny,rkind)
           write(10,"(100g17.9)") real(this%nz,rkind)
           close(10)
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
       class(igrid), intent(in) :: this
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

   subroutine readRestartFile(this, tid, rid)
       use decomp_2d_io
       use mpi
       use exits, only: message
       use kind_parameters, only: mpirkind
       class(igrid), intent(inout) :: this
       integer, intent(in) :: tid, rid
       character(len=clen) :: tempname, fname
       integer :: ierr, fid 

       write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_u.",tid
       fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
       call decomp_2d_read_one(1,this%u,fname, this%gpC)

       write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_v.",tid
       fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
       call decomp_2d_read_one(1,this%v,fname, this%gpC)

       write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_w.",tid
       fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
       call decomp_2d_read_one(1,this%w,fname, this%gpE)

       if (this%isStratified) then
           write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_T.",tid
           fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
           call decomp_2d_read_one(1,this%T,fname, this%gpC)
       end if 

       if (nrank == 0) then
           write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",rid, "_info.",tid
           fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
           fid = 10
           open(unit=fid,file=trim(fname),status="old",action="read")
           read (fid, "(100g15.5)")  this%tsim
           if (this%useControl) then
              read(fid,"(100g15.5)") this%restartPhi
           end if
           close(fid)
       end if 

       call mpi_barrier(mpi_comm_world, ierr)
       call mpi_bcast(this%tsim,1,mpirkind,0,mpi_comm_world,ierr)
       call mpi_barrier(mpi_comm_world, ierr)
       call message("================= RESTART FILE USED ======================")
       call message(0, "Simulation Time at restart:", this%tsim)
       call message("=================================== ======================")

   end subroutine

   subroutine initialize_hdf5_io(this)
       class(igrid), intent(inout) :: this
       character(len=5) :: filename_prefix

       write(filename_prefix,"(A3,I2.2)") "Run", this%runID
       if (this%ioType == 1) then
           call this%viz_hdf5%init(MPI_COMM_WORLD,this%gpC, "x",this%outputdir, filename_prefix, &
              reduce_precision=.false.,write_xdmf=.true., read_only=.false., wider_time_format=.true.) 
       elseif (this%ioType == 2) then
           call this%viz_hdf5%init(MPI_COMM_WORLD,this%gpC, "x",this%outputdir, filename_prefix, &
              reduce_precision=.true.,write_xdmf=.true., read_only=.false., wider_time_format=.true.)  
       else
           call GracefulExit("Invalid choice for IOTYPE", 312)
       end if 

       call this%viz_hdf5%write_coords(this%mesh)
   end subroutine 

   subroutine destroy_hdf5_io(this)
       class(igrid), intent(inout) :: this

       call this%viz_hdf5%destroy()
   end subroutine 

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
              call this%scalars(idx)%dumpScalarField(this%step)
           else
              call this%scalars(idx)%dumpScalarField(this%step, this%viz_hdf5)
           end if 
        end do 
     end if 
   end subroutine 

   subroutine dump_visualization_files(this)
       class(igrid), intent(inout) :: this

       select case (this%ioType) 
       case(0)
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
              call this%dumpFullField(this%kappa_bounding,'kBND')
           end if 
           if (this%computeRapidSlowPressure) then
               call this%dumpFullField(this%prapid,'prap')
               call this%dumpFullField(this%pslow,'pslo')
           end if 
           if (this%computevorticity) then
               call this%dumpFullField(this%ox,'omgX')
               call this%dumpFullField(this%oy,'omgY')
               call this%dumpFullField(this%oz,'omgZ')
           end if
           if (this%WriteTurbineForce) then 
                call this%dumpFullField(this%WindTurbineArr%fx, "TrbX")
                call this%dumpFullField(this%WindTurbineArr%fy, "TrbY")
                call this%dumpFullField(this%WindTurbineArr%fz, "TrbZ")
           end if 
       case default
           call this%viz_hdf5%update_vizcount(this%step)
           call this%viz_hdf5%start_viz(this%tsim)
           call this%viz_hdf5%write_variable(this%u, "uVel")
           call this%viz_hdf5%write_variable(this%v, "vVel")
           call this%viz_hdf5%write_variable(this%wC, "wVel")
           call this%dump_scalar_fields()           
           if (this%isStratified .or. this%initspinup) call this%viz_hdf5%write_variable(this%T,'potT')
           if (this%fastCalcPressure) call this%viz_hdf5%write_variable(this%pressure,'prss')
           if (this%computeDNSpressure) call this%viz_hdf5%write_variable(this%pressure_dns,'pdns')
           if (this%computeturbinepressure) call this%viz_hdf5%write_variable(this%pressure_turbine,'ptrb')
           if (this%computefringepressure) call this%viz_hdf5%write_variable(this%pressure_fringe,'pfrn')
           if ((this%useSGS) .and. (this%dump_NU_SGS)) call this%viz_hdf5%write_variable(this%nu_SGS,'nSGS')
           if ((this%useSGS) .and. (this%dump_KAPPA_SGS) .and. (this%isStratified)) call this%viz_hdf5%write_variable(this%kappaSGS,'kSGS')
           if ((this%useSGS) .and. (this%dump_KAPPA_SGS) .and. (this%isStratified) .and. associated(this%kappa_bounding)) then
              call this%viz_hdf5%write_variable(this%kappa_bounding,'kBND')
           end if 
           if (this%computeRapidSlowPressure) then
               call this%viz_hdf5%write_variable(this%prapid,'prap')
               call this%viz_hdf5%write_variable(this%pslow,'pslo')
           end if 
           if (this%computevorticity) then
               call this%viz_hdf5%write_variable(this%ox,'omgX')
               call this%viz_hdf5%write_variable(this%oy,'omgY')
               call this%viz_hdf5%write_variable(this%oz,'omgZ')
           end if 
           if (this%WriteTurbineForce) then 
                call this%viz_hdf5%write_variable(this%WindTurbineArr%fx, "TrbX")
                call this%viz_hdf5%write_variable(this%WindTurbineArr%fy, "TrbY")
                call this%viz_hdf5%write_variable(this%WindTurbineArr%fz, "TrbZ")
           end if

           call this%viz_hdf5%end_viz()
       end select 
       if (this%useProbes) then
            call this%dumpProbes()    
            call message(0,"Performed a scheduled dump for probes.")
       end if
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
               !call this%dumpFullField(this%u,'uVel')
               !call this%dumpFullField(this%v,'vVel')
               !call this%dumpFullField(this%wC,'wVel')
               !call this%dump_scalar_fields()
               !call this%dumpVisualizationInfo()
               !if (this%isStratified .or. this%initspinup) call this%dumpFullField(this%T,'potT')
               !if (this%fastCalcPressure) call this%dumpFullField(this%pressure,'prss')
               !if (this%computeDNSpressure) call this%dumpFullField(this%pressure_dns,'pdns')
               !if (this%computeturbinepressure) call this%dumpFullField(this%pressure_turbine,'ptrn')
               !if (this%computefringepressure) call this%dumpFullField(this%pressure_fringe,'pfrn')
               !if (this%useWindTurbines) then
               !    this%WindTurbineArr%dumpTurbField = .true.
               !    this%WindTurbineArr%step = this%step-1
               !endif
               call this%dump_visualization_files()
               call message(0,"Done with the initialization data dump.")
           end if
       end if
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
       use basic_io, only: write_2d_ascii
       class(igrid), intent(in) :: this
       character(len=clen) :: tempname, fname
       integer :: pid, idx

       do idx = 1,this%nprobes
           pid = this%probes(4,idx)
           write(tempname,"(A3,I2.2,A6,I3.3,A4,I6.6,A4,I6.6,A4)") "Run",this%runID, "_PROBE",pid,"_tst",this%probeStartStep,"_ten",this%step,".out"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call write_2d_ascii(transpose(this%probe_data(:,idx,this%probeStartStep:this%step)), fname)
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
