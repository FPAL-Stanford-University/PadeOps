program HIT_AD_get3Dmeans
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use timer, only: tic, toc
   use exits, only: message, gracefulExit
   implicit none
   real(rkind), dimension(:,:,:), allocatable :: usum, vsum, wsum 
   real(rkind), dimension(:,:,:), allocatable :: R11sum, R12sum, R13sum, R22sum, R23sum, R33sum
   real(rkind), dimension(:,:,:), allocatable :: u, v, w
   real(rkind) :: dx, dy, dz
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 10.d0*pi, Ly = 2.d0*pi, Lz = 2.d0
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 2
   logical :: isZPeriodic = .true. 
   integer :: nt 
   logical :: file_found_u, file_found_v, file_found_w

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz
   
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)          
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)

   dx = Lx/real(nx,rkind) 
   dy = Ly/real(ny,rkind) 
   dz = Lz/real(nz,rkind)

   ! Initialize the operator class
   call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NUmericalSchemeVert)

   ! Allocate all the needed memory 
   call ops%allocate3DField(u)
   call ops%allocate3DField(v)
   call ops%allocate3DField(w)
   
   call ops%allocate3DField(usum)
   call ops%allocate3DField(vsum)
   call ops%allocate3DField(wsum)
   
   call ops%allocate3DField(R11sum)
   call ops%allocate3DField(R12sum)
   call ops%allocate3DField(R13sum)
   call ops%allocate3DField(R22sum)
   call ops%allocate3DField(R23sum)
   call ops%allocate3DField(R33sum)
  
   usum = 0.d0 
   vsum = 0.d0 
   wsum = 0.d0 
   R11sum = 0.d0 
   R12sum = 0.d0 
   R13sum = 0.d0 
   R22sum = 0.d0 
   R23sum = 0.d0 
   R33sum = 0.d0 

   nt = (tstop - tstart)/tstep + 1

   call message(0,"Number of snapshots to read:", nt)

   tidx = tstart
   !do while(tidx <= tstop) 
   !  file_found_u = ops%check_dump_existence("uVel",TIDX)
   !  file_found_v = ops%check_dump_existence("uVel",TIDX)
   !  file_found_w = ops%check_dump_existence("uVel",TIDX)
   !   
   !  if (file_found_u .and. file_found_v .and. file_found_w) then
   !      call message(0, "File succcessfully found for:", TIDX)
   !  else
   !      call message(0, "File missing for tid:", TIDX)
   !      call mpi_barrier(mpi_comm_world, ierr) 
   !      call gracefulExit("Some file was missing.", 324)
   !  end if 

   !  tidx = tidx + tstep
   !end do 

   tidx = tstart
   idx = 0
   do while(tidx <= tstop)
      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(u,"uVel",TIDX)
      call ops%ReadField3D(v,"vVel",TIDX)
      call ops%ReadField3D(w,"wVel",TIDX)

      usum = usum + u
      vsum = vsum + v
      wsum = wsum + w
      
      R11sum = R11sum + u*u
      R12sum = R12sum + u*v
      R13sum = R13sum + u*w
      R22sum = R22sum + v*v
      R23sum = R23sum + v*w
      R33sum = R33sum + w*w

      tidx = tidx + tstep
      idx = idx + 1
      call toc()
   end do 

   usum = usum/real(idx,rkind)
   vsum = vsum/real(idx,rkind)
   wsum = wsum/real(idx,rkind)

    
   R11sum = R11sum/real(idx,rkind)
   R12sum = R12sum/real(idx,rkind)
   R13sum = R13sum/real(idx,rkind)
   R22sum = R22sum/real(idx,rkind)
   R23sum = R23sum/real(idx,rkind)
   R33sum = R33sum/real(idx,rkind)

   call ops%WriteField3D(usum,"uMmn",tidx)
   call ops%WriteField3D(vsum,"vMmn",tidx)
   call ops%WriteField3D(wsum,"wMmn",tidx)
  
   R11sum = R11sum - usum*usum
   R12sum = R12sum - usum*vsum
   R13sum = R13sum - usum*wsum
   R22sum = R22sum - vsum*vsum
   R23sum = R23sum - vsum*wsum
   R33sum = R33sum - wsum*wsum


   call ops%writeField3D(R11sum,"R11m",tidx)
   call ops%writeField3D(R12sum,"R12m",tidx)
   call ops%writeField3D(R13sum,"R13m",tidx)
   call ops%writeField3D(R22sum,"R22m",tidx)
   call ops%writeField3D(R23sum,"R23m",tidx)
   call ops%writeField3D(R33sum,"R33m",tidx)
   
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

