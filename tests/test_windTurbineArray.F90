program test_wtarray
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero, half, kappa
    use reductions, only: p_maxval, p_sum 
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use actuatorLineMod, only: actuatorLine
    use exits, only: message
    use mpi
    use spectralMod, only: spectral
    use TurbineMod, only: TurbineArray

    implicit none 

    type(turbineArray), allocatable, target :: WindTurbineArr

    integer, parameter :: nx = 96, ny = 96, nz = 64
    character(len=clen) :: inputFile = "/home/nghaisas/PadeOps/tests/test_actuatorLine_files/input.in"
    real(rkind), dimension(:,:,:,:), allocatable, target :: mesh
    real(rkind), dimension(:,:,:), pointer :: xG, yG, zG
    real(rkind), dimension(:,:,:), allocatable :: u, v, wC
    complex(rkind), dimension(:,:,:), pointer :: rhsu, rhsv, rhsw
    complex(rkind), dimension(:,:,:,:), allocatable, target :: rhsE, rhsC
    real(rkind), dimension(:), allocatable :: inst_horz_avg_turb
    real(rkind), parameter :: Lx = pi, Ly = pi, Lz = one
    real(rkind) :: dx, dy, dz

    complex(rkind), dimension(:,:,:,:), allocatable, target :: cbuffyC, cbuffzC, cbuffyE, cbuffzE
    real(rkind), dimension(:,:,:,:), allocatable, target :: rbuffxC


    type(decomp_info), allocatable :: gp, gpE
    type(decomp_info), pointer :: sp_gp, sp_gpE 
    type(spectral), allocatable, target :: spectE, spectC
    character(len=clen) :: filter_x="cf90" ! What filter to use: "cf90", "gaussian", "lstsq", "spectral"
    character(len=clen) :: fname
    type(actuatorLine), pointer :: actLine

    integer :: ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 1, pcol = 4, turbID, blID, ptID, tID
    real(rkind) :: xPeriods = 2.d0, yPeriods = 2.d0, zpeak = 0.3d0, epsnd = 5.d0, z0init = 1.d-4 
    real(rkind) :: dt = 0.1_rkind
    logical     :: periodicbcs(3)

    namelist /DECOMP/ prow, pcol

    periodicbcs(1) = .true.;   periodicbcs(2) = .true.;   periodicbcs(3) = .false.

    call MPI_Init(ierr)

    allocate(gp); allocate(gpE)

    open(unit=10,file=trim(inputFile),status='old',action='read')
    read(unit=10,NML=DECOMP)    
    close(10)

    call decomp_2d_init(nx, ny, nz, prow, pcol, periodicbcs)
    call get_decomp_info(gp)
    call decomp_info_init(nx,ny,nz+1,gpE)

    dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
    !filter_x = "cf90"
    allocate(spectC);     call spectC%init("x", nx, ny, nz, dx, dy, dz, "four", filter_x, 2 , .false.)
    allocate(spectE);     call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", filter_x, 2 , .false.)

    sp_gp => spectC%spectdecomp
    sp_gpE => spectE%spectdecomp


    allocate(mesh(gp%xsz(1),gp%xsz(2),gp%xsz(3),3));
    !allocate(yG(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    !allocate(zG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); 
    allocate(u(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(v(gp%xsz(1),gp%xsz(2),gp%xsz(3))); 
    allocate(wC(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    call spectC%alloc_r2c_out(rhsC,2);    rhsu => rhsC(:,:,:,1); rhsv => rhsC(:,:,:,2)
    call spectE%alloc_r2c_out(rhsE,1);    rhsw => rhsE(:,:,:,1)
    allocate(rbuffxC(gp%xsz(1),gp%xsz(2),gp%xsz(3),2)) 
    allocate(cbuffyC(sp_gp%ysz(1),sp_gp%ysz(2),sp_gp%ysz(3),2)) 
    allocate(cbuffzC(sp_gp%zsz(1),sp_gp%zsz(2),sp_gp%zsz(3),2)) 
    allocate(cbuffyE(sp_gpE%ysz(1),sp_gpE%ysz(2),sp_gpE%ysz(3),2)) 
    allocate(cbuffzE(sp_gpE%zsz(1),sp_gpE%zsz(2),sp_gpE%zsz(3),2)) 

    call mpi_barrier(mpi_comm_world, ierr); print *, '--',0;     call mpi_barrier(mpi_comm_world, ierr)

    xG => mesh(:,:,:,1);   yG => mesh(:,:,:,2);   zG => mesh(:,:,:,3)

    ix1 = gp%xst(1); iy1 = gp%xst(2); iz1 = gp%xst(3)
    ixn = gp%xen(1); iyn = gp%xen(2); izn = gp%xen(3)
    do k=1,size(mesh,3)
        do j=1,size(mesh,2)
            do i=1,size(mesh,1)
                xG(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                yG(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                zG(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
            end do
        end do
    end do
    xG = xG - dx; yG = yG - dy; zG = zG - dz 

    call mpi_barrier(mpi_comm_world, ierr); print *, '--',1, size(mesh,1), size(mesh,2), size(mesh,3), size(mesh,4);     call mpi_barrier(mpi_comm_world, ierr)
    allocate(WindTurbineArr)
    call WindTurbineArr%init(inputFile, gp, gpE, spectC, spectE, cbuffyC, cbuffyE, cbuffzC, cbuffzE, mesh, dx, dy, dz)
    allocate(inst_horz_avg_turb (8*WindTurbineArr%nTurbines))
    call mpi_barrier(mpi_comm_world, ierr); print *, '--',2;     call mpi_barrier(mpi_comm_world, ierr)

    u = (one/kappa)*log(zG/z0init) + epsnd*cos(Yperiods*two*pi*yG/Ly)*exp(-half*(zG/zpeak/Lz)**2)
    v = epsnd*(zG/Lz)*cos(Xperiods*two*pi*xG/Lx)*exp(-half*(zG/zpeak/Lz)**2)
    wC= zero  
    u = one; v = zero; wC = zero

    tid = 0
    write(fname,'(a,i4.4,a)') 'tec_global_',nrank,'.dat'
    open(10,file=fname,status='replace')
    write(10,*) 'VARIABLES="X","Y","Z","u","v","w","fx","fy","fz"'
    write(10,*) 'ZONE I=', gp%xsz(1), 'J=', gp%xsz(2), 'K=', gp%xsz(3), 'ZONETYPE=ORDERED'
    write(10,*) 'DATAPACKING=POINT, SOLUTIONTIME=',real(tid,rkind)
    do k=1,gp%xsz(3)
     do j=1,gp%xsz(2)
      do i=1,gp%xsz(1)
        write(10,'(9(e19.12,1x))') xG(i,j,k), yG(i,j,k), zG(i,j,k), u(i,j,k), v(i,j,k), wC(i,j,k), &
                                   WindTurbineArr%fx(i,j,k), WindTurbineArr%fy(i,j,k), WindTurbineArr%fz(i,j,k)
      enddo
     enddo
    enddo
    close(10)
    do turbID = 1, WindTurbineArr%nTurbines
      actLine => WindTurbineArr%turbArrayALM(turbID)
      write(fname,'(a1,i3.3,a11)') 'T',turbID,'_points.dat'
      if(nrank==0) then
        open(10,FILE=fname,status='replace')
        write(10,*) 'VARIABLES="BladeID","PointID","X","Y","Z"'
        write(10,*) 'ZONE T="SOLUTIONTIME=',real(tid,rkind), ' ", F=POINT, I=', actLine%num_blades*actLine%num_blade_points
        do blID = 1, actLine%num_blades
          do ptID = 1, actLine%num_blade_points
            write(10,'(2(i4,1x),3(e19.12,1x))') blID, ptID, actLine%blade_points(:, ptID, blID)
          enddo
        enddo
        close(10)
      endif
    enddo


    do tid = 1, 1!0!num_time_steps

      call mpi_barrier(mpi_comm_world, ierr)
      if(nrank==0) write(*,*) 'Time Step No.', tid
      call tic()
      call WindTurbineArr%getForceRHS(dt, u, v, wC, rhsu, rhsv, rhsw, .true., inst_horz_avg_turb)
      call mpi_barrier(mpi_comm_world, ierr)
      call toc()
  
      write(fname,'(a,i4.4,a)') 'tec_global_',nrank,'.dat'
      open(10,file=fname,status='old',action='write',position='append')
      write(10,*) 'ZONE I=', gp%xsz(1), 'J=', gp%xsz(2), 'K=', gp%xsz(3), 'ZONETYPE=ORDERED'
      write(10,*) 'DATAPACKING=POINT, SOLUTIONTIME=',real(tid, rkind)
      write(10,*) 'VARSHARELIST=([1,2,3]=1)'
      do k=1,gp%xsz(3)
       do j=1,gp%xsz(2)
        do i=1,gp%xsz(1)
          write(10,'(9(e19.12,1x))') u(i,j,k), v(i,j,k), wC(i,j,k), WindTurbineArr%fx(i,j,k), &
                                          WindTurbineArr%fy(i,j,k), WindTurbineArr%fz(i,j,k)
        enddo
       enddo
      enddo
      close(10)
      do turbID = 1, WindTurbineArr%nTurbines
        actLine => WindTurbineArr%turbArrayALM(turbID)
        write(fname,'(a1,i3.3,a11)') 'T',turbID,'_points.dat'
        if(nrank==0) then
          !open(10,FILE=fname,status='replace')
          open(10,FILE=fname,status='old',action='write',position='append')
          !write(10,*) 'VARIABLES="BladeID","PointID","X","Y","Z"'
          write(10,*) 'ZONE T="SOLUTIONTIME=',real(tid,rkind), ' ", F=POINT, I=', actLine%num_blades*actLine%num_blade_points
          do blID = 1, actLine%num_blades
            do ptID = 1, actLine%num_blade_points
              write(10,'(2(i4,1x),3(e19.12,1x))') blID, ptID, actLine%blade_points(:, ptID, blID)
            enddo
          enddo
          close(10)
        endif
      enddo

    enddo

    !call message(2,"Computed Source:", p_sum(sum(rhs)) * dx*dy*dz)
    !call message(2,"Expected Source:", (6.d0*0.5d0*(pi/4.d0)*(0.08d0**2)*0.65d0))

    !call message(2,"error:", 100.d0*abs(p_sum(sum(rhs)) * dx*dy*dz - 6.d0*0.5d0*(pi/4.d0)*(0.08d0**2)*0.65d0) / (6.d0*0.5d0*(pi/4.d0)*(0.08d0**2)*0.65d0))

    !------write out all diagnostics----------
    ! blade-setup interpolated values
    if(nrank==0) then
      do turbID = 1, WindTurbineArr%nTurbines
       actLine => WindTurbineArr%turbArrayALM(turbID)
       do blID = 1, actLine%num_blades
        write(fname,'(a1,i3.3,a3,i3.3,a15)') 'T',turbID,'_Bl',blID,'_interpvals.dat'
        open(10,File=fname,status='replace',action='write')
        do k = 1, actLine%num_blade_points
           write(10,'(i4,3(e19.12,1x),i4)') k, actLine%radial_dist(k,blID), actLine%twistAng(k,blID), actLine%chord(k,blID), actLine%airfoilIDIndex(k,blID)
        enddo
        close(10)
       enddo
      nullify(actLine)
      enddo
    endif

    do turbID = 1, WindTurbineArr%nTurbines
      actLine => WindTurbineArr%turbArrayALM(turbID)
      write(fname,'(a1,i3.3,a8)') 'T',turbID,'_log.dat'
      if(nrank==0) then
        open(10,FILE=fname,status='replace')
        write(10,'(a,3(e19.12,1x))') 'xLoc,    yLoc,    zLoc   = ', actLine%xLoc, actLine%yLoc, actLine%zLoc
        write(10,'(a,3(e19.12,1x))') 'tip_rad, hub_rad, hub_ht = ', actLine%tip_radius, actLine%hub_radius, actLine%hub_height
        write(10,'(a,3(e19.12,1x))') 'yaw_ang, azimuth, rotspd = ', actLine%yaw_angle, actLine%blade_azimuth, actLine%rotSpeed
        write(10,'(a,3(e19.12,1x))') 'nac_width, delta, nrmfac = ', actLine%nacelle_width, actLine%delta, actLine%normfactor
        write(10,*)
        write(10,'(a,3(e19.12,1x))') 'turbLoc                  = ', actLine%turbLoc
        write(10,'(a,3(e19.12,1x))') 'rotor_center             = ', actLine%rotor_center
        write(10,'(a,3(e19.12,1x))') 'rotor_shaft              = ', actLine%rotor_shaft 
        write(10,'(a,3(e19.12,1x))') 'turb_torque, turb_thrust = ', actLine%turb_torque, actLine%turb_thrust
        write(10,*)
        write(10,*)
        close(10)
      endif
      call mpi_barrier(mpi_comm_world,ierr)
      do k = 0, prow*pcol-1
        if(k==nrank) then
          open(10,FILE=fname,status='old',action='write',position='append')
          write(10,'(a,i3,a)') '---nrank = ', k, '----'
          write(10,*) actLine%Am_I_Active, actLine%Am_I_Split
          write(10,'(a,9(i4,1x))') 'ist, iend, jst, jend, kst, kend = ', actLine%ist, actLine%iend, actLine%jst, actLine%jend, actLine%kst, actLine%kend, actLine%xlen, actLine%ylen, actLine%zlen
          if(actLine%Am_I_Active) then
            write(10,'(a,6(e19.12,1x))') 'xst, xend, yst, yend, zst, zend = ', xG(actLine%ist,1,1), xG(actLine%iend,1,1), yG(1,actLine%jst,1), yG(1,actLine%jend,1), zG(1,1,actLine%kst), zG(1,1,actLine%kend)
          endif
          write(10,'(a,4(e19.12,1x))')   'xR, yR, zL, zR                  = ', actLine%xRightPad, actLine%yRightPad, actLine%zLeftPad, actLine%zRightPad
          write(10,'(a,3(e19.12,1x))')   'distr_thrust                    = ', actLine%distr_thrust
          close(10)
        endif
        call mpi_barrier(mpi_comm_world,ierr)
      enddo
      nullify(actLine)
    enddo

    
    ! test turbine blade rotation
    ! test yaw motion

    call WindTurbineArr%destroy()
    call mpi_barrier(mpi_comm_world, ierr); print *, '--',5;     call mpi_barrier(mpi_comm_world, ierr)

    nullify(xG, yG, zG, rhsw, rhsv, rhsu)
    deallocate(mesh, u, v, wC, rhsE, rhsC)
    call mpi_barrier(mpi_comm_world, ierr); print *, '--',6;     call mpi_barrier(mpi_comm_world, ierr)
    deallocate(rbuffxC)
    call mpi_barrier(mpi_comm_world, ierr); print *, '1-',6;     call mpi_barrier(mpi_comm_world, ierr)
    deallocate(cbuffzC)
    call mpi_barrier(mpi_comm_world, ierr); print *, '3-',6;     call mpi_barrier(mpi_comm_world, ierr)
    deallocate(cbuffyE)
    call mpi_barrier(mpi_comm_world, ierr); print *, '4-',6;     call mpi_barrier(mpi_comm_world, ierr)
    deallocate(cbuffzE)
    call mpi_barrier(mpi_comm_world, ierr); print *, '5-',6;     call mpi_barrier(mpi_comm_world, ierr)
    deallocate(inst_horz_avg_turb)
    call mpi_barrier(mpi_comm_world, ierr); print *, '--',7;     call mpi_barrier(mpi_comm_world, ierr)
    deallocate(cbuffyC)
    call mpi_barrier(mpi_comm_world, ierr); print *, '2-',6;     call mpi_barrier(mpi_comm_world, ierr)

    call MPI_Finalize(ierr)
end program 
