program test_BEMTheory
    use kind_parameters, only: rkind, clen, mpirkind
    use constants, only: pi, two, one, imi, zero, half, kappa
    use reductions, only: p_maxval, p_sum, p_minval 
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use actuatorDisk_RotMod, only: actuatorDisk_Rot
    use exits, only: message
    use mpi

    implicit none 

    type(actuatorDisk_Rot), dimension(:), allocatable :: hawts_Rot
    integer :: nx = 64, ny = 16, nz = 16, ntryparam = 2
    character(len=clen) :: inputDir = "/home/lele/nghaisas/Codes/PadeOps2/PadeOps/tests/test_BEMTheory_files"
    !character(len=clen) :: inputDir = "/fastscratch/nghaisas/runs/PadeOps/wupa/run3/bemt/4rot"
    character(len=clen) :: tempname, fname 
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG
    real(rkind), dimension(:,:,:), allocatable :: u, v, w, rhs1, rhsv, rhsw, rhs2
    real(rkind), dimension(:),     allocatable :: windspeed, angularspeed,pitch,power,thrustcoeff,shaft_torque,turbrotspeed
    real(rkind), parameter :: Lx = 28.8d0, Ly = 4.8d0, Lz = 3.0d0
    real(rkind) :: dx, dy, dz, diam, CT
    type(decomp_info) :: gp 
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0, num_turbines
    integer :: num_windspeeds, ioUnit 
    real(rkind) :: xPeriods = 2.d0, yPeriods = 2.d0, zpeak = 0.3d0, epsnd = 5.d0, z0init = 1.d-4 
    real(rkind) :: inst_val(8), inst_val_p(8), rhoRef, Uref, Lref, comp1, comp2, maxdiff

    call MPI_Init(ierr)

    ioUnit = 100; write(tempname,"(a)") 'input.dat'
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(ioUnit,file=fname,form='formatted',action='read',status='old')
    read(ioUnit,*) nx, ny, nz, ntryparam
    close(ioUnit)

    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    allocate(xG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(yG(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(zG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(u(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(v(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(w(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(rhs1(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhs2(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhsv(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhsw(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 


    dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
    ix1 = gp%xst(1); iy1 = gp%xst(2); iz1 = gp%xst(3)
    ixn = gp%xen(1); iyn = gp%xen(2); izn = gp%xen(3)
    do k=1,size(xG,3)
        do j=1,size(xG,2)
            do i=1,size(xG,1)
                xG(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                yG(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                zG(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
            end do
        end do
    end do
    xG = xG - dx; yG = yG - dy; zG = zG - dz 

    u = one; v = zero; w = zero

    num_windspeeds = 31
    allocate(windspeed(num_windspeeds),angularspeed(num_windspeeds),pitch(num_windspeeds))
    allocate(power(num_windspeeds),thrustcoeff(num_windspeeds),shaft_torque(num_windspeeds))
    allocate(turbrotspeed(num_windspeeds))

    ioUnit = 100; write(tempname,"(a)") 'WuPA_RenewEnergy2014_InputData.dat'
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(ioUnit,file=fname,form='formatted',action='read',status='old')
    read(ioUnit,*)       ! ignore header
    do i = 1, num_windspeeds
      read(ioUnit,*) windspeed(i), angularspeed(i), pitch(i)
    enddo
    close(ioUnit)

    num_turbines = 1
    allocate(hawts_Rot(num_turbines))
    do idx = 1,num_turbines
        call hawts_Rot(idx)%init(inputDir, idx, xG, yG, zG, ntryparam)
    end do 

    rhoRef = 1.2d0
    Uref = 10.0d0 ! m/s
    Lref = hawts_Rot(1)%Lref ! m (make sure this is cosnistent with turbine input files)

    ! non-dimensionalize all variables
    windspeed = windspeed/Uref
    angularspeed = angularspeed*2.0d0*pi/60.0d0 / (Uref/Lref)
    pitch = pitch*pi/180.0d0

    ! reset diam, CT
    diam = hawts_Rot(1)%diam
 
    rhs1 = 0.d0; rhsv = 0.0d0; rhsw = 0.0d0
    call tic()
    do j = 1, num_windspeeds
      call mpi_barrier(mpi_comm_world, ierr)
      u = windspeed(j)
      do idx = 1,num_turbines
          hawts_Rot(idx)%TSR = angularspeed(j)*diam/two/windspeed(j)
          hawts_Rot(idx)%blade_pitch = pitch(j)
          call hawts_Rot(idx)%get_RHS(u, v, w, rhs1, rhsv, rhsw, inst_val_p)
      enddo
      call MPI_reduce(inst_val_p, inst_val, 8*num_turbines, mpirkind, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      thrustcoeff(j)  = -inst_val(1)/(0.5d0*pi*diam*diam/4.0d0*windspeed(j)*windspeed(j))
      power(j)        = -inst_val(5)*rhoRef*Lref**2*Uref**3/1.0d6
      shaft_torque(j) = -inst_val(3)*rhoref*Lref**3*Uref**2/1.0d3
      turbrotspeed(j) = inst_val(4)/(2.0d0*pi/60.0d0 / (Uref/Lref))

      inst_val_p = 0.0d0; inst_val = 0.0d0

      call mpi_barrier(mpi_comm_world, ierr)
      call message(2,"CT ",thrustcoeff(j))
    enddo
    call toc()

    windspeed = windspeed*Uref
    angularspeed = angularspeed/(2.0d0*pi/60.0d0 / (Uref/Lref))
    pitch = pitch/(pi/180.0d0)

    ioUnit = 100+nrank; write(tempname,"(a,i4.4,a)") 'WuPA_RenewEnergy2014_OutputData_',nrank,'.dat'
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(ioUnit,file=fname,form='formatted',action='write',status='unknown')
    do i = 1, num_windspeeds
      write(ioUnit,'(6(e19.12,1x))') windspeed(i), turbrotspeed(i), pitch(i), thrustcoeff(i), power(i), shaft_torque(i)
    enddo
    close(ioUnit)

    deallocate(hawts_Rot,windspeed,angularspeed,pitch,thrustcoeff,power,shaft_torque,turbrotspeed)
    deallocate(xG, yG, zG, u, v, w, rhs1, rhs2)
    call MPI_Finalize(ierr)
end program
