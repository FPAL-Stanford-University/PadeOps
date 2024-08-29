! Template for PadeOps

#include "HIT_Periodic_moving_files/initialize.F90"       
#include "HIT_Periodic_moving_files/temporalHook.F90"  

program HIT_Periodic
    use mpi
    use kind_parameters,  only: clen
    use HIT_periodic_parameters, only: useBandpassFilter 
    use IncompressibleGrid, only: igrid
    use HIT_Periodic_parameters, only: k_bp_left, k_bp_right
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use exits, only: message
    use budgets_vol_avg_mod,   only: budgets_vol_avg  
    !DEBUG
    use basic_io, only: write_2d_ascii
    use decomp_2d, only: nrank
    !END DEBUG

    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr
    type(budgets_vol_avg)   :: budg_volavg
    !DEBUG
    integer :: nx, i
    character(len=clen) :: fname, fpath
    !END DEBUG

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize hit_grid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the hit_grid solver (see hit_grid.F90)
  
    call igp%start_io(.true.)                !<-- Start I/O by creating a header file (see io.F90)

    !nx = igp%gpC%xsz(1)
    !write(fpath,'(A)')'/scratch/06632/ryanhass/LES/HITdecay/shearlessMixing_Bodart/test_TPC'
    !do i = 1,nx
    !    write(fname,'(A,I2.2,A,I3.3,A)')'/rank',nrank,'_x',i,'_u.out'
    !    call write_2d_ascii(igp%u(i,:,:),trim(fpath)//trim(fname))
    !    write(fname,'(A,I2.2,A,I3.3,A)')'/rank',nrank,'_x',i,'_v.out'
    !    call write_2d_ascii(igp%v(i,:,:),trim(fpath)//trim(fname))
    !    write(fname,'(A,I2.2,A,I3.3,A)')'/rank',nrank,'_x',i,'_wC.out'
    !    call write_2d_ascii(igp%wC(i,:,:),trim(fpath)//trim(fname))
    !    write(fname,'(A,I2.2,A,I3.3,A)')'/rank',nrank,'_x',i,'_w.out'
    !    call write_2d_ascii(igp%w(i,:,:),trim(fpath)//trim(fname))
    !    write(fname,'(A,I2.2,A,I3.3,A)')'/rank',nrank,'_x',i,'_y.out'
    !    call write_2d_ascii(igp%mesh(i,:,:,2),trim(fpath)//trim(fname))
    !    write(fname,'(A,I2.2,A,I3.3,A)')'/rank',nrank,'_x',i,'_z.out'
    !    call write_2d_ascii(igp%mesh(i,:,:,3),trim(fpath)//trim(fname))
    !    write(fname,'(A,I2.2,A,I3.3,A)')'/rank',nrank,'_x',i,'_zE.out'
    !    call write_2d_ascii(igp%meshE(i,:,:,3),trim(fpath)//trim(fname))
    !end do

    call igp%printDivergence()

    ! Initialize bandpass filtering 
!    if (useBandpassFilter) then
!         call igp%spectC%init_bandpass_filter(k_bp_left, k_bp_right, igp%cbuffzC(:,:,:,1), igp%cbuffyC(:,:,:,1))
!    end if 

!    call budg_volavg%init(inputfile, igp)   !<-- Budget class initialization 
  
    call tic() 
    do while (igp%tsim < igp%tstop .and. igp%step < igp%nsteps) 
      
       call igp%timeAdvance()     !<-- Time stepping scheme + Pressure Proj. (see igrid.F90)
 !      call budg_volavg%doBudgets()       
       call doTemporalStuff(igp)     !<-- Go to the temporal hook (see temporalHook.F90)
       
    end do 
 
!    call budg_volavg%doBudgets(.true.)      !<-- Force dump budget information if active

!    call igp%finalize_io()                  !<-- Close the header file (wrap up i/o)

    call budg_volavg%destroy()          !<-- release memory taken by the budget class 

    call igp%destroy()                !<-- Destroy the IGRID derived type 
   

    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program
