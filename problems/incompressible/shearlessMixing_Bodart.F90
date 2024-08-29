! Template for PadeOps

#include "shearlessMixing_Bodart_files/initialize.F90"       
#include "shearlessMixing_Bodart_files/temporalHook.F90"  

program shearlessMixing
    use mpi
    use kind_parameters,    only: clen, rkind
    use IncompressibleGrid, only: igrid, readField3D
    use temporalhook,       only: doTemporalStuff
    use exits,              only: message, gracefulExit
    use budgets_xy_avg_mod, only: budgets_xy_avg
    use stats_xy_mod,       only: stats_xy
    use timer,              only: tic, toc 
    use fortran_assert,     only: assert
    use decomp_2d,          only: nproc, nrank
    use decomp_2d_io,       only: decomp_2d_write_one
    use basic_io,           only: write_2d_ascii

    implicit none

    type(igrid), allocatable, target :: SM
    type(stats_xy), dimension(:), allocatable :: stats
    character(len=clen) :: inputfile, tempname, stats_info_dir
    integer :: ierr, n, ioUnit, stid
    integer :: num_stats_instances = 1
    !DEBUG
    real(rkind), dimension(:,:,:), allocatable :: u, v, w
    integer :: k
    !END DEBUG

    ! Required for reading the namelist, but not used directly in the main program
    real(rkind) :: Lx, Ly, Lz, zmin, Tref
    logical :: symmetricDomain 
    
    call MPI_Init(ierr)

    call GETARG(1,inputfile)                                            

    namelist /SMinput/ Lx, Ly, Lz, symmetricDomain, zmin, Tref, stats_info_dir, num_stats_instances

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=SMinput)
    close(ioUnit)    

    allocate(SM)                                                     
    
    call SM%init(inputfile, .true.)                     
    call SM%start_io(.true.)                                   
    call SM%printDivergence()

    ! DEBUG
    !allocate(u(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    !allocate(v(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    !allocate(w(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    !call readField3D(74,999999,trim(SM%inputdir),'uVel',u,SM%gpC)
    !call readField3D(74,999999,trim(SM%inputdir),'vVel',v,SM%gpC)
    !call readField3D(74,999999,trim(SM%inputdir),'wVel',w,SM%gpC)
    !call SM%dumpFullField(u,'uVel',step=999)
    !call SM%dumpFullField(v,'vVel',step=999)
    !call SM%dumpFullField(w,'wVel',step=999)
    !call decomp_2d_write_one(1,u,trim(SM%outputdir)//'/Run74_uVel_t009999.out',SM%gpC)
    !call decomp_2d_write_one(1,v,trim(SM%outputdir)//'/Run74_vVel_t009999.out',SM%gpC)
    !call decomp_2d_write_one(1,w,trim(SM%outputdir)//'/Run74_wVel_t009999.out',SM%gpC)
    !! Dump ascii data
    !do k = 1,SM%gpC%xsz(2)
    !    write(tempname,'(A,I2.2,A,I3.3,A)')trim(SM%outputdir)//'/rank',nrank,'_u_y',k,'.out'
    !    call write_2d_ascii(u(:,k,:),trim(tempname))
    !    write(tempname,'(A,I2.2,A,I3.3,A)')trim(SM%outputdir)//'/rank',nrank,'_v_y',k,'.out'
    !    call write_2d_ascii(v(:,k,:),trim(tempname))
    !    write(tempname,'(A,I2.2,A,I3.3,A)')trim(SM%outputdir)//'/rank',nrank,'_w_y',k,'.out'
    !    call write_2d_ascii(w(:,k,:),trim(tempname))
    !end do
    !call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !call gracefulExit('stop',ierr)
    ! END DEBUG

    ! Initialize stats class instances
    if (num_stats_instances > 0) then
        allocate(stats(num_stats_instances))
        do stid = 1,size(stats)
            write(tempname,'(A,I2.2,A4)')trim(stats_info_dir)//"/STATS_",stid,".inp"
            call stats(stid)%init(trim(tempname),SM)
        end do
    end if

    call message("==========================================================")
    call message(0, "All memory allocated! Now running the simulation.")

    do while (SM%tsim < SM%tstop .and. SM%step < SM%nsteps)
       call tic() 
       call SM%timeAdvance()
       do n = 1,num_stats_instances
           call stats(n)%compute_stats()
       end do
       call doTemporalStuff(SM) 
       call toc()
    end do

    call message("==========================================================")
    call message(0,"Finalizing simulation")
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    do stid = 1,size(stats)
        call stats(stid)%destroy()
    end do
    if (allocated(stats)) deallocate(stats)
    call SM%finalize_io()
    call SM%destroy()
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call message(0,"Simulation finalized")
   
    !deallocate(SM)
    
    call MPI_Finalize(ierr)

end program
