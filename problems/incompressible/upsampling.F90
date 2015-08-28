! This program reads in the hit3d data and upsamples it on a grid with 
! Nx_up = (3/2)*Nx, Ny_up = (3/2)*Ny, Nz_up = (3/2)*Nz

#include "hit_files/meshgen.F90"         !<-- Meshgen File (USER SPECIFIED)
#include "hit_files/io.F90"              !<-- I/O Procedures (USER SPECIFIED)
#include "hit_files/initfields.F90"      !<-- Initializing Fields FIle (USER SPECIFIED)
#include "upsampling_files/hit3dwr.F90"  !<-- writing hit 3d style output 
program upsampling
    use mpi
    use kind_parameters,  only: rkind,clen,stdout,stderr
    use IncompressibleGrid, only: igrid
    use hitCD_IO, only: start_io, finalize_io
    use constants, only: half 
    use decomp_2d, only: nproc 
    use hit3dwr, only: write_upsampled_file, zeroPad
    use exits, only: GracefulExit, message
    

    implicit none

    include "fftw3.f"
    real(rkind), parameter :: upSamplingFactor = 1.5_rkind
    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr
    integer(kind=8) :: planS, planL

    complex(rkind), dimension(:,:,:), allocatable :: fhat
    real(rkind), dimension(:,:,:), allocatable :: uup, vup, wup
    complex(rkind), dimension(:,:,:), allocatable :: fuphat
    integer :: NxL, NyL, NzL, Nx, Ny, Nz

    call MPI_Init(ierr)               !<-- Begin MPI

    call GETARG(1,inputfile)          !<-- Get the location of the input file

    allocate(igp)                     !<-- Initialize igrid with defaults

    call igp%init(inputfile)          !<-- Properly initialize the igrid solver (see igrid.F90)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (nproc .ne. 1) then
        call GracefulExit("This program is intended to be run using a single processor", 1)
    end if 
    Nx = igp%nx
    Ny = igp%ny
    Nz = igp%nz

    ! STEP 1: Create the upsampled domain FFT
    NxL = floor(upSamplingFactor * Nx)
    NyL = floor(upSamplingFactor * Ny)
    NzL = floor(upSamplingFactor * Nz)

    ! STEP 2: Allocate the two arrays 
    allocate(fhat(Nx,Ny,Nz))   
    allocate(fuphat(NxL,NyL,NzL))   
    allocate(uup(NxL,NyL,NzL))   
    allocate(vup(NxL,NyL,NzL))   
    allocate(wup(NxL,NyL,NzL))   
 
    ! STEP 2: Create the fftw plans
   
    call dfftw_plan_dft_3d(planS, Nx,Ny,Nz, fhat,fhat, &
            &                       FFTW_FORWARD, FFTW_ESTIMATE)

    call dfftw_plan_dft_3d(planL, NxL,NyL,NzL, fuphat,fuphat, &
            &                       FFTW_BACKWARD, FFTW_ESTIMATE)

    ! STEP 3: Upsample u
    fhat = igp%u
    call dfftw_execute_dft(planS,fhat,fhat)
    call zeroPad(fhat,fuphat)
    call dfftw_execute_dft(planL,fuphat,fuphat)
    uup = real(fuphat/real(NxL)/real(NyL)/real(NzL))

    ! STEP 4: Upsample v
    fhat = igp%v
    call dfftw_execute_dft(planS,fhat,fhat)
    call zeroPad(fhat,fuphat)
    call dfftw_execute_dft(planL,fuphat,fuphat)
    vup = real(fuphat/real(NxL)/real(NyL)/real(NzL))

    ! STEP 5: Upsample w
    fhat = igp%w
    call dfftw_execute_dft(planS,fhat,fhat)
    call zeroPad(fhat,fuphat)
    call dfftw_execute_dft(planL,fuphat,fuphat)
    wup = real(fuphat/real(NxL)/real(NyL)/real(NzL))
    
    ! STEP 6:  Write the upsampled velocities
    call write_upsampled_file(igp%inputdir,uup, vup, wup)
    call message(1,"Upsampled files written successfully!")

    ! STEP 7: Deallocate all the memory
    deallocate(fuphat, fhat,uup, vup, wup)
    call dfftw_destroy_plan(planS)
    call dfftw_destroy_plan(planL)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call igp%destroy()                !<-- Destroy the IGRID derived type 
    
    deallocate(igp)                   !<-- Deallocate all the memory associated with scalar defaults
    
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

contains

    

end program 
