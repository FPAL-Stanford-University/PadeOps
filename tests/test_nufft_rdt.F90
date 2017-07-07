program test_nufft_rdt
   use kind_parameters, only: rkind
   !use nufftMod, only: nufft
   use rdtMod, only: rdt
   use constants, only: two, pi
   use mpi
   use timer, only: tic, toc

   implicit none

   real(rkind), dimension(:,:,:), allocatable :: u, v, w, uExact, vExact, wExact

   integer, parameter :: nx = 128, ny = 128, nz = 128
   integer, parameter :: nk = 100, ntheta = 100
   type(rdt), allocatable :: my_rdt
   logical :: createRNGterms = .false., isStratified = .false. 
   real(rkind) :: kmin = 5.d0, kmax = 50.d0
   real(rkind), parameter :: dx = two*pi/real(nx,rkind), dy = two*pi/real(ny,rkind), dz = two*pi/real(nz,rkind) 
   real(rkind), dimension(2), parameter :: xLims = [0.d0, two*pi], yLims = [0.d0, two*pi], zLims = [0.d0, two*pi]
   integer :: seed1 = 31345, seed2 = 43245, seed3 = 78647
   real(rkind), parameter :: eps_nufft = 1.d-8

   allocate(u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz))
   allocate(uExact(nx, ny, nz), vExact(nx, ny, nz), wExact(nx, ny, nz))
   allocate(my_rdt)

   call my_rdt%init(nk, ntheta, kmin, kmax, isStratified, seed1, seed2, seed3, nx, ny, nz, xLims, yLims, zLims, eps_nufft, createRNGterms)

   call tic()
   call my_rdt%getFields(u, v, w)
   call toc()
   print*, sum(u)/real(nx*ny*nz)!(23,23,23)
   print*, sum(v)/real(nx*ny*nz)!(23,23,23)
   print*, sum(w)/real(nx*ny*nz)!(23,23,23)

   deallocate(my_rdt)
   deallocate(u, v, w, uExact, vExact, wExact)



end program 
