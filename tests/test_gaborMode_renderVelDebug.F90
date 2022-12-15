#include "../problems/gabor/hitEnrich_files/initialize.F90"
program test_gaborMode_renderVelMPI
  use kind_parameters,       only: clen, rkind
  use enrichmentMod,         only: enrichmentOperator, nthreads, xDom, yDom, zDom
  use IncompressibleGrid,    only: igrid
  use auxiliary_openmp_subs, only: GetArguments
  use mpi
  use exits,                 only: message
  use decomp_2d,             only: nrank
  use reductions,            only: p_maxval
  use fortran_assert,          only: assert
  implicit none

  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  integer :: ioUnit, ierr, provided
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator) :: enrich
  real(rkind), dimension(:,:,:), allocatable :: u1, v1, w1
  real(rkind), dimension(2,2) :: xpos, ypos, zpos, kx, ky, kz, uhatR, uhatI, vhatR, vhatI, whatR, whatI
  real(rkind) :: tol = 1.d-13
  integer :: id1, id2, id3
  character(len=200) :: id1str, id2str, id3str
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  call GETARG(4,id1str)
  call GETARG(5,id2str)
  call GETARG(6,id3str)
  call GetArguments(nthreads,4)

  read(id1str,*) id1
  read(id2str,*) id2
  read(id3str,*) id3

  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM)
  allocate(u1(size(enrich%smallScales%u, 1),size(enrich%smallScales%u, 2),size(enrich%smallScales%u, 3)))
  allocate(v1(size(enrich%smallScales%v, 1),size(enrich%smallScales%v, 2),size(enrich%smallScales%v, 3)))
  allocate(w1(size(enrich%smallScales%wC, 1),size(enrich%smallScales%wC, 2),size(enrich%smallScales%wC, 3)))

  enrich%modeData(:,7:12) = 0.d0
  enrich%x(1)     =    3.11347110697574d0     
  enrich%y(1)     =   0.846489344291164d0     
  enrich%z(1)     =    1.50084108828105d0     
  enrich%kx(1)    =   7.538928350430579d-003
  enrich%ky(1)    =   9.203686525045829d-003
  enrich%kz(1)    =   -9.49999255035043d0     
  enrich%uhatR(1) =  -9.796411828818073d-002
  enrich%uhatI(1) =  -8.461212861900098d-002
  enrich%vhatR(1) =  -0.199034603443980d0     
  enrich%vhatI(1) =  -0.171907242779377d0     
  enrich%whatR(1) =  -2.705682717959693d-004
  enrich%whatI(1) =  -2.336912515874032d-004
 
  enrich%x(    2) =    3.22619062375460d0     
  enrich%y(    2) =   0.406661756287087d0     
  enrich%z(    2) =   0.929066698277791d0     
  enrich%kx(   2) =  -9.738632831059054d-003
  enrich%ky(   2) =  -4.393339142904240d-003
  enrich%kz(   2) =    9.49999399250346d0     
  enrich%uhatR(2) =   7.563858193824458d-002
  enrich%uhatI(2) =   2.096057085955508d-002
  enrich%vhatR(2) =   0.272165835790534d0     
  enrich%vhatI(2) =   7.542118242374235d-002
  enrich%whatR(2) =   2.034036230573583d-004
  enrich%whatI(2) =   5.636615527330887d-005
 
  enrich%x(    3) =    5.54192067941069d0     
  enrich%y(    3) =   0.930500340497201d0     
  enrich%z(    3) =   1.847319627107817d-002
  enrich%kx(   3) =  -5.685163495892143d-003
  enrich%ky(   3) =   9.045513488339469d-003
  enrich%kz(   3) =    9.49999399250346d0     
  enrich%uhatR(3) =   0.247769501861132d0     
  enrich%uhatI(3) =   0.113520915885509d0     
  enrich%vhatR(3) =   9.811114291935046d-002
  enrich%vhatI(3) =   4.495172617742562d-002
  enrich%whatR(3) =   5.485734634718119d-005
  enrich%whatI(3) =   2.513407079403549d-005
  
  enrich%tid = 0
  call enrich%wrapupTimeStep()
  u1 = enrich%smallScales%u
  v1 = enrich%smallScales%v
  w1 = enrich%smallScales%wC
  print*, maxval(abs(u1))

  enrich%modeData(:,7:12) = 0.d0
  enrich%x(    id1) =    3.11347110697574d0     
  enrich%y(    id1) =   0.846489344291156d0     
  enrich%z(    id1) =    1.50084108828105d0     
  enrich%kx(   id1) =   7.538928350430579d-003
  enrich%ky(   id1) =   9.203686525045829d-003
  enrich%kz(   id1) =   -9.49999255035043d0     
  enrich%uhatR(id1) =  -9.796411828818073d-002
  enrich%uhatI(id1) =  -8.461212861900098d-002
  enrich%vhatR(id1) =  -0.199034603443980d0     
  enrich%vhatI(id1) =  -0.171907242779377d0     
  enrich%whatR(id1) =  -2.705682717959693d-004
  enrich%whatI(id1) =  -2.336912515874032d-004
 
  enrich%x(    id2) =    3.22619062375460d0     
  enrich%y(    id2) =   0.406661756287078d0     
  enrich%z(    id2) =   0.929066698277791d0     
  enrich%kx(   id2) =  -9.738632831059054d-003
  enrich%ky(   id2) =  -4.393339142904240d-003
  enrich%kz(   id2) =    9.49999399250346d0     
  enrich%uhatR(id2) =   7.563858193824458d-002
  enrich%uhatI(id2) =   2.096057085955508d-002
  enrich%vhatR(id2) =   0.272165835790534d0     
  enrich%vhatI(id2) =   7.542118242374235d-002
  enrich%whatR(id2) =   2.034036230573583d-004
  enrich%whatI(id2) =   5.636615527330887d-005
  
  enrich%x(id3)     =     5.54192067941069d0      
  enrich%y(id3)     =    0.930500340497191d0      
  enrich%z(id3)     =    1.847319627107817d-002 
  enrich%kx(id3)    =   -5.685163495892143d-003
  enrich%ky(id3)    =    9.045513488339469d-003
  enrich%kz(id3)    =     9.49999399250346d0      
  enrich%uhatR(id3) =    0.247769501861132d0     
  enrich%uhatI(id3) =    0.113520915885509d0     
  enrich%vhatR(id3) =    9.811114291935046d-002 
  enrich%vhatI(id3) =    4.495172617742562d-002 
  enrich%whatR(id3) =    5.485734634718119d-005
  enrich%whatI(id3) =    2.513407079403549d-005
 
  call enrich%wrapupTimeStep()

  print*, id1, id2, id3
  print*, maxval(abs(enrich%smallScales%u))
  print*, "max u difference:", p_maxval(maxval(abs(enrich%smallScales%u - u1)))
  print*, "max v difference:", p_maxval(maxval(abs(enrich%smallScales%v - v1)))
  print*, "max w difference:", p_maxval(maxval(abs(enrich%smallScales%wC - w1)))
  call assert(p_maxval(maxval(abs(enrich%smallScales%u - u1))) < tol)
  call assert(p_maxval(maxval(abs(enrich%smallScales%v - v1))) < tol)
  call assert(p_maxval(maxval(abs(enrich%smallScales%wC - w1))) < tol)
  call message("Test PASSED!")
   
  deallocate(u1,v1,w1) 
  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  call MPI_Finalize(ierr)
end program 
