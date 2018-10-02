#ifdef __bgq__

! Copyright IBM Corp. 2011.
! All Rights Reserved.
!
! For XLF 12.1 and above, compile with:
!   xlf -qnoobject xlf_posix_bindings.F03
! If compiling with XLF 11.1, also add -qxflag=ignore_tkr
! The command above will generate xlf_posix_bindings.mod
! To use, add:
!   use xlf_posix_bindings
! to your compilation unit and ensure xlf_posix_bindings.mod is in the
! include path.  (e.g. by use -I<directory containing xlf_posix_bindings.mod>
!
! This module defines XLF bindings for the following POSIX functions:
!   alarm,
!   calloc,
!   clock,
!   errno,  (even on systems where errno is a macro)
!   exit,
!   free,
!   getgid,
!   getpid,
!   getuid,
!   malloc,
!   qsort,
!   realloc,
!   sleep,
!   time,   (tloc argument should be loaded via C_LOC from iso_c_binding)
!   umask,
!   usleep
!
! The module is designed to work with programs not compiled with -qextname,
! but will work with -qextname as well. (although it's only necessary for
! programs compiled without -qextname).
!
@process directive(ibm*)
module xlf_posix_bindings
  use, intrinsic :: iso_c_binding
  implicit none
  private

  interface
    function alarm(seconds) bind(c, name="alarm")
      integer(4), value :: seconds
      integer(4) alarm
    end function
  end interface
  public alarm

  interface
    function calloc(nelem, elsize) bind(c, name="calloc")
      import C_SIZE_T, C_INTPTR_T
      integer(C_SIZE_T), value :: nelem, elsize
      integer(C_INTPTR_T) calloc
    end function
  end interface
  public calloc

  interface
    pure function clock() bind(c, name="clock")
#if __linux__
      import C_LONG
      integer(C_LONG) clock
#else
      integer(4) clock
#endif
    end function
  end interface
  public clock

  interface errno
    function ierrno_() bind(c, name="ierrno_")
      integer(4) ierrno_
    end function
  end interface
  public errno

  interface exit
    subroutine exit_(status) bind(c, name="exit_")
      integer(4) status
    end subroutine
  end interface
  public exit

  interface
    subroutine free(ptr) bind(c, name="free")
      import C_INTPTR_T
      integer(C_INTPTR_T), value :: ptr
    end subroutine
  end interface
  public free

  interface
    pure function getgid() bind(c, name="getgid")
      integer(4) getgid
    end function
  end interface
  public getgid

  interface
    pure function getpid() bind(c, name="getpid")
      integer(4) getpid
    end function
  end interface
  public getpid

  interface
    pure function getuid() bind(c, name="getuid")
      integer(4) getuid
    end function
  end interface
  public getuid

  interface
    function malloc(size) bind(c, name="malloc")
      import C_SIZE_T, C_INTPTR_T
      integer(C_SIZE_T), value :: size
      integer(C_INTPTR_T) malloc
    end function
  end interface
  public malloc

  abstract interface
    integer(4) function compar_iface(a, b)
      integer, intent(in) :: a, b
! Until we implement TYPE(*), use ignore_tkr
!ibm* ignore_tkr a, b
    end function
  end interface

  interface
    subroutine qsort(base, nel, width, compar) bind(c, name="qsort")
      import C_SIZE_T, compar_iface
      integer :: base
! Until we implement TYPE(*), use ignore_tkr
!ibm* ignore_tkr base
      integer(C_SIZE_T), value :: nel, width
      procedure(compar_iface) compar   
    end subroutine
  end interface
  public qsort

  interface
    function realloc(ptr, size) bind(c, name="realloc")
      import C_SIZE_T, C_INTPTR_T
      integer(C_INTPTR_T), value :: ptr
      integer(C_SIZE_T), value :: size
      integer(C_INTPTR_T) realloc
    end function
  end interface
  public realloc

  interface
    function sleep(seconds) bind(c, name="sleep")
      integer(4), value :: seconds
      integer(4) :: sleep
    end function
  end interface
  public sleep

  interface
    function time(tloc) bind(c, name="time")
      import C_PTR, C_LONG
      type(C_PTR), value :: tloc
      integer(C_LONG) time
    end function
  end interface
  public time

  interface
    function umask(cmask) bind(c, name="cmask")
      integer(4), value :: cmask
      integer(4) umask
    end function
  end interface
  public umask

  interface
    function usleep(useconds) bind(c, name="usleep")
      integer(4), value :: useconds
      integer(4) usleep
    end function
  end interface
  public usleep
end module

#endif

module random
    use kind_parameters, only: rkind
    use constants, only: zero, one, two, pi
   
    implicit none
    private 
    public :: uniform_random, gaussian_random

    interface uniform_random
        module procedure unrand3R, unrand2R, unrand1R, unrand0R
    end interface 
    interface gaussian_random
        module procedure grand3R, grand2R, grand1R
    end interface 
contains

    subroutine grand3R(array,mu,sigma,seed)
        real(rkind), dimension(:,:,:), intent(inout) :: array
        real(rkind), intent(in) :: mu, sigma
        integer, intent(in), optional :: seed
        real(rkind), dimension(:,:,:), allocatable, target :: uarr
        real(rkind), dimension(:,:,:), pointer :: uarr1, uarr2      
 
        allocate(uarr(size(array,1),size(array,2),2*size(array,3))) 
        
        if (present(seed)) then
            call uniform_random(uarr,zero,one,seed)
        else
            call uniform_random(uarr,zero,one)
        end if 
        uarr1 => uarr(:,:,1:size(array,3)) 
        uarr2 => uarr(:,:,size(array,3)+1:2*size(array,3)) 
       
        ! Now just use the Box-Muller algorithm 
        array = sqrt(-two*log(uarr1))*cos(two*pi*uarr2)        

        array = mu + sigma*array
        nullify(uarr1, uarr2)
        deallocate(uarr)
    end subroutine


    subroutine grand2R(array,mu,sigma,seed)
        real(rkind), dimension(:,:), intent(inout) :: array
        real(rkind), intent(in) :: mu, sigma
        integer, intent(in), optional :: seed
        
        real(rkind), dimension(:,:), allocatable, target :: uarr
        real(rkind), dimension(:,:), pointer :: uarr1, uarr2      
 
        allocate(uarr(size(array,1),2*size(array,2))) 
        
        if (present(seed)) then
            call uniform_random(uarr,zero,one,seed)
        else
            call uniform_random(uarr,zero,one)
        end if 
        
        
        uarr1 => uarr(:,1:size(array,2)) 
        uarr2 => uarr(:,size(array,2)+1:2*size(array,2)) 
       
        ! Now just use the Box-Muller algorithm 
        array = sqrt(-two*log(uarr1))*cos(two*pi*uarr2)        

        array = mu + sigma*array
        nullify(uarr1, uarr2)
        nullify(uarr1, uarr2)
        deallocate(uarr)
    end subroutine

    subroutine grand1R(array,mu,sigma,seed)
        real(rkind), dimension(:), intent(inout) :: array
        real(rkind), intent(in) :: mu, sigma
        integer, intent(in), optional :: seed
        
        real(rkind), dimension(:), allocatable, target :: uarr
        real(rkind), dimension(:), pointer :: uarr1, uarr2      
 
        allocate(uarr(2*size(array,1))) 
        
        if (present(seed)) then
            call uniform_random(uarr,zero,one,seed)
        else
            call uniform_random(uarr,zero,one)
        end if 
        
        uarr1 => uarr(1:size(array,1)) 
        uarr2 => uarr(size(array,1)+1:2*size(array,1)) 
       
        ! Now just use the Box-Muller algorithm 
        array = sqrt(-two*log(uarr1))*cos(two*pi*uarr2)        

        array = mu + sigma*array
        nullify(uarr1, uarr2)
        nullify(uarr1, uarr2)
        deallocate(uarr)
    end subroutine
    
    subroutine unrand3R(array,left, right, seed)
        real(rkind), dimension(:,:,:), intent(inout) :: array
        real(rkind), intent(in) :: left, right
        integer, optional, intent(in) :: seed
        real(rkind) :: diff
        integer, allocatable :: iseed(:)
        integer :: n

        diff = right - left
   
        if (present(seed)) then
            call random_seed(size = n)
            allocate(iseed(n))
            iseed = seed
            call random_seed(put=iseed)
            deallocate(iseed)
        else
            call init_random_seed()
        end if 
        call random_number(array)
        array = diff*array
        array = array + left
    end subroutine

    subroutine unrand2R(array,left, right, seed)
        real(rkind), dimension(:,:), intent(inout) :: array
        real(rkind), intent(in) :: left, right
        real(rkind) :: diff
        integer, optional, intent(in) :: seed
        integer, allocatable :: iseed(:)
        integer :: n

        diff = right - left
    
        if (present(seed)) then
            call random_seed(size = n)
            allocate(iseed(n))
            iseed = seed
            call random_seed(put=iseed)
            deallocate(iseed)
        else
            call init_random_seed()
        end if 
        call random_number(array)
        array = diff*array
        array = array + left
    end subroutine

    subroutine unrand1R(array,left, right, seed)
        real(rkind), dimension(:), intent(inout) :: array
        real(rkind), intent(in) :: left, right
        real(rkind) :: diff
        integer, optional, intent(in) :: seed
        integer, allocatable :: iseed(:)
        integer :: n

        diff = right - left
    
        if (present(seed)) then
            call random_seed(size = n)
            allocate(iseed(n))
            iseed = seed
            call random_seed(put=iseed)
            deallocate(iseed)
        else
            call init_random_seed()
        end if 
        call random_number(array)
        array = diff*array
        array = array + left
    end subroutine

    subroutine unrand0R(val,left, right)
        real(rkind), intent(inout) :: val
        real(rkind), intent(in) :: left, right
        real(rkind) :: diff

        diff = right - left
    
        call init_random_seed()
        call random_number(val)
        val = diff*val
        val = val + left
    end subroutine


    subroutine init_random_seed()
        ! Taken from GNU 
        use iso_fortran_env, only: int64
#ifdef __INTEL_COMPILER
        use ifport
#endif
#ifdef __bgq__
        use xlf_posix_bindings
#endif
        implicit none
        integer, allocatable :: iseed(:)
        integer :: i, n, un, istat, dt(8), pid
        integer(int64) :: t
        
        call random_seed(size = n)
        allocate(iseed(n))
        ! First try if the OS provides a random number generator
        open(newunit=un, file="/dev/urandom", access="stream", &
             form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
           read(un) iseed
           close(un)
        else
           ! Fallback to XOR:ing the current time and pid. The PID is
           ! useful in case one launches multiple instances of the same
           ! program in parallel.
           call system_clock(t)
           if (t == 0) then
              call date_and_time(values=dt)
              t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                   + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                   + dt(3) * 24_int64 * 60 * 60 * 1000 &
                   + dt(5) * 60 * 60 * 1000 &
                   + dt(6) * 60 * 1000 + dt(7) * 1000 &
                   + dt(8)
           end if
           pid = getpid()
           t = ieor(t, int(pid, kind(t)))
           do i = 1, n
              iseed(i) = lcg(t)
           end do
        end if
        call random_seed(put=iseed)
    contains
        ! This simple PRNG might not be good enough for real work, but is
        ! sufficient for seeding a better PRNG.
        function lcg(s)
          integer :: lcg
          integer(int64) :: s
          if (s == 0) then
             s = 104729
          else
             s = mod(s, 4294967296_int64)
          end if
          s = mod(s * 279470273_int64, 4294967291_int64)
          lcg = int(mod(s, int(huge(0), int64)), kind(0))
        end function lcg
    end subroutine init_random_seed

end module 
