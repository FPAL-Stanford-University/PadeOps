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

#ifdef __bgq__
@process directive(ibm*)
#endif

module xlf_posix_bindings

#ifdef __bgq__
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

#endif
end module
