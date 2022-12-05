module auxiliary_openmp_subs
   use exits, only: warning
   implicit none

contains
   subroutine GetArguments(nthreads,nargs)
      use kind_parameters, only: clen
      integer, intent(out) :: nthreads
      integer, intent(in) :: nargs
      character(len=clen) :: args
      integer :: numargs

      numargs = command_argument_count()
      if (numargs < nargs) then
        call warning('User did not specify number of threads. Defaulting to 1.')
        nthreads = 1
      else
        call get_command_argument(nargs,args)
        read(args,*) nthreads
      end if

   end subroutine

end module 
