module auxiliary_openmp_subs
   implicit none

contains
   subroutine GetArguments(nthreads)
      use kind_parameters, only: clen
      integer, intent(out) :: nthreads
      character(len=clen) :: args
      integer :: numargs

      numargs = command_argument_count()
      if (numargs < 2) then
         print*, "Need to provide input-file and number of OMP threads as two separate command line arguments"
         stop
      end if

      call get_command_argument(2,args)

      read(args,*) nthreads

   end subroutine

end module 
