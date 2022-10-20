module auxiliary_openmp_subs
   implicit none

contains
   subroutine GetArguments(nthreads)
      use kind_parameters, only: clen
      integer, intent(out) :: nthreads
      character(len=clen) :: args
      integer :: numargs

      numargs = command_argument_count()
      if (numargs < 4) then
         print*, "Need to provide 3 input-files and the number of OMP threads as four separate command line arguments"
         print*, "i.e. <large scales input file> <small scales input file> <Gabor modes input file> <nthread>"
         stop
      end if

      call get_command_argument(4,args)

      read(args,*) nthreads

   end subroutine

end module 
