module exits

    ! use mpi_module, only : taskid
    implicit none
    private
    public :: GracefulExit, parallelRun
    logical :: parallelRun = .false. 
     
contains

    subroutine gracefulExit(message, errcode)
        character(len=*), intent(in) :: message
        integer, intent(in) :: errcode, rank

        if (parallelRun)
            print*, "Message from task: ", 0! taskid
            print*, message
            print*, "Error Code:", errcode
            stop
            !call mpi_abort(?,?)
        else
            print*, "Message from task: ", 0
            print*, message
            print*, "Error Code:", errcode
            stop
        end if 
    end subroutine
    
end module 
