module timer
        use decomp_2d, only: nrank 
        use mpi,       only: MPI_WTIME  
        implicit none
        
        double precision :: start, finish

        interface toc
            module procedure toc1,toc2
        end interface

contains 
        subroutine tic 
            start = MPI_WTIME()
        end subroutine tic

        subroutine toc1 
            finish = MPI_WTIME()
            if (nrank == 0) then    
                print*, "Elapsed time is ",finish - start, " seconds"
            end if 
        end subroutine toc1

        subroutine toc2(message) 
            character(len=*), intent(in) :: message
            finish = MPI_WTIME()
            if (nrank == 0) then    
                print*, message, finish - start, " seconds"
            end if 
        end subroutine toc2

end module timer
