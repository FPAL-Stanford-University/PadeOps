module timer
        use decomp_2d,       only: nrank 
        use mpi,             only: MPI_WTIME  
        use kind_parameters, only: rkind
        use reductions,      only: P_MAXVAL 
        implicit none
        
        double precision :: start, finish

        interface toc
            module procedure toc1,toc2,toc3,toc4
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
        end subroutine 

        subroutine toc2(message) 
            character(len=*), intent(in) :: message
            finish = MPI_WTIME()
            if (nrank == 0) then    
                print*, message, finish - start, " seconds"
            end if 
        end subroutine 

        subroutine toc3(message,val) 
            character(len=*), intent(in) :: message
            real(rkind), intent(out) :: val
            real(rkind) :: myval
            
            finish = MPI_WTIME()
            myval = real(finish - start,rkind)
            val = P_MAXVAL(myval)
            
            if (nrank == 0) then    
                print*, message, val, " seconds"
            end if 

        end subroutine 

        subroutine toc4(val) 
            real(rkind), intent(out) :: val
            real(rkind) :: myval

            finish = MPI_WTIME()
            myval = real(finish - start,rkind)
            val = P_MAXVAL(myval)
        end subroutine 

end module timer
