module timer
        use decomp_2d,       only: nrank 
        use mpi,             only: MPI_WTIME  
        use kind_parameters, only: rkind, clen
        use reductions,      only: P_MAXVAL 
        implicit none

        !external :: MPI_BARRIER

        double precision :: start, finish

        interface toc
            module procedure toc1,toc2,toc3,toc4, tocall, toc5
        end interface

contains 

        subroutine tocall(comm, val)
            integer, intent(in) :: comm
            integer :: ierr
            real(rkind), intent(out) :: val
            real(rkind) :: myval
            
            call mpi_barrier(comm, ierr)
            finish = MPI_WTIME()
            myval = real(finish - start,rkind)
            val = P_MAXVAL(myval)
            
        end subroutine 

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

        subroutine toc5(Nnow, Ntotal, whichMssg)
            integer, intent(in) :: Nnow, Ntotal
            character(len=clen) :: message
            integer, intent(in), optional :: whichMssg
            integer :: mssgID = 1

            if (present(whichMssg)) mssgID = whichMssg

            finish = MPI_WTIME()
            if (nrank == 0) then   
                if (mssgID == 1) then 
                  write(message,'(I0,A,I0,A,F0.3,A)') Nnow, ' of ', Ntotal, &
                    ' took ', finish - start, ' seconds'
                else
                  write(message,'(F0.3,A,F0.3,A)') real(Nnow)/real(Ntotal)*100.d0, '% complete '// &
                    'and took ', finish - start, ' seconds'
                end if
                print*, trim(message)
            end if 

        end subroutine
end module timer
