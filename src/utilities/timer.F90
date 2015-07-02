module timer
        implicit none
        real(kind=4) :: start, finish
        real(kind=4) :: TARRAY(2), TOE

        interface toc
            module procedure toc1,toc2
        end interface

contains 
        subroutine tic 
                CALL ETIME(TARRAY, TOE)
                start = TARRAY(1)
        end subroutine tic

        subroutine toc1 
                CALL ETIME(TARRAY, TOE)
                finish = TARRAY(1)
                print*, "Elapsed time is ",finish - start, " seconds"
        end subroutine toc1

        subroutine toc2(message) 
                character(len=*), intent(in) :: message
                CALL ETIME(TARRAY, TOE)
                finish = TARRAY(1)
                print*, message, finish - start, " seconds"
        end subroutine toc2

end module timer
