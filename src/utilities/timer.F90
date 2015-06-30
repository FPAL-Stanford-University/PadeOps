module timer
        implicit none
        real(kind=4) :: start, finish
        real(kind=4) :: TARRAY(2), TOE

contains 
        subroutine tic 
                CALL ETIME(TARRAY, TOE)
                start = TARRAY(1)
        end subroutine tic

        subroutine toc 
                CALL ETIME(TARRAY, TOE)
                finish = TARRAY(1)
                print*, "Elapsed time is ",finish - start, " seconds"
        end subroutine toc

end module timer
