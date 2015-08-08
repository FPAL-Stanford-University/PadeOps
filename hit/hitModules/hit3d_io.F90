module hit3d_io
    use decomp_2d, only: nrank, nproc 
    use kind_parameters, only: rkind
    use fields, only: Fields_phys
                
    use inputParameters, only: inputDir 

contains
    subroutine getHit3d_uvw
        character(len=64) :: uFile = "u.000000"
        character(len=64) :: vFile = "v.000000"
        character(len=64) :: wFile = "w.000000"


    end subroutine 
 

end module hit3d_io
