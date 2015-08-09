module timeStep
    use kind_parameters, only: rkind, mpirkind 
    use decomp_2d,       only: nrank
    use variables,       only: DT, CFL, FT, useCFL, fieldsPhys, fieldsSpec, buff_real

contains

subroutine getDT
    
    if (useCFL) then
        ! Compute the new time step 
    else
        return
    end if 

end subroutine

subroutine WrapUpTstep(tid)
    use mpi 
    integer, intent(in) :: tid
    real(rkind) :: mymax, globalmax
    integer :: ierr 

    ! First convert the updated fieldsSpec to fieldsPhys 
    call FT%ifft3_z2x(fieldsSpec(:,:,:,1),fieldsPhys(:,:,:,1))
    call FT%ifft3_z2x(fieldsSpec(:,:,:,2),fieldsPhys(:,:,:,2))
    call FT%ifft3_z2x(fieldsSpec(:,:,:,3),fieldsPhys(:,:,:,3))

    ! Do statistic, print shit to screen, etc. 
    
    ! Get tke 
    buff_Real = fieldsPhys(:,:,:,1)*fieldsPhys(:,:,:,1) + fieldsPhys(:,:,:,2)*fieldsPhys(:,:,:,2) + fieldsPhys(:,:,:,3)*fieldsPhys(:,:,:,3)
    mymax = maxval(abs(buff_real)) 
    call MPI_Reduce(mymax,globalmax, 1, mpirkind, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     
    if (nrank .eq. 0) then
        print*, "Time:", tid, "Max tke:", 0.5_rkind*globalmax
    end if

end subroutine 


end module 
