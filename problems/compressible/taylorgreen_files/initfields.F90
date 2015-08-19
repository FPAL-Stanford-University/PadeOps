subroutine initfields(decomp,dx,dy,dz,inpDirectory,mesh,fields)
    use kind_parameters,  only: rkind
    use constants,        only: zero,one,two,pi
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index
    use decomp_2d,          only: decomp_info
    
    implicit none
    character(len=*),                                               intent(in)    :: inpDirectory
    type(decomp_info),                                              intent(in)    :: decomp
    real(rkind),                                                    intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:),     intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    associate( rho => fields(:,:,:,rho_index), u => fields(:,:,:,u_index), &
                 v => fields(:,:,:,  v_index), w => fields(:,:,:,w_index), &
                 p => fields(:,:,:,  p_index), T => fields(:,:,:,T_index), &
                 e => fields(:,:,:,  e_index),                             &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        rho = one
        u   =  sin(x)*cos(y)*cos(z)
        v   = -cos(x)*sin(y)*cos(z)
        w   = zero
        p   = 100._rkind + ( (cos(two*z) + two)*(cos(two*x) + cos(two*y)) - two ) / 16._rkind
           
    end associate

end subroutine

