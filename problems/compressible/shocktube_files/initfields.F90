subroutine initfields(decomp,dx,dy,dz,inpDirectory,mesh,fields)
    use kind_parameters,  only: rkind
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index
    use decomp_2d,          only: decomp_info
    
    implicit none
    character(len=*),                                               intent(in)    :: inpDirectory
    type(decomp_info),                                              intent(in)    :: decomp
    real(rkind),                                                    intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:),     intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp

    real(rkind) :: rhoL = one
    real(rkind) :: rhoR = one/eight

    real(rkind) :: pL = one
    real(rkind) :: pR = 0.1_rkind

    associate( rho => fields(:,:,:,rho_index), u => fields(:,:,:,u_index), &
                 v => fields(:,:,:,  v_index), w => fields(:,:,:,w_index), &
                 p => fields(:,:,:,  p_index), T => fields(:,:,:,T_index), &
                 e => fields(:,:,:,  e_index),                             &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        tmp = half * ( one+tanh(x/(dx)) )

        rho = (one-tmp)*rhoL + tmp*rhoR
        u   = zero
        v   = zero
        w   = zero
        p   = (one-tmp)*pL + tmp*pR
           
    end associate

end subroutine

