module Sep1SolidEOS

    use kind_parameters, only: rkind
    use constants,       only: zero,half,one,two,three,fourth,sixth,third,twothird,fourthird
    use ElasticEOSMod,   only: elasticeos

    implicit none

    type, extends(elasticeos) :: sep1solid

        real(rkind) :: mu = zero                     ! Shear Modulus

    contains

        procedure :: init
        procedure :: get_finger
        procedure :: get_devstress
        procedure :: get_eelastic
        procedure :: get_sos

    end type

contains

    subroutine init(this,mu_)
        class(sep1solid), intent(inout) :: this
        real(rkind) :: mu_

        this%mu = mu_

    end subroutine

    pure subroutine get_finger(this,g,finger,fingersq,trG,trG2,detG)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:,:), intent(in)  :: g
        real(rkind), dimension(:,:,:,:), intent(out) :: finger
        real(rkind), dimension(:,:,:,:), intent(out), optional :: fingersq
        real(rkind), dimension(:,:,:),   intent(out), optional :: trG, trG2, detG

        associate ( g11 => g(:,:,:,1), g12 => g(:,:,:,2), g13 => g(:,:,:,3), &
                    g21 => g(:,:,:,4), g22 => g(:,:,:,5), g23 => g(:,:,:,6), &
                    g31 => g(:,:,:,7), g32 => g(:,:,:,8), g33 => g(:,:,:,9)  )

            finger(:,:,:,1) = g11*g11 + g21*g21 + g31*g31
            finger(:,:,:,2) = g11*g12 + g21*g22 + g31*g32
            finger(:,:,:,3) = g11*g13 + g21*g23 + g31*g33
            finger(:,:,:,4) = g12*g12 + g22*g22 + g32*g32
            finger(:,:,:,5) = g12*g13 + g22*g23 + g32*g33
            finger(:,:,:,6) = g13*g13 + g23*g23 + g33*g33

        end associate

        associate ( GG11 => finger(:,:,:,1), GG12 => finger(:,:,:,2), GG13 => finger(:,:,:,3), &
                    GG21 => finger(:,:,:,2), GG22 => finger(:,:,:,4), GG23 => finger(:,:,:,5), &
                    GG31 => finger(:,:,:,3), GG32 => finger(:,:,:,5), GG33 => finger(:,:,:,6)  )

            if(present(detG)) then
                detG = GG11*(GG22*GG33-GG23*GG32) - GG12*(GG21*GG33-GG31*GG23) + GG13*(GG21*GG32-GG31*GG22)
            endif

            if(present(fingersq)) then
                fingersq(:,:,:,1) = GG11*GG11 + GG12*GG21 + GG13*GG31
                fingersq(:,:,:,2) = GG11*GG12 + GG12*GG22 + GG13*GG32
                fingersq(:,:,:,3) = GG11*GG13 + GG12*GG23 + GG13*GG33
                fingersq(:,:,:,4) = GG21*GG12 + GG22*GG22 + GG23*GG32
                fingersq(:,:,:,5) = GG21*GG13 + GG22*GG23 + GG23*GG33
                fingersq(:,:,:,6) = GG31*GG13 + GG32*GG23 + GG33*GG33
            endif

            if(present(trG)) trG = GG11 + GG22 + GG33
            if(present(trG2) .and. present(fingersq)) trG2 = fingersq(:,:,:,1) + fingersq(:,:,:,4) + fingersq(:,:,:,6)
        end associate

    end subroutine

    pure subroutine get_devstress(this,finger,fingersq,trG,trG2,detG,devstress)
        use exits, only: GracefulExit
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:,:), intent(in)  :: finger
        real(rkind), dimension(:,:,:,:), intent(in)  :: fingersq
        real(rkind), dimension(:,:,:),   intent(in)  :: trG, trG2, detG
        real(rkind), dimension(:,:,:,:), intent(out) :: devstress

        real(rkind), dimension(size(finger,1),size(finger,2),size(finger,3)) :: devstmp
        integer :: i

        ! if(.not.present(fingersq)) call GracefulExit("fingersq required for devstress",1111)
        
        do i = 1,6
            devstress(:,:,:,i) = -this%mu*(detG**(-sixth)*fingersq(:,:,:,i) - detG**sixth*finger(:,:,:,i))
        end do
        devstmp = third*this%mu*(detG**(-sixth)*trG2 - detG**sixth*trG)
        devstress(:,:,:,1) = devstress(:,:,:,1) + devstmp
        devstress(:,:,:,4) = devstress(:,:,:,4) + devstmp
        devstress(:,:,:,6) = devstress(:,:,:,6) + devstmp

    end subroutine

    pure subroutine get_eelastic(this,rho0,trG,trG2,detG,eelastic)
        class(sep1solid), intent(in) :: this
        real(rkind),                   intent(in)  :: rho0
        real(rkind), dimension(:,:,:), intent(in)  :: trG,trG2,detG
        real(rkind), dimension(:,:,:), intent(out) :: eelastic

        eelastic = fourth*this%mu/rho0*(detG**(-twothird)*trG2 - two*detG**(-third)*trG + three)

    end subroutine

    pure subroutine get_sos(this,rho0,sos)
        class(sep1solid), intent(in) :: this
        real(rkind),                   intent(in)  :: rho0
        real(rkind), dimension(:,:,:), intent(inout) :: sos

        sos = sqrt(sos**two + fourthird*this%mu/rho0)

    end subroutine

end module
