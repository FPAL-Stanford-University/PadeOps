module Sep1SolidEOS

    use kind_parameters, only: rkind
    use constants,       only: zero,half,one,two,three,fourth,sixth,third,twothird,fourthird
    use ElasticEOSMod,   only: elasticeos

    implicit none

    type, extends(elasticeos) :: sep1solid

        real(rkind) :: rho0 = one                    ! Reference density
        real(rkind) :: mu = zero                     ! Shear Modulus
        real(rkind) :: yield = one                   ! Yield Stress
        real(rkind) :: tau0 = 1.0d-10                ! Plastic relaxation time scale

    contains

        ! procedure :: init
        procedure :: get_finger
        procedure :: get_devstress
        procedure :: get_eelastic
        procedure :: get_sos
        procedure :: plastic_deformation

    end type

    interface sep1solid
        module procedure init
    end interface

contains

    function init(rho0_,mu_,yield_,tau0_) result(this)
        type(sep1solid) :: this
        real(rkind), intent(in) :: rho0_,mu_, yield_, tau0_

        this%rho0 = rho0_
        this%mu = mu_
        this%yield = yield_
        this%tau0 = tau0_

    end function

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

    subroutine get_eelastic(this,trG,trG2,detG,eelastic)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: trG,trG2,detG
        real(rkind), dimension(:,:,:), intent(out) :: eelastic

        print *, this%mu, this%rho0, maxval(detG), minval(detG)
        print *, maxval(trG), minval(trG)
        print *, maxval(trG2), minval(trG2)
        !eelastic = fourth*this%mu/this%rho0*(detG**(-twothird)*trG2 - two*detG**(-third)*trG + three)
        eelastic = fourth*this%mu/this%rho0!*(detG**(-twothird)*trG2 - two*detG**(-third)*trG + three)
        print *, '------99888'

    end subroutine

    pure subroutine get_sos(this,sos)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(inout) :: sos

        sos = sqrt(sos**two + fourthird*this%mu/this%rho0)

    end subroutine

    subroutine plastic_deformation(this, gfull)
        use constants, only: twothird
        use decomp_2d, only: nrank
        class(sep1solid), target, intent(in) :: this
        real(rkind), dimension(:,:,:,:), intent(inout) :: gfull
        real(rkind), dimension(3,3) :: g, u, vt, gradf
        real(rkind), dimension(3)   :: sval, beta, Sa, f, f1, f2, dbeta, beta_new
        real(rkind), dimension(15)  :: work
        real(rkind) :: sqrt_om, betasum, Sabymu_sq, ycrit, C0, t
        real(rkind) :: tol = real(1.D-12,rkind), residual
        integer :: i,j,k
        integer :: iters, niters = 500
        integer :: lwork, info
        integer, dimension(3) :: ipiv
        integer :: nxp, nyp, nzp

        nxp = size(gfull,1); nyp = size(gfull,2); nzp = size(gfull,3);

        ! Get optimal lwork
        lwork = -1
        call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, work, lwork, info)
        lwork = work(1)

        do k = 1,nzp
            do j = 1,nyp
                do i = 1,nxp
                    g(1,1) = gfull(i,j,k,1); g(1,2) = gfull(i,j,k,2); g(1,3) = gfull(i,j,k,3)
                    g(2,1) = gfull(i,j,k,4); g(2,2) = gfull(i,j,k,5); g(2,3) = gfull(i,j,k,6)
                    g(3,1) = gfull(i,j,k,7); g(3,2) = gfull(i,j,k,8); g(3,3) = gfull(i,j,k,9)

                    ! Get SVD of g
                    call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, work, lwork, info)
                    if(info .ne. 0) print '(A,I6,A)', 'proc ', nrank, ': Problem with SVD. Please check.'

                    sqrt_om = sval(1)*sval(2)*sval(3)
                    beta = sval**two / sqrt_om**(two/three)

                    betasum = sum( beta*(beta-one) ) / three
                    Sa = -this%mu*sqrt_om * ( beta*(beta-one) - betasum )

                    Sabymu_sq = sum(Sa**two) / this%mu**two
                    ycrit = Sabymu_sq - (two/three)*(this%yield/this%mu)**two

                    if (ycrit .LE. zero) then
                        ! print '(A)', 'Inconsistency in plastic algorithm, ycrit < 0!'
                        cycle
                    end if

                    C0 = Sabymu_sq / ycrit
                    Sa = Sa*( sqrt(C0 - one)/sqrt(C0) )

                    ! Now get new beta
                    f = Sa / (this%mu*sqrt_om); f(3) = beta(1)*beta(2)*beta(3)     ! New function value (target to attain)
                    
                    betasum = sum( beta*(beta-one) ) / three
                    f1 = -( beta*(beta-one) - betasum ); f1(3) = beta(1)*beta(2)*beta(3)   ! Original function value
                    residual = sqrt(sum( (f1-f)**two ))                                    ! Norm of the residual
                    iters = 0
                    do while ( (iters < niters) .AND. (residual .GT. tol) )
                        ! Get newton step
                        gradf(1,1) = -twothird*(two*beta(1)-one); gradf(1,2) =     third*(two*beta(2)-one); gradf(1,3) = third*(two*beta(3)-one)
                        gradf(2,1) =     third*(two*beta(1)-one); gradf(2,2) = -twothird*(two*beta(2)-one); gradf(2,3) = third*(two*beta(3)-one)
                        gradf(3,1) = beta(2)*beta(3);             gradf(3,2) = beta(3)*beta(1);             gradf(3,3) = beta(1)*beta(2)

                        dbeta = (f-f1)
                        call dgesv(3, 3, gradf, 3, ipiv, dbeta, 3, info)
                        
                        ! Backtracking line search
                        t = 1.
                        beta_new = beta + t * dbeta
                        betasum = sum( beta_new*(beta_new-one) ) / three
                        f2 = -( beta_new*(beta_new-one) - betasum ); f2(3) = beta_new(1)*beta_new(2)*beta_new(3)
                        do while ( sqrt(sum( (f2-f)**two )) .GE. residual )
                            t = half*t
                            beta_new = beta + t * dbeta
                            betasum = sum( beta_new*(beta_new-one) ) / three
                            f2 = -( beta_new*(beta_new-one) - betasum ); f2(3) = beta_new(1)*beta_new(2)*beta_new(3)
                        end do
                        beta = beta_new
                        f1 = f2
                        residual = sqrt(sum( (f1-f)**two ))                                    ! Norm of the residual
                    end do
                    if (iters >= niters) print '(A)', 'Newton solve in plastic_deformation did not converge'
                    
                    ! Then get new svals
                    sval = sqrt(beta) * sqrt_om**(one/three)
                    
                    ! Get g = u*sval*vt
                    vt(1,:) = vt(1,:)*sval(1); vt(2,:) = vt(2,:)*sval(2); vt(3,:) = vt(3,:)*sval(3)  ! sval*vt
                    g = MATMUL(u,vt) ! u*sval*vt

                    gfull(i,j,k,1) = g(1,1); gfull(i,j,k,2) = g(1,2); gfull(i,j,k,3) = g(1,3)
                    gfull(i,j,k,4) = g(2,1); gfull(i,j,k,5) = g(2,2); gfull(i,j,k,6) = g(2,3)
                    gfull(i,j,k,7) = g(3,1); gfull(i,j,k,8) = g(3,2); gfull(i,j,k,9) = g(3,3)
                end do
            end do
        end do
    end subroutine

end module
