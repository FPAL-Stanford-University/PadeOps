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

        real(rkind), dimension(:), allocatable :: svdwork     ! Work array for SVD stuff

    contains

        ! procedure :: init
        procedure :: get_finger
        procedure :: get_devstress
        procedure :: get_eelastic
        procedure :: get_sos
        procedure :: get_sos2
        procedure :: plastic_deformation
        procedure :: make_tensor_SPD
        final     :: destroy

    end type

    interface sep1solid
        module procedure init
    end interface

contains

    function init(rho0_,mu_,yield_,tau0_) result(this)
        type(sep1solid) :: this
        real(rkind), intent(in) :: rho0_,mu_, yield_, tau0_

        integer :: info, lwork
        real(rkind), dimension(3,3) :: g, u, vt
        real(rkind), dimension(3)   :: sval

        this%rho0 = rho0_
        this%mu = mu_
        this%yield = yield_
        this%tau0 = tau0_

        if (allocated(this%svdwork)) deallocate(this%svdwork); allocate(this%svdwork(1))

        ! Get optimal lwork
        lwork = -1
        call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, this%svdwork, lwork, info)
        lwork = this%svdwork(1)

        deallocate(this%svdwork); allocate(this%svdwork(lwork))

    end function

    pure elemental subroutine destroy(this)
        type(sep1solid), intent(inout) :: this

        if (allocated(this%svdwork)) deallocate(this%svdwork)
    end subroutine

    subroutine get_finger(this,g,finger,fingersq,trG,trG2,detG)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:,:), intent(in)  :: g
        real(rkind), dimension(:,:,:,:), intent(out) :: finger
        real(rkind), dimension(:,:,:,:), intent(out) :: fingersq
        real(rkind), dimension(:,:,:),   intent(out) :: trG, trG2, detG

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

            detG = GG11*(GG22*GG33-GG23*GG32) - GG12*(GG21*GG33-GG31*GG23) + GG13*(GG21*GG32-GG31*GG22)

            fingersq(:,:,:,1) = GG11*GG11 + GG12*GG21 + GG13*GG31
            fingersq(:,:,:,2) = GG11*GG12 + GG12*GG22 + GG13*GG32
            fingersq(:,:,:,3) = GG11*GG13 + GG12*GG23 + GG13*GG33
            fingersq(:,:,:,4) = GG21*GG12 + GG22*GG22 + GG23*GG32
            fingersq(:,:,:,5) = GG21*GG13 + GG22*GG23 + GG23*GG33
            fingersq(:,:,:,6) = GG31*GG13 + GG32*GG23 + GG33*GG33

            trG = GG11 + GG22 + GG33
            trG2 = fingersq(:,:,:,1) + fingersq(:,:,:,4) + fingersq(:,:,:,6)
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

        eelastic = fourth*this%mu/this%rho0*(detG**(-twothird)*trG2 - two*detG**(-third)*trG + three)

    end subroutine

    pure subroutine get_sos(this,rhom,sos)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in) :: rhom
        real(rkind), dimension(:,:,:), intent(inout) :: sos

        sos = sqrt(sos**two + fourthird*this%mu/rhom)

    end subroutine

    pure subroutine get_sos2(this,rhom,sos2)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in) :: rhom
        real(rkind), dimension(:,:,:), intent(inout) :: sos2

        sos2 = sos2 + fourthird*this%mu/rhom

    end subroutine

    subroutine plastic_deformation(this, gfull)
        use kind_parameters, only: clen
        use constants,       only: eps, twothird
        use decomp_2d,       only: nrank
        use exits,           only: GracefulExit, nancheck, message
        class(sep1solid), target, intent(inout) :: this
        real(rkind), dimension(:,:,:,:), intent(inout) :: gfull
        real(rkind), dimension(3,3) :: g, u, vt, gradf, gradf_new
        real(rkind), dimension(3)   :: sval, beta, Sa, f, f1, f2, dbeta, beta_new, dbeta_new
        real(rkind) :: sqrt_om, betasum, Sabymu_sq, ycrit, C0, t
        real(rkind) :: tol = real(1.D-12,rkind), residual, residual_new
        integer :: i,j,k
        integer :: iters
        integer, parameter :: niters = 500
        integer :: lwork, info
        integer, dimension(3) :: ipiv
        integer :: nxp, nyp, nzp
        character(len=clen) :: charout

        nxp = size(gfull,1); nyp = size(gfull,2); nzp = size(gfull,3);

        ! Get optimal lwork
        lwork = -1
        call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, this%svdwork, lwork, info)
        lwork = this%svdwork(1)
        if (lwork .GT. size(this%svdwork)) then
            deallocate(this%svdwork); allocate(this%svdwork(lwork))
        end if

        if ( nancheck(gfull) ) then
            call message("NaN found in g during plastic relaxation.")
        end if

        do k = 1,nzp
            do j = 1,nyp
                do i = 1,nxp
                    g(1,1) = gfull(i,j,k,1); g(1,2) = gfull(i,j,k,2); g(1,3) = gfull(i,j,k,3)
                    g(2,1) = gfull(i,j,k,4); g(2,2) = gfull(i,j,k,5); g(2,3) = gfull(i,j,k,6)
                    g(3,1) = gfull(i,j,k,7); g(3,2) = gfull(i,j,k,8); g(3,3) = gfull(i,j,k,9)

                    ! Get SVD of g
                    call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, this%svdwork, lwork, info)
                    if(info .ne. 0) then
                        write(charout, '(A,I0,A)') 'proc ', nrank, ': Problem with SVD. Please check.'
                        call GracefulExit(charout,3475)
                    end if

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

                    ! Get newton step
                    gradf(1,1) = -twothird*(two*beta(1)-one); gradf(1,2) =     third*(two*beta(2)-one); gradf(1,3) = third*(two*beta(3)-one)
                    gradf(2,1) =     third*(two*beta(1)-one); gradf(2,2) = -twothird*(two*beta(2)-one); gradf(2,3) = third*(two*beta(3)-one)
                    gradf(3,1) = beta(2)*beta(3);             gradf(3,2) = beta(3)*beta(1);             gradf(3,3) = beta(1)*beta(2)

                    dbeta = (f-f1)
                    call dgesv(3, 1, gradf, 3, ipiv, dbeta, 3, info)
                   
                    ! Compute residual
                    residual = -sum( (f1-f)*dbeta )                                    ! lambda**2
                    iters = 0
                    do while ( (iters < niters) .AND. (abs(residual) .GT. tol) )
                        ! Backtracking line search
                        t = 1._rkind
                        beta_new = beta + t * dbeta

                        ! Get new residual
                        gradf_new(1,1) = -twothird*(two*beta_new(1)-one);
                        gradf_new(1,2) =     third*(two*beta_new(2)-one);
                        gradf_new(1,3) = third*(two*beta_new(3)-one)

                        gradf_new(2,1) =     third*(two*beta_new(1)-one);
                        gradf_new(2,2) = -twothird*(two*beta_new(2)-one);
                        gradf_new(2,3) = third*(two*beta_new(3)-one)

                        gradf_new(3,1) = beta_new(2)*beta_new(3);
                        gradf_new(3,2) = beta_new(3)*beta_new(1);
                        gradf_new(3,3) = beta_new(1)*beta_new(2)

                        betasum = sum( beta_new*(beta_new-one) ) / three
                        f2 = -( beta_new*(beta_new-one) - betasum ); f2(3) = beta_new(1)*beta_new(2)*beta_new(3)
                        dbeta_new = (f-f2)
                        call dgesv(3, 1, gradf_new, 3, ipiv, dbeta_new, 3, info)
                        residual_new = -sum( (f2-f)*dbeta_new )                                    ! lambda**2

                        do while ( (abs(residual_new) .GE. abs(residual)) .AND. (t > eps) )
                            if (iters .GT. (niters - 10)) then
                                print '(A,I0,3(A,ES15.5))', 'iters = ', iters, ', t = ', t, ', residual_new = ', residual_new, ', residual = ', residual
                            end if

                            t = half*t
                            beta_new = beta + t * dbeta

                            gradf_new(1,1) = -twothird*(two*beta_new(1)-one);
                            gradf_new(1,2) =     third*(two*beta_new(2)-one);
                            gradf_new(1,3) = third*(two*beta_new(3)-one)

                            gradf_new(2,1) =     third*(two*beta_new(1)-one);
                            gradf_new(2,2) = -twothird*(two*beta_new(2)-one);
                            gradf_new(2,3) = third*(two*beta_new(3)-one)

                            gradf_new(3,1) = beta_new(2)*beta_new(3);
                            gradf_new(3,2) = beta_new(3)*beta_new(1);
                            gradf_new(3,3) = beta_new(1)*beta_new(2)

                            betasum = sum( beta_new*(beta_new-one) ) / three
                            f2 = -( beta_new*(beta_new-one) - betasum ); f2(3) = beta_new(1)*beta_new(2)*beta_new(3)

                            dbeta_new = (f-f2)
                            call dgesv(3, 1, gradf_new, 3, ipiv, dbeta_new, 3, info)
                            residual_new = -sum( (f2-f)*dbeta_new )                                    ! lambda**2
                        end do
                        beta = beta_new
                        f1 = f2
                        dbeta = dbeta_new
                        residual = residual_new

                        iters = iters + 1
                        if (t <= eps) then
                            print '(A)', 'Newton solve in plastic_deformation did not converge'
                            exit
                        end if
                    end do
                    if ((iters >= niters) .OR. (t <= eps)) then
                        write(charout,'(4(A,I0))') 'Newton solve in plastic_deformation did not converge at index ',i,',',j,',',k,' of process ',nrank
                        print '(A)', charout
                        print '(A)', 'g = '
                        print '(4X,3(ES15.5))', gfull(i,j,k,1), gfull(i,j,k,2), gfull(i,j,k,3)
                        print '(4X,3(ES15.5))', gfull(i,j,k,4), gfull(i,j,k,5), gfull(i,j,k,6)
                        print '(4X,3(ES15.5))', gfull(i,j,k,7), gfull(i,j,k,8), gfull(i,j,k,9)
                        print '(A,ES15.5)', '( ||S||^2 - (2/3) sigma_Y^2 )/mu^2 = ', ycrit

                        print '(A,ES15.5)', 'Relaxation, t = ', t
                        print '(A,ES15.5)', 'Residual = ', residual
                        call GracefulExit(charout,6382)
                    end if

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

    subroutine make_tensor_SPD(this,gfull)
        use exits, only: GracefulExit
        class(sep1solid), intent(inout) :: this
        real(rkind), dimension(:,:,:,:), intent(inout) :: gfull

        real(rkind), dimension(3,3) :: g, u, vt
        real(rkind), dimension(3)   :: sval
        integer :: nx, ny, nz
        integer :: i, j, k
        integer :: info, lwork

        if (size(gfull,4) .NE. 9) call GracefulExit("Incorrect dimension for tensor in make_tensor_SPD.",2384)

        nx = size(gfull,1); ny = size(gfull,2); nz = size(gfull,3)

        ! Get optimal lwork
        lwork = -1
        call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, this%svdwork, lwork, info)
        lwork = this%svdwork(1)
        if (lwork .GT. size(this%svdwork)) then
            deallocate(this%svdwork); allocate(this%svdwork(lwork))
        end if

        do k = 1,nz
            do j = 1,ny
                do i = 1,nx
                    g(1,1) = gfull(i,j,k,1); g(1,2) = gfull(i,j,k,2); g(1,3) = gfull(i,j,k,3)
                    g(2,1) = gfull(i,j,k,4); g(2,2) = gfull(i,j,k,5); g(2,3) = gfull(i,j,k,6)
                    g(3,1) = gfull(i,j,k,7); g(3,2) = gfull(i,j,k,8); g(3,3) = gfull(i,j,k,9)
        
                    ! Get SVD of g
                    call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, this%svdwork, lwork, info)

                    ! Get projection (V * Sigma * V^T)
                    u = transpose(vt)
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
