module test_plastic
    use kind_parameters, only: rkind
    use Sep1SolidEOS,    only: sep1solid
    implicit none

    integer, parameter :: nx = 1, ny = 1, nz = 1
    type(sep1solid) :: elastic
    real(rkind) :: mu = real(77.D9,rkind)
    real(rkind) :: yield = real(2.49D9,rkind)

contains

    subroutine plastic_deformation(g11,g12,g13,g21,g22,g23,g31,g32,g33)
        use constants
        use decomp_2d, only: nrank
        real(rkind), dimension(nx,ny,nz), intent(inout) :: g11,g12,g13,g21,g22,g23,g31,g32,g33
        real(rkind), dimension(3,3) :: g, u, vt, gradf
        real(rkind), dimension(3)   :: sval, beta, Sa, f, f1, f2, dbeta, beta_new
        real(rkind), dimension(15)  :: work
        real(rkind) :: sqrt_om, betasum, Sabymu_sq, ycrit, C0, t
        real(rkind) :: tol = real(1.D-12,rkind), residual
        integer :: i,j,k
        integer :: iters, niters = 500
        integer :: info, lwork
        integer, dimension(3) :: ipiv

        ! Get optimal lwork
        lwork = -1
        call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, work, lwork, info)
        lwork = work(1)

        do k = 1,nz
            do j = 1,ny
                do i = 1,nx
                    g(1,1) = g11(i,j,k); g(1,2) = g12(i,j,k); g(1,3) = g13(i,j,k)
                    g(2,1) = g21(i,j,k); g(2,2) = g22(i,j,k); g(2,3) = g23(i,j,k)
                    g(3,1) = g31(i,j,k); g(3,2) = g32(i,j,k); g(3,3) = g33(i,j,k)

                    ! Get SVD of g
                    call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, work, lwork, info)
                    if(info .ne. 0) print '(A,I6,A)', 'proc ', nrank, ': Problem with SVD. Please check.'

                    sqrt_om = sval(1)*sval(2)*sval(3)
                    beta = sval**two / sqrt_om**(two/three)

                    betasum = sum( beta*(beta-one) ) / three
                    Sa = -mu*sqrt_om * ( beta*(beta-one) - betasum )

                    Sabymu_sq = sum(Sa**two) / mu**two
                    ycrit = Sabymu_sq - (two/three)*(yield/mu)**two

                    if (ycrit .LE. zero) then
                        print '(A)', 'Inconsistency in plastic algorithm, ycrit < 0!'
                        cycle
                    end if

                    C0 = Sabymu_sq / ycrit
                    Sa = Sa*( sqrt(C0 - one)/sqrt(C0) )

                    ! Now get new beta
                    f = Sa / (mu*sqrt_om); f(3) = one     ! New function value (target to attain)
                    
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
                            t = 0.9_rkind*t
                            beta_new = beta + t * dbeta
                            betasum = sum( beta_new*(beta_new-one) ) / three
                            f2 = -( beta_new*(beta_new-one) - betasum ); f2(3) = beta_new(1)*beta_new(2)*beta_new(3)
                        end do
                        beta = beta_new
                        f1 = f2
                        residual = sqrt(sum( (f1-f)**two ))                                    ! Norm of the residual
                        iters = iters + 1
                    end do
                    if (iters >= niters) print '(A)', 'Newton solve in plastic_deformation did not converge'
                    
                    print*, "residual = ", residual, " in ", iters, "iterations"

                    ! Then get new svals
                    sval = sqrt(beta) * sqrt_om**(one/three)
                    
                    ! Get g = u*sval*vt
                    vt(1,:) = vt(1,:)*sval(1); vt(2,:) = vt(2,:)*sval(2); vt(3,:) = vt(3,:)*sval(3)  ! sval*vt
                    g = MATMUL(u,vt) ! u*sval*vt

                    g11(i,j,k) = g(1,1); g12(i,j,k) = g(1,2); g13(i,j,k) = g(1,3)
                    g21(i,j,k) = g(2,1); g22(i,j,k) = g(2,2); g23(i,j,k) = g(2,3)
                    g31(i,j,k) = g(3,1); g32(i,j,k) = g(3,2); g33(i,j,k) = g(3,3)
                end do
            end do
        end do
    end subroutine

    subroutine get_SS(g,SS)
        use constants
        real(rkind), dimension(nx,ny,nz,9), intent(in) :: g
        real(rkind), dimension(nx,ny,nz), intent(out)  :: SS
        real(rkind), dimension(nx,ny,nz,6) :: finger,fingersq,devstress
        real(rkind), dimension(nx,ny,nz) :: trG, trG2, detG

        call elastic%get_finger(g,finger,fingersq,trG,trG2,detG)
        call elastic%get_devstress(finger,fingersq,trG,trG2,detG,devstress)

        SS = devstress(:,:,:,1)**2 + devstress(:,:,:,4)**2 + devstress(:,:,:,6)**2 + two*(devstress(:,:,:,2)**2 + devstress(:,:,:,3)**2 + devstress(:,:,:,5)**2)

    end subroutine

end module

program test_plastic_deformation
    use kind_parameters, only: rkind
    use constants,       only: one, two, three
    use test_plastic
    implicit none

    real(rkind), dimension(nx,ny,nz,9), target :: g
    real(rkind), dimension(:,:,:), pointer :: g11,g12,g13,g21,g22,g23,g31,g32,g33
    real(rkind), dimension(nx,ny,nz) :: SS

    call elastic%init(mu,yield)

    g11 => g(:,:,:,1); g12 => g(:,:,:,2); g13 => g(:,:,:,3)
    g21 => g(:,:,:,4); g22 => g(:,:,:,5); g23 => g(:,:,:,6)
    g31 => g(:,:,:,7); g32 => g(:,:,:,8); g33 => g(:,:,:,9)

    g11 = real(1.0D0 ,rkind); g12 = real(1.0D-3,rkind); g13 = real(0.0D0 ,rkind)
    g21 = real(2.0D-3,rkind); g22 = real(0.0D0 ,rkind); g23 = real(0.0D0 ,rkind)
    g31 = real(0.0D0 ,rkind); g32 = real(0.0D0 ,rkind); g33 = real(0.0D0 ,rkind)

    g11 = g11 + one; g22 = g22 + one; g33 = g33 + one

    call get_SS(g,SS)
    write(*,*) "Initial g:"
    write(*,*) g11(1,1,1), g12(1,1,1), g13(1,1,1)
    write(*,*) g21(1,1,1), g22(1,1,1), g23(1,1,1)
    write(*,*) g31(1,1,1), g32(1,1,1), g33(1,1,1)
    write(*,*) "Initial SS/mu - (2/3)(yield/mu)^2", SS(1,1,1)/mu**2 - (two/three)*(yield/mu)**2

    call plastic_deformation(g11,g12,g13,g21,g22,g23,g31,g32,g33)
    ! call elastic%plastic_deformation(g)

    call get_SS(g,SS)
    write(*,*) "Final g:"
    write(*,*) g11(1,1,1), g12(1,1,1), g13(1,1,1)
    write(*,*) g21(1,1,1), g22(1,1,1), g23(1,1,1)
    write(*,*) g31(1,1,1), g32(1,1,1), g33(1,1,1)
    write(*,*) "Final SS/mu - (2/3)(yield/mu)^2", SS(1,1,1)/mu**2 - (two/three)*(yield/mu)**2

end program
