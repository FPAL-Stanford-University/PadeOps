module test_plastic
    use kind_parameters, only: rkind
    implicit none

    integer, parameter :: nx = 1, ny = 1, nz = 1
    type(sep1solid) :: elastic
    real(rkind) :: mu = real(1.D8,rkind)
    real(rkind) :: yield = real(1.D8,rkind)

contains

    subroutine plastic_deformation(g11,g12,g13,g21,g22,g23,g31,g32,g33)
        real(rkind), dimension(nx,ny,nz), intent(inout) :: g11,g12,g13,g21,g22,g23,g31,g32,g33
        real(rkind), dimension(3,3) :: g, u, vt, gradf
        real(rkind), dimension(3)   :: sval, beta, Sa, f, f1, f2, dbeta, beta_new
        real(rkind), dimension(15)  :: work
        real(rkind) :: sqrt_om, betasum, Sabymu_sq, ycrit, C0, t
        real(rkind) :: tol = real(1.D-8,rkind)
        integer :: i,j,k
        integer :: iters, niters = 500
        integer :: info
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
                    if(info .ne. 0) print '(A,I,A)', 'proc ', nrank, ': Problem with SVD. Please check.'

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
                    f = Sa / (mu*sqrt_om); f(3) = beta(1)*beta(2)*beta(3)     ! New function value (target to attain)
                    
                    betasum = sum( beta*(beta-one) ) / three
                    f1 = -( beta*(beta-one) - betasum ); f1(3) = beta(1)*beta(2)*beta(3)   ! Original function value
                    residual = sqrt(sum( (f1-f)**two ))                                    ! Norm of the residual
                    iters = 0
                    do while ( (iters < niters) .AND. (residual .GT. tol) )
                        ! Get newton step
                        gradf(1,1) = -twothird*(two*beta(1)-one); gradf(1,2) =     third*(two*beta(2)-one); gradf(1,3) = third*(two*beta(3)-one)
                        gradf(2,1) =     third*(two*beta(1)-one); gradf(2,2) = -twothird*(two*beta(2)-onw); gradf(2,3) = third*(two*beta(3)-one)
                        gradf(3,1) = beta(2)*beta(3);             gradf(3,2) = beta(3)*beta(1);             gradf(3,3) = beta(1)*beta(2)

                        dbeta = (f-f1)
                        call dgesv(3, 3, gradf, 3, ipiv, dbeta, 1, info)
                        
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

                    g11(i,j,k) = g(1,1); g12(i,j,k) = g(1,2); g13(i,j,k) = g(1,3)
                    g21(i,j,k) = g(2,1); g22(i,j,k) = g(2,2); g23(i,j,k) = g(2,3)
                    g31(i,j,k) = g(3,1); g32(i,j,k) = g(3,2); g33(i,j,k) = g(3,3)
                end do
            end do
        end do
    end subroutine

    subroutine get_SS(g,SS)
        real(rkind), dimension(nx,ny,nz,9), intent(in) :: g
        real(rkind), dimension(nx,ny,nz), intent(out)  :: SS
        real(rkind), dimension(nx,ny,nz,6) :: fingersq,devstress
        real(rkind), dimension(nx,ny,nz) :: trG, trG2, detG

        call elastic%get_finger(g,finger,fingersq,trG,trG2,detG)
        call elastic%get_devstress(finger,fingersq,trG,trG2,detG,devstress)

        SS = devstress(:,:,:,1)**2 + devstress(:,:,:,4)**2 + devstress(:,:,:,6)**2 + two*(devstress(:,:,:,2)**2 + devstress(:,:,:,3)**2 + devstress(:,:,:,5)**2)

    end subroutine

end module

program test_plastic_deformation
    use kind_parameters, only: rkind
    use test_plastic
    use Sep1SolidEOS,    only: sep1solid
    implicit none

    real(rkind), dimension(nx,ny,nz,9) :: g
    real(rkind), dimension(:,:,:), pointer :: g11,g12,g13,g21,g22,g23,g31,g32,g33
    real(rkind), dimension(nx,ny,nz) :: SS

    call elastic%init(mu)

    g11 => g(:,:,:,1); g12 => g(:,:,:,2); g13 => g(:,:,:,3)
    g21 => g(:,:,:,4); g22 => g(:,:,:,5); g23 => g(:,:,:,6)
    g31 => g(:,:,:,7); g32 => g(:,:,:,8); g33 => g(:,:,:,9)

    g11 = 1.1_rkind; g12 = 0.2_rkind; g13 = 0.0_rkind
    g21 = 0.1_rkind; g22 = 1.1_rkind; g23 = 0.0_rkind
    g31 = 0.0_rkind; g32 = 0.0_rkind; g33 = 1.1_rkind

    call get_SS(g,SS)
    yield = sqrt( (three/four)*SS(1,1,1) ) ! set yield^2 to (3/4)*SS to guarantee SS > (2/3)*yield^2

    write(*,*) "Initial g:"
    write(*,*) g11(1,1,1), g12(1,1,1), g13(1,1,1)
    write(*,*) g21(1,1,1), g22(1,1,1), g23(1,1,1)
    write(*,*) g31(1,1,1), g32(1,1,1), g33(1,1,1)
    write(*,*) "Initial SS - (2/3)yield^2", SS(1,1,1) - (two/three)*yield**2

    call plastic_deformation(g11,g12,g13,g21,g22,g23,g31,g32,g33)

    call get_SS(g,SS)
    write(*,*) "Final g:"
    write(*,*) g11(1,1,1), g12(1,1,1), g13(1,1,1)
    write(*,*) g21(1,1,1), g22(1,1,1), g23(1,1,1)
    write(*,*) g31(1,1,1), g32(1,1,1), g33(1,1,1)
    write(*,*) "Final SS - (2/3)yield^2", SS(1,1,1) - (two/three)*yield**2

end program
