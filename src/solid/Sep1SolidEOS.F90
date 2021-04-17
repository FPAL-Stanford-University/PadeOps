module Sep1SolidEOS

    use kind_parameters, only: rkind
    use constants,       only: zero,half,one,two,three,fourth,sixth,third,twothird,fourthird
    use ElasticEOSMod,   only: elasticeos

    implicit none

    type, extends(elasticeos) :: sep1solid

        real(rkind) :: rho0 = one                    ! Reference density
        real(rkind) :: tau0 = 1.0d-10                ! Plastic relaxation time scale
        real(rkind) :: melt_t = 3.0d0                ! Melt temperature
        real(rkind) :: melt_c = 0.2d0                ! Melt model coefficient
        real(rkind) :: kos_b = 0.0d0                 ! LANL Kospall model coefficient
        real(rkind) :: kos_t = 0.0d0                 ! LANL Kospall model coefficient
        real(rkind) :: kos_h = 0.0d0                 ! LANL Kospall model coefficient
        real(rkind) :: kos_g = 0.0d0                 ! LANL Kospall model coefficient
        real(rkind) :: kos_m = 0.0d0                 ! LANL Kospall model coefficient
        real(rkind) :: kos_q = 0.0d0                 ! LANL Kospall model coefficient
        real(rkind) :: kos_f = 0.0d0                 ! LANL Kospall model coefficient
        real(rkind) :: kos_alpha = 1.0d0             ! LANL Kospall model coefficient
        real(rkind) :: kos_beta = 0.0d0              ! LANL Kospall model coefficient
        real(rkind) :: kos_e = 0.0d0                 ! LANL Kospall model coefficient
        real(rkind) :: mu0 = zero                    ! Shear Modulus initial
        real(rkind) :: yield0 = one                  ! Yield Stress initial
        real(rkind) :: eta_det_ge = 1.0d0               ! Density preserving factor in g eqns.
        real(rkind) :: eta_det_gp = 1.0d0               ! Density preserving factor in g eqns.
        real(rkind) :: eta_det_gt = 1.0d0               ! Density preserving factor in g eqns.
        real(rkind) :: diff_c_ge = 1.0d0               ! Curl preserving factor in g eqns.
        real(rkind) :: diff_c_gp = 1.0d0               ! Curl preserving factor in g eqns.
        real(rkind) :: diff_c_gt = 1.0d0               ! Curl preserving factor in g eqns.


        real(rkind), dimension(:,:,:), allocatable :: mu, yield ! Shear Modulus and Yield Stress

        real(rkind), dimension(:), allocatable :: svdwork     ! Work array for SVD stuff

        integer :: nxd,nyd,nzd,kos_sh


    contains

        ! procedure :: init
        procedure :: get_finger
        procedure :: get_devstress
        procedure :: get_eelastic
        procedure :: get_sos
        procedure :: get_sos2
        procedure :: get_sos2_mixture
        procedure :: plastic_deformation
        procedure :: make_tensor_SPD
        final     :: destroy

    end type

    interface sep1solid
        module procedure init
    end interface

contains

    function init(rho0_,mu_,yield_,tau0_,eta_det_ge_,eta_det_gp_,eta_det_gt_,diff_c_ge_,diff_c_gp_,diff_c_gt_,melt_t_,melt_c_,kos_b_,kos_t_,kos_h_,kos_g_,kos_m_,kos_q_,kos_f_,kos_alpha_,kos_beta_,kos_e_,kos_sh_,nxd_,nyd_,nzd_) result(this)
        type(sep1solid) :: this
        real(rkind), intent(in) :: rho0_, mu_, yield_, tau0_, eta_det_ge_, eta_det_gp_, eta_det_gt_, diff_c_ge_, diff_c_gp_, diff_c_gt_, melt_t_, melt_c_, kos_b_, kos_t_, kos_h_, kos_g_, kos_m_, kos_q_, kos_f_, kos_alpha_, kos_beta_, kos_e_
        integer, intent(in) :: nxd_,nyd_,nzd_,kos_sh_

        integer :: info, lwork
        real(rkind), dimension(3,3) :: g, u, vt
        real(rkind), dimension(3)   :: sval


        this%nxd = nxd_
        this%nyd = nyd_
        this%nzd = nzd_
        this%rho0 = rho0_
        this%mu0 = mu_
        this%yield0 = yield_
        this%tau0 = tau0_
        this%eta_det_ge = eta_det_ge_
        this%eta_det_gp = eta_det_gp_
        this%eta_det_gt = eta_det_gt_
        this%diff_c_ge = diff_c_ge_
        this%diff_c_gp = diff_c_gp_
        this%diff_c_gt = diff_c_gt_
        this%melt_t = melt_t_
        this%melt_c = melt_c_
        this%kos_b = kos_b_
        this%kos_t = kos_t_
        this%kos_h = kos_h_
        this%kos_g = kos_g_
        this%kos_m = kos_m_
        this%kos_q = kos_q_
        this%kos_f = kos_f_
        this%kos_alpha = kos_alpha_
        this%kos_beta = kos_beta_
        this%kos_e = kos_e_
        this%kos_sh = kos_sh_


        if (allocated(this%mu)) deallocate(this%mu); allocate(this%mu(this%nxd,this%nyd,this%nzd))
        if (allocated(this%yield)) deallocate(this%yield); allocate(this%yield(this%nxd,this%nyd,this%nzd))
        this%mu = this%mu0
        this%yield = this%yield0

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

    subroutine get_finger(this,g,finger,fingersq,trG,trG2,detG,use_gTg)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:,:), intent(in)  :: g
        real(rkind), dimension(:,:,:,:), intent(out) :: finger
        real(rkind), dimension(:,:,:,:), intent(out) :: fingersq
        real(rkind), dimension(:,:,:),   intent(out) :: trG, trG2, detG
        logical,                         intent(in), optional :: use_gTg

        associate ( g11 => g(:,:,:,1), g12 => g(:,:,:,2), g13 => g(:,:,:,3), &
                    g21 => g(:,:,:,4), g22 => g(:,:,:,5), g23 => g(:,:,:,6), &
                    g31 => g(:,:,:,7), g32 => g(:,:,:,8), g33 => g(:,:,:,9)  )

            if ( present(use_gTg) .and. use_gTg ) then
                finger(:,:,:,1) = g11
                finger(:,:,:,2) = g12
                finger(:,:,:,3) = g13
                finger(:,:,:,4) = g22
                finger(:,:,:,5) = g23
                finger(:,:,:,6) = g33
            else
                finger(:,:,:,1) = g11*g11 + g21*g21 + g31*g31
                finger(:,:,:,2) = g11*g12 + g21*g22 + g31*g32
                finger(:,:,:,3) = g11*g13 + g21*g23 + g31*g33
                finger(:,:,:,4) = g12*g12 + g22*g22 + g32*g32
                finger(:,:,:,5) = g12*g13 + g22*g23 + g32*g33
                finger(:,:,:,6) = g13*g13 + g23*g23 + g33*g33
            end if

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

    pure subroutine get_devstress(this,finger,fingersq,trG,trG2,detG,devstress,rho0mix,mumix)
        use exits, only: GracefulExit
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:,:), intent(in)  :: finger
        real(rkind), dimension(:,:,:,:), intent(in)  :: fingersq
        real(rkind), dimension(:,:,:),   intent(in)  :: trG, trG2, detG
        real(rkind), dimension(:,:,:,:), intent(out) :: devstress
        real(rkind), dimension(:,:,:), intent(in), optional :: mumix, rho0mix

        real(rkind), dimension(size(finger,1),size(finger,2),size(finger,3)) :: devstmp
        integer :: i

        ! if(.not.present(fingersq)) call GracefulExit("fingersq required for devstress",1111)
        
        do i = 1,6
            devstress(:,:,:,i) = -(detG**(-sixth)*fingersq(:,:,:,i) - detG**sixth*finger(:,:,:,i))
        end do
        devstmp = third*(detG**(-sixth)*trG2 - detG**sixth*trG)
        devstress(:,:,:,1) = devstress(:,:,:,1) + devstmp
        devstress(:,:,:,4) = devstress(:,:,:,4) + devstmp
        devstress(:,:,:,6) = devstress(:,:,:,6) + devstmp

        if(present(rho0mix) .and. present(mumix)) then
            do i = 1, 6
              devstress(:,:,:,i) = devstress(:,:,:,i)*mumix
            enddo
        else
           !devstress = devstress*this%mu
           do i = 1, 6
              devstress(:,:,:,i) = devstress(:,:,:,i)*this%mu
           enddo
        endif

    end subroutine

    subroutine get_eelastic(this,trG,trG2,detG,eelastic,rho0mix,mumix)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: trG,trG2,detG
        real(rkind), dimension(:,:,:), intent(out) :: eelastic
        real(rkind), dimension(:,:,:), intent(in), optional :: mumix, rho0mix

        if(present(rho0mix) .and. present(mumix)) then
            eelastic = fourth*mumix/rho0mix*(detG**(-twothird)*trG2 - two*detG**(-third)*trG + three)
        else
            eelastic = fourth*this%mu/this%rho0*(detG**(-twothird)*trG2 - two*detG**(-third)*trG + three)
        endif 


    end subroutine

    pure subroutine get_sos(this,rhom,sos)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in) :: rhom
        real(rkind), dimension(:,:,:), intent(inout) :: sos

        sos = sqrt(sos**two + fourthird*this%mu/rhom)

    end subroutine

    pure subroutine get_sos2_mixture(this,rhomix,mumix,sos2)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in) :: rhomix,mumix
        real(rkind), dimension(:,:,:), intent(inout) :: sos2

        sos2 = sos2 + fourthird*mumix/rhomix

    end subroutine

    pure subroutine get_sos2(this,rhom,sos2)
        class(sep1solid), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in) :: rhom
        real(rkind), dimension(:,:,:), intent(inout) :: sos2

        sos2 = sos2 + fourthird*this%mu/rhom

    end subroutine

    subroutine plastic_deformation(this, gfull, gpfull, pe, rho, Temp, sxx, sxy, sxz, syy, syz, szz,use_gTg,useOneG,strainHard,cnsrv_g,cnsrv_gt,cnsrv_gp,cnsrv_pe, mumix, yieldmix, rho0mix)
        use kind_parameters, only: clen
        use constants,       only: eps, twothird
        use decomp_2d,       only: nrank
        use exits,           only: GracefulExit, nancheck, message
        class(sep1solid), target, intent(inout) :: this
        real(rkind), dimension(:,:,:,:), intent(inout) :: gfull,gpfull
        real(rkind), dimension(:,:,:), intent(inout) :: pe
        real(rkind), dimension(:,:,:), intent(in) :: sxx,sxy,sxz,syy,syz,szz

        real(rkind), dimension(:,:,:), intent(in), optional :: mumix, yieldmix, rho0mix
        real(rkind), dimension(:,:,:), intent(in) :: rho,Temp

        real(rkind), dimension(3,3) :: g, u, vt, gradf, gradf_new,sigma,gpt1,gpt2,geinv,gedlt,gp,gpt3,gpt4!,gpdlt!,gpinv!,alms,
        real(rkind), dimension(3)   :: sval, beta, Sa, f, f1, f2, dbeta, beta_new, dbeta_new
        real(rkind) :: sqrt_om, betasum, Sabymu_sq, ycrit, C0, t
        real(rkind) :: tol = real(1.D-12,rkind), residual, residual_new
        real(rkind) :: asq,xsq,tot,adiv, detg,detgp,almsn,almsn0,almsnd,lpn,lpn2,sigman
        !real(rkind) :: dtsub,qtsub
        integer :: i,j,k!,n
        integer :: iters
        integer, parameter :: niters = 500
        integer :: lwork, info
        integer, dimension(3) :: ipiv
        integer :: nxp, nyp, nzp
        character(len=clen) :: charout
        real(rkind), dimension(size(gfull,1),size(gfull,2),size(gfull,3)) :: mulocal, yieldlocal,rho0local

        logical, intent(in) :: use_gTg,useOneG,strainHard,cnsrv_g,cnsrv_gt,cnsrv_gp,cnsrv_pe

        nxp = size(gfull,1); nyp = size(gfull,2); nzp = size(gfull,3);

        if(present(yieldmix) .and. present(mumix) .and. present(rho0mix)) then
            mulocal = mumix
            yieldlocal = yieldmix
            rho0local = rho0mix
        else
            mulocal = this%mu
            yieldlocal = this%yield
            rho0local = this%rho0
        endif

        if (use_gTg.and.(.not.strainHard)) then
            ! Get optimal lwork
            lwork = -1
            call dsyev('V', 'U', 3, G, 3, sval, this%svdwork, lwork, info)
            lwork = this%svdwork(1)
        else
            ! Get optimal lwork
            lwork = -1
            call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, this%svdwork, lwork, info)
            lwork = this%svdwork(1)
        end if

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

                    if (use_gTg.and.(.not.strainHard)) then
                        ! Get eigenvalues and eigenvectors of G
                        call dsyev('V', 'U', 3, G, 3, sval, this%svdwork, lwork, info)
                        if(info .ne. 0) print '(A,I6,A)', 'proc ', nrank, ': Problem with DSYEV. Please check.'
                        sval = sqrt(sval)  ! Get singular values of g
                    else
                        ! Get SVD of g
                        call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, this%svdwork, lwork, info)
                        if(info .ne. 0) then
                            write(charout, '(A,I0,A)') 'proc ', nrank, ': Problem with SVD. Please check.'
                            call GracefulExit(charout,3475)
                        end if
                    end if

                    sqrt_om = sval(1)*sval(2)*sval(3)
                    beta = sval**two / sqrt_om**(two/three)

                    betasum = sum( beta*(beta-one) ) / three

                    ! mca: old
                    ! Sa = -mulocal(i,j,k)*sqrt_om * ( beta*(beta-one) - betasum )
                    ! Sabymu_sq = sum(Sa**two) / mulocal(i,j,k)**two
                    ! ycrit = Sabymu_sq - (two/three)*(yieldlocal(i,j,k)/mulocal(i,j,k))**two
                    ! C0 = Sabymu_sq / ycrit
                    ! Sa = Sa*( sqrt(C0 - one)/sqrt(C0) )
                    ! f = Sa / (mulocal(i,j,k)*sqrt_om); f(3) = beta(1)*beta(2)*beta(3) 

                    Sabymu_sq = sqrt( sum( (beta*(beta-one) - betasum)**two ) )

                    !if (ycrit .LE. zero) then !mca: old
                    if (mulocal(i,j,k)*Sabymu_sq*sqrt_om .LE. sqrt(two/three)*yieldlocal(i,j,k)) then
                       ! print '(A)', 'Inconsistency in plastic algorithm, ycrit < 0!'
                       cycle
                    end if

                    !if( (yieldlocal(i,j,k).LE.eps) .OR. (mulocal(i,j,k).LE.eps) ) then !with condition: 4.48313318030560; without condition:8.00000000 why less?? 
                    !   C0 = 1.0d0
                    !else
                       C0 = yieldlocal(i,j,k)/mulocal(i,j,k)
                       !if(use_gTg) C0 = 2.0*C0
                    !endif
                    f = -C0*sqrt(two/three)/(sqrt_om*Sabymu_sq)*( beta*(beta-one) - betasum ); f(3) = beta(1)*beta(2)*beta(3)     ! New function value (target to attain)

                    ! mca: old
                    ! ! Now get new beta
                    ! f = Sa / (mulocal(i,j,k)*sqrt_om); f(3) = beta(1)*beta(2)*beta(3)     ! New function value (target to attain)
                    ! betasum = sum( beta*(beta-one) ) / three

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
                    t = 1._rkind
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
                    !if ((iters >= niters) .OR. ((t <= eps).and.(iters.ge.50))) then !mca
                        write(charout,'(4(A,I0))') 'Newton solve in plastic_deformation did not converge at index ',i,',',j,',',k,' of process ',nrank
                        print '(A)', charout
                        print '(A)', 'g = '
                        print '(4X,3(ES15.5))', gfull(i,j,k,1), gfull(i,j,k,2), gfull(i,j,k,3)
                        print '(4X,3(ES15.5))', gfull(i,j,k,4), gfull(i,j,k,5), gfull(i,j,k,6)
                        print '(4X,3(ES15.5))', gfull(i,j,k,7), gfull(i,j,k,8), gfull(i,j,k,9)
                        !print '(A,ES15.5)', '( ||S||^2 - (2/3) sigma_Y^2 )/mu^2 = ', ycrit

                        print '(A,ES15.5)', 'Relaxation, t = ', t
                        print '(A,ES15.5)', 'Residual = ', residual
                        print '(A,2(I6))', 'iters, niters = ', iters,niters
                        call GracefulExit(charout,6382)
                    end if

                    ! Then get new svals
                    sval = sqrt(beta) * sqrt_om**(one/three)

                    if (use_gTg.and.(.not.strainHard)) then
                        sval = sval*sval ! New eigenvalues of G
                        
                        ! Get g = v*sval*vt
                        u = G; vt = transpose(u)
                        vt(1,:) = vt(1,:)*sval(1); vt(2,:) = vt(2,:)*sval(2); vt(3,:) = vt(3,:)*sval(3)  ! eigval*vt
                        G = MATMUL(u,vt) ! v*eigval*vt
                    else
                        ! Get g = u*sval*vt
                        vt(1,:) = vt(1,:)*sval(1); vt(2,:) = vt(2,:)*sval(2); vt(3,:) = vt(3,:)*sval(3)  ! sval*vt
                        g = MATMUL(u,vt) ! u*sval*vt
                    end if
                    
                    
                    
                    if(strainHard) then
                       !update plastic g: delta g^p_{im} = - g^p_{ik}*(delta g^e_{kj})*(g^e_{jm})^(-1) 
                       ! !using pre-updated values
                       ! detg = -gfull(i,j,k,3)*gfull(i,j,k,5)*gfull(i,j,k,7) + gfull(i,j,k,2)*gfull(i,j,k,6)*gfull(i,j,k,7) + gfull(i,j,k,3)*gfull(i,j,k,4)*gfull(i,j,k,8) - gfull(i,j,k,1)*gfull(i,j,k,6)*gfull(i,j,k,8) - gfull(i,j,k,2)*gfull(i,j,k,4)*gfull(i,j,k,9) + gfull(i,j,k,1)*gfull(i,j,k,5)*gfull(i,j,k,9)

                       ! geinv(1,1) = (gfull(i,j,k,5)*gfull(i,j,k,9) - gfull(i,j,k,6)*gfull(i,j,k,8))/detg
                       ! geinv(1,2) = (gfull(i,j,k,3)*gfull(i,j,k,8) - gfull(i,j,k,2)*gfull(i,j,k,9))/detg
                       ! geinv(1,3) = (gfull(i,j,k,2)*gfull(i,j,k,6) - gfull(i,j,k,3)*gfull(i,j,k,5))/detg
                       ! geinv(2,1) = (gfull(i,j,k,6)*gfull(i,j,k,7) - gfull(i,j,k,4)*gfull(i,j,k,9))/detg
                       ! geinv(2,2) = (gfull(i,j,k,1)*gfull(i,j,k,9) - gfull(i,j,k,3)*gfull(i,j,k,7))/detg
                       ! geinv(2,3) = (gfull(i,j,k,3)*gfull(i,j,k,4) - gfull(i,j,k,1)*gfull(i,j,k,6))/detg
                       ! geinv(3,1) = (gfull(i,j,k,4)*gfull(i,j,k,8) - gfull(i,j,k,5)*gfull(i,j,k,7))/detg
                       ! geinv(3,2) = (gfull(i,j,k,2)*gfull(i,j,k,7) - gfull(i,j,k,1)*gfull(i,j,k,8))/detg
                       ! geinv(3,3) = (gfull(i,j,k,1)*gfull(i,j,k,5) - gfull(i,j,k,2)*gfull(i,j,k,4))/detg


                       !using post-updated values
                       detg = -g(1,3)*g(2,2)*g(3,1) + g(1,2)*g(2,3)*g(3,1) + g(1,3)*g(2,1)*g(3,2) - g(1,1)*g(2,3)*g(3,2) - g(1,2)*g(2,1)*g(3,3) + g(1,1)*g(2,2)*g(3,3)

                       geinv(1,1) = (g(2,2)*g(3,3) - g(2,3)*g(3,2))/detg
                       geinv(1,2) = (g(1,3)*g(3,2) - g(1,2)*g(3,3))/detg
                       geinv(1,3) = (g(1,2)*g(2,3) - g(1,3)*g(2,2))/detg
                       geinv(2,1) = (g(2,3)*g(3,1) - g(2,1)*g(3,3))/detg
                       geinv(2,2) = (g(1,1)*g(3,3) - g(1,3)*g(3,1))/detg
                       geinv(2,3) = (g(1,3)*g(2,1) - g(1,1)*g(2,3))/detg
                       geinv(3,1) = (g(2,1)*g(3,2) - g(2,2)*g(3,1))/detg
                       geinv(3,2) = (g(1,2)*g(3,1) - g(1,1)*g(3,2))/detg
                       geinv(3,3) = (g(1,1)*g(2,2) - g(1,2)*g(2,1))/detg


                       gedlt(1,1) = g(1,1) - gfull(i,j,k,1)
                       gedlt(1,2) = g(1,2) - gfull(i,j,k,2)
                       gedlt(1,3) = g(1,3) - gfull(i,j,k,3)
                       gedlt(2,1) = g(2,1) - gfull(i,j,k,4)
                       gedlt(2,2) = g(2,2) - gfull(i,j,k,5)
                       gedlt(2,3) = g(2,3) - gfull(i,j,k,6)
                       gedlt(3,1) = g(3,1) - gfull(i,j,k,7)
                       gedlt(3,2) = g(3,2) - gfull(i,j,k,8)
                       gedlt(3,3) = g(3,3) - gfull(i,j,k,9)


                       if (use_gTg) then
                          ! !old: close but needs improvment
                          ! !(j,m) or (i,m)
                          ! gpt1(1,1) = gpfull(i,j,k,1)*gedlt(1,1) + gpfull(i,j,k,4)*gedlt(2,1) + gpfull(i,j,k,7)*gedlt(3,1) 
                          ! gpt1(1,2) = gpfull(i,j,k,1)*gedlt(1,2) + gpfull(i,j,k,4)*gedlt(2,2) + gpfull(i,j,k,7)*gedlt(3,2) 
                          ! gpt1(1,3) = gpfull(i,j,k,1)*gedlt(1,3) + gpfull(i,j,k,4)*gedlt(2,3) + gpfull(i,j,k,7)*gedlt(3,3) 

                          ! gpt1(2,1) = gpfull(i,j,k,2)*gedlt(1,1) + gpfull(i,j,k,5)*gedlt(2,1) + gpfull(i,j,k,8)*gedlt(3,1) 
                          ! gpt1(2,2) = gpfull(i,j,k,2)*gedlt(1,2) + gpfull(i,j,k,5)*gedlt(2,2) + gpfull(i,j,k,8)*gedlt(3,2) 
                          ! gpt1(2,3) = gpfull(i,j,k,2)*gedlt(1,3) + gpfull(i,j,k,5)*gedlt(2,3) + gpfull(i,j,k,8)*gedlt(3,3) 

                          ! gpt1(3,1) = gpfull(i,j,k,3)*gedlt(1,1) + gpfull(i,j,k,6)*gedlt(2,1) + gpfull(i,j,k,9)*gedlt(3,1) 
                          ! gpt1(3,2) = gpfull(i,j,k,3)*gedlt(1,2) + gpfull(i,j,k,6)*gedlt(2,2) + gpfull(i,j,k,9)*gedlt(3,2) 
                          ! gpt1(3,3) = gpfull(i,j,k,3)*gedlt(1,3) + gpfull(i,j,k,6)*gedlt(2,3) + gpfull(i,j,k,9)*gedlt(3,3) 


                          ! gpt2(1,1) =    2.0*(gpt1(1,1)*geinv(1,1) + gpt1(1,2)*geinv(2,1) + gpt1(1,3)*geinv(3,1))

                          ! gpt2(1,2) =    (gpt1(1,1)*geinv(1,2) + gpt1(1,2)*geinv(2,2) + gpt1(1,3)*geinv(3,2)) &
                          !              + (gpt1(2,1)*geinv(1,1) + gpt1(2,2)*geinv(2,1) + gpt1(2,3)*geinv(3,1))

                          ! gpt2(1,3) =    (gpt1(1,1)*geinv(1,3) + gpt1(1,2)*geinv(2,3) + gpt1(1,3)*geinv(3,3)) &
                          !              + (gpt1(3,1)*geinv(1,1) + gpt1(3,2)*geinv(2,1) + gpt1(3,3)*geinv(3,1))

                          ! gpt2(2,1) =    gpt2(1,2)

                          ! gpt2(2,2) =    2.0*(gpt1(2,1)*geinv(1,2) + gpt1(2,2)*geinv(2,2) + gpt1(2,3)*geinv(3,2))

                          ! gpt2(2,3) =    (gpt1(2,1)*geinv(1,3) + gpt1(2,2)*geinv(2,3) + gpt1(2,3)*geinv(3,3)) &
                          !              + (gpt1(3,1)*geinv(1,2) + gpt1(3,2)*geinv(2,2) + gpt1(3,3)*geinv(3,2))

                          ! gpt2(3,1) =    gpt2(1,3)                               

                          ! gpt2(3,2) =    gpt2(2,3)                               

                          ! gpt2(3,3) =    2.0*(gpt1(3,1)*geinv(1,3) + gpt1(3,2)*geinv(2,3) + gpt1(3,3)*geinv(3,3))



                          !new
                          gpt1(1,1) = gpfull(i,j,k,1)*gfull(i,j,k,1) + gpfull(i,j,k,4)*gfull(i,j,k,4) + gpfull(i,j,k,7)*gfull(i,j,k,7) 
                          gpt1(1,2) = gpfull(i,j,k,1)*gfull(i,j,k,2) + gpfull(i,j,k,4)*gfull(i,j,k,5) + gpfull(i,j,k,7)*gfull(i,j,k,8) 
                          gpt1(1,3) = gpfull(i,j,k,1)*gfull(i,j,k,3) + gpfull(i,j,k,4)*gfull(i,j,k,6) + gpfull(i,j,k,7)*gfull(i,j,k,9) 

                          gpt1(2,1) = gpfull(i,j,k,2)*gfull(i,j,k,1) + gpfull(i,j,k,5)*gfull(i,j,k,4) + gpfull(i,j,k,8)*gfull(i,j,k,7) 
                          gpt1(2,2) = gpfull(i,j,k,2)*gfull(i,j,k,2) + gpfull(i,j,k,5)*gfull(i,j,k,5) + gpfull(i,j,k,8)*gfull(i,j,k,8) 
                          gpt1(2,3) = gpfull(i,j,k,2)*gfull(i,j,k,3) + gpfull(i,j,k,5)*gfull(i,j,k,6) + gpfull(i,j,k,8)*gfull(i,j,k,9) 

                          gpt1(3,1) = gpfull(i,j,k,3)*gfull(i,j,k,1) + gpfull(i,j,k,6)*gfull(i,j,k,4) + gpfull(i,j,k,9)*gfull(i,j,k,7) 
                          gpt1(3,2) = gpfull(i,j,k,3)*gfull(i,j,k,2) + gpfull(i,j,k,6)*gfull(i,j,k,5) + gpfull(i,j,k,9)*gfull(i,j,k,8) 
                          gpt1(3,3) = gpfull(i,j,k,3)*gfull(i,j,k,3) + gpfull(i,j,k,6)*gfull(i,j,k,6) + gpfull(i,j,k,9)*gfull(i,j,k,9) 


                          gpt2(1,1) = gpt1(1,1)*gfull(i,j,k,1) + gpt1(2,1)*gfull(i,j,k,4) + gpt1(3,1)*gfull(i,j,k,7) 
                          gpt2(1,2) = gpt1(1,1)*gfull(i,j,k,2) + gpt1(2,1)*gfull(i,j,k,5) + gpt1(3,1)*gfull(i,j,k,8) 
                          gpt2(1,3) = gpt1(1,1)*gfull(i,j,k,3) + gpt1(2,1)*gfull(i,j,k,6) + gpt1(3,1)*gfull(i,j,k,9)
 
                          gpt2(2,1) = gpt1(1,2)*gfull(i,j,k,1) + gpt1(2,2)*gfull(i,j,k,4) + gpt1(3,2)*gfull(i,j,k,7) 
                          gpt2(2,2) = gpt1(1,2)*gfull(i,j,k,2) + gpt1(2,2)*gfull(i,j,k,5) + gpt1(3,2)*gfull(i,j,k,8) 
                          gpt2(2,3) = gpt1(1,2)*gfull(i,j,k,3) + gpt1(2,2)*gfull(i,j,k,6) + gpt1(3,2)*gfull(i,j,k,9) 

                          gpt2(3,1) = gpt1(1,3)*gfull(i,j,k,1) + gpt1(2,3)*gfull(i,j,k,4) + gpt1(3,3)*gfull(i,j,k,7) 
                          gpt2(3,2) = gpt1(1,3)*gfull(i,j,k,2) + gpt1(2,3)*gfull(i,j,k,5) + gpt1(3,3)*gfull(i,j,k,8) 
                          gpt2(3,3) = gpt1(1,3)*gfull(i,j,k,3) + gpt1(2,3)*gfull(i,j,k,6) + gpt1(3,3)*gfull(i,j,k,9) 


                          gpt3(1,1) = geinv(1,1)*gpt2(1,1) + geinv(2,1)*gpt2(2,1) + geinv(3,1)*gpt2(3,1) 
                          gpt3(1,2) = geinv(1,1)*gpt2(1,2) + geinv(2,1)*gpt2(2,2) + geinv(3,1)*gpt2(3,2) 
                          gpt3(1,3) = geinv(1,1)*gpt2(1,3) + geinv(2,1)*gpt2(2,3) + geinv(3,1)*gpt2(3,3)
 
                          gpt3(2,1) = geinv(1,2)*gpt2(1,1) + geinv(2,2)*gpt2(2,1) + geinv(3,2)*gpt2(3,1) 
                          gpt3(2,2) = geinv(1,2)*gpt2(1,2) + geinv(2,2)*gpt2(2,2) + geinv(3,2)*gpt2(3,2) 
                          gpt3(2,3) = geinv(1,2)*gpt2(1,3) + geinv(2,2)*gpt2(2,3) + geinv(3,2)*gpt2(3,3) 

                          gpt3(3,1) = geinv(1,3)*gpt2(1,1) + geinv(2,3)*gpt2(2,1) + geinv(3,3)*gpt2(3,1) 
                          gpt3(3,2) = geinv(1,3)*gpt2(1,2) + geinv(2,3)*gpt2(2,2) + geinv(3,3)*gpt2(3,2) 
                          gpt3(3,3) = geinv(1,3)*gpt2(1,3) + geinv(2,3)*gpt2(2,3) + geinv(3,3)*gpt2(3,3) 


                          gpt4(1,1) = gpt3(1,1)*geinv(1,1) + gpt3(1,2)*geinv(2,1) + gpt3(1,3)*geinv(3,1) 
                          gpt4(1,2) = gpt3(1,1)*geinv(1,2) + gpt3(1,2)*geinv(2,2) + gpt3(1,3)*geinv(3,2) 
                          gpt4(1,3) = gpt3(1,1)*geinv(1,3) + gpt3(1,2)*geinv(2,3) + gpt3(1,3)*geinv(3,3)
 
                          gpt4(2,1) = gpt3(2,1)*geinv(1,1) + gpt3(2,2)*geinv(2,1) + gpt3(2,3)*geinv(3,1) 
                          gpt4(2,2) = gpt3(2,1)*geinv(1,2) + gpt3(2,2)*geinv(2,2) + gpt3(2,3)*geinv(3,2) 
                          gpt4(2,3) = gpt3(2,1)*geinv(1,3) + gpt3(2,2)*geinv(2,3) + gpt3(2,3)*geinv(3,3) 

                          gpt4(3,1) = gpt3(3,1)*geinv(1,1) + gpt3(3,2)*geinv(2,1) + gpt3(3,3)*geinv(3,1) 
                          gpt4(3,2) = gpt3(3,1)*geinv(1,2) + gpt3(3,2)*geinv(2,2) + gpt3(3,3)*geinv(3,2) 
                          gpt4(3,3) = gpt3(3,1)*geinv(1,3) + gpt3(3,2)*geinv(2,3) + gpt3(3,3)*geinv(3,3) 


                          ! gpt4(1,1) = gpt3(1,1)*geinv(1,1) + gpt3(2,1)*geinv(2,1) + gpt3(3,1)*geinv(3,1) 
                          ! gpt4(1,2) = gpt3(1,1)*geinv(1,2) + gpt3(2,1)*geinv(2,2) + gpt3(3,1)*geinv(3,2) 
                          ! gpt4(1,3) = gpt3(1,1)*geinv(1,3) + gpt3(2,1)*geinv(2,3) + gpt3(3,1)*geinv(3,3)
 
                          ! gpt4(2,1) = gpt3(1,2)*geinv(1,1) + gpt3(2,2)*geinv(2,1) + gpt3(3,2)*geinv(3,1) 
                          ! gpt4(2,2) = gpt3(1,2)*geinv(1,2) + gpt3(2,2)*geinv(2,2) + gpt3(3,2)*geinv(3,2) 
                          ! gpt4(2,3) = gpt3(1,2)*geinv(1,3) + gpt3(2,2)*geinv(2,3) + gpt3(3,2)*geinv(3,3) 

                          ! gpt4(3,1) = gpt3(1,3)*geinv(1,1) + gpt3(2,3)*geinv(2,1) + gpt3(3,3)*geinv(3,1) 
                          ! gpt4(3,2) = gpt3(1,3)*geinv(1,2) + gpt3(2,3)*geinv(2,2) + gpt3(3,3)*geinv(3,2) 
                          ! gpt4(3,3) = gpt3(1,3)*geinv(1,3) + gpt3(2,3)*geinv(2,3) + gpt3(3,3)*geinv(3,3) 


                          gpt2(1,1) = gpt4(1,1)
                          gpt2(1,2) = gpt4(1,2)
                          gpt2(1,3) = gpt4(1,3)
                          gpt2(2,1) = gpt4(2,1)
                          gpt2(2,2) = gpt4(2,2)
                          gpt2(2,3) = gpt4(2,3)
                          gpt2(3,1) = gpt4(3,1)
                          gpt2(3,2) = gpt4(3,2)
                          gpt2(3,3) = gpt4(3,3)

                       else
                          ! !old good
                          ! gpt1(1,1) = gpfull(i,j,k,1)*gedlt(1,1) + gpfull(i,j,k,2)*gedlt(2,1) + gpfull(i,j,k,3)*gedlt(3,1) 
                          ! gpt1(1,2) = gpfull(i,j,k,1)*gedlt(1,2) + gpfull(i,j,k,2)*gedlt(2,2) + gpfull(i,j,k,3)*gedlt(3,2) 
                          ! gpt1(1,3) = gpfull(i,j,k,1)*gedlt(1,3) + gpfull(i,j,k,2)*gedlt(2,3) + gpfull(i,j,k,3)*gedlt(3,3) 

                          ! gpt1(2,1) = gpfull(i,j,k,4)*gedlt(1,1) + gpfull(i,j,k,5)*gedlt(2,1) + gpfull(i,j,k,6)*gedlt(3,1) 
                          ! gpt1(2,2) = gpfull(i,j,k,4)*gedlt(1,2) + gpfull(i,j,k,5)*gedlt(2,2) + gpfull(i,j,k,6)*gedlt(3,2) 
                          ! gpt1(2,3) = gpfull(i,j,k,4)*gedlt(1,3) + gpfull(i,j,k,5)*gedlt(2,3) + gpfull(i,j,k,6)*gedlt(3,3) 

                          ! gpt1(3,1) = gpfull(i,j,k,7)*gedlt(1,1) + gpfull(i,j,k,8)*gedlt(2,1) + gpfull(i,j,k,9)*gedlt(3,1) 
                          ! gpt1(3,2) = gpfull(i,j,k,7)*gedlt(1,2) + gpfull(i,j,k,8)*gedlt(2,2) + gpfull(i,j,k,9)*gedlt(3,2) 
                          ! gpt1(3,3) = gpfull(i,j,k,7)*gedlt(1,3) + gpfull(i,j,k,8)*gedlt(2,3) + gpfull(i,j,k,9)*gedlt(3,3) 


                          ! gpt2(1,1) = gpt1(1,1)*geinv(1,1) + gpt1(1,2)*geinv(2,1) + gpt1(1,3)*geinv(3,1)
                          ! gpt2(1,2) = gpt1(1,1)*geinv(1,2) + gpt1(1,2)*geinv(2,2) + gpt1(1,3)*geinv(3,2)
                          ! gpt2(1,3) = gpt1(1,1)*geinv(1,3) + gpt1(1,2)*geinv(2,3) + gpt1(1,3)*geinv(3,3)

                          ! gpt2(2,1) = gpt1(2,1)*geinv(1,1) + gpt1(2,2)*geinv(2,1) + gpt1(2,3)*geinv(3,1)
                          ! gpt2(2,2) = gpt1(2,1)*geinv(1,2) + gpt1(2,2)*geinv(2,2) + gpt1(2,3)*geinv(3,2)
                          ! gpt2(2,3) = gpt1(2,1)*geinv(1,3) + gpt1(2,2)*geinv(2,3) + gpt1(2,3)*geinv(3,3)

                          ! gpt2(3,1) = gpt1(3,1)*geinv(1,1) + gpt1(3,2)*geinv(2,1) + gpt1(3,3)*geinv(3,1)
                          ! gpt2(3,2) = gpt1(3,1)*geinv(1,2) + gpt1(3,2)*geinv(2,2) + gpt1(3,3)*geinv(3,2)
                          ! gpt2(3,3) = gpt1(3,1)*geinv(1,3) + gpt1(3,2)*geinv(2,3) + gpt1(3,3)*geinv(3,3)


                          !new
                          gpt1(1,1) = gpfull(i,j,k,1)*gfull(i,j,k,1) + gpfull(i,j,k,2)*gfull(i,j,k,4) + gpfull(i,j,k,3)*gfull(i,j,k,7) 
                          gpt1(1,2) = gpfull(i,j,k,1)*gfull(i,j,k,2) + gpfull(i,j,k,2)*gfull(i,j,k,5) + gpfull(i,j,k,3)*gfull(i,j,k,8) 
                          gpt1(1,3) = gpfull(i,j,k,1)*gfull(i,j,k,3) + gpfull(i,j,k,2)*gfull(i,j,k,6) + gpfull(i,j,k,3)*gfull(i,j,k,9) 

                          gpt1(2,1) = gpfull(i,j,k,4)*gfull(i,j,k,1) + gpfull(i,j,k,5)*gfull(i,j,k,4) + gpfull(i,j,k,6)*gfull(i,j,k,7) 
                          gpt1(2,2) = gpfull(i,j,k,4)*gfull(i,j,k,2) + gpfull(i,j,k,5)*gfull(i,j,k,5) + gpfull(i,j,k,6)*gfull(i,j,k,8) 
                          gpt1(2,3) = gpfull(i,j,k,4)*gfull(i,j,k,3) + gpfull(i,j,k,5)*gfull(i,j,k,6) + gpfull(i,j,k,6)*gfull(i,j,k,9) 

                          gpt1(3,1) = gpfull(i,j,k,7)*gfull(i,j,k,1) + gpfull(i,j,k,8)*gfull(i,j,k,4) + gpfull(i,j,k,9)*gfull(i,j,k,7) 
                          gpt1(3,2) = gpfull(i,j,k,7)*gfull(i,j,k,2) + gpfull(i,j,k,8)*gfull(i,j,k,5) + gpfull(i,j,k,9)*gfull(i,j,k,8) 
                          gpt1(3,3) = gpfull(i,j,k,7)*gfull(i,j,k,3) + gpfull(i,j,k,8)*gfull(i,j,k,6) + gpfull(i,j,k,9)*gfull(i,j,k,9) 


                          gpt2(1,1) = gpt1(1,1)*geinv(1,1) + gpt1(1,2)*geinv(2,1) + gpt1(1,3)*geinv(3,1)
                          gpt2(1,2) = gpt1(1,1)*geinv(1,2) + gpt1(1,2)*geinv(2,2) + gpt1(1,3)*geinv(3,2)
                          gpt2(1,3) = gpt1(1,1)*geinv(1,3) + gpt1(1,2)*geinv(2,3) + gpt1(1,3)*geinv(3,3)

                          gpt2(2,1) = gpt1(2,1)*geinv(1,1) + gpt1(2,2)*geinv(2,1) + gpt1(2,3)*geinv(3,1)
                          gpt2(2,2) = gpt1(2,1)*geinv(1,2) + gpt1(2,2)*geinv(2,2) + gpt1(2,3)*geinv(3,2)
                          gpt2(2,3) = gpt1(2,1)*geinv(1,3) + gpt1(2,2)*geinv(2,3) + gpt1(2,3)*geinv(3,3)

                          gpt2(3,1) = gpt1(3,1)*geinv(1,1) + gpt1(3,2)*geinv(2,1) + gpt1(3,3)*geinv(3,1)
                          gpt2(3,2) = gpt1(3,1)*geinv(1,2) + gpt1(3,2)*geinv(2,2) + gpt1(3,3)*geinv(3,2)
                          gpt2(3,3) = gpt1(3,1)*geinv(1,3) + gpt1(3,2)*geinv(2,3) + gpt1(3,3)*geinv(3,3)

                       endif



                       !gp -- before implicit update
                       gp(1,1) = gpfull(i,j,k,1)
                       gp(1,2) = gpfull(i,j,k,2)
                       gp(1,3) = gpfull(i,j,k,3)
                       gp(2,1) = gpfull(i,j,k,4)
                       gp(2,2) = gpfull(i,j,k,5)
                       gp(2,3) = gpfull(i,j,k,6)
                       gp(3,1) = gpfull(i,j,k,7)
                       gp(3,2) = gpfull(i,j,k,8)
                       gp(3,3) = gpfull(i,j,k,9)

                       ! !Almansi plastic strain -- before implicit update
                       ! alms(1,1) = 0.5*(1.0 - (gp(1,1)*gp(1,1) + gp(2,1)*gp(2,1) + gp(3,1)*gp(3,1)) )
                       ! alms(1,2) = 0.5*(0.0 - (gp(1,1)*gp(1,2) + gp(2,1)*gp(2,2) + gp(3,1)*gp(3,2)) )
                       ! alms(1,3) = 0.5*(0.0 - (gp(1,1)*gp(1,3) + gp(2,1)*gp(2,3) + gp(3,1)*gp(3,3)) )
                       ! alms(2,1) = 0.5*(0.0 - (gp(1,2)*gp(1,1) + gp(2,2)*gp(2,1) + gp(3,2)*gp(3,1)) )
                       ! alms(2,2) = 0.5*(1.0 - (gp(1,2)*gp(1,2) + gp(2,2)*gp(2,2) + gp(3,2)*gp(3,2)) )
                       ! alms(2,3) = 0.5*(0.0 - (gp(1,2)*gp(1,3) + gp(2,2)*gp(2,3) + gp(3,2)*gp(3,3)) )
                       ! alms(3,1) = 0.5*(0.0 - (gp(1,3)*gp(1,1) + gp(2,3)*gp(2,1) + gp(3,3)*gp(3,1)) )
                       ! alms(3,2) = 0.5*(0.0 - (gp(1,3)*gp(1,2) + gp(2,3)*gp(2,2) + gp(3,3)*gp(3,2)) )
                       ! alms(3,3) = 0.5*(1.0 - (gp(1,3)*gp(1,3) + gp(2,3)*gp(2,3) + gp(3,3)*gp(3,3)) )

                       !Eulerian-Almansi strain tensor: ea = (I-(g_p)^T.g_p)/2 --- not explicitly calculated
                       !strain norm: e_p = sqrt(2/3*ea_ij*ea_ij)
                       if (use_gTg) then
                          almsn0 = sqrt(( (1.0d0 - gp(1,1))**2 + gp(1,2)**2 + gp(1,3)**2 + gp(2,1)**2 + (1.0d0 - gp(2,2))**2 + gp(2,3)**2 + gp(3,1)**2 + gp(3,2)**2 + (1.0d0 - gp(3,3))**2 )/6.0d0)
                       else
                          almsn0 = sqrt( ( (1.0 - gp(1,1)**2 - gp(2,1)**2 - gp(3,1)**2)**2 + (1.0 - gp(1,2)**2 - gp(2,2)**2 - gp(3,2)**2)**2 + (1.0 - gp(1,3)**2 - gp(2,3)**2 - gp(3,3)**2)**2 )/6.0d0 + ( (-gp(1,1)*gp(1,2) - gp(2,1)*gp(2,2) - gp(3,1)*gp(3,2))**2 + (-gp(1,1)*gp(1,3) - gp(2,1)*gp(2,3) - gp(3,1)*gp(3,3))**2 + (-gp(1,2)*gp(1,3) - gp(2,2)*gp(2,3) - gp(3,2)*gp(3,3))**2 )/3.0d0 ) !same as solidmod
                       endif

                       !update plastic g
                       ! !old
                       ! gpfull(i,j,k,1) = gpfull(i,j,k,1) - gpt2(1,1)
                       ! gpfull(i,j,k,2) = gpfull(i,j,k,2) - gpt2(1,2)
                       ! gpfull(i,j,k,3) = gpfull(i,j,k,3) - gpt2(1,3)
                       ! gpfull(i,j,k,4) = gpfull(i,j,k,4) - gpt2(2,1)
                       ! gpfull(i,j,k,5) = gpfull(i,j,k,5) - gpt2(2,2)
                       ! gpfull(i,j,k,6) = gpfull(i,j,k,6) - gpt2(2,3)
                       ! gpfull(i,j,k,7) = gpfull(i,j,k,7) - gpt2(3,1)
                       ! gpfull(i,j,k,8) = gpfull(i,j,k,8) - gpt2(3,2)
                       ! gpfull(i,j,k,9) = gpfull(i,j,k,9) - gpt2(3,3)

                       !new
                       gpfull(i,j,k,1) = gpt2(1,1)
                       gpfull(i,j,k,2) = gpt2(1,2)
                       gpfull(i,j,k,3) = gpt2(1,3)
                       gpfull(i,j,k,4) = gpt2(2,1)
                       gpfull(i,j,k,5) = gpt2(2,2)
                       gpfull(i,j,k,6) = gpt2(2,3)
                       gpfull(i,j,k,7) = gpt2(3,1)
                       gpfull(i,j,k,8) = gpt2(3,2)
                       gpfull(i,j,k,9) = gpt2(3,3)


                       !for plastic entropy
                       ! !delta gp
                       ! gpdlt(1,1) = gpfull(i,j,k,1) - gp(1,1)
                       ! gpdlt(1,2) = gpfull(i,j,k,2) - gp(1,2)
                       ! gpdlt(1,3) = gpfull(i,j,k,3) - gp(1,3)
                       ! gpdlt(2,1) = gpfull(i,j,k,4) - gp(2,1)
                       ! gpdlt(2,2) = gpfull(i,j,k,5) - gp(2,2)
                       ! gpdlt(2,3) = gpfull(i,j,k,6) - gp(2,3)
                       ! gpdlt(3,1) = gpfull(i,j,k,7) - gp(3,1)
                       ! gpdlt(3,2) = gpfull(i,j,k,8) - gp(3,2)
                       ! gpdlt(3,3) = gpfull(i,j,k,9) - gp(3,3)

                       !gp -- after implicit update
                       gp(1,1) = gpfull(i,j,k,1)
                       gp(1,2) = gpfull(i,j,k,2)
                       gp(1,3) = gpfull(i,j,k,3)
                       gp(2,1) = gpfull(i,j,k,4)
                       gp(2,2) = gpfull(i,j,k,5)
                       gp(2,3) = gpfull(i,j,k,6)
                       gp(3,1) = gpfull(i,j,k,7)
                       gp(3,2) = gpfull(i,j,k,8)
                       gp(3,3) = gpfull(i,j,k,9)

                       ! !using post-updated values
                       ! detgp = -gp(1,3)*gp(2,2)*gp(3,1) + gp(1,2)*gp(2,3)*gp(3,1) + gp(1,3)*gp(2,1)*gp(3,2) - gp(1,1)*gp(2,3)*gp(3,2) - gp(1,2)*gp(2,1)*gp(3,3) + gp(1,1)*gp(2,2)*gp(3,3)

                       ! gpinv(1,1) = (gp(2,2)*gp(3,3) - gp(2,3)*gp(3,2))/detgp
                       ! gpinv(1,2) = (gp(1,3)*gp(3,2) - gp(1,2)*gp(3,3))/detgp
                       ! gpinv(1,3) = (gp(1,2)*gp(2,3) - gp(1,3)*gp(2,2))/detgp
                       ! gpinv(2,1) = (gp(2,3)*gp(3,1) - gp(2,1)*gp(3,3))/detgp
                       ! gpinv(2,2) = (gp(1,1)*gp(3,3) - gp(1,3)*gp(3,1))/detgp
                       ! gpinv(2,3) = (gp(1,3)*gp(2,1) - gp(1,1)*gp(2,3))/detgp
                       ! gpinv(3,1) = (gp(2,1)*gp(3,2) - gp(2,2)*gp(3,1))/detgp
                       ! gpinv(3,2) = (gp(1,2)*gp(3,1) - gp(1,1)*gp(3,2))/detgp
                       ! gpinv(3,3) = (gp(1,1)*gp(2,2) - gp(1,2)*gp(2,1))/detgp


                       ! !Almansi plastic strain -- after implicit update
                       ! alms(1,1) = 0.5*(1.0 - (gp(1,1)*gp(1,1) + gp(2,1)*gp(2,1) + gp(3,1)*gp(3,1)) )
                       ! alms(1,2) = 0.5*(0.0 - (gp(1,1)*gp(1,2) + gp(2,1)*gp(2,2) + gp(3,1)*gp(3,2)) )
                       ! alms(1,3) = 0.5*(0.0 - (gp(1,1)*gp(1,3) + gp(2,1)*gp(2,3) + gp(3,1)*gp(3,3)) )
                       ! alms(2,1) = 0.5*(0.0 - (gp(1,2)*gp(1,1) + gp(2,2)*gp(2,1) + gp(3,2)*gp(3,1)) )
                       ! alms(2,2) = 0.5*(1.0 - (gp(1,2)*gp(1,2) + gp(2,2)*gp(2,2) + gp(3,2)*gp(3,2)) )
                       ! alms(2,3) = 0.5*(0.0 - (gp(1,2)*gp(1,3) + gp(2,2)*gp(2,3) + gp(3,2)*gp(3,3)) )
                       ! alms(3,1) = 0.5*(0.0 - (gp(1,3)*gp(1,1) + gp(2,3)*gp(2,1) + gp(3,3)*gp(3,1)) )
                       ! alms(3,2) = 0.5*(0.0 - (gp(1,3)*gp(1,2) + gp(2,3)*gp(2,2) + gp(3,3)*gp(3,2)) )
                       ! alms(3,3) = 0.5*(1.0 - (gp(1,3)*gp(1,3) + gp(2,3)*gp(2,3) + gp(3,3)*gp(3,3)) )

                       !Eulerian-Almansi strain tensor: ea = (I-(g_p)^T.g_p)/2 --- not explicitly calculated
                       !strain norm: e_p = sqrt(2/3*ea_ij*ea_ij)
                       if (use_gTg) then 
                          almsn = sqrt(( (1.0d0 - gp(1,1))**2 + gp(1,2)**2 + gp(1,3)**2 + gp(2,1)**2 + (1.0d0 - gp(2,2))**2 + gp(2,3)**2 + gp(3,1)**2 + gp(3,2)**2 + (1.0d0 - gp(3,3))**2 )/6.0d0)
                       else
                          almsn = sqrt( ( (1.0 - gp(1,1)**2 - gp(2,1)**2 - gp(3,1)**2)**2 + (1.0 - gp(1,2)**2 - gp(2,2)**2 - gp(3,2)**2)**2 + (1.0 - gp(1,3)**2 - gp(2,3)**2 - gp(3,3)**2)**2 )/6.0d0 + ( (-gp(1,1)*gp(1,2) - gp(2,1)*gp(2,2) - gp(3,1)*gp(3,2))**2 + (-gp(1,1)*gp(1,3) - gp(2,1)*gp(2,3) - gp(3,1)*gp(3,3))**2 + (-gp(1,2)*gp(1,3) - gp(2,2)*gp(2,3) - gp(3,2)*gp(3,3))**2 )/3.0d0 ) !same as solidmod
                       endif

                       !update pe -- note this is called after g and pe are converted back to primitive
                       !pe(i,j,k) =  pe(i,j,k) + (almsn-almsn0) !correct dimensinos for strain - compare this with e_p and e_pp
                       pe(i,j,k) =  pe(i,j,k) + (almsn-almsn0)*yieldlocal(i,j,k) !correct dimensions for energy
                       !pe(i,j,k) =  pe(i,j,k) + (almsn-almsn0)*yieldlocal(i,j,k)/Temp(i,j,k) !correct dimensions for entropy

                    endif


                    !update elastic g
                    gfull(i,j,k,1) = g(1,1); gfull(i,j,k,2) = g(1,2); gfull(i,j,k,3) = g(1,3)
                    gfull(i,j,k,4) = g(2,1); gfull(i,j,k,5) = g(2,2); gfull(i,j,k,6) = g(2,3)
                    gfull(i,j,k,7) = g(3,1); gfull(i,j,k,8) = g(3,2); gfull(i,j,k,9) = g(3,3)



                    ! ! old
                    ! !plastic entropy RHS
 
                    ! ! Get eigenvalues and eigenvectors of deviatoric stress
                    ! sigma(1,1) = sxx(i,j,k); sigma(1,2) = sxy(i,j,k); sigma(1,3) = sxz(i,j,k)
                    ! sigma(2,1) = sxy(i,j,k); sigma(2,2) = syy(i,j,k); sigma(2,3) = syz(i,j,k)
                    ! sigma(3,1) = sxz(i,j,k); sigma(3,2) = syz(i,j,k); sigma(3,3) = szz(i,j,k)

                    ! call dsyev('V', 'U', 3, sigma, 3, sval, this%svdwork, lwork, info)
                    ! if(info .ne. 0) then
                    !    print '(A,I6,A)', 'proc ', nrank, ': Problem with DSYEV. Please check.'
                    !    call GracefulExit(charout,3475)
                    ! end if

                    ! asq = (sval(1)**2+sval(2)**2+sval(3)**2)*mulocal(i,j,k)**2
                    ! xsq = 2.0d0/3.0d0*(yieldlocal(i,j,k)/mulocal(i,j,k))**2
                    ! tot = tdel/this%tau0
                    
                    ! !if((asq.lt.1.e-10).or.(xsq.lt.1.e-10)) then
                    ! !   print*,asq,xsq
                    ! !endif

                    ! adiv = (xsq-asq)*exp(-2.0d0*xsq*tot)+asq
                    ! adiv = sign(max(abs(adiv),eps),adiv)

                    ! !!!!!!pe = pe + max( (mulocal(i,j,k)*asq*(asq-xsq)) / (4.0d0*rho(i,j,k)**2/rho0local(i,j,k)*Temp(i,j,k)) * (1.0d0 - exp(-2.0d0*xsq*tot)) / (asq-(asq-xsq)*exp(-2.0d0*xsq*tot)) , 0.0d0)
                    
                    ! !!old
                    ! !pe = pe + max( ( (mulocal(i,j,k)*asq*(asq-xsq)) / (4.0d0*rho(i,j,k)**2/rho0local(i,j,k)*Temp(i,j,k)) ) * ( (1.0d0 - exp(-2.0d0*xsq*tot)) /  adiv ) , 0.0d0)


                    
                    ! !!!tmp  = (rho*this%Ys + this%elastic%rho0*detg *epssmall)/(this%VF + epssmall)   
                    ! !need to check correct rho below
                    ! pe(i,j,k) = pe(i,j,k) + max( ( (mulocal(i,j,k)*asq*(asq-xsq)) / (4.0d0*rho(i,j,k)/rho0local(i,j,k)*Temp(i,j,k)) ) * ( (1.0d0 - exp(-2.0d0*xsq*tot)) /  adiv ) , 0.0d0)
                    
                    

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
