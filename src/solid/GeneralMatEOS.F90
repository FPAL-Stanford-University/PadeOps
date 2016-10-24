module GeneralMatEOS

    use kind_parameters, only: rkind
    use constants,       only: half,one,zero,third,twothird,sixth,three,two,six
    use EOSMod,          only: eos
    use decomp_2d,       only: decomp_info, nrank

    implicit none

    !type, extends(GenSolFluidEos) :: generaleos
    type :: generaleos

        ! EOS (77)-(84) in Hill et al., JCP 229 (2010).
        real(rkind) :: K0   != 138.0d9
        real(rkind) :: Kp   != 4.96_rkind
        real(rkind) :: G0   != 46.9d9
        real(rkind) :: Gp   != 0.57_rkind
        real(rkind) :: beta != 0.0_rkind
        real(rkind) :: T0   != 300.0_rkind
        real(rkind) :: Cv   != 3.9d-4
        real(rkind) :: gam0 != 1.96_rkind
        real(rkind) :: qpar != 1.0_rkind
        

        ! Godunov-Romenskii EOS, (3.27)-(3.29) in Lopez-Ortega PhD Thesis.
        ! Default values: Copper in Lopez-Ortega et al., JCP, 257 (2014).
        real(rkind) :: mu0   != 39.38d9
        !real(rkind) :: beta  != 3.0_rkind
        !real(rkind) :: Kp     != 15.28d6
        real(rkind) :: alp   != 1.0_rkind
        !real(rkind) :: Cv    != 390.0_rkind
        !real(rkind) :: T0    != 300.0_rkind
        real(rkind) :: gams  != 2.0_rkind

        ! EOS (2.83)-(2.84) in Barton PhD Thesis (2010)
        ! Default values given in Table (2.1) in Barton PhD Thesis (2010)
        real(rkind) :: B0 = 2.1d3**2
        !real(rkind) :: beta  != 3.0_rkind
        !real(rkind) :: Kp    != c0^2-(4/3)B0^2 = 15.28d6
        !real(rkind) :: alp   != 1.0_rkind
        !real(rkind) :: Cv    != 390.0_rkind
        !real(rkind) :: T0    != 300.0_rkind
        !real(rkind) :: gams  != 2.0_rkind
        

        real(rkind) :: invCv

        real(rkind), allocatable, dimension(:,:,:,:) :: finger
        real(rkind), allocatable, dimension(:,:,:,:) :: fingersq
        real(rkind), allocatable, dimension(:,:,:)   :: Inv1,Inv2,Inv3,dedI1fac,dedI2fac,dedI3fac,GI3,GpI3
        real(rkind), allocatable, dimension(:,:,:)   :: Inv1G,Inv2Gfac,Inv3G

        real(rkind), allocatable, dimension(:) :: work
        integer :: lwork

        integer :: eostype

    contains

        procedure :: init
        procedure :: get_p_devstress_T_sos
        procedure :: get_e_from_rhoT
        !procedure :: get_T
        !procedure :: get_sos
        procedure :: destroy
 
    end type

contains

    subroutine init(this,decomp,eostype,eosparams)
        class(generaleos),         intent(inout) :: this
        type(decomp_info),         intent(in)    :: decomp
        integer,                   intent(in)    :: eostype
        real(rkind), dimension(:), intent(in)    :: eosparams

        real(rkind) :: gloc(3,3), eigval(3)
        integer :: info

        this%eostype = eostype

        if(this%eostype==2) then

            ! EOS (77)-(84) in Hill et al., JCP 229 (2010). Cauchy stress
            ! derived using Cayley-Hamilton theorem, general form written in eq.
            ! (21) in Hill et al. (2010).

            this%K0   = eosparams(1);     this%Kp   = eosparams(2)
            this%G0   = eosparams(3);     this%Gp   = eosparams(4)
            this%beta = eosparams(5);     this%T0   = eosparams(6)
            this%Cv   = eosparams(7);     this%gam0 = eosparams(8)
            this%qpar = eosparams(9)
    
            this%invCv = one/this%Cv
    
            allocate( this%finger  (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),6) )
            allocate( this%fingersq(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),6) )
            allocate( this%Inv1    (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
            allocate( this%Inv2    (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
            allocate( this%Inv3    (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
            allocate( this%dedI1fac(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
            allocate( this%dedI2fac(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
            allocate( this%dedI3fac(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
            allocate( this%GI3     (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
            allocate( this%GpI3    (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )

        elseif(this%eostype==3) then

            ! Godunov-Romenskii EOS (3.27)-(3.29), (6.3) in Lopez-Ortega PhD Thesis. Very similar
            ! (but not exactly identical) to eqs. (36)-(43) in Barton et al., JCP, 240 (2013); Also
            ! found in Lopez-Ortega et al., JCP, 257 (2014), eq. (41).

            this%mu0  = eosparams(1);     this%beta = eosparams(2); 
            this%Kp   = eosparams(3);     this%alp  = eosparams(4);
            this%Cv   = eosparams(5);     this%T0   = eosparams(6);
            this%gams = eosparams(7);

            ! Get optimal lwork for eigenvalue computation
            this%lwork = -1;         allocate(this%work(100))   ! arbitrary large size
            call dsyev('V', 'U', 3, gloc, 3, eigval, this%work, this%lwork, info)
            this%lwork = this%work(1)
            deallocate(this%work); allocate(this%work(this%lwork))

            allocate( this%finger  (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),6) )

        elseif(this%eostype==4) then

            ! EOS (2.82)-(2.83), in Barton PhD Thesis.

            this%B0   = eosparams(1);     this%beta = eosparams(2); 
            this%Kp   = eosparams(3);     this%alp  = eosparams(4);
            this%Cv   = eosparams(5);     this%T0   = eosparams(6);
            this%gams = eosparams(7);

            !! Get optimal lwork for eigenvalue computation
            !this%lwork = -1;         allocate(this%work(100))   ! arbitrary large size
            !call dsyev('V', 'U', 3, gloc, 3, eigval, this%work, this%lwork, info)
            !this%lwork = this%work(1)
            !deallocate(this%work); allocate(this%work(this%lwork))

            allocate( this%finger  (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),6) )
            allocate( this%fingersq(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),6) )
            allocate( this%Inv1G   (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
            allocate( this%Inv2Gfac(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
            allocate( this%Inv3G   (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
        endif

    end subroutine

    subroutine destroy(this)
        class(generaleos),         intent(inout) :: this

         
        if(this%eostype==2) then
            deallocate(this%finger,this%fingersq,this%Inv1,this%Inv2,this%Inv3,this%dedI1fac,this%dedI2fac,this%dedI3fac,this%GI3,this%GpI3)
        elseif(this%eostype==3) then
            deallocate(this%finger, this%work)
        elseif(this%eostype==4) then
            deallocate(this%Inv3G,this%Inv2Gfac,this%Inv1G,this%fingersq,this%finger)
        endif
    end subroutine

    subroutine get_p_devstress_T_sos(this,rho0,g,rho,e,entr,p,T,devstress,sos)
        use decomp_2d, only: nrank
        class(generaleos), intent(inout) :: this
        real(rkind),                     intent(in)  :: rho0
        real(rkind), dimension(:,:,:,:), intent(in)  :: g
        real(rkind), dimension(:,:,:),   intent(in)  :: rho, e
        real(rkind), dimension(:,:,:),   intent(out) :: entr, p, T, sos
        real(rkind), dimension(:,:,:,:), intent(out) :: devstress

        integer :: i, j, k, nxp, nyp, nzp, info, idebug = 0
        real(rkind) :: gloc(3,3), eigval(3), lam(3,3), hencky(3,3), mu0Byrho0, OneByAlpP2, OneByAlpP1, BetP1
        real(rkind) :: KbyAlp, Kby2Alp2, CvT0, gamFac,  ploc, I1He, expmI1He, cs2, mu, eshear, ethermal, dehdI1He, expmI1Healp, expmI1Hegam
        real(rkind) :: Alpby2, Betby2, Gamby2, B0by2, invT0, gamCv, B0fac, betFac, betFac2

        if(this%eostype==2) then
          !p = zero!(this%gam-one)*rho*e - this%gam*this%PInf
          !devstress = zero!(this%gam-one)*rho*e - this%gam*this%PInf

          !associate ( g11 => g(:,:,:,1), g12 => g(:,:,:,2), g13 => g(:,:,:,3), &
          !            g21 => g(:,:,:,4), g22 => g(:,:,:,5), g23 => g(:,:,:,6), &
          !            g31 => g(:,:,:,7), g32 => g(:,:,:,8), g33 => g(:,:,:,9)  )

          !    ! compute g*gT first
          !    this%finger(:,:,:,1) = g11*g11 + g12*g12 + g13*g13
          !    this%finger(:,:,:,2) = g11*g21 + g12*g22 + g13*g23
          !    this%finger(:,:,:,3) = g11*g31 + g12*g32 + g13*g33
          !    this%finger(:,:,:,4) = g21*g21 + g22*g22 + g23*g23
          !    this%finger(:,:,:,5) = g21*g31 + g22*g32 + g23*g33
          !    this%finger(:,:,:,6) = g31*g31 + g32*g32 + g33*g33

          !    ! store tr(finger) in Inv2
          !    this%Inv2 = this%finger(:,:,:,1) + this%finger(:,:,:,4) + this%finger(:,:,:,6)
          !    

          !    ! compute inverse of (g*gT) and store in finger. This is right Cauchy-Green tensor, C ( = g^(-T)*g^(-1) = FT*F) )

          !    ! compute invariants of C
          !    this%Inv1 = this%finger(:,:,:,1) + this%finger(:,:,:,4) + this%finger(:,:,:,6)

          !    this%Inv3 = this%finger(:,:,:,1) * (this%finger(:,:,:,4)*this%finger(:,:,:,6) - this%finger(:,:,:,5)*this%finger(:,:,:,5)) &
          !              + this%finger(:,:,:,2) * (this%finger(:,:,:,3)*this%finger(:,:,:,5) - this%finger(:,:,:,6)*this%finger(:,:,:,2)) &
          !              + this%finger(:,:,:,3) * (this%finger(:,:,:,5)*this%finger(:,:,:,2) - this%finger(:,:,:,3)*this%finger(:,:,:,4)) 

          !    this%Inv2 = this%Inv2*this%Inv3

          !    ! compute gT*g first
          !    this%finger(:,:,:,1) = g11*g11 + g21*g21 + g31*g31
          !    this%finger(:,:,:,2) = g11*g12 + g21*g22 + g31*g32
          !    this%finger(:,:,:,3) = g11*g13 + g21*g23 + g31*g33
          !    this%finger(:,:,:,4) = g12*g12 + g22*g22 + g32*g32
          !    this%finger(:,:,:,5) = g12*g13 + g22*g23 + g32*g33
          !    this%finger(:,:,:,6) = g13*g13 + g23*g23 + g33*g33

          !end associate

          !! compute inverse of (gT*g) and store in finger. This is left Cauchy-Green tensor, b ( = g^(-1)*g^(-T) = F*FT)
          !    

          !associate ( GG11 => this%finger(:,:,:,1), GG12 => this%finger(:,:,:,2), GG13 => this%finger(:,:,:,3), &
          !            GG21 => this%finger(:,:,:,2), GG22 => this%finger(:,:,:,4), GG23 => this%finger(:,:,:,5), &
          !            GG31 => this%finger(:,:,:,3), GG32 => this%finger(:,:,:,5), GG33 => this%finger(:,:,:,6)  )

          !        ! compute square of b and store in fingersq
          !        this%fingersq(:,:,:,1) = GG11*GG11 + GG12*GG21 + GG13*GG31
          !        this%fingersq(:,:,:,2) = GG11*GG12 + GG12*GG22 + GG13*GG32
          !        this%fingersq(:,:,:,3) = GG11*GG13 + GG12*GG23 + GG13*GG33
          !        this%fingersq(:,:,:,4) = GG21*GG12 + GG22*GG22 + GG23*GG32
          !        this%fingersq(:,:,:,5) = GG21*GG13 + GG22*GG23 + GG23*GG33
          !        this%fingersq(:,:,:,6) = GG31*GG13 + GG32*GG23 + GG33*GG33
          !end associate

          !this%dedI1fac = this%beta * this%GI3 * this%Inv3**(-third)
          !this%dedI2fac = (one-this%beta) * this%GI3 * this%Inv3**(-twothird)

          !! 2*rho*dedI3*I3
          !this%dedI3fac = three*this%K0*exp(-1.5D0*(this%Kp-one)*(this%Inv3**(sixth)-one))*(this%Inv3**(-sixth)-this%Inv3**(-third)) &
          !         - (two*rho0*this%Cv*this%T0*this%gam0**2/this%qpar) * this%Inv3**(half*(one-this%qpar)) * (one - this%Inv3**this%qpar) * &
          !           (exp(entr)-one) * exp(this%gam0/this%qpar*(one-this%Inv3**this%qpar))      &
          !         + this%GpI3 * (this%beta*this%Inv1*this%Inv3**(-third) + (one-this%beta)*this%Inv2*this%Inv3**(-twothird) - three) &
          !         + this%GI3/this%Inv3 * (this%beta*sixth*this%Inv1*this%Inv3**(-third) - 1.5_rkind*(one-this%beta)*this%Inv2*this%Inv3**(-twothird)-three)

          !do i = 1, 6
          !    devstress(:,:,:,i) = -this%dedI2fac * this%fingersq(:,:,:,i) + (this%dedI1fac + this%Inv1*this%dedI2fac)*this%finger(:,:,:,i)
          !enddo
          !p = -third*(devstress(:,:,:,1) + devstress(:,:,:,4) + devstress(:,:,:,6))
          !devstress(:,:,:,1) = devstress(:,:,:,1) + p 
          !devstress(:,:,:,4) = devstress(:,:,:,4) + p 
          !devstress(:,:,:,6) = devstress(:,:,:,6) + p 

          !p = p + this%Inv3*this%dedI3fac


        elseif(this%eostype==3) then

            nxp = size(g,1); nyp = size(g,2); nzp = size(g,3);

            associate ( g11 => g(:,:,:,1), g12 => g(:,:,:,2), g13 => g(:,:,:,3), &
                        g21 => g(:,:,:,4), g22 => g(:,:,:,5), g23 => g(:,:,:,6), &
                        g31 => g(:,:,:,7), g32 => g(:,:,:,8), g33 => g(:,:,:,9)  )

                ! compute gT*g first
                this%finger(:,:,:,1) = g11*g11 + g21*g21 + g31*g31
                this%finger(:,:,:,2) = g11*g12 + g21*g22 + g31*g32
                this%finger(:,:,:,3) = g11*g13 + g21*g23 + g31*g33
                this%finger(:,:,:,4) = g12*g12 + g22*g22 + g32*g32
                this%finger(:,:,:,5) = g12*g13 + g22*g23 + g32*g33
                this%finger(:,:,:,6) = g13*g13 + g23*g23 + g33*g33

            end associate

            KbyAlp = this%Kp/this%alp
            Kby2Alp2 = half*this%Kp/this%alp**2
            CvT0 = this%Cv*this%T0
            mu0Byrho0 = this%mu0/rho0
            OneByAlpP2 = one/this%alp + two
            OneByAlpP1 = one/this%alp + one
            gamFac = this%gams * (this%gams + one)
            BetP1 = this%beta + one

            do k = 1,nzp
              do j = 1,nyp
                do i = 1,nxp

                    ! Get eigenvalues and eigenvectors of (gT*g)
                    gloc(1,1) = this%finger(i,j,k,1); gloc(1,2) = this%finger(i,j,k,2); gloc(1,3) = this%finger(i,j,k,3)
                                                      gloc(2,2) = this%finger(i,j,k,4); gloc(2,3) = this%finger(i,j,k,5)
                                                                                        gloc(3,3) = this%finger(i,j,k,6)
                    call dsyev('V', 'U', 3, gloc, 3, eigval, this%work, this%lwork, info)
                    if(info .ne. 0) print '(A,I6,A)', 'proc ', nrank, ': Problem with EV. Please check.'
                    if(minval(eigval)<1.0d-12) then
                        print '(A,I6,A)', 'proc ', nrank, ': Matrix not SPD. Please check.'
                    endif

                    if(idebug==1) write(*,*) '---eigval---'
                    if(idebug==1) write(*,'(3(e19.12,1x))') eigval(1:3)
                    if(idebug==1) write(*,*) '---eigvev---'
                    if(idebug==1) write(*,'(3(e19.12,1x))') gloc(1,1:3)
                    if(idebug==1) write(*,'(3(e19.12,1x))') gloc(2,1:3)
                    if(idebug==1) write(*,'(3(e19.12,1x))') gloc(3,1:3)

                    ! compute Hencky strain matrix: log(F*FT) = log((gT*g)^(-1)) = P*(-log(Lam^))*PT
                    lam = zero; lam(1,1) = -log(eigval(1)); lam(2,2) = -log(eigval(2)); lam(3,3) = -log(eigval(3))
                    hencky = matmul(gloc,lam)
                    if(idebug==1) write(*,*) '---hencky1---'
                    if(idebug==1) write(*,'(3(e19.12,1x))') hencky(1,1:3)
                    if(idebug==1) write(*,'(3(e19.12,1x))') hencky(2,1:3)
                    if(idebug==1) write(*,'(3(e19.12,1x))') hencky(3,1:3)
                    hencky = matmul(hencky,transpose(gloc))
                    hencky = half*hencky
                    if(idebug==1) write(*,*) '---hencky ---'
                    if(idebug==1) write(*,'(3(e19.12,1x))') hencky(1,1:3)
                    if(idebug==1) write(*,'(3(e19.12,1x))') hencky(2,1:3)
                    if(idebug==1) write(*,'(3(e19.12,1x))') hencky(3,1:3)

                    ! make this traceless
                    ploc = third*(hencky(1,1)+hencky(2,2)+hencky(3,3))
                    hencky(1,1) = hencky(1,1) - ploc
                    hencky(2,2) = hencky(2,2) - ploc
                    hencky(3,3) = hencky(3,3) - ploc
                    if(idebug==1) write(*,*) '---traceless hencky ---'
                    if(idebug==1) write(*,'(3(e19.12,1x))') hencky(1,1:3)
                    if(idebug==1) write(*,'(3(e19.12,1x))') hencky(2,1:3)
                    if(idebug==1) write(*,'(3(e19.12,1x))') hencky(3,1:3)

                    ! compute entropy from energy
                    expmI1He = rho(i,j,k)/rho0
                    I1He = -log(expmI1He)
                    !I1He = three*ploc
                    !expmI1He = exp(-I1He)
                    expmI1Healp = expmI1He**this%alp
                    expmI1Hegam = expmI1He**this%gams
                    cs2 = mu0Byrho0*expmI1He**this%beta
                    mu = rho(i,j,k)*cs2
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'expmI1He   :', expmI1He
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'I1He       :', I1He
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'expmI1Healp:', expmI1Healp
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'expmI1Hegam:', expmI1Hegam
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'cs2         :', cs2
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'mu          :', mu

                    eshear = cs2*sum(hencky**2)
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'eshear     :', eshear
                    ethermal = e(i,j,k) - eshear - Kby2Alp2*(expmI1Healp-one)**2
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'ethermal   :', ethermal
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'ehydro     :', Kby2Alp2*(expmI1Healp-one)**2
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'etotal     :', e(i,j,k)                     
                    !entr(i,j,k) = this%Cv*log(ethermal/(CvT0*expmI1Hegam) + one)
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'entr       :', entr(i,j,k)
                    dehdI1He = -KbyAlp*expmI1Healp*(expmI1Healp - one) - this%gams * ethermal
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'dehdI1He   :', dehdI1He

                    devstress(i,j,k,1) = two*mu*hencky(1,1)
                    devstress(i,j,k,2) = two*mu*hencky(1,2)
                    devstress(i,j,k,3) = two*mu*hencky(1,3)
                    devstress(i,j,k,4) = two*mu*hencky(2,2)
                    devstress(i,j,k,5) = two*mu*hencky(2,3)
                    devstress(i,j,k,6) = two*mu*hencky(3,3)

                    !p(i,j,k) = -(two*mu*ploc + rho(i,j,k) * dehdI1He - twothird*mu*I1He)
                    p(i,j,k) = -rho(i,j,k) * dehdI1He
                    T(i,j,k) = ethermal/this%Cv + this%T0*expmI1Hegam
                    sos(i,j,k) = this%Kp*expmI1Healp*(OneByAlpP2*expmI1Healp-OneByAlpP1) + &
                                      gamFac*ethermal + two*cs2*(one-BetP1*I1He)
                    if(sos(i,j,k) < 1.0d-14) then
                    else
                      sos(i,j,k) = sqrt(sos(i,j,k))
                    endif
                    if(idebug==1) write(*,'(a,1x,e19.12)') 'p      :', p(i,j,k)
                enddo
              enddo
            enddo
        elseif(this%eostype==4) then

          associate ( g11 => g(:,:,:,1), g12 => g(:,:,:,2), g13 => g(:,:,:,3), &
                      g21 => g(:,:,:,4), g22 => g(:,:,:,5), g23 => g(:,:,:,6), &
                      g31 => g(:,:,:,7), g32 => g(:,:,:,8), g33 => g(:,:,:,9)  )
  
              this%finger(:,:,:,1) = g11*g11 + g21*g21 + g31*g31
              this%finger(:,:,:,2) = g11*g12 + g21*g22 + g31*g32
              this%finger(:,:,:,3) = g11*g13 + g21*g23 + g31*g33
              this%finger(:,:,:,4) = g12*g12 + g22*g22 + g32*g32
              this%finger(:,:,:,5) = g12*g13 + g22*g23 + g32*g33
              this%finger(:,:,:,6) = g13*g13 + g23*g23 + g33*g33
  
          end associate
  
          associate ( GG11 => this%finger(:,:,:,1), GG12 => this%finger(:,:,:,2), GG13 => this%finger(:,:,:,3), &
                      GG21 => this%finger(:,:,:,2), GG22 => this%finger(:,:,:,4), GG23 => this%finger(:,:,:,5), &
                      GG31 => this%finger(:,:,:,3), GG32 => this%finger(:,:,:,5), GG33 => this%finger(:,:,:,6)  )
  
              this%fingersq(:,:,:,1) = GG11*GG11 + GG12*GG21 + GG13*GG31
              this%fingersq(:,:,:,2) = GG11*GG12 + GG12*GG22 + GG13*GG32
              this%fingersq(:,:,:,3) = GG11*GG13 + GG12*GG23 + GG13*GG33
              this%fingersq(:,:,:,4) = GG21*GG12 + GG22*GG22 + GG23*GG32
              this%fingersq(:,:,:,5) = GG21*GG13 + GG22*GG23 + GG23*GG33
              this%fingersq(:,:,:,6) = GG31*GG13 + GG32*GG23 + GG33*GG33
  
              this%Inv1G = GG11 + GG22 + GG33
              this%Inv3G = GG11*(GG22*GG33-GG23*GG32) - GG12*(GG21*GG33-GG31*GG23) + GG13*(GG21*GG32-GG31*GG22)
              this%Inv2Gfac = -sixth*this%Inv1G*this%Inv1G + half*(this%fingersq(:,:,:,1) + this%fingersq(:,:,:,4) + this%fingersq(:,:,:,6))
          end associate
  
          Kby2Alp2 = half*this%Kp/this%alp**2
          Alpby2 = half*this%alp
          Gamby2 = half*this%gams
          B0by2  = half*this%B0
          Betby2 = half*this%beta
          invT0  = one/this%T0
          gamCv  = this%gams*this%Cv
          B0fac  = (half*this%beta + twothird)*this%B0 
          gamFac = this%gams**2*this%Cv
          betFac = half*this%beta+twothird
          betFac2 = this%beta*betfac
          KbyAlp = this%Kp/this%alp

          T = (e - (Kby2Alp2 * (this%Inv3G**Alpby2 - one)**2 + B0by2*this%Inv3G**Betby2*this%Inv2Gfac)) / this%Cv + this%T0*this%Inv3G**Gamby2
          entr = this%Cv*log(invT0*T/this%Inv3G**Gamby2)
          p = rho * (KbyAlp * (this%Inv3G**Alpby2 - one)*this%Inv3G**Alpby2 + gamCv*(T - this%T0*this%Inv3G**Gamby2) + B0fac*this%Inv3G**Betby2*this%Inv2Gfac)
          devstress(:,:,:,1) = two*this%B0*rho*this%Inv3G**Betby2 * (third*this%Inv2Gfac + sixth*this%Inv1G*this%finger(:,:,:,1) - half*this%fingersq(:,:,:,1))
          devstress(:,:,:,2) = two*this%B0*rho*this%Inv3G**Betby2 * (                      sixth*this%Inv1G*this%finger(:,:,:,2) - half*this%fingersq(:,:,:,2))
          devstress(:,:,:,3) = two*this%B0*rho*this%Inv3G**Betby2 * (                      sixth*this%Inv1G*this%finger(:,:,:,3) - half*this%fingersq(:,:,:,3))
          devstress(:,:,:,4) = two*this%B0*rho*this%Inv3G**Betby2 * (third*this%Inv2Gfac + sixth*this%Inv1G*this%finger(:,:,:,4) - half*this%fingersq(:,:,:,4))
          devstress(:,:,:,5) = two*this%B0*rho*this%Inv3G**Betby2 * (                      sixth*this%Inv1G*this%finger(:,:,:,5) - half*this%fingersq(:,:,:,5))
          devstress(:,:,:,6) = two*this%B0*rho*this%Inv3G**Betby2 * (third*this%Inv2Gfac + sixth*this%Inv1G*this%finger(:,:,:,6) - half*this%fingersq(:,:,:,6))
          sos = (p-devstress(:,:,:,1))/rho + this%Kp*this%Inv3G**Alpby2*(two*this%Inv3G**Alpby2 - one) + gamFac*(T - this%T0*this%Inv3G**Gamby2) + &
                third*this%B0*this%Inv3G**Betby2 * (betFac2*(this%Inv3G - one)**2 - two*(this%Inv3G - one) + 4*this%Inv3G*(betFac*sqrt(this%Inv3G) - two))
          if(minval(sos) < 1.0d-14) then
          else
            sos = sqrt(sos)
          endif
        endif

    end subroutine

    subroutine get_e_from_rhoT(this,rho0,g,rho,T,e)
        class(generaleos), intent(inout) :: this
        real(rkind),                     intent(in)   :: rho0
        real(rkind), dimension(:,:,:,:), intent(in)   :: g
        real(rkind), dimension(:,:,:),   intent(in)   :: T, rho
        real(rkind), dimension(:,:,:),   intent(out)  :: e

        integer :: i, j, k, nxp, nyp, nzp, info, idebug = 0
        real(rkind) :: I1He, expmI1He, expmI1Healp, expmI1Hegam, cs2, ehydro, ethermal, eshear
        real(rkind) :: Kby2Alp2, gloc(3,3), eigval(3), lam(3,3), hencky(3,3), ploc, mu0Byrho0
        real(rkind) :: Alpby2, Betby2, Gamby2, B0by2

        nxp = size(g,1); nyp = size(g,2); nzp = size(g,3);

        if(this%eostype==2) then
        elseif(this%eostype==3) then
          
          associate ( g11 => g(:,:,:,1), g12 => g(:,:,:,2), g13 => g(:,:,:,3), &
                      g21 => g(:,:,:,4), g22 => g(:,:,:,5), g23 => g(:,:,:,6), &
                      g31 => g(:,:,:,7), g32 => g(:,:,:,8), g33 => g(:,:,:,9)  )

              ! compute gT*g first
              this%finger(:,:,:,1) = g11*g11 + g21*g21 + g31*g31
              this%finger(:,:,:,2) = g11*g12 + g21*g22 + g31*g32
              this%finger(:,:,:,3) = g11*g13 + g21*g23 + g31*g33
              this%finger(:,:,:,4) = g12*g12 + g22*g22 + g32*g32
              this%finger(:,:,:,5) = g12*g13 + g22*g23 + g32*g33
              this%finger(:,:,:,6) = g13*g13 + g23*g23 + g33*g33

          end associate

          Kby2Alp2 = half*this%Kp/this%alp**2
          mu0Byrho0 = this%mu0/rho0

          do k = 1,nzp
            do j = 1,nyp
              do i = 1,nxp

                  ! Get eigenvalues and eigenvectors of (gT*g)
                  gloc(1,1) = this%finger(i,j,k,1); gloc(1,2) = this%finger(i,j,k,2); gloc(1,3) = this%finger(i,j,k,3)
                                                    gloc(2,2) = this%finger(i,j,k,4); gloc(2,3) = this%finger(i,j,k,5)
                                                                                      gloc(3,3) = this%finger(i,j,k,6)
                  call dsyev('V', 'U', 3, gloc, 3, eigval, this%work, this%lwork, info)
                  if(info .ne. 0) print '(A,I6,A)', 'proc ', nrank, ': Problem with EV. Please check.'
                  if(minval(eigval)<1.0d-12) then
                      print '(A,I6,A)', 'proc ', nrank, ': Matrix not SPD. Please check.'
                  endif

                  ! compute Hencky strain matrix: log(F*FT) = log((gT*g)^(-1)) = P*(Lam^(-1))*PT
                  lam = zero; lam(1,1) = -log(eigval(1)); lam(2,2) = -log(eigval(2)); lam(3,3) = -log(eigval(3))
                  hencky = matmul(gloc,lam)
                  hencky = matmul(hencky,transpose(gloc))
                  hencky = half*hencky

                  ! make this traceless
                  ploc = third*(hencky(1,1)+hencky(2,2)+hencky(3,3))
                  hencky(1,1) = hencky(1,1) - ploc
                  hencky(2,2) = hencky(2,2) - ploc
                  hencky(3,3) = hencky(3,3) - ploc

                  !I1He = three*ploc
                  !expmI1He = exp(-I1He)
                  expmI1He = rho(i,j,k)/rho0
                  I1He = -log(expmI1He)
                  expmI1Healp = expmI1He**this%alp
                  expmI1Hegam = expmI1He**this%gams
                  cs2 = mu0Byrho0*expmI1He**this%beta

                  ehydro = Kby2Alp2*(expmI1Healp-one)**2
                  ethermal = this%Cv*(T(i,j,k) - this%T0*expmI1Hegam)
                  eshear = cs2*sum(hencky**2)
                  e(i,j,k) = ehydro + ethermal + eshear
                  if(idebug==1) write(*,'(a,1x,e19.12)') 'eshear     :', eshear
                  if(idebug==1) write(*,'(a,1x,e19.12)') 'ethermal   :', ethermal
                  if(idebug==1) write(*,'(a,1x,e19.12)') 'ehydro     :', ehydro
              enddo
            enddo
          enddo
        elseif(this%eostype==4) then

          associate ( g11 => g(:,:,:,1), g12 => g(:,:,:,2), g13 => g(:,:,:,3), &
                      g21 => g(:,:,:,4), g22 => g(:,:,:,5), g23 => g(:,:,:,6), &
                      g31 => g(:,:,:,7), g32 => g(:,:,:,8), g33 => g(:,:,:,9)  )
  
              this%finger(:,:,:,1) = g11*g11 + g21*g21 + g31*g31
              this%finger(:,:,:,2) = g11*g12 + g21*g22 + g31*g32
              this%finger(:,:,:,3) = g11*g13 + g21*g23 + g31*g33
              this%finger(:,:,:,4) = g12*g12 + g22*g22 + g32*g32
              this%finger(:,:,:,5) = g12*g13 + g22*g23 + g32*g33
              this%finger(:,:,:,6) = g13*g13 + g23*g23 + g33*g33
  
          end associate
  
          associate ( GG11 => this%finger(:,:,:,1), GG12 => this%finger(:,:,:,2), GG13 => this%finger(:,:,:,3), &
                      GG21 => this%finger(:,:,:,2), GG22 => this%finger(:,:,:,4), GG23 => this%finger(:,:,:,5), &
                      GG31 => this%finger(:,:,:,3), GG32 => this%finger(:,:,:,5), GG33 => this%finger(:,:,:,6)  )
  
              this%fingersq(:,:,:,1) = GG11*GG11 + GG12*GG21 + GG13*GG31
              this%fingersq(:,:,:,2) = GG11*GG12 + GG12*GG22 + GG13*GG32
              this%fingersq(:,:,:,3) = GG11*GG13 + GG12*GG23 + GG13*GG33
              this%fingersq(:,:,:,4) = GG21*GG12 + GG22*GG22 + GG23*GG32
              this%fingersq(:,:,:,5) = GG21*GG13 + GG22*GG23 + GG23*GG33
              this%fingersq(:,:,:,6) = GG31*GG13 + GG32*GG23 + GG33*GG33
  
              this%Inv1G = GG11 + GG22 + GG33
              this%Inv3G = GG11*(GG22*GG33-GG23*GG32) - GG12*(GG21*GG33-GG31*GG23) + GG13*(GG21*GG32-GG31*GG22)
              this%Inv2Gfac = -sixth*this%Inv1G*this%Inv1G + half*(this%fingersq(:,:,:,1) + this%fingersq(:,:,:,4) + this%fingersq(:,:,:,6))
          end associate
  
          Kby2Alp2 = half*this%Kp/this%alp**2
          Alpby2 = half*this%alp
          Gamby2 = half*this%gams
          B0by2  = half*this%B0
          Betby2 = half*this%beta
  
          e = Kby2Alp2 * (this%Inv3G**Alpby2 - one)**2 + this%Cv*(T - this%T0*this%Inv3G**Gamby2) + B0by2*this%Inv3G**Betby2*this%Inv2Gfac
        endif

    end subroutine

    !pure subroutine get_T(this,e,T)
    !    class(generaleos), intent(in) :: this
    !    real(rkind), dimension(:,:,:), intent(in)  :: e
    !    real(rkind), dimension(:,:,:), intent(out) :: T

    !    T = this%invCv*e

    !end subroutine

    !pure subroutine get_sos(this,rho0,rho,devstress,sos)
    !    class(generaleos), intent(in) :: this
    !    real(rkind),                     intent(in)  :: rho0
    !    real(rkind), dimension(:,:,:),   intent(in)  :: rho
    !    real(rkind), dimension(:,:,:,:), intent(in)  :: devstress
    !    real(rkind), dimension(:,:,:),   intent(out) :: sos

    !    sos = zero!sqrt(this%gam*(p+this%PInf)/rho)

    !end subroutine

!    pure subroutine get_e_from_p(this,rho,p,e)
!        class(generaleos), intent(in) :: this
!        real(rkind), dimension(:,:,:), intent(in)  :: rho,p
!        real(rkind), dimension(:,:,:), intent(out) :: e
!
!        e = zero!(p + this%gam*this%PInf) * this%onebygam_m1 / rho
!
!    end subroutine

end module
