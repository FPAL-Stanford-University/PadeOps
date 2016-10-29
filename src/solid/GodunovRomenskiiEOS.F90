module GodunovRomenskiiEOSMod

    use kind_parameters, only: rkind,clen
    use constants,       only: one
    use exits,           only: GracefulExit
    use FingerEOSMod,    only: fingereos

    implicit none

    ! This class extends the fingereos class and implements the Godunov Romenskii EOS as defined in Barton (2009), Section 2.6.
    type, extends(fingereos) :: godromeos

        real(rkind) :: K0
        real(rkind) :: Cv
        real(rkind) :: T0
        real(rkind) :: B0
        real(rkind) :: alpha
        real(rkind) :: beta
        real(rkind) :: gam

    contains

        procedure :: get_energy_derivatives
        procedure :: get_e_from_rho_invariants_T

    end type

    interface godromeos
        module procedure init
    end interface

contains

    function init(rho0,K0,Cv,T0,B0,alpha,beta,gam,usegTg) result(this)
        type(godromeos)           :: this
        real(rkind),   intent(in) :: rho0, K0, Cv, T0, B0, alpha, beta, gam
        logical,       intent(in) :: usegTg

        this%rho0 = rho0
        this%K0 = K0
        this%Cv = Cv
        this%T0 = T0
        this%B0 = B0
        this%alpha = alpha
        this%beta = beta
        this%gam = gam
        this%usegTg = usegTg

    end function

    ! Subroutine to compute the dervatives of the internal energy with respect to the invariants of the finger tensor and
    ! entropy (in the form of the temperature)
    subroutine get_energy_derivatives(this,rho,e,I1,I2,I3,dedI1,dedI2,dedI3,T)
        use constants, only: half, one, two, three, four
        class(godromeos),                                            intent(in)  :: this
        real(rkind), dimension(:,:,:),                               intent(in)  :: rho
        real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)), intent(in)  :: e, I1, I2, I3
        real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)), intent(out) :: dedI1, dedI2, dedI3, T
        
        real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)) :: entrexp

        entrexp = e
        entrexp = entrexp - half*this%B0*(I3**(half*this%beta))*(I1**2/three - I2) ! U(I3,entr)
        entrexp = one + ( entrexp - ( this%K0/(two*this%alpha**2) ) * ( I3**(half*this%alpha) - one )**2 ) / ( this%Cv * this%T0 * I3**(half*this%gam) ) ! e^(entropy/Cv)

        T = this%T0 * I3**(half*this%gam) * entrexp ! d e / d entropy

        dedI1 = (this%B0/three) * I3**(half*this%beta) * I1
        dedI2 =  - half*this%B0 * I3**(half*this%beta)
        dedI3 = (this%B0 * this%beta / four) * I3**(half*this%beta - one) * ( I1**2 / three - I2 ) &
              + (half*this%K0/this%alpha) * (I3**(half*this%alpha) - one) * I3**(half*this%alpha - one) &
              + (half*this%gam) * this%Cv * this%T0 * I3**(half*this%gam - one) * ( entrexp - one )

    end subroutine

    ! Subroutine to get the internal energy from the temperature and material density and invariants of the finger tensor
    subroutine get_e_from_rho_invariants_T(this,rho,I1,I2,I3,T,e)
        use constants, only: half, two, three
        class(godromeos),                                              intent(in)  :: this
        real(rkind), dimension(:,:,:),                                 intent(in)  :: rho
        real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)),   intent(in)  :: I1,I2,I3,T
        real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)),   intent(out) :: e

        e = T / ( this%T0 * I3**(half*this%gam) ) ! e^(entropy / Cv)

        e = ( this%K0/(two*this%alpha**2) ) * ( I3**(half*this%alpha) - one )**2 + this%Cv*this%T0 * I3**(half*this%gam) * (e-one) & ! U(I3,entropy)
          + half*this%B0 * I3**(half*this%beta) * ( I1**2/three - I2 )                                                               ! W(I1,I2,I3)

    end subroutine

end module
