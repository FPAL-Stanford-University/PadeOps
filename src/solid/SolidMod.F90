module SolidMod

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,one,two
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit
    use EOSMod,          only: eos
    use StiffGasEOS,     only: stiffgas
    use Sep1SolidEOS,    only: sep1solid

    implicit none

    type :: solid
        integer :: nxp, nyp, nzp

        class(stiffgas ), allocatable :: hydro
        class(sep1solid), allocatable :: elastic

        type(decomp_info), pointer :: decomp
        type(derivatives), pointer :: der

        real(rkind), dimension(:,:,:), allocatable :: Ys
        real(rkind), dimension(:,:,:), allocatable :: VF
        real(rkind), dimension(:,:,:), allocatable :: eh
        real(rkind), dimension(:,:,:), allocatable :: eel

        real(rkind), dimension(:,:,:,:), allocatable :: g
        real(rkind), dimension(:,:,:),   pointer     :: g11
        real(rkind), dimension(:,:,:),   pointer     :: g12
        real(rkind), dimension(:,:,:),   pointer     :: g13
        real(rkind), dimension(:,:,:),   pointer     :: g21
        real(rkind), dimension(:,:,:),   pointer     :: g22
        real(rkind), dimension(:,:,:),   pointer     :: g23
        real(rkind), dimension(:,:,:),   pointer     :: g31
        real(rkind), dimension(:,:,:),   pointer     :: g32
        real(rkind), dimension(:,:,:),   pointer     :: g33
        
        real(rkind), dimension(:,:,:,:), allocatable :: devstress
        real(rkind), dimension(:,:,:),   pointer     :: sxx
        real(rkind), dimension(:,:,:),   pointer     :: sxy
        real(rkind), dimension(:,:,:),   pointer     :: sxz
        real(rkind), dimension(:,:,:),   pointer     :: syy
        real(rkind), dimension(:,:,:),   pointer     :: syz
        real(rkind), dimension(:,:,:),   pointer     :: szz
        
        real(rkind), dimension(:,:,:),   allocatable :: p

        ! species-specific artificial properties
        real(rkind), dimension(:,:,:),   allocatable :: kap
        real(rkind), dimension(:,:,:,:), allocatable :: Ji

        ! species-specific conserved variables
        real(rkind), dimension(:,:,:,:), allocatable :: consrv

        ! work arrays
        real(rkind), dimension(:,:,:,:), allocatable :: Qtmpg
        real(rkind), dimension(:,:,:),   allocatable :: QtmpYs
        real(rkind), dimension(:,:,:),   allocatable :: Qtmpeh
        real(rkind), dimension(:,:,:),   allocatable :: QtmpVF

    contains

        procedure :: getRHS_g
        procedure :: getPlasticSources
        procedure :: getRHS_VF
        procedure :: getRHS_eh
        procedure :: update_g
        procedure :: update_Ys
        procedure :: update_eh
        procedure :: update_VF
        procedure :: filter
        final     :: destroy

    end type

    interface solid
        module procedure init
    end interface

contains

    function init(decomp,der) result(this)
        type(solid) :: this
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der

        this%decomp => decomp
        this%der => der
       
        ! Assume everything is in Y decomposition
        this%nxp = decomp%ysz(1)
        this%nyp = decomp%ysz(2)
        this%nzp = decomp%ysz(3)

        ! if (allocated(this%hydro)) deallocate(this%hydro)
        ! allocate( this%hydro, source=hydro )
        ! 
        ! if (allocated(this%elastic)) deallocate(this%elastic)
        ! allocate( this%elastic, source=elastic )
        
        ! Allocate material massfraction
        if( allocated( this%Ys ) ) deallocate( this%Ys )
        allocate( this%Ys(this%nxp,this%nyp,this%nzp) )
        
        ! Allocate material volume fraction
        if( allocated( this%VF ) ) deallocate( this%VF )
        allocate( this%VF(this%nxp,this%nyp,this%nzp) )
        
        ! Allocate material hydrodynamic energy
        if( allocated( this%eh ) ) deallocate( this%eh )
        allocate( this%eh(this%nxp,this%nyp,this%nzp) )
        
        ! Allocate material elastic energy
        if( allocated( this%eel ) ) deallocate( this%eel )
        allocate( this%eel(this%nxp,this%nyp,this%nzp) )
        
        ! Allocate material inverse elastic deformation gradients and associate pointers
        if( allocated( this%g ) ) deallocate( this%g )
        allocate( this%g(this%nxp,this%nyp,this%nzp,9) )
        this%g11 => this%g(:,:,:,1)   
        this%g12 => this%g(:,:,:,2)   
        this%g13 => this%g(:,:,:,3)   
        this%g21 => this%g(:,:,:,4)   
        this%g22 => this%g(:,:,:,5)   
        this%g23 => this%g(:,:,:,6)   
        this%g31 => this%g(:,:,:,7)   
        this%g32 => this%g(:,:,:,8)   
        this%g33 => this%g(:,:,:,9)   
        
        ! Allocate material deviatoric stress array
        if( allocated( this%devstress ) ) deallocate( this%devstress )
        allocate( this%devstress(this%nxp,this%nyp,this%nzp,6) )
        this%sxx  => this%devstress(:,:,:,1)   
        this%sxy  => this%devstress(:,:,:,2)   
        this%sxz  => this%devstress(:,:,:,3)   
        this%syy  => this%devstress(:,:,:,4)   
        this%syz  => this%devstress(:,:,:,5)   
        this%szz  => this%devstress(:,:,:,6)   
        
        ! Allocate material pressure array
        if( allocated( this%p ) ) deallocate( this%p )
        allocate( this%p(this%nxp,this%nyp,this%nzp) )

        ! Allocate material conserved variables
        if( allocated( this%consrv ) ) deallocate( this%consrv )
        allocate( this%consrv(this%nxp,this%nyp,this%nzp,2) )

        ! Allocate work arrays
        ! g tensor equation
        if( allocated( this%Qtmpg ) ) deallocate( this%Qtmpg )
        allocate( this%Qtmpg(this%nxp,this%nyp,this%nzp,9) )

        ! Ys equation
        if( allocated( this%QtmpYs ) ) deallocate( this%QtmpYs )
        allocate( this%QtmpYs(this%nxp,this%nyp,this%nzp) )

        ! eh equation
        if( allocated( this%Qtmpeh ) ) deallocate( this%Qtmpeh )
        allocate( this%Qtmpeh(this%nxp,this%nyp,this%nzp) )

        ! VF equation
        if( allocated( this%QtmpVF ) ) deallocate( this%QtmpVF )
        allocate( this%QtmpVF(this%nxp,this%nyp,this%nzp) )

    end function

    pure elemental subroutine destroy(this)
        class(solid), intent(inout) :: this

        ! First deallocate all the arrays
        if( allocated( this%Ys )        ) deallocate( this%Ys )
        if( allocated( this%VF )        ) deallocate( this%VF )
        if( allocated( this%eh )        ) deallocate( this%eh )
        if( allocated( this%g )         ) deallocate( this%g )
        if( allocated( this%devstress ) ) deallocate( this%devstress )
        if( allocated( this%p )         ) deallocate( this%p )

        nullify( this%g11 ); nullify( this%g12 ); nullify( this%g13 )
        nullify( this%g21 ); nullify( this%g22 ); nullify( this%g23 )
        nullify( this%g31 ); nullify( this%g32 ); nullify( this%g33 )

        nullify( this%sxx ); nullify( this%sxy ); nullify( this%sxz )
                             nullify( this%syy ); nullify( this%syz )
                                                  nullify( this%szz )

        if( allocated( this%consrv ) ) deallocate( this%consrv )
        if( allocated( this%eel )        ) deallocate( this%eel )

        if( allocated( this%Qtmpg ) ) deallocate( this%Qtmpg )
        if( allocated( this%QtmpYs ) ) deallocate( this%QtmpYs )
        if( allocated( this%Qtmpeh ) ) deallocate( this%Qtmpeh )
        if( allocated( this%QtmpVF ) ) deallocate( this%QtmpVF )

        ! Now deallocate the EOS objects
        if ( allocated(this%hydro)      ) deallocate(this%hydro)
        if ( allocated(this%elastic)    ) deallocate(this%elastic)

    end subroutine

    subroutine update_g(this,isub,dt,rho,u,v,w)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: rhsg  ! RHS for g tensor equation

        call this%getRHS_g(rho,u,v,w,rhsg)

        ! advance sub-step
        if(isub==1) this%Qtmpg = zero                   ! not really needed, since RK45_A(1) = 0
        this%Qtmpg  = dt*rhsg + RK45_A(isub)*this%Qtmpg
        this%g = this%g  + RK45_B(isub)*this%Qtmpg

    end subroutine

    subroutine getRHS_g(this,rho,u,v,w,rhsg)
        use operators, only: gradient, curl
        class(solid),                                         intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg

        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg
        real(rkind) :: etafac = one/6._rkind

        rhsg = zero

        detg = this%g11*(this%g22*this%g33-this%g23*this%g32) &
             - this%g12*(this%g21*this%g33-this%g31*this%g23) &
             + this%g13*(this%g21*this%g32-this%g31*this%g22)

        tmp = rho*this%Ys/(this%VF + real(1.0D-32,rkind))  ! Get the species density = rho*Y/VF
        penalty = etafac*( tmp/detg/this%elastic%rho0-one) ! Penalty term to keep g consistent with species density

        tmp = -u*this%g11-v*this%g12-w*this%g13
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,1),rhsg(:,:,:,2),rhsg(:,:,:,3))
        
        call curl(this%decomp, this%der, this%g11, this%g12, this%g13, curlg)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + v*curlg(:,:,:,3) - w*curlg(:,:,:,2) + penalty*this%g11
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + w*curlg(:,:,:,1) - u*curlg(:,:,:,3) + penalty*this%g12
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + u*curlg(:,:,:,2) - v*curlg(:,:,:,1) + penalty*this%g13
 
        tmp = -u*this%g21-v*this%g22-w*this%g23
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,4),rhsg(:,:,:,5),rhsg(:,:,:,6))
        
        call curl(this%decomp, this%der, this%g21, this%g22, this%g23, curlg)
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + v*curlg(:,:,:,3) - w*curlg(:,:,:,2) + penalty*this%g21
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + w*curlg(:,:,:,1) - u*curlg(:,:,:,3) + penalty*this%g22
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + u*curlg(:,:,:,2) - v*curlg(:,:,:,1) + penalty*this%g23
 
        tmp = -u*this%g31-v*this%g32-w*this%g33
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,7),rhsg(:,:,:,8),rhsg(:,:,:,9))

        call curl(this%decomp, this%der, this%g31, this%g32, this%g33, curlg)
        rhsg(:,:,:,7) = rhsg(:,:,:,7) + v*curlg(:,:,:,3) - w*curlg(:,:,:,2) + penalty*this%g31
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + w*curlg(:,:,:,1) - u*curlg(:,:,:,3) + penalty*this%g32
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + u*curlg(:,:,:,2) - v*curlg(:,:,:,1) + penalty*this%g33

        if (this%explPlast) then
            call this%getPlasticSources(detg,rhsg)
        end if

    end subroutine

    subroutine getPlasticSources(this,detg,rhsg)
        use constants, only: twothird
        class(solid),                                         intent(in)    :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)    :: detg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(inout) :: rhsg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: invtaurel

        ! Get S'S'
        invtaurel = this%sxx*this%sxx + two*this%sxy*this%sxy + two*this%sxz*this%sxz &
                                      +     this%syy*this%syy + two*this%syz*this%syz &
                                                              +     this%szz*this%szz

        ! 1/tau_rel
        invtaurel = this%invtau0 * ( invtaurel - (twothird)*this%elastic%yield**2 ) / this%elastic%mu**2
        where (invtaurel .LE. zero)
            invtaurel = zero
        end where
        invtaurel = invtaurel / (two * this%elastic%mu * detg)

        ! Add (1/tau_rel)*g*S to the rhsg (explicit plastic source terms)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + invtaurel * ( this%g11*this%sxx + this%g12*this%sxy + this%g13*this%sxz ) ! g11 
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + invtaurel * ( this%g11*this%sxy + this%g12*this%syy + this%g13*this%syz ) ! g12 
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + invtaurel * ( this%g11*this%sxz + this%g12*this%syz + this%g13*this%szz ) ! g13 
 
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + invtaurel * ( this%g21*this%sxx + this%g22*this%sxy + this%g23*this%sxz ) ! g21 
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + invtaurel * ( this%g21*this%sxy + this%g22*this%syy + this%g23*this%syz ) ! g22 
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + invtaurel * ( this%g21*this%sxz + this%g22*this%syz + this%g23*this%szz ) ! g23 

        rhsg(:,:,:,7) = rhsg(:,:,:,7) + invtaurel * ( this%g31*this%sxx + this%g32*this%sxy + this%g33*this%sxz ) ! g31 
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + invtaurel * ( this%g31*this%sxy + this%g32*this%syy + this%g33*this%syz ) ! g32 
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + invtaurel * ( this%g31*this%sxz + this%g32*this%syz + this%g33*this%szz ) ! g33 
    end subroutine

    subroutine update_Ys(this,isub,dt,rho,u,v,w)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhsYs  ! RHS for mass fraction equation

        call this%getRHS_Ys(rho,u,v,w,rhsYs)

        ! advance sub-step
        if(isub==1) this%QtmpYs = zero                   ! not really needed, since RK45_A(1) = 0
        this%QtmpYs  = dt*rhsYs + RK45_A(isub)*this%QtmpYs
        this%consrv(:,:,:,1) = this%consrv(:,:,:,1)  + RK45_B(isub)*this%QtmpYs

    end subroutine

    subroutine getRHS_Ys(this,rho,u,v,w,rhsYs)
        use operators, only: gradient, curl
        class(solid),                                         intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(out) :: rhsYs

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp1, tmp2, tmp3

        rhsYs = -rho*this%Ys
        tmp1 = rhsYs*u - Ji(:,:,:,1)   
        tmp2 = rhsYs*v - Ji(:,:,:,2)
        tmp3 = rhsYs*w - Ji(:,:,:,3)

        call divergence(this%decomp,this%der,tmp1,tmp2,tmp3,rhsYs)

    end subroutine

    subroutine update_eh(this,isub,dt,rho,u,v,w,tauiiart,divu)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,tauiiart,divu

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhseh  ! RHS for eh equation

        call this%getRHS_eh(rho,u,v,w,tauiiart,divu,rhseh)

        ! advance sub-step
        if(isub==1) this%Qtmpeh = zero                   ! not really needed, since RK45_A(1) = 0
        this%Qtmpeh  = dt*rhseh + RK45_A(isub)*this%Qtmpeh
        this%consrv(:,:,:,2) = this%consrv(:,:,:,2) + RK45_B(isub)*this%Qtmpeh

    end subroutine

    subroutine getRHS_eh(this,rho,u,v,w,tauiiart,divu,rhseh)
        use operators, only: gradient, divergence
        class(solid),                                       intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: divu
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhseh

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp1, tmp2, tmp3

        ! artificial conductivity terms
        call gradient(this%decomp,this%der,-this%eh,tmp1,tmp2,tmp3)
        tmp1 = this%VF * this%kap * tmp1   
        tmp2 = this%VF * this%kap * tmp2   
        tmp3 = this%VF * this%kap * tmp3

        ! add artificial interdiffusional enthalpy terms
        rhseh = -(this%VF*this%eh + this%VF*this%VF*this%p/(this%Ys*rho + epssmall))
        tmp1 = tmp1 + rhseh * this%Ji(:,:,:,1)
        tmp2 = tmp2 + rhseh * this%Ji(:,:,:,2)
        tmp3 = tmp3 + rhseh * this%Ji(:,:,:,3)

        ! add convective terms 
        rhseh = -rho*this%Ys*this%eh
        tmp1 = tmp1 + rhseh*u
        tmp2 = tmp2 + rhseh*v
        tmp3 = tmp3 + rhseh*w

        call divergence(this%decomp,this%der,tmp1,tmp2,tmp3,rhseh)

        ! this would be required for species total energy equation. for
        ! hydrodynamic energy, only divu is required
        !dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        !dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        !dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

        !tmp1 = this%VF*(-this%p + this%sxx); tmp2 = this%VF*(-this%p + this%syy); tmp3 = this%VF*(-this%p + this%szz)
        !rhseh = rhseh  + tmp1*dudx + tmp2*dvdy + tmp3*dwdz &
        !               + this%sxy * (dudy + dvdx) &
        !               + this%sxz * (dudz + dwdx) &
        !               + this%syz * (dvdz + dwdy)

        !nullify(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)

        rhseh = rhseh - this%VF * (this%p + tauiiart) * divu        ! only diagonal terms of artificial stress tensor used here; off-diagonal assumed to affect e_elastic, not e_hydro

    end subroutine

    subroutine update_VF(this,isub,dt,rho,u,v,w)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhsVF  ! RHS for mass fraction equation

        call this%getRHS_VF(rho,u,v,w,rhsVF)

        ! advance sub-step
        if(isub==1) this%QtmpVF = zero                   ! not really needed, since RK45_A(1) = 0
        this%QtmpVF  = dt*rhsVF + RK45_A(isub)*this%QtmpVF
        this%VF = this%VF  + RK45_B(isub)*this%QtmpVF

    end subroutine

    subroutine getRHS_VF(this,u,v,w,rhsVF)
        use operators, only: gradient
        class(solid),                                       intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhsVF

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp1, tmp2, tmp3

        call gradient(this%decomp,this%der,-this%VF,tmp1,tmp2,tmp3)
        rhsVF = u*tmp1 + v*tmp2 + w*tmp3

        ! any penalty term to keep VF between 0 and 1???

    end subroutine

    subroutine filter(this, fil, 1)
        class(solid),  intent(inout) :: this
        type(filters), intent(in)    :: fil
        integer,       intent(in)    :: iflag

        ! filter g
        call this%filter(this%g11, fil, iflag)
        call this%filter(this%g12, fil, iflag)
        call this%filter(this%g13, fil, iflag)
        call this%filter(this%g21, fil, iflag)
        call this%filter(this%g22, fil, iflag)
        call this%filter(this%g23, fil, iflag)
        call this%filter(this%g31, fil, iflag)
        call this%filter(this%g32, fil, iflag)
        call this%filter(this%g33, fil, iflag)

        ! filter Ys
        call this%filter(consrv(:,:,:,1), fil, iflag)

        ! filter eh
        call this%filter(consrv(:,:,:,2), fil, iflag)

        ! filter VF
        call this%filter(this%VF, fil, iflag)

    end subroutine


    ! note ::directly copied from sgrid. Need to go through this subroutine and
    ! make necessary modifications.
    subroutine filter(this,arr,myfil,numtimes)
        class(solid), target, intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(inout) :: arr
        type(filters), target, optional, intent(in) :: myfil
        integer, optional, intent(in) :: numtimes

        type(filters), pointer :: fil2use
        integer :: times2fil
        real(rkind), dimension(:,:,:), pointer :: tmp_in_y, tmp1_in_x, tmp1_in_z, tmp2_in_x, tmp2_in_z
        integer :: lastx, lasty, lastz, idx


        if (present(myfil)) then
            fil2use => myfil
        else
            fil2use => this%fil
        end if

        if (present(numtimes)) then
            times2fil = numtimes
        else
            times2fil = 1
        end if

        ! Allocate pointers for the needed buffers 
        ! Atleast 2 buffers in x and z are assumed
        ! Last two buffers are occupied

        lastx = size(this%xbuf,4)
        lasty = size(this%ybuf,4)
        lastz = size(this%zbuf,4)

        tmp1_in_x => this%xbuf(:,:,:,lastx)
        tmp2_in_x => this%xbuf(:,:,:,lastx-1)
        tmp_in_y => this%ybuf(:,:,:,lasty)
        tmp1_in_z => this%zbuf(:,:,:,lastz)
        tmp2_in_z => this%zbuf(:,:,:,lastz-1)


        ! First filter in y
        call fil2use%filtery(arr,tmp_in_y)
        ! Subsequent refilters 
        do idx = 1,times2fil-1
            arr = tmp_in_y
            call fil2use%filtery(arr,tmp_in_y)
        end do

        ! Then transpose to x
        call transpose_y_to_x(tmp_in_y,tmp1_in_x,this%decomp)

        ! First filter in x
        call fil2use%filterx(tmp1_in_x,tmp2_in_x)
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_x = tmp2_in_x
            call fil2use%filterx(tmp1_in_x,tmp2_in_x)
        end do

        ! Now transpose back to y
        call transpose_x_to_y(tmp2_in_x,tmp_in_y,this%decomp)

        ! Now transpose to z
        call transpose_y_to_z(tmp_in_y,tmp1_in_z,this%decomp)

        !First filter in z
        call fil2use%filterz(tmp1_in_z,tmp2_in_z)
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_z = tmp2_in_z
            call fil2use%filterz(tmp1_in_z,tmp2_in_z)
        end do

        ! Now transpose back to y
        call transpose_z_to_y(tmp2_in_z,arr,this%decomp)

        ! Finished

    end subroutine

    subroutine get_eelastic_devstress(this)
        class(solid), intent(inout) :: this

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6)  :: finger,fingersq
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)    :: trG, trG2, detG

        call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG)
        call this%elastic%get_eelastic(this%rho0,trG,trG2,detG,this%eel)
        call this%elastic%get_devstress(finger, fingersq, trG, trG2, detG, this%devstress)

    end subroutine

    subroutine get_conserved(this,rho)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho

        this%consrv(:,:,:,1) = rho * this%Ys
        this%consrv(:,:,:,2) = this%consrv(:,:,:,1) * this%eh

    end subroutine

    subroutine get_primitive(this,onebyrho)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: onebyrho

        this%Ys = this%consrv(:,:,:,1) * onebyrho
        this%eh = this%consrv(:,:,:,2) / this%consrv(:,:,:,1)

    end subroutine

    !subroutine get_kap(this)
    !    class(solid),                                         intent(in)    :: this
    !    real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: VFeh, tmp1, tmp2, tmp3

    !    ! Step 0: Artificial conductivity based on this quantity
    !    VFeh = this%VF*this%eh

    !    ! Step 1: Get components of grad(e) squared individually
    !    call this%gradient(VFeh,ytmp1,ytmp2,ytmp3) ! Does not use any Y buffers
    !    ytmp1 = ytmp1*ytmp1
    !    ytmp2 = ytmp2*ytmp2
    !    ytmp3 = ytmp3*ytmp3

    !    ! Step 2: Get 4th derivative in X
    !    call transpose_y_to_x(VFeh,xtmp1,this%decomp)
    !    call this%der%d2dx2(xtmp1,xtmp2)
    !    call this%der%d2dx2(xtmp2,xtmp1)
    !    xtmp2 = xtmp1*this%dx**6
    !    call transpose_x_to_y(xtmp2,ytmp4,this%decomp)
    !    kapstar = ytmp4 * ytmp1 / (ytmp1 + ytmp2 + ytmp3)

    !    ! Step 3: Get 4th derivative in Z
    !    call transpose_y_to_z(VFeh,ztmp1,this%decomp)
    !    call this%der%d2dz2(ztmp1,ztmp2)
    !    call this%der%d2dz2(ztmp2,ztmp1)
    !    ztmp2 = ztmp1*this%dz**6
    !    call transpose_z_to_y(ztmp2,ytmp4,this%decomp)
    !    kapstar = kapstar + ytmp4 * ytmp3 / (ytmp1 + ytmp2 + ytmp3)

    !    ! Step 4: Get 4th derivative in Y
    !    call this%der%d2dy2(VFeh,ytmp4)
    !    call this%der%d2dy2(ytmp4,ytmp5)
    !    ytmp4 = ytmp5*this%dy**6
    !    kapstar = kapstar + ytmp4 * ytmp2 / (ytmp1 + ytmp2 + ytmp3)

    !    ! Now, all ytmps are free to use
    !    call this%sgas%get_sos(this%rho,this%p,ytmp1)  ! Speed of sound - hydrodynamic part
    !    call this%elastic%get_sos(this%rho0,ytmp1)     ! Speed of sound - elastic part

    !    kapstar = this%Ckap*this%rho*ytmp1*abs(kapstar)/this%T

    !    ! Filter kapstar
    !    call this%filter(kapstar, this%gfil, 2)

    !    ! Now, add to physical fluid properties
    !    this%mu   = this%mu   + mustar
    !    this%bulk = this%bulk + bulkstar
    !    this%kap  = this%kap  + kapstar

    !end subroutine

    !subroutine get_q(this,duidxj)
    !    class(sgrid), target, intent(inout) :: this
    !    real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(inout) :: duidxj

    !    real(rkind), dimension(:,:,:), pointer :: tmp1_in_x, tmp2_in_x, tmp1_in_y, tmp1_in_z, tmp2_in_z
    !    type(derivatives), pointer :: der

    !    der => this%der

    !    tmp1_in_x => this%xbuf(:,:,:,1)
    !    tmp2_in_x => this%xbuf(:,:,:,2)

    !    tmp1_in_z => this%zbuf(:,:,:,1)
    !    tmp2_in_z => this%zbuf(:,:,:,2)

    !    tmp1_in_y => this%ybuf(:,:,:,1)

    !    ! Step 1: Get qy (dvdy is destroyed)
    !    call der%ddy(this%T,tmp1_in_y)
    !    duidxj(:,:,:,qyidx) = -this%kap*tmp1_in_y

    !    ! Step 2: Get qx (dudx is destroyed)
    !    call transpose_y_to_x(this%T,tmp1_in_x,this%decomp)
    !    call der%ddx(tmp1_in_x,tmp2_in_x)
    !    call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
    !    duidxj(:,:,:,qxidx) = -this%kap*tmp1_in_y

    !    ! Step 3: Get qz (dwdz is destroyed)
    !    call transpose_y_to_z(this%T,tmp1_in_z,this%decomp)
    !    call der%ddz(tmp1_in_z,tmp2_in_z)
    !    call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
    !    duidxj(:,:,:,qzidx) = -this%kap*tmp1_in_y

    !    ! Done
    !end subroutine

end module
