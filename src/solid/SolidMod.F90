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

    contains

        procedure :: getRHS_g
        procedure :: getPlasticSources
        procedure :: getRHS_VF
        procedure :: getRHS_eh
        procedure :: update_g
        procedure :: update_VF
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

        ! Now deallocate the EOS objects
        if ( allocated(this%hydro)      ) deallocate(this%hydro)
        if ( allocated(this%elastic)    ) deallocate(this%elastic)

    end subroutine

    subroutine getRHS_g(this,rho,u,v,w,rhsg)
        use operators, only: gradient, curl
        class(solid),                                         intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg

        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg
        real(rkind) :: etafac = one/32._rkind

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

end module
