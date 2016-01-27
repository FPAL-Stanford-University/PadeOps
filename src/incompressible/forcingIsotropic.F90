module forcingMod
    use kind_parameters, only: rkind
    use decomp_2d, only: decomp_info, decomp_info_init, &
                    transpose_x_to_y, transpose_y_to_x, &
                    transpose_y_to_z, transpose_z_to_y 
    use spectralMod, only: spectral 
    use reductions, only: P_SUM
    use exits, only: message 
    use constants, only: zero

    implicit none

    private
    public :: forcing

    real(rkind) :: dissConst = 0.1_rkind
    type :: forcing
        private
        real(rkind) :: dissRate, KEforcing
        integer, dimension(:), allocatable :: shell_i, shell_j, shell_k
        integer :: nForceShellsG,nForceShellsL
        real(rkind) :: multFact
        real(rkind) :: nu
        contains
            procedure :: init
            procedure :: addForcing
            procedure :: destroy
    end type  

contains

    subroutine init(this,spect,kfmax,nxg,nyg,nzg,nu)
        class(forcing), intent(inout) :: this
        real(rkind), intent(in) :: kfmax
        class(spectral), intent(in), target :: spect
        real(rkind) :: kfmax_sq
        real(rkind), dimension(:,:,:), pointer :: kabssq
        integer :: counter, i, j, k 
        integer, intent(in) :: nxg, nyg, nzg
        real(rkind), intent(in) :: nu

        counter = 0
        kfmax_sq = kfmax**2
        kabssq => spect%kabs_sq
        do k = lbound(kabssq,3),ubound(kabssq,3)
            do j = lbound(kabssq,2),ubound(kabssq,2)
                do i = lbound(kabssq,1),ubound(kabssq,1)
                    if (kabssq(i,j,k) .le. kfmax_sq) then
                        counter = counter + 1
                    end if 
                end do 
            end do 
        end do 

        allocate(this%shell_i(counter),this%shell_j(counter),this%shell_k(counter))
        counter = 0
        do k = lbound(kabssq,3),ubound(kabssq,3)
            do j = lbound(kabssq,2),ubound(kabssq,2)
                do i = lbound(kabssq,1),ubound(kabssq,1)
                    if (kabssq(i,j,k) .le. kfmax_sq) then
                        counter = counter + 1
                        this%shell_i(counter) = i
                        this%shell_j(counter) = j
                        this%shell_k(counter) = k
                    end if 
                end do 
            end do 
        end do 

        this%nu = nu
        this%nForceShellsL = counter
        this%nForceShellsG = P_SUM(this%nForceShellsL)
        this%multFact = (real(nxg,rkind)*real(nyg,rkind)*real(nzg,rkind))**2
        call message(0,"Initiazed the FORCING derived type successfully")
        call message(1,"Total number of forced shells (global)", this%nForceShellsG)
    end subroutine

    subroutine addForcing(this,Sfield,duidxj,rhs)
        class(forcing), intent(in) :: this
        complex(rkind), dimension(:,:,:,:), intent(in), target :: Sfield
        complex(rkind), dimension(:,:,:,:), intent(inout)      :: rhs
        real(rkind) :: Elocal, Eglobal
        integer :: shellid, i, j, k
        real(rkind), dimension(:,:,:,:), intent(in), target :: duidxj
        complex(rkind), dimension(:,:,:), pointer :: uhat, vhat, what
        real(rkind) :: mconst
        real(rkind), dimension(:,:,:), pointer :: dudy, dudz, dvdx, dvdz, dwdx, dwdy
        real(rkind) :: enstrophy, dissipation

        uhat => Sfield(:,:,:,1)
        vhat => Sfield(:,:,:,1)
        what => Sfield(:,:,:,1)
        dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3)
        dvdx => duidxj(:,:,:,4); dvdz => duidxj(:,:,:,6)
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8)
        
        Elocal = zero
        do shellid = 1,this%nForceShellsL
            i = this%shell_i(shellid) 
            j = this%shell_j(shellid) 
            k = this%shell_k(shellid) 

            Elocal = Elocal + real(uhat(i,j,k),rkind)**2 + aimag(uhat(i,j,k))**2
            Elocal = Elocal + real(vhat(i,j,k),rkind)**2 + aimag(vhat(i,j,k))**2
            Elocal = Elocal + real(what(i,j,k),rkind)**2 + aimag(what(i,j,k))**2
        end do 
        Eglobal = P_SUM(Elocal)
        
        enstrophy = zero
        do k = lbound(dudy,3),ubound(dudy,3)
            do j = lbound(dudy,2),ubound(dudy,2)
                do i = lbound(dudy,1),ubound(dudy,1)
                    enstrophy = enstrophy + (dwdy(i,j,k) - dvdz(i,j,k))**2 + & 
                        (dudz(i,j,k) - dwdx(i,j,k))**2 + (dvdx(i,j,k) - dudy(i,j,k))**2
                end do 
            end do 
        end do 
        enstrophy = enstrophy/sqrt(this%multFact)
        dissipation = this%nu*P_SUM(enstrophy)
        mconst = dissConst*dissipation*this%multFact/Eglobal
        call message(0,"Mult Const:", mconst)
       
        do shellid = 1,this%nForceShellsL
            i = this%shell_i(shellid) 
            j = this%shell_j(shellid) 
            k = this%shell_k(shellid) 
            
            rhs(i,j,k,1) = rhs(i,j,k,1) + mconst*uhat(i,j,k)
            rhs(i,j,k,2) = rhs(i,j,k,2) + mconst*vhat(i,j,k)
            rhs(i,j,k,3) = rhs(i,j,k,3) + mconst*what(i,j,k)
        end do

    end subroutine

    subroutine destroy(this)
        class(forcing), intent(inout) :: this
    
        if (allocated(this%shell_i)) deallocate(this%shell_i)
        if (allocated(this%shell_j)) deallocate(this%shell_j)
        if (allocated(this%shell_k)) deallocate(this%shell_k)

    end subroutine
end module 
