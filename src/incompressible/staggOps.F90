module staggOpsMod
    use kind_parameters, only: rkind
    use exits, only: GracefulExit
    use decomp_2d
    use constants, only: half,one,two,three, zero

    implicit none

    private

    public :: staggOps

    type :: staggOps
        private
        integer :: nxC, nyC, nzC
        integer :: nxE, nyE, nzE
        integer :: nxC_cmplx, nyC_cmplx, nzC_cmplx
        integer :: nxE_cmplx, nyE_cmplx, nzE_cmplx
        type(decomp_info), pointer :: cellDecomp
        type(decomp_info), pointer :: edgeDecomp
        integer :: stagg_scheme
        real(rkind) :: dx, dy, dz
        logical :: isBotSided = .false. 
        logical :: isTopSided = .false. 
        contains
            procedure :: init
            procedure :: destroy
            procedure, private :: InterpZ_Edge2Cell_CMPLX
            procedure, private :: InterpZ_Edge2Cell_REAL
            procedure, private :: InterpZ_Cell2Edge_CMPLX
            procedure, private :: InterpZ_Cell2Edge_REAL
            procedure, private :: ddz_E2C_cmplx
            procedure, private :: ddz_E2C_real
            procedure, private :: ddz_C2E_cmplx
            procedure, private :: ddz_C2E_real
            procedure, private :: ddz_C2C_cmplx
            procedure, private :: ddz_C2C_real
            procedure, private :: d2dz2_C2C_real 
            procedure, private :: d2dz2_E2E_real
            procedure, private :: d2dz2_C2C_cmplx 
            procedure, private :: d2dz2_E2E_cmplx
            generic :: InterpZ_Cell2Edge => InterpZ_Cell2Edge_CMPLX, InterpZ_Cell2Edge_REAL
            generic :: InterpZ_Edge2Cell => InterpZ_Edge2Cell_CMPLX, InterpZ_Edge2Cell_REAL
            generic :: ddz_C2E => ddz_C2E_real, ddz_C2E_cmplx
            generic :: ddz_E2C => ddz_E2C_real, ddz_E2C_cmplx
            generic :: ddz_C2C => ddz_C2C_real, ddz_C2C_cmplx
            generic :: d2dz2_E2E => d2dz2_E2E_cmplx, d2dz2_E2E_real
            generic :: d2dz2_C2C => d2dz2_C2C_cmplx, d2dz2_C2C_real
    end type

contains

    pure subroutine ddz_E2C_Real(this,fE,dfdzC)
        class(staggOps), intent(in) :: this
        real(rkind), dimension(this%nxE,this%nyE,this%nzE), intent(in) :: fE
        real(rkind), dimension(this%nxC,this%nyC,this%nzC), intent(out):: dfdzC
        real(rkind) :: OneByDz

        OneByDz = one/this%dz
        dfdzC(:,:,1:this%nzC) = fE(:,:,2:this%nzE) - fE(:,:,1:this%nzE-1)
        dfdzC = OneByDz*dfdzC

    end subroutine

    pure subroutine ddz_E2C_Cmplx(this,fE,dfdzC)
        class(staggOps), intent(in) :: this
        complex(rkind), dimension(this%nxE_cmplx,this%nyE_cmplx,this%nzE_cmplx), intent(in) :: fE
        complex(rkind), dimension(this%nxC_cmplx,this%nyC_cmplx,this%nzC_cmplx), intent(out):: dfdzC
        real(rkind) :: OneByDz

        OneByDz = one/this%dz
        dfdzC(:,:,1:this%nzC) = fE(:,:,2:this%nzE) - fE(:,:,1:this%nzE-1)
        dfdzC = OneByDz*dfdzC

    end subroutine

    pure subroutine ddz_C2E_Real(this,fC,dfdzE, isTopEven, isBotEven)
        class(staggOps), intent(in) :: this
        real(rkind), dimension(this%nxC,this%nyC,this%nzC), intent(in) :: fC
        real(rkind), dimension(this%nxE,this%nyE,this%nzE), intent(out):: dfdzE
        real(rkind) :: OneByDz
        logical, intent(in) :: isTopEven, isBotEven

        OneByDz = one/(this%dz)
        dfdzE(:,:,2:this%nzE-1) =  fC(:,:,2:this%nzC) - fC(:,:,1:this%nzC-1) 

        if (isBotEven) then
            dfdzE(:,:,1) = zero
        else 
            dfdzE(:,:,1) = two*fc(:,:,1)
        end if

        if (isTopEven) then
            dfdzE(:,:,this%nzE) = zero
        else 
            dfdzE(:,:,this%nzE) = -two*fC(:,:,this%nzC)
        end if

        dfdzE = dfdzE*OneByDz

    end subroutine

    pure subroutine ddz_C2E_Cmplx(this,fC,dfdzE, isTopEven, isBotEven)
        class(staggOps), intent(in) :: this
        complex(rkind), dimension(this%nxC_cmplx,this%nyC_cmplx,this%nzC_cmplx), intent(in) :: fC
        complex(rkind), dimension(this%nxE_cmplx,this%nyE_cmplx,this%nzE_cmplx), intent(out):: dfdzE
        real(rkind) :: OneByDz
        logical, intent(in) :: isTopEven, isBotEven

        OneByDz = one/(this%dz)
        dfdzE(:,:,2:this%nzE-1) =  fC(:,:,2:this%nzC) - fC(:,:,1:this%nzC-1) 

        if (isBotEven) then
            dfdzE(:,:,1) = zero
        else 
            dfdzE(:,:,1) = two*fc(:,:,1)
        end if

        if (isTopEven) then
            dfdzE(:,:,this%nzE) = zero
        else 
            dfdzE(:,:,this%nzE) = -two*fC(:,:,this%nzC)
        end if

        dfdzE = dfdzE*OneByDz

    end subroutine

    pure subroutine ddz_C2C_Real(this,fC,dfdzC, isTopEven, isBotEven)
        class(staggOps), intent(in) :: this
        real(rkind), dimension(this%nxC,this%nyC,this%nzC), intent(in) :: fC
        real(rkind), dimension(this%nxC,this%nyC,this%nzC), intent(out):: dfdzC
        real(rkind) :: OneBy2Dz
        logical, intent(in) :: isTopEven, isBotEven

        OneBy2Dz = one/(two*this%dz)
        dfdzC(:,:,2:this%nzC-1) =  fC(:,:,3:this%nzC) - fC(:,:,1:this%nzC-2) 
        dfdzC(:,:,2:this%nzC-1) = OneBy2Dz*dfdzC(:,:,2:this%nzC-1)

        if (isBotEven) then
            dfdzC(:,:,1) = OneBy2Dz*(fC(:,:,2) - fC(:,:,1))
        else 
            dfdzC(:,:,1) = OneBy2Dz*(fC(:,:,2) + fC(:,:,1))
        end if

        if (isTopEven) then
            dfdzC(:,:,this%nzC) = OneBy2Dz*(fC(:,:,this%nzC) - fC(:,:,this%nzC-1))
        else 
            dfdzC(:,:,this%nzC) = -OneBy2Dz*(fC(:,:,this%nzC) + fC(:,:,this%nzC-1))
        end if

    end subroutine
  
    pure subroutine ddz_C2C_cmplx(this,fC,dfdzC, isTopEven, isBotEven)
        class(staggOps), intent(in) :: this
        complex(rkind), dimension(this%nxC_cmplx,this%nyC_cmplx,this%nzC_cmplx), intent(in) :: fC
        complex(rkind), dimension(this%nxC_cmplx,this%nyC_cmplx,this%nzC_cmplx), intent(out):: dfdzC
        real(rkind) :: OneBy2Dz
        logical, intent(in) :: isTopEven, isBotEven

        OneBy2Dz = one/(two*this%dz)
        dfdzC(:,:,2:this%nzC-1) =  fC(:,:,3:this%nzC) - fC(:,:,1:this%nzC-2) 
        dfdzC(:,:,2:this%nzC-1) = OneBy2Dz*dfdzC(:,:,2:this%nzC-1)

        if (this%isBotSided) then
            !dfdzC(:,:,1) = (-half*fC(:,:,3) + two*fC(:,:,2) - three/two*fC(:,:,1))/this%dz 
            dfdzC(:,:,1) = (fC(:,:,2) - fC(:,:,1))/this%dz 
        else
            if (isBotEven) then
                dfdzC(:,:,1) = OneBy2Dz*(fC(:,:,2) - fC(:,:,1))
            else 
                dfdzC(:,:,1) = OneBy2Dz*(fC(:,:,2) + fC(:,:,1))
            end if
        end if 

        if (this%isTopSided) then
            !dfdzC(:,:,this%nzC) = (half*fC(:,:,this%nzC-2) - two*fC(:,:,this%nzC-1) &
            !            + three/two*fC(:,:,this%nzC))/this%dz 
            dfdzC(:,:,this%nzC) = (fC(:,:,this%nzC) - fC(:,:,this%nzC-1))/this%dz 
        else
            if (isTopEven) then
                dfdzC(:,:,this%nzC) = OneBy2Dz*(fC(:,:,this%nzC) - fC(:,:,this%nzC-1))
            else 
                dfdzC(:,:,this%nzC) = -OneBy2Dz*(fC(:,:,this%nzC) + fC(:,:,this%nzC-1))
            end if
        end if 

    end subroutine


    pure subroutine d2dz2_C2C_cmplx(this,fC,d2fdz2C, isTopEven, isBotEven)
        class(staggOps), intent(in) :: this
        complex(rkind), dimension(this%nxC_cmplx,this%nyC_cmplx,this%nzC_cmplx), intent(in) :: fC
        complex(rkind), dimension(this%nxC_cmplx,this%nyC_cmplx,this%nzC_cmplx), intent(out):: d2fdz2C
        real(rkind) :: OneByDzSq
        logical, intent(in) :: isTopEven, isBotEven

        OneByDzSq = one/(this%dz**2)
        d2fdz2C(:,:,2:this%nzC-1) =  fC(:,:,3:this%nzE) - two*fC(:,:,2:this%nzC-1) + fC(:,:,1:this%nzE-2) 
        
        if (isBotEven) then
            d2fdz2C(:,:,1) = fC(:,:,2) - fC(:,:,1)
        else 
            d2fdz2C(:,:,1) = fC(:,:,2) - three*fC(:,:,1)
        end if

        if (isTopEven) then
            d2fdz2C(:,:,this%nzC) = fC(:,:,this%nzC-1) - fC(:,:,this%nzC)
        else
            d2fdz2C(:,:,this%nzC) = fC(:,:,this%nzC-1) - three*fC(:,:,this%nzC)
        end if 

        d2fdz2C = OneByDzSq*d2fdz2C

    end subroutine



    pure subroutine d2dz2_C2C_real(this,fC,d2fdz2C, isTopEven, isBotEven)
        class(staggOps), intent(in) :: this
        real(rkind), dimension(this%nxC,this%nyC,this%nzC), intent(in) :: fC
        real(rkind), dimension(this%nxC,this%nyC,this%nzC), intent(out):: d2fdz2C
        real(rkind) :: OneByDzSq
        logical, intent(in) :: isTopEven, isBotEven

        OneByDzSq = one/(this%dz**2)
        d2fdz2C(:,:,2:this%nzC-1) =  fC(:,:,3:this%nzE) - two*fC(:,:,2:this%nzC-1) + fC(:,:,1:this%nzE-2) 
        
        if (isBotEven) then
            d2fdz2C(:,:,1) = fC(:,:,2) - fC(:,:,1)
        else 
            d2fdz2C(:,:,1) = fC(:,:,2) - three*fC(:,:,1)
        end if

        if (isTopEven) then
            d2fdz2C(:,:,this%nzC) = fC(:,:,this%nzC-1) - fC(:,:,this%nzC)
        else
            d2fdz2C(:,:,this%nzC) = fC(:,:,this%nzC-1) - three*fC(:,:,this%nzC)
        end if 

        d2fdz2C = OneByDzSq*d2fdz2C

    end subroutine

    pure subroutine d2dz2_E2E_cmplx(this,fE,d2fdz2E, isTopEven, isBotEven)
        class(staggOps), intent(in) :: this
        complex(rkind), dimension(this%nxE_cmplx,this%nyE_cmplx,this%nzE_cmplx), intent(in) :: fE
        complex(rkind), dimension(this%nxE_cmplx,this%nyE_cmplx,this%nzE_cmplx), intent(out):: d2fdz2E
        real(rkind) :: OneByDzSq
        logical, intent(in) :: isTopEven, isBotEven

        OneByDzSq = one/(this%dz**2)
        d2fdz2E(:,:,2:this%nzE-1) =  fE(:,:,3:this%nzE) - two*fE(:,:,2:this%nzE-1) + fE(:,:,1:this%nzE-2) 
        
        if (isBotEven) then
            d2fdz2E(:,:,1) = two*(fE(:,:,2) - fE(:,:,1)) 
        else 
            d2fdz2E(:,:,1) = zero  
        end if

        if (isTopEven) then
            d2fdz2E(:,:,this%nzE) = two*(fE(:,:,this%nzE-1) - fE(:,:,this%nzE))
        else
            d2fdz2E(:,:,this%nzE) = zero  
        end if 

        d2fdz2E = OneByDzSq*d2fdz2E

    end subroutine

    pure subroutine d2dz2_E2E_real(this,fE,d2fdz2E, isTopEven, isBotEven)
        class(staggOps), intent(in) :: this
        real(rkind), dimension(this%nxE,this%nyE,this%nzE), intent(in) :: fE
        real(rkind), dimension(this%nxE,this%nyE,this%nzE), intent(out):: d2fdz2E
        real(rkind) :: OneByDzSq
        logical, intent(in) :: isTopEven, isBotEven

        OneByDzSq = one/(this%dz**2)
        d2fdz2E(:,:,2:this%nzE-1) =  fE(:,:,3:this%nzE) - two*fE(:,:,2:this%nzE-1) + fE(:,:,1:this%nzE-2) 
        
        if (isBotEven) then
            d2fdz2E(:,:,1) = two*(fE(:,:,2) - fE(:,:,1)) 
        else 
            d2fdz2E(:,:,1) = zero  
        end if

        if (isTopEven) then
            d2fdz2E(:,:,this%nzE) = two*(fE(:,:,this%nzE-1) - fE(:,:,this%nzE))
        else
            d2fdz2E(:,:,this%nzE) = zero  
        end if 

        d2fdz2E = OneByDzSq*d2fdz2E

    end subroutine


    subroutine init(this, gpC, gpE, stagg_scheme , dx, dy, dz, gpCspect, gpEspect, isTopSided, isBotSided)
        class(staggOps), intent(inout) :: this
        class(decomp_info), intent(in), target:: gpC, gpE
        integer, intent(in) :: stagg_scheme
        real(rkind), intent(in) :: dx, dy, dz
        class(decomp_info), intent(in), optional, target:: gpCspect, gpEspect
        logical, intent(in), optional :: isTopSided, isBotSided

        if ((present(isTopSided)) .and. (present(isBotSided))) then
            this%isTopSided = isTopSided; this%isBotSided = isBotSided
        end if 
        this%nxC = gpC%zsz(1)
        this%nyC = gpC%zsz(2)
        this%nzC = gpC%zsz(3)

        this%nxE = gpE%zsz(1) 
        this%nyE = gpE%zsz(2) 
        this%nzE = gpE%zsz(3)
        
        this%nxC_cmplx = 0 
        this%nyC_cmplx = 0 
        this%nzC_cmplx = 0 

        this%nxE_cmplx = 0 
        this%nyE_cmplx = 0 
        this%nzE_cmplx = 0 

        if (present(gpCspect)) then
            this%nxC_cmplx = gpCspect%zsz(1)
            this%nyC_cmplx = gpCspect%zsz(2)
            this%nzC_cmplx = gpCspect%zsz(3)
        end if 

        if (present(gpEspect)) then
            this%nxE_cmplx = gpEspect%zsz(1) 
            this%nyE_cmplx = gpEspect%zsz(2) 
            this%nzE_cmplx = gpEspect%zsz(3)
        end if 

        this%stagg_scheme = stagg_scheme
        this%cellDecomp => gpC
        this%edgeDecomp => gpE

        this%dx = dx
        this%dy = dy
        this%dz = dz

    end subroutine


    pure subroutine destroy(this)
        class(staggOps), intent(inout) :: this

        nullify(this%edgeDecomp)
        nullify(this%cellDecomp)

    end subroutine

    pure subroutine InterpZ_Edge2Cell_REAL(this, edgeArr, cellArr)
        class(staggOps), intent(in) :: this
        real(rkind), intent(in), dimension(this%nxE, this%nyE, this%nzE) :: edgeArr
        real(rkind), intent(out), dimension(this%nxC, this%nyC, this%nzC) :: cellArr

        cellArr = edgeArr(1:this%nxC,1:this%nyC,1:this%nzC)
        cellArr = cellArr + edgeArr(:,:,2:this%nzC+1)
        cellArr = half*cellArr   

    end subroutine
    
    pure subroutine InterpZ_Edge2Cell_CMPLX(this, edgeArr, cellArr)
        class(staggOps), intent(in) :: this
        complex(rkind), intent(in), dimension(this%nxE_cmplx, this%nyE_cmplx, this%nzE_cmplx) :: edgeArr
        complex(rkind), intent(out), dimension(this%nxC_cmplx, this%nyC_cmplx, this%nzC_cmplx) :: cellArr

        cellArr = edgeArr(1:this%nxC,1:this%nyC,1:this%nzC)
        cellArr = cellArr + edgeArr(:,:,2:this%nzC+1)
        cellArr = half*cellArr   

    end subroutine 

    pure subroutine InterpZ_Cell2Edge_CMPLX(this, cellArr, edgeArr, BotVal, TopVal)
        class(staggOps), intent(in) :: this
        complex(rkind), intent(in), dimension(this%nxC_cmplx, this%nyC_cmplx, this%nzC_cmplx) :: cellArr
        complex(rkind), intent(out), dimension(this%nxE_cmplx, this%nyE_cmplx, this%nzE_cmplx) :: edgeArr
        complex(rkind), intent(in) :: BotVal, TopVal

        edgeArr(:,:,this%nzE) = TopVal
        edgeArr(:,:,1          ) = BotVal
        edgeArr(:,:,2:this%nzE-1) = cellArr(:,:,1:this%nzC-1) 
        edgeArr(:,:,2:this%nzE-1) = edgeArr(:,:,2:this%nzE-1) + cellArr(:,:,2:this%nzC)
        edgeArr = half*edgeArr
    end subroutine

    pure subroutine InterpZ_Cell2Edge_REAL(this, cellArr, edgeArr, BotVal, TopVal)
        class(staggOps), intent(in) :: this
        real(rkind), intent(in), dimension(this%nxC, this%nyC, this%nzC) :: cellArr
        real(rkind), intent(out), dimension(this%nxE, this%nyE, this%nzE) :: edgeArr
        real(rkind), intent(in) :: BotVal, TopVal

        edgeArr(:,:,this%nzE) = TopVal
        edgeArr(:,:,1          ) = BotVal
        edgeArr(:,:,2:this%nzE-1) = cellArr(:,:,1:this%nzC-1) 
        edgeArr(:,:,2:this%nzE-1) = edgeArr(:,:,2:this%nzE-1) + cellArr(:,:, 2:this%nzC)
        edgeArr(:,:,2:this%nzE-1) = half*edgeArr(:,:,2:this%nzE)
    end subroutine

end module 
