    subroutine init(this, ModelID, spectC, spectE, gpC, gpE, dx, dy, dz, useDynamicProcedure, useClipping, zmesh, nCwall, z0, useWallModel, TfilterZ)
        class(sgs), intent(inout), target :: this
        type(spectral), intent(in), target :: spectC, spectE
        type(decomp_info), intent(in), target :: gpC, gpE
        real(rkind), intent(in) :: dx, dy, dz
        logical, intent(in) :: useDynamicProcedure, useClipping
        integer, intent(in) :: modelID
        integer, intent(in), optional :: nCwall
        real(rkind), dimension(:,:,:), intent(in), optional :: zMesh
        real(rkind), intent(in), optional :: z0
        logical, intent(in), optional :: useWallModel, TfilterZ
        integer :: ierr

        this%SGSmodel = modelID 
        this%useDynamicProcedure = useDynamicProcedure
        this%useClipping = useClipping

        allocate(this%rbuff(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),8))
        this%rbuff = zero  
        this%spectC => spectC
        this%spectE => spectE
        this%deltaFilter = ((1.5*dx)*(1.5*dy)*dz)**(one/three)
        this%deltaTFilter = deltaRatio*this%deltaFilter
        if (present(TfilterZ))  useVerticalTfilter = TfilterZ

        if (present(useWallmodel)) then
            this%useWallModel = useWallModel
            allocate(this%rbuffE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3),3))
            this%rbuffE = zero
            this%nuSGSE => this%rbuffE(:,:,:,3)
            if (this%useWallModel) then
                if (.not. present(z0)) then
                    call GracefulExit("Need to provide z0 if wall model is being used", 442)
                else
                    this%z0 = z0
                end if 
            end if 
        end if 

        select case (this%SGSmodel)
        case(0)
            this%mconst = (this%deltaFilter*c_smag)**2
            call message(1,"SMAGORINSKY SGS model initialized")
            this%useWallFunction = .false. 
        case(1)
            this%mconst = (this%deltaFilter*c_sigma)**2
            allocate(this%SIGMAbuffs(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),14))
            if (this%useDynamicProcedure) then
                call GracefulExit("The standard dynamic procedure is not &
                & available with the SIGMA model.", 321)
            end if
            this%useWallFunction = .false. 
            call message(1,"SIGMA SGS model initialized")
        case(2)
            this%cSMAG_WALL => this%rbuff(:,:,:,8)
            if (this%UseDynamicProcedure) then
                call GracefulExit("Dynamic Procedure cannot be used if damping &
                    & function is being used",3213)
            endif 
            if (.not. present(zMesh)) then
                call GracefulExit("Need to pass in the zCell, z0 values if Dynamic &
                & Procedure is used. If you intend to use the Wall model, use &
                & igridWallM instead of igrid.", 43)
            end if 
            this%cSMAG_WALL = ( c_smag**(-real(ncWall,rkind)) + (kappa*(zMesh/this%deltaFilter + &
                & z0/this%deltaFilter))**(-real(ncWall,rkind))  )**(-one/real(ncWall,rkind))
            this%cSMAG_WALL = (this%deltaFilter*this%cSMAG_WALL)**2    
            this%useWallFunction = .true. 
            call message(1,"SMAGORINSKY (w/ Wall function) SGS model initialized")
        case default 
            call GracefulExit("Invalid choice for SGS model.",2013)
        end select

        this%nuSGS => this%rbuff(:,:,:,7)
        this%nuSGSfil => this%rbuff(:,:,:,8)
        this%gpC => gpC
        this%gpE => gpE
        this%sp_gp => this%spectC%spectdecomp
        this%sp_gpE => this%spectE%spectdecomp
        
        allocate(this%cbuff(this%sp_gp%ysz(1), this%sp_gp%ysz(2), this%sp_gp%ysz(3),2))
        this%nuSGShat => this%cbuff(:,:,:,1)

        allocate(this%ctmpCz(this%sp_gp %zsz(1), this%sp_gp %zsz(2), this%sp_gp %zsz(3)))
        allocate(this%ctmpEz(this%sp_gpE%zsz(1), this%sp_gpE%zsz(2), this%sp_gpE%zsz(3)))
        allocate(this%ctmpEy(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)))
        
        allocate(this%ctmpCz2(this%sp_gp%zsz(1), this%sp_gp%zsz(2), this%sp_gp%zsz(3)))

        if ((useCompactFD) .and. (.not. this%useWallModel)) then
            allocate(this%derZ_EE, this%derZ_OO)
            call this%derZ_EE%init( this%sp_gp%zsz(3), dz, isTopEven = .true., isBotEven = .true., & 
                             isTopSided = .true., isBotSided = .true.) 
            call this%derZ_OO%init( this%sp_gp%zsz(3), dz, isTopEven = .false., isBotEven = .false., & 
                             isTopSided = .true., isBotSided = .true.) 
        else
            allocate(this%Ops2ndOrder)
            call this%Ops2ndOrder%init(gpC,gpE,0,dx,dy,dz,spectC%spectdecomp,spectE%spectdecomp, .true., .true.)
        end if 
        
        allocate(this%Lij(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),6))
        allocate(this%Mij(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),6))
        if (this%useDynamicProcedure) then
            call message(1,"Dynamic Procedure initialized")
        end if 

        this%meanFact = one/(real(gpC%xsz(1))*real(gpC%ysz(2)))
        allocate(this%rtmpY(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
        allocate(this%rtmpZ(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))
        allocate(this%rtmpZ2(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))

        allocate(this%rtmpYE(gpE%ysz(1),gpE%ysz(2),gpE%ysz(3)))
        allocate(this%rtmpZE(gpE%zsz(1),gpE%zsz(2),gpE%zsz(3)))
        allocate(this%rtmpZE2(gpE%zsz(1),gpE%zsz(2),gpE%zsz(3)))

        ierr = this%Tfilz%init(this%sp_gp%zsz(3),.false.)
        ierr = this%Gfilz%init(this%sp_gp%zsz(3),.false.)
    end subroutine

    subroutine destroy(this)
        class(sgs), intent(inout) :: this

        if ((useCompactFD) .and. (.not. this%useWallModel)) then
            call this%derZ_OO%destroy()
            call this%derZ_EE%destroy()
            deallocate(this%derZ_OO, this%derZ_EE)
        else
            call this%Ops2ndOrder%destroy()
            deallocate(this%Ops2ndOrder)
        end if   
         
        if (allocated(this%SIGMAbuffs)) deallocate(this%SIGMAbuffs) 
        nullify( this%nuSGS)
        if(this%useDynamicProcedure) deallocate(this%Lij, this%Mij)
        nullify(this%sp_gp, this%gpC, this%spectC)
        deallocate(this%rbuff, this%cbuff)
        deallocate(this%rtmpZ, this%rtmpY)
    end subroutine
