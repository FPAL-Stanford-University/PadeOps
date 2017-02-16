    subroutine init(this, ModelID, spectC, spectE, gpC, gpE, dx, dy, dz, useDynamicProcedure, useClipping, zmesh,  z0, useWallModel, wallMtype,  TfilterZ, Pr, useWallDamping, nCWall, Cs, CompStokesP , useGlobDynProcedure, GlobDynProcFreq, isStratified, isInviscid, useCoriolis, useWindTurbines, Re)
        class(sgs), intent(inout), target :: this
        type(spectral), intent(in), target :: spectC, spectE
        type(decomp_info), intent(in), target :: gpC, gpE
        real(rkind), intent(in) :: dx, dy, dz
        logical, intent(in) :: useDynamicProcedure, useClipping
        integer, intent(in) :: modelID
        real(rkind), dimension(:,:,:), intent(in), optional :: zMesh
        real(rkind), intent(in), optional :: z0, Pr, Cs, nCWall
        integer, intent(in), optional :: wallMtype, GlobDynProcFreq
        logical, intent(in), optional :: useWallModel, TfilterZ, CompStokesP, useWallDamping, useGlobDynProcedure, isStratified, isInviscid, useCoriolis, useWindTurbines
        real(rkind), intent(in), optional :: Re
        integer :: ierr

        this%SGSmodel = modelID 
        this%useDynamicProcedure = useDynamicProcedure
        this%useClipping = useClipping
        if(present(useGlobDynProcedure)) this%useGlobDynProcedure = useGlobDynProcedure
        if(present(GlobDynProcFreq)) this%GlobDynProcFreq = GlobDynProcFreq

        allocate(this%rbuff(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),11))
        this%rbuff = zero  
        this%spectC => spectC
        this%spectE => spectE
        this%deltaFilter = ((1.5*dx)*(1.5*dy)*dz)**(one/three)
        this%deltaTFilter = deltaRatio*this%deltaFilter
        this%dxsq = dx**2
        this%dysq = dy**2
        this%dzsq = dz**2
        this%dz = dz
        this%nz = gpC%zsz(3)
        if (present(TfilterZ))  useVerticalTfilter = TfilterZ
        if (present(Pr)) this%Pr = Pr
        if (present(CompStokesP)) this%CompStokesP = CompStokesP
        if (present(useWallDamping)) this%useWallFunction = useWallDamping   
        if (present(isStratified)) this%isStratified = isStratified
        if (this%isStratified) this%ntimeAvgQs = 5 ! Update the size of the inst_horz_avg vector from 3 to 5
        if (present(isInviscid)) this%isInviscid = isInviscid
        if (present(useCoriolis)) this%useCoriolis = useCoriolis
        if (present(useWindTurbines)) this%useWindTurbines = useWindTurbines
        if (present(Re)) this%invRe = one/Re

        if (present(useWallmodel)) then
            this%useWallModel = useWallModel
            allocate(this%rbuffE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3),4))
            this%rbuffE = zero
            this%nuSGSE => this%rbuffE(:,:,:,3)
            this%nuSCAE => this%rbuffE(:,:,:,4)
            if (this%useWallModel) then
                if (.not. present(z0)) then
                    call GracefulExit("Need to provide z0 if wall model is being used", 442)
                else
                    this%z0 = z0
                end if
                if (.not. present(wallMtype)) then
                    this%WallModel = 0
                else
                    this%WallModel = wallMtype
                end if  
            end if 
        end if 

        select case (this%SGSmodel)
        case(0)
            this%mconst = (this%deltaFilter*c_smag)**2
            this%eddyViscModel = .true.
            call message(1,"SMAGORINSKY SGS model initialized")
        case(1)
            this%mconst = (this%deltaFilter*c_sigma)**2
            allocate(this%SIGMAbuffs(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),14))
            if (this%useDynamicProcedure) then
                call GracefulExit("The standard dynamic procedure is not &
                & available with the SIGMA model.", 321)
            end if
            this%eddyViscModel = .true.
            call message(1,"SIGMA SGS model initialized")
        case(2)
            this%mconst = (this%deltaFilter*c_smag)**2
            allocate(this%SSI_xbuff1(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
            allocate(this%SSI_xbuff2(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))

            this%eddyViscModel = .true.
            call message(1,"SHEAR-IMPROVED SMAGORINSKY SGS model initialized")
        case(3)
            this%mconst = eight*(this%deltaFilter/c_mgm)**2
            allocate(this%MGMbuffsE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3),6))
            if (this%useDynamicProcedure) then
                call GracefulExit("The standard dynamic procedure is not &
                & available with the MGM model for now.", 321)
            end if
            this%eddyViscModel = .false.
            call message(1,"MGM SGS model initialized")
        case(4)
            this%mconst = one  ! Model constants are absorbed in cx, cy and c
            this%cx = cx_amd*(1.5d0*dx); this%cy = cy_amd*(1.5d0*dy); this%cz = cz_amd*dz
            allocate(this%AMD_Den(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
            if (this%isStratified) then
                allocate(this%dTdzC_diff(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
            end if
            call message(1,"AMD SGS model initialized")
        case default 
            call GracefulExit("Invalid choice for SGS model.",2013)
        end select

        if (this%useWallFunction) then
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
            this%cSMAG_WALL = ( Cs**(-real(ncWall,rkind)) + (kappa*(zMesh/this%deltaFilter + &
                & z0/this%deltaFilter))**(-real(ncWall,rkind))  )**(-one/real(ncWall,rkind))
            this%cSMAG_WALL = (this%deltaFilter*this%cSMAG_WALL)**2    
            call message(1,"Wall damping function for the SGS model initialized")
        end if 

        this%nuSGS => this%rbuff(:,:,:,7)
        this%nuSGSfil => this%rbuff(:,:,:,8)
        this%nuSCA => this%rbuff(:,:,:,11)
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

        if (useCompactFD) then
            allocate(this%derZ_SS, this%derTAU33)
            call this%derTAU33%init( this%sp_gp%zsz(3), dz, isTopEven = .true., isBotEven = .true., & 
                             isTopSided = .false., isBotSided = .true.) 
            call this%derZ_SS%init( this%sp_gp%zsz(3), dz, isTopEven = .false., isBotEven = .false., & 
                             isTopSided = .true., isBotSided = .true.) 
        else
            allocate(this%Ops2ndOrder)
            call this%Ops2ndOrder%init(gpC,gpE,0,dx,dy,dz,spectC%spectdecomp,spectE%spectdecomp, .true., .true.)
        end if 
        allocate(this%OpsNU)
        call this%OpsNU%init(gpC,gpE,0,dx,dy,dz,spectC%spectdecomp,spectE%spectdecomp, .true., .true.)
        
        allocate(this%Lij(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),6))
        allocate(this%Mij(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),6))
        if (this%useDynamicProcedure) then
            allocate(this%uf(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
            allocate(this%vf(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
            allocate(this%wCf(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
            call message(1,"Dynamic Procedure initialized")
        end if 

        if(this%useGlobDynProcedure) then
          !if(this%useCoriolis .or. this%isStratified .or. this%useWindTurbines) then
            allocate(this%GlobDynrbuffxC(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),15))
            allocate(this%GlobDynrbuffxE(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),2))
          !endif

          if(gpC%xst(3)==1) then
            this%zEdgeLo = 1
          else
            this%zEdgeLo = 0
          endif
          !call MPI_COMM_SPLIT(mpi_comm_world, this%zEdgeLo, nrank, this%mycomm, ierr)
          !call MPI_COMM_RANK( this%mycomm, this%myComm_nrank, ierr )
          !call MPI_COMM_SIZE( this%mycomm, this%myComm_nproc, ierr )
        endif

        this%meanFact = one/(real(gpC%xsz(1))*real(gpC%ysz(2)))
        allocate(this%rtmpY(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
        allocate(this%rtmpZ(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))
        allocate(this%rtmpZ2(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))

        allocate(this%rtmpYE(gpE%ysz(1),gpE%ysz(2),gpE%ysz(3)))
        allocate(this%rtmpZE(gpE%zsz(1),gpE%zsz(2),gpE%zsz(3)))
        allocate(this%rtmpZE2(gpE%zsz(1),gpE%zsz(2),gpE%zsz(3)))
        
        if (this%isStratified) then
            allocate(this%q1(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
            allocate(this%q2(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
            allocate(this%q3(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
        end if 

        ierr = this%Tfilz%init(this%sp_gp%zsz(3),.false.)
        ierr = this%Gfilz%init(this%sp_gp%zsz(3),.false.)
    end subroutine

    subroutine destroy(this)
        class(sgs), intent(inout) :: this

        if (useCompactFD) then
            call this%derZ_SS%destroy()
            call this%derTAU33%destroy()
            deallocate(this%derZ_SS, this%derTAU33)
        else
            call this%Ops2ndOrder%destroy()
            deallocate(this%Ops2ndOrder)
        end if   
         
        if (allocated(this%SIGMAbuffs)) deallocate(this%SIGMAbuffs) 
        if (allocated(this%MGMbuffsE)) deallocate(this%MGMbuffsE)
        nullify( this%nuSGS)
        if(this%useDynamicProcedure) deallocate(this%Lij, this%Mij)
        if(this%useGlobDynProcedure)  then
          deallocate(this%GlobDynrbuffxE)
          deallocate(this%GlobDynrbuffxC)
        endif
        nullify(this%sp_gp, this%gpC, this%spectC)
        deallocate(this%rbuff, this%cbuff)
        deallocate(this%rtmpZ, this%rtmpY)
    end subroutine
