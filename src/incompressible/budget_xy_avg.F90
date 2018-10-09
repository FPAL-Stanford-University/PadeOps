module budgets_xy_avg_mod
   use kind_parameters, only: rkind, clen
   use decomp_2d
   use cd06staggstuff, only: cd06stagg
   use reductions, only: p_sum
   use incompressibleGrid, only: igrid  

   implicit none 

   private
   public :: budgets_xy_avg

   ! BUDGET TYPE: 
   ! BUDGET_0: 6 Reynolds stress terms + 3 temp fluxes + meanU + meanV + meanT
   ! BUDGET_1: momentum equation terms (Budget0 also computed) 
   ! BUDGET_2: TKE budget (Budget0 and Budget1 also computed) 


   ! BUDGET_0 term indices:
   ! 1:  <U> 
   ! 2:  <V>
   ! 3:  <T>
   ! 4:  <uu>
   ! 5:  <uv> 
   ! 6:  <uw>
   ! 7:  <vv>
   ! 8:  <vw>
   ! 9:  <ww>
   ! 10: <uT>
   ! 11: <vT>
   ! 12: <wT>
   ! 13: <TT>
   ! 14: <tau13>
   ! 15: <tau23> 
   ! 16: <q3>


   ! BUDGET_1 term indices:  
   ! 1:  X eqn - Advection/convection term
   ! 2:  X eqn - SGS term
   ! 3:  X eqn - Viscous term 
   ! 4:  X eqn - Coriolis Term 
   ! 5:  X eqn - Geostrophic Forcing Term 
   ! 6:  X eqn - Actuator disk/turbine term 
   
   ! 7:  Y eqn - Advection/convection term
   ! 8:  Y eqn - SGS term
   ! 9:  Y eqn - Viscous term 
   ! 10: Y eqn - Coriolis Term 
   ! 11: Y eqn - Geostrophic Forcing term 

   ! 12: Z eqn - Advection term
   ! 13: Z eqn - SGS term 
   ! 14: Z eqn - Pressure gradient term (the mean coriolis, and buoyancy terms are removed)  


   type :: budgets_xy_avg
        private
        integer :: budgetType = 1, run_id, nz

        complex(rkind), dimension(:,:,:), allocatable :: uc, vc, wc, usgs, vsgs, wsgs, uvisc, vvisc, wvisc, px, py, pz, wb, ucor, vcor, wcor, uturb 
        type(igrid), pointer :: igrid_sim 
        
        real(rkind), dimension(:,:), allocatable :: budget_0, budget_1, budget_2
        integer :: counter, rid 
        real(rkind) :: avgFact 
        character(len=clen) :: budgets_dir
        real(rkind), dimension(:), allocatable :: tmp_meanC, tmp_meanE
        real(rkind), dimension(:,:,:), allocatable :: tmpC_1d, tmpE_1d
        real(rkind), dimension(:), allocatable :: buoy_hydrostatic

        type(cd06stagg) :: cd06op_z_mom_budget 

    contains
        procedure           :: init
        procedure           :: updateBudget
        procedure           :: DumpBudget
        procedure           :: destroy
        procedure           :: ResetBudget
        procedure, private  :: restartBudget
        procedure, private  :: AssembleBudget0
        procedure, private  :: AssembleBudget1
        procedure, private  :: get_xy_fluctE_from_fhatE
        procedure, private  :: get_xy_fluctC_from_fhatC
        procedure, private  :: get_xy_meanC_from_fhatC 
        procedure, private  :: get_xy_meanE_from_fhatE 
        procedure, private  :: get_xy_meanC_from_fC 
        procedure, private  :: get_xy_meanE_from_fE 
        procedure, private  :: interp_1d_Edge2Cell 
        !procedure, private  :: AssembleBudget2 (TKE budget is incomplete)
   end type 


contains 

    subroutine init(this, inputfile, igrid_sim) 
        class(budgets_xy_avg), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 
        type(igrid), intent(inout), target :: igrid_sim 
        
        character(len=clen) :: budgets_dir = "NULL"
        integer :: ioUnit, ierr,  budgetType = 1, restart_tid = 0, restart_rid = 0, restart_counter = 0
        logical :: restart_budgets = .false. 

        namelist /BUDGET_XY_AVG/ budgetType, budgets_dir, restart_budgets, restart_rid, restart_tid, restart_counter  
        
        ! STEP 1: Read in inputs, link pointers and allocate budget vectors
        ioUnit = 534
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=BUDGET_XY_AVG)
        close(ioUnit)

        this%igrid_sim => igrid_sim 
        this%run_id = igrid_sim%runid
        this%nz = igrid_sim%nz

        this%budgets_dir = budgets_dir
        this%budgetType = budgetType 
        this%avgFact = 1.d0/(real(igrid_sim%nx,rkind)*real(igrid_sim%ny,rkind))
        allocate(this%Budget_0(this%nz,16))
        allocate(this%Budget_1(this%nz,14))
        if ((trim(budgets_dir) .eq. "null") .or.(trim(budgets_dir) .eq. "NULL")) then
            this%budgets_dir = igrid_sim%outputDir
        end if 

        if (restart_budgets) then
            call this%RestartBudget(restart_rid, restart_tid, restart_counter)
        else
            this%counter = 0
            this%budget_0 = 0.d0 
            this%budget_1 = 0.d0 
        end if

        call this%cd06op_z_mom_budget%init(igrid_sim%nx, igrid_sim%dx, .false., .false., .true., .true.)
        allocate(this%tmp_meanC(this%nz)) 
        allocate(this%tmp_meanE(this%nz+1)) 
        allocate(this%tmpC_1d(1,1,this%nz)) 
        allocate(this%tmpE_1d(1,1,this%nz+1)) 

        allocate(this%buoy_hydrostatic(this%nz)) 


        ! STEP 2: Allocate memory (massive amount of memory needed)
        call igrid_sim%spectC%alloc_r2c_out(this%uc)
        call igrid_sim%spectC%alloc_r2c_out(this%usgs)
        call igrid_sim%spectC%alloc_r2c_out(this%uvisc)
        call igrid_sim%spectC%alloc_r2c_out(this%ucor)
        call igrid_sim%spectC%alloc_r2c_out(this%px)
        call igrid_sim%spectC%alloc_r2c_out(this%uturb)

        call igrid_sim%spectC%alloc_r2c_out(this%vc)
        call igrid_sim%spectC%alloc_r2c_out(this%vsgs)
        call igrid_sim%spectC%alloc_r2c_out(this%vvisc)
        call igrid_sim%spectC%alloc_r2c_out(this%vcor)
        call igrid_sim%spectC%alloc_r2c_out(this%py)

        call igrid_sim%spectE%alloc_r2c_out(this%wc)
        call igrid_sim%spectE%alloc_r2c_out(this%wsgs)
        call igrid_sim%spectE%alloc_r2c_out(this%wvisc)
        call igrid_sim%spectE%alloc_r2c_out(this%wcor)
        call igrid_sim%spectE%alloc_r2c_out(this%pz)
        call igrid_sim%spectE%alloc_r2c_out(this%wb)


        ! STEP 3: Now instrument igrid 
        call igrid_sim%instrumentForBudgets(this%uc, this%vc, this%wc, this%usgs, this%vsgs, this%wsgs, &
                   & this%uvisc, this%vvisc, this%wvisc, this%px, this%py, this%pz, this%wb, this%ucor, &
                   & this%vcor, this%wcor, this%uturb)  

    end subroutine 

    subroutine ResetBudget(this)
        class(budgets_xy_avg), intent(inout) :: this
        
        this%counter = 0
        this%budget_0 = 0.d0 
        this%budget_1 = 0.d0 

    end subroutine 
    
    subroutine destroy(this)
        class(budgets_xy_avg), intent(inout) :: this

        nullify(this%igrid_sim)
        deallocate(this%uc, this%vc, this%wc, this%usgs, this%vsgs, this%wsgs, &
                   & this%uvisc, this%vvisc, this%wvisc, this%px, this%py, this%pz, this%wb, this%ucor, &
                   & this%vcor, this%wcor, this%uturb)  
        deallocate(this%budget_0, this%budget_1)

    end subroutine 


    subroutine updateBudget(this)
        class(budgets_xy_avg), intent(inout) :: this
   
        call this%igrid_sim%getMomentumTerms()  

        select case (this%budgetType)
        case(0)
            call this%AssembleBudget0()
        case(1)
            call this%AssembleBudget0()
            call this%AssembleBudget1()
        case(2)
            call this%AssembleBudget0()
            call this%AssembleBudget1()
            !call this%AssembleBudget2()
        end select

        this%counter = this%counter + 1

    end subroutine 

    subroutine DumpBudget(this)
        class(budgets_xy_avg), intent(inout) :: this
        character(len=clen) :: fname, tempname 

        if (nrank == 0) then
            ! Budget 0: 
            write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",this%rid,"_budget0","_t",this%igrid_sim%step,"_n",this%counter,".stt"
            fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
            
            this%budget_0 = this%budget_0/(real(this%counter,rkind) + 1.d-18)
            this%budget_0(:,4)  = this%budget_0(:,4)  - this%budget_0(:,1)*this%budget_0(:,1) ! R11
            this%budget_0(:,5)  = this%budget_0(:,5)  - this%budget_0(:,1)*this%budget_0(:,2) ! R12
            this%budget_0(:,7)  = this%budget_0(:,7)  - this%budget_0(:,2)*this%budget_0(:,2) ! R22
            this%budget_0(:,10) = this%budget_0(:,10) - this%budget_0(:,1)*this%budget_0(:,3) ! <uT>
            this%budget_0(:,11) = this%budget_0(:,11) - this%budget_0(:,2)*this%budget_0(:,3) ! <vT>
            this%budget_0(:,13) = this%budget_0(:,13) - this%budget_0(:,3)*this%budget_0(:,3) ! <TT>

            call write_2d_ascii(this%budget_0, fname) 
            
            this%budget_0(:,4)  = this%budget_0(:,4)  + this%budget_0(:,1)*this%budget_0(:,1) ! R11
            this%budget_0(:,5)  = this%budget_0(:,5)  + this%budget_0(:,1)*this%budget_0(:,2) ! R12
            this%budget_0(:,7)  = this%budget_0(:,7)  + this%budget_0(:,2)*this%budget_0(:,2) ! R22
            this%budget_0(:,10) = this%budget_0(:,10) + this%budget_0(:,1)*this%budget_0(:,3) ! <uT>
            this%budget_0(:,11) = this%budget_0(:,11) + this%budget_0(:,2)*this%budget_0(:,3) ! <vT>
            this%budget_0(:,13) = this%budget_0(:,13) + this%budget_0(:,3)*this%budget_0(:,3) ! <TT>
            this%budget_0 = this%budget_0*(real(this%counter,rkind) + 1.d-18)

            ! Budget 1: 
            write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",this%rid,"_budget1","_t",this%igrid_sim%step,"_n",this%counter,".stt"
            fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
            this%budget_1 = this%budget_1/(real(this%counter,rkind) + 1.d-18)
            call write_2d_ascii(this%budget_1, fname) 
            this%budget_1 = this%budget_1*(real(this%counter,rkind) + 1.d-18)
        end if 
    end subroutine 


    subroutine restartBudget(this, rid, tid, cid)
        class(budgets_xy_avg), intent(inout) :: this
        integer, intent(in) :: rid, cid, tid
        character(len=clen) :: fname, tempname 

        ! Budget 0: 
        write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",rid,"_budget0","_t",tid,"_n",cid,".stt"
        fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
        if (allocated(this%budget_0)) deallocate(this%budget_0)
        call read_2d_ascii(this%budget_0, fname)
        this%budget_0 = this%budget_0*(real(cid,rkind) + 1.d-18)

        ! Budget 1: 
        write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",rid,"_budget1","_t",tid,"_n",cid,".stt"
        fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
        if (allocated(this%budget_1)) deallocate(this%budget_1)
        call read_2d_ascii(this%budget_1, fname)
        this%budget_1 = this%budget_1*(real(cid,rkind) + 1.d-18)

        this%counter = cid
    end subroutine 


    subroutine get_xy_meanC_from_fhatC(this, fhat, fmean)
        class(budgets_xy_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%sp_gpC%ysz(1), this%igrid_sim%sp_gpC%ysz(2), this%igrid_sim%sp_gpC%ysz(3)), intent(inout) :: fhat
        real(rkind), dimension(this%nz), intent(out) :: fmean

        call transpose_y_to_z(fhat, this%igrid_sim%cbuffzC(:,:,:,1), this%igrid_sim%sp_gpC)
        if (nrank == 0) then
            fmean = this%igrid_sim%cbuffzC(1,1,:,1)*this%avgFact
        else
            fmean = 0.d0 ! Only 0 processor has the actual mean  
        end if 

    end subroutine 

    subroutine get_xy_meanE_from_fhatE(this, fhat, fmean)
        class(budgets_xy_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%sp_gpE%ysz(1), this%igrid_sim%sp_gpE%ysz(2), this%igrid_sim%sp_gpE%ysz(3)), intent(inout) :: fhat
        real(rkind), dimension(this%nz+1), intent(out) :: fmean

        call transpose_y_to_z(fhat, this%igrid_sim%cbuffzE(:,:,:,1), this%igrid_sim%sp_gpE)
        if (nrank == 0) then
            fmean = this%igrid_sim%cbuffzE(1,1,:,1)*this%avgFact
        else
            fmean = 0.d0 ! Only 0 processor has the actual mean  
        end if 

    end subroutine 
    
    subroutine get_xy_meanC_from_fC(this, f, fmean)
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1), this%igrid_sim%gpC%xsz(2), this%igrid_sim%gpC%xsz(3)), intent(inout) :: f
        real(rkind), dimension(this%nz), intent(out) :: fmean

        call this%igrid_sim%spectC%fft(f, this%igrid_sim%cbuffyC(:,:,:,1))
        call transpose_y_to_z(this%igrid_sim%cbuffyC(:,:,:,1), this%igrid_sim%cbuffzC(:,:,:,1), this%igrid_sim%sp_gpC)
        if (nrank == 0) then
            fmean = this%igrid_sim%cbuffzC(1,1,:,1)*this%avgFact
        else
            fmean = 0.d0 ! Only 0 processor has the actual mean  
        end if 

    end subroutine 
    
    subroutine get_xy_meanE_from_fE(this, f, fmean)
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpE%xsz(1), this%igrid_sim%gpE%xsz(2), this%igrid_sim%gpE%xsz(3)), intent(inout) :: f
        real(rkind), dimension(this%nz+1), intent(out) :: fmean

        call this%igrid_sim%spectE%fft(f, this%igrid_sim%cbuffyE(:,:,:,1))
        call transpose_y_to_z(this%igrid_sim%cbuffyE(:,:,:,1), this%igrid_sim%cbuffzE(:,:,:,1), this%igrid_sim%sp_gpE)
        if (nrank == 0) then
            fmean = this%igrid_sim%cbuffzE(1,1,:,1)*this%avgFact
        else
            fmean = 0.d0 ! Only 0 processor has the actual mean  
        end if 

    end subroutine 
    
    subroutine get_xy_fluctC_from_fhatC(this, fhat, ffluct)
        class(budgets_xy_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%sp_gpC%ysz(1), this%igrid_sim%sp_gpC%ysz(2), this%igrid_sim%sp_gpC%ysz(3)), intent(inout) :: fhat
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(out) :: ffluct 

        this%igrid_sim%cbuffyC(:,:,:,1) = fhat
        if (this%igrid_sim%spectC%carryingZeroK) then
            this%igrid_sim%cbuffyC(1,1,:,1) = 0.d0 
        end if 
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1),ffluct)

    end subroutine 

    subroutine get_xy_fluctE_from_fhatE(this, fhat, ffluct)
        class(budgets_xy_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%sp_gpE%ysz(1), this%igrid_sim%sp_gpE%ysz(2), this%igrid_sim%sp_gpE%ysz(3)), intent(inout) :: fhat
        real(rkind), dimension(this%igrid_sim%gpE%xsz(1),this%igrid_sim%gpE%xsz(2),this%igrid_sim%gpE%xsz(3)), intent(out) :: ffluct 

        this%igrid_sim%cbuffyE(:,:,:,1) = fhat
        if (this%igrid_sim%spectE%carryingZeroK) then
            this%igrid_sim%cbuffyE(1,1,:,1) = 0.d0 
        end if 
        call this%igrid_sim%spectE%ifft(this%igrid_sim%cbuffyE(:,:,:,1),ffluct)

    end subroutine


    subroutine interp_1d_Edge2Cell(this, fE,fC)
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%nz+1), intent(in)  :: fE
        real(rkind), dimension(this%nz  ), intent(out) :: fC


        this%tmpE_1d(1,1,:) = fE
        call this%cd06op_z_mom_budget%interpz_E2C(this%tmpE_1d, this%tmpC_1d, 1, 1)
        fC = this%tmpC_1d(1,1,:)

    end subroutine 


    subroutine AssembleBudget0(this)
        class(budgets_xy_avg), intent(inout) :: this

        ! STEP 1: Compute mean U, V and T
        call this%get_xy_meanC_from_fhatC(this%igrid_sim%uhat, this%tmp_meanC)
        this%budget_0(:,1) = this%budget_0(:,1) + this%tmp_meanC

        call this%get_xy_meanC_from_fhatC(this%igrid_sim%vhat, this%tmp_meanC)
        this%budget_0(:,2) = this%budget_0(:,2) + this%tmp_meanC
        
        call this%get_xy_meanC_from_fhatC(this%igrid_sim%That, this%tmp_meanC)
        this%budget_0(:,3) = this%budget_0(:,3) + this%tmp_meanC

        ! STEP 2: Get Reynolds stresses (IMPORTANT: need to correct for fluctuation before dumping)
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%u*this%igrid_sim%u
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
        this%budget_0(:,4) = this%budget_0(:,4) + this%tmp_meanC

        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%u*this%igrid_sim%v
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
        this%budget_0(:,5) = this%budget_0(:,5) + this%tmp_meanC
        
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%uE*this%igrid_sim%w
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_0(:,6) = this%budget_0(:,6) + this%tmp_meanC
        
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%v*this%igrid_sim%v
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
        this%budget_0(:,7) = this%budget_0(:,7) + this%tmp_meanC

        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%vE*this%igrid_sim%w
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_0(:,8) = this%budget_0(:,8) + this%tmp_meanC
        
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%w
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_0(:,9) = this%budget_0(:,9) + this%tmp_meanC
        
        ! STEP 2: Get Temperature fluxes and variances (IMPORTANT: need to correct for fluctuation before dumping)
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%u*this%igrid_sim%T
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
        this%budget_0(:,10) = this%budget_0(:,10) + this%tmp_meanC
        
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%v*this%igrid_sim%T
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
        this%budget_0(:,11) = this%budget_0(:,11) + this%tmp_meanC

        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%TE*this%igrid_sim%w
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_0(:,12) = this%budget_0(:,12) + this%tmp_meanC

        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%T*this%igrid_sim%T
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
        this%budget_0(:,13) = this%budget_0(:,13) + this%tmp_meanC

        ! STEP 3: SGS stress (also viscous stress if finite reynolds number is being used)
        call this%get_xy_meanE_from_fE(this%igrid_sim%tau13, this%tmp_meanE) 
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_0(:,14) = this%budget_0(:,14) + this%tmp_meanC

        call this%get_xy_meanE_from_fE(this%igrid_sim%tau23, this%tmp_meanE) 
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_0(:,15) = this%budget_0(:,15) + this%tmp_meanC

        if (associated(this%igrid_Sim%q3)) then
            call this%get_xy_meanE_from_fE(this%igrid_sim%q3, this%tmp_meanE) 
            call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
            this%budget_0(:,16) = this%budget_0(:,16) + this%tmp_meanC
        end if 

    end subroutine 


    subroutine AssembleBudget1(this)
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind) :: tmp1, tmp2

        ! Get the geostrophic forcing 
        call this%igrid_sim%get_geostrophic_forcing(tmp1, tmp2)         ! Forcing in x and y directions respectively 

        ! Get u- momentum budget
        call this%get_xy_meanC_from_fhatC(this%uc, this%tmp_meanC)
        this%budget_1(:,1) = this%budget_1(:,1) + this%tmp_meanC

        call this%get_xy_meanC_from_fhatC(this%usgs, this%tmp_meanC)
        this%budget_1(:,2) = this%budget_1(:,2) + this%tmp_meanC
        
        call this%get_xy_meanC_from_fhatC(this%uvisc, this%tmp_meanC)
        this%budget_1(:,3) = this%budget_1(:,3) + this%tmp_meanC

        call this%get_xy_meanC_from_fhatC(this%ucor, this%tmp_meanC)
        this%budget_1(:,4) = this%budget_1(:,4) + this%tmp_meanC - tmp1 ! Remove the geostrophic forcing term

        this%budget_1(:,5) = this%budget_1(:,5) + tmp1

        call this%get_xy_meanC_from_fhatC(this%uturb, this%tmp_meanC)
        this%budget_1(:,6) = this%budget_1(:,6) + this%tmp_meanC
       
        ! Get y- momentum budget 
        call this%get_xy_meanC_from_fhatC(this%vc, this%tmp_meanC)
        this%budget_1(:,7) = this%budget_1(:,7) + this%tmp_meanC

        call this%get_xy_meanC_from_fhatC(this%vsgs, this%tmp_meanC)
        this%budget_1(:,8) = this%budget_1(:,8) + this%tmp_meanC

        call this%get_xy_meanC_from_fhatC(this%vvisc, this%tmp_meanC)
        this%budget_1(:,9) = this%budget_1(:,9) + this%tmp_meanC

        call this%get_xy_meanC_from_fhatC(this%vcor, this%tmp_meanC)
        this%budget_1(:,10) = this%budget_1(:,10) + this%tmp_meanC - tmp2 ! Remove the geostrophic forcing term

        this%budget_1(:,11) = this%budget_1(:,11) + tmp2

        ! Get z- momentum equation 
        call this%get_xy_meanE_from_fhatE(this%wc, this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_1(:,12) = this%budget_1(:,12) + this%tmp_meanC

        call this%get_xy_meanE_from_fhatE(this%wsgs, this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_1(:,13) = this%budget_1(:,13) + this%tmp_meanC

        call this%get_xy_meanE_from_fhatE(this%pz, this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_1(:,14) = this%budget_1(:,14) + this%tmp_meanC

    end subroutine 
end module 
