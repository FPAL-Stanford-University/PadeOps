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


   ! BUDGET_1 term indices:  
   ! 1:  X eqn - Advection/convection term
   ! 2:  X eqn - SGS term
   ! 3:  X eqn - Viscous term 
   ! 4:  X eqn - Coriolis Term 
   ! 5:  X eqn - Geostrophic Forcing Term 
   ! 7:  X eqn - Actuator disk/turbine term 
   
   ! 7:  Y eqn - Advection/convection term
   ! 8:  Y eqn - SGS term
   ! 9:  Y eqn - Viscous term 
   ! 10: Y eqn - Coriolis Term 
   ! 11: Y eqn - Geostrophic Forcing term 

   ! 12: Z eqn - Advection term 
   ! 13: Z eqn - Pressure gradient term 
   ! 14: Z eqn - Coriolis term 
   ! 15: Z eqn - Buoyancy term 


   type :: budgets_xy_avg
        private
        integer :: budgetType = 1, run_id, nz

        complex(rkind), dimension(:,:,:), allocatable :: uc, vc, wc, usgs, vsgs, wsgs, uvisc, vvisc, wvisc, px, py, pz, wb, ucor, vcor, wcor, uturb 
        type(igrid), pointer :: igrid_sim 
        
        real(rkind), dimension(:,:), allocatable :: budget_0, budget_1, budget_2
        integer :: counter, tid, rid 
        real(rkind) :: avgFact 
        character(len=clen) :: budgets_dir

        type(cd06stagg) :: cd06op_z_mom_budget 

    contains
        procedure           :: init
        procedure           :: updateBudget
        procedure           :: DumpBudget
        procedure           :: destroy
        procedure, private  :: restartBudget
        procedure           :: ResetBudget
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
        allocate(this%Budget_1(this%nz,15))
        allocate(this%Budget_0(this%nz,12))
        if ((trim(budgets_dir) .eq. "null") .or.(trim(budgets_dir) .eq. "NULL")) then
            this%budgets_dir = igrid_sim%outputDir
        end if 

        if (restart_budgets) then
            this%counter = restart_counter
            call this%RestartBudget(restart_rid, restart_tid, restart_counter)
        else
            this%counter = 0
            this%budget_0 = 0.d0 
            this%budget_1 = 0.d0 
        end if

        call this%cd06op_z_mom_budget%init(igrid_sim%nx, igrid_sim%dx, .false., .false., .true., .true.)


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
    end subroutine 


    subroutine DumpBudget(this)
        class(budgets_xy_avg), intent(inout) :: this
    end subroutine 


    subroutine restartBudget(this, rid, tid, cid)
        class(budgets_xy_avg), intent(inout) :: this
        integer, intent(in) :: rid, cid, tid

    end subroutine 


end module 
