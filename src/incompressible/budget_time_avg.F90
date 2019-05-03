module budgets_time_avg_mod
   use kind_parameters, only: rkind, clen, mpirkind
   use decomp_2d
   use reductions, only: p_sum
   use incompressibleGrid, only: igrid  
   use exits, only: message, GracefulExit
   use basic_io, only: read_2d_ascii, write_2d_ascii
   use constants, only: half
   use mpi 

   implicit none 

   private
   public :: budgets_time_avg

   ! BUDGET TYPE: 
   ! BUDGET_0: 6 Reynolds stress terms + 3 temp fluxes + meanU + meanV + meanT
   ! BUDGET_1: momentum equation terms (Budget0 also computed) 
   ! BUDGET_2: MKE budget (Budget0 and Budget1 also computed) 
   ! BUDGET 3: TKE budget (Budget 0, 1 and 2 also included)


   ! BUDGET_0 term indices:
   ! 1:  <U> 
   ! 2:  <V>
   ! 3:  <W>
   ! 4:  <uu>
   ! 5:  <uv> 
   ! 6:  <uw>
   ! 7:  <vv>
   ! 8:  <vw>
   ! 9:  <ww>
   ! 10: <P>
   ! 11: <tau11> 
   ! 12: <tau12> 
   ! 13: <tau13>
   ! 14: <tau22> 
   ! 15: <tau23> 
   ! 16: <tau33> 


   ! BUDGET_1 term indices:  
   ! 1:  X eqn - Advection/convection term
   ! 2:  X eqn - Pressure gradient term
   ! 3:  X eqn - SGS term
   ! 4:  X eqn - Actuator disk/turbine term 
   
   ! 5:  Y eqn - Advection/convection term
   ! 6:  Y eqn - Pressure gradient term
   ! 7:  Y eqn - SGS term

   ! 8:  Z eqn - Advection term
   ! 9:  Z eqn - Pressure gradient term 
   ! 10: Z eqn - SGS term 


   ! BUDGET_2 term indices: 
   ! 1:  Loss to Resolved TKE      (G)
   ! 2:  Advective transport       (B)
   ! 3:  Reynolds stress transport (E)
   ! 4:  Pressure transport        (C)
   ! 5:  SGS + viscous transport   (D+F)
   ! 6:  Loss to SGS TKE + viscous dissipation (H+I)
   ! 7:  Actuator disk sink        (J)


   ! BUDGET_3 term indices:
   ! 1. TKE production              (G)
   ! 2. convective transport        (B)
   ! 3. turbulent transport         (C)
   ! 4. Pressure transport          (D)
   ! 5. SGS + viscous transport     (E+F)
   ! 6. SGS + viscous dissipation   (H+I)
   ! 7. Actuator disk/Turbine sink  (J)


   ! BUDGET_4_ij term indices: 
   ! <incomplete/not-needed for now>

   type :: budgets_time_avg
        private
        integer :: budgetType = 1, run_id, nz

        complex(rkind), dimension(:,:,:), allocatable :: uc, vc, wc, usgs, vsgs, wsgs, px, py, pz, uturb 
        type(igrid), pointer :: igrid_sim 
        
        real(rkind), dimension(:,:,:,:), allocatable :: budget_0, budget_1, budget_2, budget_3
        integer :: counter
        character(len=clen) :: budgets_dir

        integer :: tidx_dump 
        integer :: tidx_compute
        integer :: tidx_budget_start 
        logical :: do_budgets

    contains
        procedure           :: init
        procedure           :: destroy
        procedure           :: ResetBudget
        procedure           :: DoBudgets
        
        procedure, private  :: updateBudget
        procedure, private  :: DumpBudget
        procedure, private  :: restartBudget
        procedure, private  :: dump_budget_field 
        
        procedure, private  :: AssembleBudget0
        procedure, private  :: DumpBudget0 
        
        procedure, private  :: AssembleBudget1
        procedure, private  :: DumpBudget1
        
        procedure, private  :: AssembleBudget2
        procedure, private  :: DumpBudget2
        
        procedure, private  :: AssembleBudget3
        procedure, private  :: DumpBudget3
        
        procedure, private :: ddx_R2R
        procedure, private :: ddy_R2R
        procedure, private :: ddz_R2R
        procedure, private :: ddx_C2R
        procedure, private :: ddy_C2R
        procedure, private :: ddz_C2R
        procedure, private :: interp_Edge2Cell
    end type 


contains 

    subroutine init(this, inputfile, igrid_sim) 
        class(budgets_time_avg), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 
        type(igrid), intent(inout), target :: igrid_sim 
        
        character(len=clen) :: budgets_dir = "NULL"
        integer :: ioUnit, ierr,  budgetType = 1, restart_tid = 0, restart_rid = 0, restart_counter = 0
        logical :: restart_budgets = .false. 
        integer :: tidx_compute = 1000000, tidx_dump = 1000000, tidx_budget_start = 0
        logical :: do_budgets = .false. 
        namelist /BUDGET_TIME_AVG/ budgetType, budgets_dir, restart_budgets, restart_rid, restart_tid, restart_counter, tidx_dump, tidx_compute, do_budgets, tidx_budget_start 
        
        ! STEP 1: Read in inputs, link pointers and allocate budget vectors
        ioUnit = 534
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=BUDGET_TIME_AVG)
        close(ioUnit)

        this%igrid_sim => igrid_sim 
        this%run_id = igrid_sim%runid
        this%nz = igrid_sim%nz
        this%do_budgets = do_budgets
        this%tidx_dump = tidx_dump
        this%tidx_compute = tidx_compute
        this%tidx_budget_start = tidx_budget_start  

        this%budgets_dir = budgets_dir
        this%budgetType = budgetType 
       
        if(this%do_budgets) then 
            allocate(this%budget_0(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3),16))
            allocate(this%budget_1(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3),10))
            allocate(this%budget_2(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3),07))
            allocate(this%budget_3(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3),08))

            if ((trim(budgets_dir) .eq. "null") .or.(trim(budgets_dir) .eq. "NULL")) then
                this%budgets_dir = igrid_sim%outputDir
            end if 

            if (restart_budgets) then
                call GracefulExit("To be done",1234)
                call this%RestartBudget(restart_rid, restart_tid, restart_counter)
            else
                call this%resetBudget()
            end if

            ! STEP 2: Allocate memory (massive amount of memory needed)
            call igrid_sim%spectC%alloc_r2c_out(this%uc)
            call igrid_sim%spectC%alloc_r2c_out(this%usgs)
            call igrid_sim%spectC%alloc_r2c_out(this%px)
            call igrid_sim%spectC%alloc_r2c_out(this%uturb)

            call igrid_sim%spectC%alloc_r2c_out(this%vc)
            call igrid_sim%spectC%alloc_r2c_out(this%vsgs)
            call igrid_sim%spectC%alloc_r2c_out(this%py)

            call igrid_sim%spectE%alloc_r2c_out(this%wc)
            call igrid_sim%spectE%alloc_r2c_out(this%wsgs)
            call igrid_sim%spectE%alloc_r2c_out(this%pz)


            ! STEP 3: Now instrument igrid 
            call igrid_sim%instrumentForBudgets_TimeAvg(this%uc, this%vc, this%wc, this%usgs, this%vsgs, this%wsgs, &
                       & this%px, this%py, this%pz, this%uturb)  

        end if

    end subroutine 


    subroutine doBudgets(this)
        class(budgets_time_avg), intent(inout) :: this

        if (this%do_budgets .and. (this%igrid_sim%step>this%tidx_budget_start)) then
        
            if (mod(this%igrid_sim%step,this%tidx_compute) .eq. 0) then
                call this%updateBudget()
            end if

            if (mod(this%igrid_sim%step,this%tidx_dump) .eq. 0) then
                call this%dumpBudget()
                call message(0,"Dumped a budget .stt file")
            end if 

        end if 
    end subroutine 

    subroutine updateBudget(this)
        class(budgets_time_avg), intent(inout) :: this
   
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
            ! Budget 2 need not be assembled now; it only needs to be assembled
            ! before writing to disk 
        case(3)
            call this%AssembleBudget0()
            call this%AssembleBudget1()
            call this%AssembleBudget3()
        end select

        this%counter = this%counter + 1

    end subroutine 


    subroutine DumpBudget(this)
        class(budgets_time_avg), intent(inout) :: this
        
        ! MKE budget is only assembled before dumping
        if (this%budgetType>1) call this%AssembleBudget2() 
       
        ! Budget 0: 
        call this%dumpbudget0()

        ! Budget 1: 
        if (this%budgetType>0) then
            call this%dumpbudget1()
        end if 
        
        ! Budget 2: 
        if (this%budgetType>1) then
            call this%dumpbudget2()
        end if 

        ! Budget 3: 
        if (this%budgetType>2) then
            call this%dumpbudget3()
        end if 
    end subroutine 

    ! ---------------------- Budget 0 ------------------------
    subroutine DumpBudget0(this)
        class(budgets_time_avg), intent(inout) :: this
        integer :: idx 
        
        ! Step 1: Get the average from sum
        this%budget_0 = this%budget_0/(real(this%counter,rkind) + 1.d-18)
        
        ! Step 2: Get the <Rij> from <ui uj>
        this%budget_0(:,:,:,4)  = this%budget_0(:,:,:,4)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
        this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
        this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
        this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
        this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
        this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  - this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
        
        ! Step 3: Dump the full budget 
        do idx = 1,size(this%budget_0,4)
            call this%dump_budget_field(this%budget_0(:,:,:,idx),idx,0)
        end do 
        
        ! Step 4: Go back to <ui uj> from <Rij>
        this%budget_0(:,:,:,4)  = this%budget_0(:,:,:,4)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
        this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
        this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
        this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
        this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
        this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
        
        ! Step 5: Go back to summing instead of averaging
        this%budget_0 = this%budget_0*(real(this%counter,rkind) + 1.d-18)

    end subroutine 

    subroutine AssembleBudget0(this)
        class(budgets_time_avg), intent(inout) :: this

        ! STEP 1: Compute mean U, V and W
        this%budget_0(:,:,:,1) = this%budget_0(:,:,:,1) + this%igrid_sim%u
        this%budget_0(:,:,:,2) = this%budget_0(:,:,:,2) + this%igrid_sim%v
        this%budget_0(:,:,:,3) = this%budget_0(:,:,:,3) + this%igrid_sim%wC

        ! STEP 2: Get Reynolds stresses (IMPORTANT: need to correct for fluctuation before dumping)
        this%budget_0(:,:,:,4) = this%budget_0(:,:,:,4) + this%igrid_sim%u*this%igrid_sim%u
        this%budget_0(:,:,:,5) = this%budget_0(:,:,:,5) + this%igrid_sim%u*this%igrid_sim%v
        this%budget_0(:,:,:,6) = this%budget_0(:,:,:,6) + this%igrid_sim%u*this%igrid_sim%w
        this%budget_0(:,:,:,7) = this%budget_0(:,:,:,7) + this%igrid_sim%v*this%igrid_sim%v
        this%budget_0(:,:,:,8) = this%budget_0(:,:,:,8) + this%igrid_sim%v*this%igrid_sim%w
        this%budget_0(:,:,:,9) = this%budget_0(:,:,:,9) + this%igrid_sim%w*this%igrid_sim%w

        ! STEP 3: Pressure
        this%budget_0(:,:,:,10) = this%budget_0(:,:,:,10) + this%igrid_sim%pressure

        ! STEP 3: SGS stresses (also viscous stress if finite reynolds number is being used)
        call this%igrid_sim%sgsmodel%populate_tauij_E_to_C()
        this%budget_0(:,:,:,11:16) = this%budget_0(:,:,:,11:16) + this%igrid_sim%tauSGS_ij 

    end subroutine 

    ! ---------------------- Budget 1 ------------------------
    subroutine AssembleBudget1(this)
        class(budgets_time_avg), intent(inout) :: this

        ! STEP 1: Get 4 terms from u-equation 
        call this%igrid_sim%spectC%ifft(this%uc,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_1(:,:,:,1) = this%budget_1(:,:,:,1) + this%igrid_sim%rbuffxC(:,:,:,1)

        call this%igrid_sim%spectC%ifft(this%px,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_1(:,:,:,2) = this%budget_1(:,:,:,2) + this%igrid_sim%rbuffxC(:,:,:,1)
        
        call this%igrid_sim%spectC%ifft(this%usgs,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_1(:,:,:,3) = this%budget_1(:,:,:,3) + this%igrid_sim%rbuffxC(:,:,:,1)
        
        call this%igrid_sim%spectC%ifft(this%uturb,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_1(:,:,:,4) = this%budget_1(:,:,:,4) + this%igrid_sim%rbuffxC(:,:,:,1)

        ! STEP 2: Get 3 terms from v-equation 
        call this%igrid_sim%spectC%ifft(this%vc,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_1(:,:,:,5) = this%budget_1(:,:,:,5) + this%igrid_sim%rbuffxC(:,:,:,1)

        call this%igrid_sim%spectC%ifft(this%py,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_1(:,:,:,6) = this%budget_1(:,:,:,6) + this%igrid_sim%rbuffxC(:,:,:,1)
        
        call this%igrid_sim%spectC%ifft(this%vsgs,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_1(:,:,:,7) = this%budget_1(:,:,:,7) + this%igrid_sim%rbuffxC(:,:,:,1)
        
        ! STEP 2: Get 3 terms from w-equation 
        call this%igrid_sim%spectE%ifft(this%wc,this%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%igrid_sim%rbuffxE(:,:,:,1), this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_1(:,:,:,8) = this%budget_1(:,:,:,8) + this%igrid_sim%rbuffxC(:,:,:,1)

        call this%igrid_sim%spectE%ifft(this%pz,this%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%igrid_sim%rbuffxE(:,:,:,1), this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_1(:,:,:,9) = this%budget_1(:,:,:,9) + this%igrid_sim%rbuffxC(:,:,:,1)
        
        call this%igrid_sim%spectE%ifft(this%wsgs,this%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%igrid_sim%rbuffxE(:,:,:,1), this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_1(:,:,:,10) = this%budget_1(:,:,:,10) + this%igrid_sim%rbuffxC(:,:,:,1)
        
        
    end subroutine

    subroutine DumpBudget1(this)
        class(budgets_time_avg), intent(inout) :: this
        integer :: idx 

        ! Step 1: Get the average from sum
        this%budget_1 = this%budget_1/(real(this%counter,rkind) + 1.d-18)
        
        ! Step 2: Dump the full budget 
        do idx = 1,size(this%budget_1,4)
            call this%dump_budget_field(this%budget_1(:,:,:,idx),idx,1)
        end do 
        
        ! Step 3: Go back to summing instead of averaging
        this%budget_1 = this%budget_1*(real(this%counter,rkind) + 1.d-18)
    end subroutine 

    ! ---------------------- Budget 2 ------------------------
    subroutine AssembleBudget2(this)
        class(budgets_time_avg), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer :: Umn, Vmn, Wmn, R11, R12, R13, R22, R23, R33
        real(rkind), dimension(:,:,:), pointer :: Pmn, tau11, tau12, tau13, tau22, tau23, tau33
        real(rkind), dimension(:,:,:), pointer :: buff, buff2

        if (this%counter > 0) then
            ! < Incomplete: Look at budget_xy_avg for reference. > 

            ! Get the average from sum
            this%budget_0 = this%budget_0/(real(this%counter,rkind) + 1.d-18)
            this%budget_0(:,:,:,4)  = this%budget_0(:,:,:,4)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
            this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
            this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
            this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
            this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
            this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  - this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
            this%budget_1 = this%budget_1/(real(this%counter,rkind) + 1.d-18)

            Umn => this%budget_0(:,:,:,1);    Vmn => this%budget_0(:,:,:,2);      Wmn => this%budget_0(:,:,:,3);
            R11 => this%budget_0(:,:,:,4);    R12 => this%budget_0(:,:,:,5);      R13 => this%budget_0(:,:,:,6)
            R22 => this%budget_0(:,:,:,7);    R23 => this%budget_0(:,:,:,8);      R33 => this%budget_0(:,:,:,9)
            Pmn => this%budget_0(:,:,:,10);   tau11 => this%budget_0(:,:,:,11);   tau12 => this%budget_0(:,:,:,12)
            tau13 => this%budget_0(:,:,:,13); tau22 => this%budget_0(:,:,:,14);   tau23 => this%budget_0(:,:,:,15)
            tau33 => this%budget_0(:,:,:,16); buff => this%igrid_sim%rbuffxC(:,:,:,1); buff2 => this%igrid_sim%rbuffxC(:,:,:,2)

            ! 1:  Loss to Resolved TKE  (G)  && 6 : Loss to SGS + viscous dissipation (H+I)
            call this%ddx_R2R(Umn, buff); 
            this%budget_2(:,:,:,1) = R11*buff;
            this%budget_2(:,:,:,6) = tau11*buff
            
            call this%ddx_R2R(Vmn, buff);
            this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + R12*buff
            this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + tau12*buff

            call this%ddx_R2R(Wmn, buff);
            this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + R13*buff
            this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + tau13*buff

            call this%ddy_R2R(Umn, buff);
            this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + R12*buff
            this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + tau12*buff

            call this%ddy_R2R(Vmn, buff);
            this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + R22*buff
            this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + tau22*buff

            call this%ddy_R2R(Wmn, buff);
            this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + R23*buff
            this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + tau23*buff

            call this%ddz_R2R(Umn, buff);
            this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + R13*buff
            this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + tau13*buff

            call this%ddz_R2R(Vmn, buff);
            this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + R23*buff
            this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + tau23*buff

            call this%ddz_R2R(Wmn, buff);
            this%budget_2(:,:,:,1) = this%budget_2(:,:,:,1) + R33*buff
            this%budget_2(:,:,:,6) = this%budget_2(:,:,:,6) + tau33*buff

            ! 2:  Advective transport       (B)
            buff2 = half*(Umn*Umn + Vmn*Vmn + Wmn*Wmn)
            call this%ddx_R2R(buff2,buff); this%budget_2(:,:,:,2) = -Umn*buff
            call this%ddy_R2R(buff2,buff); this%budget_2(:,:,:,2) = this%budget_2(:,:,:,2) - Vmn*buff
            call this%ddz_R2R(buff2,buff); this%budget_2(:,:,:,2) = this%budget_2(:,:,:,2) - Wmn*buff

            ! 3:  Reynolds stress transport (E)
            this%budget_2(:,:,:,3) = Umn*this%budget_1(:,:,:,1) + Vmn*this%budget_1(:,:,:,5) + Wmn*this%budget_1(:,:,:,8)
            this%budget_2(:,:,:,3) = this%budget_2(:,:,:,3) - this%budget_2(:,:,:,1) - this%budget_2(:,:,:,2)

            ! 4:  Pressure transport        (C)
            this%budget_2(:,:,:,4) = Umn*this%budget_1(:,:,:,2) + Vmn*this%budget_1(:,:,:,6) + Wmn*this%budget_1(:,:,:,9)

            ! 5:  SGS + viscous transport   (D+F)
            this%budget_2(:,:,:,5) = Umn*this%budget_1(:,:,:,3) + Vmn*this%budget_1(:,:,:,7) + Wmn*this%budget_1(:,:,:,10)
            this%budget_2(:,:,:,5) = this%budget_2(:,:,:,5) - this%budget_2(:,:,:,6)

            ! 7:  Actuator disk sink        (J)
            this%budget_2(:,:,:,7) = Umn*this%budget_1(:,:,:,4)


            nullify(Umn,Vmn,Wmn,R11,R12,R13,R22,R23,R33,Pmn,tau11,tau12,tau13,tau22,tau23,tau33,buff,buff2)

            ! Go back to sum
            this%budget_1 = this%budget_1*(real(this%counter,rkind) + 1.d-18)
            this%budget_0(:,:,:,4)  = this%budget_0(:,:,:,4)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
            this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
            this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
            this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
            this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
            this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
            this%budget_0 = this%budget_0*(real(this%counter,rkind) + 1.d-18)

        end if 

    end subroutine 
    
    subroutine DumpBudget2(this)
        class(budgets_time_avg), intent(inout) :: this
        integer :: idx

        ! Dump the full budget 
        do idx = 1,size(this%budget_2,4)
            call this%dump_budget_field(this%budget_2(:,:,:,idx),idx,2)
        end do 

    end subroutine 

    
    ! ---------------------- Budget 3 ------------------------
    subroutine AssembleBudget3(this)
        class(budgets_time_avg), intent(inout) :: this

        call this%igrid_sim%sgsmodel%populate_tauij_E_to_C()

        ! 3. turbulent transport         (C)
        call this%igrid_sim%spectC%ifft(this%uc,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,3) = this%budget_3(:,:,:,3) + this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,1)

        call this%igrid_sim%spectC%ifft(this%vc,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,3) = this%budget_3(:,:,:,3) + this%igrid_sim%v*this%igrid_sim%rbuffxC(:,:,:,1)

        call this%igrid_sim%spectE%ifft(this%wc,this%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%igrid_sim%rbuffxE(:,:,:,1), this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,3) = this%budget_3(:,:,:,3) + this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,1)


        ! 4. Pressure transport          (D)
        call this%igrid_sim%spectC%ifft(this%px,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,1)
        
        call this%igrid_sim%spectC%ifft(this%py,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + this%igrid_sim%v*this%igrid_sim%rbuffxC(:,:,:,1)

        call this%igrid_sim%spectE%ifft(this%pz,this%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%igrid_sim%rbuffxE(:,:,:,1), this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,1)


        ! 5. SGS + viscous transport     (E+F)
        call this%igrid_sim%spectC%ifft(this%usgs,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,1)
        
        call this%igrid_sim%spectC%ifft(this%vsgs,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + this%igrid_sim%v*this%igrid_sim%rbuffxC(:,:,:,1)

        call this%igrid_sim%spectE%ifft(this%wsgs,this%igrid_sim%rbuffxE(:,:,:,1))
        call this%interp_Edge2Cell(this%igrid_sim%rbuffxE(:,:,:,1), this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,1)
        

        ! 6. SGS + viscous dissipation   (H+I)
        call this%ddx_R2R(this%igrid_sim%u, this%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%igrid_sim%rbuffxC(:,:,:,1)*this%igrid_sim%tauSGS_ij(:,:,:,1)

        call this%ddx_R2R(this%igrid_sim%v, this%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%igrid_sim%rbuffxC(:,:,:,1)*this%igrid_sim%tauSGS_ij(:,:,:,2)

        call this%ddx_R2R(this%igrid_sim%wC, this%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%igrid_sim%rbuffxC(:,:,:,1)*this%igrid_sim%tauSGS_ij(:,:,:,3)

        call this%ddy_R2R(this%igrid_sim%u, this%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%igrid_sim%rbuffxC(:,:,:,1)*this%igrid_sim%tauSGS_ij(:,:,:,2)

        call this%ddy_R2R(this%igrid_sim%v, this%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%igrid_sim%rbuffxC(:,:,:,1)*this%igrid_sim%tauSGS_ij(:,:,:,4)

        call this%ddy_R2R(this%igrid_sim%wC, this%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%igrid_sim%rbuffxC(:,:,:,1)*this%igrid_sim%tauSGS_ij(:,:,:,5)

        call this%ddz_R2R(this%igrid_sim%u, this%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%igrid_sim%rbuffxC(:,:,:,1)*this%igrid_sim%tauSGS_ij(:,:,:,3)

        call this%ddz_R2R(this%igrid_sim%v, this%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%igrid_sim%rbuffxC(:,:,:,1)*this%igrid_sim%tauSGS_ij(:,:,:,5)

        call this%ddz_R2R(this%igrid_sim%wC, this%igrid_sim%rbuffxC(:,:,:,1)); 
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%igrid_sim%rbuffxC(:,:,:,1)*this%igrid_sim%tauSGS_ij(:,:,:,6)


        ! 7. Actuator disk/Turbine sink  (J)
        call this%igrid_sim%spectC%ifft(this%uturb,this%igrid_sim%rbuffxC(:,:,:,1))
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,1)
 
    end subroutine 
    
    subroutine DumpBudget3(this)
        class(budgets_time_avg), intent(inout), target :: this
        integer :: idx
        real(rkind), dimension(:,:,:), pointer :: Umn, Vmn, Wmn, R11, R12, R13, R22, R23, R33
        real(rkind), dimension(:,:,:), pointer :: Pmn, tau11, tau12, tau13, tau22, tau23, tau33
        real(rkind), dimension(:,:,:), pointer :: buff, buff2


        ! Get the average from sum
        this%budget_0 = this%budget_0/(real(this%counter,rkind) + 1.d-18)
        this%budget_0(:,:,:,4)  = this%budget_0(:,:,:,4)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
        this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
        this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  - this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
        this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
        this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  - this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
        this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  - this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
        this%budget_1 = this%budget_1/(real(this%counter,rkind) + 1.d-18)
        this%budget_3 = this%budget_3/(real(this%counter,rkind) + 1.d-18)

        Umn => this%budget_0(:,:,:,1);    Vmn => this%budget_0(:,:,:,2);      Wmn => this%budget_0(:,:,:,3);
        R11 => this%budget_0(:,:,:,4);    R12 => this%budget_0(:,:,:,5);      R13 => this%budget_0(:,:,:,6)
        R22 => this%budget_0(:,:,:,7);    R23 => this%budget_0(:,:,:,8);      R33 => this%budget_0(:,:,:,9)
        Pmn => this%budget_0(:,:,:,10);   tau11 => this%budget_0(:,:,:,11);   tau12 => this%budget_0(:,:,:,12)
        tau13 => this%budget_0(:,:,:,13); tau22 => this%budget_0(:,:,:,14);   tau23 => this%budget_0(:,:,:,15)
        tau33 => this%budget_0(:,:,:,16); buff => this%igrid_sim%rbuffxC(:,:,:,1); buff2 => this%igrid_sim%rbuffxC(:,:,:,2)

        ! 1. TKE production              (G)
        this%budget_3(:,:,:,1) = -this%budget_2(:,:,:,1)

        ! 2. convective transport        (B)
        buff2 = half*(R11*R11 + R22*R22 + R33*R33)
        call this%ddx_R2R(buff2,buff); this%budget_3(:,:,:,2) = -Umn*buff
        call this%ddy_R2R(buff2,buff); this%budget_3(:,:,:,2) = this%budget_3(:,:,:,2) - Vmn*buff
        call this%ddz_R2R(buff2,buff); this%budget_3(:,:,:,2) = this%budget_3(:,:,:,2) - Wmn*buff

        ! 3. turbulent transport         (C)
        this%budget_3(:,:,:,3) = this%budget_3(:,:,:,3) - this%budget_3(:,:,:,2) - this%budget_2(:,:,:,2) - this%budget_2(:,:,:,3)

        ! 4. Pressure transport          (D)
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) - this%budget_2(:,:,:,4)

        ! 5. SGS + viscous transport     (E+F)
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) - this%budget_3(:,:,:,6) - this%budget_2(:,:,:,5) - this%budget_2(:,:,:,6)

        ! 6. SGS + viscous dissipation   (H+I)
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) - this%budget_2(:,:,:,6)

        ! 7. Actuator disk/Turbine sink  (J)
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) - Umn*this%budget_1(:,:,:,4)


        ! Dump the full budget 
        do idx = 1,size(this%budget_3,4)
            call this%dump_budget_field(this%budget_3(:,:,:,idx),idx,3)
        end do 


        ! Revert arrays to the correct state for Assemble (Order is very
        ! important throughout this subroutine, particularly indices 5 and 6)
        this%budget_3(:,:,:,7) = this%budget_3(:,:,:,7) + Umn*this%budget_1(:,:,:,4)
        this%budget_3(:,:,:,6) = this%budget_3(:,:,:,6) + this%budget_2(:,:,:,6)
        this%budget_3(:,:,:,5) = this%budget_3(:,:,:,5) + this%budget_3(:,:,:,6) + this%budget_2(:,:,:,5) + this%budget_2(:,:,:,6)
        this%budget_3(:,:,:,4) = this%budget_3(:,:,:,4) + this%budget_2(:,:,:,4)
        this%budget_3(:,:,:,3) = this%budget_3(:,:,:,3) + this%budget_3(:,:,:,2) + this%budget_2(:,:,:,2) + this%budget_2(:,:,:,3)
        

        nullify(Umn,Vmn,Wmn,R11,R12,R13,R22,R23,R33,Pmn,tau11,tau12,tau13,tau22,tau23,tau33,buff,buff2)

        ! Go back to sum
        this%budget_3 = this%budget_3*(real(this%counter,rkind) + 1.d-18)
        this%budget_1 = this%budget_1*(real(this%counter,rkind) + 1.d-18)
        this%budget_0(:,:,:,4)  = this%budget_0(:,:,:,4)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,1) ! R11
        this%budget_0(:,:,:,5)  = this%budget_0(:,:,:,5)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,2) ! R12
        this%budget_0(:,:,:,6)  = this%budget_0(:,:,:,6)  + this%budget_0(:,:,:,1)*this%budget_0(:,:,:,3) ! R13
        this%budget_0(:,:,:,7)  = this%budget_0(:,:,:,7)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,2) ! R22
        this%budget_0(:,:,:,8)  = this%budget_0(:,:,:,8)  + this%budget_0(:,:,:,2)*this%budget_0(:,:,:,3) ! R23
        this%budget_0(:,:,:,9)  = this%budget_0(:,:,:,9)  + this%budget_0(:,:,:,3)*this%budget_0(:,:,:,3) ! R33
        this%budget_0 = this%budget_0*(real(this%counter,rkind) + 1.d-18)

    end subroutine 

    

    ! ----------------------supproting subroutines ------------------------
    subroutine dump_budget_field(this, field, fieldID, BudgetID)
        use decomp_2d_io
        class(budgets_time_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(in) :: field
        integer, intent(in) :: fieldID, BudgetID
        character(len=clen) :: fname, tempname 

        write(tempname,"(A3,I2.2,A7,I1.1,A5,I2.2,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_budget",BudgetID,"_term",fieldID,"_t",this%igrid_sim%step,"_n",this%counter,".s3D"
        fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)

        call decomp_2d_write_one(1,field,fname, this%igrid_sim%gpC)

    end subroutine 
    
    subroutine restartBudget(this, rid, tid, cid)
        class(budgets_time_avg), intent(inout) :: this
        integer, intent(in) :: rid, cid, tid
        !character(len=clen) :: fname, tempname 

        this%counter = cid

        ! << Incomplete for now - write after completing dumpbudget and look at
        ! budget_xy_avg for reference. >>
    end subroutine 
    
    subroutine ResetBudget(this)
        class(budgets_time_avg), intent(inout) :: this
        
        this%counter = 0
        this%budget_0 = 0.d0 
        this%budget_1 = 0.d0 
        this%budget_2 = 0.d0 
        this%budget_3 = 0.d0 
    end subroutine 
    
    subroutine destroy(this)
        class(budgets_time_avg), intent(inout) :: this

        nullify(this%igrid_sim)
        if(this%do_budgets) then
            deallocate(this%uc, this%vc, this%wc, this%usgs, this%vsgs, this%wsgs, this%px, this%py, this%pz, this%uturb)  
            deallocate(this%budget_0, this%budget_1)
        end if

    end subroutine 

    ! ----------------------private derivative operators ------------------------
    subroutine ddx_R2R(this, f, dfdx)
        class(budgets_time_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(in) :: f
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(out) :: dfdx
        
        call this%igrid_sim%spectC%fft(f,this%igrid_sim%cbuffyC(:,:,:,1))
        call this%igrid_sim%spectC%mtimes_ik1_ip(this%igrid_sim%cbuffyC(:,:,:,1))
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1), dfdx)
    end subroutine 

    subroutine ddx_C2R(this, fhat, dfdx)
        class(budgets_time_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%spectC%spectdecomp%ysz(1),this%igrid_sim%spectC%spectdecomp%ysz(2),this%igrid_sim%spectC%spectdecomp%ysz(3)), intent(in) :: fhat
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(out) :: dfdx
        
        call this%igrid_sim%spectC%mtimes_ik1_oop(fhat,this%igrid_sim%cbuffyC(:,:,:,1))
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1), dfdx)
    end subroutine 
    
    subroutine ddy_R2R(this, f, dfdy)
        class(budgets_time_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(in) :: f
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(out) :: dfdy
        
        call this%igrid_sim%spectC%fft(f,this%igrid_sim%cbuffyC(:,:,:,1))
        call this%igrid_sim%spectC%mtimes_ik2_ip(this%igrid_sim%cbuffyC(:,:,:,1))
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1), dfdy)
    end subroutine 

    subroutine ddy_C2R(this, fhat, dfdy)
        class(budgets_time_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%spectC%spectdecomp%ysz(1),this%igrid_sim%spectC%spectdecomp%ysz(2),this%igrid_sim%spectC%spectdecomp%ysz(3)), intent(in) :: fhat
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(out) :: dfdy
        
        call this%igrid_sim%spectC%mtimes_ik2_oop(fhat,this%igrid_sim%cbuffyC(:,:,:,1))
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1), dfdy)
    end subroutine 
    
    subroutine ddz_R2R(this, f, dfdz)
        class(budgets_time_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(in) :: f
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(out) :: dfdz
        
        call transpose_x_to_y(f,this%igrid_sim%rbuffyC(:,:,:,1),this%igrid_sim%gpC)
        call transpose_y_to_z(this%igrid_sim%rbuffyC(:,:,:,1),this%igrid_sim%rbuffzC(:,:,:,1),this%igrid_sim%gpC)
        call this%igrid_sim%Pade6opZ%ddz_C2C(this%igrid_sim%rbuffzC(:,:,:,1),this%igrid_sim%rbuffzC(:,:,:,2),0,0)
        call transpose_z_to_y(this%igrid_sim%rbuffzC(:,:,:,2),this%igrid_sim%rbuffyC(:,:,:,1),this%igrid_sim%gpC)
        call transpose_y_to_x(this%igrid_sim%rbuffyC(:,:,:,1),dfdz,this%igrid_sim%gpC)

    end subroutine 
    
    subroutine ddz_C2R(this, fhat, dfdz)
        class(budgets_time_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%spectC%spectdecomp%ysz(1),this%igrid_sim%spectC%spectdecomp%ysz(2),this%igrid_sim%spectC%spectdecomp%ysz(3)), intent(in) :: fhat
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(out) :: dfdz
        
        call transpose_y_to_z(fhat,this%igrid_sim%cbuffzC(:,:,:,1),this%igrid_sim%sp_gpC)
        call this%igrid_sim%Pade6opZ%ddz_C2C(this%igrid_sim%cbuffzC(:,:,:,1),this%igrid_sim%cbuffzC(:,:,:,2),0,0)
        call transpose_z_to_y(this%igrid_sim%cbuffzC(:,:,:,2),this%igrid_sim%cbuffyC(:,:,:,1),this%igrid_sim%sp_gpC)
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1), dfdz)
        
    end subroutine 

    subroutine interp_Edge2Cell(this, fE, fC)
        class(budgets_time_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpE%xsz(1),this%igrid_sim%gpE%xsz(2),this%igrid_sim%gpE%xsz(3)), intent(in) :: fE
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(out) :: fC

        call transpose_x_to_y(fE,this%igrid_sim%rbuffyE(:,:,:,1),this%igrid_sim%gpE)
        call transpose_y_to_z(this%igrid_sim%rbuffyE(:,:,:,1),this%igrid_sim%rbuffzE(:,:,:,1),this%igrid_sim%gpE)
        call this%igrid_sim%Pade6opZ%interpz_E2C(this%igrid_sim%rbuffzE(:,:,:,1),this%igrid_sim%rbuffzC(:,:,:,2),0,0)
        call transpose_z_to_y(this%igrid_sim%rbuffzC(:,:,:,2),this%igrid_sim%rbuffyC(:,:,:,1),this%igrid_sim%gpC)
        call transpose_y_to_x(this%igrid_sim%rbuffyC(:,:,:,1),fC,this%igrid_sim%gpC)
        
    end subroutine 
end module 
