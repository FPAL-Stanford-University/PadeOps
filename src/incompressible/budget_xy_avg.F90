module budgets_xy_avg_mod
   use kind_parameters, only: rkind, clen, mpirkind
   use decomp_2d
   use reductions, only: p_sum
   use incompressibleGrid, only: igrid  
   use exits, only: message, GracefulExit
   use constants, only: half,two,one,zero
   use basic_io, only: read_2d_ascii, write_2d_ascii
   use mpi 

   implicit none 

   private
   public :: budgets_xy_avg

   external :: MPI_REDUCE, MPI_BCAST
   
   ! BUDGET TYPE: 
   ! BUDGET_0: 6 Reynolds stress terms + 3 temp fluxes + meanU + meanV + meanT
   ! BUDGET_1: momentum equation terms (Budget0 also computed) 
   ! BUDGET_2: MKE budget (Budget0 and Budget1 also computed) 
   ! BUDGET 3: TKE budget (Budget 0, 1 and 2 also included)
   ! BUDGET 4: RS budgets (Budget 0, 1, 2, and 4 also included)


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
   ! 17: <P>
   ! 18: <tau11> 
   ! 19: <tau12> 
   ! 20: <tau22> 
   ! 21: <tau33> 


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


   ! BUDGET_2 term indices: 
   ! 1:  Loss to Resolved TKE
   ! 2:  Convective transport 
   ! 3:  Loss to SGS TKE + viscous dissipation 
   ! 4:  SGS + viscous transport
   ! 5:  Actuator disk sink
   ! 6:  Geostrophic forcing source
   ! 7:  Coriolis work (should be zero)


   ! BUDGET_3 term indices:
   ! 1. TKE production
   ! 2. convective transport 
   ! 3. Pressure transport
   ! 4. SGS + viscous transport
   ! 5. SGS + viscous dissipation
   ! 6. Actuator disk/Turbine sink 
   ! 7. Coriolis work (should be zero)
   ! 8. Buoyancy transfer


   ! BUDGET_4_ij term indices: 
   ! 1. Shear Production
   ! 2. Turbulent transport
   ! 3. Pressure strain
   ! 4. Pressure transport 
   ! 5. SGS + Viscous transport
   ! 6. SGS + Viscous dissipation
   ! 7. Buoyancy contribution
   ! 8. Coriolis exchange
   ! 9. Actuatir disk sink/source

   type :: budgets_xy_avg
        private
        integer :: budgetType = 1, run_id, nz, num_vars_rz

        complex(rkind), dimension(:,:,:), allocatable :: uc, vc, wc, usgs, vsgs, wsgs, uvisc, vvisc, wvisc, px, py, pz, wb, ucor, vcor, wcor, uturb
        type(igrid), pointer :: igrid_sim 
        real(rkind), dimension(:), allocatable :: U_mean, V_mean, dUdz, dVdz, uw, vw
        real(rkind), dimension(:), allocatable :: dUdzE, dVdzE
        
        real(rkind), dimension(:,:), allocatable :: budget_0, budget_1, budget_2, budget_4s
        real(rkind), dimension(:,:), allocatable :: budget_4_13, budget_4_23, budget_4_12, budget_4_11, budget_4_22, budget_4_33
        real(rkind), dimension(:,:), allocatable :: budget_3s, budget_3, budget_0s, budget_1s, mean_qty
        integer :: counter
        real(rkind) :: avgFact 
        character(len=clen) :: budgets_dir
        real(rkind), dimension(:), allocatable :: tmp_meanC, tmp_meanE
        real(rkind), dimension(:,:,:), allocatable :: tmpC_1d, tmpE_1d, ddz_tmpC_1d, xspectra_mean, RFx, RFy, RFz
        real(rkind), dimension(:), allocatable :: buoy_hydrostatic
        real(rkind), dimension(:), allocatable :: wTh, P_mean, tau_13_mean, tau_23_mean, tau_33_mean

        real(rkind), dimension(:), allocatable :: meanZ_bcast
        real(rkind), dimension(:,:), allocatable :: tmpz1, tmpz2


        integer :: tidx_dump 
        integer :: tidx_compute
        integer :: tidx_budget_start 
        real(rkind) :: time_budget_start 
        logical :: do_budgets, do_spectra, do_autocorrel, forceDump

    contains
        procedure           :: init
        procedure           :: destroy
        procedure           :: ResetBudget
        procedure           :: DoBudgets
        procedure, private  :: updateBudget
        procedure, private  :: DumpBudget
        procedure, private  :: restartBudget
        procedure, private  :: AssembleBudget0
        procedure, private  :: AssembleBudget1
        procedure, private  :: AssembleBudget2
        procedure, private  :: AssembleBudget3
        procedure, private  :: AssembleBudget4_13
        procedure, private  :: AssembleBudget4_23 
        procedure, private  :: AssembleBudget4_33 
        procedure, private  :: AssembleBudget4_11 
        procedure, private  :: get_xy_fluctE_from_fhatE
        procedure, private  :: get_xy_fluctC_from_fhatC
        procedure, private  :: get_xy_meanC_from_fhatC 
        procedure, private  :: get_xy_meanE_from_fhatE 
        procedure, private  :: get_xy_meanC_from_fC 
        procedure, private  :: get_xy_meanE_from_fE 
        procedure, private  :: interp_1d_Edge2Cell 
        procedure, private  :: ddz_1d_Cell2Cell
        procedure, private  :: Assemble_spectra
        procedure, private  :: dump_spectra
        procedure, private  :: Assemble_autocorrel_x
        procedure, private  :: Assemble_autocorrel_y
        procedure, private  :: Assemble_autocorrel_z
        procedure, private  :: dump_autocorrel_x
        procedure, private  :: dump_autocorrel_y
        procedure, private  :: dump_autocorrel_z

   end type 


contains 

    subroutine init(this, inputfile, igrid_sim) 
        class(budgets_xy_avg), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 
        type(igrid), intent(inout), target :: igrid_sim 
        
        character(len=clen) :: budgets_dir = "NULL"
        integer :: ioUnit, ierr,  budgetType = 1, restart_tid = 0, restart_rid = 0, restart_counter = 0
        logical :: restart_budgets = .false. 
        integer :: tidx_compute = 1000000, tidx_dump = 1000000, tidx_budget_start = -100
        logical :: do_budgets = .false., do_spectra = .false., do_autocorrel = .false.
        real(rkind) :: time_budget_start = -1.0d0
        namelist /BUDGET_XY_AVG/ budgetType, budgets_dir, restart_budgets, restart_rid, restart_tid, restart_counter, tidx_dump, tidx_compute, do_budgets, tidx_budget_start, time_budget_start, do_spectra, do_autocorrel
        
        ! STEP 1: Read in inputs, link pointers and allocate budget vectors
        ioUnit = 534
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=BUDGET_XY_AVG)
        close(ioUnit)

        this%igrid_sim => igrid_sim 
        this%run_id = igrid_sim%runid
        this%nz = igrid_sim%nz
        this%do_budgets = do_budgets
        this%do_spectra = do_spectra
        this%do_autocorrel = do_autocorrel
        this%tidx_dump = tidx_dump
        this%tidx_compute = tidx_compute
        this%tidx_budget_start = tidx_budget_start  
        this%time_budget_start = time_budget_start  
        this%forceDump = .false.

        this%budgets_dir = budgets_dir
        this%budgetType = budgetType 
        this%avgFact = 1.d0/(real(igrid_sim%nx,rkind)*real(igrid_sim%ny,rkind))

        if((this%tidx_budget_start > 0) .and. (this%time_budget_start > 0.0d0)) then
            call GracefulExit("Both tidx_budget_start and time_budget_start in budget_xy_avg are positive. Turn one negative", 100)
        endif

        if(this%do_budgets) then
            allocate(this%Budget_0s(this%nz,21))
            allocate(this%Budget_0(this%nz,21))
            allocate(this%Budget_1(this%nz,14))
            allocate(this%Budget_1s(this%nz,14))
            allocate(this%Budget_2(this%nz,7))
            allocate(this%Budget_3(this%nz,8))
            allocate(this%Budget_3s(this%nz,8))
            
            allocate(this%Budget_4s(this%nz,9))
            allocate(this%Budget_4_13(this%nz,9))
            allocate(this%Budget_4_23(this%nz,9))
            allocate(this%Budget_4_33(this%nz,9))
            allocate(this%Budget_4_11(this%nz,9))

            if ((trim(budgets_dir) .eq. "null") .or.(trim(budgets_dir) .eq. "NULL")) then
                this%budgets_dir = igrid_sim%outputDir
            end if 

            allocate(this%mean_qty(11,1))

            if(this%do_spectra) then
                ! note :: the number of variables for which spectra are to be calculated
                ! must be smaller than nyg
                allocate(this%xspectra_mean(this%igrid_sim%sp_gpC%ysz(1),this%igrid_sim%sp_gpC%ysz(2),this%igrid_sim%sp_gpC%ysz(3)))   
            endif

            if(this%do_autocorrel) then
                allocate(this%RFx(this%igrid_sim%gpC%ysz(1), this%igrid_sim%gpC%ysz(2), this%igrid_sim%gpC%ysz(3)))
                allocate(this%RFy(this%igrid_sim%gpC%xsz(1), this%igrid_sim%gpC%xsz(2), this%igrid_sim%gpC%xsz(3)))
                this%num_vars_rz = 1 ! only u for now
                allocate(this%RFz(this%nz,this%nz,this%num_vars_rz), this%tmpz1(this%nz,this%nz), this%tmpz2(this%nz,this%nz))

            endif

            if (restart_budgets) then
                if(this%do_spectra) then
                    call GracefulExit("restart budgets not supported with do_spectra. Set one of them to false", 100)
                endif
                if(this%do_autocorrel) then
                    call GracefulExit("restart budgets not supported with do_autocorrel. Set one of them to false", 100)
                endif
                call this%RestartBudget(restart_rid, restart_tid, restart_counter)
            else
                call this%resetBudget()
            end if

            allocate(this%tmp_meanC(this%nz)) 
            allocate(this%tmp_meanE(this%nz+1)) 
            allocate(this%tmpC_1d(1,1,this%nz)) 
            allocate(this%ddz_tmpC_1d(1,1,this%nz)) 
            allocate(this%tmpE_1d(1,1,this%nz+1)) 
            allocate(this%U_mean(this%nz), this%V_mean(this%nz))
            allocate(this%dUdz(this%nz), this%dVdz(this%nz))
            allocate(this%dUdzE(this%nz+1), this%dVdzE(this%nz+1))
            allocate(this%uw(this%nz), this%vw(this%nz))
            allocate(this%wTh(this%nz))
            allocate(this%P_mean(this%nz))
            allocate(this%tau_13_mean(this%nz), this%tau_23_mean(this%nz), this%tau_33_mean(this%nz))
            allocate(this%buoy_hydrostatic(this%nz)) 
            allocate(this%meanZ_bcast(this%nz)) 


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

        end if 

    end subroutine 


    subroutine doBudgets(this, forceDump)
        class(budgets_xy_avg), intent(inout) :: this
        logical, intent(in), optional :: forceDump

        if(present(forceDump)) then
            this%forceDump = forceDump
        endif

        if (this%do_budgets) then
            if( ( (this%tidx_budget_start>0) .and. (this%igrid_sim%step>this%tidx_budget_start) ) .or. &
                ( (this%time_budget_start>0) .and. (this%igrid_sim%tsim>this%time_budget_start) ) ) then
        
                if (mod(this%igrid_sim%step,this%tidx_compute) .eq. 0) then
                    call this%updateBudget()
                end if

                if ((mod(this%igrid_sim%step,this%tidx_dump) .eq. 0) .or. this%forceDump) then
                    call this%dumpBudget()
                    call message(0,"Dumped a budget .stt file")
                end if 
            end if 
        end if 

        this%forceDump = .false. ! reset to default value

    end subroutine 

    subroutine ResetBudget(this)
        class(budgets_xy_avg), intent(inout) :: this
        
        this%counter = 0
        this%budget_0 = 0.d0 
        this%budget_1 = 0.d0 
        this%budget_2 = 0.d0 
        this%budget_3 = 0.d0 
        this%budget_4_13 = 0.d0 
        this%budget_4_23 = 0.d0 
        this%budget_4_33 = 0.d0 
        this%budget_4_11 = 0.d0 
        this%mean_qty = 0.d0 

        if(this%do_spectra) then
            this%xspectra_mean = zero
        endif

        if(this%do_autocorrel) then
            ! initialize all allocated variables
            this%RFx = zero
            this%RFy = zero
            this%RFz = zero
        endif

    end subroutine 
    
    subroutine destroy(this)
        class(budgets_xy_avg), intent(inout) :: this

        nullify(this%igrid_sim)
        if(this%do_budgets) then    !---should we ideally check if (allocated(...)) for each array before deallocating ??
            deallocate(this%uc, this%vc, this%wc, this%usgs, this%vsgs, this%wsgs, &
                       & this%uvisc, this%vvisc, this%wvisc, this%px, this%py, this%pz, this%wb, this%ucor, &
                       & this%vcor, this%wcor, this%uturb)
            deallocate(this%budget_0, this%budget_1)
            deallocate(this%mean_qty)
            if(this%do_spectra) then
                deallocate(this%xspectra_mean)
            endif
            if(this%do_autocorrel) then
                deallocate(this%RFx, this%RFy, this%RFz, this%tmpz1, this%tmpz2)
            endif
        endif

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
            ! Budget 2 need not be assembled now; it only needs to be assembled
            ! before writing to disk 
        case(3)
            call this%AssembleBudget0()
            call this%AssembleBudget1()
            call this%AssembleBudget3()
        case(4)
            call this%AssembleBudget0()
            call this%AssembleBudget1()
            call this%AssembleBudget3()
            call this%AssembleBudget4_13()
            call this%AssembleBudget4_23()
            call this%AssembleBudget4_33()
            call this%AssembleBudget4_11()
        end select

        if(this%do_spectra) then
            call this%Assemble_spectra()
        endif

        if(this%do_autocorrel) then
            call this%Assemble_autocorrel_x()
            call this%Assemble_autocorrel_y()
            call this%Assemble_autocorrel_z()
        endif

        this%counter = this%counter + 1

    end subroutine 

    subroutine DumpBudget(this)
        class(budgets_xy_avg), intent(inout) :: this
        character(len=clen) :: fname, tempname 

        if (this%budgetType>1) call this%AssembleBudget2()
        
        if (nrank == 0) then
            ! Budget 0: 
            write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_budget0","_t",this%igrid_sim%step,"_n",this%counter,".stt"
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

            ! Mean Quantities:
            write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_meanqty","_t",this%igrid_sim%step,"_n",this%counter,".stt"
            fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
            this%mean_qty = this%mean_qty / (real(this%counter,rkind) + 1.d-18)
            call write_2d_ascii(this%mean_qty, fname)
            this%mean_qty = this%mean_qty * (real(this%counter,rkind) + 1.d-18)

            ! Budget 1: 
            if (this%budgetType>0) then
                write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_budget1","_t",this%igrid_sim%step,"_n",this%counter,".stt"
                fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
                this%budget_1 = this%budget_1/(real(this%counter,rkind) + 1.d-18)
                call write_2d_ascii(this%budget_1, fname) 
                this%budget_1 = this%budget_1*(real(this%counter,rkind) + 1.d-18)
            end if

            ! Budget 2: 
            if (this%budgetType>1) then
                write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_budget2","_t",this%igrid_sim%step,"_n",this%counter,".stt"
                fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
                call write_2d_ascii(this%budget_2, fname) 
            end if

            ! Budget 3:
            ! subtract out sgs dissipation from sgs transport 
            this%U_mean = this%budget_0(:,1)/real(this%counter,rkind)
            this%V_mean = this%budget_0(:,2)/real(this%counter,rkind)
            call this%ddz_1d_Cell2Cell(this%U_mean, this%dUdz, 0, 0)
            call this%ddz_1d_Cell2Cell(this%V_mean, this%dVdz, 0, 0)
            
            this%budget_0s = this%budget_0/real(this%counter,rkind) 
            this%budget_1s = this%budget_1/real(this%counter,rkind) 
            this%budget_3s = this%budget_3/real(this%counter,rkind) 
            
            ! Production 
            this%budget_3s(:,1) = -this%budget_2(:,1)
           
            ! Convective transport 
            ! TKE transport = TOTAL KE transport - MKE trasport  
            this%budget_3s(:,2) = this%budget_3s(:,2) - this%budget_2(:,2)

            ! Pressure transport 
            ! Done already because <dpdx> and <dpdy> are 0. 

            ! SGS + Viscous dissipation
            this%budget_3s(:,5) = this%budget_3s(:,5) - this%budget_0s(:,18)*this%dUdz - this%budget_0s(:,19)*this%dVdz  

            ! SGS + Viscous transport 
            ! First get time average
            this%budget_3s(:,4) = this%budget_3s(:,4) - this%budget_1s(:,2)*this%U_mean - this%budget_1s(:,3)*this%U_mean &
                                - this%budget_1s(:,8)*this%V_mean - this%budget_1s(:,9)*this%V_mean 
            ! Now subtract out the dissipation
            this%budget_3s(:,4) = this%budget_3s(:,4) - this%budget_3s(:,5)

            ! Turbine  
            this%budget_3s(:,6) = this%budget_3s(:,6) - this%U_mean*this%budget_1s(:,6)

            ! Coriolis
            this%budget_3s(:,7) = this%budget_3s(:,7) - this%U_mean*this%budget_1s(:,4) - this%V_mean*this%budget_1s(:,10)

            ! Buoyancy 
            ! nothing to do since <W> = 0
            if (this%budgetType>2) then
                write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_budget3","_t",this%igrid_sim%step,"_n",this%counter,".stt"
                fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
                call write_2d_ascii(this%budget_3s, fname) 
            end if

            ! Budget 4: 
            ! 13 component:
            this%budget_4s = this%budget_4_13/real(this%counter,rkind) 
            ! Shear production: 
            this%budget_4s(:,1) = -this%budget_0s(:,9)*this%dUdz 
            ! Turbulent transport 
            this%budget_4s(:,2) = this%budget_4s(:,2) - (this%budget_1s(:,12)*this%U_mean) - this%budget_4s(:,1)
            ! Pressure strain 
            this%budget_4s(:,3) = this%budget_4s(:,3) - (this%budget_0s(:,17)*this%dUdz)
            ! Pressure transport 
            this%budget_4s(:,4) = this%budget_4s(:,4) - (this%budget_1s(:,14)*this%u_mean) - this%budget_4s(:,3)
            ! SGS + viscous dissipation 
            this%budget_4s(:,5) = this%budget_4s(:,5) - (this%budget_0s(:,21)*this%dUdz)
            ! SGS + viscous transport 
            this%budget_4s(:,6) = this%budget_4s(:,6) - (this%budget_1s(:,13)*this%u_mean) - this%budget_4s(:,5)
            ! rest of the indices do not need to be corrected. 
            if (this%budgetType>3) then
                write(tempname,"(A3,I2.2,A8,A3,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_budget4","_13","_t",this%igrid_sim%step,"_n",this%counter,".stt"
                fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
                call write_2d_ascii(this%budget_4s, fname) 
            end if

            ! 23 component:
            this%budget_4s = this%budget_4_23/real(this%counter,rkind) 
            ! Shear production: 
            this%budget_4s(:,1) = -this%budget_0s(:,9)*this%dVdz 
            ! Turbulent transport 
            this%budget_4s(:,2) = this%budget_4s(:,2) - (this%budget_1s(:,12)*this%v_mean) - this%budget_4s(:,1)
            ! Pressure strain 
            this%budget_4s(:,3) = this%budget_4s(:,3) - (this%budget_0s(:,17)*this%dVdz)
            ! Pressure transport 
            this%budget_4s(:,4) = this%budget_4s(:,4) - (this%budget_1s(:,14)*this%v_mean) - this%budget_4s(:,3)
            ! SGS + viscous dissipation 
            this%budget_4s(:,5) = this%budget_4s(:,5) - (this%budget_0s(:,21)*this%dVdz)
            ! SGS + viscous transport 
            this%budget_4s(:,6) = this%budget_4s(:,6) - (this%budget_1s(:,13)*this%v_mean) - this%budget_4s(:,5)
            ! rest of the indices do not need to be corrected. 
            if (this%budgetType>3) then
                write(tempname,"(A3,I2.2,A8,A3,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_budget4","_23","_t",this%igrid_sim%step,"_n",this%counter,".stt"
                fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
                call write_2d_ascii(this%budget_4s, fname) 
            end if
            
            ! 33 component:
            this%budget_4s = this%budget_4_33/real(this%counter,rkind) 
            ! Pressure transport 
            this%budget_4s(:,4) = this%budget_4s(:,4) - this%budget_4s(:,3)
            ! Viscous + SGS transport 
            this%budget_4s(:,6) = this%budget_4s(:,6) - this%budget_4s(:,5)
            if (this%budgetType>3) then
                write(tempname,"(A3,I2.2,A8,A3,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_budget4","_33","_t",this%igrid_sim%step,"_n",this%counter,".stt"
                fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
                call write_2d_ascii(this%budget_4s, fname) 
            end if

            ! 11 component:
            this%budget_4s = this%budget_4_11/real(this%counter,rkind) 
            ! Shear production: 
            this%budget_4s(:,1) = -2.d0*this%budget_0s(:,6)*this%dUdz 
            ! Turbulent transport 
            this%budget_4s(:,2) = this%budget_4s(:,2) - 2.d0*(this%budget_1s(:,1)*this%u_mean) - this%budget_4s(:,1)
            ! Pressure strain 
            ! Do nothing 
            ! Pressure transport 
            this%budget_4s(:,4) = this%budget_4s(:,4) - this%budget_4s(:,3)
            ! SGS + viscous dissipation 
            this%budget_4s(:,5) = this%budget_4s(:,5) - 2.d0*(this%budget_0s(:,14)*this%dUdz)
            ! SGS + viscous transport 
            this%budget_4s(:,6) = this%budget_4s(:,6) - (2.d0*(this%budget_1s(:,3) + this%budget_1s(:,2))*this%u_mean) - this%budget_4s(:,5)
            ! Coriolis term 
            this%budget_4s(:,8) = this%budget_4s(:,8) - 2.d0*(this%u_mean*(this%budget_1s(:,4) + this%budget_1s(:,5)))
            ! Actuator disk term 
            this%budget_4s(:,9) = this%budget_4s(:,9) - 2.d0*this%u_mean*this%budget_1s(:,6)

            if (this%budgetType>3) then
                write(tempname,"(A3,I2.2,A8,A3,A2,I6.6,A2,I6.6,A4)") "Run",this%run_id,"_budget4","_11","_t",this%igrid_sim%step,"_n",this%counter,".stt"
                fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
                call write_2d_ascii(this%budget_4s, fname) 
            end if

        end if 

        if(this%do_spectra) then
            call this%dump_spectra()
        endif

        if(this%do_autocorrel) then
            call this%dump_autocorrel_x()
            call this%dump_autocorrel_y()
            call this%dump_autocorrel_z()
        endif

    end subroutine 

    subroutine Assemble_autocorrel_x(this)
        class(budgets_xy_avg), intent(inout) :: this

        integer :: i, j, k, i2, ir, jindx

        this%igrid_sim%rbuffxC(:,:,:,1) = zero
        jindx = 1 ! u velocity
        do k = 1, this%igrid_sim%gpC%xsz(3)
          do j = 1, this%igrid_sim%gpC%xsz(2)
            do ir = 1, this%igrid_sim%gpC%xsz(1)/2 + 1
              do i = 1, this%igrid_sim%gpC%xsz(1)
                i2 = mod(i+ir-1, this%igrid_sim%gpC%xsz(1))
                if(i2==0) i2 = this%igrid_sim%gpC%xsz(1)
                !this%RFx(ir,j,k) = this%RFx(ir,j,k) + this%igrid_sim%u(i,j,k)*this%igrid_sim%u(i2,j,k)
                !this%RFx(ir,jindx,k) = this%RFx(ir,jindx,k) + sum(this%igrid_sim%u(i,:,k)*this%igrid_sim%u(i2,:,k))
                this%igrid_sim%rbuffxC(ir,j,k,1) = this%igrid_sim%rbuffxC(ir,j,k,1) + this%igrid_sim%u(i,j,k)*this%igrid_sim%u(i2,j,k)
              enddo
            enddo
          enddo
        enddo
        call transpose_x_to_y(this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
        do k = 1, this%igrid_sim%gpC%ysz(3)
          do j = 1, this%igrid_sim%gpC%ysz(2)
              this%RFx(:,jindx,k) = this%RFx(:,jindx,k) + this%igrid_sim%rbuffyC(:,j,k,1)
          enddo
        enddo

        this%igrid_sim%rbuffxC(:,:,:,1) = zero
        jindx = 2 ! v velocity
        do k = 1, this%igrid_sim%gpC%xsz(3)
          do j = 1, this%igrid_sim%gpC%xsz(2)
            do ir = 1, this%igrid_sim%gpC%xsz(1)/2 + 1
              do i = 1, this%igrid_sim%gpC%xsz(1)
                i2 = mod(i+ir-1, this%igrid_sim%gpC%xsz(1))
                if(i2==0) i2 = this%igrid_sim%gpC%xsz(1)
                this%igrid_sim%rbuffxC(ir,j,k,1) = this%igrid_sim%rbuffxC(ir,j,k,1) + this%igrid_sim%v(i,j,k)*this%igrid_sim%v(i2,j,k)
              enddo
            enddo
          enddo
        enddo
        call transpose_x_to_y(this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
        do k = 1, this%igrid_sim%gpC%ysz(3)
          do j = 1, this%igrid_sim%gpC%ysz(2)
              this%RFx(:,jindx,k) = this%RFx(:,jindx,k) + this%igrid_sim%rbuffyC(:,j,k,1)
          enddo
        enddo

        this%igrid_sim%rbuffxC(:,:,:,1) = zero
        jindx = 3 ! w velocity
        do k = 1, this%igrid_sim%gpC%xsz(3)
          do j = 1, this%igrid_sim%gpC%xsz(2)
            do ir = 1, this%igrid_sim%gpC%xsz(1)/2 + 1
              do i = 1, this%igrid_sim%gpC%xsz(1)
                i2 = mod(i+ir-1, this%igrid_sim%gpC%xsz(1))
                if(i2==0) i2 = this%igrid_sim%gpC%xsz(1)
                this%igrid_sim%rbuffxC(ir,j,k,1) = this%igrid_sim%rbuffxC(ir,j,k,1) + this%igrid_sim%wC(i,j,k)*this%igrid_sim%wC(i2,j,k)
              enddo
            enddo
          enddo
        enddo
        call transpose_x_to_y(this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
        do k = 1, this%igrid_sim%gpC%ysz(3)
          do j = 1, this%igrid_sim%gpC%ysz(2)
              this%RFx(:,jindx,k) = this%RFx(:,jindx,k) + this%igrid_sim%rbuffyC(:,j,k,1)
          enddo
        enddo

        ! pressure
        if(this%igrid_sim%fastCalcPressure .or. this%igrid_sim%storePressure) then
          this%igrid_sim%rbuffxC(:,:,:,1) = zero
          jindx = 4
          do k = 1, this%igrid_sim%gpC%xsz(3)
            do j = 1, this%igrid_sim%gpC%xsz(2)
              do ir = 1, this%igrid_sim%gpC%xsz(1)/2 + 1
                do i = 1, this%igrid_sim%gpC%xsz(1)
                  i2 = mod(i+ir-1, this%igrid_sim%gpC%xsz(1))
                  if(i2==0) i2 = this%igrid_sim%gpC%xsz(1)
                  this%igrid_sim%rbuffxC(ir,j,k,1) = this%igrid_sim%rbuffxC(ir,j,k,1) + this%igrid_sim%pressure(i,j,k)*this%igrid_sim%pressure(i2,j,k)
                enddo
              enddo
            enddo
          enddo
          call transpose_x_to_y(this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
          do k = 1, this%igrid_sim%gpC%ysz(3)
            do j = 1, this%igrid_sim%gpC%ysz(2)
                this%RFx(:,jindx,k) = this%RFx(:,jindx,k) + this%igrid_sim%rbuffyC(:,j,k,1)
            enddo
          enddo
        endif

        if(this%igrid_sim%isStratified) then
          this%igrid_sim%rbuffxC(:,:,:,1) = zero
          jindx = jindx + 1    ! T
          do k = 1, this%igrid_sim%gpC%xsz(3)
            do j = 1, this%igrid_sim%gpC%xsz(2)
              do ir = 1, this%igrid_sim%gpC%xsz(1)/2 + 1
                do i = 1, this%igrid_sim%gpC%xsz(1)
                  i2 = mod(i+ir-1, this%igrid_sim%gpC%xsz(1))
                  if(i2==0) i2 = this%igrid_sim%gpC%xsz(1)
                  this%igrid_sim%rbuffxC(ir,j,k,1) = this%igrid_sim%rbuffxC(ir,j,k,1) + this%igrid_sim%T(i,j,k)*this%igrid_sim%T(i2,j,k)
                enddo
              enddo
            enddo
          enddo
          call transpose_x_to_y(this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
          do k = 1, this%igrid_sim%gpC%ysz(3)
            do j = 1, this%igrid_sim%gpC%ysz(2)
                this%RFx(:,jindx,k) = this%RFx(:,jindx,k) + this%igrid_sim%rbuffyC(:,j,k,1)
            enddo
          enddo
        endif

    end subroutine 

    subroutine Assemble_autocorrel_y(this)
        class(budgets_xy_avg), intent(inout) :: this

        integer :: iindx, j, j2, jr, k

        this%igrid_sim%rbuffyC(:,:,:,2) = zero
        iindx = 1 ! u velocity
        call transpose_x_to_y(this%igrid_sim%u, this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
        do k = 1, this%igrid_sim%gpC%ysz(3)
          do jr = 1, this%igrid_sim%gpC%ysz(2)/2 + 1
            do j = 1, this%igrid_sim%gpC%ysz(2)
                j2 = mod(j+jr-1, this%igrid_sim%gpC%ysz(2))
                if(j2==0) j2 = this%igrid_sim%gpC%ysz(2)
                this%igrid_sim%rbuffyC(:,jr,k,2) = this%igrid_sim%rbuffyC(:,jr,k,2) + this%igrid_sim%rbuffyC(:,j,k,1)*this%igrid_sim%rbuffyC(:,j2,k,1)
            enddo
          enddo
        enddo
        call transpose_y_to_x(this%igrid_sim%rbuffyC(:,:,:,2), this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%gpC)
        do k = 1, this%igrid_sim%gpC%xsz(3)
          do j = 1, this%igrid_sim%gpC%xsz(2)
              this%RFy(iindx,j,k) = this%RFy(iindx,j,k) + sum(this%igrid_sim%rbuffxC(:,j,k,1))
          enddo
        enddo

        this%igrid_sim%rbuffyC(:,:,:,2) = zero
        iindx = 2 ! v velocity
        call transpose_x_to_y(this%igrid_sim%v, this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
        do k = 1, this%igrid_sim%gpC%ysz(3)
          do jr = 1, this%igrid_sim%gpC%ysz(2)/2 + 1
            do j = 1, this%igrid_sim%gpC%ysz(2)
                j2 = mod(j+jr-1, this%igrid_sim%gpC%ysz(2))
                if(j2==0) j2 = this%igrid_sim%gpC%ysz(2)
                this%igrid_sim%rbuffyC(:,jr,k,2) = this%igrid_sim%rbuffyC(:,jr,k,2) + this%igrid_sim%rbuffyC(:,j,k,1)*this%igrid_sim%rbuffyC(:,j2,k,1)
            enddo
          enddo
        enddo
        call transpose_y_to_x(this%igrid_sim%rbuffyC(:,:,:,2), this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%gpC)
        do k = 1, this%igrid_sim%gpC%xsz(3)
          do j = 1, this%igrid_sim%gpC%xsz(2)
              this%RFy(iindx,j,k) = this%RFy(iindx,j,k) + sum(this%igrid_sim%rbuffxC(:,j,k,1))
          enddo
        enddo

        this%igrid_sim%rbuffyC(:,:,:,2) = zero
        iindx = 3 ! w velocity
        call transpose_x_to_y(this%igrid_sim%wC, this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
        do k = 1, this%igrid_sim%gpC%ysz(3)
          do jr = 1, this%igrid_sim%gpC%ysz(2)/2 + 1
            do j = 1, this%igrid_sim%gpC%ysz(2)
                j2 = mod(j+jr-1, this%igrid_sim%gpC%ysz(2))
                if(j2==0) j2 = this%igrid_sim%gpC%ysz(2)
                this%igrid_sim%rbuffyC(:,jr,k,2) = this%igrid_sim%rbuffyC(:,jr,k,2) + this%igrid_sim%rbuffyC(:,j,k,1)*this%igrid_sim%rbuffyC(:,j2,k,1)
            enddo
          enddo
        enddo
        call transpose_y_to_x(this%igrid_sim%rbuffyC(:,:,:,2), this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%gpC)
        do k = 1, this%igrid_sim%gpC%xsz(3)
          do j = 1, this%igrid_sim%gpC%xsz(2)
              this%RFy(iindx,j,k) = this%RFy(iindx,j,k) + sum(this%igrid_sim%rbuffxC(:,j,k,1))
          enddo
        enddo

        ! pressure
        if(this%igrid_sim%fastCalcPressure .or. this%igrid_sim%storePressure) then
          this%igrid_sim%rbuffyC(:,:,:,2) = zero
          iindx = 4 ! pressure
          call transpose_x_to_y(this%igrid_sim%pressure, this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
          do k = 1, this%igrid_sim%gpC%ysz(3)
            do jr = 1, this%igrid_sim%gpC%ysz(2)/2 + 1
              do j = 1, this%igrid_sim%gpC%ysz(2)
                  j2 = mod(j+jr-1, this%igrid_sim%gpC%ysz(2))
                  if(j2==0) j2 = this%igrid_sim%gpC%ysz(2)
                  this%igrid_sim%rbuffyC(:,jr,k,2) = this%igrid_sim%rbuffyC(:,jr,k,2) + this%igrid_sim%rbuffyC(:,j,k,1)*this%igrid_sim%rbuffyC(:,j2,k,1)
              enddo
            enddo
          enddo
          call transpose_y_to_x(this%igrid_sim%rbuffyC(:,:,:,2), this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%gpC)
          do k = 1, this%igrid_sim%gpC%xsz(3)
            do j = 1, this%igrid_sim%gpC%xsz(2)
                this%RFy(iindx,j,k) = this%RFy(iindx,j,k) + sum(this%igrid_sim%rbuffxC(:,j,k,1))
            enddo
          enddo
        endif

        if(this%igrid_sim%isStratified) then
          this%igrid_sim%rbuffyC(:,:,:,2) = zero
          iindx = iindx + 1    ! T
          call transpose_x_to_y(this%igrid_sim%pressure, this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
          do k = 1, this%igrid_sim%gpC%ysz(3)
            do jr = 1, this%igrid_sim%gpC%ysz(2)/2 + 1
              do j = 1, this%igrid_sim%gpC%ysz(2)
                  j2 = mod(j+jr-1, this%igrid_sim%gpC%ysz(2))
                  if(j2==0) j2 = this%igrid_sim%gpC%ysz(2)
                  this%igrid_sim%rbuffyC(:,jr,k,2) = this%igrid_sim%rbuffyC(:,jr,k,2) + this%igrid_sim%rbuffyC(:,j,k,1)*this%igrid_sim%rbuffyC(:,j2,k,1)
              enddo
            enddo
          enddo
          call transpose_y_to_x(this%igrid_sim%rbuffyC(:,:,:,2), this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%gpC)
          do k = 1, this%igrid_sim%gpC%xsz(3)
            do j = 1, this%igrid_sim%gpC%xsz(2)
                this%RFy(iindx,j,k) = this%RFy(iindx,j,k) + sum(this%igrid_sim%rbuffxC(:,j,k,1))
            enddo
          enddo
        endif

    end subroutine 

    subroutine Assemble_autocorrel_z(this)
        class(budgets_xy_avg), intent(inout) :: this
        integer :: iindx, k, kr, k2, nsample, ierr

        iindx = 1 ! u velocity
        call transpose_x_to_y(this%igrid_sim%u,                this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
        call transpose_y_to_z(this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%rbuffzC(:,:,:,1), this%igrid_sim%gpC)
        this%tmpz1 = 0.0d0
        do k = 1, this%igrid_sim%gpC%zsz(3)
          do kr = 1, this%igrid_sim%gpC%zsz(3)
            nsample = 0

            k2 = k+(kr-1)
            if((k2>=1) .and. (k2<=this%igrid_sim%gpC%zsz(3))) then 
              nsample = nsample + 1
              this%tmpz1(kr, k) = this%tmpz1(kr, k) + sum(this%igrid_sim%rbuffzC(:,:,k,1)*this%igrid_sim%rbuffzC(:,:,k2,1))
            endif

            k2 = k-(kr-1)
            if((k2>=1) .and. (k2<=this%igrid_sim%gpC%zsz(3))) then 
              nsample = nsample + 1
              this%tmpz1(kr, k) = this%tmpz1(kr, k) + sum(this%igrid_sim%rbuffzC(:,:,k,1)*this%igrid_sim%rbuffzC(:,:,k2,1))
            endif

            this%tmpz1(kr, k) = this%tmpz1(kr, k)/(real(nsample, rkind) + 1.0d-18)
          enddo
        enddo
        call mpi_reduce(this%tmpz1, this%tmpz2, this%nz*this%nz, mpirkind, mpi_sum, 0, mpi_comm_world, ierr)
        this%RFz(:,:,iindx) = this%RFz(:,:,iindx) + this%tmpz2
        !print *, '--1111--', nrank, maxval(abs(this%tmpz2))

    end subroutine 

    subroutine dump_autocorrel_x(this)
        use decomp_2d_io
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpC%ysz(1),this%igrid_sim%gpC%ysz(2),this%igrid_sim%gpC%ysz(3)) :: tmpvar
        real(rkind) :: normfac, utmpz
        integer :: dirid, decompdir, jindx, k
        character(len=clen) :: fname, tempname 

        ! get the squares of the means
        do k = 1, size(this%budget_0, 1)
            jindx = 1 ! u velocity
            utmpz = this%budget_0(k, jindx)/(real(this%counter,rkind) + 1.d-18)
            this%igrid_sim%rbuffzC(:,jindx,k,1) = utmpz*utmpz

            jindx = 2 ! v velocity
            utmpz = this%budget_0(k, jindx)/(real(this%counter,rkind) + 1.d-18)
            this%igrid_sim%rbuffzC(:,jindx,k,1) = utmpz*utmpz

            jindx = 3 ! w velocity
            ! mean is zero; do nothing

            if(this%igrid_sim%fastCalcPressure .or. this%igrid_sim%storePressure) then
                jindx = 4 ! pressure
                utmpz = this%budget_0(k, 17)/(real(this%counter,rkind) + 1.d-18) ! note :: the index of pressure in budget_0 is 17
                this%igrid_sim%rbuffzC(:,jindx,k,1) = utmpz*utmpz
            endif

            if(this%igrid_sim%isStratified) then
                jindx = jindx + 1 ! T
                utmpz = this%budget_0(k, 3)/(real(this%counter,rkind) + 1.d-18) ! note :: the index of Temperature in budget_0 is 3
                this%igrid_sim%rbuffzC(:,jindx,k,1) = utmpz*utmpz
            endif
        enddo
        call transpose_z_to_y(this%igrid_sim%rbuffzC(:,:,:,1), this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)

        ! normalize and subtract the squares of the means
        normfac = one/real(this%igrid_sim%gpC%xsz(1) * this%igrid_sim%gpC%ysz(2) * this%counter, rkind)
        tmpvar = normfac*this%RFx
        tmpvar = tmpvar - this%igrid_sim%rbuffyC(:,:,:,1)

        dirid = 2; decompdir = 2

        jindx = 1 ! u
        write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_u_x_t",this%igrid_sim%step,".out"
        fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
        call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%gpC)

        jindx = 2 ! v
        write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_v_x_t",this%igrid_sim%step,".out"
        fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
        call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%gpC)

        jindx = 3 ! w
        write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_w_x_t",this%igrid_sim%step,".out"
        fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
        call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%gpC)

        if(this%igrid_sim%fastCalcPressure .or. this%igrid_sim%storePressure) then
            jindx = 4 ! pressure
            write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_p_x_t",this%igrid_sim%step,".out"
            fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
            call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%gpC)
        endif

        if(this%igrid_sim%isStratified) then
            jindx = jindx + 1    ! T
            write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_T_x_t",this%igrid_sim%step,".out"
            fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
            call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%gpC)
        endif

    end subroutine 

    subroutine dump_autocorrel_y(this)
        use decomp_2d_io
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)) :: tmpvar
        real(rkind) :: normfac, utmpz
        integer :: dirid, decompdir, iindx, j, k
        character(len=clen) :: fname, tempname 

        ! get the squares of the means
        do k = 1, this%igrid_sim%gpC%zsz(3)
          do j = 1, this%igrid_sim%gpC%zsz(2)
            iindx = 1 ! u velocity
            utmpz = this%budget_0(k, 1)/(real(this%counter,rkind) + 1.d-18)
            this%igrid_sim%rbuffzC(iindx,j,k,1) = utmpz*utmpz

            iindx = 2 ! v velocity
            utmpz = this%budget_0(k, 2)/(real(this%counter,rkind) + 1.d-18)
            this%igrid_sim%rbuffzC(iindx,j,k,1) = utmpz*utmpz

            iindx = 3 ! w velocity
            ! mean is zero; do nothing

            if(this%igrid_sim%fastCalcPressure .or. this%igrid_sim%storePressure) then
                iindx = 4 ! pressure
                utmpz = this%budget_0(k, 17)/(real(this%counter,rkind) + 1.d-18) ! note :: the index of pressure in budget_0 is 17
                this%igrid_sim%rbuffzC(iindx,j,k,1) = utmpz*utmpz
            endif

            if(this%igrid_sim%isStratified) then
                iindx = iindx + 1 ! T
                utmpz = this%budget_0(k, 3)/(real(this%counter,rkind) + 1.d-18) ! note :: the index of Temperature in budget_0 is 3
                this%igrid_sim%rbuffzC(iindx,j,k,1) = utmpz*utmpz
            endif
          enddo
        enddo
        call transpose_z_to_y(this%igrid_sim%rbuffzC(:,:,:,1), this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%gpC)
        call transpose_y_to_x(this%igrid_sim%rbuffyC(:,:,:,1), this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%gpC)

        ! normalize and subtract the squares of the means
        normfac = one/real(this%igrid_sim%gpC%xsz(1) * this%igrid_sim%gpC%ysz(2) * this%counter, rkind)
        tmpvar = normfac*this%RFy
        tmpvar = tmpvar - this%igrid_sim%rbuffxC(:,:,:,1)

        dirid = 1; decompdir = 1

        iindx = 1 ! u
        write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_u_y_t",this%igrid_sim%step,".out"
        fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
        call decomp_2d_write_plane(decompdir, tmpvar, dirid, iindx, fname, this%igrid_sim%gpC)

        iindx = 2 ! v
        write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_v_y_t",this%igrid_sim%step,".out"
        fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
        call decomp_2d_write_plane(decompdir, tmpvar, dirid, iindx, fname, this%igrid_sim%gpC)

        iindx = 3 ! w
        write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_w_y_t",this%igrid_sim%step,".out"
        fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
        call decomp_2d_write_plane(decompdir, tmpvar, dirid, iindx, fname, this%igrid_sim%gpC)

        if(this%igrid_sim%fastCalcPressure .or. this%igrid_sim%storePressure) then
            iindx = 4 ! pressure
            write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_p_y_t",this%igrid_sim%step,".out"
            fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
            call decomp_2d_write_plane(decompdir, tmpvar, dirid, iindx, fname, this%igrid_sim%gpC)
        endif

        if(this%igrid_sim%isStratified) then
            iindx = iindx + 1    ! T
            write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_T_y_t",this%igrid_sim%step,".out"
            fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
            call decomp_2d_write_plane(decompdir, tmpvar, dirid, iindx, fname, this%igrid_sim%gpC)
        endif

    end subroutine 

    subroutine dump_autocorrel_z(this)
        use decomp_2d_io
        class(budgets_xy_avg), intent(inout) :: this
        integer :: iindx, k, kr, k2, nsample
        real(rkind), dimension(this%nz) :: utmpz
        real(rkind) :: normfac
        character(len=clen) :: fname, tempname 

        if(nrank==0) then

          iindx = 1 ! u velocity
          utmpz = this%budget_0(:, 1)/(real(this%counter,rkind) + 1.d-18)
          this%tmpz1 = 0.0d0
          do k = 1, this%igrid_sim%gpC%zsz(3)
            do kr = 1, this%igrid_sim%gpC%zsz(3)
              nsample = 0

              k2 = k+(kr-1)
              if((k2>=1) .and. (k2<=this%igrid_sim%gpC%zsz(3))) then 
                nsample = nsample + 1
                this%tmpz1(kr, k) = this%tmpz1(kr, k) + utmpz(k)*utmpz(k2)
              endif

              k2 = k-(kr-1)
              if((k2>=1) .and. (k2<=this%igrid_sim%gpC%zsz(3))) then 
                nsample = nsample + 1
                this%tmpz1(kr, k) = this%tmpz1(kr, k) + utmpz(k)*utmpz(k2)
              endif

              this%tmpz1(kr, k) = this%tmpz1(kr, k)/(real(nsample, rkind) + 1.0d-18)
            enddo
          enddo
          normfac = one/real(this%igrid_sim%gpC%xsz(1) * this%igrid_sim%gpC%ysz(2) * this%counter, rkind)

          this%tmpz2 = normfac*this%RFz(:,:,iindx) - this%tmpz1
          !print *, '--2222--', nrank, maxval(abs(this%tmpz2)), maxval(abs(this%RFz(:,:,iindx)))
          write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorrel_u_z_t",this%igrid_sim%step,".out"
          fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
          call write_2d_ascii(this%tmpz2, fname) 

          !write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorre1_u_z_t",this%igrid_sim%step,".out"
          !fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
          !call write_2d_ascii(this%tmpz1, fname) 

          !this%tmpz2 = normfac*this%RFz(:,:,iindx) - this%tmpz1
          !write(tempname,"(A3,I2.2,A17,I6.6,A4)") "Run", this%run_id,"_autocorre3_u_z_t",this%igrid_sim%step,".out"
          !fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
          !call write_2d_ascii(this%tmpz2, fname) 

        endif

    end subroutine 

    subroutine Assemble_spectra(this)
        class(budgets_xy_avg), intent(inout) :: this
        integer :: k, j, jindx

        ! u velocity
        jindx = 1
        call this%igrid_sim%spectC%fft1_x2y(this%igrid_sim%u,this%igrid_sim%cbuffyC(:,:,:,1))
        do k = 1, size(this%igrid_sim%cbuffyC, 3)
          do j = 1, size(this%igrid_sim%cbuffyC, 2)
            this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%igrid_sim%cbuffyC(:,j,k,1))
          end do
        end do

        ! v velocity
        jindx = 2
        call this%igrid_sim%spectC%fft1_x2y(this%igrid_sim%v,this%igrid_sim%cbuffyC(:,:,:,1))
        do k = 1, size(this%igrid_sim%cbuffyC, 3)
          do j = 1, size(this%igrid_sim%cbuffyC, 2)
            this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%igrid_sim%cbuffyC(:,j,k,1))
          end do
        end do

        ! w velocity
        jindx = 3
        call this%igrid_sim%spectC%fft1_x2y(this%igrid_sim%wC,this%igrid_sim%cbuffyC(:,:,:,1))
        do k = 1, size(this%igrid_sim%cbuffyC, 3)
          do j = 1, size(this%igrid_sim%cbuffyC, 2)
            this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%igrid_sim%cbuffyC(:,j,k,1))
          end do
        end do

        ! TKE
        jindx = 4
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%u*this%igrid_sim%u
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%rbuffxC(:,:,:,1) + this%igrid_sim%v*this%igrid_sim%v
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%rbuffxC(:,:,:,1) + this%igrid_sim%wC*this%igrid_sim%wC
        !this%igrid_sim%rbuffxC(:,:,:,1) = half*this%igrid_sim%rbuffxC(:,:,:,1) ! not required because this only changes the mean

        call this%igrid_sim%spectC%fft1_x2y(this%igrid_sim%rbuffxC(:,:,:,1), this%igrid_sim%cbuffyC(:,:,:,1))
        do k = 1, size(this%igrid_sim%cbuffyC, 3)
          do j = 1, size(this%igrid_sim%cbuffyC, 2)
            this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%igrid_sim%cbuffyC(:,j,k,1))
          end do
        end do

        ! pressure
        if(this%igrid_sim%fastCalcPressure .or. this%igrid_sim%storePressure) then
            jindx = 5
            call this%igrid_sim%spectC%fft1_x2y(this%igrid_sim%pressure,this%igrid_sim%cbuffyC(:,:,:,1))
            do k = 1, size(this%igrid_sim%cbuffyC, 3)
              do j = 1, size(this%igrid_sim%cbuffyC, 2)
                this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%igrid_sim%cbuffyC(:,j,k,1))
              end do
            end do
        endif

        if(this%igrid_sim%isStratified) then
            jindx = jindx + 1    ! T
            call this%igrid_sim%spectC%fft1_x2y(this%igrid_sim%T,this%igrid_sim%cbuffyC(:,:,:,1))
            do k = 1, size(this%igrid_sim%cbuffyC, 3)
              do j = 1, size(this%igrid_sim%cbuffyC, 2)
                this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%igrid_sim%cbuffyC(:,j,k,1))
              end do
            end do
        endif

    end subroutine 

    subroutine dump_spectra(this)
        use decomp_2d_io
        class(budgets_xy_avg), intent(inout) :: this
        integer :: dirid, decompdir, nspectra, jindx
        real(rkind) :: normfac
        character(len=clen) :: fname, tempname 
        real(rkind), dimension(this%igrid_sim%sp_gpC%ysz(1),this%igrid_sim%sp_gpC%ysz(2),this%igrid_sim%sp_gpC%ysz(3)) :: tmpvar

        ! Dump horizontally averaged x-spectra
        dirid = 2; decompdir = 2

        ! --- only 4, 5 or 6 planes of xspextra_mean in y-direction are being used
        nspectra = 4
        if(this%igrid_sim%fastCalcPressure .or. this%igrid_sim%storePressure) nspectra = nspectra + 1
        if(this%igrid_sim%isStratified)                                       nspectra = nspectra + 1

        ! --- for k1 = 1, multiplication factor is 1.0,
        ! --- for k1 = 2:Nx/2+1, multiplication factor is 2.0
        normfac = two/real(size(this%igrid_sim%cbuffyC(:,:,:,1),2),rkind)/real(this%counter, rkind)
        tmpvar(1:this%igrid_sim%sp_gpC%ysz(1),1:nspectra,:) = normfac*this%xspectra_mean(1:this%igrid_sim%sp_gpC%ysz(1),1:nspectra,:)
        !this%igrid_sim%cbuffyC(1:this%igrid_sim%sp_gpC%ysz(1),1:nspectra,:,2) = normfac*this%xspectra_mean(1:this%igrid_sim%sp_gpC%ysz(1),1:nspectra,:)
        if(this%igrid_sim%sp_gpC%yst(1)==1) then
            tmpvar(1,1:nspectra,:) = half*tmpvar(1,1:nspectra,:)
            !this%igrid_sim%cbuffyC(1,1:nspectra,:,2) = half*this%igrid_sim%cbuffyC(1,1:nspectra,:,2)
        endif

        ! u Velocity
        jindx = 1 
        write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%run_id,"_specu_t",this%igrid_sim%step,".out"
        fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
        !call decomp_2d_write_plane(decompdir, this%igrid_sim%cbuffyC(:,:,:,2), dirid, jindx, fname, this%igrid_sim%sp_gpC)
        call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%sp_gpC)

        ! v Velocity
        jindx = 2 
        write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%run_id,"_specv_t",this%igrid_sim%step,".out"
        fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
        !call decomp_2d_write_plane(decompdir, this%igrid_sim%cbuffyC(:,:,:,2), dirid, jindx, fname, this%igrid_sim%sp_gpC)
        call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%sp_gpC)

        ! w Velocity
        jindx = 3 
        write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%run_id,"_specw_t",this%igrid_sim%step,".out"
        fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
        !call decomp_2d_write_plane(decompdir, this%igrid_sim%cbuffyC(:,:,:,2), dirid, jindx, fname, this%igrid_sim%sp_gpC)
        call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%sp_gpC)

        ! TKE
        jindx = 4 
        write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%run_id,"_speck_t",this%igrid_sim%step,".out"
        fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
        !call decomp_2d_write_plane(decompdir, this%igrid_sim%cbuffyC(:,:,:,2), dirid, jindx, fname, this%igrid_sim%sp_gpC)
        call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%sp_gpC)

        if(this%igrid_sim%fastCalcPressure .or. this%igrid_sim%storePressure) then
            jindx = jindx + 1 ! p
            write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%run_id,"_specp_t",this%igrid_sim%step,".out"
            fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
            !call decomp_2d_write_plane(decompdir, this%igrid_sim%cbuffyC(:,:,:,2), dirid, jindx, fname, this%igrid_sim%sp_gpC)
            call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%sp_gpC)
        endif

        if(this%igrid_sim%isStratified) then
            jindx = jindx + 1 ! T
            write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%run_id,"_specT_t",this%igrid_sim%step,".out"
            fname = this%budgets_dir(:len_trim(this%budgets_dir))//"/"//trim(tempname)
            !call decomp_2d_write_plane(decompdir, this%igrid_sim%cbuffyC(:,:,:,2), dirid, jindx, fname, this%igrid_sim%sp_gpC)
            call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%igrid_sim%sp_gpC)
        endif

    end subroutine 

    subroutine restartBudget(this, rid, tid, cid)
        class(budgets_xy_avg), intent(inout) :: this
        integer, intent(in) :: rid, cid, tid
        character(len=clen) :: fname, tempname 

        this%counter = cid
        
        ! Budget 0: 
        write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",rid,"_budget0","_t",tid,"_n",cid,".stt"
        fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
        if (allocated(this%budget_0)) deallocate(this%budget_0)
        call read_2d_ascii(this%budget_0, fname)
        this%budget_0 = this%budget_0*(real(cid,rkind) + 1.d-18)

        ! Mean Quantities: 
        write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",rid,"_meanqty","_t",tid,"_n",cid,".stt"
        fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
        if (allocated(this%budget_0)) deallocate(this%budget_0)
        call read_2d_ascii(this%mean_qty, fname)
        this%mean_qty = this%mean_qty*(real(cid,rkind) + 1.d-18)

        ! Budget 1: 
        if (this%budgetType>0) then
            write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",rid,"_budget1","_t",tid,"_n",cid,".stt"
            fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
            if (allocated(this%budget_1)) deallocate(this%budget_1)
            call read_2d_ascii(this%budget_1, fname)
            this%budget_1 = this%budget_1*(real(cid,rkind) + 1.d-18)
        end if 

        ! Budget 2: need not be read in since it's computed using budget_0 and
        ! budget_1 entirely. 
        
        ! Budget 3:
        if (this%budgetType>2) then
            write(tempname,"(A3,I2.2,A8,A2,I6.6,A2,I6.6,A4)") "Run",rid,"_budget3","_t",tid,"_n",cid,".stt"
            fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
            if (allocated(this%budget_3)) deallocate(this%budget_3)
            call read_2d_ascii(this%budget_3, fname)
            this%budget_3 = this%budget_3*(real(cid,rkind) + 1.d-18)
        end if

        ! Budget 4
        if (this%budgetType>3) then
            write(tempname,"(A3,I2.2,A8,A3,A2,I6.6,A2,I6.6,A4)") "Run",rid,"_budget4","_13","_t",tid,"_n",cid,".stt"
            fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
            if (allocated(this%budget_4_13)) deallocate(this%budget_4_13)
            call read_2d_ascii(this%budget_4_13, fname)
            this%budget_4_13 = this%budget_4_13*(real(cid,rkind) + 1.d-18)
            
            write(tempname,"(A3,I2.2,A8,A3,A2,I6.6,A2,I6.6,A4)") "Run",rid,"_budget4","_23","_t",tid,"_n",cid,".stt"
            fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
            if (allocated(this%budget_4_23)) deallocate(this%budget_4_23)
            call read_2d_ascii(this%budget_4_23, fname)
            this%budget_4_23 = this%budget_4_23*(real(cid,rkind) + 1.d-18)
            
            write(tempname,"(A3,I2.2,A8,A3,A2,I6.6,A2,I6.6,A4)") "Run",rid,"_budget4","_33","_t",tid,"_n",cid,".stt"
            fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
            if (allocated(this%budget_4_33)) deallocate(this%budget_4_33)
            call read_2d_ascii(this%budget_4_33, fname)
            this%budget_4_33 = this%budget_4_33*(real(cid,rkind) + 1.d-18)
            
            write(tempname,"(A3,I2.2,A8,A3,A2,I6.6,A2,I6.6,A4)") "Run",rid,"_budget4","_11","_t",tid,"_n",cid,".stt"
            fname = this%budgets_Dir(:len_trim(this%budgets_Dir))//"/"//trim(tempname)
            if (allocated(this%budget_4_11)) deallocate(this%budget_4_11)
            call read_2d_ascii(this%budget_4_11, fname)
            this%budget_4_11 = this%budget_4_11*(real(cid,rkind) + 1.d-18)

        end if 

        if(this%do_spectra) then
            ! Read in previous spectra
            ! ---- to be done -------
        endif

        if(this%do_autocorrel) then
            ! Read in previous spectra
            ! ---- to be done -------
        endif

    end subroutine 


    subroutine get_xy_meanC_from_fhatC(this, fhat, fmean)
        class(budgets_xy_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%sp_gpC%ysz(1), this%igrid_sim%sp_gpC%ysz(2), this%igrid_sim%sp_gpC%ysz(3)), intent(in) :: fhat
        real(rkind), dimension(this%nz), intent(out) :: fmean
        integer :: ierr

        call transpose_y_to_z(fhat, this%igrid_sim%cbuffzC(:,:,:,1), this%igrid_sim%sp_gpC)
        if (nrank == 0) then
            fmean = real(this%igrid_sim%cbuffzC(1,1,:,1),rkind)*this%avgFact
        else
            fmean = 0.d0 ! Only 0 processor has the actual mean  
        end if 
        
        call mpi_bcast(fmean,this%nz,mpirkind,0,mpi_comm_world, ierr)

    end subroutine 

    subroutine get_xy_meanE_from_fhatE(this, fhat, fmean)
        class(budgets_xy_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%sp_gpE%ysz(1), this%igrid_sim%sp_gpE%ysz(2), this%igrid_sim%sp_gpE%ysz(3)), intent(in) :: fhat
        real(rkind), dimension(this%nz+1), intent(out) :: fmean
        integer :: ierr

        call transpose_y_to_z(fhat, this%igrid_sim%cbuffzE(:,:,:,1), this%igrid_sim%sp_gpE)
        if (nrank == 0) then
            fmean = real(this%igrid_sim%cbuffzE(1,1,:,1),rkind)*this%avgFact
        else
            fmean = 0.d0 ! Only 0 processor has the actual mean  
        end if 
        
        call mpi_bcast(fmean,this%nz+1,mpirkind,0,mpi_comm_world, ierr)

    end subroutine 
    
    subroutine get_xy_meanC_from_fC(this, f, fmean)
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1), this%igrid_sim%gpC%xsz(2), this%igrid_sim%gpC%xsz(3)), intent(in) :: f
        real(rkind), dimension(this%nz), intent(out) :: fmean
        integer :: ierr

        call this%igrid_sim%spectC%fft(f, this%igrid_sim%cbuffyC(:,:,:,1))
        call transpose_y_to_z(this%igrid_sim%cbuffyC(:,:,:,1), this%igrid_sim%cbuffzC(:,:,:,1), this%igrid_sim%sp_gpC)
        if (nrank == 0) then
            fmean = real(this%igrid_sim%cbuffzC(1,1,:,1),rkind)*this%avgFact
        else
            fmean = 0.d0 ! Only 0 processor has the actual mean  
        end if 
        call mpi_bcast(fmean,this%nz,mpirkind,0,mpi_comm_world, ierr)

    end subroutine 
    
    subroutine get_xy_meanE_from_fE(this, f, fmean)
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%igrid_sim%gpE%xsz(1), this%igrid_sim%gpE%xsz(2), this%igrid_sim%gpE%xsz(3)), intent(in) :: f
        real(rkind), dimension(this%nz+1), intent(out) :: fmean
        integer :: ierr

        call this%igrid_sim%spectE%fft(f, this%igrid_sim%cbuffyE(:,:,:,1))
        call transpose_y_to_z(this%igrid_sim%cbuffyE(:,:,:,1), this%igrid_sim%cbuffzE(:,:,:,1), this%igrid_sim%sp_gpE)
        if (nrank == 0) then
            fmean = real(this%igrid_sim%cbuffzE(1,1,:,1),rkind)*this%avgFact
        else
            fmean = 0.d0 ! Only 0 processor has the actual mean  
        end if 
        call mpi_bcast(fmean,this%nz+1,mpirkind,0,mpi_comm_world, ierr)

    end subroutine 
    
    subroutine get_xy_fluctC_from_fhatC(this, fhat, ffluct)
        class(budgets_xy_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%sp_gpC%ysz(1), this%igrid_sim%sp_gpC%ysz(2), this%igrid_sim%sp_gpC%ysz(3)), intent(in) :: fhat
        real(rkind), dimension(this%igrid_sim%gpC%xsz(1),this%igrid_sim%gpC%xsz(2),this%igrid_sim%gpC%xsz(3)), intent(out) :: ffluct 

        this%igrid_sim%cbuffyC(:,:,:,1) = fhat
        if (this%igrid_sim%spectC%carryingZeroK) then
            this%igrid_sim%cbuffyC(1,1,:,1) = 0.d0 
        end if 
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1),ffluct)

    end subroutine 

    subroutine get_xy_fluctE_from_fhatE(this, fhat, ffluct)
        class(budgets_xy_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%sp_gpE%ysz(1), this%igrid_sim%sp_gpE%ysz(2), this%igrid_sim%sp_gpE%ysz(3)), intent(in) :: fhat
        real(rkind), dimension(this%igrid_sim%gpE%xsz(1),this%igrid_sim%gpE%xsz(2),this%igrid_sim%gpE%xsz(3)), intent(out) :: ffluct 

        this%igrid_sim%cbuffyE(:,:,:,1) = fhat
        if (this%igrid_sim%spectE%carryingZeroK) then
            this%igrid_sim%cbuffyE(1,1,:,1) = 0.d0 
        end if 
        call this%igrid_sim%spectE%ifft(this%igrid_sim%cbuffyE(:,:,:,1),ffluct)

    end subroutine
    
    subroutine interp_1d_Edge2Cell(this, fE,fC, bc_bot, bc_top)
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%nz+1), intent(in)  :: fE
        real(rkind), dimension(this%nz  ), intent(out) :: fC
        integer, intent(in) :: bc_bot, bc_top 


        this%tmpE_1d(1,1,:) = fE
        call this%igrid_sim%Pade6opZ%interp_1d_E2C(this%tmpE_1d,this%tmpC_1d, bc_bot, bc_top)
        fC = this%tmpC_1d(1,1,:)

    end subroutine 

    subroutine ddz_1d_Cell2Cell(this, fC, ddz_fC, bc_bot, bc_top)
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%nz  ), intent(in)  :: fC
        real(rkind), dimension(this%nz  ), intent(out) :: ddz_fC
        integer, intent(in) :: bc_bot, bc_top 

        this%tmpC_1d(1,1,:) = fC
        call this%igrid_sim%Pade6opZ%ddz_1d_C2C(this%tmpC_1d,this%ddz_tmpC_1d, bc_bot, bc_top)
        ddz_fC = this%ddz_tmpC_1d(1,1,:)

    end subroutine 

    subroutine AssembleBudget0(this)
        class(budgets_xy_avg), intent(inout) :: this

        ! STEP 1: Compute mean U, V and T
        call this%get_xy_meanC_from_fhatC(this%igrid_sim%uhat, this%tmp_meanC)
        this%budget_0(:,1) = this%budget_0(:,1) + this%tmp_meanC

        call this%get_xy_meanC_from_fhatC(this%igrid_sim%vhat, this%tmp_meanC)
        this%budget_0(:,2) = this%budget_0(:,2) + this%tmp_meanC
        
        if(this%igrid_sim%isStratified) then
            call this%get_xy_meanC_from_fhatC(this%igrid_sim%That, this%tmp_meanC)
            this%budget_0(:,3) = this%budget_0(:,3) + this%tmp_meanC
        endif

        !! STEP 2: Get Reynolds stresses (IMPORTANT: need to correct for fluctuation before dumping)
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%u*this%igrid_sim%u
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
        this%budget_0(:,4) = this%budget_0(:,4) + this%tmp_meanC

        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%u*this%igrid_sim%v
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
        this%budget_0(:,5) = this%budget_0(:,5) + this%tmp_meanC
        
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%uE*this%igrid_sim%w
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,-1,-1)
        this%budget_0(:,6) = this%budget_0(:,6) + this%tmp_meanC
        
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%v*this%igrid_sim%v
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
        this%budget_0(:,7) = this%budget_0(:,7) + this%tmp_meanC

        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%vE*this%igrid_sim%w
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,-1,-1)
        this%budget_0(:,8) = this%budget_0(:,8) + this%tmp_meanC
        
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%w
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_0(:,9) = this%budget_0(:,9) + this%tmp_meanC
      
        ! time- and horizontally-averaged surface quantities
        this%mean_qty(1,1) = this%mean_qty(1,1) + this%igrid_sim%sgsmodel%get_ustar()
        this%mean_qty(2,1) = this%mean_qty(2,1) + this%igrid_sim%sgsmodel%get_uw_surf()
        this%mean_qty(3,1) = this%mean_qty(3,1) + this%igrid_sim%sgsmodel%get_vw_surf()
        this%mean_qty(4,1) = this%mean_qty(4,1) + this%igrid_sim%sgsmodel%get_umean()
        this%mean_qty(5,1) = this%mean_qty(5,1) + this%igrid_sim%sgsmodel%get_vmean()
        this%mean_qty(6,1) = this%mean_qty(6,1) + this%igrid_sim%sgsmodel%get_uspeedmean()
        this%mean_qty(7,1) = this%mean_qty(7,1) + this%igrid_sim%sgsmodel%get_GlobalConstant()
        this%mean_qty(8,1) = this%mean_qty(8,1) + this%igrid_sim%sgsmodel%getMax_DynSmagConst()
        this%mean_qty(9,1) = this%mean_qty(9,1) + this%igrid_sim%sgsmodel%getMax_DynPrandtl()
 
        ! STEP 2: Get Temperature fluxes and variances (IMPORTANT: need to correct for fluctuation before dumping)
        if(this%igrid_sim%isStratified) then
            this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%u*this%igrid_sim%T
            call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
            this%budget_0(:,10) = this%budget_0(:,10) + this%tmp_meanC
            
            this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%v*this%igrid_sim%T
            call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
            this%budget_0(:,11) = this%budget_0(:,11) + this%tmp_meanC

            this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%TE*this%igrid_sim%w
            call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
            call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,-1,-1)
            this%budget_0(:,12) = this%budget_0(:,12) + this%tmp_meanC

            this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%T*this%igrid_sim%T
            call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1), this%tmp_meanC)
            this%budget_0(:,13) = this%budget_0(:,13) + this%tmp_meanC

            this%mean_qty(10,1) = this%mean_qty(10,1) + this%igrid_sim%sgsmodel%get_InvObLength()
            this%mean_qty(11,1) = this%mean_qty(11,1) + this%igrid_sim%sgsmodel%get_wTh_surf()
        endif

        ! STEP 3: SGS stress (also viscous stress if finite reynolds number is being used)
        call this%get_xy_meanE_from_fE(this%igrid_sim%tau13, this%tmp_meanE) 
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_0(:,14) = this%budget_0(:,14) + this%tmp_meanC

        call this%get_xy_meanE_from_fE(this%igrid_sim%tau23, this%tmp_meanE) 
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_0(:,15) = this%budget_0(:,15) + this%tmp_meanC

        if (associated(this%igrid_Sim%q3)) then
            call this%get_xy_meanE_from_fE(this%igrid_sim%q3, this%tmp_meanE) 
            call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
            this%budget_0(:,16) = this%budget_0(:,16) + this%tmp_meanC
        end if 

        ! Get mean pressure
        ! call this%get_
        call this%get_xy_meanC_from_fC(this%igrid_sim%pressure, this%tmp_meanC)
        this%budget_0(:,17) = this%budget_0(:,17) + this%tmp_meanC

        ! Get mean tau_11
        call this%get_xy_meanC_from_fC(this%igrid_sim%tauSGS_ij(:,:,:,1), this%tmp_meanC)
        this%budget_0(:,18) = this%budget_0(:,18) + this%tmp_meanC

        ! Get mean tau_12
        call this%get_xy_meanC_from_fC(this%igrid_sim%tauSGS_ij(:,:,:,2), this%tmp_meanC)
        this%budget_0(:,19) = this%budget_0(:,19) + this%tmp_meanC

        ! Get mean tau_22
        call this%get_xy_meanC_from_fC(this%igrid_sim%tauSGS_ij(:,:,:,4), this%tmp_meanC)
        this%budget_0(:,20) = this%budget_0(:,20) + this%tmp_meanC

        ! Get mean tau_23
        call this%get_xy_meanC_from_fC(this%igrid_sim%tauSGS_ij(:,:,:,6), this%tmp_meanC)
        this%budget_0(:,21) = this%budget_0(:,21) + this%tmp_meanC

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
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_1(:,12) = this%budget_1(:,12) + this%tmp_meanC

        call this%get_xy_meanE_from_fhatE(this%wsgs, this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_1(:,13) = this%budget_1(:,13) + this%tmp_meanC

        call this%get_xy_meanE_from_fhatE(this%pz, this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_1(:,14) = this%budget_1(:,14) + this%tmp_meanC

    end subroutine


    subroutine AssembleBudget2(this)
        class(budgets_xy_avg), intent(inout) :: this

        if (this%counter > 0) then

            this%U_mean = this%budget_0(:,1)/real(this%counter,rkind)
            this%V_mean = this%budget_0(:,2)/real(this%counter,rkind)

            this%uw = this%budget_0(:,6)/real(this%counter,rkind)
            this%vw = this%budget_0(:,8)/real(this%counter,rkind)

            call this%ddz_1d_Cell2Cell(this%U_mean, this%dUdz, 0, 0)
            call this%ddz_1d_Cell2Cell(this%V_mean, this%dVdz, 0, 0)

            ! Loss of MKE to Resolved TKE
            this%budget_2(:,1) = this%uw*this%dUdz + this%vw*this%dVdz
            
            ! Convective transport: (more accurate to compute this way than to
            ! take derivative)
            this%tmp_meanC = this%U_mean*(this%budget_1(:,1)/real(this%counter,rkind))
            this%tmp_meanC = this%tmp_meanC + this%V_mean*(this%budget_1(:,7)/real(this%counter,rkind))
            this%budget_2(:,2) = this%tmp_meanC - this%budget_2(:,1)

            ! Loss of MKE to SGS TKE + viscous dissipation 
            this%budget_2(:,3) = (this%budget_0(:,14)/real(this%counter+1,rkind))*this%dUdz &
                               + (this%budget_0(:,15)/real(this%counter+1,rkind))*this%dVdz
       
            ! Viscous + SGS transport
            this%tmp_meanC = this%U_mean*((this%budget_1(:,2) + this%budget_1(:,3))/real(this%counter,rkind))
            this%tmp_meanC = this%tmp_meanC + this%V_mean*((this%budget_1(:,8) + this%budget_1(:,9))/real(this%counter,rkind))
            this%budget_2(:,4) = this%tmp_meanC - this%budget_2(:,3)

            ! Actuator disk sink
            this%budget_2(:,5) = this%U_mean*(this%budget_1(:,6)/real(this%counter,rkind))

            ! Geostrophic forcing term
            this%budget_2(:,6) = this%U_mean*(this%budget_1(:,5)/real(this%counter,rkind)) &
                               + this%V_mean*(this%budget_1(:,11)/real(this%counter,rkind))
            
            ! Coriolis forcing term (should be 0)
            this%budget_2(:,7) = this%U_mean*(this%budget_1(:,4)/real(this%counter,rkind)) &
                               + this%V_mean*(this%budget_1(:,10)/real(this%counter,rkind))

        end if 

    end subroutine 
    
    subroutine AssembleBudget3(this)
        use IncompressibleGrid, only: wBC_bottom, wBC_top 
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind) :: tmp1, tmp2
        
        ! First get the mean statistics
        !this%U_mean = this%budget_0(:,1)/real(this%counter+1,rkind)
        !this%V_mean = this%budget_0(:,2)/real(this%counter+1,rkind)
        
        
        
        ! TKE production: 
        ! this%budget_3(:,1) = this%budget_2(:,1) <- compute this before writing 

        ! Convective transport 
        ! Really just compute total kinetic energy transport: u_i* ddx_j (u_i * u_j)  = ddx_j (u_i u_i /2 u_j) = u_i * RHS_i
        ! The TKE convetive transport is just Total KE transport - MKE transport 
        call this%igrid_sim%spectC%ifft(this%uc,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%igrid_sim%spectC%ifft(this%vc,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%rbuffxC(:,:,:,4) + this%igrid_sim%v*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_3(:,2) = this%budget_3(:,2) + this%tmp_meanC
        call this%igrid_sim%spectE%ifft(this%wc,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, -wBC_bottom, -wBC_top)
        this%budget_3(:,2) = this%budget_3(:,2) + this%tmp_meanC 
    
        ! Pressure transport 
        call this%igrid_sim%spectC%ifft(this%px,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%rbuffxC(:,:,:,3)*this%igrid_sim%u
        call this%igrid_sim%spectC%ifft(this%py,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%rbuffxC(:,:,:,4) + this%igrid_sim%rbuffxC(:,:,:,3)*this%igrid_sim%v
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_3(:,3) = this%budget_3(:,3)  + this%tmp_meanC 
        call this%igrid_sim%spectE%ifft(this%pz,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%rbuffxE(:,:,:,1)*this%igrid_sim%w
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_3(:,3) = this%budget_3(:,3)  + this%tmp_meanC 

        ! SGS transport 
        this%igrid_sim%cbuffyC(:,:,:,1) = this%usgs + this%uvisc
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1),this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,3)
        this%igrid_sim%cbuffyC(:,:,:,1) = this%vsgs + this%vvisc
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1),this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%rbuffxC(:,:,:,4) + this%igrid_sim%v*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_3(:,4) = this%budget_3(:,4) + this%tmp_meanC
        this%igrid_sim%cbuffyE(:,:,:,1) = this%wsgs + this%wvisc
        call this%igrid_sim%spectE%ifft(this%igrid_sim%cbuffyE(:,:,:,1),this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_3(:,4) = this%budget_3(:,4) + this%tmp_meanC 
        
        ! SGS dissipation
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%tauSGS_ij(:,:,:,1)*this%igrid_sim%duidxjC(:,:,:,1) &
                                        + this%igrid_sim%tauSGS_ij(:,:,:,2)*this%igrid_sim%duidxjC(:,:,:,4) &
                                        + this%igrid_sim%tauSGS_ij(:,:,:,4)*this%igrid_sim%duidxjC(:,:,:,5) &
                                        + this%igrid_sim%tauSGS_ij(:,:,:,2)*this%igrid_sim%duidxjC(:,:,:,2) &
                                        + this%igrid_sim%tauSGS_ij(:,:,:,6)*this%igrid_sim%duidxjC(:,:,:,9) 

        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1),this%tmp_meanC)
        this%budget_3(:,5) = this%budget_3(:,5) + this%tmp_meanC

        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%tau13*this%igrid_sim%duidxjE(:,:,:,3) &
                                        + this%igrid_sim%tau23*this%igrid_sim%duidxjE(:,:,:,6) &
                                        + this%igrid_sim%tau13*this%igrid_sim%duidxjE(:,:,:,7) &
                                        + this%igrid_sim%tau23*this%igrid_sim%duidxjE(:,:,:,8) 

        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1),this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_3(:,5) = this%budget_3(:,5) + this%tmp_meanC 

        ! Turbine sink
        call this%igrid_sim%spectC%ifft(this%uturb,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%rbuffxC(:,:,:,3)*this%igrid_sim%u
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_3(:,6) = this%budget_3(:,6)  + this%tmp_meanC 

        ! Coriolis sink 
        ! Get the geostrophic forcing 
        call this%igrid_sim%get_geostrophic_forcing(tmp1, tmp2)         ! Forcing in x and y directions respectively 
        
        call this%igrid_sim%spectC%ifft(this%ucor,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,3) = this%igrid_sim%rbuffxC(:,:,:,3) - tmp1 ! remove the geostrphic forcing term 
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%igrid_sim%spectC%ifft(this%vcor,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,3) = this%igrid_sim%rbuffxC(:,:,:,3) - tmp2 ! remove the geostrphic forcing term 
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%rbuffxC(:,:,:,4) + this%igrid_sim%v*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_3(:,7) = this%budget_3(:,7) + this%tmp_meanC
        call this%igrid_sim%spectE%ifft(this%wcor,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, -1, -1)
        this%budget_3(:,7) = this%budget_3(:,7) + this%tmp_meanC 


        ! Buoyancy transfer
        call this%igrid_sim%spectE%ifft(this%wb,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_3(:,8) = this%budget_3(:,8) + this%tmp_meanC 

    end subroutine

   
    subroutine AssembleBudget4_13(this)
        class(budgets_xy_avg), intent(inout) :: this
        
        ! First get the mean statistics
        !this%U_mean = this%budget_0(:,1)/real(this%counter+1,rkind)
        !this%V_mean = this%budget_0(:,2)/real(this%counter+1,rkind)
       
        !Indices: i = 1, k = 3
            
        ! Shear production: -(w*w)*dU/dz 
        ! Do nothing here (compute before writing)

        ! Turbulent transport 
        ! Assemble: w*RHSc_u + u*RHSc_w
        call this%igrid_sim%spectC%ifft(this%uc,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_13(:,2) = this%budget_4_13(:,2) + this%tmp_meanC
        call this%igrid_sim%spectE%ifft(this%wc,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%uE*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_13(:,2) = this%budget_4_13(:,2) + this%tmp_meanC

        ! pressure strain rate 
        call this%igrid_sim%interpolate_cellField_to_edgeField(this%igrid_sim%pressure,this%igrid_sim%rbuffxE(:,:,:,2),0,0) 
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%duidxjE(:,:,:,3) + this%igrid_sim%duidxjE(:,:,:,7) 
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%rbuffxE(:,:,:,1)*this%igrid_sim%rbuffxE(:,:,:,2)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_13(:,3) = this%budget_4_13(:,3) + this%tmp_meanC

        ! pressure transport 
        call this%igrid_sim%spectC%ifft(this%px,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_13(:,4) = this%budget_4_13(:,4) + this%tmp_meanC
        call this%igrid_sim%spectE%ifft(this%pz,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%uE*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_13(:,4) = this%budget_4_13(:,4) + this%tmp_meanC
        
        ! sgs dissipation 
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%tauSGS_ij(:,:,:,1)*this%igrid_sim%duidxjC(:,:,:,7) &
                                        + this%igrid_sim%tauSGS_ij(:,:,:,2)*this%igrid_sim%duidxjC(:,:,:,8) &
                                        + this%igrid_sim%tauSGS_ij(:,:,:,6)*this%igrid_sim%duidxjC(:,:,:,3) 
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1),this%tmp_meanC)
        this%budget_4_13(:,5) = this%budget_4_13(:,5) + this%tmp_meanC
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%tau13*this%igrid_sim%duidxjE(:,:,:,9) &
                                        + this%igrid_sim%tau13*this%igrid_sim%duidxjE(:,:,:,1) &
                                        + this%igrid_sim%tau23*this%igrid_sim%duidxjE(:,:,:,2) 
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_13(:,5) = this%budget_4_13(:,5) + this%tmp_meanC
        
        ! SGS + viscous transport 
        this%igrid_sim%cbuffyC(:,:,:,1) = this%usgs + this%uvisc
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1),this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_13(:,6) = this%budget_4_13(:,6) + this%tmp_meanC
        this%igrid_sim%cbuffyE(:,:,:,1) = this%wsgs + this%wvisc
        call this%igrid_sim%spectE%ifft(this%igrid_sim%cbuffyE(:,:,:,1),this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%uE*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_4_13(:,6) = this%budget_4_13(:,6) + this%tmp_meanC

        ! Buoyancy transfer: 
        call this%igrid_sim%spectE%ifft(this%wb,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%uE*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_13(:,7) = this%budget_4_13(:,7) + this%tmp_meanC

        ! Coriolis transfer: 
        call this%igrid_sim%spectC%ifft(this%ucor,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_13(:,8) = this%budget_4_13(:,8) + this%tmp_meanC
        call this%igrid_sim%spectE%ifft(this%wcor,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%uE*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_13(:,8) = this%budget_4_13(:,8) + this%tmp_meanC
        
        ! Actuator disk source/sink: 
        call this%igrid_sim%spectC%ifft(this%uturb,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_13(:,9) = this%budget_4_13(:,9) + this%tmp_meanC


    end subroutine

    subroutine AssembleBudget4_23(this)
        class(budgets_xy_avg), intent(inout) :: this
        
        ! First get the mean statistics
        !this%U_mean = this%budget_0(:,1)/real(this%counter+1,rkind)
        !this%V_mean = this%budget_0(:,2)/real(this%counter+1,rkind)
       
        !Indices: i = 2, k = 3
            
        ! Shear production: -(w*w)*dU/dz 
        ! Do nothing here (compute before writing)

        ! Turbulent transport 
        ! Assemble: w*RHSc_v + v*RHSc_w
        call this%igrid_sim%spectC%ifft(this%vc,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_23(:,2) = this%budget_4_23(:,2) + this%tmp_meanC
        call this%igrid_sim%spectE%ifft(this%wc,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%vE*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_23(:,2) = this%budget_4_23(:,2) + this%tmp_meanC

        ! pressure strain rate 
        call this%igrid_sim%interpolate_cellField_to_edgeField(this%igrid_sim%pressure,this%igrid_sim%rbuffxE(:,:,:,2),0,0) 
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%duidxjE(:,:,:,6) + this%igrid_sim%duidxjE(:,:,:,8) 
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%rbuffxE(:,:,:,1)*this%igrid_sim%rbuffxE(:,:,:,2)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_23(:,3) = this%budget_4_23(:,3) + this%tmp_meanC

        ! pressure transport 
        call this%igrid_sim%spectC%ifft(this%py,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_23(:,4) = this%budget_4_23(:,4) + this%tmp_meanC
        call this%igrid_sim%spectE%ifft(this%pz,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%vE*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_23(:,4) = this%budget_4_23(:,4) + this%tmp_meanC
        
        ! sgs dissipation 
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%tauSGS_ij(:,:,:,2)*this%igrid_sim%duidxjC(:,:,:,7) &
                                        + this%igrid_sim%tauSGS_ij(:,:,:,4)*this%igrid_sim%duidxjC(:,:,:,8) &
                                        + this%igrid_sim%tauSGS_ij(:,:,:,6)*this%igrid_sim%duidxjC(:,:,:,6) 
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1),this%tmp_meanC)
        this%budget_4_23(:,5) = this%budget_4_23(:,5) + this%tmp_meanC
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%tau13*this%igrid_sim%duidxjE(:,:,:,4) &
                                        + this%igrid_sim%tau23*this%igrid_sim%duidxjE(:,:,:,5) &
                                        + this%igrid_sim%tau23*this%igrid_sim%duidxjE(:,:,:,9) 
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_23(:,5) = this%budget_4_23(:,5) + this%tmp_meanC
        
        ! SGS + viscous transport 
        this%igrid_sim%cbuffyC(:,:,:,1) = this%vsgs + this%vvisc
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1),this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_23(:,6) = this%budget_4_23(:,6) + this%tmp_meanC
        this%igrid_sim%cbuffyE(:,:,:,1) = this%wsgs + this%wvisc
        call this%igrid_sim%spectE%ifft(this%igrid_sim%cbuffyE(:,:,:,1),this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%vE*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_4_23(:,6) = this%budget_4_23(:,6) + this%tmp_meanC

        ! Buoyancy transfer: 
        call this%igrid_sim%spectE%ifft(this%wb,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%vE*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_23(:,7) = this%budget_4_23(:,7) + this%tmp_meanC

        ! Coriolis transfer: 
        call this%igrid_sim%spectC%ifft(this%vcor,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%wC*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_23(:,8) = this%budget_4_23(:,8) + this%tmp_meanC
        call this%igrid_sim%spectE%ifft(this%wcor,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%vE*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_23(:,8) = this%budget_4_23(:,8) + this%tmp_meanC
        
        ! Actuator disk source/sink: 
        ! Doesn't show up. 

    end subroutine

    subroutine AssembleBudget4_33(this)
        class(budgets_xy_avg), intent(inout) :: this
        
        ! First get the mean statistics
        !this%U_mean = this%budget_0(:,1)/real(this%counter+1,rkind)
        !this%V_mean = this%budget_0(:,2)/real(this%counter+1,rkind)
       
        !Indices: i = 3, k = 3
            
        ! Shear production: -(w*w)*dU/dz 
        ! Do nothing here (compute before writing)

        ! Turbulent transport 
        call this%igrid_sim%spectE%ifft(this%wc,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_33(:,2) = this%budget_4_33(:,2) + 2.d0*this%tmp_meanC

        ! pressure strain rate 
        call this%igrid_sim%interpolate_cellField_to_edgeField(this%igrid_sim%pressure,this%igrid_sim%rbuffxE(:,:,:,2),0,0) 
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%duidxjE(:,:,:,9)*this%igrid_sim%rbuffxE(:,:,:,2)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_33(:,3) = this%budget_4_33(:,3) + 2.d0*this%tmp_meanC

        ! pressure transport 
        call this%igrid_sim%spectE%ifft(this%pz,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_33(:,4) = this%budget_4_33(:,4) + 2.d0*this%tmp_meanC
        
        ! sgs dissipation 
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%tauSGS_ij(:,:,:,6)*this%igrid_sim%duidxjC(:,:,:,9) 
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1),this%tmp_meanC)
        this%budget_4_33(:,5) = this%budget_4_33(:,5) + 2.d0*this%tmp_meanC
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%tau13*this%igrid_sim%duidxjE(:,:,:,7) &
                                        + this%igrid_sim%tau23*this%igrid_sim%duidxjE(:,:,:,8) 
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_33(:,5) = this%budget_4_33(:,5) + 2.d0*this%tmp_meanC
        
        ! SGS + viscous transport 
        this%igrid_sim%cbuffyE(:,:,:,1) = this%wsgs + this%wvisc
        call this%igrid_sim%spectE%ifft(this%igrid_sim%cbuffyE(:,:,:,1),this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC,0,0)
        this%budget_4_33(:,6) = this%budget_4_33(:,6) + 2.d0*this%tmp_meanC

        ! Buoyancy transfer: 
        call this%igrid_sim%spectE%ifft(this%wb,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_33(:,7) = this%budget_4_33(:,7) + 2.d0*this%tmp_meanC

        ! Coriolis transfer: 
        call this%igrid_sim%spectE%ifft(this%wcor,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_33(:,8) = this%budget_4_33(:,8) + 2.d0*this%tmp_meanC
        
        ! Actuator disk source/sink: 
        ! Doesn't show up. 

    end subroutine
    
    subroutine AssembleBudget4_11(this)
        class(budgets_xy_avg), intent(inout) :: this
        
        ! First get the mean statistics
        !this%U_mean = this%budget_0(:,1)/real(this%counter+1,rkind)
        !this%V_mean = this%budget_0(:,2)/real(this%counter+1,rkind)
       
        !Indices: i = 3, k = 3
            
        ! Shear production: -(w*w)*dU/dz 
        ! Do nothing here (compute before writing)

        ! Turbulent transport 
        call this%igrid_sim%spectC%ifft(this%uc,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_11(:,2) = this%budget_4_11(:,2) + 2.d0*this%tmp_meanC


        ! pressure strain rate 
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%duidxjC(:,:,:,1)*this%igrid_sim%pressure
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_11(:,3) = this%budget_4_11(:,3) + 2.d0*this%tmp_meanC

        ! pressure transport 
        call this%igrid_sim%spectC%ifft(this%px,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_11(:,4) = this%budget_4_11(:,4) + 2.d0*this%tmp_meanC
        
        ! sgs dissipation 
        this%igrid_sim%rbuffxC(:,:,:,1) = this%igrid_sim%tauSGS_ij(:,:,:,1)*this%igrid_sim%duidxjC(:,:,:,1) & 
                                        + this%igrid_sim%tauSGS_ij(:,:,:,2)*this%igrid_sim%duidxjC(:,:,:,2) 
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,1),this%tmp_meanC)
        this%budget_4_11(:,5) = this%budget_4_11(:,5) + 2.d0*this%tmp_meanC
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%tau13*this%igrid_sim%duidxjE(:,:,:,3) 
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC, 0, 0)
        this%budget_4_11(:,5) = this%budget_4_11(:,5) + 2.d0*this%tmp_meanC
        
        ! SGS + viscous transport
        this%igrid_sim%cbuffyC(:,:,:,1) = this%usgs + this%uvisc
        call this%igrid_sim%spectC%ifft(this%igrid_sim%cbuffyC(:,:,:,1),this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_11(:,6) = this%budget_4_11(:,6) + 2.d0*this%tmp_meanC

        ! Buoyancy transfer: 
        ! Doesn't show up 

        ! Coriolis transfer: 
        call this%igrid_sim%spectC%ifft(this%ucor,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_11(:,8) = this%budget_4_11(:,8) + 2.d0*this%tmp_meanC

        ! Actuator disk source/sink: 
        call this%igrid_sim%spectC%ifft(this%uturb,this%igrid_sim%rbuffxC(:,:,:,3))
        this%igrid_sim%rbuffxC(:,:,:,4) = this%igrid_sim%u*this%igrid_sim%rbuffxC(:,:,:,3)
        call this%get_xy_meanC_from_fC(this%igrid_sim%rbuffxC(:,:,:,4),this%tmp_meanC)
        this%budget_4_11(:,9) = this%budget_4_11(:,9) + 2.d0*this%tmp_meanC

    end subroutine


end module 
