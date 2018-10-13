module budgets_xy_avg_mod
   use kind_parameters, only: rkind, clen, mpirkind
   use decomp_2d
   use cd06staggstuff, only: cd06stagg
   use reductions, only: p_sum
   use incompressibleGrid, only: igrid  
   use exits, only: message
   use basic_io, only: read_2d_ascii, write_2d_ascii
   use mpi 

   implicit none 

   private
   public :: budgets_xy_avg

   ! BUDGET TYPE: 
   ! BUDGET_0: 6 Reynolds stress terms + 3 temp fluxes + meanU + meanV + meanT
   ! BUDGET_1: momentum equation terms (Budget0 also computed) 
   ! BUDGET_2: MKE budget (Budget0 and Budget1 also computed) 
   ! BUDGET 3: TKE budget (Budget 0, 1 and 2 also included)


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



   type :: budgets_xy_avg
        private
        integer :: budgetType = 1, run_id, nz

        complex(rkind), dimension(:,:,:), allocatable :: uc, vc, wc, usgs, vsgs, wsgs, uvisc, vvisc, wvisc, px, py, pz, wb, ucor, vcor, wcor, uturb 
        type(igrid), pointer :: igrid_sim 
        real(rkind), dimension(:), allocatable :: U_mean, V_mean, dUdz, dVdz, uw, vw
        real(rkind), dimension(:), allocatable :: dUdzE, dVdzE
        
        real(rkind), dimension(:,:), allocatable :: budget_0, budget_1, budget_2
        real(rkind), dimension(:,:), allocatable :: budget_3s, budget_3, budget_0s, budget_1s
        integer :: counter
        real(rkind) :: avgFact 
        character(len=clen) :: budgets_dir
        real(rkind), dimension(:), allocatable :: tmp_meanC, tmp_meanE
        real(rkind), dimension(:,:,:), allocatable :: tmpC_1d, tmpE_1d, ddz_tmpC_1d
        real(rkind), dimension(:), allocatable :: buoy_hydrostatic
        real(rkind), dimension(:), allocatable :: wTh, P_mean, tau_13_mean, tau_23_mean, tau_33_mean

        real(rkind), dimension(:), allocatable :: meanZ_bcast

        type(cd06stagg) :: cd06op_z_mom_budget 

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
        procedure, private  :: AssembleBudget0
        procedure, private  :: AssembleBudget1
        procedure, private  :: AssembleBudget2
        procedure, private  :: AssembleBudget3
        procedure, private  :: get_xy_fluctE_from_fhatE
        procedure, private  :: get_xy_fluctC_from_fhatC
        procedure, private  :: get_xy_meanC_from_fhatC 
        procedure, private  :: get_xy_meanE_from_fhatE 
        procedure, private  :: get_xy_meanC_from_fC 
        procedure, private  :: get_xy_meanE_from_fE 
        procedure, private  :: interp_1d_Edge2Cell 
        procedure, private  :: ddz_1d_Cell2Cell

   end type 


contains 

    subroutine init(this, inputfile, igrid_sim) 
        class(budgets_xy_avg), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 
        type(igrid), intent(inout), target :: igrid_sim 
        
        character(len=clen) :: budgets_dir = "NULL"
        integer :: ioUnit, ierr,  budgetType = 1, restart_tid = 0, restart_rid = 0, restart_counter = 0
        logical :: restart_budgets = .false. 
        integer :: tidx_compute = 1000000, tidx_dump = 1000000, tidx_budget_start = 0
        logical :: do_budgets = .false. 
        namelist /BUDGET_XY_AVG/ budgetType, budgets_dir, restart_budgets, restart_rid, restart_tid, restart_counter, tidx_dump, tidx_compute, do_budgets, tidx_budget_start 
        
        ! STEP 1: Read in inputs, link pointers and allocate budget vectors
        ioUnit = 534
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=BUDGET_XY_AVG)
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
        this%avgFact = 1.d0/(real(igrid_sim%nx,rkind)*real(igrid_sim%ny,rkind))
        allocate(this%Budget_0s(this%nz,21))
        allocate(this%Budget_0(this%nz,21))
        allocate(this%Budget_1(this%nz,14))
        allocate(this%Budget_1s(this%nz,14))
        allocate(this%Budget_2(this%nz,7))
        allocate(this%Budget_3(this%nz,8))
        allocate(this%Budget_3s(this%nz,8))

        if ((trim(budgets_dir) .eq. "null") .or.(trim(budgets_dir) .eq. "NULL")) then
            this%budgets_dir = igrid_sim%outputDir
        end if 

        if (restart_budgets) then
            call this%RestartBudget(restart_rid, restart_tid, restart_counter)
        else
            call this%resetBudget()
        end if

        call this%cd06op_z_mom_budget%init(igrid_sim%nx, igrid_sim%dx, .false., .false., .true., .true.)
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

    end subroutine 


    subroutine doBudgets(this)
        class(budgets_xy_avg), intent(inout) :: this

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

    subroutine ResetBudget(this)
        class(budgets_xy_avg), intent(inout) :: this
        
        this%counter = 0
        this%budget_0 = 0.d0 
        this%budget_1 = 0.d0 
        this%budget_2 = 0.d0 
        this%budget_3 = 0.d0 

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
            call this%ddz_1d_Cell2Cell(this%U_mean, this%dUdz)
            call this%ddz_1d_Cell2Cell(this%V_mean, this%dVdz)
            
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

        end if 
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
    end subroutine 


    subroutine get_xy_meanC_from_fhatC(this, fhat, fmean)
        class(budgets_xy_avg), intent(inout) :: this
        complex(rkind), dimension(this%igrid_sim%sp_gpC%ysz(1), this%igrid_sim%sp_gpC%ysz(2), this%igrid_sim%sp_gpC%ysz(3)), intent(in) :: fhat
        real(rkind), dimension(this%nz), intent(out) :: fmean
        integer :: ierr

        call transpose_y_to_z(fhat, this%igrid_sim%cbuffzC(:,:,:,1), this%igrid_sim%sp_gpC)
        if (nrank == 0) then
            fmean = this%igrid_sim%cbuffzC(1,1,:,1)*this%avgFact
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
            fmean = this%igrid_sim%cbuffzE(1,1,:,1)*this%avgFact
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
            fmean = this%igrid_sim%cbuffzC(1,1,:,1)*this%avgFact
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
            fmean = this%igrid_sim%cbuffzE(1,1,:,1)*this%avgFact
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
    
    subroutine interp_1d_Edge2Cell(this, fE,fC)
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%nz+1), intent(in)  :: fE
        real(rkind), dimension(this%nz  ), intent(out) :: fC


        this%tmpE_1d(1,1,:) = fE
        call this%cd06op_z_mom_budget%interpz_E2C(this%tmpE_1d, this%tmpC_1d, 1, 1)
        fC = this%tmpC_1d(1,1,:)

    end subroutine 

    subroutine ddz_1d_Cell2Cell(this, fC, ddz_fC)
        class(budgets_xy_avg), intent(inout) :: this
        real(rkind), dimension(this%nz  ), intent(in)  :: fC
        real(rkind), dimension(this%nz  ), intent(out) :: ddz_fC

        this%tmpC_1d(1,1,:) = fC
        call this%cd06op_z_mom_budget%ddz_C2C(this%tmpC_1d, this%ddz_tmpC_1d,1,1)
        ddz_fC = this%ddz_tmpC_1d(1,1,:)

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

        !! STEP 2: Get Reynolds stresses (IMPORTANT: need to correct for fluctuation before dumping)
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

        ! Get mean pressure
        ! call this%get_
        call this%get_xy_meanC_from_fC(this%igrid_sim%pressure, this%tmp_meanC)
        this%budget_0(:,17) = this%budget_0(:,17) + this%tmp_meanC


        ! Get the remaining tau_sgs

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


    subroutine AssembleBudget2(this)
        class(budgets_xy_avg), intent(inout) :: this

        if (this%counter > 0) then

            this%U_mean = this%budget_0(:,1)/real(this%counter,rkind)
            this%V_mean = this%budget_0(:,2)/real(this%counter,rkind)

            this%uw = this%budget_0(:,6)/real(this%counter,rkind)
            this%vw = this%budget_0(:,8)/real(this%counter,rkind)

            call this%ddz_1d_Cell2Cell(this%U_mean, this%dUdz)
            call this%ddz_1d_Cell2Cell(this%V_mean, this%dVdz)

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
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
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
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
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
        call this%igrid_sim%spectE%ifft(this%wsgs,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
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
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
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
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_3(:,7) = this%budget_3(:,7) + this%tmp_meanC 


        ! Buoyancy transfer
        call this%igrid_sim%spectE%ifft(this%wb,this%igrid_sim%rbuffxE(:,:,:,1))
        this%igrid_sim%rbuffxE(:,:,:,1) = this%igrid_sim%w*this%igrid_sim%rbuffxE(:,:,:,1)
        call this%get_xy_meanE_from_fE(this%igrid_sim%rbuffxE(:,:,:,1), this%tmp_meanE)
        call this%interp_1d_Edge2Cell(this%tmp_meanE, this%tmp_meanC)
        this%budget_3(:,8) = this%budget_3(:,8) + this%tmp_meanC 

    end subroutine


end module 
