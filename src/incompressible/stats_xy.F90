module stats_xy_mod
    ! Module which computes instantaneous xy-averaged quantities and dumps them
    ! as scheduled. Unlike the budget_xy_avg.F90, a time-average is not
    ! accumulated so the user can specify averaging windows later. This gives a
    ! stats class valid for non-stationary problems.
    use kind_parameters,    only: rkind, clen, MPIRKIND
    use incompressibleGrid, only: igrid
    use MPI
    use decomp_2d,          only: nrank, transpose_y_to_z
    use basic_io,           only: write_2d_ascii 
    use exits,              only: message
    use fortran_assert,     only: assert
    !use stdlib_sorting,     only: sort_index
    implicit none

    private
    public :: stats_xy
    type :: stats_xy
        ! The quantities computed in this class are the following:
        ! > TKE budget
        ! > Buoyancy variance budget (referred to as "TPE budget" in this code)
        !   > Note: the terms in the buoyancy variance equation need to be
        !   scaled by 1/N^2 (1 over the buoyancy frequency squared) during post processing
        ! > Reynolds stresses
        ! > TKE
        ! > Force variance
        ! > Mean u,v,w,T,p
        ! > Temperature fluxes (<u'T'>, <v'T'>, <w'T'>)
        ! > Temperatur variance
        ! > Mean temperature gradient (note, buoyancy frequency, N ~ sqrt(dTdz),
        !   which is needed for the tpe_budget above)

        private
        real(rkind), dimension(:,:,:), allocatable :: stats, q3
        real(rkind), dimension(:,:,:,:), allocatable :: stats_sca
        real(rkind), dimension(:,:,:), pointer :: T, dTdxC, dTdyC, dTdzC, q1, q2
        real(rkind), dimension(:,:), pointer :: tke_budget, TT_budget, R11_budget, &
          R22_budget, R33_budget, wT_budget, duidxj_var, dTdxj_var, tauij_var, qj_var, &
          tke_flux, sca_var_flux
        real(rkind), dimension(:), pointer :: mke, tke, fifi, meanU, meanV, meanW, meanT, meanP, &
          meanFx, meanFy, meanFz,&
          uu, vv, ww, uv, uw, vw, uT, vT, wT, TT, dTdz, S12_var, S13_var, S23_var, &
          meanq3
        real(rkind), dimension(:), allocatable :: time
        real(rkind), dimension(:,:), allocatable :: zbuff, zbuff_, qold, qnew
        integer, dimension(:), allocatable :: tid
        real(rkind), dimension(:,:,:), allocatable :: tau11,tau12,tau13,tau22,tau23,tau33

        integer :: tid_start, compute_freq, dump_freq, tidx, tid_ddt
        integer :: nwrite, nscalars
        character(len=clen) :: outputdir
        type(igrid), pointer :: sim
        real(rkind) :: avgFact
        logical :: do_stats

        contains
          procedure          :: init
          procedure          :: destroy
          procedure          :: compute_stats
          procedure, private :: dump_stats
          procedure, private :: write_time_info
          procedure, private :: write_labels
          procedure, private :: link_pointers
          procedure, private :: get_tke_budget_rhs
          procedure, private :: get_TT_budget_rhs
          procedure, private :: get_wT_budget_rhs
          procedure, private :: compute_unsteady_terms
          procedure, private :: turbulent_transport
          procedure, private :: gradient_production
          procedure, private :: convective_transport
          procedure, private :: pseudo_dissipation
          procedure, private :: molecular_transport
          procedure, private :: SGS_transport
          procedure, private :: SGS_dissipation
          procedure, private :: covariance
          procedure, private :: mean
          procedure, private :: ddz
          procedure, private :: d2dz2
    end type
    contains

      subroutine init(this,inputfile,sim)
        class(stats_xy), intent(inout), target :: this
        character(len=*), intent(in) :: inputfile
        type(igrid), intent(in), target :: sim
        integer :: tid_start = 1000000, compute_freq = 1000000, dump_freq = 10000000
        character(len=clen) :: outputdir
        integer :: ioUnit, ierr, nstore, n
        integer, parameter :: nterms      = 86
        integer, parameter :: nterms_sca  = 41
        integer, parameter :: n_ddt_terms = 6
        real, parameter :: zero = 0.d0
        logical :: do_stats = .false.

        namelist /STATS_XY/ tid_start, compute_freq, dump_freq, outputdir, &
          do_stats

        ioUnit = 534
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=STATS_XY)
        close(ioUnit)

        if (do_stats) then
            this%tid_start    = tid_start
            this%compute_freq = compute_freq
            this%dump_freq    = dump_freq
            this%outputdir    = outputdir
            this%do_stats     = do_stats

            this%sim => sim
            this%avgFact = 1.d0/(real(sim%nx,rkind)*real(sim%ny,rkind))

            this%nscalars = 0
            if (sim%useScalars)   this%nscalars = this%sim%n_scalars 
            if (sim%isStratified) this%nscalars = this%nscalars + 1

            ! How many snapshots should be stored before dumping (don't want to have
            ! to dump to disk so often for performance reasons)
            nstore = dump_freq/compute_freq
            this%nwrite = 0
            
            ! User safeguards
            call assert(this%compute_freq > 3,'Compute frequency must be set '//&
              'to greater than 3 -- stats_xy.F90')

            ! Allocate memory and link pointers
            allocate(this%stats(sim%nz,nterms,nstore))
            if (this%nscalars > 0) then
                allocate(this%stats_sca(sim%nz,nterms_sca,nstore,this%nscalars))
                allocate(this%q3(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3)))
            end if
            allocate(this%time(nstore))
            allocate(this%tid(nstore))
            allocate(this%zbuff(sim%nz,4))
            allocate(this%zbuff_(sim%nz,2))
            allocate(this%qold(sim%nz,n_ddt_terms), this%qnew(sim%nz,n_ddt_terms))
            allocate(this%tau11(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3)))
            allocate(this%tau12(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3)))
            allocate(this%tau13(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3)))
            allocate(this%tau22(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3)))
            allocate(this%tau23(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3)))
            allocate(this%tau33(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3)))

            this%tidx = 0
            call this%link_pointers(this%tidx+1,min(1,this%nscalars))
            call this%write_labels()

            call message("==========================================")
            call message("stats_xy instance initialized successfully")
            call message(" ")
        end if
      end subroutine

      subroutine destroy(this)
          class(stats_xy), intent(inout) :: this

          nullify(this%tke_budget,this%TT_budget,this%dTdz,this%tke,this%fifi,&
            this%meanU,this%meanV,this%meanW,this%meanP,this%meanT,&
            this%uu,this%vv,this%ww,this%uv,this%uw,this%vw,this%uT,&
            this%vT,this%wT,this%TT,this%meanFx,this%meanFy,this%meanFz)
          deallocate(this%stats)
          nullify(this%sim)
          if (allocated(this%time)) deallocate(this%time)
          if (allocated(this%tid))  deallocate(this%tid)
      end subroutine

      subroutine compute_stats(this)
          class(stats_xy), intent(inout) :: this
          integer :: n

          if ( (this%do_stats)                             .and. &
               (this%sim%step >= this%tid_start          )) then

             ! Step 1: Get budget RHS terms and all other budget terms at t^n
             if (mod(this%sim%step,this%compute_freq) == 0) then
                 this%tidx = this%tidx + 1
                 call this%link_pointers(this%tidx,min(1,this%nscalars))

                 ! Get tauij_SGS
                 call this%sim%sgsModel%populate_tauij_E_to_C()
                 call this%sim%sgsModel%get_tauijC(this%tau11,this%tau12,this%tau13,&
                   this%tau22,this%tau23,this%tau33)

                 ! Get mean profiles
                 call this%mean(this%sim%u       ,this%meanU)
                 call this%mean(this%sim%v       ,this%meanV)
                 call this%mean(this%sim%wC      ,this%meanW)
                 ! TODO: This is not the true mean pressure. Need to account for
                 ! meanT component
                 call this%mean(this%sim%pressure,this%meanP) 

                 !   > mean force components
                 if (this%sim%localizedForceLayer == 2) then
                     this%sim%rbuffxC(:,:,:,1) = this%sim%spectForceLayer%ampFact_x*this%sim%spectForceLayer%fx
                     this%sim%rbuffxC(:,:,:,2) = this%sim%spectForceLayer%ampFact_y*this%sim%spectForceLayer%fy
                     this%sim%rbuffxE(:,:,:,1) = this%sim%spectForceLayer%ampFact_z*this%sim%spectForceLayer%fz
                     call this%sim%spectForceLayer%interpE2C(&
                       this%sim%rbuffxE(:,:,:,1), this%sim%rbuffxC(:,:,:,3), &
                       this%sim%rbuffyC(:,:,:,1), this%sim%rbuffzC(:,:,:,1), & 
                       this%sim%rbuffyE(:,:,:,1), this%sim%rbuffzE(:,:,:,1))

                     call this%mean(this%sim%rbuffxC(:,:,:,1), this%meanFx)
                     call this%mean(this%sim%rbuffxC(:,:,:,2), this%meanFy)
                     call this%mean(this%sim%rbuffxC(:,:,:,3), this%meanFz)

                     ! force variance
                     this%sim%rbuffxC(:,:,:,2) = &
                       this%sim%spectForceLayer%ampFact_x*this%sim%spectForceLayer%ampFact_x*&
                       this%sim%spectForceLayer%fx*this%sim%spectForceLayer%fx + & 
                       this%sim%spectForceLayer%ampFact_y*this%sim%spectForceLayer%ampFact_y*&
                       this%sim%spectForceLayer%fy*this%sim%spectForceLayer%fy + & 
                       this%sim%spectForceLayer%ampFact_z*this%sim%spectForceLayer%ampFact_z*&
                       this%sim%rbuffxC(:,:,:,1)*this%sim%rbuffxC(:,:,:,1)
                     call this%mean(this%sim%rbuffxC(:,:,:,2),this%fifi)
                     this%fifi = this%fifi - (this%meanFx*this%meanFx + &
                       this%meanFy*this%meanFy + this%meanFz*this%meanFz)
                 else
                     this%meanFx = 0.d0
                     this%meanFy = 0.d0
                     this%meanFz = 0.d0
                     this%fifi   = 0.d0
                 end if

                 ! Compute MKE
                 this%mke = 0.5d0*(this%meanU*this%meanU + this%meanV*this%meanV + &
                   this%meanW*this%meanW)

                 ! Reynolds stresses
                 !   > Normal stresses
                 call this%covariance(this%sim%u,  this%sim%u,  this%uu)
                 call this%covariance(this%sim%v,  this%sim%v,  this%vv)
                 call this%covariance(this%sim%wC, this%sim%wC, this%ww)

                 !   > Cross stresses
                 call this%covariance(this%sim%u,  this%sim%v,  this%uv)
                 call this%covariance(this%sim%u,  this%sim%wC, this%uw)
                 call this%covariance(this%sim%v,  this%sim%wC, this%vw)

                 ! Compute TKE
                 this%tke = 0.5d0*(this%uu + this%vv + this%ww)

                 ! Compute tau_SGS stats
                 call this%covariance(this%tau11,this%tau11,this%tauij_var(:,1))
                 call this%covariance(this%tau12,this%tau12,this%tauij_var(:,2))
                 call this%covariance(this%tau13,this%tau13,this%tauij_var(:,3))

                 call this%covariance(this%tau22,this%tau22,this%tauij_var(:,4))
                 call this%covariance(this%tau23,this%tau23,this%tauij_var(:,5))
                 call this%covariance(this%tau33,this%tau33,this%tauij_var(:,6))

                 ! Add time ID and simulation time to their respective vectors
                 this%tid(this%tidx) = this%sim%step
                 this%time(this%tidx) = this%sim%tsim

                 ! Assemble budgets
                 call this%get_tke_budget_rhs(this%tke_budget,this%R11_budget,&
                   this%R22_budget,this%R33_budget,&
                   this%tau11,this%tau12,this%tau13,this%tau22,this%tau23,this%tau33)

                 ! TKE fluxes
                 this%sim%rbuffxC(:,:,:,1) = this%sim%u*this%sim%u + &
                   this%sim%v*this%sim%v + this%sim%wC*this%sim%wC
                 call this%covariance(this%sim%u ,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,1))
                 call this%covariance(this%sim%v ,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,2))
                 call this%covariance(this%sim%wC,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,3))

                 this%tke_flux(:,1) = 0.5d0*this%zbuff(:,1) - &
                   (this%meanU*this%uu + this%meanV*this%uv + this%meanW*this%uw)
                 this%tke_flux(:,2) = 0.5d0*this%zbuff(:,2) - &
                   (this%meanU*this%uv + this%meanV*this%vv + this%meanW*this%vw)
                 this%tke_flux(:,3) = 0.5d0*this%zbuff(:,3) - &
                   (this%meanU*this%uw + this%meanV*this%vw + this%meanW*this%ww)


                 do n = 1,this%nscalars
                     call this%link_pointers(this%tidx,n)
                     
                     ! TODO: Mean scalar profile
                     call this%mean(this%T, this%meanT)

                     ! TODO: Mean scalar gradient
                     call this%mean(this%dTdzC,this%dTdz)

                     ! Scalar variance and fluxes
                     call this%covariance(this%T,      this%T, this%TT)
                     call this%covariance(this%sim%u,  this%T, this%uT)
                     call this%covariance(this%sim%v,  this%T, this%vT)
                     call this%covariance(this%sim%wC, this%T, this%wT)

                     ! Scalar graident variances
                     call this%covariance(this%dTdxC,this%dTdxC,this%dTdxj_var(:,1))
                     call this%covariance(this%dTdyC,this%dTdyC,this%dTdxj_var(:,2))
                     call this%covariance(this%dTdzC,this%dTdzC,this%dTdxj_var(:,3))

                     ! SGS heat flux moments
                     call this%covariance(this%q1,this%q1,this%qj_var(:,1))
                     call this%covariance(this%q2,this%q2,this%qj_var(:,2))
                     call this%covariance(this%q3,this%q3,this%qj_var(:,3))
                     call this%mean(this%q3,this%meanq3)

                     call this%get_TT_budget_rhs()
                     call this%get_wT_budget_rhs()

                     this%sim%rbuffxC(:,:,:,1) = this%T*this%T
                     call this%covariance(this%sim%u ,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,1))
                     call this%covariance(this%sim%v ,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,2))
                     call this%covariance(this%sim%wC,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,3))

                     this%sca_var_flux(:,1) = this%zbuff(:,1) - 2.d0*this%meanT*this%uT
                     this%sca_var_flux(:,2) = this%zbuff(:,2) - 2.d0*this%meanT*this%vT
                     this%sca_var_flux(:,3) = this%zbuff(:,3) - 2.d0*this%meanT*this%wT
                 end do

                 ! Increment nwrite. This is to handle the edge-case where we
                 ! started the simulation at a tid such that the first dump call is
                 ! before we've filled the stats buffers, so we only want to write
                 ! data that has been computed.
                 this%nwrite = this%nwrite + 1

             end if

             ! Step 2: compute unsteady terms
             call this%compute_unsteady_terms()

             ! Dump stats at the dump frequency depending on if dt is constant
             ! or dynamic. If dynamic, dump the data because we have ddt terms.
             ! If it is constant, we use a second order central difference for
             ! ddt and so we need to wait until we have the ddt term at the
             ! next time step to dump the data
             if (mod(this%sim%step,this%dump_freq) == 0) then
                 if (this%sim%CFL >= 0.d0) then
                     call this%dump_stats()
                     this%nwrite = 0
                     this%tidx = 0
                 end if
             elseif (mod(this%sim%step-1,this%dump_freq) == 0) then
                 if (this%sim%CFL <  0.d0) then
                     call this%dump_stats()
                     this%nwrite = 0
                     this%tidx = 0
                 end if
             end if
          end if
      end subroutine

      subroutine get_tke_budget_rhs(this,tke_budget,R11_budget,R22_budget,R33_budget,&
          tau_11,tau_12,tau_13,tau_22,tau_23,tau_33)
          class(stats_xy), intent(inout) :: this
          real(rkind), dimension(:,:), intent(inout) :: tke_budget, &
            R11_budget, R22_budget, R33_budget
          real(rkind), dimension(:,:,:), intent(in) :: tau_11, tau_12, tau_13, &
                                                       tau_22, tau_23, tau_33
          real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, &
                                                    dvdx, dvdy, dvdz, &
                                                    dwdx, dwdy, dwdz
          real(rkind) :: TwobyRe, OnebyRe

          tke_budget = 0.d0
          OnebyRe = 1.d0/this%sim%Re
          TwobyRe = 2.d0*OnebyRe

          associate( &
              dudx => this%sim%duidxjC(:,:,:,1), dudy => this%sim%duidxjC(:,:,:,2), dudz => this%sim%duidxjC(:,:,:,3), &
              dvdx => this%sim%duidxjC(:,:,:,4), dvdy => this%sim%duidxjC(:,:,:,5), dvdz => this%sim%duidxjC(:,:,:,6), &
              dwdx => this%sim%duidxjC(:,:,:,7), dwdy => this%sim%duidxjC(:,:,:,8), dwdz => this%sim%duidxjC(:,:,:,9))

              ! Shear production
              call this%gradient_production(dudx,dudy,dudz,this%uu,this%uv,this%uw,R11_budget(:,2))
              call this%gradient_production(dvdx,dvdy,dvdz,this%uv,this%vv,this%vw,R22_budget(:,2))
              call this%gradient_production(dwdx,dwdy,dwdz,this%uw,this%vw,this%ww,R33_budget(:,2))
              tke_budget(:,2) = R11_budget(:,2) + R22_budget(:,2) + R33_budget(:,3)
              R11_budget(:,2) = 2.d0*R11_budget(:,2)
              R22_budget(:,2) = 2.d0*R22_budget(:,2)
              R33_budget(:,2) = 2.d0*R33_budget(:,2)

              ! Force production
              call this%sim%spectForceLayer%interpE2C(&
                this%sim%spectForceLayer%fz, this%sim%rbuffxC(:,:,:,1), &
                this%sim%rbuffyC(:,:,:,1), this%sim%rbuffzC(:,:,:,1), & 
                this%sim%rbuffyE(:,:,:,1), this%sim%rbuffzE(:,:,:,1))
              call this%covariance(this%sim%spectForceLayer%fx, this%sim%u,  this%zbuff(:,1))
              call this%covariance(this%sim%spectForceLayer%fy, this%sim%v,  this%zbuff(:,2))
              call this%covariance(this%sim%rbuffxC(:,:,:,1)  , this%sim%wC, this%zbuff(:,3))
              tke_budget(:,3) = this%sim%spectForceLayer%ampFact_x*this%zbuff(:,1) + &
                this%sim%spectForceLayer%ampFact_y*this%zbuff(:,2) + &
                this%sim%spectForceLayer%ampFact_z*this%zbuff(:,3)
              R11_budget(:,3) = 2.d0*this%sim%spectForceLayer%ampFact_x*this%zbuff(:,1)
              R22_budget(:,3) = 2.d0*this%sim%spectForceLayer%ampFact_y*this%zbuff(:,2)
              R33_budget(:,3) = 2.d0*this%sim%spectForceLayer%ampFact_z*this%zbuff(:,3)

              ! Convective transport
              call this%convective_transport(this%uu,R11_budget(:,4))
              call this%convective_transport(this%vv,R22_budget(:,4))
              call this%convective_transport(this%ww,R33_budget(:,4))
              tke_budget(:,4) = 0.5d0*(R11_budget(:,4) + R22_budget(:,4) + R33_budget(:,4))

              ! Turbulent transport
              call this%turbulent_transport(this%sim%u,this%sim%u,&
                dudx,dudy,dudz,dudx,dudy,dudz,R11_budget(:,4),R11_budget(:,2),R11_budget(:,5))
              call this%turbulent_transport(this%sim%v,this%sim%v,&
                dvdx,dvdy,dvdz,dvdx,dvdy,dvdz,R22_budget(:,4),R22_budget(:,2),R22_budget(:,5))
              call this%turbulent_transport(this%sim%wC,this%sim%wC,&
                dwdx,dwdy,dwdz,dwdx,dwdy,dwdz,R33_budget(:,4),R33_budget(:,2),R33_budget(:,5))

              tke_budget(:,5) = 0.5d0*(R11_budget(:,5) + R22_budget(:,5) + R33_budget(:,5))

              ! Pressure transport
              call this%covariance(this%sim%wC, this%sim%pressure, this%zbuff(:,1))
              this%zbuff(:,1) = -1.d0*this%zbuff(:,1)
              call this%ddz(this%zbuff(:,1), tke_budget(:,6), 0, 0)
              R11_budget(:,6) = 0.d0
              R22_budget(:,6) = 0.d0
              R33_budget(:,6) = 2.d0*tke_budget(:,6)

              ! Viscous transport: ddz(<ui*Si3>)
              this%sim%rbuffxC(:,:,:,1) = 0.5d0*(dudz + dwdx)
              this%sim%rbuffxC(:,:,:,2) = 0.5d0*(dvdz + dwdy)
              call this%covariance(this%sim%u ,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,1))
              call this%covariance(this%sim%v ,this%sim%rbuffxC(:,:,:,2),this%zbuff(:,2))
              call this%covariance(this%sim%wC,dwdz                     ,this%zbuff(:,3))
              this%zbuff(:,4) = TwobyRe*(this%zbuff(:,1) + this%zbuff(:,2) + this%zbuff(:,3))
              call this%ddz(this%zbuff(:,4),tke_budget(:,7),0,0)

              call this%molecular_transport(this%sim%u ,dudz,OnebyRe,R11_budget(:,7))
              call this%molecular_transport(this%sim%v ,dvdz,OnebyRe,R22_budget(:,7))
              call this%molecular_transport(this%sim%wC,dwdz,OnebyRe,R33_budget(:,7))
              tke_budget(:,13) = R11_budget(:,7) + R22_budget(:,7) + R33_budget(:,7)
              R11_budget(:,7) = 2.d0*R11_budget(:,7)
              R22_budget(:,7) = 2.d0*R22_budget(:,7)
              R33_budget(:,7) = 2.d0*R33_budget(:,7)

              ! Viscous dissipation
              !   > NOTE: this is positive. Will need to plot
              !     negative of this for budget visualization
              ! TODO: There is bug in viscous dissipation
              this%sim%rbuffxC(:,:,:,3) = 0.5d0*(dudy + dvdx) ! S12
              call this%covariance(dudx,dudx,this%zbuff(:,1)) ! <S11'*S11'>
              call this%covariance(dvdy,dvdy,this%zbuff(:,2)) ! <S22'*S22'>
              call this%covariance(dwdz,dwdz,this%zbuff(:,3)) ! <S33'*S33'>
              tke_budget(:,8) = 0.d0
              this%zbuff(:,4) = this%zbuff(:,1) + this%zbuff(:,2) + this%zbuff(:,3)
              
              call this%covariance(this%sim%rbuffxC(:,:,:,3),&
                this%sim%rbuffxC(:,:,:,3),this%zbuff(:,3))    ! 0.5*(<S12'*S12'> + <S21'*S21'>)
              call this%covariance(this%sim%rbuffxC(:,:,:,1),&
                this%sim%rbuffxC(:,:,:,1),this%zbuff(:,1))    ! 0.5*(<S13'*S13'> + <S31'*S31'>)
              call this%covariance(this%sim%rbuffxC(:,:,:,2),&
                this%sim%rbuffxC(:,:,:,2),this%zbuff(:,2))    ! 0.5*(<S23'*S23'> + <S32'*S32'>)
              tke_budget(:,8) = TwobyRe*(this%zbuff(:,4) + &
                2.d0*(this%zbuff(:,3) + this%zbuff(:,1) + this%zbuff(:,2)))

              call this%pseudo_dissipation(dudx,dudy,dudz,dudx,dudy,dudz,TwobyRe,R11_budget(:,8))
              call this%pseudo_dissipation(dvdx,dvdy,dvdz,dvdx,dvdy,dvdz,TwobyRe,R22_budget(:,8))
              call this%pseudo_dissipation(dwdx,dwdy,dwdz,dwdx,dwdy,dwdz,TwobyRe,R33_budget(:,8))
              tke_budget(:,12) = 0.5d0*(R11_budget(:,8) + R22_budget(:,8) + R33_budget(:,8))

              ! Compute velocity gradient variances
              call this%covariance(dudx,dudx,this%duidxj_var(:,1))
              call this%covariance(dudy,dudy,this%duidxj_var(:,2))
              call this%covariance(dudz,dudz,this%duidxj_var(:,3))

              call this%covariance(dvdx,dvdx,this%duidxj_var(:,4))
              call this%covariance(dvdy,dvdy,this%duidxj_var(:,5))
              call this%covariance(dvdz,dvdz,this%duidxj_var(:,6))

              call this%covariance(dwdx,dwdx,this%duidxj_var(:,7))
              call this%covariance(dwdy,dwdy,this%duidxj_var(:,8))
              call this%covariance(dwdz,dwdz,this%duidxj_var(:,9))

              call this%covariance(this%sim%rbuffxC(:,:,:,3),this%sim%rbuffxC(:,:,:,3),this%S12_var)
              call this%covariance(this%sim%rbuffxC(:,:,:,1),this%sim%rbuffxC(:,:,:,1),this%S13_var)
              call this%covariance(this%sim%rbuffxC(:,:,:,2),this%sim%rbuffxC(:,:,:,2),this%S23_var)

              ! SGS transport
              call this%SGS_transport(this%sim%u ,this%tau13,R11_budget(:,9))
              call this%SGS_transport(this%sim%v ,this%tau23,R22_budget(:,9))
              call this%SGS_transport(this%sim%wC,this%tau33,R33_budget(:,9))
              tke_budget(:,9) = R11_budget(:,9) + R22_budget(:,9) + R33_budget(:,9)
              R11_budget(:,9) = 2.d0*R11_budget(:,9)
              R22_budget(:,9) = 2.d0*R22_budget(:,9)
              R33_budget(:,9) = 2.d0*R33_budget(:,9)

              ! SGS dissipation
              !   > NOTE: consistent with viscous dissipation, the sign on this
              !   is opposite of what needs to be plotted so that epsilon_total
              !   = epsilon + epsilon_SGS and -epsilon_total is on RHS in budget
              call this%SGS_dissipation(dudx,dudy,dudz,this%tau11,this%tau12,this%tau13,R11_budget(:,10))
              call this%SGS_dissipation(dvdx,dvdy,dvdz,this%tau12,this%tau22,this%tau23,R22_budget(:,10))
              call this%SGS_dissipation(dwdx,dwdy,dwdz,this%tau13,this%tau23,this%tau33,R33_budget(:,10))
              tke_budget(:,10) = R11_budget(:,10) + R22_budget(:,10) + R33_budget(:,10)
              R11_budget(:,10) = 2.d0*R11_budget(:,10)
              R22_budget(:,10) = 2.d0*R22_budget(:,10)
              R33_budget(:,10) = 2.d0*R33_budget(:,10)
              
              ! Buoyancy transfer
              if (this%sim%isStratified) then
                  call this%covariance(this%sim%wC, this%sim%T, tke_budget(:,11))
                  tke_budget(:,11) = -1.d0/this%sim%buoyancyFact*tke_budget(:,11)
                  R11_budget(:,11) = 0.d0
                  R22_budget(:,11) = 0.d0
                  R33_budget(:,11) = 2.d0*tke_budget(:,11)
              else
                  tke_budget(:,11) = 0.d0
                  R11_budget(:,11) = 0.d0
                  R22_budget(:,11) = 0.d0
                  R33_budget(:,11) = 0.d0 
              end if

              ! Pressure-strain correlations
              call this%covariance(this%sim%pressure,dudx,R11_budget(:,12))
              call this%covariance(this%sim%pressure,dvdy,R22_budget(:,12))
              call this%covariance(this%sim%pressure,dwdz,R33_budget(:,12))
              R11_budget(:,12) = 2.d0*R11_budget(:,12)
              R22_budget(:,12) = 2.d0*R22_budget(:,12)
              R33_budget(:,12) = 2.d0*R33_budget(:,12)

          end associate

      end subroutine

      subroutine get_TT_budget_rhs(this)
          class(stats_xy), intent(inout) :: this
          real(rkind) :: OnebyRePr
          
          this%TT_budget = 0.d0
          OnebyRePr = 1.d0/(this%sim%Re*this%sim%PrandtlFluid)

          ! Gradient production
          call this%gradient_production(this%dTdxC,this%dTdyC,this%dTdzC,&
            this%uT,this%vT,this%wT,this%TT_budget(:,2))
          this%TT_budget(:,2) = 2.d0*this%TT_budget(:,2)

          ! Convective transport
          call this%convective_transport(this%TT,this%TT_budget(:,3))

          ! Turbulent transport
          call this%turbulent_transport(this%T,this%T,this%dTdxC,this%dTdyC,this%dTdzC,&
            this%dTdxC,this%dTdyC,this%dTdzC,this%TT_budget(:,3), this%TT_budget(:,2), &
            this%TT_budget(:,4))

          ! Molecular transport
          call this%molecular_transport(this%T,this%dTdzC,OnebyRePr,this%TT_budget(:,5))
          this%TT_budget(:,5) = 2.d0*this%TT_budget(:,5)

          ! SGS transport
          call this%SGS_transport(this%T,this%q3,this%TT_budget(:,6))
          this%TT_budget(:,6) = 2.d0*this%TT_budget(:,6)

          ! Molecular destruction
          call this%pseudo_dissipation(this%dTdxC,this%dTdyC,this%dTdzC,this%dTdxC,this%dTdyC,this%dTdzC,OnebyRePr,this%TT_budget(:,7))
          this%TT_budget(:,7) = 2.d0*this%TT_budget(:,7)

          ! SGS destruction
          call this%SGS_dissipation(this%dTdxC,this%dTdyC,this%dTdzC,this%q1,this%q2,this%q3,this%TT_budget(:,8))
          this%TT_budget(:,8) = 2.d0*this%TT_budget(:,8)

      end subroutine

      subroutine get_wT_budget_rhs(this)
          class(stats_xy), intent(inout) :: this
          real(rkind) :: OnebyRePr, OnebyRe, OnebyFr2, viscosity
          real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, &
            dvdy, dvdz, dwdx, dwdy, dwdz

          this%wT_budget = 0.d0
          OnebyRePr = 1.d0/(this%sim%Re*this%sim%PrandtlFluid)
          OnebyRe = 1.d0/this%sim%Re
          OnebyFr2 = this%sim%buoyancyFact

          associate( &
              dudx => this%sim%duidxjC(:,:,:,1), dudy => this%sim%duidxjC(:,:,:,2), dudz => this%sim%duidxjC(:,:,:,3), &
              dvdx => this%sim%duidxjC(:,:,:,4), dvdy => this%sim%duidxjC(:,:,:,5), dvdz => this%sim%duidxjC(:,:,:,6), &
              dwdx => this%sim%duidxjC(:,:,:,7), dwdy => this%sim%duidxjC(:,:,:,8), dwdz => this%sim%duidxjC(:,:,:,9))

              ! Scalar gradient production
              call this%gradient_production(this%dTdxC,this%dTdyC,this%dTdzC,this%uw,this%vw,this%ww,this%wT_budget(:,2))

              ! Shear gradient production
              call this%gradient_production(dwdx,dwdy,dwdz,this%uT,this%vT,this%wT,this%wT_budget(:,3))

              ! Force production
              this%sim%rbuffxE(:,:,:,1) = this%sim%spectForceLayer%ampFact_z*this%sim%spectForceLayer%fz
              call this%sim%spectForceLayer%interpE2C(&
                this%sim%rbuffxE(:,:,:,1), this%sim%rbuffxC(:,:,:,1), &
                this%sim%rbuffyC(:,:,:,1), this%sim%rbuffzC(:,:,:,1), & 
                this%sim%rbuffyE(:,:,:,1), this%sim%rbuffzE(:,:,:,1))
              call this%covariance(this%T,this%sim%rbuffxC(:,:,:,1),this%wT_budget(:,4))

              ! Buoyancy production
              this%wT_budget(:,5) = OnebyFr2*this%TT

              ! Convective transport
              call this%convective_transport(this%wT,this%wT_budget(:,6))

              ! Turbulent transport
              this%wT_budget(:,2) = this%wT_budget(:,2) + this%wT_budget(:,3)
              call this%turbulent_transport(this%sim%wC,this%T,dwdx,dwdy,dwdz,&
                this%dTdxC,this%dTdyC,this%dTdzC,this%wT_budget(:,6), this%wT_budget(:,2), &
                this%wT_budget(:,7))
              this%wT_budget(:,2) = this%wT_budget(:,2) - this%wT_budget(:,3)

              ! Pressure transport
              call this%covariance(this%T,this%sim%pressure,this%zbuff(:,1))
              call this%ddz(this%zbuff(:,1),this%wT_budget(:,8),0,0)
              this%wT_budget(:,8) = -1.d0*this%wT_budget(:,8)

              ! Molecular transport
              call this%molecular_transport(this%sim%wC,this%dTdzC,OnebyRePr,this%wT_budget(:,9 ))
              call this%molecular_transport(this%T     ,dwdz      ,OnebyRe  ,this%wT_budget(:,10))

              ! SGS transport
              call this%SGS_transport(this%sim%wC,this%q3   ,this%wT_budget(:,11))
              call this%SGS_transport(this%T     ,this%tau33,this%wT_budget(:,12))

              ! Pressure scrambling
              call this%covariance(this%sim%pressure, this%dTdzC, this%wT_budget(:,13))

              ! Molecular destruction
              call this%pseudo_dissipation(dwdx,dwdy,dwdz,this%dTdxC,this%dTdyC,this%dTdzC,OnebyRePr,this%wT_budget(:,14))
              call this%pseudo_dissipation(dwdx,dwdy,dwdz,this%dTdxC,this%dTdyC,this%dTdzC,OnebyRe  ,this%wT_budget(:,15))

              ! SGS destruction
              call this%SGS_dissipation(dwdx,dwdy,dwdz,this%q1,this%q2,this%q3,this%wT_budget(:,16))
              call this%SGS_dissipation(this%dTdxC,this%dTdyC,this%dTdzC,this%tau13,this%tau23,this%tau33,this%wT_budget(:,17))
          end associate
      end subroutine

      subroutine link_pointers(this,tidx,sca)
        class(stats_xy), intent(inout), target :: this
        integer, intent(in) :: tidx, sca
        integer :: id

        id = 13
        this%tke_budget => this%stats(:,1:id    ,tidx); id = id + 1
        this%tke        => this%stats(:,id      ,tidx); id = id + 1
        this%fifi       => this%stats(:,id      ,tidx); id = id + 1
        this%meanU      => this%stats(:,id      ,tidx); id = id + 1
        this%meanV      => this%stats(:,id      ,tidx); id = id + 1
        this%meanW      => this%stats(:,id      ,tidx); id = id + 1
        this%meanP      => this%stats(:,id      ,tidx); id = id + 1
        this%uu         => this%stats(:,id      ,tidx); id = id + 1
        this%vv         => this%stats(:,id      ,tidx); id = id + 1
        this%ww         => this%stats(:,id      ,tidx); id = id + 1
        this%uv         => this%stats(:,id      ,tidx); id = id + 1
        this%uw         => this%stats(:,id      ,tidx); id = id + 1
        this%vw         => this%stats(:,id      ,tidx); id = id + 1
        this%meanFx     => this%stats(:,id      ,tidx); id = id + 1
        this%meanFy     => this%stats(:,id      ,tidx); id = id + 1
        this%meanFz     => this%stats(:,id      ,tidx); id = id + 1
        this%mke        => this%stats(:,id      ,tidx); id = id + 1
        this%R11_budget => this%stats(:,id:id+11,tidx); id = id + 12
        this%R22_budget => this%stats(:,id:id+11,tidx); id = id + 12
        this%R33_budget => this%stats(:,id:id+11,tidx); id = id + 12
        this%duidxj_var => this%stats(:,id:id+8 ,tidx); id = id + 9
        this%S12_var    => this%stats(:,id      ,tidx); id = id + 1
        this%S13_var    => this%stats(:,id      ,tidx); id = id + 1
        this%S23_var    => this%stats(:,id      ,tidx); id = id + 1
        this%tauij_var  => this%stats(:,id:id+5 ,tidx); id = id + 6
        this%tke_flux   => this%stats(:,id:id+2 ,tidx); id = id + 3

        if (this%nscalars > 0) then
            id = 1
            this%TT_budget    => this%stats_sca(:,id:id+7 ,tidx,sca); id = id + 8
            this%dTdz         => this%stats_sca(:,id      ,tidx,sca); id = id + 1
            this%meanT        => this%stats_sca(:,id      ,tidx,sca); id = id + 1
            this%uT           => this%stats_sca(:,id      ,tidx,sca); id = id + 1
            this%vT           => this%stats_sca(:,id      ,tidx,sca); id = id + 1
            this%wT           => this%stats_sca(:,id      ,tidx,sca); id = id + 1
            this%TT           => this%stats_sca(:,id      ,tidx,sca); id = id + 1
            this%wT_budget    => this%stats_sca(:,id:id+16,tidx,sca); id = id + 17
            this%dTdxj_var    => this%stats_sca(:,id:id+2 ,tidx,sca); id = id + 3
            this%qj_var       => this%stats_sca(:,id:id+2 ,tidx,sca); id = id + 3
            this%meanq3       => this%stats_sca(:,id      ,tidx,sca); id = id + 1
            this%sca_var_flux => this%stats_sca(:,id:id+2 ,tidx,sca); id = id + 3 

            if (sca == 1) then
                this%T    => this%sim%T
                this%dTdxC => this%sim%dTdxC
                this%dTdyC => this%sim%dTdyC
                this%dTdzC => this%sim%dTdzC
                this%q1 => this%sim%q1_T
                this%q2 => this%sim%q2_T
                call this%sim%spectForceLayer%interpE2C(this%sim%q3_T, this%q3, &
                  this%sim%rbuffyC(:,:,:,1), this%sim%rbuffzC(:,:,:,1), & 
                  this%sim%rbuffyE(:,:,:,1), this%sim%rbuffzE(:,:,:,1))
            else
                this%T    => this%sim%scalars(sca-1)%F
                this%dTdxC => this%sim%scalars(sca-1)%dFdxC
                this%dTdyC => this%sim%scalars(sca-1)%dFdyC
                this%dTdzC => this%sim%scalars(sca-1)%dFdzC
                this%q1 => this%sim%scalars(sca-1)%q1
                this%q2 => this%sim%scalars(sca-1)%q2
                call this%sim%spectForceLayer%interpE2C(this%sim%scalars(sca-1)%q3, this%q3, &
                  this%sim%rbuffyC(:,:,:,1), this%sim%rbuffzC(:,:,:,1), & 
                  this%sim%rbuffyE(:,:,:,1), this%sim%rbuffzE(:,:,:,1))
            end if
        end if

      end subroutine
      
      subroutine dump_stats(this)
          class(stats_xy), intent(inout) :: this
          integer :: n, sca
          character(len=clen) :: tempname, fname

          if (nrank == 0) then
              do n = 1,this%nwrite
                  write(tempname,"(A3,I2.2,A11,I6.6,A4)")"Run",this%sim%runID,"_stats_xy_t",this%tid(n),".out"
                  fname = trim(this%outputdir)//"/"//trim(tempname)
                  call write_2d_ascii(this%stats(:,:,n),trim(fname))
                  do sca = 1,this%nscalars
                      write(tempname,"(A3,I2.2,A13,I2.2,A2,I6.6,A4)")&
                        "Run",this%sim%runID,"_stats_xy_sca",sca,"_t",this%tid(n),".out"
                      fname = trim(this%outputdir)//"/"//trim(tempname)
                      call write_2d_ascii(this%stats_sca(:,:,n,sca),trim(fname))
                  end do
              end do
          end if
          call this%write_time_info()
      end subroutine

      subroutine write_labels(this)
          class(stats_xy), intent(inout) :: this
          integer :: ierr, fid, id
          character(len=clen) :: tempname, fname
          logical :: exists

          fid = 12
          
          write(tempname,"(A3,I2.2,A29)")"Run",this%sim%runID,"_momentum_stats_index_key.out"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)

          if (nrank == 0) then
              id = 1
              open(fid, file=fname)!, status="new", action="write")
              write(fid,*) "TKE budget terms"
              write(fid,"(A4,I2,A15)") "  > ",id,". Unsteady term";            id = id + 1
              write(fid,"(A4,I2,A18)") "  > ",id,". Shear production";         id = id + 1
              write(fid,"(A4,I2,A18)") "  > ",id,". Force production";         id = id + 1
              write(fid,"(A4,I2,A22)") "  > ",id,". Convective transport";     id = id + 1
              write(fid,"(A4,I2,A21)") "  > ",id,". Turbulent transport";      id = id + 1
              write(fid,"(A4,I2,A20)") "  > ",id,". Pressure transport";       id = id + 1
              write(fid,"(A4,I2,A19)") "  > ",id,". Viscous transport";        id = id + 1
              write(fid,"(A4,I2,A21)") "  > ",id,". Viscous dissipation";      id = id + 1
              write(fid,"(A4,I2,A15)") "  > ",id,". SGS transport";            id = id + 1
              write(fid,"(A4,I2,A17)") "  > ",id,". SGS dissipation";          id = id + 1
              write(fid,"(A4,I2,A19)") "  > ",id,". Buoyancy transfer";        id = id + 1
              write(fid,"(A4,I2,A20)") "  > ",id,". Pseudo dissipation";       id = id + 1
              write(fid,"(A4,I2,A26)") "  > ",id,". Pseudo viscous transport"; id = id + 1
              write(fid,*) " "
              write(fid,*) "Misc terms:"
              write(fid,"(A4,I2,A33)") "  > ",id,". Turbulence kinetic energy (TKE)"; id = id + 1
              write(fid,"(A4,I2,A16)") "  > ",id,". Force variance";                  id = id + 1
              write(fid,"(A4,I2,A17)") "  > ",id,". Mean u-velocity";                 id = id + 1
              write(fid,"(A4,I2,A17)") "  > ",id,". Mean v-velocity";                 id = id + 1 
              write(fid,"(A4,I2,A17)") "  > ",id,". Mean w-velocity";                 id = id + 1 
              write(fid,"(A4,I2,A15)") "  > ",id,". Mean pressure";                   id = id + 1 
              write(fid,"(A4,I2,A8)" ) "  > ",id,". <u'u'>";                          id = id + 1 
              write(fid,"(A4,I2,A8)" ) "  > ",id,". <v'v'>";                          id = id + 1 
              write(fid,"(A4,I2,A8)" ) "  > ",id,". <w'w'>";                          id = id + 1 
              write(fid,"(A4,I2,A8)" ) "  > ",id,". <u'v'>";                          id = id + 1 
              write(fid,"(A4,I2,A8)" ) "  > ",id,". <u'w'>";                          id = id + 1 
              write(fid,"(A4,I2,A8)" ) "  > ",id,". <v'w'>";                          id = id + 1 
              write(fid,"(A4,I2,A14)") "  > ",id,". Mean x-force";                    id = id + 1 
              write(fid,"(A4,I2,A14)") "  > ",id,". Mean y-force";                    id = id + 1 
              write(fid,"(A4,I2,A14)") "  > ",id,". Mean z-force";                    id = id + 1 
              write(fid,"(A4,I2,A27)") "  > ",id,". Mean kinetic energy (MKE)";       id = id + 1 
              write(fid,*) " "
              write(fid,*) "Rij budgets:"
              write(fid,*) "  > R11 budget"
              write(fid,"(A6,I2,A15)") "    > ",id,". Unsteady term";        id = id + 1
              write(fid,"(A6,I2,A18)") "    > ",id,". Shear production";     id = id + 1
              write(fid,"(A6,I2,A18)") "    > ",id,". Force production";     id = id + 1
              write(fid,"(A6,I2,A22)") "    > ",id,". Convective transport"; id = id + 1
              write(fid,"(A6,I2,A21)") "    > ",id,". Turbulent transport";  id = id + 1
              write(fid,"(A6,I2,A20)") "    > ",id,". Pressure transport";   id = id + 1
              write(fid,"(A6,I2,A19)") "    > ",id,". Viscous transport";    id = id + 1
              write(fid,"(A6,I2,A21)") "    > ",id,". Viscous dissipation";  id = id + 1
              write(fid,"(A6,I2,A15)") "    > ",id,". SGS transport";        id = id + 1
              write(fid,"(A6,I2,A17)") "    > ",id,". SGS dissipation";      id = id + 1
              write(fid,"(A6,I2,A19)") "    > ",id,". Buoyancy transfer";    id = id + 1
              write(fid,"(A6,I2,A29)") "    > ",id,". Pressure-strain correlation";    id = id + 1
              write(fid,*) " "
              write(fid,*) "  > R22 budget"
              write(fid,"(A6,I2,A15)") "    > ",id,". Unsteady term";        id = id + 1
              write(fid,"(A6,I2,A18)") "    > ",id,". Shear production";     id = id + 1
              write(fid,"(A6,I2,A18)") "    > ",id,". Force production";     id = id + 1
              write(fid,"(A6,I2,A22)") "    > ",id,". Convective transport"; id = id + 1
              write(fid,"(A6,I2,A21)") "    > ",id,". Turbulent transport";  id = id + 1
              write(fid,"(A6,I2,A20)") "    > ",id,". Pressure transport";   id = id + 1
              write(fid,"(A6,I2,A19)") "    > ",id,". Viscous transport";    id = id + 1
              write(fid,"(A6,I2,A21)") "    > ",id,". Viscous dissipation";  id = id + 1
              write(fid,"(A6,I2,A15)") "    > ",id,". SGS transport";        id = id + 1
              write(fid,"(A6,I2,A17)") "    > ",id,". SGS dissipation";      id = id + 1
              write(fid,"(A6,I2,A19)") "    > ",id,". Buoyancy transfer";    id = id + 1
              write(fid,"(A6,I2,A29)") "    > ",id,". Pressure-strain correlation";    id = id + 1
              write(fid,*) " "
              write(fid,*) "  > R33 budget"
              write(fid,"(A6,I2,A15)") "    > ",id,". Unsteady term";        id = id + 1
              write(fid,"(A6,I2,A18)") "    > ",id,". Shear production";     id = id + 1
              write(fid,"(A6,I2,A18)") "    > ",id,". Force production";     id = id + 1
              write(fid,"(A6,I2,A22)") "    > ",id,". Convective transport"; id = id + 1
              write(fid,"(A6,I2,A21)") "    > ",id,". Turbulent transport";  id = id + 1
              write(fid,"(A6,I2,A20)") "    > ",id,". Pressure transport";   id = id + 1
              write(fid,"(A6,I2,A19)") "    > ",id,". Viscous transport";    id = id + 1
              write(fid,"(A6,I2,A21)") "    > ",id,". Viscous dissipation";  id = id + 1
              write(fid,"(A6,I2,A15)") "    > ",id,". SGS transport";        id = id + 1
              write(fid,"(A6,I2,A17)") "    > ",id,". SGS dissipation";      id = id + 1
              write(fid,"(A6,I2,A19)") "    > ",id,". Buoyancy transfer";    id = id + 1
              write(fid,"(A6,I2,A29)") "    > ",id,". Pressure-strain correlation";    id = id + 1
              write(fid,*) " "
              write(fid,*) "Velocity gradient variances"
              write(fid,"(A4,I2,A11)") "  > ",id,". Var(dudx)"; id = id + 1
              write(fid,"(A4,I2,A11)") "  > ",id,". Var(dudy)"; id = id + 1
              write(fid,"(A4,I2,A11)") "  > ",id,". Var(dudz)"; id = id + 1
              write(fid,"(A4,I2,A11)") "  > ",id,". Var(dvdx)"; id = id + 1
              write(fid,"(A4,I2,A11)") "  > ",id,". Var(dvdy)"; id = id + 1
              write(fid,"(A4,I2,A11)") "  > ",id,". Var(dvdz)"; id = id + 1
              write(fid,"(A4,I2,A11)") "  > ",id,". Var(dwdx)"; id = id + 1
              write(fid,"(A4,I2,A11)") "  > ",id,". Var(dwdy)"; id = id + 1
              write(fid,"(A4,I2,A11)") "  > ",id,". Var(dwdz)"; id = id + 1
              write(fid,"(A4,I2,A10)") "  > ",id,". Var(S12)" ; id = id + 1
              write(fid,"(A4,I2,A10)") "  > ",id,". Var(S13)" ; id = id + 1
              write(fid,"(A4,I2,A10)") "  > ",id,". Var(S23)" ; id = id + 1
              write(fid,*) " "
              write(fid,*) "SGS stress tensor variances"
              write(fid,"(A4,I2,A12)") "  > ",id,". Var(tau11)"; id = id + 1
              write(fid,"(A4,I2,A12)") "  > ",id,". Var(tau12)"; id = id + 1
              write(fid,"(A4,I2,A12)") "  > ",id,". Var(tau13)"; id = id + 1
              write(fid,"(A4,I2,A12)") "  > ",id,". Var(tau22)"; id = id + 1
              write(fid,"(A4,I2,A12)") "  > ",id,". Var(tau23)"; id = id + 1
              write(fid,"(A4,I2,A12)") "  > ",id,". Var(tau33)"; id = id + 1
              write(fid,*) " "
              write(fid,*) "TKE fluxes"
              write(fid,"(A4,I2,A14)") "  > ",id,". <u'ui'ui'>/2"; id = id + 1
              write(fid,"(A4,I2,A14)") "  > ",id,". <v'ui'ui'>/2"; id = id + 1
              write(fid,"(A4,I2,A14)") "  > ",id,". <w'ui'ui'>/2"; id = id + 1
              
              close(fid)
          end if
          call MPI_Barrier(MPI_COMM_WORLD,ierr)

          if (this%sim%isStratified .or. this%sim%useScalars) then
              write(tempname,"(A3,I2.2,A27)")"Run",this%sim%runID,"_scalar_stats_index_key.out"
              fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)

              if (nrank == 0) then
                  id = 1
                  open(fid, file=fname)!, status="new", action="write")
                  write(fid,*) "Scalar variance budget terms"
                  write(fid,"(A4,I2,A15)") "  > ",id,". Unsteady term";            id = id + 1
                  write(fid,"(A4,I2,A26)") "  > ",id,". Mean gradient production"; id = id + 1
                  write(fid,"(A4,I2,A22)") "  > ",id,". Convective transport";     id = id + 1
                  write(fid,"(A4,I2,A21)") "  > ",id,". Turbulent transport";      id = id + 1
                  write(fid,"(A4,I2,A21)") "  > ",id,". Molecular transport";      id = id + 1
                  write(fid,"(A4,I2,A15)") "  > ",id,". SGS transport";            id = id + 1
                  write(fid,"(A4,I2,A23)") "  > ",id,". Molecular destruction";    id = id + 1
                  write(fid,"(A4,I2,A17)") "  > ",id,". SGS destruction";          id = id + 1
                  write(fid,*) " "
                  write(fid,*) "Misc terms:"
                  write(fid,"(A4,I2,A28)") "  > ",id,". Mean scalar gradient, dTdz"; id = id + 1
                  write(fid,"(A4,I2,A14)") "  > ",id,". Mean scalar ";               id = id + 1 
                  write(fid,"(A4,I2,A8)" ) "  > ",id,". <u'T'>";                     id = id + 1 
                  write(fid,"(A4,I2,A8)" ) "  > ",id,". <v'T'>";                     id = id + 1
                  write(fid,"(A4,I2,A8)" ) "  > ",id,". <w'T'>";                     id = id + 1 
                  write(fid,"(A4,I2,A8)" ) "  > ",id,". <T'T'>";                     id = id + 1 
                  write(fid,*) " "
                  write(fid,*) "Scalar flux budgets:"
                  write(fid,*) "  > <w'T'> budget"
                  write(fid,"(A6,I2,A15)") "    > ",id,". Unsteady term";                     id = id + 1
                  write(fid,"(A6,I2,A33)") "    > ",id,". Mean scalar gradient production";   id = id + 1
                  write(fid,"(A6,I2,A23)") "    > ",id,". Mean shear production";             id = id + 1
                  write(fid,"(A6,I2,A18)") "    > ",id,". Force production";                  id = id + 1
                  write(fid,"(A6,I2,A21)") "    > ",id,". Buoyancy production";               id = id + 1
                  write(fid,"(A6,I2,A22)") "    > ",id,". Convective transport";              id = id + 1
                  write(fid,"(A6,I2,A21)") "    > ",id,". Turbulent transport";               id = id + 1
                  write(fid,"(A6,I2,A20)") "    > ",id,". Pressure transport";                id = id + 1
                  write(fid,"(A6,I2,A30)") "    > ",id,". Molecular transport (w*dTdz)";      id = id + 1
                  write(fid,"(A6,I2,A30)") "    > ",id,". Molecular transport (T*dwdz)";      id = id + 1
                  write(fid,"(A6,I2,A21)") "    > ",id,". SGS transport (w*q)";               id = id + 1
                  write(fid,"(A6,I2,A23)") "    > ",id,". SGS transport (T*tau)";             id = id + 1
                  write(fid,"(A6,I2,A21)") "    > ",id,". Pressure scrambling";               id = id + 1
                  write(fid,"(A6,I2,A23)") "    > ",id,". Molecular destruction";             id = id + 1
                  write(fid,"(A6,I2,A23)") "    > ",id,". Molecular destruction";             id = id + 1
                  write(fid,"(A6,I2,A28)") "    > ",id,". SGS destruction (qi*dwdxi)";        id = id + 1
                  write(fid,"(A6,I2,A32)") "    > ",id,". SGS destruction (tau_i3*dTdxi)";    id = id + 1
                  write(fid,*) " "
                  write(fid,*) "Temperature gradient variances"
                  write(fid,"(A4,I2,A11)") "  > ",id,". Var(dTdx)"; id = id + 1
                  write(fid,"(A4,I2,A11)") "  > ",id,". Var(dTdy)"; id = id + 1
                  write(fid,"(A4,I2,A11)") "  > ",id,". Var(dTdz)"; id = id + 1
                  write(fid,*) " "
                  write(fid,*) "SGS heat flux moments"
                  write(fid,"(A4,I2,A9 )") "  > ",id,". Var(q1)" ; id = id + 1
                  write(fid,"(A4,I2,A9 )") "  > ",id,". Var(q2)" ; id = id + 1
                  write(fid,"(A4,I2,A9 )") "  > ",id,". Var(q3)" ; id = id + 1
                  write(fid,"(A4,I2,A10)") "  > ",id,". Mean(q3)"; id = id + 1
                  write(fid,*) " "
                  write(fid,*) "Scalar variance fluxes"
                  write(fid,"(A4,I2,A10)") "  > ",id,". <u'T'T'>"; id = id + 1
                  write(fid,"(A4,I2,A10)") "  > ",id,". <v'T'T'>"; id = id + 1
                  write(fid,"(A4,I2,A10)") "  > ",id,". <w'T'T'>"; id = id + 1
                  
                  close(fid)
              end if
              call MPI_Barrier(MPI_COMM_WORLD,ierr)
          end if
      end subroutine

      subroutine write_time_info(this)
          class(stats_xy), intent(inout) :: this
          character(len=clen) :: tempname, fname, fname2
          logical :: exists 
          integer :: n, ierr, fid1, fid2

          fid1 = 12
          fid2 = 13

          write(tempname,"(A3,I2.2,A18)")"Run",this%sim%runID,"_time_info_all.out"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)

          if (nrank == 0) then
              inquire(file=fname, exist=exists)
              if (.not. exists) then
                  open(fid1, file=fname, status="new", action="write")
                  close(fid1)
              end if
              do n = 1,this%nwrite
                  ! Append time info to master file
                  open(fid1, file=fname, status="old", position="append", action="write")
                  write(fid1,*) this%time(n), this%tid(n)
                  close(fid1)

                  ! Dump a single file with tid and time for the curren id
                  write(tempname,"(A3,I2.2,A12,I6.6,A4)")"Run",this%sim%runID,&
                    "_time_info_t",this%tid(n),".out"
                  fname2 = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
                  inquire(file=fname2, exist=exists)
                  if (exists) then
                      open(fid2, file=fname2, status="replace", action="write")
                  else
                      open(fid2, file=fname2, status="new", action="write")
                  end if
                  write(fid2,*) this%time(n), this%tid(n)
                  close(fid2)
              end do
          end if 
          
          call MPI_Barrier(MPI_COMM_WORLD,ierr)

      end subroutine

      subroutine mean(this, f, fmean)
          class(stats_xy), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: f
          real(rkind), dimension(:), intent(out) :: fmean
          integer :: ierr

          call this%sim%spectC%fft(f, this%sim%cbuffyC(:,:,:,1))
          call transpose_y_to_z(this%sim%cbuffyC(:,:,:,1), this%sim%cbuffzC(:,:,:,1), this%sim%sp_gpC)
          if (nrank == 0) then
              fmean = real(this%sim%cbuffzC(1,1,:,1),rkind)*this%avgFact
          else
              fmean = 0.d0 ! Only 0 processor has the actual mean  
          end if 
          call MPI_BCAST(fmean,this%sim%nz,MPIRKIND,0,MPI_COMM_WORLD, ierr)

      end subroutine 
  
      subroutine covariance(this,f, g, covar_fg)
          class(stats_xy), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: f, g
          real(rkind), dimension(:), intent(out) :: covar_fg

          call this%mean(f,  this%zbuff_(:,1))
          call this%mean(g,  this%zbuff_(:,2))
          this%sim%rbuffxC(:,:,:,4) = f*g
          call this%mean(this%sim%rbuffxC(:,:,:,4), covar_fg)

          covar_fg = covar_fg - this%zbuff_(:,1)*this%zbuff_(:,2)
        
      end subroutine

    subroutine ddz(this, f, ddz_f, bc_bot, bc_top)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:), intent(in)  :: f
        real(rkind), dimension(:), intent(out) :: ddz_f
        integer, intent(in) :: bc_bot, bc_top 

        this%zbuff_(:,1) = f
        call this%sim%Pade6opZ%ddz_1d_C2C(this%zbuff_(:,1),this%zbuff_(:,2), bc_bot, bc_top)
        ddz_f = this%zbuff_(:,2)

    end subroutine 

    subroutine d2dz2(this,f,d2dz2_f, bc_bot, bc_top)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:), intent(in) :: f
        real(rkind), dimension(:), intent(out) :: d2dz2_f
        integer, intent(in) :: bc_bot, bc_top

        call this%sim%Pade6opZ%d2dz2_C2C_1d(f,d2dz2_f,bc_bot,bc_top)

    end subroutine

    subroutine turbulent_transport(this,f,g,dfdx,dfdy,dfdz,dgdx,dgdy,dgdz,&
        ConvTrans_fg,production,TurbTrans_fg)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in) :: f, g, dfdx, dfdy, dfdz, dgdx, dgdy, dgdz
        real(rkind), dimension(:), intent(in) :: ConvTrans_fg, production
        real(rkind), dimension(:), intent(out) :: TurbTrans_fg

        this%sim%rbuffxC(:,:,:,1) = this%sim%u*dgdx + this%sim%v*dgdy + this%sim%wC*dgdz
        this%sim%rbuffxC(:,:,:,2) = this%sim%u*dfdx + this%sim%v*dfdy + this%sim%wC*dfdz
        call this%covariance(f,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,1))
        call this%covariance(g,this%sim%rbuffxC(:,:,:,2),this%zbuff(:,2))
        TurbTrans_fg = -this%zbuff(:,1) - this%zbuff(:,2) + ConvTrans_fg - production
    end subroutine

    subroutine gradient_production(this,dgdx,dgdy,dgdz,fu,fv,fw,production)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in) :: dgdx, dgdy, dgdz
        real(rkind), dimension(:), intent(in) :: fu, fv, fw
        real(rkind), dimension(:), intent(out) :: production

        call this%mean(dgdx,this%zbuff(:,1))
        call this%mean(dgdy,this%zbuff(:,2))
        call this%mean(dgdz,this%zbuff(:,3))

        production = -fu*this%zbuff(:,1) - fv*this%zbuff(:,2) - fw*this%zbuff(:,3)

    end subroutine

    subroutine convective_transport(this,fg,convTrans)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:), intent(in) :: fg
        real(rkind), dimension(:), intent(out) :: convTrans

        call this%ddz(fg,convTrans,0,0)
        convTrans = this%meanW*convTrans
    end subroutine

    subroutine pseudo_dissipation(this,dfdx,dfdy,dfdz,dgdx,dgdy,dgdz,viscosity,dissipation)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in) :: dfdx, dfdy, dfdz, dgdx, dgdy, dgdz
        real(rkind), intent(in) :: viscosity
        real(rkind), dimension(:), intent(out) :: dissipation

        call this%covariance(dfdx, dgdx, this%zbuff(:,1))
        call this%covariance(dfdy, dgdy, this%zbuff(:,2))
        call this%covariance(dfdz, dgdz, this%zbuff(:,3))

        dissipation = viscosity*(this%zbuff(:,1) + this%zbuff(:,2) + this%zbuff(:,3))

    end subroutine

    subroutine molecular_transport(this,f,dgdz,diffusivity,molecTrans)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in) :: f, dgdz
        real(rkind), intent(in) :: diffusivity ! the molecular diffusivities
        real(rkind), dimension(:), intent(out) :: molecTrans

        call this%covariance(f,dgdz,this%zbuff(:,1))
        call this%ddz(this%zbuff(:,1),this%zbuff(:,3),0,0)

        molecTrans = diffusivity*this%zbuff(:,3)

    end subroutine

    subroutine SGS_transport(this,f,SGS_g,SGStrans)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in) :: f, SGS_g
        real(rkind), dimension(:), intent(out) :: SGStrans

        call this%covariance(f,SGS_g,this%zbuff(:,1))
        this%zbuff(:,1) = -1.d0*this%zbuff(:,1)
        call this%ddz(this%zbuff(:,1),SGStrans,0,0)

    end subroutine

    subroutine SGS_dissipation(this,dfdx,dfdy,dfdz,gSGS1,gSGS2,gSGS3,dissipation)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in) :: dfdx, dfdy, dfdz
        real(rkind), dimension(:,:,:), intent(in) :: gSGS1, gSGS2, gSGS3 
        real(rkind), dimension(:), intent(out) :: dissipation

        call this%covariance(gSGS1,dfdx,this%zbuff(:,1))
        call this%covariance(gSGS2,dfdy,this%zbuff(:,2))
        call this%covariance(gSGS3,dfdz,this%zbuff(:,3))
        dissipation = -1.d0*(this%zbuff(:,1) + this%zbuff(:,2) + this%zbuff(:,3))

    end subroutine

    subroutine compute_unsteady_terms(this)
        class(stats_xy), intent(inout) :: this

        ! Step 1: store TKE and temperature variacne required for unsteady term at t^n-1
        if (mod(this%sim%step+1,this%compute_freq) == 0) then
            call this%covariance(this%sim%u,  this%sim%u,  this%qold(:,2))
            call this%covariance(this%sim%v,  this%sim%v,  this%qold(:,3))
            call this%covariance(this%sim%wC, this%sim%wC, this%qold(:,4))
            !call this%covariance(this%sim%T,  this%sim%T,  this%qold(:,5))
            !call this%covariance(this%sim%wC, this%sim%T,  this%qold(:,6))

            this%qold(:,1) = 0.5d0*(this%qold(:,2) + this%qold(:,3) + this%qold(:,4))
            
        ! Step 2: Compute ddt terms (using t^n+1 if dt=const or t^n
        ! otherwise)
        elseif (mod(this%sim%step,this%compute_freq) == 0) then
            if (this%sim%CFL >= 0.d0) then
                call this%covariance(this%sim%u,  this%sim%u,  this%qnew(:,2))
                call this%covariance(this%sim%v,  this%sim%v,  this%qnew(:,3))
                call this%covariance(this%sim%wC, this%sim%wC, this%qnew(:,4))
                !call this%covariance(this%sim%T,  this%sim%T,  this%qnew(:,5))
                !call this%covariance(this%sim%wC, this%sim%T,  this%qnew(:,6))
                this%qnew(:,1) = 0.5d0*(this%qnew(:,2) + this%qnew(:,3) + this%qnew(:,4))

                this%tke_budget(:,1) = (this%qnew(:,1) - this%qold(:,1))/(this%sim%dt)
                this%R11_budget(:,1) = (this%qnew(:,2) - this%qold(:,2))/(this%sim%dt)
                this%R22_budget(:,1) = (this%qnew(:,3) - this%qold(:,3))/(this%sim%dt)
                this%R33_budget(:,1) = (this%qnew(:,4) - this%qold(:,4))/(this%sim%dt)
                !this%TT_budget(:,1) = (this%qnew(:,2) - this%qold(:,2))/(this%sim%dt)
            end if
        elseif( (mod(this%sim%step-1,this%compute_freq) == 0) .and. &
                (this%nwrite > 0)                           ) then
            if (this%sim%CFL < 0.d0) then
                call this%covariance(this%sim%u,  this%sim%u,  this%qnew(:,2))
                call this%covariance(this%sim%v,  this%sim%v,  this%qnew(:,3))
                call this%covariance(this%sim%wC, this%sim%wC, this%qnew(:,4))
                !call this%covariance(this%sim%T,  this%sim%T,  this%qnew(:,2))
                !call this%covariance(this%sim%wC, this%sim%T,  this%qnew(:,6))
                this%qnew(:,1) = 0.5d0*(this%qnew(:,2) + this%qnew(:,3) + this%qnew(:,4))

                this%tke_budget(:,1) = (this%qnew(:,1) - this%qold(:,1))/(2.d0*this%sim%dt)
                this%R11_budget(:,1) = (this%qnew(:,2) - this%qold(:,2))/(2.d0*this%sim%dt)
                this%R22_budget(:,1) = (this%qnew(:,3) - this%qold(:,3))/(2.d0*this%sim%dt)
                this%R33_budget(:,1) = (this%qnew(:,4) - this%qold(:,4))/(2.d0*this%sim%dt)
                !this%TT_budget(:,1) = (this%qnew(:,2) - this%qold(:,2))/(2.d0*this%sim%dt)
            end if
        end if
    end subroutine

end module stats_xy_mod
