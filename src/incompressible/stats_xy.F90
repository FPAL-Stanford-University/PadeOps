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
    use spectralMod,        only: spectral
    use reductions,         only: p_maxval
    !use stdlib_sorting,     only: sort_index
    implicit none

    private
    public :: stats_xy
    type :: stats_xy

        private
        ! Stats
        real(rkind), dimension(:,:,:,:), allocatable :: stats
        real(rkind), dimension(:,:,:,:,:), allocatable :: stats_sca
        real(rkind), dimension(:,:), pointer :: tke_budget, TT_budget, R11_budget, &
          R22_budget, R33_budget, wT_budget, duidxj_var, dTdxj_var, tauij_var, qj_var, &
          tke_flux, sca_var_flux
        real(rkind), dimension(:), pointer :: mke, tke, fifi, meanU, meanV, meanW, meanT, meanP, &
          meanFx, meanFy, meanFz, meanFT, fTfT, &
          uu, vv, ww, uv, uw, vw, uT, vT, wT, TT, dTdz, S12_var, S13_var, S23_var, &
          meanq3, pp, up, vp, wp

        ! Scalar fields
        real(rkind), dimension(:,:,:), pointer :: T, dTdxC, dTdyC, dTdzC, q1, q2, q3, &
          T_all, dTdxC_all, dTdyC_all, dTdzC_all
        real(rkind), dimension(:,:,:,:), allocatable :: Tsplit
        real(rkind), dimension(:,:,:,:,:), allocatable :: dTdxj

        ! Momentum fields
        real(rkind), dimension(:,:,:,:), allocatable :: psplit, fxsplit, fysplit, fzsplit, fTsplit
        real(rkind), dimension(:,:,:,:,:), allocatable,public :: ujsplit, duidxj 
        real(rkind), dimension(:,:,:), pointer :: u, v, wC, pressure, dudx, dudy, dudz, &
          dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, fx, fy, fz, fT

        ! SGS fields
        real(rkind), dimension(:,:,:,:,:), allocatable :: qjsplit
        real(rkind), dimension(:,:,:,:,:), allocatable :: tauij
        real(rkind), dimension(:,:,:), pointer :: tau11,tau12,tau13,tau22,tau23,tau33

        ! Buffers
        real(rkind), dimension(:,:), allocatable :: zbuff, zbuff_
        real(rkind), dimension(:,:,:), allocatable :: uvarold, uvarnew
        real(rkind), dimension(:,:,:,:), allocatable :: svarold, svarnew
        integer, dimension(:), allocatable :: tid

        ! Spectral
        type(spectral), pointer :: spectC
        real(rkind), dimension(:,:,:), allocatable, public :: k2d, mask
        real(rkind), dimension(:), allocatable :: kc_vec
        real(rkind) :: kmax

        ! Misc
        real(rkind), dimension(:), allocatable :: time
        integer :: tid_start, compute_freq, dump_freq, tidx, tid_ddt
        integer :: nwrite, nscalars, nscales
        character(len=clen) :: outputdir
        type(igrid), pointer :: sim
        real(rkind) :: avgFact
        logical :: do_stats, scale_split

        contains
          procedure          :: init
          procedure          :: destroy
          procedure          :: compute_stats

          ! Memory-mangament and IO
          procedure, private :: dump_stats
          procedure, private :: write_time_info
          procedure, private :: write_labels
          procedure, private :: link_pointers

          ! Budgets
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
          procedure, private :: interscale_transfer
          procedure, private :: mixed_scale_turbulent_transport

          ! Stats
          procedure, private :: covariance
          procedure, private :: mean

          ! Derivatives
          procedure, private :: ddz
          procedure, private :: d2dz2

          ! Scale-splitting
          procedure, private :: do_scale_splitting_multiple
          procedure, private :: do_scale_splitting_single
          generic,   private :: do_scale_splitting => do_scale_splitting_multiple, do_scale_splitting_single
          procedure, private :: copy_and_scale_split
          procedure, private :: low_pass_filter

          ! Accessors
          procedure          :: copy_stats
          procedure          :: set_tidx
          procedure          :: set_nwrite
          procedure          :: print_var
    end type
    contains

      subroutine init(this,inputfile,sim)
        class(stats_xy), intent(inout), target :: this
        character(len=*), intent(in) :: inputfile
        type(igrid), intent(in), target :: sim
        integer :: nscales = 1
        real(rkind), dimension(:), allocatable :: kc_for_scale_splitting
        integer :: tid_start = 1000000, compute_freq = 1000000, dump_freq = 10000000
        character(len=clen) :: outputdir
        integer :: ioUnit, ierr, nstore, n
        integer, parameter :: nterms      = 102
        integer, parameter :: nterms_sca  = 51
        integer, parameter :: n_ddt_terms = 6
        real, parameter :: zero = 0.d0
        logical :: do_stats = .false.

        namelist /STATS_XY/ tid_start, compute_freq, dump_freq, outputdir, &
          do_stats, nscales
        namelist /SCALE_SPLITTING/ kc_for_scale_splitting

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
            if (nscales > 1) then
                allocate(kc_for_scale_splitting(nscales))
                open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
                read(unit=ioUnit, NML=SCALE_SPLITTING)
                close(ioUnit)

                this%scale_split = .true.
                this%nscales = nscales
            else 
                this%scale_split = .false.
                this%nscales = 1
            end if
            call message(0,"Initializing stats class")
            call message(1,"tid_start    ", tid_start)
            call message(1,"compute_freq ",compute_freq)
            call message(1,"dump_freq    ",dump_freq)
            call message(1,"outputdir    "//trim(this%outputdir))
            call message(1,"do_stats     ",do_stats)
            call message(1,"scale_split  ",this%scale_split)

            this%sim    => sim
            this%spectC => this%sim%spectC
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

            ! Set things up for scale-splitting
            allocate(this%kc_vec(this%nscales + 1))
            allocate(this%k2d( this%sim%sp_gpC%ysz(1),this%sim%sp_gpC%ysz(2),this%sim%sp_gpC%ysz(3)))
            allocate(this%mask(this%sim%sp_gpC%ysz(1),this%sim%sp_gpC%ysz(2),this%sim%sp_gpC%ysz(3)))
            this%k2d = sqrt(this%spectC%k1**2.d0 + this%spectC%k2**2.d0)
            this%kmax = p_maxval(maxval(this%k2d))
            this%kc_vec(1) = 0.d0
            this%kc_vec(this%nscales+1) = this%kmax
            if (this%scale_split) then
                this%kc_vec(2:this%nscales) = kc_for_scale_splitting
            end if
            call assert(maxval(this%kc_vec) <= this%kmax,"kc is larger than the Nyquist wavenumber -- stats_xy.F90")

            ! Allocate memory and link pointers
            if (this%nscalars > 0) then
                allocate(this%stats_sca(sim%nz,nterms_sca,nstore,this%nscalars,this%nscales))
                allocate(this%qjsplit(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales,3))
                allocate(this%Tsplit(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales))
                allocate(this%dTdxj(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales,3))
            end if
            allocate(this%stats(sim%nz,nterms,nstore,this%nscales))
            allocate(this%time(nstore))
            allocate(this%tid(nstore))
            allocate(this%zbuff(sim%nz,4))
            allocate(this%zbuff_(sim%nz,2))
            allocate(this%uvarold(sim%nz,this%nscales,3), this%uvarnew(sim%nz,this%nscales,3))
            allocate(this%svarold(sim%nz,this%nscalars,this%nscales,2), this%svarnew(sim%nz,this%nscalars,this%nscales,2))

            ! These arrays are required to have all scales defined simultaneously for computing interscale transfer
            allocate(this%duidxj(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales,9))
            allocate(this%ujsplit(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales,3))

            ! These arrays do not require it, but for now, we store all scales. If memory becomes an issue, we can change this
            allocate(this%tauij(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales,6))
            allocate(this%psplit(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales))
            allocate(this%fxsplit(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales))
            allocate(this%fysplit(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales))
            allocate(this%fzsplit(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales))
            allocate(this%fTsplit(this%sim%gpC%xsz(1),this%sim%gpC%xsz(2),this%sim%gpC%xsz(3),this%nscales))
            
            this%tidx = 0
            call this%link_pointers(this%tidx+1,min(1,this%nscalars),1)
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
            this%vT,this%wT,this%TT,this%meanFx,this%meanFy,this%meanFz,&
            this%pp,this%up,this%vp,this%wp,this%meanfT,this%fTfT)
          if (allocated(this%stats)) deallocate(this%stats)
          if (allocated(this%stats_sca)) deallocate(this%stats_sca)
          if (allocated(this%Tsplit)) deallocate(this%Tsplit)
          if (allocated(this%dTdxj)) deallocate(this%dTdxj)
          if (allocated(this%psplit)) deallocate(this%psplit)
          if (allocated(this%fxsplit)) deallocate(this%fxsplit)
          if (allocated(this%fysplit)) deallocate(this%fysplit)
          if (allocated(this%fzsplit)) deallocate(this%fzsplit)
          if (allocated(this%fTsplit)) deallocate(this%fTsplit)
          if (allocated(this%ujsplit)) deallocate(this%ujsplit)
          if (allocated(this%duidxj)) deallocate(this%duidxj)
          if (allocated(this%dTdxj)) deallocate(this%dTdxj)
          if (allocated(this%qjsplit)) deallocate(this%qjsplit)
          if (allocated(this%tauij)) deallocate(this%tauij)
          if (allocated(this%zbuff)) deallocate(this%zbuff)
          if (allocated(this%zbuff_)) deallocate(this%zbuff_)
          if (allocated(this%uvarold)) deallocate(this%uvarold)
          if (allocated(this%uvarnew)) deallocate(this%uvarnew)
          if (allocated(this%svarold)) deallocate(this%svarold)
          if (allocated(this%svarnew)) deallocate(this%svarnew)
          if (allocated(this%tid)) deallocate(this%tid)
          if (allocated(this%k2d)) deallocate(this%k2d)
          if (allocated(this%mask)) deallocate(this%mask)
          if (allocated(this%kc_vec)) deallocate(this%kc_vec)
          if (allocated(this%time)) deallocate(this%time)

          if (associated(this%sim)) nullify(this%sim)
          if (associated(this%q1)) nullify(this%q1)
          if (associated(this%q2)) nullify(this%q2)
          if (associated(this%q3)) nullify(this%q3)
          if (associated(this%tke_budget)) nullify(this%tke_budget)
          if (associated(this%TT_budget)) nullify(this%TT_budget)
          if (associated(this%R11_budget)) nullify(this%R11_budget)
          if (associated(this%R22_budget)) nullify(this%R22_budget)
          if (associated(this%R33_budget)) nullify(this%R33_budget)
          if (associated(this%wT_budget)) nullify(this%wT_budget)
          if (associated(this%duidxj_var)) nullify(this%duidxj_var)
          if (associated(this%dTdxj_var)) nullify(this%dTdxj_var)
          if (associated(this%tauij_var)) nullify(this%tauij_var)
          if (associated(this%qj_var)) nullify(this%qj_var)
          if (associated(this%tke_flux)) nullify(this%tke_flux)
          if (associated(this%sca_var_flux)) nullify(this%sca_var_flux)
          if (associated(this%mke)) nullify(this%mke)
          if (associated(this%tke)) nullify(this%tke)
          if (associated(this%fifi)) nullify(this%fifi)
          if (associated(this%meanU)) nullify(this%meanU)
          if (associated(this%meanV)) nullify(this%meanV)
          if (associated(this%meanW)) nullify(this%meanW)
          if (associated(this%meanT)) nullify(this%meanT)
          if (associated(this%meanP)) nullify(this%meanP)
          if (associated(this%meanFx)) nullify(this%meanFx)
          if (associated(this%meanFy)) nullify(this%meanFy)
          if (associated(this%meanFz)) nullify(this%meanFz)
          if (associated(this%meanFT)) nullify(this%meanFT)
          if (associated(this%uu)) nullify(this%uu)
          if (associated(this%vv)) nullify(this%vv)
          if (associated(this%ww)) nullify(this%ww)
          if (associated(this%uv)) nullify(this%uv)
          if (associated(this%uw)) nullify(this%uw)
          if (associated(this%vw)) nullify(this%vw)
          if (associated(this%uT)) nullify(this%uT)
          if (associated(this%vT)) nullify(this%vT)
          if (associated(this%wT)) nullify(this%wT)
          if (associated(this%TT)) nullify(this%TT)
          if (associated(this%dTdz)) nullify(this%dTdz)
          if (associated(this%S12_var)) nullify(this%S12_var)
          if (associated(this%S13_var)) nullify(this%S13_var)
          if (associated(this%S23_var)) nullify(this%S23_var)
          if (associated(this%meanq3)) nullify(this%meanq3)
          if (associated(this%pp)) nullify(this%pp)
          if (associated(this%up)) nullify(this%up)
          if (associated(this%vp)) nullify(this%vp)
          if (associated(this%wp)) nullify(this%wp)
          if (associated(this%T)) nullify(this%T)
          if (associated(this%dTdxC)) nullify(this%dTdxC)
          if (associated(this%dTdyC)) nullify(this%dTdyC)
          if (associated(this%dTdzC)) nullify(this%dTdzC)
          if (associated(this%T_all)) nullify(this%T_all)
          if (associated(this%dTdxC_all)) nullify(this%dTdxC_all)
          if (associated(this%dTdyC_all)) nullify(this%dTdyC_all)
          if (associated(this%dTdzC_all)) nullify(this%dTdzC_all)
          if (associated(this%u)) nullify(this%u)
          if (associated(this%v)) nullify(this%v)
          if (associated(this%wC)) nullify(this%wC)
          if (associated(this%pressure)) nullify(this%pressure)
          if (associated(this%dudx)) nullify(this%dudx)
          if (associated(this%dudy)) nullify(this%dudy)
          if (associated(this%dudz)) nullify(this%dudz)
          if (associated(this%dvdx)) nullify(this%dvdx)
          if (associated(this%dvdy)) nullify(this%dvdy)
          if (associated(this%dvdz)) nullify(this%dvdz)
          if (associated(this%dwdx)) nullify(this%dwdx)
          if (associated(this%dwdy)) nullify(this%dwdy)
          if (associated(this%dwdz)) nullify(this%dwdz)
          if (associated(this%fx)) nullify(this%fx)
          if (associated(this%fy)) nullify(this%fy)
          if (associated(this%fz)) nullify(this%fz)
          if (associated(this%fT)) nullify(this%fT)
          if (associated(this%tau11)) nullify(this%tau11)
          if (associated(this%tau12)) nullify(this%tau12)
          if (associated(this%tau13)) nullify(this%tau13)
          if (associated(this%tau22)) nullify(this%tau22)
          if (associated(this%tau23)) nullify(this%tau23)
          if (associated(this%spectC)) nullify(this%spectC)
      end subroutine

      subroutine copy_and_scale_split(this,kc,sca,scalar_only)
          class(stats_xy), intent(inout) :: this
          real(rkind), dimension(:), intent(in) :: kc
          integer, intent(in) :: sca
          real(rkind), dimension(:,:,:,:), pointer :: rbuff
          complex(rkind), dimension(:,:,:), pointer :: cbuff
          real(rkind), dimension(:,:,:), pointer :: F, dFdx, dFdy, dFdz, q1, q2, q3
          logical, intent(in) :: scalar_only
          integer :: ij, n

          ! nullify pointers
          F    => null()
          dFdx => null()
          dFdy => null()
          dFdz => null()
          q1   => null()
          q2   => null()
          q3   => null()
          rbuff => null()
          cbuff => null()

          associate( rbuff => this%sim%rbuffxC(:,:,:,1:2), &
                     cbuff => this%sim%cbuffyC(:,:,:,1))

              if (this%nscalars > 0) then
                  if (this%sim%isStratified) then
                      if (sca == 1) then
                          F    => this%sim%T
                          dFdx => this%sim%dTdxC
                          dFdy => this%sim%dTdyC
                          dFdz => this%sim%dTdzC
                          q1   => this%sim%q1_T
                          q2   => this%sim%q2_T
                          q3   => this%sim%q3_T
                      else
                          F    => this%sim%scalars(sca-1)%F
                          dFdx => this%sim%scalars(sca-1)%dFdxC
                          dFdy => this%sim%scalars(sca-1)%dFdyC
                          dFdz => this%sim%scalars(sca-1)%dFdzC
                          q1   => this%sim%scalars(sca-1)%q1
                          q2   => this%sim%scalars(sca-1)%q2
                          q3   => this%sim%scalars(sca-1)%q3
                      end if
                  else
                      F    => this%sim%scalars(sca)%F
                      dFdx => this%sim%scalars(sca)%dFdxC
                      dFdy => this%sim%scalars(sca)%dFdyC
                      dFdz => this%sim%scalars(sca)%dFdzC
                      q1   => this%sim%scalars(sca)%q1
                      q2   => this%sim%scalars(sca)%q2
                      q3   => this%sim%scalars(sca)%q3
                  end if

                  ! Scalar
                  call this%do_scale_splitting(F,this%Tsplit,kc,rbuff,cbuff)

                  ! Scalar gradient
                  call this%do_scale_splitting(dFdx,this%dTdxj(:,:,:,:,1),kc,rbuff,cbuff)
                  call this%do_scale_splitting(dFdy,this%dTdxj(:,:,:,:,2),kc,rbuff,cbuff)
                  call this%do_scale_splitting(dFdz,this%dTdxj(:,:,:,:,3),kc,rbuff,cbuff)
                  
                  ! SGS flux
                  call this%do_scale_splitting(q1,this%qjsplit(:,:,:,:,1),kc,rbuff,cbuff)
                  call this%do_scale_splitting(q2,this%qjsplit(:,:,:,:,2),kc,rbuff,cbuff)
                  call this%sim%spectForceLayer%interpE2C(q3, this%sim%rbuffxC(:,:,:,3), &
                    this%sim%rbuffyC(:,:,:,1), this%sim%rbuffzC(:,:,:,1), & 
                    this%sim%rbuffyE(:,:,:,1), this%sim%rbuffzE(:,:,:,1))
                  call this%do_scale_splitting(this%sim%rbuffxC(:,:,:,3),this%qjsplit(:,:,:,:,3),kc,rbuff,cbuff)

                  ! Forcing
                  call this%do_scale_splitting(this%sim%spectForceLayer%fT,this%fTsplit,kc,rbuff,cbuff)

                  nullify(F, dFdx, dFdy, dFdz, q1, q2, q3)
              endif

              if (scalar_only) then
                  continue
              else
                  ! Velocity
                  call this%do_scale_splitting(this%sim%u ,this%ujsplit(:,:,:,:,1),kc,rbuff,cbuff)
                  call this%do_scale_splitting(this%sim%v ,this%ujsplit(:,:,:,:,2),kc,rbuff,cbuff)
                  call this%do_scale_splitting(this%sim%wC,this%ujsplit(:,:,:,:,3),kc,rbuff,cbuff)

                  ! Velocity gradient
                  do n = 1,9
                      call this%do_scale_splitting(this%sim%duidxjC(:,:,:,n),this%duidxj(:,:,:,:,n),kc,rbuff,cbuff)
                  end do

                  ! Pressure
                  call this%do_scale_splitting(this%sim%pressure,this%psplit,kc,rbuff,cbuff)

                  ! Forcing
                  call this%do_scale_splitting(this%sim%spectForceLayer%fx,this%fxsplit,kc,rbuff,cbuff)
                  call this%do_scale_splitting(this%sim%spectForceLayer%fy,this%fysplit,kc,rbuff,cbuff)
                  call this%sim%spectForceLayer%interpE2C(&
                    this%sim%spectForceLayer%fz  , this%sim%rbuffxC(:,:,:,3), &
                    this%sim%rbuffyC(:,:,:,1), this%sim%rbuffzC(:,:,:,1), & 
                    this%sim%rbuffyE(:,:,:,1), this%sim%rbuffzE(:,:,:,1))
                  call this%do_scale_splitting(this%sim%rbuffxC(:,:,:,3),this%fzsplit,kc,rbuff,cbuff)

                  ! tauij_SGS
                  call this%sim%sgsModel%populate_tauij_E_to_C()
                  call this%sim%sgsModel%get_tauijC(this%tauij(:,:,:,1,1),&
                    this%tauij(:,:,:,1,2),this%tauij(:,:,:,1,3),this%tauij(:,:,:,1,4),&
                    this%tauij(:,:,:,1,5),this%tauij(:,:,:,1,6))
                  do ij = 1,6
                      this%sim%rbuffxC(:,:,:,3) = this%tauij(:,:,:,1,ij)
                      call this%do_scale_splitting(this%sim%rbuffxC(:,:,:,3),this%tauij(:,:,:,:,ij),kc,rbuff,cbuff)
                  end do

              end if

          end associate

          ! nullify pointers
          if (associated(F    )) nullify(F) 
          if (associated(dFdx )) nullify(dFdx)
          if (associated(dFdy )) nullify(dFdy)
          if (associated(dFdz )) nullify(dFdz)
          if (associated(q1   )) nullify(q1)
          if (associated(q2   )) nullify(q2)
          if (associated(q3   )) nullify(q3)
          if (associated(rbuff)) nullify(rbuff)
          if (associated(cbuff)) nullify(cbuff)

      end subroutine

      subroutine do_scale_splitting_multiple(this,fin,fout,kc,rbuff,cbuff)
          class(stats_xy), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in) :: fin
          real(rkind), dimension(:,:,:,:), intent(out) :: fout
          real(rkind), dimension(:), intent(in) :: kc
          real(rkind), dimension(:,:,:,:), intent(inout) :: rbuff
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuff
          integer :: scl

          ! Copy the input array so if size(kc) = 1, this routine simply copies the data from fin
          fout(:,:,:,1) = fin

          ! First, subtract the mean
          call this%do_scale_splitting_single(fin,rbuff(:,:,:,1),0.d0,1.d-14,rbuff(:,:,:,2),cbuff)
          rbuff(:,:,:,1) = fin - rbuff(:,:,:,1)

          ! Keep the mean in u1:
          !rbuff(:,:,:,1) = fin

          ! Split the scales along the k2d boundaries defined by kc
          do scl = 1,size(kc)-1
              call this%do_scale_splitting_single(rbuff(:,:,:,1),fout(:,:,:,scl),kc(scl),kc(scl+1),rbuff(:,:,:,2),cbuff)
          end do
      end subroutine

      subroutine do_scale_splitting_single(this,fin,fout,kleft,kright,rbuff,cbuff)
          class(stats_xy), intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in)  :: fin
          real(rkind), dimension(:,:,:), intent(out) :: fout  
          real(rkind), intent(in) :: kleft, kright
          real(rkind), dimension(:,:,:), intent(inout) :: rbuff
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuff

          ! Band-pass filter via two low-pass filters
          call this%low_pass_filter(fin,rbuff,kleft ,cbuff)
          call this%low_pass_filter(fin,fout ,kright,cbuff)

          fout = fout - rbuff

      end subroutine

      subroutine low_pass_filter(this,f,ffil,kc,cbuff)
          class(stats_xy), intent(inout) :: this 
          real(rkind), dimension(:,:,:), intent(in) :: f
          real(rkind), dimension(:,:,:), intent(out) :: ffil
          real(rkind), intent(in) :: kc
          complex(rkind), dimension(:,:,:), intent(inout) :: cbuff

          this%mask = 0.d0
          where (this%k2d < kc) this%mask = 1.d0

          call this%spectC%fft(f,cbuff)
          cbuff = this%mask*cbuff
          call this%spectC%ifft(cbuff,ffil)
      end subroutine

      subroutine compute_stats(this)
          class(stats_xy), intent(inout) :: this
          integer :: n, scl, sca

          if ( (this%do_stats)                             .and. &
               (this%sim%step >= this%tid_start          )) then

             ! Step 1: Get budget RHS terms and all other budget terms at t^n
             if (mod(this%sim%step,this%compute_freq) == 0) then
                 call message(1,"Computing stats")

                 ! If we aren't scale-splitting, this routine will simply copy the full fields
                 call this%copy_and_scale_split(this%kc_vec,sca=1,scalar_only=.false.)

                 this%tidx = this%tidx + 1
                 do scl = 1,this%nscales
                     call this%link_pointers(this%tidx,min(1,this%nscalars),scl)

                     ! Get mean profiles
                     call this%mean(this%sim%u  ,this%meanU)
                     call this%mean(this%sim%v  ,this%meanV)
                     call this%mean(this%sim%wC ,this%meanW)

                     ! TODO: This is not the true mean pressure. Need to account for
                     ! meanT component
                     call this%mean(this%sim%pressure,this%meanP) 

                     !   > mean force components
                     if (this%sim%localizedForceLayer == 2) then
                         this%sim%rbuffxC(:,:,:,1) = this%sim%spectForceLayer%ampFact_x*this%sim%spectForceLayer%fx
                         this%sim%rbuffxC(:,:,:,2) = this%sim%spectForceLayer%ampFact_y*this%sim%spectForceLayer%fy
                         call this%sim%spectForceLayer%interpE2C(&
                           this%sim%spectForceLayer%fz, this%sim%rbuffxC(:,:,:,3), &
                           this%sim%rbuffyC(:,:,:,1), this%sim%rbuffzC(:,:,:,1), & 
                           this%sim%rbuffyE(:,:,:,1), this%sim%rbuffzE(:,:,:,1))
                         this%sim%rbuffxC(:,:,:,3) = this%sim%spectForceLayer%ampFact_z*this%sim%rbuffxC(:,:,:,3)

                         call this%mean(this%sim%rbuffxC(:,:,:,1), this%meanFx)
                         call this%mean(this%sim%rbuffxC(:,:,:,2), this%meanFy)
                         call this%mean(this%sim%rbuffxC(:,:,:,3), this%meanFz)

                         ! force variance
                         this%sim%rbuffxC(:,:,:,2) = &
                           this%sim%spectForceLayer%ampFact_x*this%sim%spectForceLayer%ampFact_x*&
                           this%fx*this%fx + & 
                           this%sim%spectForceLayer%ampFact_y*this%sim%spectForceLayer%ampFact_y*&
                           this%fy*this%fy + & 
                           this%sim%spectForceLayer%ampFact_z*this%sim%spectForceLayer%ampFact_z*&
                           this%fz*this%fz
                         call this%mean(this%sim%rbuffxC(:,:,:,2),this%fifi)
                         this%fifi = this%fifi - (this%meanFx*this%meanFx + &
                           this%meanFy*this%meanFy + this%meanFz*this%meanFz)

                         if (this%nscalars>0) then
                             call this%mean(this%sim%spectForceLayer%fT,this%meanfT)
                             call this%covariance(this%fT,this%fT,this%fTfT)
                         end if
                     else
                         this%meanFx = 0.d0
                         this%meanFy = 0.d0
                         this%meanFz = 0.d0
                         this%fifi   = 0.d0
                         if (this%nscalars>0) then
                             this%meanFT = 0.d0
                             this%fTfT = 0.d0
                         end if
                     end if

                     ! Compute MKE
                     this%mke = 0.5d0*(this%meanU*this%meanU + this%meanV*this%meanV + &
                       this%meanW*this%meanW)

                     ! Reynolds stresses
                     !   > Normal stresses
                     call this%covariance(this%u,  this%u,  this%uu)
                     call this%covariance(this%v,  this%v,  this%vv)
                     call this%covariance(this%wC, this%wC, this%ww)

                     !   > Cross stresses
                     call this%covariance(this%u,  this%v,  this%uv)
                     call this%covariance(this%u,  this%wC, this%uw)
                     call this%covariance(this%v,  this%wC, this%vw)

                     ! Pressure covariances
                     call this%covariance(this%pressure, this%pressure, this%pp)
                     call this%covariance(this%u       , this%pressure, this%up)
                     call this%covariance(this%v       , this%pressure, this%vp)
                     call this%covariance(this%wC      , this%pressure, this%wp)

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
                       this%tau11,this%tau12,this%tau13,this%tau22,this%tau23,this%tau33,scl)

                     ! TKE fluxes
                     this%sim%rbuffxC(:,:,:,1) = this%u*this%u + this%v*this%v + this%wC*this%wC
                     call this%covariance(this%u ,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,1))
                     call this%covariance(this%v ,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,2))
                     call this%covariance(this%wC,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,3))

                     this%tke_flux(:,1) = 0.5d0*this%zbuff(:,1)! - &
                       !(this%meanU*this%uu + this%meanV*this%uv + this%meanW*this%uw)
                     this%tke_flux(:,2) = 0.5d0*this%zbuff(:,2)! - &
                       !(this%meanU*this%uv + this%meanV*this%vv + this%meanW*this%vw)
                     this%tke_flux(:,3) = 0.5d0*this%zbuff(:,3)! - &
                       !(this%meanU*this%uw + this%meanV*this%vw + this%meanW*this%ww)

                 end do

                 ! Scalar stats
                 do n = 1,this%nscalars
                     ! Run scale-splitting on each scalar so we don't have to store fields for all scales and all scalars
                     call this%copy_and_scale_split(this%kc_vec,sca=n,scalar_only=.true.)

                     do scl = 1,this%nscales
                         call this%link_pointers(this%tidx,n,scl)
                         
                         ! Mean scalar profile
                         call this%mean(this%T_all, this%meanT)

                         ! Mean scalar gradient
                         call this%mean(this%dTdzC_all,this%dTdz)

                         ! Scalar variance and fluxes
                         call this%covariance(this%T,  this%T, this%TT)
                         call this%covariance(this%u,  this%T, this%uT)
                         call this%covariance(this%v,  this%T, this%vT)
                         call this%covariance(this%wC, this%T, this%wT)

                         ! Scalar graident variances
                         call this%covariance(this%dTdxC,this%dTdxC,this%dTdxj_var(:,1))
                         call this%covariance(this%dTdyC,this%dTdyC,this%dTdxj_var(:,2))
                         call this%covariance(this%dTdzC,this%dTdzC,this%dTdxj_var(:,3))

                         ! SGS heat flux moments
                         call this%covariance(this%q1,this%q1,this%qj_var(:,1))
                         call this%covariance(this%q2,this%q2,this%qj_var(:,2))
                         call this%covariance(this%q3,this%q3,this%qj_var(:,3))
                         call this%mean(this%q3,this%meanq3)

                         call this%get_TT_budget_rhs(scl)
                         call this%get_wT_budget_rhs(scl)

                         this%sim%rbuffxC(:,:,:,1) = this%T*this%T
                         call this%covariance(this%u ,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,1))
                         call this%covariance(this%v ,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,2))
                         call this%covariance(this%wC,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,3))

                         this%sca_var_flux(:,1) = this%zbuff(:,1)! - 2.d0*this%meanT*this%uT
                         this%sca_var_flux(:,2) = this%zbuff(:,2)! - 2.d0*this%meanT*this%vT
                         this%sca_var_flux(:,3) = this%zbuff(:,3)! - 2.d0*this%meanT*this%wT
                     end do
                 end do

                 ! Increment nwrite. This is to handle the edge-case where we
                 ! started the simulation at a tid such that the first dump call is
                 ! before we've filled the stats buffers, so we only want to write
                 ! data that has been computed.
                 this%nwrite = this%nwrite + 1

             end if

             ! Step 2: compute unsteady terms
             ! Step 2.a: If we haven't done scaling splitting yet, do it now
             if ( (mod(this%sim%step+1,this%compute_freq) == 0) .or. &
                  (mod(this%sim%step-1,this%compute_freq) == 0) ) then
                  call this%copy_and_scale_split(this%kc_vec, sca=1,scalar_only=.false.)
             end if
             ! Step 2.b: Compute unsteady terms for uu, vv, and ww
             ! TODO: Need to loop over scales, link pointers, and compute the unsteady term. But we need a dimension in uvarold/new
             ! for the different scales
             do scl = 1,this%nscales
                 call this%link_pointers(max(this%tidx,1),1,scl)
                 call this%compute_unsteady_terms(terms='velocity',sca=1,scl=scl)
             end do

             ! Step 2.c: Loop over scalars, do scale splitting, link pointers, and compute/store unsteady terms
             do sca = 1,this%nscalars
                 call this%copy_and_scale_split(this%kc_vec, sca=sca,scalar_only=.true.)
                 do scl = 1,this%nscales
                     call this%link_pointers(max(this%tidx,1),sca,scl)
                     call this%compute_unsteady_terms(terms='scalar',sca=sca,scl=scl)
                 end do
             end do

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
          tau_11,tau_12,tau_13,tau_22,tau_23,tau_33,scl)
          class(stats_xy), intent(inout) :: this
          real(rkind), dimension(:,:), intent(inout) :: tke_budget, &
            R11_budget, R22_budget, R33_budget
          real(rkind), dimension(:,:,:), intent(in) :: tau_11, tau_12, tau_13, &
                                                       tau_22, tau_23, tau_33
          integer, intent(in) :: scl                                           
          real(rkind), dimension(:,:,:), pointer :: dudx_all, dudy_all, dudz_all, &
                                                    dvdx_all, dvdy_all, dvdz_all, &
                                                    dwdx_all, dwdy_all, dwdz_all
          real(rkind) :: TwobyRe, OnebyRe

          tke_budget = 0.d0
          OnebyRe = 1.d0/this%sim%Re
          TwobyRe = 2.d0*OnebyRe

          associate( &
              dudx_all => this%sim%duidxjC(:,:,:,1), dudy_all => this%sim%duidxjC(:,:,:,2), dudz_all => this%sim%duidxjC(:,:,:,3), &
              dvdx_all => this%sim%duidxjC(:,:,:,4), dvdy_all => this%sim%duidxjC(:,:,:,5), dvdz_all => this%sim%duidxjC(:,:,:,6), &
              dwdx_all => this%sim%duidxjC(:,:,:,7), dwdy_all => this%sim%duidxjC(:,:,:,8), dwdz_all => this%sim%duidxjC(:,:,:,9))

              ! Shear production
              call this%gradient_production(dudx_all,dudy_all,dudz_all,this%uu,this%uv,this%uw,R11_budget(:,2))!,"R11 production calculation")
              call this%gradient_production(dvdx_all,dvdy_all,dvdz_all,this%uv,this%vv,this%vw,R22_budget(:,2))
              call this%gradient_production(dwdx_all,dwdy_all,dwdz_all,this%uw,this%vw,this%ww,R33_budget(:,2))
              tke_budget(:,2) = R11_budget(:,2) + R22_budget(:,2) + R33_budget(:,2)
              R11_budget(:,2) = 2.d0*R11_budget(:,2)
              R22_budget(:,2) = 2.d0*R22_budget(:,2)
              R33_budget(:,2) = 2.d0*R33_budget(:,2)
          end associate

          ! Force production
          call this%covariance(this%fx, this%u,  this%zbuff(:,1))
          call this%covariance(this%fy, this%v,  this%zbuff(:,2))
          call this%covariance(this%fz, this%wC, this%zbuff(:,3))
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
          call this%turbulent_transport(this%u,this%u,&
            this%dudx,this%dudy,this%dudz,this%dudx,this%dudy,this%dudz,&
            R11_budget(:,4),R11_budget(:,5))
          call this%turbulent_transport(this%v,this%v,&
            this%dvdx,this%dvdy,this%dvdz,this%dvdx,this%dvdy,this%dvdz,&
            R22_budget(:,4),R22_budget(:,5))
          call this%turbulent_transport(this%wC,this%wC,&
            this%dwdx,this%dwdy,this%dwdz,this%dwdx,this%dwdy,this%dwdz,&
            R33_budget(:,4),R33_budget(:,5))

          tke_budget(:,5) = 0.5d0*(R11_budget(:,5) + R22_budget(:,5) + R33_budget(:,5))

          ! Pressure transport
          call this%covariance(this%wC, this%pressure, this%zbuff(:,1))
          this%zbuff(:,1) = -1.d0*this%zbuff(:,1)
          call this%ddz(this%zbuff(:,1), tke_budget(:,6), 0, 0)
          R11_budget(:,6) = 0.d0
          R22_budget(:,6) = 0.d0
          R33_budget(:,6) = 2.d0*tke_budget(:,6)

          ! Viscous transport: ddz(<ui*Si3>)
          this%sim%rbuffxC(:,:,:,1) = 0.5d0*(this%dudz + this%dwdx)
          this%sim%rbuffxC(:,:,:,2) = 0.5d0*(this%dvdz + this%dwdy)
          call this%covariance(this%u ,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,1))
          call this%covariance(this%v ,this%sim%rbuffxC(:,:,:,2),this%zbuff(:,2))
          call this%covariance(this%wC,this%dwdz                ,this%zbuff(:,3))
          this%zbuff(:,4) = TwobyRe*(this%zbuff(:,1) + this%zbuff(:,2) + this%zbuff(:,3))
          call this%ddz(this%zbuff(:,4),tke_budget(:,7),0,0)

          call this%molecular_transport(this%u ,this%dudz,OnebyRe,R11_budget(:,7))
          call this%molecular_transport(this%v ,this%dvdz,OnebyRe,R22_budget(:,7))
          call this%molecular_transport(this%wC,this%dwdz,OnebyRe,R33_budget(:,7))
          tke_budget(:,13) = R11_budget(:,7) + R22_budget(:,7) + R33_budget(:,7)
          R11_budget(:,7) = 2.d0*R11_budget(:,7)
          R22_budget(:,7) = 2.d0*R22_budget(:,7)
          R33_budget(:,7) = 2.d0*R33_budget(:,7)

          ! Viscous dissipation
          !   > NOTE: this is positive. Will need to plot
          !     negative of this for budget visualization
          this%sim%rbuffxC(:,:,:,3) = 0.5d0*(this%dudy + this%dvdx) ! S12
          call this%covariance(this%dudx,this%dudx,this%zbuff(:,1)) ! <S11'*S11'>
          call this%covariance(this%dvdy,this%dvdy,this%zbuff(:,2)) ! <S22'*S22'>
          call this%covariance(this%dwdz,this%dwdz,this%zbuff(:,3)) ! <S33'*S33'>
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

          call this%pseudo_dissipation(this%dudx,this%dudy,this%dudz,this%dudx,&
            this%dudy,this%dudz,TwobyRe,R11_budget(:,8))
          call this%pseudo_dissipation(this%dvdx,this%dvdy,this%dvdz,this%dvdx,&
            this%dvdy,this%dvdz,TwobyRe,R22_budget(:,8))
          call this%pseudo_dissipation(this%dwdx,this%dwdy,this%dwdz,this%dwdx,&
            this%dwdy,this%dwdz,TwobyRe,R33_budget(:,8))
          tke_budget(:,12) = 0.5d0*(R11_budget(:,8) + R22_budget(:,8) + R33_budget(:,8))

          ! Compute velocity gradient variances
          call this%covariance(this%dudx,this%dudx,this%duidxj_var(:,1))
          call this%covariance(this%dudy,this%dudy,this%duidxj_var(:,2))
          call this%covariance(this%dudz,this%dudz,this%duidxj_var(:,3))

          call this%covariance(this%dvdx,this%dvdx,this%duidxj_var(:,4))
          call this%covariance(this%dvdy,this%dvdy,this%duidxj_var(:,5))
          call this%covariance(this%dvdz,this%dvdz,this%duidxj_var(:,6))

          call this%covariance(this%dwdx,this%dwdx,this%duidxj_var(:,7))
          call this%covariance(this%dwdy,this%dwdy,this%duidxj_var(:,8))
          call this%covariance(this%dwdz,this%dwdz,this%duidxj_var(:,9))

          call this%covariance(this%sim%rbuffxC(:,:,:,3),this%sim%rbuffxC(:,:,:,3),this%S12_var)
          call this%covariance(this%sim%rbuffxC(:,:,:,1),this%sim%rbuffxC(:,:,:,1),this%S13_var)
          call this%covariance(this%sim%rbuffxC(:,:,:,2),this%sim%rbuffxC(:,:,:,2),this%S23_var)

          ! SGS transport
          call this%SGS_transport(this%u ,this%tau13,R11_budget(:,9))
          call this%SGS_transport(this%v ,this%tau23,R22_budget(:,9))
          call this%SGS_transport(this%wC,this%tau33,R33_budget(:,9))
          tke_budget(:,9) = R11_budget(:,9) + R22_budget(:,9) + R33_budget(:,9)
          R11_budget(:,9) = 2.d0*R11_budget(:,9)
          R22_budget(:,9) = 2.d0*R22_budget(:,9)
          R33_budget(:,9) = 2.d0*R33_budget(:,9)

          ! SGS dissipation
          !   > NOTE: consistent with viscous dissipation, the sign on this
          !   is opposite of what needs to be plotted so that epsilon_total
          !   = epsilon + epsilon_SGS and -epsilon_total is on RHS in budget
          call this%SGS_dissipation(this%dudx,this%dudy,this%dudz,this%tau11,this%tau12,this%tau13,R11_budget(:,10))
          call this%SGS_dissipation(this%dvdx,this%dvdy,this%dvdz,this%tau12,this%tau22,this%tau23,R22_budget(:,10))
          call this%SGS_dissipation(this%dwdx,this%dwdy,this%dwdz,this%tau13,this%tau23,this%tau33,R33_budget(:,10))
          tke_budget(:,10) = R11_budget(:,10) + R22_budget(:,10) + R33_budget(:,10)
          R11_budget(:,10) = 2.d0*R11_budget(:,10)
          R22_budget(:,10) = 2.d0*R22_budget(:,10)
          R33_budget(:,10) = 2.d0*R33_budget(:,10)
          
          ! Buoyancy transfer
          if (this%sim%isStratified) then
              call this%covariance(this%wC, this%T, tke_budget(:,11))
              tke_budget(:,11) = -this%sim%buoyancyFact*tke_budget(:,11)
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
          call this%covariance(this%pressure,this%dudx,R11_budget(:,12))
          call this%covariance(this%pressure,this%dvdy,R22_budget(:,12))
          call this%covariance(this%pressure,this%dwdz,R33_budget(:,12))
          R11_budget(:,12) = 2.d0*R11_budget(:,12)
          R22_budget(:,12) = 2.d0*R22_budget(:,12)
          R33_budget(:,12) = 2.d0*R33_budget(:,12)

          ! Interscale transfer
          call this%interscale_transfer(this%ujsplit(:,:,:,:,1),this%ujsplit(:,:,:,:,1),&
            this%ujsplit,this%duidxj(:,:,:,:,1:3),this%duidxj(:,:,:,:,1:3),scl,&
            R11_budget(:,13),R11_budget(:,14))
          call this%interscale_transfer(this%ujsplit(:,:,:,:,2),this%ujsplit(:,:,:,:,2),&
            this%ujsplit,this%duidxj(:,:,:,:,4:6),this%duidxj(:,:,:,:,4:6),scl,&
            R22_budget(:,13),R22_budget(:,14))
          call this%interscale_transfer(this%ujsplit(:,:,:,:,3),this%ujsplit(:,:,:,:,3),&
            this%ujsplit,this%duidxj(:,:,:,:,7:9),this%duidxj(:,:,:,:,7:9),scl,&
            R33_budget(:,13),R33_budget(:,14))
          R11_budget(:,13:14) = 0.5d0*R11_budget(:,13:14)
          R22_budget(:,13:14) = 0.5d0*R22_budget(:,13:14)
          R33_budget(:,13:14) = 0.5d0*R33_budget(:,13:14)
          tke_budget(:,14) = 0.5d0*(R11_budget(:,13) + R22_budget(:,13) + R33_budget(:,13))
          tke_budget(:,15) = 0.5d0*(R11_budget(:,14) + R22_budget(:,14) + R33_budget(:,14))

          ! Mixed scale turbulent transport
          call this%mixed_scale_turbulent_transport(this%ujsplit(:,:,:,:,1),this%ujsplit(:,:,:,:,1),this%ujsplit,scl,R11_budget(:,15))
          call this%mixed_scale_turbulent_transport(this%ujsplit(:,:,:,:,2),this%ujsplit(:,:,:,:,2),this%ujsplit,scl,R22_budget(:,15))
          call this%mixed_scale_turbulent_transport(this%ujsplit(:,:,:,:,3),this%ujsplit(:,:,:,:,3),this%ujsplit,scl,R33_budget(:,15))
          R11_budget(:,15) = R11_budget(:,15)
          R22_budget(:,15) = R22_budget(:,15)
          R33_budget(:,15) = R33_budget(:,15)
          tke_budget(:,16) = 0.5d0*(R11_budget(:,15) + R22_budget(:,15) + R33_budget(:,15))


      end subroutine

      subroutine interscale_transfer(this,f,g,uj,dfdxj,dgdxj,s,I_d,I_id)
          class(stats_xy), intent(inout) :: this
          real(rkind), dimension(:,:,:,:), intent(in) :: f, g
          real(rkind), dimension(:,:,:,:,:), intent(in) :: uj, dfdxj, dgdxj
          integer, intent(in) :: s
          real(rkind), dimension(:), intent(out) :: I_d, I_id
          integer :: sig, gam, nscales, j
          real(rkind), dimension(:,:,:), pointer :: rbuff
          ! Outputs:
          !   I_d  --> Interscale transfer (Direct): transfer between scales sigma and s
          !   I_id --> Interscale transfer (InDirect): transfer between scales sigma, gamma, and s

          nscales = size(f,4)

          associate(rbuff => this%sim%rbuffxC(:,:,:,1))
              ! "Direct" interscale transfer
              rbuff = 0.d0
              do sig = 1,nscales
                  if (sig .ne. s) then
                      do j = 1,3
                        rbuff = rbuff + f(:,:,:,s  )*uj(:,:,:,s  ,j)*dgdxj(:,:,:,sig,j) + &
                                        g(:,:,:,s  )*uj(:,:,:,s  ,j)*dfdxj(:,:,:,sig,j) - &
                                        f(:,:,:,sig)*uj(:,:,:,sig,j)*dgdxj(:,:,:,s  ,j) - &
                                        g(:,:,:,sig)*uj(:,:,:,sig,j)*dfdxj(:,:,:,s  ,j)

                      end do
                  end if
              end do
              call this%mean(rbuff,I_d)

              ! "Indirect" interscale transfer
              rbuff = 0.d0
              do sig = 1,nscales
                  if (sig .ne. s) then
                      do gam = 1,nscales
                          if ((gam .ne. s) .and. (gam .ne. sig)) then
                              do j = 1,3
                                  rbuff = rbuff + g(:,:,:,sig)*uj(:,:,:,gam,j)*dfdxj(:,:,:,s,j) + &
                                                  f(:,:,:,sig)*uj(:,:,:,gam,j)*dgdxj(:,:,:,s,j)
                              end do
                          end if
                      end do
                  end if
              end do
              call this%mean(rbuff,I_id)
          end associate

      end subroutine

      subroutine mixed_scale_turbulent_transport(this,f,g,uj,s,budget)
          class(stats_xy), intent(inout) :: this
          real(rkind), dimension(:,:,:,:), intent(in) :: f, g
          real(rkind), dimension(:,:,:,:,:), intent(in) :: uj
          integer, intent(in) :: s
          real(rkind), dimension(:), intent(out) :: budget
          integer :: sig, gam, j, nscales
          real(rkind), dimension(:,:,:), pointer :: rbuff

          nscales = size(f,4)

          associate(rbuff => this%sim%rbuffxC(:,:,:,1))
              rbuff = 0.d0
              do sig = 1,nscales
                  if (sig .ne. s) then
                      do gam = 1,nscales
                          if (gam .ne. s) then
                              rbuff = rbuff - &
                                f(:,:,:,s)*g(:,:,:,sig)*uj(:,:,:,gam,3) - &
                                g(:,:,:,s)*f(:,:,:,sig)*uj(:,:,:,gam,3)
                          end if
                      end do
                  end if
              end do

              call this%mean(rbuff,this%zbuff(:,1))
              call this%ddz(this%zbuff(:,1),budget,0,0)

          end associate
      end subroutine

      subroutine get_TT_budget_rhs(this,scl)
          class(stats_xy), intent(inout) :: this
          integer, intent(in) :: scl
          real(rkind) :: OnebyRePr
          real(rkind), dimension(:,:,:), pointer :: dTdx, dTdy, dTdz
          
          this%TT_budget = 0.d0
          OnebyRePr = 1.d0/(this%sim%Re*this%sim%PrandtlFluid)

          ! Gradient production
          call this%gradient_production(this%dTdxC_all,this%dTdyC_all,this%dTdzC_all,&
            this%uT,this%vT,this%wT,this%TT_budget(:,2))
          this%TT_budget(:,2) = 2.d0*this%TT_budget(:,2)

          ! Convective transport
          call this%convective_transport(this%TT,this%TT_budget(:,3))

          ! Turbulent transport
          call this%turbulent_transport(this%T,this%T,this%dTdxC,this%dTdyC,this%dTdzC,&
            this%dTdxC,this%dTdyC,this%dTdzC,this%TT_budget(:,3), this%TT_budget(:,4))

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

          ! Interscale transfer
          call this%interscale_transfer(this%Tsplit,this%Tsplit,&
            this%ujsplit,this%dTdxj,this%dTdxj,scl,&
            this%TT_budget(:,9),this%TT_budget(:,10))
          this%TT_budget(:,9:10) = 0.5d0*this%TT_budget(:,9:10)

          ! Mixed scale turbulent transport
          call this%mixed_scale_turbulent_transport(this%Tsplit,this%Tsplit,this%ujsplit,scl,this%TT_budget(:,11))
          this%TT_budget(:,11) = this%TT_budget(:,11)

          ! Force production
          call this%covariance(this%T,this%fT,this%TT_budget(:,12))


      end subroutine

      subroutine get_wT_budget_rhs(this,scl)
          class(stats_xy), intent(inout) :: this
          integer, intent(in) :: scl
          real(rkind) :: OnebyRePr, OnebyRe, OnebyFr2, viscosity
          real(rkind), dimension(:,:,:), pointer :: dwdx, dwdy, dwdz

          this%wT_budget = 0.d0
          OnebyRePr = 1.d0/(this%sim%Re*this%sim%PrandtlFluid)
          OnebyRe = 1.d0/this%sim%Re
          OnebyFr2 = this%sim%buoyancyFact

          ! Scalar gradient production
          call this%gradient_production(this%dTdxC_all,this%dTdyC_all,this%dTdzC_all,this%uw,this%vw,this%ww,this%wT_budget(:,2))

          ! Shear gradient production
          associate(dwdx_all => this%sim%duidxjC(:,:,:,7), dwdy_all => this%sim%duidxjC(:,:,:,8), dwdz_all => this%sim%duidxjC(:,:,:,9))
              call this%gradient_production(dwdx_all,dwdy_all,dwdz_all,this%uT,this%vT,this%wT,this%wT_budget(:,3))
          end associate

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
          call this%turbulent_transport(this%wC,this%T,this%dwdx,this%dwdy,this%dwdz,&
            this%dTdxC,this%dTdyC,this%dTdzC,this%wT_budget(:,6), this%wT_budget(:,7))

          ! Pressure transport
          call this%covariance(this%T,this%sim%pressure,this%zbuff(:,1))
          call this%ddz(this%zbuff(:,1),this%wT_budget(:,8),0,0)
          this%wT_budget(:,8) = -1.d0*this%wT_budget(:,8)

          ! Molecular transport
          call this%molecular_transport(this%sim%wC,this%dTdzC,OnebyRePr,this%wT_budget(:,9 ))
          call this%molecular_transport(this%T     ,this%dwdz ,OnebyRe  ,this%wT_budget(:,10))

          ! SGS transport
          call this%SGS_transport(this%sim%wC,this%q3   ,this%wT_budget(:,11))
          call this%SGS_transport(this%T     ,this%tau33,this%wT_budget(:,12))

          ! Pressure scrambling
          call this%covariance(this%sim%pressure, this%dTdzC, this%wT_budget(:,13))

          ! Molecular destruction
          call this%pseudo_dissipation(this%dwdx,this%dwdy,this%dwdz,this%dTdxC,this%dTdyC,this%dTdzC,OnebyRePr,this%wT_budget(:,14))
          call this%pseudo_dissipation(this%dwdx,this%dwdy,this%dwdz,this%dTdxC,this%dTdyC,this%dTdzC,OnebyRe  ,this%wT_budget(:,15))

          ! SGS destruction
          call this%SGS_dissipation(this%dwdx,this%dwdy,this%dwdz,this%q1,this%q2,this%q3,this%wT_budget(:,16))
          call this%SGS_dissipation(this%dTdxC,this%dTdyC,this%dTdzC,this%tau13,this%tau23,this%tau33,this%wT_budget(:,17))

          ! Interscale transfer
          call this%interscale_transfer(this%Tsplit,this%ujsplit(:,:,:,:,3),&
            this%ujsplit,this%dTdxj,this%duidxj(:,:,:,:,7:9),scl,&
            this%wT_budget(:,18),this%wT_budget(:,19))

          ! Mixed scale turbulent transport
          call this%mixed_scale_turbulent_transport(this%Tsplit,this%ujsplit(:,:,:,:,3),this%ujsplit,scl,this%wT_budget(:,20))

          ! Scalar force production
          call this%covariance(this%wC,this%fT,this%wT_budget(:,21))
      end subroutine

      subroutine link_pointers(this,tidx,sca,scl)
        class(stats_xy), intent(inout), target :: this
        integer, intent(in) :: tidx, sca, scl
        integer :: id

        id = 16
        this%tke_budget => this%stats(:,1:id    ,tidx,scl); id = id + 1
        this%tke        => this%stats(:,id      ,tidx,scl); id = id + 1
        this%fifi       => this%stats(:,id      ,tidx,scl); id = id + 1
        this%meanU      => this%stats(:,id      ,tidx,scl); id = id + 1
        this%meanV      => this%stats(:,id      ,tidx,scl); id = id + 1
        this%meanW      => this%stats(:,id      ,tidx,scl); id = id + 1
        this%meanP      => this%stats(:,id      ,tidx,scl); id = id + 1
        this%uu         => this%stats(:,id      ,tidx,scl); id = id + 1
        this%vv         => this%stats(:,id      ,tidx,scl); id = id + 1
        this%ww         => this%stats(:,id      ,tidx,scl); id = id + 1
        this%uv         => this%stats(:,id      ,tidx,scl); id = id + 1
        this%uw         => this%stats(:,id      ,tidx,scl); id = id + 1
        this%vw         => this%stats(:,id      ,tidx,scl); id = id + 1
        this%pp         => this%stats(:,id      ,tidx,scl); id = id + 1
        this%up         => this%stats(:,id      ,tidx,scl); id = id + 1
        this%vp         => this%stats(:,id      ,tidx,scl); id = id + 1
        this%wp         => this%stats(:,id      ,tidx,scl); id = id + 1
        this%meanFx     => this%stats(:,id      ,tidx,scl); id = id + 1
        this%meanFy     => this%stats(:,id      ,tidx,scl); id = id + 1
        this%meanFz     => this%stats(:,id      ,tidx,scl); id = id + 1
        this%mke        => this%stats(:,id      ,tidx,scl); id = id + 1
        this%R11_budget => this%stats(:,id:id+14,tidx,scl); id = id + 15
        this%R22_budget => this%stats(:,id:id+14,tidx,scl); id = id + 15
        this%R33_budget => this%stats(:,id:id+14,tidx,scl); id = id + 15
        this%duidxj_var => this%stats(:,id:id+8 ,tidx,scl); id = id + 9
        this%S12_var    => this%stats(:,id      ,tidx,scl); id = id + 1
        this%S13_var    => this%stats(:,id      ,tidx,scl); id = id + 1
        this%S23_var    => this%stats(:,id      ,tidx,scl); id = id + 1
        this%tauij_var  => this%stats(:,id:id+5 ,tidx,scl); id = id + 6
        this%tke_flux   => this%stats(:,id:id+2 ,tidx,scl); id = id + 3

        this%tau11 => this%tauij(:,:,:,scl,1)
        this%tau12 => this%tauij(:,:,:,scl,2)
        this%tau13 => this%tauij(:,:,:,scl,3)
        this%tau22 => this%tauij(:,:,:,scl,4)
        this%tau23 => this%tauij(:,:,:,scl,5)
        this%tau33 => this%tauij(:,:,:,scl,6)

        this%dudx => this%duidxj(:,:,:,scl,1)
        this%dudy => this%duidxj(:,:,:,scl,2)
        this%dudz => this%duidxj(:,:,:,scl,3)
        this%dvdx => this%duidxj(:,:,:,scl,4)
        this%dvdy => this%duidxj(:,:,:,scl,5)
        this%dvdz => this%duidxj(:,:,:,scl,6)
        this%dwdx => this%duidxj(:,:,:,scl,7)
        this%dwdy => this%duidxj(:,:,:,scl,8)
        this%dwdz => this%duidxj(:,:,:,scl,9)

        this%u        => this%ujsplit(:,:,:,scl,1)
        this%v        => this%ujsplit(:,:,:,scl,2)
        this%wC       => this%ujsplit(:,:,:,scl,3)
        this%pressure => this%psplit(:,:,:,scl)

        this%fx => this%fxsplit(:,:,:,scl)
        this%fy => this%fysplit(:,:,:,scl)
        this%fz => this%fzsplit(:,:,:,scl)

        if (this%nscalars > 0) then
            id = 1
            this%TT_budget    => this%stats_sca(:,id:id+11,tidx,sca,scl); id = id + 12
            this%dTdz         => this%stats_sca(:,id      ,tidx,sca,scl); id = id + 1
            this%meanT        => this%stats_sca(:,id      ,tidx,sca,scl); id = id + 1
            this%meanfT       => this%stats_sca(:,id      ,tidx,sca,scl); id = id + 1
            this%fTfT         => this%stats_sca(:,id      ,tidx,sca,scl); id = id + 1
            this%uT           => this%stats_sca(:,id      ,tidx,sca,scl); id = id + 1
            this%vT           => this%stats_sca(:,id      ,tidx,sca,scl); id = id + 1
            this%wT           => this%stats_sca(:,id      ,tidx,sca,scl); id = id + 1
            this%TT           => this%stats_sca(:,id      ,tidx,sca,scl); id = id + 1
            this%wT_budget    => this%stats_sca(:,id:id+20,tidx,sca,scl); id = id + 21
            this%dTdxj_var    => this%stats_sca(:,id:id+2 ,tidx,sca,scl); id = id + 3
            this%qj_var       => this%stats_sca(:,id:id+2 ,tidx,sca,scl); id = id + 3
            this%meanq3       => this%stats_sca(:,id      ,tidx,sca,scl); id = id + 1
            this%sca_var_flux => this%stats_sca(:,id:id+2 ,tidx,sca,scl); id = id + 3 

            this%T => this%Tsplit(:,:,:,scl)
            this%q1 => this%qjsplit(:,:,:,scl,1)
            this%q2 => this%qjsplit(:,:,:,scl,2)
            this%q3 => this%qjsplit(:,:,:,scl,3)
            this%dTdxC => this%dTdxj(:,:,:,scl,1)
            this%dTdyC => this%dTdxj(:,:,:,scl,2)
            this%dTdzC => this%dTdxj(:,:,:,scl,3)
            this%fT => this%fTsplit(:,:,:,scl)
            if (this%sim%isStratified) then
                if (sca == 1) then
                    this%T_all     => this%sim%T
                    this%dTdxC_all => this%sim%dTdxC
                    this%dTdyC_all => this%sim%dTdyC
                    this%dTdzC_all => this%sim%dTdzC
                else
                    this%T_all     => this%sim%scalars(sca-1)%F
                    this%dTdxC_all => this%sim%scalars(sca-1)%dFdxC
                    this%dTdyC_all => this%sim%scalars(sca-1)%dFdyC
                    this%dTdzC_all => this%sim%scalars(sca-1)%dFdzC
                end if
            else
                this%T_all     => this%sim%scalars(sca)%F
                this%dTdxC_all => this%sim%scalars(sca)%dFdxC
                this%dTdyC_all => this%sim%scalars(sca)%dFdyC
                this%dTdzC_all => this%sim%scalars(sca)%dFdzC
            end if
        end if

      end subroutine
      
      subroutine dump_stats(this)
          class(stats_xy), intent(inout) :: this
          integer :: n, sca, scl
          character(len=clen) :: tempname, fname

          if (nrank == 0) then
              do scl = 1,this%nscales
                  do n = 1,this%nwrite
                      write(tempname,"(A3,I2.2,A11,I6.6,A6,I2.2,A4)")&
                        "Run",this%sim%runID,"_stats_xy_t",this%tid(n),"_scale",scl,".out"
                      fname = trim(this%outputdir)//"/"//trim(tempname)
                      call write_2d_ascii(this%stats(:,:,n,scl),trim(fname))
                      do sca = 1,this%nscalars
                          write(tempname,"(A3,I2.2,A13,I2.2,A2,I6.6,A6,I2.2,A4)")&
                            "Run",this%sim%runID,"_stats_xy_sca",sca,"_t",this%tid(n),"_scale",scl,".out"
                          fname = trim(this%outputdir)//"/"//trim(tempname)
                          call write_2d_ascii(this%stats_sca(:,:,n,sca,scl),trim(fname))
                      end do
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
              write(fid,"(A4,I3,A15)") "  > ",id,". Unsteady term";                     id = id + 1
              write(fid,"(A4,I3,A18)") "  > ",id,". Shear production";                  id = id + 1
              write(fid,"(A4,I3,A18)") "  > ",id,". Force production";                  id = id + 1
              write(fid,"(A4,I3,A22)") "  > ",id,". Convective transport";              id = id + 1
              write(fid,"(A4,I3,A21)") "  > ",id,". Turbulent transport";               id = id + 1
              write(fid,"(A4,I3,A20)") "  > ",id,". Pressure transport";                id = id + 1
              write(fid,"(A4,I3,A19)") "  > ",id,". Viscous transport";                 id = id + 1
              write(fid,"(A4,I3,A21)") "  > ",id,". Viscous dissipation";               id = id + 1
              write(fid,"(A4,I3,A15)") "  > ",id,". SGS transport";                     id = id + 1
              write(fid,"(A4,I3,A17)") "  > ",id,". SGS dissipation";                   id = id + 1
              write(fid,"(A4,I3,A19)") "  > ",id,". Buoyancy transfer";                 id = id + 1
              write(fid,"(A4,I3,A20)") "  > ",id,". Pseudo dissipation";                id = id + 1
              write(fid,"(A4,I3,A26)") "  > ",id,". Pseudo viscous transport";          id = id + 1
              write(fid,"(A4,I3,A30)") "  > ",id,". Interscale transfer (direct)";      id = id + 1
              write(fid,"(A4,I3,A32)") "  > ",id,". Interscale transfer (indirect)";    id = id + 1
              write(fid,"(A4,I3,A35)") "  > ",id,". Turbulent transport (mixed scale)"; id = id + 1
              write(fid,*) " "
              write(fid,*) "Misc terms:"
              write(fid,"(A4,I3,A33)") "  > ",id,". Turbulence kinetic energy (TKE)"; id = id + 1
              write(fid,"(A4,I3,A16)") "  > ",id,". Force variance";                  id = id + 1
              write(fid,"(A4,I3,A17)") "  > ",id,". Mean u-velocity";                 id = id + 1
              write(fid,"(A4,I3,A17)") "  > ",id,". Mean v-velocity";                 id = id + 1 
              write(fid,"(A4,I3,A17)") "  > ",id,". Mean w-velocity";                 id = id + 1 
              write(fid,"(A4,I3,A15)") "  > ",id,". Mean pressure";                   id = id + 1 
              write(fid,"(A4,I3,A8)" ) "  > ",id,". <u'u'>";                          id = id + 1 
              write(fid,"(A4,I3,A8)" ) "  > ",id,". <v'v'>";                          id = id + 1 
              write(fid,"(A4,I3,A8)" ) "  > ",id,". <w'w'>";                          id = id + 1 
              write(fid,"(A4,I3,A8)" ) "  > ",id,". <u'v'>";                          id = id + 1 
              write(fid,"(A4,I3,A8)" ) "  > ",id,". <u'w'>";                          id = id + 1 
              write(fid,"(A4,I3,A8)" ) "  > ",id,". <v'w'>";                          id = id + 1 
              write(fid,"(A4,I3,A8)" ) "  > ",id,". <p'p'>";                          id = id + 1 
              write(fid,"(A4,I3,A8)" ) "  > ",id,". <u'p'>";                          id = id + 1 
              write(fid,"(A4,I3,A8)" ) "  > ",id,". <v'p'>";                          id = id + 1 
              write(fid,"(A4,I3,A8)" ) "  > ",id,". <w'p'>";                          id = id + 1 
              write(fid,"(A4,I3,A14)") "  > ",id,". Mean x-force";                    id = id + 1 
              write(fid,"(A4,I3,A14)") "  > ",id,". Mean y-force";                    id = id + 1 
              write(fid,"(A4,I3,A14)") "  > ",id,". Mean z-force";                    id = id + 1 
              write(fid,"(A4,I3,A27)") "  > ",id,". Mean kinetic energy (MKE)";       id = id + 1 
              write(fid,*) " "
              write(fid,*) "Rij budgets:"
              write(fid,*) "  > R11 budget"
              write(fid,"(A6,I3,A15)") "    > ",id,". Unsteady term";                     id = id + 1
              write(fid,"(A6,I3,A18)") "    > ",id,". Shear production";                  id = id + 1
              write(fid,"(A6,I3,A18)") "    > ",id,". Force production";                  id = id + 1
              write(fid,"(A6,I3,A22)") "    > ",id,". Convective transport";              id = id + 1
              write(fid,"(A6,I3,A21)") "    > ",id,". Turbulent transport";               id = id + 1
              write(fid,"(A6,I3,A20)") "    > ",id,". Pressure transport";                id = id + 1
              write(fid,"(A6,I3,A19)") "    > ",id,". Viscous transport";                 id = id + 1
              write(fid,"(A6,I3,A21)") "    > ",id,". Viscous dissipation";               id = id + 1
              write(fid,"(A6,I3,A15)") "    > ",id,". SGS transport";                     id = id + 1
              write(fid,"(A6,I3,A17)") "    > ",id,". SGS dissipation";                   id = id + 1
              write(fid,"(A6,I3,A19)") "    > ",id,". Buoyancy transfer";                 id = id + 1
              write(fid,"(A6,I3,A29)") "    > ",id,". Pressure-strain correlation";       id = id + 1
              write(fid,"(A6,I3,A30)") "    > ",id,". Interscale transfer (direct)";      id = id + 1
              write(fid,"(A6,I3,A32)") "    > ",id,". Interscale transfer (indirect)";    id = id + 1
              write(fid,"(A6,I3,A35)") "    > ",id,". Turbulent transport (mixed scale)"; id = id + 1
              write(fid,*) " "
              write(fid,*) "  > R22 budget"
              write(fid,"(A6,I3,A15)") "    > ",id,". Unsteady term";                     id = id + 1
              write(fid,"(A6,I3,A18)") "    > ",id,". Shear production";                  id = id + 1
              write(fid,"(A6,I3,A18)") "    > ",id,". Force production";                  id = id + 1
              write(fid,"(A6,I3,A22)") "    > ",id,". Convective transport";              id = id + 1
              write(fid,"(A6,I3,A21)") "    > ",id,". Turbulent transport";               id = id + 1
              write(fid,"(A6,I3,A20)") "    > ",id,". Pressure transport";                id = id + 1
              write(fid,"(A6,I3,A19)") "    > ",id,". Viscous transport";                 id = id + 1
              write(fid,"(A6,I3,A21)") "    > ",id,". Viscous dissipation";               id = id + 1
              write(fid,"(A6,I3,A15)") "    > ",id,". SGS transport";                     id = id + 1
              write(fid,"(A6,I3,A17)") "    > ",id,". SGS dissipation";                   id = id + 1
              write(fid,"(A6,I3,A19)") "    > ",id,". Buoyancy transfer";                 id = id + 1
              write(fid,"(A6,I3,A29)") "    > ",id,". Pressure-strain correlation";       id = id + 1
              write(fid,"(A6,I3,A30)") "    > ",id,". Interscale transfer (direct)";      id = id + 1
              write(fid,"(A6,I3,A32)") "    > ",id,". Interscale transfer (indirect)";    id = id + 1
              write(fid,"(A6,I3,A35)") "    > ",id,". Turbulent transport (mixed scale)"; id = id + 1
              write(fid,*) " "
              write(fid,*) "  > R33 budget"
              write(fid,"(A6,I3,A15)") "    > ",id,". Unsteady term";                     id = id + 1
              write(fid,"(A6,I3,A18)") "    > ",id,". Shear production";                  id = id + 1
              write(fid,"(A6,I3,A18)") "    > ",id,". Force production";                  id = id + 1
              write(fid,"(A6,I3,A22)") "    > ",id,". Convective transport";              id = id + 1
              write(fid,"(A6,I3,A21)") "    > ",id,". Turbulent transport";               id = id + 1
              write(fid,"(A6,I3,A20)") "    > ",id,". Pressure transport";                id = id + 1
              write(fid,"(A6,I3,A19)") "    > ",id,". Viscous transport";                 id = id + 1
              write(fid,"(A6,I3,A21)") "    > ",id,". Viscous dissipation";               id = id + 1
              write(fid,"(A6,I3,A15)") "    > ",id,". SGS transport";                     id = id + 1
              write(fid,"(A6,I3,A17)") "    > ",id,". SGS dissipation";                   id = id + 1
              write(fid,"(A6,I3,A19)") "    > ",id,". Buoyancy transfer";                 id = id + 1
              write(fid,"(A6,I3,A29)") "    > ",id,". Pressure-strain correlation";       id = id + 1
              write(fid,"(A6,I3,A30)") "    > ",id,". Interscale transfer (direct)";      id = id + 1
              write(fid,"(A6,I3,A32)") "    > ",id,". Interscale transfer (indirect)";    id = id + 1
              write(fid,"(A6,I3,A35)") "    > ",id,". Turbulent transport (mixed scale)"; id = id + 1
              write(fid,*) " "
              write(fid,*) "Velocity gradient variances"
              write(fid,"(A4,I3,A11)") "  > ",id,". Var(dudx)"; id = id + 1
              write(fid,"(A4,I3,A11)") "  > ",id,". Var(dudy)"; id = id + 1
              write(fid,"(A4,I3,A11)") "  > ",id,". Var(dudz)"; id = id + 1
              write(fid,"(A4,I3,A11)") "  > ",id,". Var(dvdx)"; id = id + 1
              write(fid,"(A4,I3,A11)") "  > ",id,". Var(dvdy)"; id = id + 1
              write(fid,"(A4,I3,A11)") "  > ",id,". Var(dvdz)"; id = id + 1
              write(fid,"(A4,I3,A11)") "  > ",id,". Var(dwdx)"; id = id + 1
              write(fid,"(A4,I3,A11)") "  > ",id,". Var(dwdy)"; id = id + 1
              write(fid,"(A4,I3,A11)") "  > ",id,". Var(dwdz)"; id = id + 1
              write(fid,"(A4,I3,A10)") "  > ",id,". Var(S12)" ; id = id + 1
              write(fid,"(A4,I3,A10)") "  > ",id,". Var(S13)" ; id = id + 1
              write(fid,"(A4,I3,A10)") "  > ",id,". Var(S23)" ; id = id + 1
              write(fid,*) " "
              write(fid,*) "SGS stress tensor variances"
              write(fid,"(A4,I3,A12)") "  > ",id,". Var(tau11)"; id = id + 1
              write(fid,"(A4,I3,A12)") "  > ",id,". Var(tau12)"; id = id + 1
              write(fid,"(A4,I3,A12)") "  > ",id,". Var(tau13)"; id = id + 1
              write(fid,"(A4,I3,A12)") "  > ",id,". Var(tau22)"; id = id + 1
              write(fid,"(A4,I3,A12)") "  > ",id,". Var(tau23)"; id = id + 1
              write(fid,"(A4,I3,A12)") "  > ",id,". Var(tau33)"; id = id + 1
              write(fid,*) " "
              write(fid,*) "TKE fluxes"
              write(fid,"(A4,I3,A14)") "  > ",id,". <u'ui'ui'>/2"; id = id + 1
              write(fid,"(A4,I3,A14)") "  > ",id,". <v'ui'ui'>/2"; id = id + 1
              write(fid,"(A4,I3,A14)") "  > ",id,". <w'ui'ui'>/2"; id = id + 1
              
              close(fid)
          end if
          call MPI_Barrier(MPI_COMM_WORLD,ierr)

          if (this%sim%isStratified .or. this%sim%useScalars) then
              write(tempname,"(A3,I2.2,A27)")"Run",this%sim%runID,"_scalar_stats_index_key.out"
              fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)

              if (nrank == 0) then
                  id = 1
                  open(fid, file=fname)
                  write(fid,*) "Scalar variance budget terms"
                  write(fid,"(A4,I2,A15)") "  > ",id,". Unsteady term";                     id = id + 1
                  write(fid,"(A4,I2,A26)") "  > ",id,". Mean gradient production";          id = id + 1
                  write(fid,"(A4,I2,A22)") "  > ",id,". Convective transport";              id = id + 1
                  write(fid,"(A4,I2,A21)") "  > ",id,". Turbulent transport";               id = id + 1
                  write(fid,"(A4,I2,A21)") "  > ",id,". Molecular transport";               id = id + 1
                  write(fid,"(A4,I2,A15)") "  > ",id,". SGS transport";                     id = id + 1
                  write(fid,"(A4,I2,A23)") "  > ",id,". Molecular destruction";             id = id + 1
                  write(fid,"(A4,I2,A17)") "  > ",id,". SGS destruction";                   id = id + 1
                  write(fid,"(A4,I2,A30)") "  > ",id,". Interscale transfer (direct)";      id = id + 1
                  write(fid,"(A4,I2,A32)") "  > ",id,". Interscale transfer (indirect)";    id = id + 1
                  write(fid,"(A4,I2,A35)") "  > ",id,". Turbulent transport (mixed scale)"; id = id + 1
                  write(fid,"(A4,I2,A18)") "  > ",id,". Force production";                  id = id + 1
                  write(fid,*) " "
                  write(fid,*) "Misc terms:"
                  write(fid,"(A4,I2,A28)") "  > ",id,". Mean scalar gradient, dTdz"; id = id + 1
                  write(fid,"(A4,I2,A14)") "  > ",id,". Mean scalar ";               id = id + 1 
                  write(fid,"(A4,I2,A20)") "  > ",id,". Mean scalar source";         id = id + 1 
                  write(fid,"(A4,I2,A10)") "  > ",id,". <fT'fT'>";                   id = id + 1 
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
                  write(fid,"(A6,I2,A30)") "    > ",id,". Interscale transfer (direct)";      id = id + 1
                  write(fid,"(A6,I2,A32)") "    > ",id,". Interscale transfer (indirect)";    id = id + 1
                  write(fid,"(A6,I2,A35)") "    > ",id,". Turbulent transport (mixed scale)"; id = id + 1
                  write(fid,"(A6,I2,A25)") "    > ",id,". Scalar force production";           id = id + 1
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
        ConvTrans_fg,TurbTrans_fg)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in) :: f, g, dfdx, dfdy, dfdz, dgdx, dgdy, dgdz
        real(rkind), dimension(:), intent(in) :: ConvTrans_fg
        real(rkind), dimension(:), intent(out) :: TurbTrans_fg
        ! Inputs:
        !   > f, g, ddxj(f,g) --> These are turbulent quantities (zero-mean)
        ! In the calculation below we use the full u,v,w fields which could potentially have mean components, hence we need to
        ! subtract out the convective transport

        this%sim%rbuffxC(:,:,:,1) = this%sim%u*dgdx + this%sim%v*dgdy + this%sim%wC*dgdz
        this%sim%rbuffxC(:,:,:,2) = this%sim%u*dfdx + this%sim%v*dfdy + this%sim%wC*dfdz
        call this%covariance(f,this%sim%rbuffxC(:,:,:,1),this%zbuff(:,1))
        call this%covariance(g,this%sim%rbuffxC(:,:,:,2),this%zbuff(:,2))
        TurbTrans_fg = -this%zbuff(:,1) - this%zbuff(:,2) + ConvTrans_fg
    end subroutine

    subroutine gradient_production(this,dgdx,dgdy,dgdz,fu,fv,fw,production,desc)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in) :: dgdx, dgdy, dgdz
        real(rkind), dimension(:), intent(in) :: fu, fv, fw
        real(rkind), dimension(:), intent(out) :: production
        character(len=*), intent(in), optional :: desc

        call this%mean(dgdx,this%zbuff(:,1))
        call this%mean(dgdy,this%zbuff(:,2))
        call this%mean(dgdz,this%zbuff(:,3))

        if (present(desc)) then
            print*, "max_abs of mean(dgdx) for "//trim(desc)//" = ",maxval(abs(this%zbuff(:,1)))
            print*, "max_abs of mean(dgdy) for "//trim(desc)//" = ",maxval(abs(this%zbuff(:,2)))
            print*, "max_abs of mean(dgdz) for "//trim(desc)//" = ",maxval(abs(this%zbuff(:,3)))
            print*, "max_abs of fu         for "//trim(desc)//" = ",maxval(abs(fu))
            print*, "max_abs of fv         for "//trim(desc)//" = ",maxval(abs(fv))
            print*, "max_abs of fw         for "//trim(desc)//" = ",maxval(abs(fw))

            print*, "mean(dgdx) for "//trim(desc)//" = ",this%zbuff(:,1)
            print*, "mean(dgdy) for "//trim(desc)//" = ",this%zbuff(:,2)
            print*, "mean(dgdz) for "//trim(desc)//" = ",this%zbuff(:,3)
            print*, "fu         for "//trim(desc)//" = ",fu
            print*, "fv         for "//trim(desc)//" = ",fv
            print*, "fw         for "//trim(desc)//" = ",fw
        end if

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

    subroutine compute_unsteady_terms(this,terms,sca,scl)
        class(stats_xy), intent(inout) :: this
        character(len=*), intent(in) :: terms
        integer, intent(in) :: sca, scl

        ! Step 1: store TKE and temperature variance required for unsteady term at t^n-1
        if (mod(this%sim%step+1,this%compute_freq) == 0) then
            if (terms == 'velocity') then
                call this%covariance(this%u,  this%u,  this%uvarold(:,scl,1))
                call this%covariance(this%v,  this%v,  this%uvarold(:,scl,2))
                call this%covariance(this%wC, this%wC, this%uvarold(:,scl,3))
            else if (terms == 'scalar') then
                call this%covariance(this%T,  this%T,  this%svarold(:,sca,scl,1))
                call this%covariance(this%wC, this%T,  this%svarold(:,sca,scl,2))
            end if
            
        ! Step 2: Compute ddt terms (using t^n+1 if dt=const or t^n
        ! otherwise)
        elseif (mod(this%sim%step,this%compute_freq) == 0) then
            if (this%sim%CFL >= 0.d0) then
                if (terms == 'velocity') then
                    call this%covariance(this%u,  this%u,  this%uvarnew(:,scl,1))
                    call this%covariance(this%v,  this%v,  this%uvarnew(:,scl,2))
                    call this%covariance(this%wC, this%wC, this%uvarnew(:,scl,3))

                    this%R11_budget(:,1) = (this%uvarnew(:,scl,1) - this%uvarold(:,scl,1))/(this%sim%dt)
                    this%R22_budget(:,1) = (this%uvarnew(:,scl,2) - this%uvarold(:,scl,2))/(this%sim%dt)
                    this%R33_budget(:,1) = (this%uvarnew(:,scl,3) - this%uvarold(:,scl,3))/(this%sim%dt)
                    this%tke_budget(:,1) = 0.5d0*(this%R11_budget(:,1) + this%R22_budget(:,1) + this%R33_budget(:,1))
                else if (terms == 'scalar') then
                    call this%covariance(this%T,  this%T,  this%svarnew(:,sca,scl,1))
                    call this%covariance(this%wC, this%T,  this%svarnew(:,sca,scl,2))
                    this%TT_budget(:,1) = (this%svarnew(:,sca,scl,1) - this%svarold(:,sca,scl,1))/(this%sim%dt)
                    this%wT_budget(:,1) = (this%svarnew(:,sca,scl,2) - this%svarold(:,sca,scl,2))/(this%sim%dt)
                end if
            end if
        elseif( (mod(this%sim%step-1,this%compute_freq) == 0) .and. &
                (this%nwrite > 0)                           ) then
            if (this%sim%CFL < 0.d0) then
                if (terms == 'velocity') then
                    call this%covariance(this%u,  this%u,  this%uvarnew(:,scl,1))
                    call this%covariance(this%v,  this%v,  this%uvarnew(:,scl,2))
                    call this%covariance(this%wC, this%wC, this%uvarnew(:,scl,3))

                    this%R11_budget(:,1) = (this%uvarnew(:,scl,1) - this%uvarold(:,scl,1))/(2.d0*this%sim%dt)
                    this%R22_budget(:,1) = (this%uvarnew(:,scl,2) - this%uvarold(:,scl,2))/(2.d0*this%sim%dt)
                    this%R33_budget(:,1) = (this%uvarnew(:,scl,3) - this%uvarold(:,scl,3))/(2.d0*this%sim%dt)
                    this%tke_budget(:,1) = 0.5d0*(this%R11_budget(:,1) + this%R22_budget(:,1) + this%R33_budget(:,1))
                else if (terms == 'scalar') then
                    call this%covariance(this%T,  this%T,  this%svarnew(:,sca,scl,1))
                    call this%covariance(this%wC, this%T,  this%svarnew(:,sca,scl,2))
                    this%TT_budget(:,1) = (this%svarnew(:,sca,scl,1) - this%svarold(:,sca,scl,1))/(2.d0*this%sim%dt)
                    this%wT_budget(:,1) = (this%svarnew(:,sca,scl,2) - this%svarold(:,sca,scl,2))/(2.d0*this%sim%dt)
                end if
            end if
        end if
    end subroutine

    subroutine copy_stats(this,stats,stats_sca)
        class(stats_xy), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(inout) :: stats
        real(rkind), dimension(:,:,:,:), intent(inout), optional :: stats_sca

        call assert(size(stats,1) == size(this%stats,1),'size(stats,1) == size(this%stats,1)')
        call assert(size(stats,2) == size(this%stats,2),'size(stats,2) == size(this%stats,2)')
        call assert(size(stats,3) == size(this%stats,4),'size(stats,3) == size(this%stats,4)')

        stats = this%stats(:,:,this%tidx,:)

        if (present(stats_sca)) then
            call assert(size(stats_sca,1) == size(this%stats_sca,1),'size(stats_sca,1) == size(this%stats_sca,1)')
            call assert(size(stats_sca,2) == size(this%stats_sca,2),'size(stats_sca,2) == size(this%stats_sca,2)')
            call assert(size(stats_sca,3) == size(this%stats_sca,4),'size(stats_sca,3) == size(this%stats_sca,4)')
            call assert(size(stats_sca,4) == size(this%stats_sca,5),'size(stats_sca,4) == size(this%stats_sca,5)')

            stats_sca = this%stats_sca(:,:,this%tidx,:,:)
        end if

    end subroutine

    subroutine set_tidx(this,tid_new)
        class(stats_xy), intent(inout) :: this
        integer, intent(in) :: tid_new
        this%tidx = tid_new
    end subroutine
    subroutine set_nwrite(this,nwrite_new)
        class(stats_xy), intent(inout) :: this
        integer, intent(in) :: nwrite_new
        this%nwrite = nwrite_new
    end subroutine
    
    subroutine print_var(this,xst,xen,yst,yen,zst,zen,scl,var)
        class(stats_xy), intent(inout) :: this
        integer, intent(in) :: xst, xen, yst, yen, zst, zen, scl, var

        if (var < 4) then
            print*, this%ujsplit(xst:xen,yst:yen,zst:zen,scl,var)
        else if (var >= 4 .and. var < 11) then
            print*, this%tauij(xst:xen,yst:yen,zst:zen,scl,var-3)
        end if
    end subroutine

end module stats_xy_mod
