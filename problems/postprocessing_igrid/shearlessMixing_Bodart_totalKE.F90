#include "../incompressible/shearlessMixing_Bodart_files/initialize.F90"

module KEmod
  use incompressibleGrid, only: igrid
  use kind_parameters, only: rkind
  use decomp_2d, only: transpose_y_to_z, transpose_z_to_y

  contains
    subroutine interp_E2C(sim,fE,fC)
        type(igrid), intent(inout) :: sim
        real(rkind), dimension(:,:,:), intent(in) :: fE
        real(rkind), dimension(:,:,:), intent(out) :: fC

        call transpose_y_to_z(fE,sim%rbuffzE(:,:,:,1),sim%sp_gpE)
        call sim%Pade6opZ%interpz_E2C(sim%rbuffzE(:,:,:,1),sim%rbuffzC(:,:,:,1),0,0)
        call transpose_z_to_y(sim%rbuffzC(:,:,:,1),fC,sim%sp_gpC)
    end subroutine

end module

program totalKE
    use mpi
    use kind_parameters,    only: clen, rkind
    use IncompressibleGrid, only: igrid
    use exits,              only: message
    use reductions,         only: p_sum, p_maxval
    use KEmod,              only: interp_E2C

    implicit none

    type(igrid), allocatable, target :: SM
    character(len=clen) :: inputfile
    integer :: ierr, ioUnit, tidst, tiden, tstep, tid
    real(rkind), dimension(:,:,:), allocatable :: vcon, wcon, uvisc, &
      vvisc, wvisc, usgs, vsgs, wsgs, px, py, pz, Ef, Econ, Evisc, Esgs, &
      Eprss, Esp, spx, spy, spz
    complex(rkind), dimension(:,:,:), allocatable :: ucon_hat, vcon_hat, wcon_hat,&
      uvisc_hat, vvisc_hat, wvisc_hat, usgs_hat, vsgs_hat, wsgs_hat, spx_hat, &
      spy_hat, spz_hat
    real(rkind), dimension(:,:,:), allocatable, target :: ucon
    real(rkind), dimension(:,:,:), pointer :: KE
   
    ! Initialize MPI 
    call MPI_Init(ierr)                                               
    
    ! Read input file to get parameters for post processing
    call GETARG(1,inputfile)                                            

    namelist/postProcess/ tidst, tiden, tstep    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=postProcess)
    close(ioUnit)    

    allocate(SM)                                                     
    
    call SM%init(inputfile, .true.)                              
    call SM%start_io(.true.)                                          
    call SM%printDivergence()

    ! Allocate memory
    allocate(ucon(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(vcon(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(wcon(SM%gpE%xsz(1),SM%gpE%xsz(2),SM%gpE%xsz(3)))

    allocate(uvisc(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(vvisc(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(wvisc(SM%gpE%xsz(1),SM%gpE%xsz(2),SM%gpE%xsz(3)))

    allocate(usgs(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(vsgs(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(wsgs(SM%gpE%xsz(1),SM%gpE%xsz(2),SM%gpE%xsz(3)))

    allocate(px(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(py(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(pz(SM%gpE%xsz(1),SM%gpE%xsz(2),SM%gpE%xsz(3)))

    allocate(spx(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(spy(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(spz(SM%gpE%xsz(1),SM%gpE%xsz(2),SM%gpE%xsz(3)))

    allocate(Ef(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(Econ(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(Evisc(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(Esgs(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(Eprss(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
    allocate(Esp(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))

    call SM%spectC%alloc_r2c_out(ucon_hat)
    call SM%spectC%alloc_r2c_out(vcon_hat)
    call SM%spectE%alloc_r2c_out(wcon_hat)

    call SM%spectC%alloc_r2c_out(uvisc_hat)
    call SM%spectC%alloc_r2c_out(vvisc_hat)
    call SM%spectE%alloc_r2c_out(wvisc_hat)

    call SM%spectC%alloc_r2c_out(usgs_hat)
    call SM%spectC%alloc_r2c_out(vsgs_hat)
    call SM%spectE%alloc_r2c_out(wsgs_hat)

    call SM%spectC%alloc_r2c_out(spx_hat)
    call SM%spectC%alloc_r2c_out(spy_hat)
    call SM%spectE%alloc_r2c_out(spz_hat)

    KE => ucon

! Test the time stepping
!call SM%timeAdvance()
!call SM%reinit(tidst)
!call SM%timeAdvance()
!call SM%reinit(tidst + tstep)
!call SM%timeAdvance()
!call MPI_Barrier(MPI_COMM_WORLD,ierr)
!stop
    do tid = tidst,tiden,tstep
        call message("TID", tid)
        call SM%reinit(tid)
        
        ! Get momentum RHS terms
        call SM%getConvectiveTerms(ucon_hat,vcon_hat,wcon_hat,ucon,vcon,wcon)
        call SM%getViscousTerms(uvisc_hat,vvisc_hat,wvisc_hat,uvisc,vvisc,wvisc)
        call SM%getSGSterms(usgs_hat,vsgs_hat,wsgs_hat,usgs,vsgs,wsgs)
        call SM%getPressureGradient(px,py,pz)
        call SM%getSpongeTerms(spx_hat,spy_hat,spz_hat,spx,spy,spz)

        usgs = usgs - uvisc
        vsgs = vsgs - vvisc
        wsgs = wsgs - wvisc
       
        ! Viscous
        SM%rbuffxE(:,:,:,1) = SM%w*wvisc
        call interp_E2C(SM,SM%rbuffxE(:,:,:,1),Evisc)
        Evisc = 1.d0/SM%Re*(SM%u*uvisc + SM%v*vvisc + Evisc)
        call message("Evisc ",p_sum(sum(Evisc)))

        ! Convection
        SM%rbuffxE(:,:,:,1) = SM%w*wsgs
        call interp_E2C(SM,SM%rbuffxE(:,:,:,1),Esgs)
        Esgs = Esgs + SM%u*usgs + SM%v*vsgs 
        call message("Esgs ",p_sum(sum(Esgs)))

        ! Pressure gradient work
        SM%rbuffxE(:,:,:,1) = -SM%w*pz
        call interp_E2C(SM,SM%rbuffxE(:,:,:,1),Eprss)
        Eprss = Eprss - SM%u*px - SM%v*py
        call message("Eprss ",p_sum(sum(Eprss)))

        ! Force layer work
        SM%rbuffxE(:,:,:,1) = SM%w*SM%forceLayer%fz
        call interp_E2C(SM,SM%rbuffxE(:,:,:,1),Ef)
        Ef = Ef + SM%u*SM%forceLayer%fx + SM%v*SM%forceLayer%fy
        call message("Ef ",p_sum(sum(Ef)))

        ! Convection
        SM%rbuffxE(:,:,:,1) = -SM%w*wcon
        call interp_E2C(SM,SM%rbuffxE(:,:,:,1),Econ)
        Econ = Econ - SM%u*ucon - SM%v*vcon
        call message("Econ ",p_sum(sum(Econ)))

        ! Sponge dissipation
        SM%rbuffxE(:,:,:,1) = SM%w*spz
        call interp_E2C(SM,SM%rbuffxE(:,:,:,1),Esp)
        Esp = Esp + SM%u*spx + SM%v*spy
        call message("Esp ",p_sum(sum(Esp)))

        ! Kinetic energy
        KE = 0.5d0*(SM%u*SM%u + SM%v*SM%v + SM%wC + SM%wC)
        call message("KEold", p_sum(sum(KE)))

        ! Compute E^{n+1} for ddt(E)
        !call SM%timeAdvance()
        !KE = 0.5d0*(SM%u*SM%u + SM%v*SM%v + SM%wC + SM%wC)
        !call message("KEnew", p_sum(sum(KE)))

        ! Sum over domain
        !residual(n) = p_sum(sum(Econ + Evisc + Esgs + Eprss + Ef))
        !call message("residual ",residual)

        !! Get first half of budget
        !budget(1) = budget(1) + 0.5d0*p_sum(sum(this%u*RHSx + this%v*RHSy + this%wC*RHSz))

        !! Get u_n+1
        !call SM%timeAdvance()

        !! Compute second half of budget terms
        !budget(1) = budget(1) + 0.5d0*p_sum(sum(this%u*RHSx + this%v*RHSy + this%wC*RHSz))

    end do

    deallocate(ucon, vcon, wcon, uvisc, vvisc, wvisc, usgs, vsgs, wsgs, &
      px, py, pz, Ef, Econ, Evisc, Esgs, Eprss, ucon_hat, vcon_hat, wcon_hat,&
      uvisc_hat, vvisc_hat, wvisc_hat, usgs_hat, vsgs_hat, wsgs_hat, Esp, &
      spx, spy, spz, spx_hat, spy_hat, spz_hat)
    nullify(KE)

    call SM%finalize_io()
    call SM%destroy()

    call MPI_Finalize(ierr)
end program totalKE
