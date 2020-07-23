module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max, GracefulExit
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi

    implicit none 

    integer :: nt_print2screen = 1
    real(rkind) ::  maxDiv, DomMaxDiv
    integer :: ierr 

contains

    subroutine doTemporalStuff(igp)
        class(igrid), intent(inout) :: igp 
        real(rkind) :: maxnusgs, maxkappasgs
      
        if (mod(igp%step,nt_print2screen) == 0) then
            maxDiv = maxval(igp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",igp%tsim)
            call message(1,"TIDX:",igp%step)
            !call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(igp%u)), p_maxval(maxval(igp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(igp%v)), p_maxval(maxval(igp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(igp%w)), p_maxval(maxval(igp%w)))
            call message_min_max(1,"Bounds for T:", p_minval(minval(igp%T)), p_maxval(maxval(igp%T)))
            call message(1, "Time step limit:" // trim(igp%dtlimit))
            call message(1,"T_surf:",igp%Tsurf)
            call message(1,"u_star:",igp%sgsmodel%get_ustar())
            call message(1,"Inv. Ob. Length:",igp%sgsmodel%get_InvObLength())
            call message(1,"wTh_surf:",igp%sgsmodel%get_wTh_surf())
            call message(1,"G_geostrophic:",igp%G_geostrophic)
            maxnusgs = p_maxval(igp%nu_SGS)
            maxkappasgs = p_maxval(igp%kappaSGS)
            call message(1,"Maximum SGS viscosity:", maxnusgs)
            call message(1,"Maximum SGS scalar kappa:", maxkappasgs)
            if (igp%useCFL) then
                call message(1,"Current dt:",igp%dt)
            end if
            call message("==========================================================")
            call toc()
            call tic()
        end if 

    end subroutine

    !subroutine set_bottomwall_heatflux(igp)
    !    class(igrid), intent(inout) :: igp

    !    real(rkind) :: tloc, Lscale, Uscale, Tscale, Tempscale

    !    Lscale = 2000.0_rkind
    !    Uscale = sqrt(3.0_rkind**2+9.0_rkind**2)
    !    Tscale = Lscale/Uscale
    !    Tempscale =1.0_rkind

    !    tloc = igp%tsim*Tscale/3600._rkind + 16.0_rkind
    !    if(tloc .le. 16.5_rkind) then
    !        igp%wTh_surf = -0.0338d0 + 0.2d0*cos(0.44d0*tloc + 0.5d0)
    !    elseif(tloc .le. 30.5_rkind) then
    !        igp%wTh_surf = -0.015d0
    !    elseif(tloc .le. 40.5_rkind) then
    !        igp%wTh_surf = -0.2094d0 - 0.395d0*cos(0.2098*tloc + 1.97d0)
    !    elseif(tloc .le. 54.5_rkind) then
    !        igp%wTh_surf = -0.0102d0
    !    elseif(tloc .le. 65.0_rkind) then
    !        igp%wTh_surf = -0.471d0 - 0.656d0*cos(0.15107d0*tloc + 0.4d0)
    !    else
    !        igp%wTh_surf = -0.002286d0*tloc + 0.13712d0   !0.1346d0
    !    endif
    !    igp%wTh_surf = igp%wTh_surf/(Uscale*Tempscale)
    !end subroutine

    !subroutine set_bottomwall_temperature(igp)
    !    class(igrid), intent(inout) :: igp

    !    real(rkind) :: tloc, Lscale, Uscale, Tscale, Tempscale

    !    Lscale = 2000.0_rkind
    !    Uscale = sqrt(3.0_rkind**2+9.0_rkind**2)
    !    Tscale = Lscale/Uscale
    !    Tempscale =1.0_rkind

    !    tloc = igp%tsim*Tscale/3600._rkind + 16.0_rkind
    !    if(tloc .le. 17.4_rkind) then
    !        igp%Tsurf = -10.0d0 - 25.0d0*cos(0.22d0*tloc + 0.2d0) + 273.15d0
    !    elseif(tloc .le. 30.0_rkind) then
    !        igp%Tsurf = -0.54d0*tloc + 15.2d0 + 273.15d0
    !    elseif(tloc .le. 41.9_rkind) then
    !        igp%Tsurf = -7.0d0 - 25.0d0*cos(0.21*tloc + 1.8d0) + 273.15d0
    !    elseif(tloc .le. 53.3_rkind) then
    !        igp%Tsurf = -0.37d0*tloc + 18.0d0 + 273.15d0
    !    elseif(tloc .le. 65.6_rkind) then
    !        igp%Tsurf = -4.0d0 - 25.0d0*cos(0.22d0*tloc + 2.5d0) + 273.15d0
    !    else
    !        igp%Tsurf = 4.4d0 + 273.15d0
    !    endif
    !    igp%Tsurf = igp%Tsurf/Tempscale
    !end subroutine

end module 
