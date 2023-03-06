module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval, p_sum
    use exits,              only: message, message_min_max
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    use decomp_2d
    use HIT_Periodic_parameters, only: uTarget, vTarget, wTarget, x_shift, &
      uadvect, useBandpassFilter, confirmFFTresult, checkEnergyInjectionRate, &
      checkEnergyInjectionRateFreq, dumpFieldsForDebugging
    use fortran_assert,     only: assert

    implicit none 

    integer :: nt_print2screen = 1
    real(rkind) :: DomMaxDiv
    integer :: ierr 

contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout), target :: gp 
        real(rkind), dimension(:,:,:), pointer ::x, y, z
        real(rkind) :: Einjection
        complex(rkind), dimension(gp%sp_gpC%ysz(1), gp%sp_gpC%ysz(2), &
          gp%sp_gpC%ysz(3)) :: uhat_xy, vhat_xy, whatC_xy
        complex(rkind), dimension(gp%sp_gpE%ysz(1), gp%sp_gpE%ysz(2), &
          gp%sp_gpE%ysz(3)) :: what_xy
        complex(rkind), dimension(gp%sp_gpC%zsz(1), gp%sp_gpC%zsz(2), &
          gp%sp_gpC%zsz(3)) :: uhat, vhat, whatC_z
        complex(rkind), dimension(gp%sp_gpE%zsz(1), gp%sp_gpE%zsz(2), &
          gp%sp_gpE%zsz(3)) :: what_z
        !complex(rkind), dimension(:,:,:), pointer :: uhat_xy, vhat_xy, &
        !  whatC_xy, what_xy!, uhat, vhat, whatC_z, what_z
        real(rkind) :: uXoddMaxI, uYoddMaxI, &
                       uXoddMaxR, uYoddMaxR, &
                       vXoddMaxI, vYoddMaxI, &
                       vXoddMaxR, vYoddMaxR, &
                       wXoddMaxI, wYoddMaxI, &
                       wXoddMaxR, wYoddMaxR

        !uhat_xy => gp%cbuffyC(:,:,:,1)
        !vhat_xy => gp%cbuffyC(:,:,:,2)
        !whatC_xy => gp%cbuffyC(:,:,:,3)
        !what_xy => gp%cbuffyE(:,:,:,1)
         
         z => gp%mesh(:,:,:,3)
         y => gp%mesh(:,:,:,2)
         x => gp%mesh(:,:,:,1)
   
         x_shift = x_shift + uadvect*gp%get_dt()

         if (mod(gp%step,gp%t_dataDump) == 0) then
            if (useBandpassFilter) then
               call message(0,"Shifting solution by", x_shift)
               call gp%spectC%bandpassFilter_and_PhaseShift(gp%uhat , uTarget, x_shift)
               call gp%spectC%bandpassFilter_and_PhaseShift(gp%vhat , vTarget, x_shift)
               call gp%spectC%bandpassFilter_and_PhaseShift(gp%whatC, wTarget, x_shift)
               call toc()
               call gp%dumpFullField(uTarget,'uBPF')
               call gp%dumpFullField(vTarget,'vBPF')
               call gp%dumpFullField(wTarget,'wBPF')
               call message(0, "Just dumped bandpass filtered fields")
            end if 
         end if 
        
         if (mod(gp%step,nt_print2screen) == 0) then
            DomMaxDiv = p_maxval(gp%divergence)
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX",gp%step)
            call message(1,"MaxDiv",DomMaxDiv)
            call message(1,"Mean TKE",gp%getMeanKE())
            call message(1,"Mean  uu",gp%getMeanuu())
            call message(1,"Mean  uv",gp%getMeanuv())
            call message(1,"Mean  uw",gp%getMeanuw())
            call message(1,"Mean  vv",gp%getMeanvv())
            call message(1,"Mean  vw",gp%getMeanvw())
            call message(1,"Mean  ww",gp%getMeanww())

            call message(1,"Mean  u", p_sum(sum(gp%u) )/real(gp%nx*gp%ny*gp%nz,rkind))
            call message(1,"Mean  v", p_sum(sum(gp%v) )/real(gp%nx*gp%ny*gp%nz,rkind))
            call message(1,"Mean  w", p_sum(sum(gp%wC))/real(gp%nx*gp%ny*gp%nz,rkind))
            
            call message_min_max(1,"Bounds for u", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
           
            if (confirmFFTresult)  then
              call gp%spectC%fft(gp%u, uhat_xy)
              call gp%spectC%fft(gp%v, vhat_xy)
              call gp%spectC%fft(gp%wC, whatC_xy)
              call gp%spectE%fft(gp%w, what_xy)

              call message(1,"max real(uhat) difference",p_maxval(maxval(real(uhat_xy - gp%uhat,rkind))))
              call message(1,"max dimag(uhat) difference",p_maxval(maxval(dimag(uhat_xy - gp%uhat))))
              call message(1,"max real(vhat) difference",p_maxval(maxval(real(vhat_xy - gp%vhat,rkind))))
              call message(1,"max dimag(vhat) difference",p_maxval(maxval(dimag(vhat_xy - gp%vhat))))
              call message(1,"max real(what) difference",p_maxval(maxval(real(what_xy - gp%what,rkind))))
              call message(1,"max dimag(what) difference",p_maxval(maxval(dimag(what_xy - gp%what))))
              call message(1,"max real(whatC) difference",p_maxval(maxval(real(whatC_xy - gp%whatC,rkind))))
              call message(1,"max dimag(whatC) difference",p_maxval(maxval(dimag(whatC_xy - gp%whatC))))
            end if

            if (checkEnergyInjectionRate .and. allocated(gp%hitforce)) then
              if (mod(gp%step,checkEnergyInjectionRateFreq) == 0) then
                call transpose_y_to_z(gp%uhat,uhat,gp%sp_gpC)
                call gp%spectC%take_fft1d_z2z_ip(uhat)

                call transpose_y_to_z(gp%vhat,vhat,gp%sp_gpC)
                call gp%spectC%take_fft1d_z2z_ip(vhat)

                call transpose_y_to_z(gp%whatC,whatC_z,gp%sp_gpC)
                call gp%spectC%take_fft1d_z2z_ip(whatC_z)
                
                call gp%hitforce%computeEnergyInjectionRate(uhat, vhat, whatC_z, &
                  gp%hitforce%fxhatwrite, gp%hitforce%fyhatwrite, gp%hitforce%fzhatwrite,&
                  Einjection)
                call message(1,"HIT energy injection rate from spectral velocity:", Einjection)
                if (abs(gp%hitforce%EpsAmplitude - Einjection) > 1.e-11) then
                  call message(0,'Energy injection rate does not match input file. Difference:',&
                    abs(gp%hitforce%EpsAmplitude - Einjection))
                end if

                call transpose_y_to_z(uhat_xy, uhat, gp%sp_gpC)
                call gp%spectC%take_fft1d_z2z_ip(uhat)

                call transpose_y_to_z(vhat_xy, vhat, gp%sp_gpC)
                call gp%spectC%take_fft1d_z2z_ip(vhat)

                call transpose_y_to_z(whatC_xy, whatC_z, gp%sp_gpC)
                call gp%Pade6opZ%interpz_C2E(whatC_z, what_z, 0, 0)
                !call transpose_z_to_y(what_z,what_xy,gp%sp_gpE)
                whatC_z = what_z(:,:,1:gp%sp_gpC%zsz(3))
                call gp%spectC%take_fft1d_z2z_ip(whatC_z)
                call gp%spectC%shiftz_E2C(whatC_z)
                
                call gp%hitforce%computeEnergyInjectionRate(uhat, vhat, whatC_z, &
                  gp%hitforce%fxhatwrite, gp%hitforce%fyhatwrite, gp%hitforce%fzhatwrite,&
                  Einjection)
                call message(1,"HIT energy injection rate from physical velocity:", Einjection)

                if (abs(gp%hitforce%EpsAmplitude - Einjection) > 1.e-11) then
                  call message(0,'Energy injection rate does not match input file. Difference:',&
                    abs(gp%hitforce%EpsAmplitude - Einjection))
                  if (dumpFieldsForDebugging) then
                    call gp%dumpFullField(gp%u,'uphC')
                    call gp%dumpFullField(gp%v,'vphC')
                    call gp%dumpFullField(gp%wC,'wphC')

                    call gp%dumpFullField(gp%uE,'uphE',gp%gpE)
                    call gp%dumpFullField(gp%vE,'vphE',gp%gpE)
                    call gp%dumpFullField(gp%w,'wphE',gp%gpE)

                    call gp%dumpSpectralField(dimag(gp%uhat),'uI_C',gp%sp_gpC)
                    call gp%dumpSpectralField(dimag(gp%vhat),'vI_C',gp%sp_gpC)
                    call gp%dumpSpectralField(dimag(gp%whatC),'wI_C',gp%sp_gpC)
                    call gp%dumpSpectralField(real(gp%uhat,rkind),'uR_C',gp%sp_gpC)
                    call gp%dumpSpectralField(real(gp%vhat,rkind),'vR_C',gp%sp_gpC)
                    call gp%dumpSpectralField(real(gp%whatC,rkind),'wR_C',gp%sp_gpC)
                    
                    call gp%dumpSpectralField(dimag(gp%uEhat),'uI_E',gp%sp_gpE)
                    call gp%dumpSpectralField(dimag(gp%vEhat),'vI_E',gp%sp_gpE)
                    call gp%dumpSpectralField(dimag(gp%what),'wI_E',gp%sp_gpE)
                    call gp%dumpSpectralField(real(gp%uEhat,rkind),'uR_E',gp%sp_gpE)
                    call gp%dumpSpectralField(real(gp%vEhat,rkind),'vR_E',gp%sp_gpE)
                    call gp%dumpSpectralField(real(gp%what,rkind),'wR_E',gp%sp_gpE)

                    call gp%hitforce%dumpForcing(gp%outputdir,gp%RunID,gp%step,dumpSpec=.true.)
                  
                  end if
                  call assert(abs(gp%hitforce%EpsAmplitude - Einjection) < 1.e-11)
                end if
              end if
            end if
            if (allocated(gp%pressure)) then
               call message_min_max(1,"Bounds for P:", p_minval(minval(gp%pressure)), p_maxval(maxval(gp%pressure)))
            end if
            if (gp%useSGS) then
                call message_min_max(1,"Bounds for nuSGS:", p_minval(minval(gp%nu_sgs)), p_maxval(maxval(gp%nu_sgs)))
                if (gp%sgsModel%usingDynProc()) then
                    call message(2,"Model multiplier, (c*delta)^2:",gp%sgsModel%getMax_DynSmagConst())
                end if 
            end if
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            call toc()
            call tic()
        end if 

    end subroutine

end module 
