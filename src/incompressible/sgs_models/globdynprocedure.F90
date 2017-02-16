    ! -- in dynamicprocedure.F90 ---
    !subroutine testFilter_ip(this, field)
    !subroutine testFilter_oop(this, fin, fout)
    !subroutine testFilter_oop_C2R(this, fhat, fout)

    subroutine GlobDynProcedure(this, u, v, wC, uhat, vhat, wChat, S13, S23, S13C, S23C, duidxjC, duidxjE, duidxjChat, totBodyForce)
        class(sgs),                                                                         intent(inout), target :: this
        real(rkind),    dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),         intent(in)            :: u, v, wC
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)),   intent(in)            :: uhat, vhat, wChat
        real(rkind),    dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)),         intent(in)            :: S13, S23
        real(rkind),    dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),         intent(in)            :: S13C, S23C
        real(rkind),    dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9),       intent(in),    target :: duidxjC
        real(rkind),    dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),4),       intent(in),    target :: duidxjE
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3),9), intent(in)            :: duidxjChat
        real(rkind),    dimension(:,:,:,:),                                                               pointer :: totBodyForce


        real(rkind),    dimension(:,:,:), pointer :: dudx,  dudy,  dudzC, dvdx,  dvdy,  dvdzC, dwdxC, dwdyC, dwdz
        complex(rkind), dimension(:,:,:), pointer :: ctmpY, ctmpZ 
        real(rkind),    dimension(:,:,:), pointer :: uf, vf, wcf, numerator, denominator, numerTF, denomTF, tmpbuff, S13f, S23f
        real(rkind),    dimension(:,:),   pointer :: tauWM13, tauWM23, S13CTF, S23CTF
        real(rkind) :: sum_numerator, sum_denominator
        integer :: idx, interpOption = 2

        if(nrank==0) write(*,*) 'In Global Dynamic Procedure'
        ! denominator first
        denomTF     => this%GlobDynrbuffxC(:,:,:,13)
        denominator => this%GlobDynrbuffxC(:,:,:,14)
        tmpbuff     => this%GlobDynrbuffxC(:,:,:,15)
        !------------Compute denominator :: f(testFilter-quantities) part ---------------
        do idx = 1, 9
          call this%testFilter_oop_C2R(duidxjChat(:,:,:,idx),this%GlobDynrbuffxC(:,:,:,3+idx))
        enddo
        dudx  => this%GlobDynrbuffxC(:,:,:,4);  dudy  => this%GlobDynrbuffxC(:,:,:,5);  dudzC => this%GlobDynrbuffxC(:,:,:,6); 
        dvdx  => this%GlobDynrbuffxC(:,:,:,7);  dvdy  => this%GlobDynrbuffxC(:,:,:,8);  dvdzC => this%GlobDynrbuffxC(:,:,:,9); 
        dwdxC => this%GlobDynrbuffxC(:,:,:,10); dwdyC => this%GlobDynrbuffxC(:,:,:,11); dwdz  => this%GlobDynrbuffxC(:,:,:,12); 
        call this%get_SIGMA_Op(denomTF, dudx, dudy, dudzC, dvdx, dvdy, dvdzC, dwdxC, dwdyC, dwdz)
        !write(*,*) 'nuSGSTF: ', maxval(denomTF), minval(denomTF)
        if(nrank==0) write(*,*) 'here 1'

        ! multiply by SijTF
        denominator = dudx*dudx + dvdy*dvdy + dwdz*dwdz
        tmpbuff = half*(dudy  + dvdx);  denominator = denominator + two*tmpbuff*tmpbuff
        tmpbuff = half*(dudzC + dwdxC); denominator = denominator + two*tmpbuff*tmpbuff
        tmpbuff = half*(dvdzC + dwdyC); denominator = denominator + two*tmpbuff*tmpbuff
        denomTF = denomTF * denominator
        !write(*,*) 'denomTF: ', maxval(denomTF), minval(denomTF)
        !------------Done denominator :: f(testFilter-quantities) part ---------------
        if(nrank==0) write(*,*) 'here 2'

        if(interpOption==1) then
          ! store S13CTF(:,:,1) and S23CTF(:,:,1) in rbuffxC(:,:,1:2,13) for later use
          S13CTF => this%GlobDynrbuffxC(:,:,1,13);         S23CTF => this%GlobDynrbuffxC(:,:,2,13)
          S13CTF = half*(dudzC(:,:,1) + dwdxC(:,:,1));     S23CTF = half*(dvdzC(:,:,1) + dwdyC(:,:,1))
        endif
        if(nrank==0) write(*,*) 'here 3'

        !------------Compute numerator :: f(testFilter-quantities) part involving TF(duidxj)-------
        numerTF => this%GlobDynrbuffxC(:,:,:,4)
        if(this%isInviscid) then
          numerTF = zero
        else
          numerTF = this%GlobDynrbuffxC(:,:,:,4) * this%GlobDynrbuffxC(:,:,:,4)
          do idx = 2, 9
            numerTF = numerTF + this%GlobDynrbuffxC(:,:,:,3+idx) * this%GlobDynrbuffxC(:,:,:,3+idx)
          enddo
          numerTF = this%invRe*numerTF
        endif
        !------------Done numerator :: f(testFilter-quantities) part involving TF(duidxj)-------
        !---- Now TF(duidxj) has been used and this%GlobDynrbuffxC(:,:,:,5:12) may be overwritten----- 
        !---- 1:3 :: totBodyForce; 4 :: numerTF; 13 :: S13CTF, S23CTF; 14 :: denominator ------------- 
        !---- should be preserved, all else may be overwritten----------------------------------------
        if(nrank==0) write(*,*) 'here 4'

        ! Allocate pointers
        dudx  => duidxjC(:,:,:,1); dudy  => duidxjC(:,:,:,2); dudzC => duidxjC(:,:,:,3); 
        dvdx  => duidxjC(:,:,:,4); dvdy  => duidxjC(:,:,:,5); dvdzC => duidxjC(:,:,:,6); 
        dwdxC => duidxjC(:,:,:,7); dwdyC => duidxjC(:,:,:,8); dwdz  => duidxjC(:,:,:,9); 

        ctmpZ => this%ctmpCz; ctmpY => this%cbuff(:,:,:,1)

        uf  => this%GlobDynrbuffxC(:,:,:,5); vf  => this%GlobDynrbuffxC(:,:,:,6); wcf => this%GlobDynrbuffxC(:,:,:,7)
        numerator => this%GlobDynrbuffxC(:,:,:,8)
        tauWM13 => this%GlobDynrbuffxC(:,:,3,13); tauWM23 => this%GlobDynrbuffxC(:,:,4,13)
        S13f => this%GlobDynrbuffxE(:,:,:,1); S23f => this%GlobDynrbuffxE(:,:,:,2)


        ! compute remaining part of denominator
        !------------Compute denominator :: testFilter of f(grid-quantities) part ---------------
        denominator = dudx*dudx + dvdy*dvdy + dwdz*dwdz
        tmpbuff = half*(dudy  + dvdx);  denominator = denominator + two*tmpbuff*tmpbuff
        tmpbuff = half*(dudzC + dwdxC); denominator = denominator + two*tmpbuff*tmpbuff
        tmpbuff = half*(dvdzC + dwdyC); denominator = denominator + two*tmpbuff*tmpbuff
        !write(*,*) 'dudx   : ', maxval(dudx), minval(dudx)
        !write(*,*) 'dudy   : ', maxval(dudy), minval(dudy)
        !write(*,*) 'dudzC  : ', maxval(dudzC), minval(dudzC)
        !write(*,*) 'dvdx   : ', maxval(dvdx), minval(dvdx)
        !write(*,*) 'dvdy   : ', maxval(dvdy), minval(dvdy)
        !write(*,*) 'dvdzC  : ', maxval(dvdzC), minval(dvdzC)
        !write(*,*) 'dwdxC  : ', maxval(dwdxC), minval(dwdxC)
        !write(*,*) 'dwdyC  : ', maxval(dwdyC), minval(dwdyC)
        !write(*,*) 'dwdyz  : ', maxval(dwdz), minval(dwdz)
        !write(*,*) 'denom  : ', maxval(denominator), minval(denominator)
        denominator = denominator*this%nuSGS
        !write(*,*) 'nuSGS  : ', maxval(this%nuSGS), minval(this%nuSGS)
        !write(*,*) 'nuSGSij: ', maxval(denominator), minval(denominator)
        call this%testFilter_ip(denominator)
        !write(*,*) 'TFnuSGS: ', maxval(denominator), minval(denominator)
        denominator = denominator - four*denomTF
        !write(*,*) 'denomin: ', maxval(denominator), minval(denominator)
        !------------Done denominator :: testFilter of f(grid-quantities) part ---------------

        !------------Compute numerator :: testFilter of f(grid-quantities) part ---------------
        if(this%isInviscid) then
          numerator = zero
        else
          numerator = -duidxjC(:,:,:,1) * duidxjC(:,:,:,1)
          do idx = 2,9
              numerator = numerator - duidxjC(:,:,:,idx) * duidxjC(:,:,:,idx)
          end do
          numerator = this%invRe*numerator 
        endif

        if(associated(totBodyForce)) then
          numerator = numerator + u  * totBodyForce(:,:,:,1) & 
                                + v  * totBodyForce(:,:,:,2) & 
                                + wC * totBodyForce(:,:,:,3)
        endif

        if(nrank==0) write(*,*) 'here 5'
        select case(this%WallModel)
        case(0)
          ! Compute wall shear stress - Moeng's model
        case(1)
          ! Compute wall shear stress - Bou-Zeid (2005) model
          call transpose_y_to_z(uhat, ctmpZ, this%sp_gp)
          call this%spectC%SurfaceFilter_ip(ctmpZ(:,:,1))
          call transpose_z_to_y(ctmpZ, ctmpY, this%sp_gp)
          call this%spectC%ifft(ctmpY, uf)
          !write(*,*) 'u1     : ', maxval(u(:,:,1)), minval(u(:,:,1))
          !write(*,*) 'u1f    : ', maxval(uf(:,:,1)), minval(uf(:,:,1))

          call transpose_y_to_z(vhat, ctmpZ, this%sp_gp)
          call this%spectC%SurfaceFilter_ip(ctmpZ(:,:,1))
          call transpose_z_to_y(ctmpZ, ctmpY, this%sp_gp)
          call this%spectC%ifft(ctmpY, vf)
          !write(*,*) 'u2     : ', maxval(v(:,:,1)), minval(v(:,:,1))
          !write(*,*) 'u2f    : ', maxval(vf(:,:,1)), minval(vf(:,:,1))

          if(this%zEdgeLo==1) then
            ! wall model is active here
            tauWM13(:,:) = this%WallMFactor*sqrt(uf(:,:,1)*uf(:,:,1) + vf(:,:,1)*vf(:,:,1))
            tauWM23(:,:) = tauWM13 * vf(:,:,1)
            tauWM13(:,:) = tauWM13 * uf(:,:,1)
            !write(*,*) 'numerat: ', maxval(numerator), minval(numerator)
            !write(*,*) 'taWM13 : ', maxval(tauWM13), minval(tauWM13)
            !write(*,*) 'taWM23 : ', maxval(tauWM23), minval(tauWM23)
            !write(*,*) 'S13C   : ', maxval(S13C  ), minval(S13C  )
            !write(*,*) 'S23C   : ', maxval(S23C  ), minval(S23C  )
            !write(*,*) '------------------gridfilter-----------------'
            if(interpOption==1) then
              numerator(:,:,1) = numerator(:,:,1) + tauWM13*S13C(:,:,1) + tauWM23*S23C(:,:,1)
              !write(*,*) 'S13    : ', maxval(S13C(:,:,1)), minval(S13C(:,:,1))
              !write(*,*) 'S23    : ', maxval(S23C(:,:,1)), minval(S23C(:,:,1))
            else
              numerator(:,:,1) = numerator(:,:,1) + tauWM13*S13(:,:,1) + tauWM23*S23(:,:,1)
              !write(*,*) 'S13    : ', maxval(S13(:,:,1)), minval(S13(:,:,1))
              !write(*,*) 'S23    : ', maxval(S23(:,:,1)), minval(S23(:,:,1))
            endif
            !write(*,*) 'taWM13 : ', maxval(tauWM13), minval(tauWM13)
            !write(*,*) 'taWM23 : ', maxval(tauWM23), minval(tauWM23)
            !write(*,*) 'numerat: ', maxval(numerator), minval(numerator), p_sum(sum(numerator))
          endif 
        end select
        if(nrank==0) write(*,*) 'here 6'

        call this%testFilter_ip(numerator)
        if(nrank==0) write(*,*) 'here 7'
        !write(*,*) 'numefil: ', maxval(numerator), minval(numerator), p_sum(sum(numerator))
        !------------Done numerator :: testFilter of f(grid-quantities) part ---------------

        !------------Compute numerator :: f(testFilter-quantities) part ---------------
        if(associated(totBodyForce)) then
          do idx = 1, 3
            call this%testFilter_ip(totBodyForce(:,:,:,idx))
          enddo
          call this%testFilter_oop_C2R(uhat,uf)     !rbuufx4
          call this%testFilter_oop_C2R(vhat,vf)     !rbuffx5
          call this%testFilter_oop_C2R(wChat,wcf)   !rbuffx6
       
          numerTF = numerTF - uf  * totBodyForce(:,:,:,1) & 
                            - vf  * totBodyForce(:,:,:,2) & 
                            - wcf * totBodyForce(:,:,:,3)
        endif
        if(nrank==0) write(*,*) 'here 8'

        select case(this%WallModel)
        case(0)
        case(1)
          ! Compute wall shear stress - Bou-Zeid (2005) model
          !write(*,*) 'u1TF   : ', maxval(uf(:,:,1)), minval(uf(:,:,1))
          call this%spectC%fft(uf, ctmpY)
          call transpose_y_to_z(ctmpY, ctmpZ, this%sp_gp)
          call this%spectC%SurfaceFilter_ip(ctmpZ(:,:,1))
          call transpose_z_to_y(ctmpZ, ctmpY, this%sp_gp)
          call this%spectC%ifft(ctmpY, uf)  ! uf now contains twice-filtered velocity
          !write(*,*) 'u1fTF  : ', maxval(uf(:,:,1)), minval(uf(:,:,1))
           
          !write(*,*) 'u2TF   : ', maxval(vf(:,:,1)), minval(vf(:,:,1))
          call this%spectC%fft(vf, ctmpY)
          call transpose_y_to_z(ctmpY, ctmpZ, this%sp_gp)
          call this%spectC%SurfaceFilter_ip(ctmpZ(:,:,1))
          call transpose_z_to_y(ctmpZ, ctmpY, this%sp_gp)
          call this%spectC%ifft(ctmpY, vf)  ! vf now contains twice-filtered velocity
          !write(*,*) 'u2fTF  : ', maxval(vf(:,:,1)), minval(vf(:,:,1))
          if(nrank==0) write(*,*) 'here 9'

          !write(*,*) '------------------testfilter-----------------'
          if(interpOption==1) then
            !write(*,*) 'S13    : ', maxval(S13CTF(:,:)), minval(S13CTF(:,:))
            !write(*,*) 'S23    : ', maxval(S23CTF(:,:)), minval(S23CTF(:,:))
          else
            ! test-filter S13(:,:,1)
            call this%spectE%fft(S13, this%ctmpEy)
            call transpose_y_to_z(this%ctmpEy, this%ctmpEz, this%sp_gpE)
            call this%spectE%TestFilter_ip(this%ctmpEz(:,:,1))
            call transpose_z_to_y(this%ctmpEz, this%ctmpEy, this%sp_gpE)
            call this%spectE%ifft(this%ctmpEy, S13f)
            if(nrank==0) write(*,*) 'here 10'

            ! test-filter S23(:,:,1)
            call this%spectE%fft(S23, this%ctmpEy)
            call transpose_y_to_z(this%ctmpEy, this%ctmpEz, this%sp_gpE)
            call this%spectE%TestFilter_ip(this%ctmpEz(:,:,1))
            call transpose_z_to_y(this%ctmpEz, this%ctmpEy, this%sp_gpE)
            call this%spectE%ifft(this%ctmpEy, S23f)
            if(nrank==0) write(*,*) 'here 11'

            !write(*,*) 'S13    : ', maxval(S13f(:,:,1)), minval(S13f(:,:,1))
            !write(*,*) 'S23    : ', maxval(S23f(:,:,1)), minval(S23f(:,:,1))
          endif
          if(this%zEdgeLo==1) then
            ! wall model is active here
            tauWM13(:,:) = this%WallMFactor*sqrt(uf(:,:,1)*uf(:,:,1) + vf(:,:,1)*vf(:,:,1))
            tauWM23(:,:) = tauWM13 * vf(:,:,1)
            tauWM13(:,:) = tauWM13 * uf(:,:,1)

            if(interpOption==1) then
              numerTF(:,:,1) = numerTF(:,:,1) - tauWM13*S13CTF - tauWM23*S23CTF
            else
              numerTF(:,:,1) = numerTF(:,:,1) - tauWM13*S13f(:,:,1) - tauWM23*S23f(:,:,1)
            endif 
            !write(*,*) 'taWM13 : ', maxval(tauWM13), minval(tauWM13)
            !write(*,*) 'taWM23 : ', maxval(tauWM23), minval(tauWM23)
            !write(*,*) 'numTFf : ', maxval(numerTF), minval(numerTF), p_sum(sum(numerTF))
          endif
        end select
        numerator = numerator + numerTF
        if(nrank==0) write(*,*) 'here 12'
        !write(*,*) 'numerat: ', maxval(numerator), minval(numerator), p_sum(sum(numerator))
        !------------Done numerator :: f(testFilter-quantities) part ---------------

        ! compute volume averages and ensure no backscatter
        sum_numerator = p_sum(sum(numerator));   sum_denominator = p_sum(sum(denominator))
        if(abs(sum_denominator) > 1.0d-9) this%GlobDynCoeff = sum_numerator/sum_denominator
        if(this%GlobDynCoeff < zero) then
          this%GlobDynCoeff = zero
        else
          this%GlobDynCoeff = this%GlobDynCoeff/(two*this%deltaFilter**2)
        endif
        if(nrank==0) then
          open(77, file='gdcoeff.dat', status='unknown', action='write', position='append')
          write(77,'(3(e19.12,1x))') sum_numerator, sum_denominator, this%GlobDynCoeff
          close(77)
        endif
        if(nrank==0) write(*,*) 'here 13'
    end subroutine
