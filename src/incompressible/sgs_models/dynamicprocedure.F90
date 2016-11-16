    subroutine testFilter_ip(this, field)
        class(sgs), intent(inout), target :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(inout) :: field
        complex(rkind), dimension(:,:,:), pointer :: ctmpY, ctmpZ1, ctmpZ2
        
        ctmpY => this%cbuff(:,:,:,1)
        ctmpZ1 => this%ctmpCz
        ctmpZ2 => this%ctmpCz2
        call this%spectC%fft(field,ctmpY)
        call this%spectC%testFilter_ip(ctmpY)
        if (useVerticalTfilter) then
            call transpose_y_to_z(ctmpY,ctmpZ1,this%sp_gp)
            call this%Gfilz%filter3(ctmpZ1,ctmpZ2,size(ctmpZ1,1),size(ctmpZ1,2))
            call transpose_z_to_y(ctmpZ2,ctmpY,this%sp_gp)
        end if 
        call this%spectC%ifft(ctmpY,field)
    end subroutine
   
    subroutine testFilter_oop(this, fin, fout)
        class(sgs), intent(inout), target :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in)  :: fin
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: fout
        complex(rkind), dimension(:,:,:), pointer :: ctmpY, ctmpZ1, ctmpZ2
        
        ctmpY => this%cbuff(:,:,:,1)
        ctmpZ1 => this%ctmpCz
        ctmpZ2 => this%ctmpCz2
        call this%spectC%fft(fin,ctmpY)
        call this%spectC%testFilter_ip(ctmpY)
        if (useVerticalTfilter) then
            call transpose_y_to_z(ctmpY,ctmpZ1,this%sp_gp)
            call this%Gfilz%filter3(ctmpZ1,ctmpZ2,size(ctmpZ1,1),size(ctmpZ1,2))
            call transpose_z_to_y(ctmpZ2,ctmpY,this%sp_gp)
        end if 
        call this%spectC%ifft(ctmpY,fout)
    end subroutine
   

    subroutine testFilter_oop_C2R(this, fhat, fout)
        class(sgs), intent(in), target :: this
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: fhat 
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: fout
        complex(rkind), dimension(:,:,:), pointer :: ctmpY, ctmpZ1, ctmpZ2

        ctmpY => this%cbuff(:,:,:,1)
        ctmpZ1 => this%ctmpCz
        ctmpZ2 => this%ctmpCz2
        call this%spectC%TestFilter_oop(fhat, ctmpY)
        if (useVerticalTfilter) then
            call transpose_y_to_z(ctmpY,ctmpZ1,this%sp_gp)
            call this%Gfilz%filter3(ctmpZ1,ctmpZ2,size(ctmpZ1,1),size(ctmpZ1,2))
            call transpose_z_to_y(ctmpZ2,ctmpY,this%sp_gp)
        end if 
        call this%spectC%ifft(ctmpY,fout)

    end subroutine

    subroutine DynamicProcedure(this,u,v,wC,uhat,vhat, wChat, duidxjhat)
        class(sgs), intent(inout), target :: this
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3),9), intent(inout), target :: duidxjhat
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(inout) :: u, v, wC
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: uhat, vhat, wChat

        real(rkind), dimension(:,:,:), pointer :: L11, L12, L13, L22, L23, L33
        real(rkind), dimension(:,:,:), pointer :: M11, M12, M13, M22, M23, M33
        complex(rkind), dimension(:,:,:), pointer :: dudxH, dudyH, dudzH 
        complex(rkind), dimension(:,:,:), pointer :: dvdxH, dvdyH, dvdzH
        complex(rkind), dimension(:,:,:), pointer :: dwdxH, dwdyH, dwdzH
        complex(rkind), dimension(:,:,:), pointer :: ctmpY 
        real(rkind), dimension(:,:,:), pointer :: rtmpX
        real(rkind), dimension(:,:,:), pointer :: numerator, denominator        
        integer :: idx
        real(rkind), dimension(:,:,:), pointer :: C_DYN 


        ! STEP 0: Allocate pointers
        C_DYN => this%rbuff(:,:,:,8); numerator => this%Lij(:,:,:,1); denominator => this%Mij(:,:,:,1)
        
        M11 => this%Mij(:,:,:,1); M12 => this%Mij(:,:,:,2); M13 => this%Mij(:,:,:,3);
        M22 => this%Mij(:,:,:,4); M23 => this%Mij(:,:,:,5); M33 => this%Mij(:,:,:,6);
        
        L11 => this%Lij(:,:,:,1); L12 => this%Lij(:,:,:,2); L13 => this%Lij(:,:,:,3);
        L22 => this%Lij(:,:,:,4); L23 => this%Lij(:,:,:,5); L33 => this%Lij(:,:,:,6);
        
        rtmpX => this%rbuff(:,:,:,8); ctmpY => this%cbuff(:,:,:,2)
        
        dudxH => duidxjhat(:,:,:,1); dudyH => duidxjhat(:,:,:,2); dudzH => duidxjhat(:,:,:,3); 
        dvdxH => duidxjhat(:,:,:,4); dvdyH => duidxjhat(:,:,:,5); dvdzH => duidxjhat(:,:,:,6); 
        dwdxH => duidxjhat(:,:,:,7); dwdyH => duidxjhat(:,:,:,8); dwdzH => duidxjhat(:,:,:,9); 


        ! STEP 1: Compute Mij = delta^2 * \testF{Ssmag*Sij)
        do idx = 1,6
            this%Mij(:,:,:,idx) = this%nuSGS * this%rbuff(:,:,:,idx)
            call this%testFilter_ip(this%Mij(:,:,:,idx))
        end do 
        this%Mij = ( this%deltaFilter * this%deltaFilter ) * this%Mij 
     
        ! STEP 2: Compute Mij = Mij - deltaF^2 * \testF{Ssmag}*\testF{Sij} 
        !         NOTE: Lij is used as a temporary storage
         
        ! Step 2a: Compute \testF{Sij} - Stored in Lij
        call this%testFilter_oop_C2R(dudxH,L11)
        call this%testFilter_oop_C2R(dvdyH,L22)
        call this%testFilter_oop_C2R(dwdzH,L33)
        ctmpY = half*(dudyH + dvdxH)
        call this%testFilter_oop_C2R(ctmpY,L12)
        ctmpY = half*(dudzH + dwdxH)
        call this%testFilter_oop_C2R(ctmpY,L13)
        ctmpY = half*(dvdzH + dwdyH)
        call this%testFilter_oop_C2R(ctmpY,L23)

        ! Step 2b: Compute \testF{Ssmag}
        call this%get_SMAG_Op(rtmpX,L11,L22,L33,L12,L13,L23)
       
        ! Step 2c: multiply and add!
        rtmpX = (this%deltaTfilter * this%deltaTfilter) * rtmpX
        do idx = 1,6
            this%Lij(:,:,:,idx) = rtmpx*this%Lij(:,:,:,idx)
        end do 
        this%Mij = this%Mij - this%Lij
        this%Mij = two*this%Mij 


        ! STEP 3: Compute Lij = \testF{ui*uj}
        L11 = u*u ; L12 = u*v ; L13 = u *wC
        L22 = v*v ; L23 = v*wC; L33 = wC*wC
        do idx = 1,6
            call this%testFilter_ip(this%Lij(:,:,:,idx))
        end do 


        ! STEP 4: Compute Lij = Lij - \testF{ui}*\testF{uj}
        call this%testFilter_oop_C2R(uhat,u)
        call this%testFilter_oop_C2R(vhat,v)
        call this%testFilter_oop_C2R(wChat,wC)
        
        L11 = L11 - u*u ; L12 = L12 - u*v  ; L13 = L13 - u *wC
        L22 = L22 - v*v ; L23 = L23 - v*wC ; L33 = L33 - wC*wC

        ! STEP 5: Compute Numerator = Mij * Lij 
        this%Lij = this%Lij * this%Mij 
        do idx = 2,6
            numerator = numerator + this%Lij(:,:,:,idx)
        end do 
        numerator = numerator + L12 + L13 + L23

        ! STEP 6: Compute Denominator = Mij * Mij 
        this%Mij = this%Mij * this%Mij 
        do idx = 2,6
            denominator = denominator + this%Mij(:,:,:,idx)
        end do 
        denominator = denominator + M12 + M13 + M23

        ! STEP 7: Filter and clip 
        call this%planarAverage(numerator,this%useClipping)
        call this%planarAverage(denominator,.false.)
        denominator = denominator + 1d-14
        
        ! STEP 8: Compute the SMAG constant and the nuSGS
        numerator = numerator/denominator
        !call this%spectC%fft(numerator,this%nuSGShat)
        !call this%spectC%dealias(this%nuSGShat)
        !call this%spectC%ifft(this%nuSGShat,numerator)
        !call this%FiltCoeff(numerator)
        !this%nuSGS =  numerator*(this%deltaFilter * this%deltaFilter) * this%nuSGS
    end subroutine

    subroutine filtCoeff(this, csmag)
        class(sgs), intent(inout) :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(inout) :: csmag


        call transpose_x_to_y(csmag,this%rtmpY,this%gpC)
        call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)

        call this%Gfilz%filter3(this%rtmpZ,this%rtmpZ2,size(this%rtmpZ,1),size(this%rtmpZ,2))

        call transpose_z_to_y(this%rtmpZ2,this%rtmpY,this%gpC)
        call transpose_y_to_x(this%rtmpY,csmag,this%gpC)

    end subroutine
    subroutine planarAverage(this,f, useClipping)
        class(sgs), intent(inout) :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(inout) :: f
        logical, intent(in) :: useClipping
        integer :: k
        real(rkind) :: mnVal

        call transpose_x_to_y(f,this%rtmpY,this%gpC)
        call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)

        do k = 1,this%gpC%zsz(3)
            mnVal = P_SUM(sum(this%rtmpZ(:,:,k)))*this%meanFact
            if (useClipping) then
                this%rtmpZ(:,:,k) = max(mnVal, zero)
            else
                this%rtmpZ(:,:,k) = mnVal
            end if 
        end do 
       
        call transpose_z_to_y(this%rtmpZ,this%rtmpY,this%gpC)
        call transpose_y_to_x(this%rtmpY,f,this%gpC)

    end subroutine


    subroutine planarAverage_oop(this,f, fout)
        class(sgs), intent(inout) :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: f
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: fout
        integer :: k
        real(rkind) :: mnVal

        call transpose_x_to_y(f,this%rtmpY,this%gpC)
        call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)

        do k = 1,this%gpC%zsz(3)
            mnVal = P_SUM(sum(this%rtmpZ(:,:,k)))*this%meanFact
            this%rtmpZ(:,:,k) = mnVal
        end do

        call transpose_z_to_y(this%rtmpZ,this%rtmpY,this%gpC)
        call transpose_y_to_x(this%rtmpY,fout,this%gpC)

    end subroutine

