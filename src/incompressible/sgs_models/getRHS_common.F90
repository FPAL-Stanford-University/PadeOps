        call this%spect%fft(tau11,tauhat)
        call this%spect%mtimes_ik1_ip(tauhat)
        urhs = urhs + tauhat

        call this%spect%fft(tau22,tauhat)
        call this%spect%mtimes_ik2_ip(tauhat)
        vrhs = vrhs + tauhat
        
        call this%spect%fft(tau33,tauhat)
        call transpose_y_to_z(tauhat,this%ctmpCz,this%sp_gp)
        call this%derZ%ddz_C2E(this%ctmpCz,this%ctmpEz,size(this%ctmpCz,1),size(this%ctmpCz,2))
        call transpose_z_to_y(this%ctmpEz,this%ctmpEy,this%sp_gpE)
        wrhs = wrhs + this%ctmpEy

        call this%spect%fft(tau12,tauhat)
        call this%spect%mtimes_ik2_oop(tauhat,tauhat2)
        urhs = urhs + tauhat2
        call this%spect%mtimes_ik1_ip(tauhat)
        vrhs = vrhs + tauhat

        call this%spect%fft(tau13,tauhat)
        call transpose_y_to_z(tauhat,this%ctmpCz,this%sp_gp)
        call this%derZ%InterpZ_C2E(this%ctmpCz,this%ctmpEz,size(this%ctmpCz,1),size(this%ctmpCz,2))
        call transpose_z_to_y(this%ctmpEz,this%ctmpEy,this%sp_gpE)
        call this%spectE%mtimes_ik1_ip(this%ctmpEy)
        wrhs = wrhs + this%ctmpEy
       
        call this%derZ%ddz_C2C(this%ctmpCz,this%ctmpCz2, size(this%ctmpCz,1),size(this%ctmpCz,2)) 
        call transpose_z_to_y(this%ctmpCz2,tauhat,this%sp_gp)
        urhs = urhs + tauhat
    
        
        call this%spect%fft(tau23,tauhat)
        call transpose_y_to_z(tauhat,this%ctmpCz,this%sp_gp)
        call this%derZ%InterpZ_C2E(this%ctmpCz,this%ctmpEz,size(this%ctmpCz,1),size(this%ctmpCz,2))
        call transpose_z_to_y(this%ctmpEz,this%ctmpEy,this%sp_gpE)
        call this%spectE%mtimes_ik2_ip(this%ctmpEy)
        wrhs = wrhs + this%ctmpEy

        call this%derZ%ddz_C2C(this%ctmpCz,this%ctmpCz2, size(this%ctmpCz,1),size(this%ctmpCz,2)) 
        call transpose_z_to_y(this%ctmpCz2,tauhat,this%sp_gp)
        vrhs = vrhs + tauhat
