!! FFTPACK version of NUFFT adapted from the NYU NUFFT Code at: 
!! http://www.cims.nyu.edu/cmcl/nufft/nufft.html. Avoid using this version, and
!! use the FFTW version instead. 


module nufftMod
    use kind_parameters, only: rkind
    use constants, only: pi, zero
    use timer, only: tic, toc
    implicit none
    private 
    public :: nufft 

    type :: nufft
        private
        integer :: nx, ny, nz
        integer :: nw1, nw2, nw3, iw10, iw11, iw12, iw13, iw14, iw15, iw16
        real(rkind) :: hx, hy, hz
        real(rkind) :: t1, t2, t3
        real(rkind), dimension(:), allocatable :: fw, fwstore
        integer :: nspread, nf1, nf2, nf3
        
        contains
            procedure :: init
            procedure :: destroy
            procedure :: fft 
    end type

contains
    subroutine init(this,eps,ms,mt,mu)
        class(nufft), intent(inout) :: this
        real(rkind), intent(in) :: eps
        integer, intent(in) :: ms, mt, mu
        real(rkind) :: rat, r2lamb

        real(rkind) :: cross1  
        integer :: ii, k1

        if ( (eps < 1.d-13) .or. (eps > 1.d-1)) then
            print*, "Cannot initialize nufft"
            stop 
        end if 

        if (eps.le.1d-11) then
           rat = 3.0d0
        else
           rat = 2.0d0
        endif

        this%nx = ms
        this%ny = mt
        this%nz = mu
        this%nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0) + 1
        this%nf1 = rat*ms
        this%nf2 = rat*mt
        this%nf3 = rat*mu
        
         
        if (2*this%nspread.gt.min(this%nf1,this%nf2,this%nf3)) then
           ii = next235(2d0*this%nspread)
           this%nf1 = max(this%nf1, ii)
           this%nf2 = max(this%nf2, ii)
           this%nf3 = max(this%nf3, ii)
        endif
        
        r2lamb = rat*rat * this%nspread / (rat*(rat-.5d0))
        this%hx = 2*pi/this%nf1
        this%hy = 2*pi/this%nf2
        this%hz = 2*pi/this%nf3
        this%nw1 = this%nf1*this%nf2*this%nf3 
        this%nw2 = this%nw1 + max(this%nf2,this%nf3)
        this%nw3 = this%nw2 + max(this%nf1,this%nf2,this%nf3)
        this%iw10 = 2*this%nw3
        this%iw11 = this%iw10 + ms
        this%iw12 = this%iw11 + mt
        this%iw13 = this%iw12 + mu/2 + 1
        this%iw14 = this%iw13 + 48+16+2*this%nf1
        this%iw15 = this%iw14 + 48+16+2*this%nf2
        this%iw16 = this%iw15 + 48+16+2*this%nf3
        allocate( this%fw(0:this%iw16-1) )
        allocate( this%fwstore(0:this%iw16-1) )
        this%t1 = pi/r2lamb
        this%t2 = pi/r2lamb
        this%t3 = pi/r2lamb
        do k1 = 1, this%nspread
           this%fw(this%iw13+k1) = exp(-this%t1*k1**2)
           this%fw(this%iw14+k1) = exp(-this%t2*k1**2)
           this%fw(this%iw15+k1) = exp(-this%t3*k1**2)
        enddo

        call dcfti1(this%nf1,this%fw(this%iw13+64),this%fw(this%iw13+48))
        call dcfti1(this%nf2,this%fw(this%iw14+64),this%fw(this%iw14+48))
        call dcfti1(this%nf3,this%fw(this%iw15+64),this%fw(this%iw15+48))
        this%t1 = pi * r2lamb / dble(this%nf1)**2
        cross1 = 1d0-2d0*mod(ms/2,2)
        do k1 = -ms/2, (ms-1)/2
           this%fw(this%iw10+k1+ms/2)=cross1*exp(this%t1*dble(k1)**2)
           cross1 = -cross1
        enddo
        this%t2 = pi * r2lamb / dble(this%nf2)**2
        cross1 = 1d0-2d0*mod(mt/2,2)
        do k1 = -mt/2, (mt-1)/2
           this%fw(this%iw11+k1+mt/2)=cross1*exp(this%t2*dble(k1)**2)
           cross1 = -cross1
        enddo
        this%t3 = pi * r2lamb / dble(this%nf3)**2
        cross1 = 1d0 / (sqrt(r2lamb)*r2lamb)
        do k1 = 0, mu/2
           this%fw(this%iw12+k1)=cross1*exp(this%t3*dble(k1)**2)
           cross1 = -cross1
        enddo
        do k1 = 0, 2*this%nf1*this%nf2*this%nf3-1
           this%fw(k1) = 0d0
        enddo
        this%t1 = pi/r2lamb
        this%t2 = pi/r2lamb
        this%t3 = pi/r2lamb

        this%fwstore = this%fw
    end subroutine 
   

    subroutine fft(this,nj,xj,yj,zj,cj,fk)
        class(nufft), intent(inout) :: this
        integer, intent(in) :: nj
        real(rkind), dimension(nj), intent(in) :: xj, yj, zj
        complex(rkind), dimension(nj), intent(in) :: cj
        complex(rkind), dimension(-this%nx/2:(this%nx-1)/2,-this%ny/2:(this%ny-1)/2,-this%nz/2:(this%nz-1)/2), intent(out) :: fk 
        
        real(rkind) :: xc(-47:47),yc(-47:47),zc(-47:47)
        integer :: j, k1, k2, k3, i2(-47:47)
        real(rkind) :: cross1, cross, diff1,diff2,diff3
        complex(rkind) :: ccj, cc,c2,zz
        integer :: jb1,jb1u,jb1d,jb2,jb3, ii, istart,i3, is2 
        integer, parameter :: iflag = 1

        this%fw = this%fwstore
        do j = 1, nj
           ccj = cj(j)/dble(nj)
           jb1 = int((xj(j)+pi)/this%hx)
           diff1 = (xj(j)+pi)/this%hx - jb1
           jb1 = mod(jb1, this%nf1)
           if (jb1.lt.0) jb1=jb1+this%nf1
           jb2 = int((yj(j)+pi)/this%hy)
           diff2 = (yj(j)+pi)/this%hy - jb2
           jb2 = mod(jb2, this%nf2)
           if (jb2.lt.0) jb2=jb2+this%nf2
           jb3 = int((zj(j)+pi)/this%hz)
           diff3 = (zj(j)+pi)/this%hz - jb3
           jb3 = mod(jb3, this%nf3)
           if (jb3.lt.0) jb3=jb3+this%nf3
           xc(0) = exp(-this%t1*diff1**2-this%t2*diff2**2-this%t3*diff3**2)
           cross = xc(0)
           cross1 = exp(2d0*this%t1 * diff1)
           do k1 = 1, this%nspread
              cross = cross * cross1
              xc(k1) = this%fw(this%iw13+k1)*cross
           enddo
           cross = xc(0)
           cross1 = 1d0/cross1
           do k1 = 1, this%nspread-1
              cross = cross * cross1
              xc(-k1) = this%fw(this%iw13+k1)*cross
           enddo
           yc(0) = 1d0
           cross = exp(2d0*this%t2 * diff2)
           cross1 = cross
           do k2 = 1, this%nspread-1
              yc(k2) = this%fw(this%iw14+k2)*cross
              yc(-k2) = this%fw(this%iw14+k2)/cross
              cross = cross * cross1
           enddo
           yc(this%nspread) = this%fw(this%iw14+this%nspread)*cross
           zc(0) = 1d0
           cross = exp(2d0*this%t3 * diff3)
           cross1 = cross
           do k3 = 1, this%nspread-1
              zc(k3) = this%fw(this%iw15+k3)*cross
              zc(-k3) = this%fw(this%iw15+k3)/cross
              cross = cross * cross1
           enddo
           zc(this%nspread) = this%fw(this%iw15+this%nspread)*cross
           jb1d = min(this%nspread-1, jb1)
           jb1u = min(this%nspread, this%nf1-jb1-1)
           do k2 = -this%nspread+1, this%nspread
              i2(k2) = jb2+k2
              if (i2(k2).lt.0) then
                 i2(k2) = i2(k2) + this%nf2
              elseif (i2(k2).ge.this%nf2) then
                 i2(k2) = i2(k2) - this%nf2
              endif
           enddo
           do k3 = -this%nspread+1, this%nspread
              i3 = jb3+k3
              if (i3.lt.0) then
                 i3 = i3 + this%nf3
              elseif (i3.ge.this%nf3) then
                 i3 = i3 - this%nf3
              endif
              c2 = zc(k3)*ccj
              do k2 = -this%nspread+1, this%nspread
                 cc = yc(k2)*c2
                 ii = jb1 + i2(k2)*this%nf1 + i3*this%nf1*this%nf2 ! cfine(ib,jb+k2,kb+k3)
                 do k1 = -this%nspread+1, -jb1d-1
	                istart = 2*(k1+this%nf1+ii)
	                zz = xc(k1)*cc
                    this%fw(istart) = this%fw(istart) + dreal(zz)
                    this%fw(istart+1) = this%fw(istart+1) + dimag(zz)
                 enddo
                 do k1 = -jb1d, jb1u
	                istart = 2*(k1+ii)
	                zz = xc(k1)*cc
                    this%fw(istart) = this%fw(istart) + dreal(zz)
                    this%fw(istart+1) = this%fw(istart+1) + dimag(zz)
                 enddo
                 do k1 = jb1u+1, this%nspread
	                istart = 2*(k1-this%nf1+ii)
	                zz = xc(k1)*cc
                    this%fw(istart) = this%fw(istart) + dreal(zz)
                    this%fw(istart+1) = this%fw(istart+1) + dimag(zz)
                 enddo
              enddo
           enddo
        enddo


        i3 = this%iw13 + 48
        do k3 = 0, this%nf3-1
           do k2 = 0, this%nf2-1
              ii = this%nf1 * (k2+k3*this%nf2) 
!              if (iflag .ge. 0) then
                 call dcftb1(this%nf1,this%fw(2*ii),this%fw(2*this%nw2),this%fw(i3+16),this%fw(i3))
!              else
!                 call dcftf1(this%nf1,this%fw(2*ii),this%fw(2*this%nw2),this%fw(i3+16),this%fw(i3))
!              endif
           enddo
        enddo
        i3 = this%iw14 + 48
        do k3 = 0, this%nf3-1
           do k1 = 0, this%nf1-1
              ii = k1 + k3 * this%nf1*this%nf2 
              do k2 = 0, this%nf2-1
                 istart = 2*(this%nw1+k2)
                 is2 = 2*(ii + k2*this%nf1)
                 this%fw(istart) = this%fw(is2)
                 this%fw(istart+1) = this%fw(is2+1)
              enddo
              !if (iflag .ge. 0) then
                 call dcftb1(this%nf2,this%fw(2*this%nw1),this%fw(2*this%nw2),this%fw(i3+16),this%fw(i3))
              !else
              !   call dcftf1(this%nf2,this%fw(2*this%nw1),this%fw(2*this%nw2),this%fw(i3+16),this%fw(i3))
              !endif
              do k2 = 0, this%nf2-1
                 istart = 2*(ii + k2*this%nf1)
                 is2 = 2*(this%nw1+k2)
                 this%fw(istart) = this%fw(is2)
                 this%fw(istart+1) = this%fw(is2+1)
              enddo
           enddo
        enddo

        i3 = this%iw15+48
        do k2 = -this%ny/2,  (this%ny-1)/2
           do k1 = -this%nx/2, (this%nx-1)/2
              ii = k1
              if (k1.lt.0) ii = ii+this%nf1
              ii = ii + k2*this%nf1
              if (k2.lt.0) ii = ii+this%nf1*this%nf2
              do k3 = 0, this%nf3-1
                 istart = 2*(this%nw1+k3)
                 is2 = 2*(ii + k3*this%nf1*this%nf2)
                 this%fw(istart) = this%fw(is2)
                 this%fw(istart+1) = this%fw(is2+1)
              enddo
              !if (iflag .ge. 0) then
                 call dcftb1(this%nf3,this%fw(2*this%nw1),this%fw(2*this%nw2),this%fw(i3+16),this%fw(i3))
              !else
              !   call dcftf1(this%nf3,this%fw(2*this%nw1),this%fw(2*this%nw2),this%fw(i3+16),this%fw(i3))
              !endif
              cross = this%fw(this%iw10+k1+this%nx/2) * this%fw(this%iw11+k2+this%ny/2)
	          zz = dcmplx(this%fw(2*this%nw1),this%fw(2*this%nw1+1))
              fk(k1, k2, 0) = (cross*this%fw(this%iw12))*zz
              do k3 = 1, (this%nz-1)/2
                 istart = 2*(this%nw1+k3)
	             zz = dcmplx(this%fw(istart),this%fw(istart+1))
                 fk(k1,k2,k3) = (cross*this%fw(this%iw12+k3))*zz
                 istart = 2*(this%nw1+this%nf3-k3)
	             zz = dcmplx(this%fw(istart),this%fw(istart+1))
                 fk(k1,k2,-k3) = (cross*this%fw(this%iw12+k3))*zz
              enddo
              if (this%nz/2*2.eq.this%nz) then 
                 istart = 2*(this%nw1+this%nf3-this%nz/2)
	             zz = dcmplx(this%fw(istart),this%fw(istart+1))
                 fk(k1,k2,-this%nz/2) =  (cross*this%fw(this%iw12+this%nz/2))*zz
	          endif
           enddo
        enddo

    end subroutine 


    subroutine destroy(this)
        class(nufft), intent(inout) :: this
        deallocate(this%fw, this%fwstore)
    end subroutine 
    
    pure function next235(base)
        implicit none
        integer next235, numdiv
        real(rkind), intent(in) ::  base
        
        next235 = 2 * int(base/2d0+.9999d0)
        if (next235.le.0) next235 = 2

100     numdiv = next235
        do while (numdiv/2*2 .eq. numdiv)
           numdiv = numdiv /2
        enddo
        do while (numdiv/3*3 .eq. numdiv)
           numdiv = numdiv /3
        enddo
        do while (numdiv/5*5 .eq. numdiv)
           numdiv = numdiv /5
        enddo
        if (numdiv .eq. 1) return
        next235 = next235 + 2
        goto 100
    end function

end module 
