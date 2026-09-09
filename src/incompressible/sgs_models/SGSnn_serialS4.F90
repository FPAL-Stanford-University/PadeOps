subroutine init_SGSnn(this,NNtype)
   use fortran_assert, only: assert
   class(sgs_igrid), intent(inout) :: this
   integer, intent(in) :: NNtype
   !real(rkind), dimension(:, :, :, :), intent(out), allocatable :: flg1
   !real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),1), intent(out) :: flg1
   allocate(this%flg1(1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))

   this%isEddyViscosityModel = .false.
     
   allocate(this%strain2(       1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%strain2_m(     1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%rot2(          1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%strainrot_m(   1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%rotstrainrot_m(1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%strain2rot_m(  1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%rotstrain2_m(  1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   !allocate(this%finalout( 1,6,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))

   allocate(this%strain2E(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),3,3))
   allocate(this%strain2_mE(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),3,3))
   allocate(this%rot2E(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),3,3))
   allocate(this%strainrot_mE(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),3,3))
   allocate(this%rotstrainrot_mE(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),3,3))
   allocate(this%strain2rot_mE(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),3,3))
   allocate(this%rotstrain2_mE(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),3,3))
   !allocate(this%finaloutE(1,6,this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4)))
   !real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),1), intent(out) :: flag
   

   select case (NNtype)
   case (1) ! Unet
     allocate(this%invariants(1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),7))
     allocate(this%invariantsnorm(1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),7))
     allocate(this%delta(1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),8))
     allocate(this%delta1(1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),8))
     allocate(this%invariantsE(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),7))
     allocate(this%invariantsEnorm(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),7))
     allocate(this%deltaE(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),8))
     allocate(this%delta1E(1,this%gpE%xsz(1),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),this%gpE%xsz(2),8))
     this%flg1(1, 1, 1, 1)=1.0
     call writeInvariants2(0, 0, this%datadir, this%gpC, this%flg1, 'V')
     !call execute_command_line("python interfacePadeOpsCont.py", wait=.false.)
     !call sleep(45)

 case (2) ! LSTM
     allocate(this%invariantsLSTM(1,10,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),7))
     allocate(this%deltaLSTM(     1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),8))
     allocate(this%delta1LSTM(1,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),8))
     allocate(this%invariantsLSTME(1,10,this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),7))
     allocate(this%deltaLSTM(1,this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),8))
     allocate(this%delta1LSTM(1,this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),8))

 case default
       call assert(.false.,'NNtype must be 1 or 2')
   end select

   !call getSijRijforNNmod(duidxjC, this%S_ij_C, this%R_ij_C, this%gpC%xsz(1), & 
   !  this%gpC%xsz(2), this%gpC%xsz(3))
   !call getSijRijforNNmod(duidxjE, this%S_ij_E, this%R_ij_E, this%gpE%xsz(1), & 
   !  this%gpE%xsz(2), this%gpE%xsz(3))
   call message(1,"NN model initialized")
end subroutine

subroutine getSijRijforNNmod(duidxj, Sij, Rij, nxL, nyL, nzL)
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), dimension(nxL,nyL,nzL,9), intent(in) :: duidxj
   real(rkind), dimension(nxL,nyL,nzL,9), intent(out) :: Rij
   real(rkind), dimension(nxL,nyL,nzL,6), intent(out) :: Sij

   integer :: i, j, k

   do k = 1,nzL
      do j = 1,nyL
         !$omp simd
         do i = 1,nxL
            Sij(i,j,k,1) = duidxj(i,j,k,1) ! S11 = dudx
            Sij(i,j,k,2) = 0.5d0*(duidxj(i,j,k,2) + duidxj(i,j,k,4)) ! S12 = 0.5*(dudy + dvdx)
            Sij(i,j,k,3) = 0.5d0*(duidxj(i,j,k,3) + duidxj(i,j,k,7)) ! S13 = 0.5*(dudz + dwdx)
            Sij(i,j,k,4) = duidxj(i,j,k,5) ! S22 = dvdy
            Sij(i,j,k,5) = 0.5d0*(duidxj(i,j,k,6) + duidxj(i,j,k,8)) ! S23 = 0.5*(dvdz + dwdy)
            Sij(i,j,k,6) = duidxj(i,j,k,9) ! S33 = dwdz


            !Sij(i,j,k,5) = duidxj(i,j,k,5) ! S22 = dvdy
            !Sij(i,j,k,6) = 0.5d0*(duidxj(i,j,k,6) + duidxj(i,j,k,8)) ! S23 = 0.5*(dvdz + dwdy)
            
            !Sij(i,j,k,9) = duidxj(i,j,k,9) ! S33 = dwdz

            Rij(i,j,k,2) = 0.5d0*(duidxj(i,j,k,2) - duidxj(i,j,k,4)) ! R12 = 0.5*(dudy - dvdx)
            Rij(i,j,k,3) = 0.5d0*(duidxj(i,j,k,3) - duidxj(i,j,k,7)) ! R13 = 0.5*(dudz - dwdx)
            Rij(i,j,k,6) = 0.5d0*(duidxj(i,j,k,6) - duidxj(i,j,k,8)) ! R23 = 0.5*(dvdz - dwdy)
         end do 
      end do 
   end do 
   !Sij(:,:,:,4) = Sij(:,:,:,2) ! S21 = S12
   !Sij(:,:,:,7) = Sij(:,:,:,3) ! S31 = S13
   !Sij(:,:,:,8) = Sij(:,:,:,6) ! S32 = S23

   Rij(:,:,:,1) = 0.d0
   Rij(:,:,:,5) = 0.d0
   Rij(:,:,:,9) = 0.d0

   Rij(:,:,:,4) = -Rij(:,:,:,2) ! R21 = -R12
   Rij(:,:,:,7) = -Rij(:,:,:,3) ! R31 = -R13
   Rij(:,:,:,8) = -Rij(:,:,:,6) ! R32 = -R23
end subroutine

subroutine getVijforNN(duidxj_C1, duidxj_E1, a_duidxj_C, a_duidxj_E, nxLC, nyLC, nzLC, nxLE, nyLE, nzLE)
    integer, intent(in) :: nxLC, nyLC, nzLC, nxLE, nyLE, nzLE
    real(rkind), dimension(nxLC,nyLC,nzLC,9), intent(in) :: duidxj_C1
    real(rkind), dimension(nxLE,nyLE,nzLE,9), intent(in) :: duidxj_E1
    real(rkind), dimension(nxLC,nyLC,nzLC,9), intent(out) :: a_duidxj_C
    real(rkind), dimension(nxLE,nyLE,nzLE,9), intent(out) :: a_duidxj_E
    !integer :: i, j, k, l
    
    a_duidxj_C = duidxj_C1
    a_duidxj_E = duidxj_E1

    !do l = 1, 9
    !  do k = 1, nzLC
    !    do j = 1, nyLC
    !      do i = 1, nxLC
    !        a_duidxj_C(i, j, k, l) = duidxj_C1(i, j, k, l)
    !      end do
    !    end do
    !  end do
    !end do
    !do l = 1, 9
    !  do k = 1, nzLE
    !    do j = 1, nyLE
    !      do i = 1, nxLE
    !        print*, a_duidxj_E
    !      end do
    !    end do
    !  end do
    !end do
end subroutine

subroutine readInvariants(tid,runID,datadir,gp,dat,arrType)
    integer, intent(in) :: tid, runID
    character(len=*), intent(in) :: datadir
    class(decomp_info), intent(in) :: gp
    real(rkind), dimension(:,:,:,:), intent(out) :: dat
    character(len=1), intent(in) :: arrType
    character(len=clen) :: fname
    integer :: n
   
    do n = 1,size(dat,4) 
        write(fname,'(A,I2.2,A5,I1,A2,I6.6,A4)') trim(datadir)//'/Run',runID,'_inv'//arrType,n,'_t',tid,'.out'
        call decomp_2d_read_one(1,dat(n,:,:,:),trim(fname),gp)
    end do
end subroutine

subroutine writeInvariants(tid,runID,datadir,gp,dat,arrType)
    integer, intent(in) :: tid, runID
    character(len=*), intent(in) :: datadir
    class(decomp_info), intent(in) :: gp
    real(rkind), dimension(:,:,:,:), intent(out) :: dat
    character(len=1), intent(in) :: arrType
    character(len=clen) :: fname
    integer :: n
   
    do n = 1,size(dat,4) 
        write(fname,'(A,I2.2,A5,I1,A2,I6.6,A4)') trim(datadir)//'/Run',runID,'_inv'//arrType,n,'_t',tid,'.out'
        call decomp_2d_write_one(1,dat(:,:,:,n),trim(fname),gp)
    end do
end subroutine


subroutine writeInvariants2(tid,runID,datadir,gp,dat,arrType)
    integer, intent(in) :: tid, runID
    character(len=*), intent(in) :: datadir
    class(decomp_info), intent(in) :: gp
    real(rkind), dimension(:,:,:,:), intent(out) :: dat
    character(len=1), intent(in) :: arrType
    character(len=clen) :: fname
    integer :: n
    do n = 1,size(dat,1)
        write(fname,'(A,I2.2,A5,I1,A2,I6.6,A4)') trim(datadir)//'/Run',runID,'_inv'//arrType,n,'_t',tid,'.out'
        call decomp_2d_write_one(1,dat(n,:,:,:),trim(fname),gp)
    end do
end subroutine

subroutine writeInvariants3(tid,runID,datadir,gp,dat,arrType)
    integer, intent(in) :: tid, runID
    character(len=*), intent(in) :: datadir
    class(decomp_info), intent(in) :: gp
    real(rkind), dimension(:,:,:,:,:), intent(out) :: dat
    character(len=1), intent(in) :: arrType
    character(len=clen) :: fname
    integer :: n, j, k
    k = 1
    do n = 1,size(dat,4)
        do j = 1,size(dat,5)
            write(fname,'(A,I2.2,A5,I1,A2,I6.6,A4)') trim(datadir)//'/Run',runID,'_inv'//arrType,k,'_t',tid,'.out'
            call decomp_2d_write_one(1,dat(:,:,:,n,j),trim(fname),gp)
            k = k+1
        end do
    end do
end subroutine



subroutine compute_tauij_NN(this,tidNow)
  !use, intrinsic :: iso_c_binding
  !implicit none
  !type(c_ptr) :: c_pstruc
  !real(c_float), pointer :: f_pstruc(:) => NULL()
  !!real(c_float), pointer, dimension(:) :: f_pstruc
  !type(c_ptr) :: c_pmag
  !real(c_float), pointer :: f_pmag(:) => NULL()
  !!real(c_float), pointer, dimension(:) :: f_pmag
  !type(c_ptr) :: c_pEstruc
  !real(c_float), pointer :: f_pEstruc(:) => NULL()
  !!real(c_float), pointer, dimension(:) :: f_pEstruc
  !type(c_ptr) :: c_pEmag
  !real(c_float), pointer :: f_pEmag(:) => NULL()
  !!real(c_float), pointer, dimension(:) :: f_pEmag


  
  class(sgs_igrid), intent(inout) :: this
  integer, intent(in) :: tidNow
  integer :: tidPast, n
  integer :: i, j, k, l, m
  integer :: count
  real(rkind), dimension(3, 3) :: arr_temp1
  real(rkind), dimension(3, 3) :: arr_temp2
  real(rkind), dimension(3, 3) :: arr_temp3
  real(rkind), dimension(3, 3) :: arr_temp4
  real(rkind), dimension(3, 3) :: arr_temp5
  real(rkind), dimension(3, 3) :: arr_temp6
  real(rkind), dimension(9) :: S_ij_a
  real(rkind), dimension(9) :: R_ij_a
  integer, dimension(9) :: timelist
  real(rkind) :: del
  real(rkind) :: tscale
  real(rkind) :: a_1
  real(rkind) :: a_2
  real(rkind) :: f
  !real(rkind) :: trace

  real(rkind), dimension(6, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: finalout
  !real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: finalout1
  !real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: finalout2
  !real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: finalout3
  !real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: finalout4
  !real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: finalout5
  !real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: finalout6
  !real(rkind), dimension(1, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: magout
  real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: trace

  real(rkind), dimension(6, this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)) :: finaloutE
  !real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4)) :: finaloutE1
  !real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4)) :: finaloutE2
  !real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4)) :: finaloutE3
  !real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4)) :: finaloutE4
  !real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4)) :: finaloutE5
  !real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4)) :: finaloutE6
  real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4),9) :: velgradE
  !real(rkind), dimension(1, this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4)) :: magoutE
  real(rkind), dimension(1,this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: flag_f


  !a_1 = this%dy/this%dz
  !a_2 = this%dz/this%dx
  !f = cosh(sqrt(4 / 27 * (log(a_1)**2 - log(a_1)*log(a_2) + log(a_2)**2)))
  !del = (this%dx*this%dy*this%dz)**(2.0/3.0)*f**(2) 
  
  ! TODO: Andy, compute your invariants, call your NN, compute tauij
  call writeInvariants(0, 0, this%datadir, this%gpC, this%duidxj_C, 'W')
    
  !print*, "Fortran C End"
  !print*, this%duidxj_C(1, 1, 1, 1)
  !print*, this%duidxj_C(1, 1, 2, 1)
  !print*, this%duidxj_C(1, 1, 3, 1)
  !print*, this%duidxj_C(20, 32, 31, 2)
  !print*, this%duidxj_C(22, 38, 18, 3)
  !print*, this%duidxj_C(26, 44, 21, 4)

  !print*, "Fortran E End"
  !print*, this%duidxj_E(1, 1, 1, 1)
  !print*, this%duidxj_E(1, 1, 2, 1)
  !print*, this%duidxj_E(1, 1, 3, 1)
  !print*, this%duidxj_E(20, 32, 31, 2)
  !print*, this%duidxj_E(22, 38, 18, 3)
  !print*, this%duidxj_E(26, 44, 21, 4)

  !do i = 61,this%gpE%xsz(1)
  !  do j = 61, this%gpE%xsz(2)
  !    do k = 61, this%gpE%xsz(3)
  !      do l = 1, 9
  !         print*, this%duidxj_E(i, j, k, l)
  !      end do
  !    end do
  !  end do
  !end do

  
  call writeInvariants(0, 0, this%datadir, this%gpE, this%duidxj_E, 'X')
  !flag_f(1, 1, 1, 1)=0.0
  !call writeInvariants2(0, 0, this%datadir, this%gpC, flag_f, 'V')
  
  !print*, "invars written" 
  call execute_command_line("python interfacePadeOps.py", wait=.true.)
  !print*, "python end"
  call readInvariants(0, 0, this%datadir, this%gpC, finalout, 'Y')
  call readInvariants(0, 0, this%datadir, this%gpE, finaloutE, 'Z')

  !do
  !  !print*, "entered and waiting"
  !  !print*, int(flag_f(1, 1, 1, 1))
  !  !call sleep(2)
  !  !print*, "done sleep"
  !  call readInvariants(0, 0, this%datadir, this%gpC, flag_f, 'V')
  !  !print*, int(flag_f(1, 1, 1, 1))
  !  IF (int(flag_f(1, 1, 1, 1))==1) EXIT
  !end do
  !! ensure that what we read was correct just to be sure
  !call sleep(1)
  !do
  !  call readInvariants(0, 0, this%datadir, this%gpC, flag_f, 'V')
  !  IF (int(flag_f(1, 1, 1, 1))==1) EXIT
  !end do
  !call readInvariants(0, 0, this%datadir, this%gpC, finalout, 'Y')
  !call readInvariants(0, 0, this%datadir, this%gpE, finaloutE, 'Z') 
  
  !do i = 62,this%gpC%xsz(1)
  !  do j = 62, this%gpC%xsz(2)
  !    do k = 62, this%gpC%xsz(3)
  !      do l = 1, 6
  !         print*, finalout(l, i, j, k)
  !      end do
  !    end do
  !  end do
  !end do

  !print*, "Fortran C Read"
  !print*, finalout(1, 1, 1, 1)
  !print*, finalout(1, 1, 2, 1)
  !print*, finalout(1, 1, 3, 1)
  !print*, finalout(2, 20, 32, 31)
  !print*, finalout(3, 22, 38, 18)
  !print*, finalout(4, 26, 44, 21)

  !print*, "Fortran E Read"
  !print*, finaloutE(1, 1, 1, 1)
  !print*, finaloutE(1, 1, 2, 1)
  !print*, finaloutE(1, 1, 3, 1)
  !print*, finaloutE(2, 20, 32, 31)
  !print*, finaloutE(3, 22, 38, 18)
  !print*, finaloutE(4, 26, 44, 21)

  
  trace = finalout(1,:,:,:)+finalout(4,:,:,:)+finalout(6,:,:,:)
  this%tau_11 = finalout(1,:,:,:)-1/3*(trace)
  this%tau_12 = finalout(2,:,:,:)
  this%tau_13 = finaloutE(3,1:this%gpE%xsz(1),1:this%gpE%xsz(2),1:this%gpE%xsz(3))
  this%tau_22 = finalout(4,:,:,:)-1/3*(trace)
  this%tau_23 = finaloutE(5,1:this%gpE%xsz(1),1:this%gpE%xsz(2),1:this%gpE%xsz(3))
  this%tau_33 = finalout(6,:,:,:)-1/3*(trace)

  !call writeInvariants2(tidNow,this%runID,this%datadir,this%gpC,finalout(1,:,:,:,:),'C')
  !call writeInvariants2(tidNow,this%runID,this%datadir,this%gpE,finaloutE(1,:,:,:,:),'E')
  !call writeInvariants(tidNow,this%runID,this%datadir,this%gpC,this%S_ij_C(:,:,:,:),'A')
  !call writeInvariants(tidNow,this%runID,this%datadir,this%gpE,this%S_ij_E(:,:,:,:),'B')
  !call writeInvariants(tidNow,this%runID,this%datadir,this%gpC,this%R_ij_C(:,:,:,:),'F') 
  !call writeInvariants(tidNow,this%runID,this%datadir,this%gpE,this%R_ij_E(:,:,:,:),'G')

  !call writeInvariants(tidNow,this%runID,this%datadir,this%gpC,this%invariants(1,:,:,:,:),'H')
  !call writeInvariants(tidNow,this%runID,this%datadir,this%gpC,this%invariantsnorm(1,:,:,:,:),'P')
  !call writeInvariants3(tidNow,this%runID,this%datadir,this%gpC,this%strain2(1,:,:,:,:,:),'I')
  !call writeInvariants3(tidNow,this%runID,this%datadir,this%gpC,this%strain2_m(1,:,:,:,:,:),'J')
  !call writeInvariants3(tidNow,this%runID,this%datadir,this%gpC,this%rot2(1,:,:,:,:,:),'K')
  !call writeInvariants3(tidNow,this%runID,this%datadir,this%gpC,this%strainrot_m(1,:,:,:,:,:),'L')
  !call writeInvariants3(tidNow,this%runID,this%datadir,this%gpC,this%rotstrainrot_m(1,:,:,:,:,:),'M')
  !call writeInvariants3(tidNow,this%runID,this%datadir,this%gpC,this%strain2rot_m(1,:,:,:,:,:),'N')
  !call writeInvariants3(tidNow,this%runID,this%datadir,this%gpC,this%rotstrain2_m(1,:,:,:,:,:),'O')



  !call this%viz_hdf5%write_variable(this%tau_12, "tau_12")
  !call this%viz_hdf5%write_variable(this%tau_13, "tau_13")
  !call this%viz_hdf5%write_variable(this%tau_23, "tau_23")
  !call this%viz_hdf5%write_variable(this%tau_22, "tau_22")
  !call this%viz_hdf5%write_variable(this%tau_33, "tau_33")
  !call this%viz_hdf5%write_variable(this%S_ij_C(:,:,:,1), "SijC_11")
  !call this%viz_hdf5%write_variable(this%S_ij_C(:,:,:,2), "SijC_12")
  !call this%viz_hdf5%write_variable(this%S_ij_C(:,:,:,5), "SijC_22")
  !call this%viz_hdf5%write_variable(this%S_ij_C(:,:,:,9), "SijC_33")
  !call this%viz_hdf5%write_variable(this%S_ij_E(:,:,:,3), "SijE_13")
  !call this%viz_hdf5%write_variable(this%S_ij_E(:,:,:,6), "SijE_23")


  !do i = 1,this%gpC%xsz(1)
  !  do j = 1,this%gpC%xsz(2)
  !    print*, this%strain2(1,i,j,1,1,1)
  !  end do
  !end do

  !print *, "mags"
  !do i = 1,this%gpC%xsz(1)
  !  do j = 1,this%gpC%xsz(2)
  !    print*, magout(1,i,j,1)
  !  end do
  !end do

  !print *, "struc"
  !do i = 1,this%gpC%xsz(1)
  !  do j = 1,this%gpC%xsz(2)
  !    print*, finalout(1,1,i,j,1)
  !  end do
  !end do


end subroutine
