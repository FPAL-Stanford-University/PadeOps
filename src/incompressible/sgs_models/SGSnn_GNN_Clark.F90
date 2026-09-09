subroutine init_SGSnn(this,NNtype)
   use fortran_assert, only: assert
   class(sgs_igrid), intent(inout) :: this
   integer, intent(in) :: NNtype

   this%isEddyViscosityModel = .false.
     
   allocate(this%strain2(       1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%strain2_m(     1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%rot2(          1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%rot(           1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%strainrot_m(   1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%rotstrainrot_m(1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%strain2rot_m(  1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   allocate(this%rotstrain2_m(  1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),3,3))
   !allocate(this%finalout( 1,6,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))

   allocate(this%strain2E(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),3,3))
   allocate(this%strain2_mE(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),3,3))
   allocate(this%rot2E(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),3,3))
   allocate(this%rotE(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),3,3))
   allocate(this%strainrot_mE(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),3,3))
   allocate(this%rotstrainrot_mE(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),3,3))
   allocate(this%strain2rot_mE(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),3,3))
   allocate(this%rotstrain2_mE(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),3,3))
   !allocate(this%finaloutE(1,6,this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4)))
   

   select case (NNtype)
   case (1) ! Unet
     allocate(this%invariants(1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),5))
     allocate(this%invariantsnorm(1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),5))
     allocate(this%delta(1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),8))
     allocate(this%delta1(1,this%gpC%xsz(1),this%gpC%xsz(3),this%gpC%xsz(2),8))
     allocate(this%invariantsE(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),5))
     allocate(this%invariantsEnorm(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),5))
     allocate(this%deltaE(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),8))
     allocate(this%delta1E(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),8))

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
        call decomp_2d_read_one(1,dat(:,:,:,n),trim(fname),gp)
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
  use, intrinsic :: iso_c_binding
  implicit none
  type(c_ptr) :: c_pstruc
  real(c_float), pointer :: f_pstruc(:) => NULL()
  type(c_ptr) :: c_pEstruc
  real(c_float), pointer :: f_pEstruc(:) => NULL()
  
  
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

  interface
      function loadpytorchgpu(train_loadf, strain2f, strain2_mf, rot2f, strainrot_mf, rotstrainrot_mf, strain2rot_mf, &
        rotstrain2_mf, rotf, dom1, dom2, dom3, spac1f, spac2f, spac3f) bind(c)
        import
        type(c_ptr) :: loadpytorchgpu
          real(rkind), dimension(2), intent(in) :: train_loadf
          real(rkind), dimension(2), intent(in) :: strain2f
          real(rkind), dimension(2), intent(in) :: strain2_mf
          real(rkind), dimension(2), intent(in) :: rot2f
          real(rkind), dimension(2), intent(in) :: rotf
          real(rkind), dimension(2), intent(in) :: strainrot_mf
          real(rkind), dimension(2), intent(in) :: rotstrainrot_mf
          real(rkind), dimension(2), intent(in) :: strain2rot_mf
          real(rkind), dimension(2), intent(in) :: rotstrain2_mf
          real(rkind), intent(in) :: spac1f
          real(rkind), intent(in) :: spac2f
          real(rkind), intent(in) :: spac3f
          integer, intent(in) :: dom1
          integer, intent(in) :: dom2
          integer, intent(in) :: dom3

      end function loadpytorchgpu
  

      subroutine c_free(ptr) bind(c,name="free")
          import
          type(c_ptr), value :: ptr
      end subroutine c_free
  end interface


  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*11) :: train_loadin
  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*11) :: train_loadinnorm
  real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9) :: tgrad
  !real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*8) :: deltain1
  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*7*10) :: train_loadinLSTM
  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*11) :: train_loadinE
  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*11) :: train_loadinEnorm
  real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),9) :: tgradE
  !real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*8) :: deltain1E
  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*7*10) :: train_loadinLSTME


  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: strain2in
  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: strain2inmag
  real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: strain2norm 
  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: strain2_min
  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: rot2in
  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: rotin
  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: strainrot_min
  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: rotstrainrot_min
  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: strain2rot_min
  real(rkind), dimension(this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9) :: rotstrain2_min
  real(rkind), dimension(1, 6, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: finalout
  real(rkind), dimension(1, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: magout
  real(rkind), dimension(1, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: tgradmag
  real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: trace

  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*9) :: strain2inE
  real(rkind), dimension(1,this%gpE%xsz(1),this%gpE%xsz(3),this%gpE%xsz(2),3,3) :: strain2modE
  real(rkind),dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*9) :: strain2inmagE
  real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)) :: strain2normE
  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*9) :: strain2_minE
  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*9) :: rot2inE
  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*9) :: rotinE
  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*9) :: strainrot_minE
  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*9) :: rotstrainrot_minE
  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*9) :: strain2rot_minE
  real(rkind), dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)*9) :: rotstrain2_minE
  !real, dimension(this%gpE%xsz(1)*this%gpE%xsz(2)*(this%gpE%xsz(3)+4-MOD(this%gpE%xsz(3),4))*7) :: f_pEstruc2

  real(rkind), dimension(1, 6, this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)) :: finaloutE
  real(rkind), dimension(1, this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)) :: magoutE
  real(rkind), dimension(1, this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)) :: tgradmagE




  !a_1 = this%dy/this%dz
  !a_2 = this%dz/this%dx
  a_1 = this%dz/this%dx
  a_2 = this%dy/this%dx
  f = cosh(sqrt(4.0 / 27.0 * (log(a_1)**2 - log(a_1)*log(a_2) + log(a_2)**2)))
  del = 1.0/12.0*((this%dx*this%dy*this%dz)**(1.0/3.0)*f)**2.0
  
  do k = 1, this%gpC%xsz(3)
    do j = 1, this%gpC%xsz(2)
      do i = 1, this%gpC%xsz(1)
        S_ij_a(1) = this%S_ij_C(i, j, k, 1)
        S_ij_a(2) = this%S_ij_C(i, j, k, 2)
        S_ij_a(3) = this%S_ij_C(i, j, k, 3)
        S_ij_a(4) = this%S_ij_C(i, j, k, 2)
        S_ij_a(5) = this%S_ij_C(i, j, k, 4)
        S_ij_a(6) = this%S_ij_C(i, j, k, 5)
        S_ij_a(7) = this%S_ij_C(i, j, k, 3)
        S_ij_a(8) = this%S_ij_C(i, j, k, 5)
        S_ij_a(9) = this%S_ij_C(i, j, k, 6)

        R_ij_a(1) = this%R_ij_C(i, j, k, 1)
        R_ij_a(2) = this%R_ij_C(i, j, k, 2)
        R_ij_a(3) = this%R_ij_C(i, j, k, 3)
        R_ij_a(4) = this%R_ij_C(i, j, k, 4)
        R_ij_a(5) = this%R_ij_C(i, j, k, 5)
        R_ij_a(6) = this%R_ij_C(i, j, k, 6)
        R_ij_a(7) = this%R_ij_C(i, j, k, 7)
        R_ij_a(8) = this%R_ij_C(i, j, k, 8)
        R_ij_a(9) = this%R_ij_C(i, j, k, 9)

        arr_temp1 = matmul(reshape(this%duidxj_C(i, j, k, :),[3,3],order=[2,1]),reshape(this%duidxj_C(i, j, k, :),[3, 3],order=[2,1]))

        if (this%NNtype == 1) then
          tgrad(i, j, k, 1) = del*(this%duidxj_C(i, j, k, 1)**2+this%duidxj_C(i, j, k, 2)**2+this%duidxj_C(i, j, k, 3)**2)
          tgrad(i, j, k, 2) = del*(this%duidxj_C(i, j, k, 1)*this%duidxj_C(i, j, k, 4)+this%duidxj_C(i, j, k, 2)* &
                                        this%duidxj_C(i, j, k, 5)+this%duidxj_C(i, j, k, 3)*this%duidxj_C(i, j, k, 6))
          tgrad(i, j, k, 3) = del*(this%duidxj_C(i, j, k, 1)*this%duidxj_C(i, j, k, 7)+this%duidxj_C(i, j, k, 2)* &
                                        this%duidxj_C(i, j, k, 8)+this%duidxj_C(i, j, k, 3)*this%duidxj_C(i, j, k, 9))
          tgrad(i, j, k, 4) = tgrad(i, j, k, 2)
          tgrad(i, j, k, 5) = del*(this%duidxj_C(i, j, k, 4)**2+this%duidxj_C(i, j, k, 5)**2+this%duidxj_C(i, j, k, 6)**2)
          tgrad(i, j, k, 6) = del*(this%duidxj_C(i, j, k, 4)*this%duidxj_C(i, j, k, 7)+this%duidxj_C(i, j, k, 5)* &
                                        this%duidxj_C(i, j, k, 8)+this%duidxj_C(i, j, k, 6)*this%duidxj_C(i, j, k, 9))
          tgrad(i, j, k, 7) = tgrad(i, j, k, 3)
          tgrad(i, j, k, 8) = tgrad(i, j, k, 6)
          tgrad(i, j, k, 9) = del*(this%duidxj_C(i, j, k, 7)**2+this%duidxj_C(i, j, k, 8)**2+this%duidxj_C(i, j, k, 9)**2)
          tgradmag(1, i, j, k) = norm2(tgrad(i, j, k, :))

          magout(1, i, j, k) = norm2(this%duidxj_C(i, j, k, :))+1.0*10**(-8)
          this%invariants(1, i, k, j, 1) = 0.0
          this%invariants(1, i, k, j, 2) = 0.5*(0.0**2-(arr_temp1(1,1)+arr_temp1(2,2)+arr_temp1(3,3)))
          this%invariants(1, i, k, j, 3) = -1*(this%duidxj_C(i, j, k, 1)*this%duidxj_C(i, j, k, 5)*this%duidxj_C(i, j, k, 9)+ &
                                               this%duidxj_C(i, j, k, 2)*this%duidxj_C(i, j, k, 6)*this%duidxj_C(i, j, k, 7)+ &
                                               this%duidxj_C(i, j, k, 3)*this%duidxj_C(i, j, k, 4)*this%duidxj_C(i, j, k, 8)- &
                                               this%duidxj_C(i, j, k, 3)*this%duidxj_C(i, j, k, 5)*this%duidxj_C(i, j, k, 7)- &
                                               this%duidxj_C(i, j, k, 1)*this%duidxj_C(i, j, k, 6)*this%duidxj_C(i, j, k, 8)- &
                                               this%duidxj_C(i, j, k, 2)*this%duidxj_C(i, j, k, 4)*this%duidxj_C(i, j, k, 9))
          this%invariants(1, i, k, j, 4) = norm2(S_ij_a(:))
          this%invariants(1, i, k, j, 5) = norm2(R_ij_a(:))

          this%invariantsnorm(1, i, k, j, 1) = 0.0
          this%invariantsnorm(1, i, k, j, 2) = this%invariants(1, i, k, j, 2)/(magout(1, i, j, k)**2)
          this%invariantsnorm(1, i, k, j, 3) = this%invariants(1, i, k, j, 3)/(magout(1, i, j, k)**3)
          this%invariantsnorm(1, i, k, j, 4) = this%invariants(1, i, k, j, 4)/magout(1, i, j, k)
          this%invariantsnorm(1, i, k, j, 5) = this%invariants(1, i, k, j, 5)/magout(1, i, j, k)
          
          this%strain2(1, i, k, j, :, :) = reshape(S_ij_a(:),[3, 3],order=[2,1])
          this%rot(1, i, k, j, :, :) = reshape(R_ij_a(:),[3, 3],order=[2,1])
          this%strain2_m(1, i, k, j, :, :) = matmul(reshape(S_ij_a(:),[3, 3],order=[2,1]),reshape(S_ij_a(:),[3, 3],order=[2,1]))
          this%rot2(1, i, k, j, :, :) =  matmul(reshape(R_ij_a(:),[3, 3],order=[2,1]),reshape(R_ij_a(:),[3, 3],order=[2,1]))
          this%strainrot_m(1, i, k, j, :, :) = matmul(reshape(S_ij_a(:),[3, 3],order=[2,1]),reshape(R_ij_a(:),[3, 3],order=[2,1])) &
                         - matmul(reshape(R_ij_a(:),[3,3], order=[2,1]),reshape(S_ij_a(:),[3, 3],order=[2,1]))
          this%rotstrainrot_m(1, i, k, j, :, :) = matmul(matmul(reshape(R_ij_a(:),[3, 3], order=[2,1]), & 
                         reshape(S_ij_a(:),[3, 3],order=[2,1])),reshape(R_ij_a(:),[3, 3],order=[2,1]))
          this%strain2rot_m(1, i, k, j, :, :) = matmul(matmul(reshape(S_ij_a(:),[3, 3], order=[2,1]), &
                         reshape(S_ij_a(:),[3, 3], order=[2,1])),reshape(R_ij_a(:),[3, 3], order=[2,1])) & 
                         - matmul(reshape(R_ij_a(:),[3, 3], order=[2,1]), &
                         matmul(reshape(S_ij_a(:),[3, 3], order=[2,1]),reshape(S_ij_a(:),[3, 3], order=[2,1])))
          this%rotstrain2_m(1, i, k, j, :, :) = matmul(reshape(R_ij_a(:), [3, 3], order=[2,1]), matmul(reshape(S_ij_a(:),[3, 3], & 
                         order=[2,1]), matmul(reshape(R_ij_a(:),[3, 3], order=[2,1]), &
                         reshape(R_ij_a(:),[3, 3], order=[2,1])))) - matmul(matmul(matmul(reshape(R_ij_a(:), [3, 3], & 
                         order=[2,1]), reshape(R_ij_a(:), [3, 3], order=[2,1])), reshape(S_ij_a(:),[3, 3], & 
                         order=[2,1])),reshape(R_ij_a(:),[3, 3], order=[2,1]))
        end if

        if (this%NNtype == 2) then
          this%invariantsLSTM(1, 10, i, j, k, 1) = S_ij_a(1)+S_ij_a(5)+S_ij_a(9)
          this%invariantsLSTM(1, 10, i, j, k, 2) = arr_temp1(1,1)+arr_temp1(2,2)+arr_temp1(3,3)
          this%invariantsLSTM(1, 10, i, j, k, 3) = arr_temp2(1,1)+arr_temp2(2,2)+arr_temp2(3,3)
          this%invariantsLSTM(1, 10, i, j, k, 4) = arr_temp3(1,1)+arr_temp3(2,2)+arr_temp3(3,3)
          this%invariantsLSTM(1, 10, i, j, k, 5) = arr_temp4(1,1)+arr_temp4(2,2)+arr_temp4(3,3)
          this%invariantsLSTM(1, 10, i, j, k, 6) = arr_temp5(1,1)+arr_temp5(2,2)+arr_temp5(3,3)
          this%invariantsLSTM(1, 10, i, j, k, 7) = arr_temp6(1,1)+arr_temp6(2,2)+arr_temp6(3,3)
        end if
      end do
    end do
  end do

  if (this%NNtype == 1) then
    
  end if

  do k = 1, this%gpC%xsz(3)
     do j = 1, this%gpC%xsz(2)
       do i = 1, this%gpC%xsz(1)
          this%strain2(1, i, k, j, :, :) =  this%strain2(1, i, k, j, :, :)/magout(1, i, j, k)
          this%rot(1, i, k, j, :, :) =  this%rot(1, i, k, j, :, :)/magout(1, i, j, k)
          this%strain2_m(1, i, k, j, :, :) =  this%strain2_m(1, i, k, j, :, :)/(magout(1, i, j, k)**2)
          this%rot2(1, i, k, j, :, :) =  this%rot2(1, i, k, j, :, :)/(magout(1, i, j, k)**2)
          this%strainrot_m(1, i, k, j, :, :) =  this%strainrot_m(1, i, k, j, :, :)/(magout(1, i, j, k)**2)
          this%rotstrainrot_m(1, i, k, j, :, :) =  this%rotstrainrot_m(1, i, k, j, :, :)/(magout(1, i, j, k)**3)
          this%strain2rot_m(1, i, k, j, :, :) =  this%strain2rot_m(1, i, k, j, :, :)/(magout(1, i, j, k)**3)
          this%rotstrain2_m(1, i, k, j, :, :) =  this%rotstrain2_m(1, i, k, j, :, :)/(magout(1, i, j, k)**4)
       end do
     end do
  end do
  

  if (this%NNtype == 2) then
      if (tidNow < 512) then
        ! call Unetand onlu feed in LSTM at time step "10"
      end if
      if (tidNow >= 512) then
         timelist(1) = tidNow-511
         timelist(2) = tidNow-256
         timelist(3) = tidNow-64
         timelist(4) = tidNow-32
         timelist(5) = tidNow-16
         timelist(6) = tidNow-8
         timelist(7) = tidNow-4
         timelist(8) = tidNow-2
         timelist(9) = tidNow-1
         do n = 1,size(this%invariantsLSTM,2)-1
           call readInvariants(timelist(n),this%runID,this%datadir,this%gpC,this%invariantsLSTM(1,n,:,:,:,:),'C')
         end do
      end if
      call writeInvariants(tidNow,this%runID,this%datadir,this%gpC,this%invariantsLSTM(1,10,:,:,:,:),'C')
  end if


  if (this%NNtype == 1) then
    count = 1
    do i = 1,this%gpC%xsz(1)
      do k = 1,this%gpC%xsz(3)
        do j = 1,this%gpC%xsz(2)
          do l = 1,5
            train_loadin(count) = this%invariants(1, i, k, j, l)
            train_loadinnorm(count) = this%invariantsnorm(1, i, k, j, l)
            count = count + 1
          end do
        end do
      end do
    end do

    count = 1
    do i = 1,this%gpC%xsz(1)
      do k = 1,this%gpC%xsz(3)
        do j = 1,this%gpC%xsz(2)
          do l = 1,3
            do m = 1,3
              strain2in(count) = this%strain2(1, i, k, j, l, m)
              strain2_min(count) = this%strain2_m(1, i, k, j, l, m)
              rot2in(count) = this%rot2(1, i, k, j, l, m)
              rotin(count) = this%rot(1, i, k, j, l, m)
              strainrot_min(count) = this%strainrot_m(1, i, k, j, l, m)
              rotstrainrot_min(count) = this%rotstrainrot_m(1, i, k, j, l, m)
              strain2rot_min(count) = this%strain2rot_m(1, i, k, j, l, m)
              rotstrain2_min(count) = this%rotstrain2_m(1, i, k, j, l,  m)
              count = count + 1
            end do
          end do
        end do
      end do
    end do
        
    c_pstruc = loadpytorchgpu(train_loadinnorm, strain2in, strain2_min, rot2in, strainrot_min, &
        rotstrainrot_min, strain2rot_min, rotstrain2_min, rotin, this%gpC%xsz(1), this%gpC%xsz(3), &
        this%gpC%xsz(2), this%dx, this%dz, this%dy)
    call c_f_pointer(c_pstruc, f_pstruc, [7*this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)])
    
    count = this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)+1
    do l = 1,6
      do i = 1,this%gpC%xsz(1)
        do k = 1,this%gpC%xsz(3)
          do j = 1,this%gpC%xsz(2)
            !this%finalout(1,l,i,j,k) = f_pstruc(count)
            finalout(1,l,i,j,k) = f_pstruc(count)
            count = count + 1
          end do
        end do
      end do
    end do

    nullify(f_pstruc) 
    call c_free(c_pstruc)
    c_pstruc = c_null_ptr
  end if
 

  if (this%NNtype == 2) then
    train_loadin = reshape(this%invariantsLSTM,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*7*10])
    !deltain = reshape(this%delta,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*8])
    strain2in = reshape(this%strain2,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    strain2_min = reshape(this%strain2_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    rot2in = reshape(this%rot2,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    strainrot_min = reshape(this%strainrot_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    rotstrainrot_min = reshape(this%rotstrainrot_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    strain2rot_min = reshape(this%strain2rot_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
    rotstrain2_min = reshape(this%rotstrain2_m,[this%gpC%xsz(1)*this%gpC%xsz(2)*this%gpC%xsz(3)*9])
  end if



  do k = 1, this%gpE%xsz(3)
     do j = 1, this%gpE%xsz(2)
       do i = 1, this%gpE%xsz(1)
         S_ij_a(1) = this%S_ij_E(i, j, k, 1)
         S_ij_a(2) = this%S_ij_E(i, j, k, 2)
         S_ij_a(3) = this%S_ij_E(i, j, k, 3)
         S_ij_a(4) = this%S_ij_E(i, j, k, 2)
         S_ij_a(5) = this%S_ij_E(i, j, k, 4)
         S_ij_a(6) = this%S_ij_E(i, j, k, 5)
         S_ij_a(7) = this%S_ij_E(i, j, k, 3)
         S_ij_a(8) = this%S_ij_E(i, j, k, 5)
         S_ij_a(9) = this%S_ij_E(i, j, k, 6)
         
         R_ij_a(1) = this%R_ij_E(i, j, k, 1)
         R_ij_a(2) = this%R_ij_E(i, j, k, 2)
         R_ij_a(3) = this%R_ij_E(i, j, k, 3)
         R_ij_a(4) = this%R_ij_E(i, j, k, 4)
         R_ij_a(5) = this%R_ij_E(i, j, k, 5)
         R_ij_a(6) = this%R_ij_E(i, j, k, 6)
         R_ij_a(7) = this%R_ij_E(i, j, k, 7)
         R_ij_a(8) = this%R_ij_E(i, j, k, 8)
         R_ij_a(9) = this%R_ij_E(i, j, k, 9)

         
         arr_temp1 = matmul(reshape(this%duidxj_E(i, j, k, :),[3,3],order=[2,1]),reshape(this%duidxj_E(i, j, k, :),[3,3],order=[2,1])) 
         if (this%NNtype == 1) then
           tgradE(i, j, k, 1) = del*(this%duidxj_E(i, j, k, 1)**2+this%duidxj_E(i, j, k, 2)**2+this%duidxj_E(i, j, k, 3)**2)
           tgradE(i, j, k, 2) = del*(this%duidxj_E(i, j, k, 1)*this%duidxj_E(i, j, k, 4)+this%duidxj_E(i, j, k, 2)* &
                                        this%duidxj_E(i, j, k, 5)+this%duidxj_E(i, j, k, 3)*this%duidxj_E(i, j, k, 6))
           tgradE(i, j, k, 3) = del*(this%duidxj_E(i, j, k, 1)*this%duidxj_E(i, j, k, 7)+this%duidxj_E(i, j, k, 2)* &
                                        this%duidxj_E(i, j, k, 8)+this%duidxj_E(i, j, k, 3)*this%duidxj_E(i, j, k, 9))
           tgradE(i, j, k, 4) = tgradE(i, j, k, 2)
           tgradE(i, j, k, 5) = del*(this%duidxj_E(i, j, k, 4)**2+this%duidxj_E(i, j, k, 5)**2+this%duidxj_E(i, j, k, 6)**2)
           tgradE(i, j, k, 6) = del*(this%duidxj_E(i, j, k, 4)*this%duidxj_E(i, j, k, 7)+this%duidxj_E(i, j, k, 5)* &
                                        this%duidxj_E(i, j, k, 8)+this%duidxj_E(i, j, k, 6)*this%duidxj_E(i, j, k, 9))
           tgradE(i, j, k, 7) = tgradE(i, j, k, 3)
           tgradE(i, j, k, 8) = tgradE(i, j, k, 6)
           tgradE(i, j, k, 9) = del*(this%duidxj_E(i, j, k, 7)**2+this%duidxj_E(i, j, k, 8)**2+this%duidxj_E(i, j, k, 9)**2)
           tgradmagE(1, i, j, k) = norm2(tgradE(i, j, k, :))

           magoutE(1, i, j, k) = norm2(this%duidxj_E(i, j, k, :))+1.0*10**(-8)
           this%invariantsE(1, i, k, j, 1) = 0.0
           this%invariantsE(1, i, k, j, 2) = 0.5*(0.0**2-(arr_temp1(1,1)+arr_temp1(2,2)+arr_temp1(3,3)))
           this%invariantsE(1, i, k, j, 3) = -1*(this%duidxj_E(i, j, k, 1)*this%duidxj_E(i, j, k, 5)*this%duidxj_E(i, j, k, 9)+ &
                                               this%duidxj_E(i, j, k, 2)*this%duidxj_E(i, j, k, 6)*this%duidxj_E(i, j, k, 7)+ &
                                               this%duidxj_E(i, j, k, 3)*this%duidxj_E(i, j, k, 4)*this%duidxj_E(i, j, k, 8)- &
                                               this%duidxj_E(i, j, k, 3)*this%duidxj_E(i, j, k, 5)*this%duidxj_E(i, j, k, 7)- &
                                               this%duidxj_E(i, j, k, 1)*this%duidxj_E(i, j, k, 6)*this%duidxj_E(i, j, k, 8)- &
                                               this%duidxj_E(i, j, k, 2)*this%duidxj_E(i, j, k, 4)*this%duidxj_E(i, j, k, 9))
           this%invariantsE(1, i, k, j, 4) = norm2(S_ij_a(:))
           this%invariantsE(1, i, k, j, 5) = norm2(R_ij_a(:))

           this%invariantsEnorm(1, i, k, j, 1) = 0.0
           this%invariantsEnorm(1, i, k, j, 2) = this%invariantsE(1, i, k, j, 2)/(magoutE(1, i, j, k)**2)
           this%invariantsEnorm(1, i, k, j, 3) = this%invariantsE(1, i, k, j, 3)/(magoutE(1, i, j, k)**3)
           this%invariantsEnorm(1, i, k, j, 4) = this%invariantsE(1, i, k, j, 4)/magoutE(1, i, j, k)
           this%invariantsEnorm(1, i, k, j, 5) = this%invariantsE(1, i, k, j, 5)/magoutE(1, i, j, k)
           
           this%strain2E(1, i, k, j, :, :) = reshape(S_ij_a(:),[3,3],order=[2,1])
           this%rotE(1, i, k, j, :, :) = reshape(R_ij_a(:),[3,3],order=[2,1])
           this%strain2_mE(1, i, k, j, :, :) = matmul(reshape(S_ij_a(:),[3, 3],order=[2,1]),reshape(S_ij_a(:),[3, 3],order=[2,1]))
           this%rot2E(1, i, k, j, :, :) = matmul(reshape(R_ij_a(:),[3, 3],order=[2,1]),reshape(R_ij_a(:),[3, 3],order=[2,1]))
           this%strainrot_mE(1, i, k, j, :, :) = matmul(reshape(S_ij_a(:),[3,3],order=[2,1]),reshape(R_ij_a(:),[3,3],order=[2,1])) &
                         - matmul(reshape(R_ij_a(:),[3,3], order=[2,1]),reshape(S_ij_a(:),[3, 3],order=[2,1]))
           this%rotstrainrot_mE(1, i, k, j, :, :) = matmul(matmul(reshape(R_ij_a(:),[3, 3], order=[2,1]), &
                         reshape(S_ij_a(:),[3, 3],order=[2,1])),reshape(R_ij_a(:),[3, 3],order=[2,1]))
           this%strain2rot_mE(1, i, k, j, :, :) = matmul(matmul(reshape(S_ij_a(:),[3, 3], order=[2,1]), &
                         reshape(S_ij_a(:),[3, 3],order=[2,1])),reshape(R_ij_a(:),[3, 3], order=[2,1])) &
                         - matmul(reshape(R_ij_a(:),[3, 3],  order=[2,1]), &
                         matmul(reshape(S_ij_a(:),[3, 3], order=[2,1]),reshape(S_ij_a(:),[3, 3], order=[2,1])))
           this%rotstrain2_mE(1, i, k, j, :, :) = matmul(reshape(R_ij_a(:), [3,3], order=[2,1]), matmul(reshape(S_ij_a(:),[3, 3], &
                         order=[2,1]), matmul(reshape(R_ij_a(:),[3, 3], order=[2,1]), &
                         reshape(R_ij_a(:),[3, 3], order=[2,1])))) - matmul(matmul(matmul(reshape(R_ij_a(:), [3, 3], &
                         order=[2,1]), reshape(R_ij_a(:), [3, 3], order=[2,1])), reshape(S_ij_a(:),[3, 3], &
                         order=[2,1])),reshape(R_ij_a(:),[3, 3], order=[2,1]))
         end if

       end do
     end do
  end do


  do k = 1, this%gpE%xsz(3)
    do j = 1, this%gpE%xsz(2)
      do i = 1, this%gpE%xsz(1)
         this%strain2E(1, i, k, j, :, :) =  this%strain2E(1, i, k, j, :, :)/(magoutE(1, i, j, k))
         this%rotE(1, i, k, j, :, :) =  this%rotE(1, i, k, j, :, :)/(magoutE(1, i, j, k))
         this%strain2_mE(1, i, k, j, :, :) =  this%strain2_mE(1, i, k, j, :, :)/(magoutE(1, i, j, k)**2)
         this%rot2E(1, i, k, j, :, :) =  this%rot2E(1, i, k, j, :, :)/(magoutE(1, i, j, k)**2)
         this%strainrot_mE(1, i, k, j, :, :) =  this%strainrot_mE(1, i, k, j, :, :)/(magoutE(1, i, j, k)**2)
         this%rotstrainrot_mE(1, i, k, j, :, :) =  this%rotstrainrot_mE(1, i, k, j, :, :)/(magoutE(1, i, j, k)**3)
         this%strain2rot_mE(1, i, k, j, :, :) =  this%strain2rot_mE(1, i, k, j, :, :)/(magoutE(1, i, j, k)**3)
         this%rotstrain2_mE(1, i, k, j, :, :) =  this%rotstrain2_mE(1, i, k, j, :, :)/(magoutE(1, i, j, k)**4)
         
      end do
    end do
  end do

  if (this%NNtype == 1) then
    count = 1
    do i = 1,this%gpE%xsz(1)
      do k = 1,this%gpE%xsz(3)
        do j = 1,this%gpE%xsz(2)
          do l = 1,5
            train_loadinE(count) = this%invariantsE(1, i, k, j, l)
            train_loadinEnorm(count) = this%invariantsEnorm(1, i, k, j, l)
            count = count + 1
          end do
        end do
      end do
    end do
    
    count = 1
    do i = 1,this%gpE%xsz(1)
      do k = 1,this%gpE%xsz(3)
        do j = 1,this%gpE%xsz(2)
          do l = 1,3
            do m = 1,3
              strain2inE(count) = this%strain2E(1, i, k, j, l, m)
              rotinE(count) = this%rotE(1, i, k, j, l, m)
              strain2_minE(count) = this%strain2_mE(1, i, k, j, l, m)
              rot2inE(count) = this%rot2E(1, i, k, j, l, m)
              strainrot_minE(count) = this%strainrot_mE(1, i, k, j, l, m)
              rotstrainrot_minE(count) = this%rotstrainrot_mE(1, i, k, j, l, m)
              strain2rot_minE(count) = this%strain2rot_mE(1, i, k, j, l, m)
              rotstrain2_minE(count) = this%rotstrain2_mE(1, i, k, j, l, m)
              count = count + 1
            end do
          end do
        end do
      end do
    end do
    
    c_pEstruc = loadpytorchgpu(train_loadinEnorm, strain2inE, strain2_minE, rot2inE, strainrot_minE, &
       rotstrainrot_minE, strain2rot_minE, rotstrain2_minE, rotinE, this%gpE%xsz(1), this%gpE%xsz(3), &
       this%gpE%xsz(2), this%dx, this%dz, this%dy)
    call c_f_pointer(c_pEstruc, f_pEstruc, [7*this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)])
          
   
    count = this%gpE%xsz(1)*this%gpE%xsz(2)*this%gpE%xsz(3)+1
    do l = 1,6
      do i = 1,this%gpE%xsz(1)
        do k = 1,this%gpE%xsz(3)
          do j = 1,this%gpE%xsz(2)
            !this%finaloutE(1,l,i,j,k) = f_pEstruc(count)
            finaloutE(1,l,i,j,k) = f_pEstruc(count)
            count = count + 1
          end do
         end do
      end do
    end do

    nullify(f_pEstruc)
    call c_free(c_pEstruc)
    c_pEstruc = c_null_ptr
       
  end if

  trace = finalout(1,1,:,:,:)*tgradmag(1,:,:,:)+finalout(1,4,:,:,:)*tgradmag(1,:,:,:)+finalout(1,6,:,:,:)*tgradmag(1,:,:,:)
  this%tau_11 = finalout(1,1,:,:,:)*tgradmag(1,:,:,:)-1.0/3.0*(trace)
  this%tau_12 = finalout(1,2,:,:,:)*tgradmag(1,:,:,:)
  this%tau_13 = finaloutE(1,3,:,:,:)*tgradmagE(1,:,:,:)
  this%tau_22 = finalout(1,4,:,:,:)*tgradmag(1,:,:,:)-1.0/3.0*(trace)
  this%tau_23 = finaloutE(1,5,:,:,:)*tgradmagE(1,:,:,:)
  this%tau_33 = finalout(1,6,:,:,:)*tgradmag(1,:,:,:)-1.0/3.0*(trace)


end subroutine
