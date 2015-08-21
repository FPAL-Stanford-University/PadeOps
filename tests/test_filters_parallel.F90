module transfer_funcs

    use kind_parameters, only: rkind
    use constants,       only: two
    implicit none

contains
    
    function GetTransferFunctionCF90(k) result(T)
        use cf90stuff, only: alpha90, beta90, a90, b90, c90, d90, e90
        real(rkind), intent(in) :: k
        real(rkind) :: T

        T = (a90 + two* b90*COS(k) + two*c90*COS(two*k) +two*d90*COS(3._rkind*k) + two*e90*COS(4._rkind*k) ) &
          / (1._rkind + two*alpha90*COS(k) + two*beta90*COS(two*k) )

    end function
    
    function GetTransferFunctionGaussian(k) result(T)
        use gaussianstuff, only: agf, bgf, cgf, dgf, egf
        real(rkind), intent(in) :: k
        real(rkind) :: T

        T = (agf + two* bgf*COS(k) + two*cgf*COS(two*k) +two*dgf*COS(3._rkind*k) + two*egf*COS(4._rkind*k) )

    end function

    function GetTransferFunctionLstsq(k) result(T)
        use lstsqstuff, only: als, bls, cls, dls, els
        real(rkind), intent(in) :: k
        real(rkind) :: T

        T = (als + two* bls*COS(k) + two*cls*COS(two*k) +two*dls*COS(3._rkind*k) + two*els*COS(4._rkind*k) )

    end function

end module


program test_filters

    use mpi
    use decomp_2d, only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use kind_parameters, only: rkind
    use constants,       only: two,pi
    use exits,           only: message
    use reductions,      only: P_MAXVAL
    use timer,           only: tic,toc
    use FiltersMod,      only: filters
    use transfer_funcs,  only: GetTransferFunctionCF90,GetTransferFunctionGaussian,GetTransferFunctionLstsq
    implicit none

    integer :: nx = 256, ny=256, nz=256


    type( decomp_info ) :: gp
    type( filters ) :: method1, method2, method3
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df
    real(rkind), dimension(:,:,:), allocatable :: xtmp1,xtmp2,ztmp1,ztmp2
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 64._rkind

    real(rkind) :: k_norm_x,k_norm_y,k_norm_z,TF1_x,TF1_y,TF1_z,TF2_x,TF2_y,TF2_z,TF3_x,TF3_y,TF3_z

    integer :: prow = 0, pcol = 0
    integer :: i,j,k,ierr

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)
    
    allocate( x     ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( y     ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( z     ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( f     ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( df    ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )

    allocate( xtmp1 ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( xtmp2 ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    
    allocate( ztmp1 ( gp%zsz(1), gp%zsz(2), gp%zsz(3) ) )
    allocate( ztmp2 ( gp%zsz(1), gp%zsz(2), gp%zsz(3) ) )
    
    dx = two*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    dz = two*pi/real(nz,rkind)

    do k=1,gp%ysz(3)
        do j=1,gp%ysz(2)
            do i=1,gp%ysz(1)
                x(i,j,k) = real( gp%yst(1) - 1 + i - 1, rkind )*dx
                y(i,j,k) = real( gp%yst(2) - 1 + j - 1, rkind )*dy
                z(i,j,k) = real( gp%yst(3) - 1 + k - 1, rkind )*dz
                f(i,j,k) = cos(omega * x(i,j,k)) * cos(omega * y(i,j,k)) * cos(omega * z(i,j,k))
            end do
        end do 
    end do 

    k_norm_x = omega * dx
    TF1_x = GetTransferFunctionCF90    (k_norm_x)
    TF2_x = GetTransferFunctionGaussian(k_norm_x)
    TF3_x = GetTransferFunctionLstsq   (k_norm_x)
    
    k_norm_y = omega * dy
    TF1_y = GetTransferFunctionCF90    (k_norm_y)
    TF2_y = GetTransferFunctionGaussian(k_norm_y)
    TF3_y = GetTransferFunctionLstsq   (k_norm_y)
    
    k_norm_z = omega * dz
    TF1_z = GetTransferFunctionCF90    (k_norm_z)
    TF2_z = GetTransferFunctionGaussian(k_norm_z)
    TF3_z = GetTransferFunctionLstsq   (k_norm_z)

    call message("Created initial data")

    call method1%init(                          gp, &
                            .TRUE., .TRUE., .TRUE., &
                            "cf90", "cf90", "cf90" )

    call method2%init(                           gp, &
                            .TRUE., .TRUE.,  .TRUE., &
                            "gaussian", "gaussian", "gaussian" )
    
    call method3%init(                           gp, &
                            .TRUE., .TRUE.,  .TRUE., &
                            "lstsq", "lstsq", "lstsq" )

    call message("Initialized all methods")
   
   
    call message("===========================================")
    call message("Now trying METHOD 1: CF90"                  )
    call message("===========================================")
    call tic()
    call transpose_y_to_x(f,xtmp1,gp)
    call method1 % filterx(xtmp1,xtmp2)
    call transpose_x_to_y(xtmp2,df,gp)
    call toc ("Time to get the x filter:")
    call message("Maximum error", P_MAXVAL(df - TF1_x*f))

    call tic() 
    call method1 % filtery(f,df)
    call toc ("Time to get the y filter:")
    call message("Maximum error", P_MAXVAL(df - TF1_y*f))

    call tic()
    call transpose_y_to_z(f,ztmp1,gp)
    call method1 % filterz(ztmp1,ztmp2)
    call transpose_z_to_y(ztmp2,df,gp)
    call toc ("Time to get the z filter:")
    call message("Maximum error", P_MAXVAL(df - TF1_z*f))


    call message("===========================================")
    call message("Now trying METHOD 2: GAUSSIAN"              )
    call message("===========================================")
    call tic() 
    call transpose_y_to_x(f,xtmp1,gp)
    call method2 % filterx(xtmp1,xtmp2)
    call transpose_x_to_y(xtmp2,df,gp)
    call toc ("Time to get the x filter:")
    call message("Maximum error", P_MAXVAL(df - TF2_x*f))

    call tic() 
    call method2 % filtery(f,df)
    call toc ("Time to get the y filter:")
    call message("Maximum error", P_MAXVAL(df - TF2_y*f))

    call tic() 
    call transpose_y_to_z(f,ztmp1,gp)
    call method2 % filterz(ztmp1,ztmp2)
    call transpose_z_to_y(ztmp2,df,gp)
    call toc ("Time to get the z filter:")
    call message("Maximum error", P_MAXVAL(df - TF2_z*f))


    call message("===========================================")
    call message("Now trying METHOD 3: LSTSQ"                 )
    call message("===========================================")
    call tic() 
    call transpose_y_to_x(f,xtmp1,gp)
    call method3 % filterx(xtmp1,xtmp2)
    call transpose_x_to_y(xtmp2,df,gp)
    call toc ("Time to get the x filter:")
    call message("Maximum error", P_MAXVAL(df - TF3_x*f))

    call tic() 
    call method3 % filtery(f,df)
    call toc ("Time to get the y filter:")
    call message("Maximum error", P_MAXVAL(df - TF3_y*f))

    call tic() 
    call transpose_y_to_z(f,ztmp1,gp)
    call method3 % filterz(ztmp1,ztmp2)
    call transpose_z_to_y(ztmp2,df,gp)
    call toc ("Time to get the z filter:")
    call message("Maximum error", P_MAXVAL(df - TF3_z*f))

    call message("==========================================")
    call message("Destroying everything"                     )
    call message("==========================================")

    call method1%destroy
    call method2%destroy
    call method3%destroy
    deallocate( x )
    deallocate( y )
    deallocate( z )
    deallocate( f )
    deallocate( df )

    deallocate( xtmp1 )
    deallocate( xtmp2 )
    
    deallocate( ztmp1 )
    deallocate( ztmp2 )
    
    call decomp_2d_finalize
    call MPI_Finalize(ierr)
    

end program
