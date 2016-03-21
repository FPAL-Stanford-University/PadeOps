program test_projection
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero
    use reductions, only: p_maxval
    use timer, only: tic, toc
    use basic_io, only: read_2d_ascii 
    use decomp_2d
    use spectralMod, only: spectral
    use cd06staggstuff, only: cd06stagg 
    use PadePoissonMod, only: padepoisson
    use exits, only: message
    use cf90stuff, only: cf90
    use basic_io, only: write_binary

    implicit none
   
    real(rkind), dimension(:,:,:), allocatable :: x, y, z, u, v, w, wE
    real(rkind), dimension(:,:,:), allocatable :: tmp
    type(decomp_info) :: gpC, gpE
    integer :: nx = 100, ny = 100, nz = 60
    integer :: ierr, prow = 0, pcol = 0
    real(rkind) :: dx, dy, dz
    integer :: dimTransform = 2
    character(len=clen) :: filename = "/home/aditya90/Codes/PadeOps/data/OpenFoam_AllData.txt" 
    real(rkind), dimension(:,:), allocatable :: temp 
    type(spectral), allocatable :: spect
    type(spectral), allocatable :: spectE
    type(padepoisson), allocatable :: poiss 
    real(rkind), dimension(:,:,:), allocatable :: divergence, rtmp1z, rtmp2z, rtmp1x, rtmp1y, rtmp1yE, rtmp1zE
    complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what
    type(cd06stagg), allocatable :: der
    complex(rkind), dimension(:,:,:), allocatable :: ctmp1, ctmp2, ctmp3, ctmp4
    type(cf90) :: filzC, filzE

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol)

    call get_decomp_info(gpC)
    call decomp_info_init(nx, ny, nz+1, gpE)

    call read_2d_ascii(temp,filename)

    allocate(tmp(nx,ny,nz))
    allocate( x ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( y ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( z ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( u ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( v ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( w ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( divergence ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    
   

    ! Set x values
    tmp = reshape(temp(:,1),[nx,ny,nz])
    x = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))

    ! Set y values
    tmp = reshape(temp(:,2),[nx,ny,nz])
    y = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))
     
    ! Set z values
    tmp = reshape(temp(:,3),[nx,ny,nz])
    z = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))
    
    ! Set u values
    tmp = reshape(temp(:,4),[nx,ny,nz])
    u = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))
    
    ! Set v values
    tmp = reshape(temp(:,5),[nx,ny,nz])
    v = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))
    
    ! Set w values
    tmp = reshape(temp(:,6),[nx,ny,nz])
    w = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))
   
    dx = x(2,1,1) - x(1,1,1)
    dy = y(1,2,1) - y(1,1,1)
    dz = z(1,1,2) - z(1,1,1)

    deallocate(temp,tmp)

    !! Allocate the operator, spectral and poisson classes
    allocate(spect)
    allocate(spectE)
    allocate(poiss)
    allocate(der)

    call der%init(nz, dz, .false., .false.) 
    call spect%init("x", nx, ny, nz, dx, dy, dz, "four", "2/3rd", dimTransform,.false.)
    call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", "2/3rd", dimTransform,.false.,.false.,.false.,.false.)
    call poiss%init(dx, dy, dz, spect, spectE, der)
    ierr = filzC%init(nz,.false.)
    print*, ierr
    ierr = filzE%init(nz+1,.false.)
    print*, ierr

    allocate(ctmp1(spect%spectdecomp%zsz(1),spect%spectdecomp%zsz(2),spect%spectdecomp%zsz(3)))
    allocate(ctmp2(spect%spectdecomp%zsz(1),spect%spectdecomp%zsz(2),spect%spectdecomp%zsz(3)))
    allocate(ctmp3(spectE%spectdecomp%zsz(1),spectE%spectdecomp%zsz(2),spectE%spectdecomp%zsz(3)))
    allocate(ctmp4(spectE%spectdecomp%zsz(1),spectE%spectdecomp%zsz(2),spectE%spectdecomp%zsz(3)))

    !! Generate wE
    allocate(wE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    allocate(rtmp1y(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
    allocate(rtmp1z(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))
    allocate(rtmp1zE(gpE%zsz(1),gpE%zsz(2),gpE%zsz(3)))
    allocate(rtmp1yE(gpE%ysz(1),gpE%ysz(2),gpE%ysz(3)))

    ! Step 1: take w from x -> z
    call transpose_x_to_y(w,rtmp1y,gpC)
    call transpose_y_to_z(rtmp1y,rtmp1z,gpC)

    ! Step 2: Interpolate onto rtmp1zE
    call der%InterpZ_C2E(rtmp1z,rtmp1zE,size(rtmp1z,1),size(rtmp1z,2))

    ! Step 3: Transpose back from z -> x
    call transpose_z_to_y(rtmp1zE,rtmp1yE,gpE)
    call transpose_y_to_x(rtmp1yE,wE,gpE)

    
    deallocate(rtmp1y,rtmp1z,rtmp1zE,rtmp1yE)


    !! Allocate storage for ffts
    call spect%alloc_r2c_out(uhat)
    call spect%alloc_r2c_out(vhat)
    call spectE%alloc_r2c_out(what)
    
    !! Take FFTs
    call spect%fft(u,uhat)
    call spect%fft(v,vhat)
    call spectE%fft(wE,what)

    !! Set Oddballs to zero
    uhat(:,ny/2+1,:) = zero + imi*zero
    vhat(:,ny/2+1,:) = zero + imi*zero
    what(:,ny/2+1,:) = zero + imi*zero
    call spectE%ifft(what, wE, .true.) 
    call spect%ifft(uhat, u, .true.) 
    call spect%ifft(vhat, v, .true.) 
    call spect%fft(u,uhat)
    call spect%fft(v,vhat)
    call spectE%fft(wE,what)

    call write_binary(u,"u_preF.dat")
    call write_binary(v,"v_preF.dat")
    call write_binary(wE,"w_preF.dat")

    !! Filter in z
    call tic()
    call transpose_y_to_z(uhat,ctmp1,spect%spectdecomp)
    call filzC%filter3(ctmp1,ctmp2,size(ctmp1,1),size(ctmp1,2))
    call transpose_z_to_y(ctmp2,uhat,spect%spectdecomp)
    call transpose_y_to_z(vhat,ctmp1,spect%spectdecomp)
    call filzC%filter3(ctmp1,ctmp2,size(ctmp1,1),size(ctmp1,2))
    call transpose_z_to_y(ctmp2,vhat,spect%spectdecomp)
    call transpose_y_to_z(what,ctmp3,spectE%spectdecomp)
    call filzE%filter3(ctmp3,ctmp4,size(ctmp3,1),size(ctmp3,2))
    call transpose_z_to_y(ctmp4,what,spectE%spectdecomp)
    call toc()
    
    
    call spect%ifft(uhat,u)
    call spect%ifft(vhat,v)
    call spect%ifft(what,wE)

    call write_binary(u,"u_postF.dat")
    call write_binary(v,"v_postF.dat")
    call write_binary(wE,"w_postF.dat")

    !! Divergence Check 
    call poiss%DivergenceCheck(uhat,vhat,what,divergence) 
    call message(1,"Max Divergence before projection:", p_maxval(divergence))

    call tic()
    !! Call Poisson Projection
    call poiss%PressureProjection(uhat,vhat,what) 
    call toc()   


    !! Divergence Check 
    call poiss%DivergenceCheck(uhat,vhat,what,divergence) 
    call message(1,"Max Divergence after projection:", p_maxval(divergence))

    !! IFFT back to x 
    call spect%ifft(uhat,u)
    call spect%ifft(vhat,v)
    call spect%ifft(what,wE)
    
    call write_binary(u,"u_postP.dat")
    call write_binary(v,"v_postP.dat")
    call write_binary(wE,"w_postP.dat")

    !! Deallocate storage
    deallocate(uhat, vhat, what)

    deallocate(wE,w,u,v, divergence)

    !! Destroy classes
    call poiss%destroy()
    call der%destroy()
    call spect%destroy()
    call spectE%destroy()
    deallocate(spect)
    deallocate(der)
    deallocate(poiss)


    call MPI_Finalize(ierr)

end program 
