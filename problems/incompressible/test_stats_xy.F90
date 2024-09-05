! Template for PadeOps

#include "test_stats_xy_files/initialize.F90"       

program test_stats_xy
    use mpi
    use kind_parameters,    only: clen, rkind
    use IncompressibleGrid, only: igrid
    use stats_xy_mod,       only: stats_xy
    use constants,          only: pi
    use exits,              only: message
    use fortran_assert,     only: assert
    use decomp_2d,          only: transpose_x_to_y, transpose_y_to_z, nrank
    use test_stats_xy_parameters, only: coefs, wavenums, taucoefs  
    use reductions,         only: p_maxval
    use basic_io,           only: read_1d_ascii
    implicit none

    type(igrid), allocatable, target :: SM, SM1
    type(stats_xy) :: stats
    character(len=clen) :: inputfile1, inputfile, tempname, stats_info_dir
    character(len=clen) :: stmt
    integer :: ierr, n, ioUnit, stid, i
    real(rkind), dimension(:,:,:), allocatable :: tau11, tau12, tau13, tau22, tau23, tau33, tau13C, tau23C
    real(rkind), dimension(:,:,:), pointer :: x, y, z, xE, yE, zE
    real(rkind), dimension(:,:,:), allocatable :: stats_cpy
    real(rkind), dimension(:,:,:,:), allocatable :: stats_sca_cpy
    real(rkind), dimension(:), allocatable :: analyticalX, analyticalY, analyticalZ, analyticalT, analyticalWT, z1d
    integer, dimension(3) :: zids, steps
    real(rkind) :: omega
    integer :: argc

    ! Required for reading the namelist, but not used directly in the main program
    real(rkind) :: Lx, Ly, Lz, zmin, Tref
    logical :: symmetricDomain = .true.
    integer :: num_stats_instances = 1
    
    call MPI_Init(ierr)
    namelist /SMinput/ Lx, Ly, Lz, symmetricDomain, zmin, Tref, stats_info_dir, num_stats_instances

    argc = command_argument_count()
    if (argc > 1) then
      call GETARG(1,inputfile1)                                            
      call GETARG(2,stmt) 

      ioUnit = 11
      open(unit=ioUnit, file=trim(inputfile1), form='FORMATTED')
      read(unit=ioUnit, NML=SMinput)
      close(ioUnit)    

      !===================================================================================================
      !==================================================================================================
      ! Read in fields and compute stats and compare to the stats computed in MATLAB

      allocate(SM1)
      call SM1%init(inputfile1,.true.)
      SM1%step = 4

      write(tempname,'(A,A4)')trim(stats_info_dir)//"/STATS_test_interscale.inp"
      call stats%init(trim(tempname),SM1)
      call stats%compute_stats()

      ! Allocate memory
      allocate(stats_cpy(SM1%gpC%zsz(3),102,2))
      allocate(stats_sca_cpy(SM1%gpC%zsz(3),47,1,2))

      ! Copy stats
      call stats%set_tidx(1)
      call stats%copy_stats(stats_cpy,stats_sca_cpy)
      call message(1,"Stats copied")
      call message(" ")

      !   Turbulent transport:
      call message("Turbulent transport")
      ! Read in the MATLAB data
      call read_1d_ascii(analyticalX ,trim(SM1%inputdir)//'/turb_trans_X_scl01.dat')
      call read_1d_ascii(analyticalY ,trim(SM1%inputdir)//'/turb_trans_Y_scl01.dat')
      call read_1d_ascii(analyticalZ ,trim(SM1%inputdir)//'/turb_trans_Z_scl01.dat')
      call read_1d_ascii(analyticalT ,trim(SM1%inputdir)//'/turb_trans_T_scl01.dat')
      call read_1d_ascii(analyticalwT,trim(SM1%inputdir)//'/turb_trans_wT_scl01.dat')
      call message("        > max absolute difference R11: ",maxval(abs(stats_cpy(    :,41  ,1) - analyticalX )))
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(    :,56  ,1) - analyticalY )))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(    :,71  ,1) - analyticalZ )))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,4 ,1,1) - analyticalT )))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,24,1,1) - analyticalWT)))
      call message("    > Scale2:")
      call read_1d_ascii(analyticalX ,trim(SM1%inputdir)//'/turb_trans_X_scl02.dat')
      call read_1d_ascii(analyticalY ,trim(SM1%inputdir)//'/turb_trans_Y_scl02.dat')
      call read_1d_ascii(analyticalZ ,trim(SM1%inputdir)//'/turb_trans_Z_scl02.dat')
      call read_1d_ascii(analyticalT ,trim(SM1%inputdir)//'/turb_trans_T_scl02.dat')
      call read_1d_ascii(analyticalwT,trim(SM1%inputdir)//'/turb_trans_wT_scl02.dat')
      call message("        > max absolute difference R11: ",maxval(abs(stats_cpy(    :,41  ,2) - analyticalX )))
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(    :,56  ,2) - analyticalY )))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(    :,71  ,2) - analyticalZ )))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,4 ,1,2) - analyticalT )))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,24,1,2) - analyticalWT)))
      call message("*** *** ***")

      call message("Turbulent-transport mixed-scale")
      call read_1d_ascii(analyticalX ,trim(SM1%inputdir)//'/turb_trans_mixed_X_scl01.dat')
      call read_1d_ascii(analyticalY ,trim(SM1%inputdir)//'/turb_trans_mixed_Y_scl01.dat')
      call read_1d_ascii(analyticalZ ,trim(SM1%inputdir)//'/turb_trans_mixed_Z_scl01.dat')
      call read_1d_ascii(analyticalT ,trim(SM1%inputdir)//'/turb_trans_mixed_T_scl01.dat')
      call read_1d_ascii(analyticalwT,trim(SM1%inputdir)//'/turb_trans_mixed_wT_scl01.dat')
      call message("        > max absolute difference R11: ",maxval(abs(stats_cpy(    :,51  ,1) - analyticalX )))
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(    :,66  ,1) - analyticalY )))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(    :,81  ,1) - analyticalZ )))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,11,1,1) - analyticalT )))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,37,1,1) - analyticalWT)))
      call message("    > Scale2:")
      call read_1d_ascii(analyticalX ,trim(SM1%inputdir)//'/turb_trans_mixed_X_scl02.dat')
      call read_1d_ascii(analyticalY ,trim(SM1%inputdir)//'/turb_trans_mixed_Y_scl02.dat')
      call read_1d_ascii(analyticalZ ,trim(SM1%inputdir)//'/turb_trans_mixed_Z_scl02.dat')
      call read_1d_ascii(analyticalT ,trim(SM1%inputdir)//'/turb_trans_mixed_T_scl02.dat')
      call read_1d_ascii(analyticalwT,trim(SM1%inputdir)//'/turb_trans_mixed_wT_scl02.dat')
      call message("        > max absolute difference R11: ",maxval(abs(stats_cpy(    :,51  ,2) - analyticalX )))
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(    :,66  ,2) - analyticalY )))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(    :,81  ,2) - analyticalZ )))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,11,1,2) - analyticalT )))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,37,1,2) - analyticalWT)))
      call message("*** *** ***")

      ! Interscale transfer
      call message("Interscale transfer")
      call read_1d_ascii(analyticalX ,trim(SM1%inputdir)//'/interscale_X.dat')
      call read_1d_ascii(analyticalY ,trim(SM1%inputdir)//'/interscale_Y.dat')
      call read_1d_ascii(analyticalZ ,trim(SM1%inputdir)//'/interscale_Z.dat')
      call read_1d_ascii(analyticalT ,trim(SM1%inputdir)//'/interscale_T.dat')
      call read_1d_ascii(analyticalwT,trim(SM1%inputdir)//'/interscale_wT.dat')
      call message("        > max absolute difference R11: ",maxval(abs(stats_cpy(    :,49  ,1) + analyticalX )))
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(    :,64  ,1) + analyticalY )))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(    :,79  ,1) + analyticalZ )))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,9 ,1,1) + analyticalT )))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,35,1,1) + analyticalWT)))
      call message("    > Scale2:")
      call message("        > max absolute difference R11: ",maxval(abs(stats_cpy(    :,49  ,2) - analyticalX )))
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(    :,64  ,2) - analyticalY )))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(    :,79  ,2) - analyticalZ )))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,9 ,1,2) - analyticalT )))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,35,1,2) - analyticalWT)))
      call message("*** *** ***")

      deallocate(stats_cpy,stats_sca_cpy,analyticalX,analyticalY,analyticalZ,analyticalT,analyticalwT)
      call stats%destroy()
      call SM1%destroy()
      deallocate(SM1)

      !===================================================================================================
      !==================================================================================================
    else
      call getarg(1,inputfile)

      ioUnit = 11
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
      read(unit=ioUnit, NML=SMinput)
      close(ioUnit)    

      allocate(SM)                                                     
      
      call SM%init(inputfile, .true.)                     
      
      ! Set the time step to trigger the stats calculation
      SM%step = 4

      ! Initialize stats class instances
      write(tempname,'(A,A4)')trim(stats_info_dir)//"/STATS_test.inp"
      call stats%init(trim(tempname),SM)

      call message(" ")
      call message("==========================================================")
      call message(0, "stats and igrid initialized! Now running the test.")

      ! Allocate memory
      allocate(tau11(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
      allocate(tau12(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
      allocate(tau13C(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
      allocate(tau13(SM%gpE%xsz(1),SM%gpE%xsz(2),SM%gpE%xsz(3)))
      allocate(tau22(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
      allocate(tau23C(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
      allocate(tau23(SM%gpE%xsz(1),SM%gpE%xsz(2),SM%gpE%xsz(3)))
      allocate(tau33(SM%gpC%xsz(1),SM%gpC%xsz(2),SM%gpC%xsz(3)))
      allocate(stats_cpy(SM%gpC%zsz(3),102,2))
      allocate(stats_sca_cpy(SM%gpC%zsz(3),47,1,2))
      allocate(analyticalX(SM%gpC%zsz(3)))
      allocate(analyticalY(SM%gpC%zsz(3)))
      allocate(analyticalZ(SM%gpC%zsz(3)))
      allocate(analyticalT(SM%gpC%zsz(3)))
      allocate(analyticalWT(SM%gpC%zsz(3)))
      allocate(z1d(SM%gpC%zsz(3)))
      call message(1,"memory allocated")
      
      ! Mesh
      x  => SM%mesh(:,:,:,1)
      y  => SM%mesh(:,:,:,2)
      z  => SM%mesh(:,:,:,3)
      xE => SM%meshE(:,:,:,1)
      yE => SM%meshE(:,:,:,2)
      zE => SM%meshE(:,:,:,3)
      call transpose_x_to_y(z,SM%rbuffyC(:,:,:,1),SM%gpC)
      call transpose_y_to_z(SM%rbuffyC(:,:,:,1),SM%rbuffzC(:,:,:,1),SM%gpC)
      z1d = SM%rbuffzC(1,1,:,1)
      call message(1,"Finished setting up the mesh")

      ! Redefine Lz ("Lz" used in the input file is really zmax)
      Lz = Lz - zmin

      ! Set the value for the fields
      SM%u  = 3.0d0 + coefs(1,1)*z  + &
        coefs(1,2)*cos(wavenums(1)*x )*cos(wavenums(1)*y )*cos(wavenums(1)*(z -zmin-Lz/2)) + &
        coefs(1,3)*cos(wavenums(2)*x )*cos(wavenums(2)*y )*cos(wavenums(2)*(z -zmin-Lz/2))
      SM%v  = 2.5d0 + coefs(2,1)*z  + &
        coefs(2,2)*sin(wavenums(1)*x )*cos(wavenums(1)*y )*cos(wavenums(1)*(z -zmin-Lz/2)) + &
        coefs(2,3)*sin(wavenums(2)*x )*cos(wavenums(2)*y )*cos(wavenums(2)*(z -zmin-Lz/2))
      SM%wC = 2.6d0 + coefs(3,1)*z  + &
        coefs(3,2)*sin(wavenums(1)*x )*cos(wavenums(1)*y )*sin(wavenums(1)*(z -zmin-Lz/2)) + &
        coefs(3,3)*sin(wavenums(2)*x )*cos(wavenums(2)*y )*sin(wavenums(2)*(z -zmin-Lz/2))
      SM%w  = 2.6d0 + coefs(3,1)*zE + &
        coefs(3,2)*sin(wavenums(1)*xE)*cos(wavenums(1)*yE)*sin(wavenums(1)*(zE-zmin-Lz/2)) + &
        coefs(3,3)*sin(wavenums(2)*xE)*cos(wavenums(2)*yE)*sin(wavenums(2)*(zE-zmin-Lz/2))
      SM%pressure = 1.1d0 + coefs(4,1)*z + &
        coefs(4,2)*sin(wavenums(1)*x )*cos(wavenums(1)*y )*sin(wavenums(1)*(z -zmin-Lz/2)) + &
        coefs(4,3)*sin(wavenums(2)*x )*cos(wavenums(2)*y )*sin(wavenums(2)*(z -zmin-Lz/2))
      SM%T = 3.4d0 + coefs(5,1)*z + &
         coefs(5,2)*sin(wavenums(1)*x )*cos(wavenums(1)*y )*sin(wavenums(1)*(z -zmin-Lz/2)) + &
         coefs(5,3)*sin(wavenums(2)*x )*cos(wavenums(2)*y )*sin(wavenums(2)*(z -zmin-Lz/2)) 
      ! dudx
      SM%duidxjC(:,:,:,1) =            -wavenums(1)*coefs(1,2)*sin(wavenums(1)*x )*cos(wavenums(1)*y )*cos(wavenums(1)*(z -zmin-Lz/2)) + &
                                       -wavenums(2)*coefs(1,3)*sin(wavenums(2)*x )*cos(wavenums(2)*y )*cos(wavenums(2)*(z -zmin-Lz/2)) 
      ! dudy
      SM%duidxjC(:,:,:,2) =            -wavenums(1)*coefs(1,2)*cos(wavenums(1)*x )*sin(wavenums(1)*y )*cos(wavenums(1)*(z -zmin-Lz/2)) + &
                                       -wavenums(2)*coefs(1,3)*cos(wavenums(2)*x )*sin(wavenums(2)*y )*cos(wavenums(2)*(z -zmin-Lz/2))
      ! dudz
      SM%duidxjC(:,:,:,3) = coefs(1,1) -wavenums(1)*coefs(1,2)*cos(wavenums(1)*x )*cos(wavenums(1)*y )*sin(wavenums(1)*(z -zmin-Lz/2)) + &
                                       -wavenums(2)*coefs(1,3)*cos(wavenums(2)*x )*cos(wavenums(2)*y )*sin(wavenums(2)*(z -zmin-Lz/2)) 
      ! dvdx
      SM%duidxjC(:,:,:,4) =             wavenums(1)*coefs(2,2)*cos(wavenums(1)*x )*cos(wavenums(1)*y )*cos(wavenums(1)*(z -zmin-Lz/2)) + &
                                        wavenums(2)*coefs(2,3)*cos(wavenums(2)*x )*cos(wavenums(2)*y )*cos(wavenums(2)*(z -zmin-Lz/2)) 
      ! dvdy
      SM%duidxjC(:,:,:,5) =            -wavenums(1)*coefs(2,2)*sin(wavenums(1)*x )*sin(wavenums(1)*y )*cos(wavenums(1)*(z -zmin-Lz/2)) + &
                                       -wavenums(2)*coefs(2,3)*sin(wavenums(2)*x )*sin(wavenums(2)*y )*cos(wavenums(2)*(z -zmin-Lz/2)) 
      ! dvdz
      SM%duidxjC(:,:,:,6) = coefs(2,1) -wavenums(1)*coefs(2,2)*sin(wavenums(1)*x )*cos(wavenums(1)*y )*sin(wavenums(1)*(z -zmin-Lz/2)) + &
                                       -wavenums(2)*coefs(2,3)*sin(wavenums(2)*x )*cos(wavenums(2)*y )*sin(wavenums(2)*(z -zmin-Lz/2))
      ! dwdx
      SM%duidxjC(:,:,:,7) =             wavenums(1)*coefs(3,2)*cos(wavenums(1)*x )*cos(wavenums(1)*y )*sin(wavenums(1)*(z -zmin-Lz/2)) + &
                                        wavenums(2)*coefs(3,3)*cos(wavenums(2)*x )*cos(wavenums(2)*y )*sin(wavenums(2)*(z -zmin-Lz/2)) 
      ! dwdy
      SM%duidxjC(:,:,:,8) =            -wavenums(1)*coefs(3,2)*sin(wavenums(1)*x )*sin(wavenums(1)*y )*sin(wavenums(1)*(z -zmin-Lz/2)) + &
                                       -wavenums(2)*coefs(3,3)*sin(wavenums(2)*x )*sin(wavenums(2)*y )*sin(wavenums(2)*(z -zmin-Lz/2))
      ! dwdz
      SM%duidxjC(:,:,:,9) = coefs(3,1) +wavenums(1)*coefs(3,2)*sin(wavenums(1)*x )*cos(wavenums(1)*y )*cos(wavenums(1)*(z -zmin-Lz/2)) + & 
                                        wavenums(2)*coefs(3,3)*sin(wavenums(2)*x )*cos(wavenums(2)*y )*cos(wavenums(2)*(z -zmin-Lz/2))

      ! SGS fields
      tau11  = tauCoefs(1)*sin(wavenums(3)*x )*cos(wavenums(3)*y )*cos(wavenums(3)*(z -zmin-Lz/2))
      tau12  = tauCoefs(2)*sin(wavenums(3)*x )*cos(wavenums(3)*y )*cos(wavenums(3)*(z -zmin-Lz/2))
      tau13C = tauCoefs(3)*sin(wavenums(3)*x )*cos(wavenums(3)*y )*cos(wavenums(3)*(z -zmin-Lz/2))
      tau13  = tauCoefs(3)*sin(wavenums(3)*xE)*cos(wavenums(3)*yE)*cos(wavenums(3)*(zE-zmin-Lz/2))
      tau22  = tauCoefs(4)*sin(wavenums(3)*x )*cos(wavenums(3)*y )*cos(wavenums(3)*(z -zmin-Lz/2))
      tau23  = tauCoefs(5)*sin(wavenums(3)*xE)*cos(wavenums(3)*yE)*cos(wavenums(3)*(zE-zmin-Lz/2))
      tau23C = tauCoefs(5)*sin(wavenums(3)*x )*cos(wavenums(3)*y )*cos(wavenums(3)*(z -zmin-Lz/2))
      tau33  = tauCoefs(6)*sin(wavenums(3)*x )*cos(wavenums(3)*y )*cos(wavenums(3)*(z -zmin-Lz/2))
      call SM%sgsmodel%set_tauij(tau11,tau12,tau13,tau22,tau23,tau33,tau13C,tau23C)

      SM%q1_T  = tauCoefs(7)*sin(wavenums(3)*x )*cos(wavenums(3)*y )*cos(wavenums(3)*(z -zmin-Lz/2)) 
      SM%q2_T  = tauCoefs(8)*sin(wavenums(3)*x )*cos(wavenums(3)*y )*cos(wavenums(3)*(z -zmin-Lz/2))
      SM%q3_T  = tauCoefs(9)*sin(wavenums(3)*xE)*cos(wavenums(3)*yE)*cos(wavenums(3)*(zE-zmin-Lz/2))
      SM%q3_TC = tauCoefs(9)*sin(wavenums(3)*x )*cos(wavenums(3)*y )*cos(wavenums(3)*(z -zmin-Lz/2))

      ! dTdx
      SM%dTdxC =               wavenums(1)*coefs(5,2)*cos(wavenums(1)*x )*cos(wavenums(1)*y )*sin(wavenums(1)*(z -zmin-Lz/2)) + &
                               wavenums(2)*coefs(5,3)*cos(wavenums(2)*x )*cos(wavenums(2)*y )*sin(wavenums(2)*(z -zmin-Lz/2)) 
      ! dTdy
      SM%dTdyC =              -wavenums(1)*coefs(5,2)*sin(wavenums(1)*x )*sin(wavenums(1)*y )*sin(wavenums(1)*(z -zmin-Lz/2)) + &
                              -wavenums(2)*coefs(5,3)*sin(wavenums(2)*x )*sin(wavenums(2)*y )*sin(wavenums(2)*(z -zmin-Lz/2)) 
      ! dTdz
      SM%dTdzC = coefs(5,1) +  wavenums(1)*coefs(5,2)*sin(wavenums(1)*x )*cos(wavenums(1)*y )*cos(wavenums(1)*(z -zmin-Lz/2)) + &
                               wavenums(2)*coefs(5,3)*sin(wavenums(2)*x )*cos(wavenums(2)*y )*cos(wavenums(2)*(z -zmin-Lz/2)) 

      ! Force components
      SM%spectForceLayer%fx = & 
        coefs(6,2)*cos(wavenums(1)*x )*cos(wavenums(1)*y )*cos(wavenums(1)*(z -zmin-Lz/2)) + &
        coefs(6,3)*cos(wavenums(2)*x )*cos(wavenums(2)*y )*cos(wavenums(2)*(z -zmin-Lz/2))
      SM%spectForceLayer%fy = &
        coefs(7,2)*sin(wavenums(1)*x )*cos(wavenums(1)*y )*cos(wavenums(1)*(z -zmin-Lz/2)) + &
        coefs(7,3)*sin(wavenums(2)*x )*cos(wavenums(2)*y )*cos(wavenums(2)*(z -zmin-Lz/2))
      SM%spectForceLayer%fz = &
        coefs(8,2)*sin(wavenums(1)*xE)*cos(wavenums(1)*yE)*sin(wavenums(1)*(zE-zmin-Lz/2)) + &
        coefs(8,3)*sin(wavenums(2)*xE)*cos(wavenums(2)*yE)*sin(wavenums(2)*(zE-zmin-Lz/2))
      SM%spectForceLayer%ampFact_x = 1.d0
      SM%spectForceLayer%ampFact_y = 1.d0
      SM%spectForceLayer%ampFact_z = 1.d0

      call message(1,"Finished defining the analytical fields")

      ! Compute the budgets
      call stats%compute_stats()
      call message(1,"Stats computed")

      ! Copy the data
      call stats%copy_stats(stats_cpy,stats_sca_cpy)
      call message(1,"Stats copied")
      call message(" ")

      ! TODO: Confirm the computed stats terms match the analytical answers
      ! Mean velocity:
      analyticalX = 3.0d0 + coefs(1,1)*z1d
      analyticalY = 2.5d0 + coefs(2,1)*z1d
      analyticalZ = 2.6d0 + coefs(3,1)*z1d
      analyticalT = 0.5d0*(analyticalX**2.d0 + analyticalY**2.d0 + analyticalZ**2.d0)
      write(stmt,'(F)')maxval(abs(stats_cpy(    :,36  ,1) - analyticalT))
      call assert(maxval(abs(stats_cpy(    :,36  ,1) - analyticalT)) < 1.d-14,'MKE,r -- '//trim(stmt))
      call assert(maxval(abs(stats_cpy(    :,19  ,1) - analyticalX)) < 1.d-14,'meanU')
      call assert(maxval(abs(stats_cpy(    :,20  ,1) - analyticalY)) < 1.d-14,'meanV')
      call assert(maxval(abs(stats_cpy(    :,21  ,1) - analyticalZ)) < 1.d-14,'meanW')
      call assert(maxval(abs(stats_cpy(    :,36  ,2) - analyticalT)) < 1.d-14,'MKE,s')
      call assert(maxval(abs(stats_cpy(    :,19  ,2) - analyticalX)) < 1.d-14,'meanU')
      call assert(maxval(abs(stats_cpy(    :,20  ,2) - analyticalY)) < 1.d-14,'meanV')
      call assert(maxval(abs(stats_cpy(    :,21  ,2) - analyticalZ)) < 1.d-14,'meanW')
      analyticalT = 3.4d0 + coefs(5,1)*z1d 
      call assert(maxval(abs(stats_sca_cpy(:,13,1,1) - analyticalT)) < 1.d-14,'meanT')
      call assert(maxval(abs(stats_sca_cpy(:,13,1,2) - analyticalT)) < 1.d-14,'meanT')
      analyticalT = coefs(5,1) 
      call assert(maxval(abs(stats_sca_cpy(:,12,1,1) - analyticalT)) < 1.d-14,'meandTdz')
      call assert(maxval(abs(stats_sca_cpy(:,12,1,2) - analyticalT)) < 1.d-14,'meandTdz')
      call message("Mean fields test PASSED!")

      ! Velocity variance
      analyticalX = 0.25d0*(coefs(1,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 
      analyticalY = 0.25d0*(coefs(2,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 
      analyticalZ = 0.25d0*(coefs(3,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 
      analyticalT = 0.5d0*(analyticalX + analyticalY + analyticalZ)
      call assert(maxval(abs(stats_cpy(    :,17  ,1) - analyticalT)) < 1.d-13,'tke,r')
      call assert(maxval(abs(stats_cpy(    :,23  ,1) - analyticalX)) < 1.d-13,'uu')
      call assert(maxval(abs(stats_cpy(    :,24  ,1) - analyticalY)) < 1.d-13,'vv')
      call assert(maxval(abs(stats_cpy(    :,25  ,1) - analyticalZ)) < 1.d-13,'ww')
      call message("  > Velocity variance test for ur PASSED!")
      analyticalX = 0.25d0*(coefs(1,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 
      analyticalY = 0.25d0*(coefs(2,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 
      analyticalZ = 0.25d0*(coefs(3,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 
      analyticalT = 0.5d0*(analyticalX + analyticalY + analyticalZ)
      call assert(maxval(abs(stats_cpy(    :,17  ,2) - analyticalT)) < 1.d-13,'tke,s')
      call assert(maxval(abs(stats_cpy(    :,23  ,2) - analyticalX)) < 1.d-13,'uu')
      call assert(maxval(abs(stats_cpy(    :,24  ,2) - analyticalY)) < 1.d-13,'vv')
      call assert(maxval(abs(stats_cpy(    :,25  ,2) - analyticalZ)) < 1.d-13,'ww')
      call message("  > Velocity variance test for us PASSED!")
      analyticalT = 0.25d0*(coefs(5,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0
      call assert(maxval(abs(stats_sca_cpy(:,17,1,1) - analyticalT)) < 1.d-13,'TT')
      analyticalT = 0.25d0*(coefs(5,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0
      call assert(maxval(abs(stats_sca_cpy(:,17,1,2) - analyticalT)) < 1.d-13,'TT')
      call message("Velocity variance test PASSED!")

      ! Velocity covariance
      analyticalX = 0.d0 
      analyticalY = 0.d0 
      analyticalZ = 0.25d0*coefs(2,2)*coefs(3,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0))*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)) 
      call assert(maxval(abs(stats_cpy(:,26,1) - analyticalX)) < 1.d-14,'uv')
      call assert(maxval(abs(stats_cpy(:,27,1) - analyticalY)) < 1.d-14,'uw')
      call assert(maxval(abs(stats_cpy(:,28,1) - analyticalZ)) < 1.d-14,'vw')
      analyticalX = 0.d0 
      analyticalY = 0.d0 
      analyticalZ = 0.25d0*coefs(2,3)*coefs(3,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0))*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)) 
      call assert(maxval(abs(stats_cpy(:,26,2) - analyticalX)) < 1.d-14,'uv')
      call assert(maxval(abs(stats_cpy(:,27,2) - analyticalY)) < 1.d-14,'uw')
      call assert(maxval(abs(stats_cpy(:,28,2) - analyticalZ)) < 1.d-14,'vw')
      call message("Velocity covariance test PASSED!")

      ! Pressure stuff
      analyticalT = (0.5d0*coefs(4,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0
      analyticalX = 0.d0
      analyticalY = 0.25d0*coefs(2,2)*coefs(4,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0))*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)) 
      analyticalZ = 0.25d0*coefs(3,2)*coefs(4,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0))**2.d0
      call assert(maxval(abs(stats_cpy(:,29,1) - analyticalT)) < 1.d-14,'<pp>,r')
      call assert(maxval(abs(stats_cpy(:,30,1) - analyticalX)) < 1.d-14,'<up>,r')
      call assert(maxval(abs(stats_cpy(:,31,1) - analyticalY)) < 1.d-14,'<vp>,r')
      call assert(maxval(abs(stats_cpy(:,32,1) - analyticalZ)) < 1.d-14,'<wp>,r')
      analyticalT = (0.5d0*coefs(4,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0
      analyticalX = 0.d0
      analyticalY = 0.25d0*coefs(2,3)*coefs(4,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0))*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)) 
      analyticalZ = 0.25d0*coefs(3,3)*coefs(4,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0))**2.d0
      call assert(maxval(abs(stats_cpy(:,29,2) - analyticalT)) < 1.d-14,'<pp>,s')
      call assert(maxval(abs(stats_cpy(:,30,2) - analyticalX)) < 1.d-14,'<up>,s')
      call assert(maxval(abs(stats_cpy(:,31,2) - analyticalY)) < 1.d-14,'<vp>,s')
      call assert(maxval(abs(stats_cpy(:,32,2) - analyticalZ)) < 1.d-14,'<wp>,s')
      analyticalT = 1.1d0 + coefs(4,1)*z1d
      call assert(maxval(abs(stats_cpy(:,22,1) - analyticalT)) < 1.d-14,'meanP,r') 
      call assert(maxval(abs(stats_cpy(:,22,2) - analyticalT)) < 1.d-14,'meanP,s') 


      ! Force stuff
      call message("Force energy (fluctuations)")
      analyticalX = (0.5d0*coefs(6,2)*cos(wavenums(1)*(z1d-zmin-Lz/2)))**2.d0 + &
                    (0.5d0*coefs(7,2)*cos(wavenums(1)*(z1d-zmin-Lz/2)))**2.d0 + & 
                    (0.5d0*coefs(8,2)*sin(wavenums(1)*(z1d-zmin-Lz/2)))**2.d0 
      analyticalY = (0.5d0*coefs(6,3)*cos(wavenums(2)*(z1d-zmin-Lz/2)))**2.d0 + &
                    (0.5d0*coefs(7,3)*cos(wavenums(2)*(z1d-zmin-Lz/2)))**2.d0 + & 
                    (0.5d0*coefs(8,3)*sin(wavenums(2)*(z1d-zmin-Lz/2)))**2.d0 
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call message("        > max absolute difference for <fifi>,r: ",p_maxval(maxval(abs(stats_cpy(:,18,1) - analyticalX))))
      call message("        > max absolute difference for <fifi>,s: ",p_maxval(maxval(abs(stats_cpy(:,18,2) - analyticalY))))
      !call assert(maxval(abs(stats_cpy(:,18,1) - analyticalX)) < 1.d-12,'<fifi>,r')
      !call assert(maxval(abs(stats_cpy(:,18,2) - analyticalY)) < 1.d-12,'<fifi>,s -- '//trim(stmt))
      call assert(maxval(abs(stats_cpy(:,33,1))) < 1.d-12,'meanFx,r')
      call assert(maxval(abs(stats_cpy(:,34,1))) < 1.d-12,'meanFy,r')
      call assert(maxval(abs(stats_cpy(:,35,1))) < 1.d-12,'meanFz,r')
      call assert(maxval(abs(stats_cpy(:,33,2))) < 1.d-12,'meanFx,s')
      call assert(maxval(abs(stats_cpy(:,34,2))) < 1.d-12,'meanFy,s')
      call assert(maxval(abs(stats_cpy(:,35,2))) < 1.d-12,'meanFz,s')
      call message("*** *** ***")


      ! Temperature flux
      analyticalX = 0.d0
      analyticalY = 0.25d0*coefs(2,2)*coefs(5,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0))*sin(wavenums(1)*(z1d-zmin-Lz/2.d0))
      analyticalZ = 0.25d0*coefs(3,2)*coefs(5,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0))**2.d0
      write(stmt,'(E)')p_maxval(maxval(abs(stats_sca_cpy(:,15,1,1) - analyticalZ)))
      call assert(maxval(abs(stats_sca_cpy(:,14,1,1) - analyticalX)) < 1.d-14,'uT,r')
      call assert(maxval(abs(stats_sca_cpy(:,15,1,1) - analyticalY)) < 1.d-14,'vT,r -- '//trim(stmt))
      call assert(maxval(abs(stats_sca_cpy(:,16,1,1) - analyticalZ)) < 1.d-14,'wT,r')
      analyticalY = 0.25d0*coefs(2,3)*coefs(5,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0))*sin(wavenums(2)*(z1d-zmin-Lz/2.d0))
      analyticalZ = 0.25d0*coefs(3,3)*coefs(5,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0))**2.d0
      write(stmt,'(E)')p_maxval(maxval(abs(stats_sca_cpy(:,15,1,2) - analyticalZ)))
      call assert(maxval(abs(stats_sca_cpy(:,14,1,2) - analyticalX)) < 1.d-14,'uT,s')
      call assert(maxval(abs(stats_sca_cpy(:,15,1,2) - analyticalY)) < 1.d-14,'vT,s -- '//trim(stmt))
      call assert(maxval(abs(stats_sca_cpy(:,16,1,2) - analyticalZ)) < 1.d-14,'wT,s')

      ! Scalar variance fluxes
      write(stmt,'(E)')p_maxval(maxval(abs(stats_sca_cpy(:,46,1,1))))
      call assert(maxval(abs(stats_sca_cpy(:,45,1,1))) < 1.d-14,'<uTT>r')
      call assert(maxval(abs(stats_sca_cpy(:,46,1,1))) < 1.d-14,'<vTT>r -- '//trim(stmt))
      call assert(maxval(abs(stats_sca_cpy(:,47,1,1))) < 1.d-14,'<wTT>r')
      call assert(maxval(abs(stats_sca_cpy(:,45,1,2))) < 1.d-14,'<uTT>s')
      call assert(maxval(abs(stats_sca_cpy(:,46,1,2))) < 1.d-14,'<vTT>s')
      call assert(maxval(abs(stats_sca_cpy(:,47,1,2))) < 1.d-14,'<wTT>s')

      ! Kinetic energy flux
      call assert(maxval(abs(stats_cpy(:,100,1))) < 1.d-14,'<u*uiui>/2,r')
      call assert(maxval(abs(stats_cpy(:,101,1))) < 1.d-14,'<v*uiui>/2,r')
      call assert(maxval(abs(stats_cpy(:,102,1))) < 1.d-14,'<w*uiui>/2,r')
      call assert(maxval(abs(stats_cpy(:,100,2))) < 1.d-14,'<u*uiui>/2,r')
      call assert(maxval(abs(stats_cpy(:,101,2))) < 1.d-14,'<v*uiui>/2,r')
      call assert(maxval(abs(stats_cpy(:,102,2))) < 1.d-14,'<w*uiui>/2,r')

      ! Velocity gradient variances
      analyticalX = 0.25d0*(wavenums(1)*coefs(1,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dudx^2>
      analyticalY = analyticalX !<dudy^2>
      analyticalZ = 0.25d0*(wavenums(1)*coefs(1,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dudz^2>
      call assert(maxval(abs(stats_cpy(:,82,1) - analyticalX)) < 1.d-14,'Var(dudx),r')
      call assert(maxval(abs(stats_cpy(:,83,1) - analyticalY)) < 1.d-14,'Var(dudy),r')
      call assert(maxval(abs(stats_cpy(:,84,1) - analyticalZ)) < 1.d-14,'Var(dudz),r')
      analyticalX = 0.25d0*(wavenums(2)*coefs(1,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dudx^2>
      analyticalY = analyticalX !<dudy^2>
      analyticalZ = 0.25d0*(wavenums(2)*coefs(1,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dudz^2>
      call assert(maxval(abs(stats_cpy(:,82,2) - analyticalX)) < 1.d-14,'Var(dudx),s')
      call assert(maxval(abs(stats_cpy(:,83,2) - analyticalY)) < 1.d-14,'Var(dudy),s')
      call assert(maxval(abs(stats_cpy(:,84,2) - analyticalZ)) < 1.d-14,'Var(dudz),s')
      
      analyticalX = 0.25d0*(wavenums(1)*coefs(2,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dvdx^2>
      analyticalY = analyticalX !<dvdy^2>
      analyticalZ = 0.25d0*(wavenums(1)*coefs(2,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dvdz^2>
      call assert(maxval(abs(stats_cpy(:,85,1) - analyticalX)) < 1.d-14,'Var(dvdx),r')
      call assert(maxval(abs(stats_cpy(:,86,1) - analyticalY)) < 1.d-14,'Var(dvdy),r')
      call assert(maxval(abs(stats_cpy(:,87,1) - analyticalZ)) < 1.d-14,'Var(dvdz),r')
      analyticalX = 0.25d0*(wavenums(2)*coefs(2,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dvdx^2>
      analyticalY = analyticalX !<dvdy^2>
      analyticalZ = 0.25d0*(wavenums(2)*coefs(2,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dvdz^2>
      call assert(maxval(abs(stats_cpy(:,85,2) - analyticalX)) < 1.d-14,'Var(dvdx),s')
      call assert(maxval(abs(stats_cpy(:,86,2) - analyticalY)) < 1.d-14,'Var(dvdy),s')
      call assert(maxval(abs(stats_cpy(:,87,2) - analyticalZ)) < 1.d-14,'Var(dvdz),s')
      
      analyticalX = 0.25d0*(wavenums(1)*coefs(3,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dwdx^2>
      analyticalY = analyticalX !<dwdy^2>
      analyticalZ = 0.25d0*(wavenums(1)*coefs(3,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dwdz^2>
      call assert(maxval(abs(stats_cpy(:,88,1) - analyticalX)) < 1.d-14,'Var(dwdx),r')
      call assert(maxval(abs(stats_cpy(:,89,1) - analyticalY)) < 1.d-14,'Var(dwdy),r')
      call assert(maxval(abs(stats_cpy(:,90,1) - analyticalZ)) < 1.d-14,'Var(dwdz),r')
      analyticalX = 0.25d0*(wavenums(2)*coefs(3,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dwdx^2>
      analyticalY = analyticalX !<dwdy^2>
      analyticalZ = 0.25d0*(wavenums(2)*coefs(3,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 !<dwdz^2>
      call assert(maxval(abs(stats_cpy(:,88,2) - analyticalX)) < 1.d-14,'Var(dwdx),s')
      call assert(maxval(abs(stats_cpy(:,89,2) - analyticalY)) < 1.d-14,'Var(dwdy),s')
      call assert(maxval(abs(stats_cpy(:,90,2) - analyticalZ)) < 1.d-14,'Var(dwdz),s')
      
      ! Temperature gradient variances
      analyticalX = (0.5d0*coefs(5,2)*wavenums(1)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 
      analyticalY = (0.5d0*coefs(5,2)*wavenums(1)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 
      analyticalZ = (0.5d0*coefs(5,2)*wavenums(1)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0 
      call assert(maxval(abs(stats_sca_cpy(:,38,1,1) - analyticalX)) < 1.d-14,'Var(dTdx),r')
      call assert(maxval(abs(stats_sca_cpy(:,39,1,1) - analyticalY)) < 1.d-14,'Var(dTdy),r')
      call assert(maxval(abs(stats_sca_cpy(:,40,1,1) - analyticalZ)) < 1.d-14,'Var(dTdz),r')
      analyticalX = (0.5d0*coefs(5,3)*wavenums(2)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 
      analyticalY = (0.5d0*coefs(5,3)*wavenums(2)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 
      analyticalZ = (0.5d0*coefs(5,3)*wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0 
      call assert(maxval(abs(stats_sca_cpy(:,38,1,2) - analyticalX)) < 1.d-14,'Var(dTdx),s')
      call assert(maxval(abs(stats_sca_cpy(:,39,1,2) - analyticalY)) < 1.d-14,'Var(dTdy),s')
      call assert(maxval(abs(stats_sca_cpy(:,40,1,2) - analyticalZ)) < 1.d-14,'Var(dTdz),s')

      ! SGS tensor variances
      call message("SGS tensor variance")
      call assert(maxval(abs(stats_cpy(:,94,1))) < 1.d-14,'Var(tau11),r')
      call assert(maxval(abs(stats_cpy(:,95,1))) < 1.d-14,'Var(tau12),r')
      call assert(maxval(abs(stats_cpy(:,96,1))) < 1.d-14,'Var(tau13),r')
      call assert(maxval(abs(stats_cpy(:,97,1))) < 1.d-14,'Var(tau22),r')
      call assert(maxval(abs(stats_cpy(:,98,1))) < 1.d-14,'Var(tau23),r')
      call assert(maxval(abs(stats_cpy(:,99,1))) < 1.d-14,'Var(tau33),r')
      analyticalX = (0.5d0*taucoefs(1)*cos(wavenums(3)*(z1d-zmin-Lz/2.d0)))**2.d0
      analyticalY = (0.5d0*taucoefs(2)*cos(wavenums(3)*(z1d-zmin-Lz/2.d0)))**2.d0
      analyticalZ = (0.5d0*taucoefs(3)*cos(wavenums(3)*(z1d-zmin-Lz/2.d0)))**2.d0
      call assert(maxval(abs(stats_cpy(:,94,2) - analyticalX)) < 1.d-14,'Var(tau11),s')
      call assert(maxval(abs(stats_cpy(:,95,2) - analyticalY)) < 1.d-14,'Var(tau12),s')
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call message("        > max absolute difference Var(tau13s): ",p_maxval(maxval(abs(stats_cpy(:,96,2) - analyticalZ))))
      analyticalX = (0.5d0*taucoefs(4)*cos(wavenums(3)*(z1d-zmin-Lz/2.d0)))**2.d0
      analyticalY = (0.5d0*taucoefs(5)*cos(wavenums(3)*(z1d-zmin-Lz/2.d0)))**2.d0
      analyticalZ = (0.5d0*taucoefs(6)*cos(wavenums(3)*(z1d-zmin-Lz/2.d0)))**2.d0
      call assert(maxval(abs(stats_cpy(:,97,2) - analyticalX)) < 1.d-14,'Var(tau22),s')
      call message("        > max absolute difference Var(tau23s): ",p_maxval(maxval(abs(stats_cpy(:,98,2) - analyticalY))))
      call assert(maxval(abs(stats_cpy(:,99,2) - analyticalZ)) < 1.d-14,'Var(tau33),s')
      call message("*** *** ***")

      ! SGS heat flux moments
      call message("SGS heat flux moments")
      call assert(maxval(abs(stats_sca_cpy(:,41,1,1))) < 1.d-14,'Var(q1),r')
      call assert(maxval(abs(stats_sca_cpy(:,42,1,1))) < 1.d-14,'Var(q2),r')
      call assert(maxval(abs(stats_sca_cpy(:,43,1,1))) < 1.d-14,'Var(q3),r')
      call assert(maxval(abs(stats_sca_cpy(:,44,1,1))) < 1.d-14,'Mean(q3),r')
      call assert(maxval(abs(stats_sca_cpy(:,44,1,2))) < 1.d-14,'Mean(q3),s')
      analyticalX = (0.5d0*taucoefs(7)*cos(wavenums(3)*(z1d-zmin-Lz/2.d0)))**2.d0
      analyticalY = (0.5d0*taucoefs(8)*cos(wavenums(3)*(z1d-zmin-Lz/2.d0)))**2.d0
      analyticalZ = (0.5d0*taucoefs(9)*cos(wavenums(3)*(z1d-zmin-Lz/2.d0)))**2.d0
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call assert(maxval(abs(stats_sca_cpy(:,41,1,2) - analyticalX)) < 1.d-14,'Var(q1),s')
      call assert(maxval(abs(stats_sca_cpy(:,42,1,2) - analyticalY)) < 1.d-14,'Var(q2),s')
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call message("        > max absolute difference Var(q3s): ",p_maxval(&
        maxval(abs(stats_sca_cpy(:,43,1,2) - analyticalZ))))
      call message("*** *** ***")

      ! Strain rate variances (S12, S13, S23 only)
      analyticalX = (0.25d0*wavenums(1)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0*(coefs(1,2)**2.d0 + coefs(2,2)**2.d0) !<S12^2>
      analyticalY = (0.25d0*wavenums(1)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0*(coefs(1,2)**2.d0 + coefs(3,2)**2.d0 - 2.d0*coefs(1,2)*coefs(3,2))
      analyticalZ = (0.25d0*wavenums(1)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0*(coefs(2,2)**2.d0 + coefs(3,2)**2.d0) !<S23^2>
      write(stmt,'(E)')maxval(abs(stats_cpy(:,92,1) - analyticalY))
      call assert(maxval(abs(stats_cpy(:,91,1) - analyticalX)) < 1.d-14,'Var(S12),r')
      call assert(maxval(abs(stats_cpy(:,92,1) - analyticalY)) < 1.d-14,'Var(S13),r -- '//trim(stmt))
      call assert(maxval(abs(stats_cpy(:,93,1) - analyticalZ)) < 1.d-14,'Var(S23),r')
      analyticalX = (0.25d0*wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0*(coefs(1,3)**2.d0 + coefs(2,3)**2.d0) !<S12^2>
      analyticalY = (0.25d0*wavenums(2)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0*(coefs(1,3)**2.d0 + coefs(3,3)**2.d0 - 2.d0*coefs(1,3)*coefs(3,3))
      analyticalZ = (0.25d0*wavenums(2)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0*(coefs(2,3)**2.d0 + coefs(3,3)**2.d0) !<S23^2>
      call assert(maxval(abs(stats_cpy(:,91,2) - analyticalX)) < 1.d-14,'Var(S12),s')
      call assert(maxval(abs(stats_cpy(:,92,2) - analyticalY)) < 1.d-14,'Var(S13),s')
      call assert(maxval(abs(stats_cpy(:,93,2) - analyticalZ)) < 1.d-14,'Var(S23),s')


      ! TKE and Rii budget terms
      !   Production:
      call message("Shear production")
      analyticalX = 0.d0 ! R11 production
      analyticalY = -0.5d0*coefs(2,2)*coefs(3,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0))*sin(wavenums(1)*(z1d-zmin-Lz/2.d0))*coefs(2,1)
      analyticalZ = -0.5d0*(coefs(3,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0*coefs(3,1)
      analyticalWT = -0.25d0*coefs(3,2)*coefs(5,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0))**2.d0*coefs(3,1)! <wT> budget
      call assert(maxval(abs(stats_cpy(    :, 2  ,1) - (analyticalX+analyticalY+analyticalZ)/2.d0)) < 1.d-14,'TKE shear production')
      call assert(maxval(abs(stats_cpy(    :,38  ,1) - analyticalX)) < 1.d-14,'R11 shear production')
      call assert(maxval(abs(stats_cpy(    :,53  ,1) - analyticalY)) < 1.d-14,'R22 shear production')
      call assert(maxval(abs(stats_cpy(    :,68  ,1) - analyticalZ)) < 1.d-14,'R33 shear production')
      call assert(maxval(abs(stats_sca_cpy(:,20,1,1) - analyticalWT)) < 1.d-14,'wT,r shear production')
      call message("    R11, R22, R33, and TKE shear production for 'r' fields test PASSED!")
      analyticalX = 0.d0 ! R11 production
      analyticalY = -0.5d0*coefs(2,3)*coefs(3,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0))*sin(wavenums(2)*(z1d-zmin-Lz/2.d0))*coefs(2,1)
      analyticalZ = -0.5d0*(coefs(3,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0*coefs(3,1)
      analyticalWT = -0.25d0*coefs(3,3)*coefs(5,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0))**2.d0*coefs(3,1)! <wT> budget
      call assert(maxval(abs(stats_cpy(    :, 2  ,2) - (analyticalX+analyticalY+analyticalZ)/2.d0)) < 1.d-14,'TKE shear production')
      call assert(maxval(abs(stats_cpy(    :,38  ,2) - analyticalX)) < 1.d-14,'R11 shear production')
      call assert(maxval(abs(stats_cpy(    :,53  ,2) - analyticalY)) < 1.d-14,'R22 shear production')
      call assert(maxval(abs(stats_cpy(    :,68  ,2) - analyticalZ)) < 1.d-14,'R33 shear production')
      call assert(maxval(abs(stats_sca_cpy(:,20,1,2) - analyticalWT)) < 1.d-14,'wT,s shear production')
      call message("    R11, R22, R33, and TKE shear production for 's' fields test PASSED!")

      call message("Gradient production")
      analyticalT  = -2.d0*0.25d0*coefs(3,2)*coefs(5,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0))**2.d0*coefs(5,1) ! <TT> budget
      analyticalWT = -0.25d0*(coefs(3,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)))**2.d0*coefs(5,1) ! <wT> budget
      call assert(maxval(abs(stats_sca_cpy(:,2 ,1,1) - analyticalT)) < 1.d-12,'TT,r gradient production')
      call assert(maxval(abs(stats_sca_cpy(:,19,1,1) - analyticalWT)) < 1.d-12,'wT,r gradient production')
      analyticalT  = -2.d0*0.25d0*coefs(3,3)*coefs(5,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0))**2.d0*coefs(5,1) ! <TT> budget
      analyticalWT = -0.25d0*(coefs(3,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)))**2.d0*coefs(5,1) ! <wT> budget
      call assert(maxval(abs(stats_sca_cpy(:,2 ,1,2) - analyticalT)) < 1.d-12,'TT,s gradient production')
      call assert(maxval(abs(stats_sca_cpy(:,19,1,2) - analyticalWT)) < 1.d-12,'wT,s gradient production')


      !   TODO: Force production
      !   Convective transport
      call message("Convective transport")
      analyticalX  = -(2.6d0 + coefs(3,1)*z1d) * &
        0.5d0*coefs(1,2)**2.d0*cos(wavenums(1)*(z1d-zmin-Lz/2.d0)) * &
        wavenums(1)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0))
      analyticalY  = -(2.6d0 + coefs(3,1)*z1d) * &
        0.5d0*coefs(2,2)**2.d0*cos(wavenums(1)*(z1d-zmin-Lz/2.d0)) * &
        wavenums(1)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0))
      analyticalZ  = (2.6d0 + coefs(3,1)*z1d) * &
        0.5d0*coefs(3,2)**2.d0*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)) * &
        wavenums(1)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0))
      analyticalT  = (2.6d0 + coefs(3,1)*z1d) * &
        0.5d0*coefs(5,2)**2.d0*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)) * &
        wavenums(1)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0))
      analyticalWT = (2.6d0 + coefs(3,1)*z1d) * &
        0.5d0*coefs(3,2)*coefs(5,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.d0)) * &
        wavenums(1)*cos(wavenums(1)*(z1d-zmin-Lz/2.d0))
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call message("        > max absolute difference TKE: ",p_maxval(maxval(abs(stats_cpy(    :,4   ,1) - 0.5d0*(analyticalX+analyticalY+analyticalZ)))))
      call message("        > max absolute difference R11: ",p_maxval(maxval(abs(stats_cpy(    :,40  ,1) - analyticalX ))))
      call message("        > max absolute difference R22: ",p_maxval(maxval(abs(stats_cpy(    :,55  ,1) - analyticalY ))))
      call message("        > max absolute difference R33: ",p_maxval(maxval(abs(stats_cpy(    :,70  ,1) - analyticalZ ))))
      call message("        > max absolute difference TT:  ",p_maxval(maxval(abs(stats_sca_cpy(:,3 ,1,1) - analyticalT ))))
      call message("        > max absolute difference wT1:  ",p_maxval(maxval(abs(stats_sca_cpy(:,23,1,1) - analyticalWT))))
      call message("    > Scale 2:")
      analyticalX = -(2.6d0 + coefs(3,1)*z1d) * &
        0.5d0*coefs(1,3)**2.d0*cos(wavenums(2)*(z1d-zmin-Lz/2.d0)) * &
        wavenums(2)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0))
      analyticalY = -(2.6d0 + coefs(3,1)*z1d) * &
        0.5d0*coefs(2,3)**2.d0*cos(wavenums(2)*(z1d-zmin-Lz/2.d0)) * &
        wavenums(2)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0))
      analyticalZ =  (2.6d0 + coefs(3,1)*z1d) * &
        0.5d0*coefs(3,3)**2.d0*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)) * &
        wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0))
      analyticalT  = (2.6d0 + coefs(3,1)*z1d) * &
        0.5d0*coefs(5,3)**2.d0*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)) * &
        wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0))
      analyticalWT = (2.6d0 + coefs(3,1)*z1d) * &
        0.5d0*coefs(3,3)*coefs(5,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.d0)) * &
        wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.d0))
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call message("        > max absolute difference TKE: ",p_maxval(maxval(abs(stats_cpy(    :,4   ,2) - 0.5d0*(analyticalX+analyticalY+analyticalZ)))))
      call message("        > max absolute difference R11: ",p_maxval(maxval(abs(stats_cpy(    :,40  ,2) - analyticalX ))))
      call message("        > max absolute difference R22: ",p_maxval(maxval(abs(stats_cpy(    :,55  ,2) - analyticalY ))))
      call message("        > max absolute difference R33: ",p_maxval(maxval(abs(stats_cpy(    :,70  ,2) - analyticalZ ))))
      call message("        > max absolute difference TT:  ",p_maxval(maxval(abs(stats_sca_cpy(:,3 ,1,2) - analyticalT ))))
      call message("        > max absolute difference wT1:  ",p_maxval(maxval(abs(stats_sca_cpy(:,23,1,2) - analyticalWT))))
      call message('    Convective transport test finished. Verify pass/fail separately.')
      call message("*** *** ***")

      !    Pressure transport
      call message("Pressure transport")
      analyticalX = 0.d0 
      analyticalY = 0.d0 
      analyticalZ = -2.d0*0.5d0*coefs(3,2)*coefs(4,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.0))*wavenums(1)*cos(wavenums(1)*(z1d-zmin-Lz/2.0))
      analyticalwT = -0.5d0*coefs(4,2)*coefs(5,2)*wavenums(1)*sin(wavenums(1)*(z1d-zmin-Lz/2.0))*cos(wavenums(1)*(z1d-zmin-Lz/2.0))
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call message("        > max absolute difference TKE: ",p_maxval(maxval(abs(stats_cpy(:,6 ,1) - 0.5d0*(analyticalX+analyticalY+analyticalZ)))))
      call message("        > max absolute difference R33: ",p_maxval(maxval(abs(stats_cpy(:,72,1) - analyticalZ))))
      call message("        > max absolute difference wT1:  ",p_maxval(maxval(abs(stats_sca_cpy(:,25,1,1) - analyticalwT))))
      call assert(maxval(abs(stats_cpy(:,42,1) - analyticalX)) < 1.d-4,'R11 pressure transport')
      call assert(maxval(abs(stats_cpy(:,57,1) - analyticalY)) < 1.d-4,'R22 pressure transport')
      call message("    > Scale 2:")
      analyticalZ = -2.d0*0.5d0*coefs(3,3)*coefs(4,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))*wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))
      analyticalwT = -0.5d0*coefs(4,3)*coefs(5,3)*wavenums(2)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(2)*(z1d-zmin-Lz/2.0))
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call message("        > max absolute difference TKE: ",p_maxval(maxval(abs(stats_cpy(:,6 ,2) - 0.5d0*(analyticalX+analyticalY+analyticalZ)))))
      call message("        > max absolute difference R33: ",p_maxval(maxval(abs(stats_cpy(:,72,2) - analyticalZ))))
      call message("        > max absolute difference wT1:  ",p_maxval(maxval(abs(stats_sca_cpy(:,25,1,2) - analyticalwT))))
      call assert(maxval(abs(stats_cpy(:,42,2) - analyticalX)) < 1.d-4,'R11 pressure transport')
      call assert(maxval(abs(stats_cpy(:,57,2) - analyticalY)) < 1.d-4,'R22 pressure transport')
      call message('    Pressure transport test finished. Verify pass/fail separately')
      call message("*** *** ***")

      ! Viscous transport
      call message("Molecular transport")
      ! stats_cpy(:,13,1) is ddz(ddz(<uiui>)), but viscous transport is also computed as ddz(<uisi3>):
      analyticalX = 1.d0/SM%Re*(&
        0.25d0*coefs(1,2)*coefs(1,2)*wavenums(1)**2.d0*(sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 - cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0) + & !ddz<u*dudz>
        0.25d0*coefs(1,2)*coefs(3,2)*wavenums(1)**2.d0*(cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0) + & !ddz<u*dwdx>
        0.25d0*coefs(2,2)*coefs(2,2)*wavenums(1)**2.d0*(sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 - cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0) + & !ddz<v*dvdz>
        0.d0                                                                                                                                 + & !ddz<v*dwdy>
        0.50d0*coefs(3,2)*coefs(3,2)*wavenums(1)**2.d0*(cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0))    !ddz<w*dwdz>
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call message("        > max absolute difference TKE: ",maxval(abs(stats_cpy(:,7 ,1) - analyticalX)))
      analyticalX = 0.5d0/SM%Re*(coefs(1,2)*wavenums(1))**2.d0*(sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 - cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0)
      analyticalY = 0.5d0/SM%Re*(coefs(2,2)*wavenums(1))**2.d0*(sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 - cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0)
      analyticalZ = 0.5d0/SM%Re*(coefs(3,2)*wavenums(1))**2.d0*(cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0)
      analyticalT = 0.5d0/(SM%PrandtlFluid*SM%Re)*(coefs(5,2)*wavenums(1))**2.d0*(cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0)
      analyticalwT = 0.25d0/SM%Re*coefs(3,2)*coefs(5,2)*wavenums(1)**2.d0*(cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0) 
      call message("        > max absolute difference TKE_b: ",maxval(abs(stats_cpy(:,13,1) - 0.5d0*(analyticalX + analyticalY + analyticalZ))))
      call message("        > max absolute difference R11: ",maxval(abs(stats_cpy(:,43,1) - analyticalX)))
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(:,58,1) - analyticalY)))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(:,73,1) - analyticalZ)))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,5 ,1,1) - analyticalT)))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,27,1,1) - analyticalwT)))
      analyticalwT = 0.25d0/(SM%PrandtlFluid*SM%Re)*coefs(3,2)*coefs(5,2)*wavenums(1)**2.d0*(cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0)
      call message("        > max absolute difference wT2: ",maxval(abs(stats_sca_cpy(:,26,1,1) - analyticalwT)))
      call message("    > Scale 2:")
      analyticalX = 1.d0/SM%Re*(&
        0.25d0*coefs(1,3)*coefs(1,3)*wavenums(2)**2.d0*(sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 - cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0) + & !ddz<u*dudz>
        0.25d0*coefs(1,3)*coefs(3,3)*wavenums(2)**2.d0*(cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0) + & !ddz<u*dwdx>
        0.25d0*coefs(2,3)*coefs(2,3)*wavenums(2)**2.d0*(sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 - cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0) + & !ddz<v*dvdz>
        0.d0                                                                                                                                 + & !ddz<v*dwdy>
        0.50d0*coefs(3,3)*coefs(3,3)*wavenums(2)**2.d0*(cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0))    !ddz<w*dwdz>
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call message("        > max absolute difference TKE: ",maxval(abs(stats_cpy(:,7 ,2) - analyticalX)))
      analyticalX = 0.5d0/SM%Re*coefs(1,3)**2.d0*wavenums(2)**2.d0*(sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 - cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0)
      analyticalY = 0.5d0/SM%Re*coefs(2,3)**2.d0*wavenums(2)**2.d0*(sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 - cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0)
      analyticalZ = 0.5d0/SM%Re*coefs(3,3)**2.d0*wavenums(2)**2.d0*(cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0)
      analyticalT = 0.5d0/(SM%PrandtlFluid*SM%Re)*(coefs(5,3)*wavenums(2))**2.d0*(cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0)
      analyticalwT = 0.25d0/SM%Re*coefs(3,3)*coefs(5,3)*wavenums(2)**2.d0*(cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0) 
      call message("        > max absolute difference TKE_b: ",maxval(abs(stats_cpy(:,13,2) - 0.5d0*(analyticalX + analyticalY + analyticalZ))))
      call message("        > max absolute difference R11: ",maxval(abs(stats_cpy(:,43,2) - analyticalX)))
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(:,58,2) - analyticalY)))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(:,73,2) - analyticalZ)))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,5 ,1,2) - analyticalT)))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,27,1,2) - analyticalwT)))
      analyticalwT = 0.25d0/(SM%PrandtlFluid*SM%Re)*coefs(3,3)*coefs(5,3)*wavenums(2)**2.d0*(cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 - sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0)
      call message("        > max absolute difference wT2: ",maxval(abs(stats_sca_cpy(:,26,1,2) - analyticalwT)))
      call message("*** *** ***")

      ! Viscous dissipation
      call message("Molecular destruction")
      analyticalX  = 2.d0/SM%Re*(0.5d0*coefs(1,2)*wavenums(1))**2.d0*(2.d0*cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 +      sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0) 
      analyticalY  = 2.d0/SM%Re*(0.5d0*coefs(2,2)*wavenums(1))**2.d0*(2.d0*cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 +      sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0) 
      analyticalZ  = 2.d0/SM%Re*(0.5d0*coefs(3,2)*wavenums(1))**2.d0*(     cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 + 2.d0*sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0) 
      analyticalT  = 2.d0/(SM%Re*SM%PrandtlFluid)*(0.5d0*coefs(5,2)*wavenums(1))**2.d0*(cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 + 2.d0*sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0) 
      analyticalwT = (0.5d0*wavenums(1))**2.d0*coefs(3,2)*coefs(5,2)*(cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 + 2.d0*sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0) 
      call assert(maxval(abs(stats_cpy(:,12,1) - 0.5d0*(analyticalX + analyticalY + analyticalZ))) < 1.d-14,'molecular destruction, TKEbr')
      call assert(maxval(abs(stats_cpy(:,44,1) - analyticalX)) < 1.d-14,'molecular destruction, R11r')
      call assert(maxval(abs(stats_cpy(:,59,1) - analyticalY)) < 1.d-14,'molecular destruction, R22r')
      call assert(maxval(abs(stats_cpy(:,74,1) - analyticalZ)) < 1.d-14,'molecular destruction, R33r')
      call assert(maxval(abs(stats_sca_cpy(:,7 ,1,1) - analyticalT)) < 1.d-14,'molecular destruction, RTTr')
      call assert(maxval(abs(stats_sca_cpy(:,31,1,1) - analyticalwT/(SM%Re*SM%PrandtlFluid))) < 1.d-14,'molecular destruction, wT1r')
      call assert(maxval(abs(stats_sca_cpy(:,32,1,1) - analyticalwT/SM%Re)) < 1.d-14,'molecular destruction, wT2r')
      analyticalX = 0.5d0*(analyticalX + analyticalY + analyticalZ) + &
        1.d0/SM%Re*((0.5d0*wavenums(1))**2.d0*(coefs(1,2)**2.d0 + coefs(2,2)**2.d0 + coefs(3,2)**2.d0)*cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 + & ! (dudx)^2
                    -2.d0*(0.5d0*wavenums(1))**2.d0*coefs(1,2)*coefs(3,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 ) ! 2*dudz*dwdx
      call assert(maxval(abs(stats_cpy(:,8 ,1) - analyticalX)) < 1.d-14,'molecular destruction, TKEr')
      call message("    > Scale 2:")
      analyticalX  = 2.d0/SM%Re*(0.5d0*coefs(1,3)*wavenums(2))**2.d0*(2.d0*cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 +      sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0) 
      analyticalY  = 2.d0/SM%Re*(0.5d0*coefs(2,3)*wavenums(2))**2.d0*(2.d0*cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 +      sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0) 
      analyticalZ  = 2.d0/SM%Re*(0.5d0*coefs(3,3)*wavenums(2))**2.d0*(     cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 + 2.d0*sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0) 
      analyticalT  = 2.d0/(SM%Re*SM%PrandtlFluid)*(0.5d0*coefs(5,3)*wavenums(2))**2.d0*(cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 + 2.d0*sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0) 
      analyticalwT = (0.5d0*wavenums(2))**2.d0*coefs(3,3)*coefs(5,3)*(cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 + 2.d0*sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0) 
      call assert(maxval(abs(stats_cpy(:,12,2) - 0.5d0*(analyticalX + analyticalY + analyticalZ))) < 1.d-14,'molecular destruction, TKEbr')
      call assert(maxval(abs(stats_cpy(:,44,2) - analyticalX)) < 1.d-14,'molecular destruction, R11r')
      call assert(maxval(abs(stats_cpy(:,59,2) - analyticalY)) < 1.d-14,'molecular destruction, R22r')
      call assert(maxval(abs(stats_cpy(:,74,2) - analyticalZ)) < 1.d-14,'molecular destruction, R33r')
      call assert(maxval(abs(stats_sca_cpy(:,7 ,1,2) - analyticalT)) < 1.d-14,'molecular destruction, RTTr')
      call assert(maxval(abs(stats_sca_cpy(:,31,1,2) - analyticalwT/(SM%Re*SM%PrandtlFluid))) < 1.d-14,'molecular destruction, wT1r')
      call assert(maxval(abs(stats_sca_cpy(:,32,1,2) - analyticalwT/SM%Re)) < 1.d-14,'molecular destruction, wT2r')
      analyticalX = 0.5d0*(analyticalX + analyticalY + analyticalZ) + &
        1.d0/SM%Re*((0.5d0*wavenums(2))**2.d0*(coefs(1,3)**2.d0 + coefs(2,3)**2.d0 + coefs(3,3)**2.d0)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 + & ! (dudx)^2
                    -2.d0*(0.5d0*wavenums(2))**2.d0*coefs(1,3)*coefs(3,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 ) ! 2*dudz*dwdx
      call assert(maxval(abs(stats_cpy(:,8 ,2) - analyticalX)) < 1.d-14,'molecular destruction, TKEs')
      call message("*** *** ***")

      ! SGS transport
      call message("SGS transport")
      call assert(maxval(abs(stats_cpy(:,45,1))) < 1.d-14,'SGS transport, R11r')
      call assert(maxval(abs(stats_cpy(:,60,1))) < 1.d-14,'SGS transport, R22r')
      call assert(maxval(abs(stats_cpy(:,75,1))) < 1.d-14,'SGS transport, R33r')
      call assert(maxval(abs(stats_cpy(:,9 ,1))) < 1.d-14,'SGS transport, TKEr')
      call assert(maxval(abs(stats_sca_cpy(:,6 ,1,1))) < 1.d-14,'SGS transport, TTr')
      call assert(maxval(abs(stats_sca_cpy(:,28,1,1))) < 1.d-14,'SGS transport, wT1r')
      call assert(maxval(abs(stats_sca_cpy(:,29,1,1))) < 1.d-14,'SGS transport, wT2r')
      call message("    > Scale 2:")
      ! There is only transport if wavenums(2) == wavenums(3). Check for this
      call assert(abs(wavenums(2) - wavenums(3)) < 1.d-14,'Must set tau wavenumbers equal to scale 2 wavenumbers to test SGS transport')
      analyticalX = 0.d0
      analyticalY =  0.5d0*coefs(2,3)*taucoefs(5)*(&
        wavenums(2)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0)) + &
        wavenums(3)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))*sin(wavenums(3)*(z1d-zmin-Lz/2.0)))
      analyticalZ = -0.5d0*coefs(3,3)*taucoefs(6)*(&
        wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0)) - &
        wavenums(3)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))*sin(wavenums(3)*(z1d-zmin-Lz/2.0)))
      analyticalT = -0.5d0*coefs(5,3)*taucoefs(9)*(&
        wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0)) - &
        wavenums(3)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))*sin(wavenums(3)*(z1d-zmin-Lz/2.0)))
      analyticalwT = -0.25d0*coefs(3,3)*taucoefs(9)*(&
        wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0)) - &
        wavenums(3)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))*sin(wavenums(3)*(z1d-zmin-Lz/2.0)))
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call assert(maxval(abs(stats_cpy(:,45,2))) < 1.d-14,'SGS transport, R11r')
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(:,60,2) - analyticalY)))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(:,75,2) - analyticalZ)))
      call message("        > max absolute difference TKE: ",maxval(abs(stats_cpy(:,9 ,2) - 0.5d0*(analyticalX + analyticalY + analyticalZ))))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,6 ,1,2) - analyticalT)))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,28,1,2) - analyticalwT)))
      analyticalwT = -0.25d0*coefs(5,3)*taucoefs(6)*(&
        wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0)) - &
        wavenums(3)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))*sin(wavenums(3)*(z1d-zmin-Lz/2.0)))
      call message("        > max absolute difference wT2: ",maxval(abs(stats_sca_cpy(:,29,1,2) - analyticalwT)))
      call message("*** *** ***")

      ! SGS dissipation
      call message("SGS dissipation")
      call assert(maxval(abs(stats_cpy(:,10,1))) < 1.d-14,'SGS dissipation, TKEr')
      call assert(maxval(abs(stats_cpy(:,46,1))) < 1.d-14,'SGS dissipation, R11r')
      call assert(maxval(abs(stats_cpy(:,61,1))) < 1.d-14,'SGS dissipation, R22r')
      call assert(maxval(abs(stats_cpy(:,76,1))) < 1.d-14,'SGS dissipation, R33r')
      call assert(maxval(abs(stats_sca_cpy(:,8 ,1,1))) < 1.d-14,'SGS dissipation, TTr')
      call assert(maxval(abs(stats_sca_cpy(:,33,1,1))) < 1.d-14,'SGS dissipation, wT1r')
      call assert(maxval(abs(stats_sca_cpy(:,34,1,1))) < 1.d-14,'SGS dissipation, wT2r')
      call message("    > Scale 2:")
      analyticalX = 0.5d0*coefs(1,3)*taucoefs(1)*wavenums(2)*&
        cos(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0))
      analyticalY = 0.5d0*coefs(2,3)*taucoefs(5)*wavenums(2)*&
        sin(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0))
      analyticalZ =  -0.5d0*coefs(3,3)*taucoefs(6)*wavenums(2)*&
        cos(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0))
      analyticalT = -0.5d0*coefs(5,3)*taucoefs(9)*wavenums(2)*&
        cos(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0))
      analyticalwT = -0.25d0*coefs(3,3)*taucoefs(9)*wavenums(2)*&
        cos(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0))
      call message("    > Make sure this scales as 6th order to confirm this is a numerical method error and not a bug:")
      call assert(maxval(abs(stats_cpy(:,46,2) - analyticalX)) < 1.d-14,'SGS dissipation, R11s')
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(:,61,2) - analyticalY)))
      call assert(maxval(abs(stats_cpy(:,76,2) - analyticalZ)) < 1.d-14,'SGS dissipation, R33s')
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,8 ,1,2) - analyticalT)))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,33,1,2) - analyticalwT)))
      analyticalwT = -0.25d0*coefs(5,3)*taucoefs(6)*wavenums(2)*&
        cos(wavenums(2)*(z1d-zmin-Lz/2.0))*cos(wavenums(3)*(z1d-zmin-Lz/2.0))
      call assert(maxval(abs(stats_sca_cpy(:,34,1,2) - analyticalwT)) < 1.d-14, 'SGS dissipation, wT2s')
      analyticalX = 0.5d0*(analyticalX + analyticalY + analyticalZ)
      call message("        > max absolute difference TKE: ",maxval(abs(stats_cpy(:,10 ,2) - analyticalX)))
      call message("*** *** ***")

      ! Force production
      call message("Force production")
      analyticalX = 0.5d0*coefs(1,2)*coefs(6,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0
      analyticalY = 0.5d0*coefs(2,2)*coefs(7,2)*cos(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0
      analyticalZ = 0.5d0*coefs(3,2)*coefs(8,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0
      analyticalwT = 0.25d0*coefs(5,2)*coefs(8,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0
      call assert(maxval(abs(stats_cpy(:,39,1) - analyticalX)) < 1.d-14,'Force production, R11r')
      call assert(maxval(abs(stats_cpy(:,54,1) - analyticalY)) < 1.d-14,'Force production, R22r')
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(:,69,1) - analyticalZ)))
      call message("        > max absolute difference TKE: ",maxval(abs(stats_cpy(:,3 ,1) - 0.5d0*(analyticalX+analyticalY+analyticalZ))))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,21,1,1) - analyticalwT)))
      call message("    > Scale 2:")
      analyticalX = 0.5d0*coefs(1,3)*coefs(6,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0
      analyticalY = 0.5d0*coefs(2,3)*coefs(7,3)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0
      analyticalZ = 0.5d0*coefs(3,3)*coefs(8,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0
      analyticalwT = 0.25d0*coefs(5,3)*coefs(8,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0
      call assert(maxval(abs(stats_cpy(:,39,2) - analyticalX)) < 1.d-14,'Force production, R11s')
      call assert(maxval(abs(stats_cpy(:,54,2) - analyticalY)) < 1.d-14,'Force production, R22s')
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(:,69,2) - analyticalZ)))
      call message("        > max absolute difference TKE: ",maxval(abs(stats_cpy(:,3 ,2) - 0.5d0*(analyticalX+analyticalY+analyticalZ))))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,21,1,2) - analyticalwT)))
      call message("*** *** ***")

      ! Buoyancy transfer
      call message("Buoyancy transfer")
      analyticalZ = -0.5d0*SM%buoyancyFact*coefs(3,2)*coefs(5,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.0))**2.d0 
      call assert(maxval(abs(stats_cpy(:,47,1)              )) < 1.d-14,'Buoyancy transfer, R11r')
      call assert(maxval(abs(stats_cpy(:,62,1)              )) < 1.d-14,'Buoyancy transfer, R22r')
      call assert(maxval(abs(stats_cpy(:,77,1) - analyticalZ)) < 1.d-14,'Buoyancy transfer, R33r')
      call assert(maxval(abs(stats_cpy(:,11,1) - 0.5d0*analyticalZ)) < 1.d-14,'Buoyancy transfer, R33r')
      call message("    > Scale 2:")
      analyticalZ = -0.5d0*SM%buoyancyFact*coefs(3,3)*coefs(5,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.0))**2.d0 
      call assert(maxval(abs(stats_cpy(:,47,2)              )) < 1.d-14,'Buoyancy transfer, R11s')
      call assert(maxval(abs(stats_cpy(:,62,2)              )) < 1.d-14,'Buoyancy transfer, R22s')
      call assert(maxval(abs(stats_cpy(:,77,2) - analyticalZ)) < 1.d-14,'Buoyancy transfer, R33r')
      call assert(maxval(abs(stats_cpy(:,11,2) - 0.5d0*analyticalZ)) < 1.d-14,'Buoyancy transfer, R33r')
      call message("*** *** ***")

      ! Buoyancy production
      call message("Buoyancy production")
      analyticalwT = 0.25d0*SM%buoyancyFact*(coefs(5,2)*sin(wavenums(1)*(z1d-zmin-Lz/2.0)))**2.d0
      call assert(maxval(abs(stats_sca_cpy(:,22,1,1) - analyticalwT)) < 1.d-14,'Buoyancy production, wTr')
      analyticalwT = 0.25d0*SM%buoyancyFact*(coefs(5,3)*sin(wavenums(2)*(z1d-zmin-Lz/2.0)))**2.d0
      call assert(maxval(abs(stats_sca_cpy(:,22,1,2) - analyticalwT)) < 1.d-14,'Buoyancy production, wTs')
      call message("*** *** ***")

      ! Pressure strain correlations
      call message("Pressure-strain correlations")
      analyticalX = -0.5d0*coefs(1,2)*coefs(4,2)*wavenums(1)*cos(wavenums(1)*(z1d-zmin-Lz/2.0))*sin(wavenums(1)*(z1d-zmin-Lz/2.0))
      analyticalZ =  0.5d0*coefs(3,2)*coefs(4,2)*wavenums(1)*cos(wavenums(1)*(z1d-zmin-Lz/2.0))*sin(wavenums(1)*(z1d-zmin-Lz/2.0))
      call assert(maxval(abs(stats_cpy(:,48,1) - analyticalX)) < 1.d-14,'Pressure-strain, R11r')
      call assert(maxval(abs(stats_cpy(:,63,1)              )) < 1.d-14,'Pressure-strain, R22r')
      call assert(maxval(abs(stats_cpy(:,78,1) - analyticalZ)) < 1.d-14,'Pressure-strain, R33r')
      analyticalX = -0.5d0*coefs(1,3)*coefs(4,3)*wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))*sin(wavenums(2)*(z1d-zmin-Lz/2.0))
      analyticalZ =  0.5d0*coefs(3,3)*coefs(4,3)*wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))*sin(wavenums(2)*(z1d-zmin-Lz/2.0))
      call assert(maxval(abs(stats_cpy(:,48,2) - analyticalX)) < 1.d-14,'Pressure-strain, R11s')
      call assert(maxval(abs(stats_cpy(:,63,2)              )) < 1.d-14,'Pressure-strain, R22s')
      call assert(maxval(abs(stats_cpy(:,78,2) - analyticalZ)) < 1.d-14,'Pressure-strain, R33s')
      call message("*** *** ***")

      ! Pressure scrambling
      call message("Pressure scrambling")
      analyticalwT = 0.25d0*coefs(4,2)*coefs(5,2)*wavenums(1)*cos(wavenums(1)*(z1d-zmin-Lz/2.0))*sin(wavenums(1)*(z1d-zmin-Lz/2.0))
      call assert(maxval(abs(stats_sca_cpy(:,30,1,1) - analyticalwT)) < 1.d-14,'Pressure scrambling, wTr')
      analyticalwT = 0.25d0*coefs(4,3)*coefs(5,3)*wavenums(2)*cos(wavenums(2)*(z1d-zmin-Lz/2.0))*sin(wavenums(2)*(z1d-zmin-Lz/2.0))
      call assert(maxval(abs(stats_sca_cpy(:,30,1,2) - analyticalwT)) < 1.d-14,'Pressure scrambling, wTs')
      call message("*** *** ***")
      
      ! DEFINE NEW FIELDS TO TEST UNSTEADY TERMS 
      call message("Unsteady terms")
      zids = [33,34,35]
      steps = [3,4,5]
      omega = pi
      SM%dt = z1d(2)-z1d(1)
      do i = 1,3,2
          SM%step = steps(i)
          SM%u  = sin(omega*z1d(zids(i)))*(coefs(1,2)*cos(wavenums(1)*x ) + coefs(1,3)*cos(wavenums(2)*x ))
          SM%v  = sin(omega*z1d(zids(i)))*(coefs(2,2)*cos(wavenums(1)*x ) + coefs(2,3)*cos(wavenums(2)*x ))
          SM%w  = sin(omega*z1d(zids(i)))*(coefs(3,2)*cos(wavenums(1)*xE) + coefs(3,3)*cos(wavenums(2)*xE))
          SM%wC = sin(omega*z1d(zids(i)))*(coefs(3,2)*cos(wavenums(1)*x ) + coefs(3,3)*cos(wavenums(2)*x ))
          SM%T  = sin(omega*z1d(zids(i)))*(coefs(5,2)*cos(wavenums(1)*x ) + coefs(5,3)*cos(wavenums(2)*x ))
          call stats%set_tidx(1)
          call stats%compute_stats()
          !if (nrank == 0) then
          !    call stats%print_var(1,1,1,1,1,10,i,'u')
          !    call stats%print_var(1,1,1,1,1,10,i,'v')
          !    call stats%print_var(1,1,1,1,1,10,i,'w')
          !end if
          !call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end do
      analyticalX  = coefs(1,2)**2.d0*omega*sin(omega*z1d(zids(2)))*cos(omega*z1d(zids(2)))
      analyticalY  = coefs(2,2)**2.d0*omega*sin(omega*z1d(zids(2)))*cos(omega*z1d(zids(2)))
      analyticalZ  = coefs(3,2)**2.d0*omega*sin(omega*z1d(zids(2)))*cos(omega*z1d(zids(2)))
      analyticalT  = coefs(5,2)**2.d0*omega*sin(omega*z1d(zids(2)))*cos(omega*z1d(zids(2)))
      analyticalwT = coefs(3,2)*coefs(5,2)*omega*sin(omega*z1d(zids(2)))*cos(omega*z1d(zids(2)))
      call stats%set_tidx(1)
      call stats%copy_stats(stats_cpy,stats_sca_cpy)
      if (nrank == 0) then
          print*, "analytical ddt<TT>:", analyticalT(1:10)
          print*, "computed ddt<TT>:", stats_sca_cpy(1:10,1,1,1)
      end if
      call message("        > max absolute difference TKE: ",maxval(abs(stats_cpy(:,1 ,1) - 0.5d0*(analyticalX+analyticalY+analyticalZ))))
      call message("        > max absolute difference R11: ",maxval(abs(stats_cpy(:,37,1) - analyticalX)))
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(:,52,1) - analyticalY)))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(:,67,1) - analyticalZ)))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,1 ,1,1) - analyticalT)))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,18,1,1) - analyticalwT)))
      call message("    > Scale 2:")
      analyticalX  = coefs(1,3)**2.d0*omega*sin(omega*z1d(zids(2)))*cos(omega*z1d(zids(2)))
      analyticalY  = coefs(2,3)**2.d0*omega*sin(omega*z1d(zids(2)))*cos(omega*z1d(zids(2)))
      analyticalZ  = coefs(3,3)**2.d0*omega*sin(omega*z1d(zids(2)))*cos(omega*z1d(zids(2)))
      analyticalT  = coefs(5,3)**2.d0*omega*sin(omega*z1d(zids(2)))*cos(omega*z1d(zids(2)))
      analyticalwT = coefs(3,3)*coefs(5,3)*omega*sin(omega*z1d(zids(2)))*cos(omega*z1d(zids(2)))
      call message("        > max absolute difference TKE: ",maxval(abs(stats_cpy(:,1 ,2) - 0.5d0*(analyticalX+analyticalY+analyticalZ))))
      call message("        > max absolute difference R11: ",maxval(abs(stats_cpy(:,37,2) - analyticalX)))
      call message("        > max absolute difference R22: ",maxval(abs(stats_cpy(:,52,2) - analyticalY)))
      call message("        > max absolute difference R33: ",maxval(abs(stats_cpy(:,67,2) - analyticalZ)))
      call message("        > max absolute difference TT:  ",maxval(abs(stats_sca_cpy(:,1 ,1,2) - analyticalT)))
      call message("        > max absolute difference wT1: ",maxval(abs(stats_sca_cpy(:,18,1,2) - analyticalwT)))
      call message("*** *** ***")

      call message("==========================================================")
      call message(0,"Finalizing simulation")
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call stats%destroy()
      call SM%destroy()
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call message(0,"Simulation finalized")
   
      !deallocate(SM)
      deallocate(tau11,tau12,tau13,tau22,tau23,tau33)
      deallocate(stats_cpy,stats_sca_cpy)
      deallocate(analyticalX,analyticalY,analyticalZ)
      deallocate(z1d)
      nullify(x,y,z,xE,yE,zE)
    end if
    
    call MPI_Finalize(ierr)

end program test_stats_xy
