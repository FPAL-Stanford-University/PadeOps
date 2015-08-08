module fields
    use kind_parameters, only: rkind
    use inputParameters 
    use fft_3d_Stuff, only: fft_3d 

    real(rkind)   , dimension(:,:,:), allocatable :: fieldsPhys
    complex(rkind), dimension(:,:,:), allocatable :: fieldsSpec





end module
