module RK_coefficients
       use kind_parameters, only: rkind
       real(rkind), parameter :: b01 = 0.39175222657189d0 , b12 = 0.368410593050371d0, b23 = 0.25189177427169d0, b34 = 0.54497475022852d0
       real(rkind), parameter :: b35 = 0.06369246866629d0 , b45 = 0.22600748323690d0
       
       real(rkind), parameter :: a20 = 0.444370493651235d0, a21 = 0.555629506348765d0
       real(rkind), parameter :: a30 = 0.620101851488403d0, a32 = 0.379898148511597d0
       real(rkind), parameter :: a40 = 0.17807995439313d0 , a43 = 0.821920045606868d0
       real(rkind), parameter :: a52 = 0.517231671970585d0, a53 = 0.096059710526147d0, a54 = 0.386708617503269d0
end module
