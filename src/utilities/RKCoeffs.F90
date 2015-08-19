module RKCoeffs

    use kind_parameters, only: rkind
    implicit none

    integer,                   parameter :: RK45_steps = 5
    real(rkind), dimension(5), parameter :: RK45_A = [0.0_rkind, &
                                                      -6234157559845.0_rkind/12983515589748.0_rkind, &
                                                      -6194124222391.0_rkind/4410992767914.0_rkind,  &
                                                     -31623096876824.0_rkind/15682348800105.0_rkind, &
                                                     -12251185447671.0_rkind/11596622555746.0_rkind ]

    real(rkind), dimension(5), parameter :: RK45_B = [ 494393426753.0_rkind/4806282396855.0_rkind,  &
                                                      4047970641027.0_rkind/5463924506627.0_rkind,  &
                                                      9795748752853.0_rkind/13190207949281.0_rkind, &
                                                      4009051133189.0_rkind/8539092990294.0_rkind,  &
                                                      1348533437543.0_rkind/7166442652324.0_rkind   ]
end module
