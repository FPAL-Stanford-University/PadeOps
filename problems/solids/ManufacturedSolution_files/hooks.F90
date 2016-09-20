module ManufacturedSolution_data
    use kind_parameters,  only: rkind
    use constants,        only: one,eight
    implicit none

    real(rkind) :: sigma_0 = one
    real(rkind) :: gamma, p_infty, rho_0
    real(rkind) :: muL, lamL, cL
    real(rkind) :: dtprob

contains

    !******************************************************************************
    !*                      Code generated with sympy 0.7.5                       *
    !*                                                                            *
    !*              See http://www.sympy.org/ for more information.               *
    !*                                                                            *
    !*                       This file is part of 'project'                       *
    !******************************************************************************
    
    REAL*8 function g_source(gamma, mu, p_infty, rho_0, sigma_0, t, x)
    implicit none
    REAL*8, intent(in) :: gamma
    REAL*8, intent(in) :: mu
    REAL*8, intent(in) :: p_infty
    REAL*8, intent(in) :: rho_0
    REAL*8, intent(in) :: sigma_0
    REAL*8, intent(in) :: t
    REAL*8, intent(in) :: x
    
    g_source = sigma_0**2*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0)*( &
          200.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) - &
          200.0d0*x + 100.0d0)*exp(-2.0d0*(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)**2*(-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**2) + sigma_0*sqrt(( &
          gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0)*(200.0d0*t*sqrt((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)*( &
          -sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0 &
          )*mu) + 1)) + 20.0d0*sigma_0*sqrt((gamma*p_infty + (4.0d0/3.0d0)* &
          mu)/rho_0)*(-10.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma &
          *p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          (gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**2)
    
    end function
    
    !******************************************************************************
    !*                      Code generated with sympy 0.7.5                       *
    !*                                                                            *
    !*              See http://www.sympy.org/ for more information.               *
    !*                                                                            *
    !*                       This file is part of 'project'                       *
    !******************************************************************************
    
    REAL*8 function momentum_source(gamma, mu, p_infty, rho_0, sigma_0, t, x)
    implicit none
    REAL*8, intent(in) :: gamma
    REAL*8, intent(in) :: mu
    REAL*8, intent(in) :: p_infty
    REAL*8, intent(in) :: rho_0
    REAL*8, intent(in) :: sigma_0
    REAL*8, intent(in) :: t
    REAL*8, intent(in) :: x
    
    momentum_source = gamma*p_infty*sigma_0*(200.0d0*t*sqrt((gamma*p_infty + &
          (4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*(1.0/(-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1))**gamma*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma &
          *p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          (gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)) + mu*sigma_0*((2.0d0/ &
          3.0d0)*(-1 + (-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-4))/((-sigma_0*exp(-(-10.0d0* &
          t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x &
          - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**( &
          2.0d0/3.0d0) - 2.0d0/3.0d0*(-1 + (-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))/((-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1)**(-2))**(1.0d0/3.0d0))*(200.0d0*t*sqrt((gamma*p_infty + ( &
          4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*exp(-(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**2 &
          ) + mu*(-8.0d0/9.0d0*sigma_0*(-1 + (-sigma_0*exp(-(-10.0d0*t*sqrt &
          ((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-4))*( &
          200.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) - &
          200.0d0*x + 100.0d0)*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma &
          *p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          (gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)*((-sigma_0*exp(-(-10.0d0* &
          t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x &
          - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**( &
          2.0d0/3.0d0)) + (4.0d0/9.0d0)*sigma_0*(-1 + (-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**( &
          -2))*(200.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) - &
          200.0d0*x + 100.0d0)*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma &
          *p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          (gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)*((-sigma_0*exp(-(-10.0d0* &
          t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x &
          - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**( &
          1.0d0/3.0d0)) - 4.0d0/3.0d0*sigma_0*(200.0d0*t*sqrt((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)*( &
          -sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0 &
          )*mu) + 1)**3*((-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(1.0d0/3.0d0)) + (8.0d0/ &
          3.0d0)*sigma_0*(200.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu) &
          /rho_0) - 200.0d0*x + 100.0d0)*exp(-(-10.0d0*t*sqrt((gamma* &
          p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          ((gamma*p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**5*((-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1)**(-2))**(2.0d0/3.0d0)))/(-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1) + 20.0d0*rho_0* &
          sigma_0**2*((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0)*(-10.0d0*t* &
          sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) + 10.0d0*x - 5.0d0 &
          )*exp(-2.0d0*(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0* &
          mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0 &
          )*mu)**2*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**2) + 20.0d0*rho_0*sigma_0*(( &
          gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0)*(-10.0d0*t*sqrt((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)/rho_0) + 10.0d0*x - 5.0d0)*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)*( &
          -sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0 &
          )*mu) + 1)) + sigma_0**3*(200.0d0*t*sqrt((gamma*p_infty + (4.0d0/ &
          3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*exp(-3.0d0*(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)**2*(-sigma_0*exp(- &
          (-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**2 &
          ) + sigma_0**2*(400.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu) &
          /rho_0) - 400.0d0*x + 200.0d0)*exp(-2.0d0*(-10.0d0*t*sqrt((gamma* &
          p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          ((gamma*p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1))
    
    end function
    
    !******************************************************************************
    !*                      Code generated with sympy 0.7.5                       *
    !*                                                                            *
    !*              See http://www.sympy.org/ for more information.               *
    !*                                                                            *
    !*                       This file is part of 'project'                       *
    !******************************************************************************
    
    REAL*8 function energy_source(gamma, mu, p_infty, rho_0, sigma_0, t, x)
    implicit none
    REAL*8, intent(in) :: gamma
    REAL*8, intent(in) :: mu
    REAL*8, intent(in) :: p_infty
    REAL*8, intent(in) :: rho_0
    REAL*8, intent(in) :: sigma_0
    REAL*8, intent(in) :: t
    REAL*8, intent(in) :: x
    
    energy_source = rho_0*sigma_0**2*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu) &
          /rho_0)*(200.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0 &
          ) - 200.0d0*x + 100.0d0)*((1.0d0/4.0d0)*mu*((2 + (-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**( &
          -4))/((-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(2.0d0/3.0d0) - 2*(2 + ( &
          -sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0 &
          )*mu) + 1)**(-2))/((-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty &
          + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(1.0d0/3.0d0) + 3)/rho_0 &
          + (1.0d0/2.0d0)*sigma_0**2*exp(-2.0d0*(-10.0d0*t*sqrt((gamma* &
          p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/( &
          rho_0*(gamma*p_infty + (4.0d0/3.0d0)*mu)) + (gamma*p_infty + &
          p_infty*((1.0/(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1))**gamma - 1))*(-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)/( &
          rho_0*(gamma - 1)))*exp(-2.0d0*(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)**2*(-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**2) + rho_0*sigma_0* &
          sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0)*(200.0d0*t*sqrt(( &
          gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)* &
          ((1.0d0/4.0d0)*mu*((2 + (-sigma_0*exp(-(-10.0d0*t*sqrt((gamma* &
          p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/( &
          gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-4))/((-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**( &
          -2))**(2.0d0/3.0d0) - 2*(2 + (-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))/((-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1)**(-2))**(1.0d0/3.0d0) + 3)/rho_0 + (1.0d0/2.0d0)*sigma_0** &
          2*exp(-2.0d0*(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0* &
          mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(rho_0*(gamma*p_infty + (4.0d0 &
          /3.0d0)*mu)) + (gamma*p_infty + p_infty*((1.0/(-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1))** &
          gamma - 1))*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)/(rho_0*(gamma - 1)))*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)*( &
          -sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0 &
          )*mu) + 1)) + rho_0*sigma_0*sqrt((gamma*p_infty + (4.0d0/3.0d0)* &
          mu)/rho_0)*(gamma*p_infty*sigma_0*(200.0d0*t*sqrt((gamma*p_infty &
          + (4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*(1.0/(-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1))**gamma*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(rho_0*( &
          gamma - 1)*(gamma*p_infty + (4.0d0/3.0d0)*mu)) + (1.0d0/4.0d0)*mu &
          *(-4.0d0/3.0d0*sigma_0*(2 + (-sigma_0*exp(-(-10.0d0*t*sqrt((gamma &
          *p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          (gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-4))*(200.0d0*t*sqrt(( &
          gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu &
          )*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)*((-sigma_0*exp(-(-10.0d0*t*sqrt &
          ((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(2.0d0 &
          /3.0d0)) + (4.0d0/3.0d0)*sigma_0*(2 + (-sigma_0*exp(-(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))*( &
          200.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) - &
          200.0d0*x + 100.0d0)*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma &
          *p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          (gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)*((-sigma_0*exp(-(-10.0d0* &
          t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x &
          - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**( &
          1.0d0/3.0d0)) - 4*sigma_0*(200.0d0*t*sqrt((gamma*p_infty + (4.0d0 &
          /3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0* &
          t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x &
          - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**3*(( &
          -sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0 &
          )*mu) + 1)**(-2))**(1.0d0/3.0d0)) + 4*sigma_0*(200.0d0*t*sqrt(( &
          gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu &
          )*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**5*((-sigma_0*exp(-(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(2.0d0 &
          /3.0d0)))/rho_0 + (1.0d0/2.0d0)*sigma_0**2*(400.0d0*t*sqrt((gamma &
          *p_infty + (4.0d0/3.0d0)*mu)/rho_0) - 400.0d0*x + 200.0d0)*exp( &
          -2.0d0*(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(rho_0*(gamma*p_infty + (4.0d0/ &
          3.0d0)*mu)) - sigma_0*(gamma*p_infty + p_infty*((1.0/(-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1))**gamma - 1))*(200.0d0*t*sqrt((gamma*p_infty + (4.0d0/ &
          3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(rho_0*(gamma - 1)*(gamma*p_infty + (4.0d0/3.0d0)*mu)))*exp( &
          -(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)*( &
          -sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0 &
          )*mu) + 1)) + 20.0d0*rho_0*sigma_0*sqrt((gamma*p_infty + (4.0d0/ &
          3.0d0)*mu)/rho_0)*(-10.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)* &
          mu)/rho_0) + 10.0d0*x - 5.0d0)*((1.0d0/4.0d0)*mu*((2 + (-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1)**(-4))/((-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(2.0d0/3.0d0) - 2*(2 + ( &
          -sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0 &
          )*mu) + 1)**(-2))/((-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty &
          + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(1.0d0/3.0d0) + 3)/rho_0 &
          + (1.0d0/2.0d0)*sigma_0**2*exp(-2.0d0*(-10.0d0*t*sqrt((gamma* &
          p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/( &
          rho_0*(gamma*p_infty + (4.0d0/3.0d0)*mu)) + (gamma*p_infty + &
          p_infty*((1.0/(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1))**gamma - 1))*(-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)/( &
          rho_0*(gamma - 1)))*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma &
          *p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          (gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**2) + rho_0*(20.0d0*gamma &
          *p_infty*sigma_0*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0)*( &
          -10.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) + 10.0d0 &
          *x - 5.0d0)*(1.0/(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1))**gamma*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(rho_0*(gamma - 1)*(gamma*p_infty + (4.0d0/3.0d0)*mu)) + ( &
          1.0d0/4.0d0)*mu*(-26.6666666666667d0*sigma_0*sqrt((gamma*p_infty &
          + (4.0d0/3.0d0)*mu)/rho_0)*(2 + (-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-4))*(-10.0d0*t* &
          sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) + 10.0d0*x - 5.0d0 &
          )*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu &
          )*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)*((-sigma_0*exp(-(-10.0d0*t*sqrt &
          ((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(2.0d0 &
          /3.0d0)) + 26.6666666666667d0*sigma_0*sqrt((gamma*p_infty + ( &
          4.0d0/3.0d0)*mu)/rho_0)*(2 + (-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))*(-10.0d0*t* &
          sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) + 10.0d0*x - 5.0d0 &
          )*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu &
          )*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)*((-sigma_0*exp(-(-10.0d0*t*sqrt &
          ((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(1.0d0 &
          /3.0d0)) - 80.0d0*sigma_0*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu) &
          /rho_0)*(-10.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0 &
          ) + 10.0d0*x - 5.0d0)*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma &
          *p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          (gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**3*((-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**( &
          -2))**(1.0d0/3.0d0)) + 80.0d0*sigma_0*sqrt((gamma*p_infty + ( &
          4.0d0/3.0d0)*mu)/rho_0)*(-10.0d0*t*sqrt((gamma*p_infty + (4.0d0/ &
          3.0d0)*mu)/rho_0) + 10.0d0*x - 5.0d0)*exp(-(-10.0d0*t*sqrt((gamma &
          *p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          ((gamma*p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**5*((-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1)**(-2))**(2.0d0/3.0d0)))/rho_0 + 20.0d0*sigma_0**2*sqrt(( &
          gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0)*(-10.0d0*t*sqrt((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)/rho_0) + 10.0d0*x - 5.0d0)*exp(-2.0d0 &
          *(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(rho_0*(gamma*p_infty + (4.0d0/3.0d0)*mu)) &
          - 20.0d0*sigma_0*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0)*( &
          gamma*p_infty + p_infty*((1.0/(-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1))**gamma - 1))*( &
          -10.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) + 10.0d0 &
          *x - 5.0d0)*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(rho_0*( &
          gamma - 1)*(gamma*p_infty + (4.0d0/3.0d0)*mu)))/(-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1) - &
          sigma_0*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0)*(-mu*(( &
          2.0d0/3.0d0)*(-1 + (-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty &
          + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-4))/((-sigma_0*exp(-(-10.0d0* &
          t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x &
          - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**( &
          2.0d0/3.0d0) - 2.0d0/3.0d0*(-1 + (-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))/((-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1)**(-2))**(1.0d0/3.0d0))/(-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1) - p_infty*((1.0/( &
          -sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0 &
          )*mu) + 1))**gamma - 1))*(200.0d0*t*sqrt((gamma*p_infty + (4.0d0/ &
          3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) - sigma_0*sqrt((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)/rho_0)*(-gamma*p_infty*sigma_0*( &
          200.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) - &
          200.0d0*x + 100.0d0)*(1.0/(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma* &
          p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/( &
          gamma*p_infty + (4.0d0/3.0d0)*mu) + 1))**gamma*exp(-(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)) - &
          mu*sigma_0*((2.0d0/3.0d0)*(-1 + (-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-4))/((-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1)**(-2))**(2.0d0/3.0d0) - 2.0d0/3.0d0*(-1 + (-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**( &
          -2))/((-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(1.0d0/3.0d0))*(200.0d0* &
          t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + &
          100.0d0)*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/ &
          3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**2) - mu*(-8.0d0/9.0d0*sigma_0*( &
          -1 + (-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-4))*(200.0d0*t*sqrt((gamma* &
          p_infty + (4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)*( &
          -sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0 &
          *mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0 &
          )*mu) + 1)*((-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(2.0d0/3.0d0)) + (4.0d0/ &
          9.0d0)*sigma_0*(-1 + (-sigma_0*exp(-(-10.0d0*t*sqrt((gamma* &
          p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/( &
          gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))*(200.0d0*t*sqrt(( &
          gamma*p_infty + (4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu &
          )*(-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)*((-sigma_0*exp(-(-10.0d0*t*sqrt &
          ((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(1.0d0 &
          /3.0d0)) - 4.0d0/3.0d0*sigma_0*(200.0d0*t*sqrt((gamma*p_infty + ( &
          4.0d0/3.0d0)*mu)/rho_0) - 200.0d0*x + 100.0d0)*exp(-(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/((gamma*p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-( &
          -10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + &
          10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**3* &
          ((-sigma_0*exp(-(-10.0d0*t*sqrt((gamma*p_infty + &
          1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma* &
          p_infty + (4.0d0/3.0d0)*mu) + 1)**(-2))**(1.0d0/3.0d0)) + (8.0d0/ &
          3.0d0)*sigma_0*(200.0d0*t*sqrt((gamma*p_infty + (4.0d0/3.0d0)*mu) &
          /rho_0) - 200.0d0*x + 100.0d0)*exp(-(-10.0d0*t*sqrt((gamma* &
          p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0)**2)/ &
          ((gamma*p_infty + (4.0d0/3.0d0)*mu)*(-sigma_0*exp(-(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1)**5*((-sigma_0* &
          exp(-(-10.0d0*t*sqrt((gamma*p_infty + 1.33333333333333d0*mu)/ &
          rho_0) + 10.0d0*x - 5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu &
          ) + 1)**(-2))**(2.0d0/3.0d0)))/(-sigma_0*exp(-(-10.0d0*t*sqrt(( &
          gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - 5.0d0) &
          **2)/(gamma*p_infty + (4.0d0/3.0d0)*mu) + 1))*exp(-(-10.0d0*t* &
          sqrt((gamma*p_infty + 1.33333333333333d0*mu)/rho_0) + 10.0d0*x - &
          5.0d0)**2)/(gamma*p_infty + (4.0d0/3.0d0)*mu)
    
    end function

end module

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info

    use ManufacturedSolution_data

    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh

    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)
    
    ! Create mesh from [0,2*pi)x[0,2*pi)x[0,2*pi) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = one/real(nx-1,rkind)
        dy = dx
        dz = dx

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 - 1 + i - 1, rkind ) * dx
                    y(i,j,k) = real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,rho0,mu,yield,gam,PInf,tau0,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,third,half,one,two,three,pi,four,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,&
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info
    use StiffGasEOS,      only: stiffgas
    use Sep1SolidEOS,     only: sep1solid
    
    use ManufacturedSolution_data

    implicit none
    character(len=*),                                               intent(in)    :: inputfile
    type(decomp_info),                                              intent(in)    :: decomp
    real(rkind),                                                    intent(in)    :: dx,dy,dz
    real(rkind),                                          optional, intent(inout) :: rho0, mu, gam, PInf, tstop, dt, tviz, yield, tau0
    real(rkind), dimension(:,:,:,:),     intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields

    integer :: ioUnit
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: sig11, eps11
    real(rkind) :: t = zero, p_star, rho_star, c_star

    namelist /PROBINPUT/  sigma_0
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    associate( rho => fields(:,:,:,rho_index),   u => fields(:,:,:,  u_index), &
                 v => fields(:,:,:,  v_index),   w => fields(:,:,:,  w_index), &
                 p => fields(:,:,:,  p_index),   T => fields(:,:,:,  T_index), &
                 e => fields(:,:,:,  e_index), g11 => fields(:,:,:,g11_index), &
               g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), & 
               g23 => fields(:,:,:,g23_index), g31 => fields(:,:,:,g31_index), & 
               g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        
        p_star = gam*PInf + (four/three)*mu
        rho_star = rho0
        c_star = sqrt(p_star/rho_star)
        sigma_0 = sigma_0 / p_star
        tstop = tstop * c_star
        dt = dt * c_star
        tviz = tviz * c_star

        dtprob = dt

        print*, "p_star = ", p_star
        print*, "rho_star = ", rho_star
        print*, "c_star = ", c_star
        print*, "sigma_0 = ", sigma_0
        print*, "tstop = ", tstop
        print*, "dt = ", dt
        print*, "tviz = ", tviz

        ! Non-dimensionalize problem parameters
        rho0 = rho0 / rho_star
        mu = mu / p_star
        PInf = PInf / p_star

        gamma = gam
        p_infty = PInf
        rho_0 = rho0

        muL = mu
        lamL = gam * PInf - (two/three)*mu
        cL = sqrt( (lamL + two*muL)/ rho0 )

        sig11 = -sigma_0 * exp( -((x - cL*t - half)/(0.1_rkind))**2 )
        eps11 = sig11 / (lamL + two*muL)

        u   = - (cL * sig11) / (lamL + two*muL)
        v   = zero
        w   = zero

        g11 = one;  g12 = zero; g13 = zero
        g21 = zero; g22 = one;  g23 = zero
        g31 = zero; g32 = zero; g33 = one

        g11 = one / (one + eps11)

        p = PInf*(g11**gam - one)

        ! Get rho compatible with det(g) and rho0
        eps11 = g11*(g22*g33-g23*g32) - g12*(g21*g33-g31*g23) + g13*(g21*g32-g31*g22)
        rho = rho0 * eps11

    end associate

end subroutine

subroutine hook_output(decomp,der,fil,dx,dy,dz,outputdir,mesh,fields,tsim,vizcount,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,eps,half,one,two,pi,eight
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use FiltersMod,       only: filters

    use ManufacturedSolution_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(filters),                   intent(in) :: fil
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    real(rkind), dimension(2),       intent(in) :: x_bc, y_bc, z_bc
    integer                                     :: outputunit=229

    character(len=clen) :: outputfile, str
    integer :: i

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => fields(:,:,:,g11_index), g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), g23 => fields(:,:,:,g23_index), &
               g31 => fields(:,:,:,g31_index), g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        write(str,'(ES7.1E2,A1,I4.4,A1,ES7.1E2)') sigma_0, "_", decomp%ysz(1), "_", dtprob
        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ManufacturedSolution_"//trim(str)//"_", vizcount, ".dat"

        open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        write(outputunit,'(6ES26.16)') tsim, sigma_0, muL, lamL, cL, p_infty
        do i=1,decomp%ysz(1)
            write(outputunit,'(14ES26.16)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
                                           g11(i,1,1), g21(i,1,1), mu(i,1,1), bulk(i,1,1), kap(i,1,1), &
                                           momentum_source(gamma, muL, p_infty, rho_0, sigma_0, tsim, x(i,1,1)), &
                                           energy_source(gamma, muL, p_infty, rho_0, sigma_0, tsim, x(i,1,1)), &
                                           g_source(gamma, muL, p_infty, rho_0, sigma_0, tsim, x(i,1,1)), &
                                           rho_0 * g_source(gamma, muL, p_infty, rho_0, sigma_0, tsim, x(i,1,1))
        end do
        close(outputunit)

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info

    use ManufacturedSolution_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind), dimension(2),       intent(in) :: x_bc, y_bc, z_bc

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
    end associate
end subroutine

subroutine hook_timestep(decomp,mesh,fields,step,tsim)
    use kind_parameters,  only: rkind
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index
    use decomp_2d,        only: decomp_info
    use exits,            only: message
    use reductions,       only: P_MAXVAL

    use ManufacturedSolution_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        

    end associate
end subroutine

subroutine hook_source(decomp,mesh,fields,tsim,rhs,rhsg)
    use kind_parameters,  only: rkind
    use constants,        only: zero
    use SolidGrid,        only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index, &
                                g11_index,g12_index,g13_index,g21_index,g22_index,g23_index,g31_index,g32_index,g33_index
    use decomp_2d,        only: decomp_info

    use ManufacturedSolution_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhsg

    integer :: i,j,k

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
               g11 => fields(:,:,:,g11_index), g12 => fields(:,:,:,g12_index), g13 => fields(:,:,:,g13_index), & 
               g21 => fields(:,:,:,g21_index), g22 => fields(:,:,:,g22_index), g23 => fields(:,:,:,g23_index), &
               g31 => fields(:,:,:,g31_index), g32 => fields(:,:,:,g32_index), g33 => fields(:,:,:,g33_index), & 
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        do k = 1,decomp%ysz(3)
            do j = 1,decomp%ysz(2)
                do i = 1,decomp%ysz(1)
                    rhs (i,j,k,2) = rhs (i,j,k,2) + momentum_source(gamma, muL, p_infty, rho_0, sigma_0, tsim, x(i,j,k))
                    rhs (i,j,k,5) = rhs (i,j,k,5) + energy_source(gamma, muL, p_infty, rho_0, sigma_0, tsim, x(i,j,k))
                    rhsg(i,j,k,1) = rhsg(i,j,k,1) + g_source(gamma, muL, p_infty, rho_0, sigma_0, tsim, x(i,j,k))
                    rhs (i,j,k,1) = rhs (i,j,k,1) + rho_0 * g_source(gamma, muL, p_infty, rho_0, sigma_0, tsim, x(i,j,k))
                end do
            end do
        end do

    end associate

end subroutine

subroutine hook_finalize()

end subroutine
