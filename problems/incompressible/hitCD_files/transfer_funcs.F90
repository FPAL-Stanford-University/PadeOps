pure elemental function GetFilterTransferFunction(kin,dx) result(T)
    use constants, only: two
    use cf90stuff, only: alpha90, beta90, a90, b90, c90, d90, e90
    real(rkind), intent(in) :: kin,dx
    real(rkind) :: k
    real(rkind) :: T

    k = kin*dx
    T = (a90 + two* b90*COS(k) + two*c90*COS(two*k) +two*d90*COS(3._rkind*k) + two*e90*COS(4._rkind*k) ) &
      / (1._rkind + two*alpha90*COS(k) + two*beta90*COS(two*k) )

end function

pure elemental function GetCD10ModWaveNum(kin,dx) result(kp)
    use constants, only: two
    use cd10stuff, only: alpha10d1, beta10d1, a10d1, b10d1, c10d1
    real(rkind), intent(in) :: kin,dx
    real(rkind) :: k
    real(rkind) :: kp

    k = kin*dx
    kp = ( two*a10d1*sin(k) + two*b10d1*sin(two*k) + two*c10d1*sin(3._rkind*k) ) / (1._rkind + two*alpha10d1*cos(k) + two*beta10d1*cos(two*k))

end function

pure elemental function GetCD06ModWaveNum(kin,dx) result(kp)
    use constants, only: two
    use cd06stuff, only: alpha06d1, a06d1, b06d1
    real(rkind), intent(in) :: kin,dx
    real(rkind) :: k
    real(rkind) :: kp

    k = kin*dx
    kp = ( two*a06d1*sin(k) + two*b06d1*sin(two*k) ) / (1._rkind + two*alpha06d1*cos(k))

end function
