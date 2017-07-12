module ElasticEOSMod

    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit

    type, abstract :: elasticeos

    contains

        procedure(get_finger_interface),    deferred :: get_finger
        procedure(get_devstress_interface), deferred :: get_devstress
        procedure(get_eelastic_interface),  deferred :: get_eelastic
        procedure(get_sos_interface),       deferred :: get_sos

    end type

    abstract interface

        pure subroutine get_finger_interface(this,g,finger,fingersq,trG,trG2,detG)
            import :: elasticeos
            import :: rkind
            class(elasticeos), intent(in) :: this
            real(rkind), dimension(:,:,:,:), intent(in)  :: g
            real(rkind), dimension(:,:,:,:), intent(out) :: finger
            real(rkind), dimension(:,:,:,:), intent(out), optional :: fingersq
            real(rkind), dimension(:,:,:),   intent(out), optional :: trG, trG2, detG
        end subroutine

        pure subroutine get_devstress_interface(this,finger,fingersq,trG,trG2,detG,devstress)
            import :: elasticeos
            import :: rkind
            class(elasticeos), intent(in) :: this
            real(rkind), dimension(:,:,:,:), intent(in)  :: finger
            real(rkind), dimension(:,:,:,:), intent(in)  :: fingersq
            real(rkind), dimension(:,:,:),   intent(in)  :: trG, trG2, detG
            real(rkind), dimension(:,:,:,:), intent(out) :: devstress
        end subroutine

        pure subroutine get_eelastic_interface(this,rho0,trG,trG2,detG,eelastic)
            import :: elasticeos
            import :: rkind
            class(elasticeos), intent(in) :: this
            real(rkind),                   intent(in)  :: rho0
            real(rkind), dimension(:,:,:), intent(in)  :: trG,trG2,detG
            real(rkind), dimension(:,:,:), intent(out) :: eelastic
        end subroutine

        pure subroutine get_sos_interface(this,rho0,sos)
            import :: elasticeos
            import :: rkind
            class(elasticeos), intent(in) :: this
            real(rkind),                   intent(in)  :: rho0
            real(rkind), dimension(:,:,:), intent(inout) :: sos
        end subroutine

    end interface

end module
