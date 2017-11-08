module ElasticEOSMod

    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit

    implicit none

    type, abstract :: elasticeos

    contains

        procedure(get_finger_interface),    deferred :: get_finger
        procedure(get_devstress_interface), deferred :: get_devstress
        procedure(get_eelastic_interface),  deferred :: get_eelastic
        procedure(get_sos_interface),       deferred :: get_sos

    end type

    abstract interface

        subroutine get_finger_interface(this,g,finger,fingersq,trG,trG2,detG,use_gTg)
            import :: elasticeos
            import :: rkind
            class(elasticeos), intent(in) :: this
            real(rkind), dimension(:,:,:,:), intent(in)  :: g
            real(rkind), dimension(:,:,:,:), intent(out) :: finger
            real(rkind), dimension(:,:,:,:), intent(out) :: fingersq
            real(rkind), dimension(:,:,:),   intent(out) :: trG, trG2, detG
            logical,                         intent(in), optional :: use_gTg
        end subroutine

        pure subroutine get_devstress_interface(this,finger,fingersq,trG,trG2,detG,devstress,rho0mix,mumix)
            import :: elasticeos
            import :: rkind
            class(elasticeos), intent(in) :: this
            real(rkind), dimension(:,:,:,:), intent(in)  :: finger
            real(rkind), dimension(:,:,:,:), intent(in)  :: fingersq
            real(rkind), dimension(:,:,:),   intent(in)  :: trG, trG2, detG
            real(rkind), dimension(:,:,:,:), intent(out) :: devstress
            real(rkind), dimension(:,:,:),   intent(in),optional  :: rho0mix, mumix
        end subroutine

        subroutine get_eelastic_interface(this,trG,trG2,detG,eelastic,rho0mix,mumix)
            import :: elasticeos
            import :: rkind
            class(elasticeos), intent(in) :: this
            real(rkind), dimension(:,:,:), intent(in)  :: trG,trG2,detG
            real(rkind), dimension(:,:,:), intent(out) :: eelastic
            real(rkind), dimension(:,:,:), intent(in),optional  :: rho0mix, mumix
        end subroutine

        pure subroutine get_sos_interface(this,rhom,sos)
            import :: elasticeos
            import :: rkind
            class(elasticeos), intent(in) :: this
            real(rkind), dimension(:,:,:), intent(in) :: rhom
            real(rkind), dimension(:,:,:), intent(inout) :: sos
        end subroutine

    end interface

end module
