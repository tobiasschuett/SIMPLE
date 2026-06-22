! Shared helpers for the collision tests (test_profiles, test_collision_subcycling).
module collision_test_utils
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple_profiles
    implicit none

contains

    ! Flat (constant) two-power profiles, so the splined collision grid reproduces
    ! the scalar coefficients. Amplitudes come from the *_scale module variables.
    subroutine set_flat_two_power()
        profile_type = "two_power"
        active_profile = TWO_POWER
        Te_p1 = 1.0d0; Te_p2 = 0.0d0
        Ti1_p1 = 1.0d0; Ti1_p2 = 0.0d0
        Ti2_p1 = 1.0d0; Ti2_p2 = 0.0d0
        ni1_p1 = 1.0d0; ni1_p2 = 0.0d0
        ni2_p1 = 1.0d0; ni2_p2 = 0.0d0
    end subroutine set_flat_two_power

    ! Representative orbit microstep dtaumin for the tests, from the same formula
    ! params_init uses: dtaumin = 2*pi*rbig/npoiper2 (see src/params.f90, params_init).
    ! Test values: npoiper2 = 256, rbig = 659 cm (~6.6 m, reactor-scale major radius);
    ! gives dtaumin ~ 16.2, matching a real run.
    function reference_dtaumin() result(dtaumin)
        real(dp) :: dtaumin
        integer, parameter :: npoiper2 = 256
        real(dp), parameter :: rbig = 659.0d0
        real(dp) :: pi
        pi = 4.0d0*atan(1.0d0)
        dtaumin = 2.0d0*pi*rbig/dble(npoiper2)
    end function reference_dtaumin

    ! Deterministic RNG seed so collisional tests are reproducible.
    subroutine seed_rng()
        integer :: n, i
        integer, allocatable :: seed(:)
        call random_seed(size=n)
        allocate (seed(n))
        do i = 1, n
            seed(i) = 4242 + i*7
        end do
        call random_seed(put=seed)
    end subroutine seed_rng

end module collision_test_utils
