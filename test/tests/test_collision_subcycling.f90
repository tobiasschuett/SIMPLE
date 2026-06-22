! Tests for the collision-step sub-cycling in collide() (src/simple_main.f90).
!
! The explicit collision step in stost becomes unstable for slow (near-thermal)
! markers, where the pitch-diffusion coefficient dhh ~ 1/p**2 grows large and a
! single Langevin increment overshoots its |lambda| <= 1 bound ("Error in stost").
! collide() sub-cycles the step, re-evaluating the coefficients each sub-step.
!
! Tests 1 and 2 call the real collide() and observe what it did through the
! diag_counters event counters it bumps (EVT_COLLIDE_SUBSTEP, etc.), the same
! way test_axis_crossing reads EVT_R_NEGATIVE.
!
! Test 1: at the standard background density/temperature a fresh marker needs only
!         a single sub-step, so collide() behaves exactly as before.
! Test 2: at the high background density that triggered the failures, a single
!         step trips the pitch overshoot but collide() sub-cycles and avoids it.
! Test 3: the safety ordering - as a marker slows, sub-cycling switches on (nsub > 1)
!         at a larger p than the p where a single step first overshoots, so the step
!         is already being resolved before the instability is reached.
program test_collision_subcycling
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use simple_profiles
    use collis_alp
    use simple_main, only: collide, collision_nsub_max
    use diag_counters, only: diag_counters_init, diag_counters_total, &
                             EVT_STOST_PITCH_OVERSHOOT, EVT_COLLIDE_SUBSTEP
    use collision_test_utils, only: set_flat_two_power, seed_rng, reference_dtaumin
    implicit none

    logical :: all_passed
    integer :: n_failed

    all_passed = .true.
    n_failed = 0

    call test_standard_density_single_substep(all_passed, n_failed)
    call test_substep_prevents_stost_error(all_passed, n_failed)
    call test_subcycling_engages_before_instability(all_passed, n_failed)

    if (all_passed) then
        print *, 'All collision tests passed'
    else
        print *, 'FAILED: ', n_failed, ' collision tests failed'
        error stop 1
    end if

contains

    ! Standard D-T reactor background (the params.f90 defaults): a fresh, fast
    ! marker (p = 1) is weakly collisional over one step, so collide() takes a
    ! single sub-step - i.e. the original single stost call, behavior unchanged.
    subroutine test_standard_density_single_substep(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail
        real(dp) :: am1, am2, Z1, Z2, ealpha, v0
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe
        real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: z(5), dtaumin
        integer(int64) :: nsub

        print *, 'Testing standard density: collide() takes a single sub-step...'

        am1 = 2.0d0; am2 = 3.0d0; Z1 = 1.0d0; Z2 = 1.0d0; ealpha = 3.5d6
        ni1_scale = 0.5d20; ni2_scale = 0.5d20
        Te_scale = 1.0d4; Ti1_scale = 1.0d4; Ti2_scale = 1.0d4
        densi1 = ni1_scale*1.0d-6; densi2 = ni2_scale*1.0d-6
        tempi1 = Ti1_scale; tempi2 = Ti2_scale; tempe = Te_scale

        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
                          ealpha, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        call set_flat_two_power()
        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)
        call seed_rng()

        ! Fresh birth-energy marker (p = 1) and a representative collision step.
        ! Set the cap the way init_collisions does; diag_counters_init() both
        ! allocates the counters and zeroes them, so the totals below count exactly
        ! this collide() call.
        dtaumin = reference_dtaumin()
        collision_nsub_max = collision_substep_cap(dtaumin)
        call diag_counters_init()
        z = fresh_marker(1.0d0)
        call collide(z, dtaumin)
        nsub = diag_counters_total(EVT_COLLIDE_SUBSTEP)

        if (nsub == 1) then
            print *, '  PASS: nsub = 1, collide reduces to a single stost call'
        else
            print *, '  FAIL: nsub =', nsub, ' (expected 1; behavior would differ from before)'
            passed = .false.; nfail = nfail + 1
        end if
    end subroutine test_standard_density_single_substep

    ! Regression: at the reported high-density parameters a single big collision
    ! step trips the pitch overshoot (ierr=2); collide() sub-cycles and removes it.
    subroutine test_substep_prevents_stost_error(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail
        real(dp) :: am1, am2, Z1, Z2, ealpha, v0
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe
        real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: dpp, dhh, fpeff, p, p_thermal, dtauc, z(5)
        integer, parameter :: max_cool = 10000   ! bound the cooling search (no infinite loop)
        integer :: ierr, icool
        integer(int64) :: nsub, novershoot
        logical :: single_overshoot

        print *, 'Testing sub-cycling prevents the stost pitch blow-up...'

        ! High density, hot background as in the reported failing run.
        am1 = 2.0d0; am2 = 3.0d0; Z1 = 1.0d0; Z2 = 1.0d0; ealpha = 3.5d6
        ni1_scale = 1.0d21; ni2_scale = 1.0d21
        Te_scale = 1.0d4; Ti1_scale = 1.0d4; Ti2_scale = 1.0d4
        densi1 = ni1_scale*1.0d-6; densi2 = ni2_scale*1.0d-6
        tempi1 = Ti1_scale; tempi2 = Ti2_scale; tempe = Te_scale

        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
                          ealpha, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        call set_flat_two_power()
        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)
        call seed_rng()

        ! Realistic collision step: dtaumin from npoiper2 and rbig (reference_dtaumin).
        dtauc = reference_dtaumin()

        ! Locate a marker that triggers the blow-up. At the realistic dtaumin only
        ! cold, sub-thermal markers overshoot in a single step; a real run reaches
        ! them when the explicit energy step drives p far below thermal. Cool p from
        ! thermal until a single step would overshoot (sqrt(2*dhh*dtauc) > sqrt(2)).
        p_thermal = 1.0d0/maxval(velrat(1:2))
        p = p_thermal
        call coleff(p, dpp, dhh, fpeff)
        do icool = 1, max_cool
            if (2.0d0*dhh*dtauc > 2.0d0) exit
            p = 0.95d0*p
            call coleff(p, dpp, dhh, fpeff)
        end do
        if (2.0d0*dhh*dtauc <= 2.0d0) then
            print *, '  FAIL: cooling did not reach the overshoot threshold in', max_cool, 'steps'
            passed = .false.; nfail = nfail + 1
            return
        end if
        print *, '  single-step overshoot requires p =', p, &
                 ', p/p_thermal =', p/p_thermal
        collision_nsub_max = collision_substep_cap(dtauc)

        ! A: one un-sub-cycled step at this marker (lambda = 0, maximal pitch
        ! diffusion) blows up - the pitch overshoot, ierr = 2. (mod(ierr,10)==2
        ! catches it with or without the +10 energy-clamp flag.)
        z = fresh_marker(p)
        call stost(z, dtauc, 1, ierr)
        single_overshoot = (mod(ierr, 10) == 2)

        ! B: the real collide() on the same marker and step sub-cycles and prevents
        ! it (nsub > 1, no overshoot). Counters are zeroed first to count only this call.
        call diag_counters_init()
        z = fresh_marker(p)
        call collide(z, dtauc)
        nsub = diag_counters_total(EVT_COLLIDE_SUBSTEP)
        novershoot = diag_counters_total(EVT_STOST_PITCH_OVERSHOOT)

        if (single_overshoot .and. nsub > 1) then
            print *, '  PASS: single step overshoots; collide sub-cycles, nsub =', nsub
        else
            print *, '  FAIL: did not reproduce setup; single_overshoot =', single_overshoot, &
                     ' nsub =', nsub
            passed = .false.; nfail = nfail + 1
        end if

        if (novershoot == 0) then
            print *, '  PASS: collide produced no pitch overshoot'
        else
            print *, '  FAIL: collide still overshot', novershoot, 'times'
            passed = .false.; nfail = nfail + 1
        end if
    end subroutine test_substep_prevents_stost_error

    ! Safety ordering: as a marker slows, sub-cycling must switch on (nsub > 1) at a
    ! larger p than the p where a single step first overshoots, so the step is already
    ! being resolved before the instability is reached. This holds because sub-cycling
    ! engages at 2*dhh*dtaumin = collision_eps**2 while the overshoot is at
    ! 2*dhh*dtaumin = 1, and collision_eps < 1. Guards against eps being set too large.
    subroutine test_subcycling_engages_before_instability(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail
        real(dp) :: am1, am2, Z1, Z2, ealpha, v0
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe
        real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: dpp, dhh, fpeff, p, p_thermal, p_subcycle_on, p_unstable, dtauc
        integer, parameter :: max_cool = 10000   ! bound the cooling searches (no infinite loop)
        integer :: icool

        print *, 'Testing sub-cycling engages before the single-step instability...'

        am1 = 2.0d0; am2 = 3.0d0; Z1 = 1.0d0; Z2 = 1.0d0; ealpha = 3.5d6
        ni1_scale = 1.0d21; ni2_scale = 1.0d21
        Te_scale = 1.0d4; Ti1_scale = 1.0d4; Ti2_scale = 1.0d4
        densi1 = ni1_scale*1.0d-6; densi2 = ni2_scale*1.0d-6
        tempi1 = Ti1_scale; tempi2 = Ti2_scale; tempe = Te_scale
        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
                          ealpha, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        call set_flat_two_power()
        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)

        dtauc = reference_dtaumin()
        p_thermal = 1.0d0/maxval(velrat(1:2))

        ! Cool from thermal to where sub-cycling first engages (nsub > 1).
        p = p_thermal
        call coleff(p, dpp, dhh, fpeff)
        do icool = 1, max_cool
            if (collision_substeps(dpp, dhh, fpeff, p, dtauc) > 1) exit
            p = 0.95d0*p
            call coleff(p, dpp, dhh, fpeff)
        end do
        if (collision_substeps(dpp, dhh, fpeff, p, dtauc) == 1) then
            print *, '  FAIL: sub-cycling never engaged in', max_cool, 'cooling steps'
            passed = .false.; nfail = nfail + 1
            return
        end if
        p_subcycle_on = p

        ! Cool further to where a single step first overshoots (sqrt(2*dhh*dtauc) > 1).
        do icool = 1, max_cool
            if (2.0d0*dhh*dtauc > 1.0d0) exit
            p = 0.95d0*p
            call coleff(p, dpp, dhh, fpeff)
        end do
        if (2.0d0*dhh*dtauc <= 1.0d0) then
            print *, '  FAIL: single-step instability not reached in', max_cool, 'cooling steps'
            passed = .false.; nfail = nfail + 1
            return
        end if
        p_unstable = p

        print *, '  p(sub-cycling on) =', p_subcycle_on, &
                 ', p(single-step unstable) =', p_unstable

        if (p_subcycle_on > p_unstable) then
            print *, '  PASS: sub-cycling engages above the instability, ratio', &
                     p_subcycle_on/p_unstable
        else
            print *, '  FAIL: sub-cycling engages at p =', p_subcycle_on, &
                     ' not above instability p =', p_unstable
            passed = .false.; nfail = nfail + 1
        end if
    end subroutine test_subcycling_engages_before_instability

    ! Fresh marker at mid-radius (s = 0.5), normalized speed p, lambda = 0
    ! (maximal pitch diffusion).
    function fresh_marker(p) result(z)
        real(dp), intent(in) :: p
        real(dp) :: z(5)
        z = 0.0d0
        z(1) = 0.5d0
        z(4) = p
    end function fresh_marker

end program test_collision_subcycling
