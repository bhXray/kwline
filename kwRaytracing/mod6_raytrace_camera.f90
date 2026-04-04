! kerrwh v1 21/1/2026
module mod6_raytrace_camera
    use constants
    use blcoordinate, only: metricg, p_total, YNOGK, lambdaq, Pemdisk, Pemdisk_all
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: raytrace_camera, trace_single_ray
    ! module-level parameters (shared by all procedures)
    double precision, parameter :: PTOL_MIN = 1.0d-250
    double precision, parameter :: ESCAPE_FACTOR = 3.0d0

    contains

    subroutine camera_setup(sinobs, muobs, a_spin, lambdaBar, robs, velocity, rhorizon, rthroat)
        implicit none
            double precision, intent(in)  :: sinobs, muobs, a_spin, lambdaBar, robs
            double precision, intent(out) :: velocity(3), rhorizon, rthroat

            double precision :: somiga, expnu, exppsi, expmu1, expmu2

        call metricg(robs, sinobs, muobs, a_spin, lambdaBar, somiga, expnu, exppsi, expmu1, expmu2)

        velocity(1) = 0.0d0
        velocity(2) = 0.0d0
        velocity(3) = exppsi/expnu * ( one/(robs**(three/two) + a_spin) - somiga )

        rthroat = one + lambdaBar**two + sqrt((one + lambdaBar**two)**two - a_spin**two)
        rhorizon = (one + sqrt(one - a_spin**two))
    end subroutine camera_setup


    subroutine ray_init(alpha, beta, sinobs, muobs, a_spin, lambdaBar, robs, scal, velocity, &
                                            f1234, lambda, q, ptotal, debug)
        implicit none
            double precision, intent(in)  :: alpha, beta, sinobs, muobs, a_spin, lambdaBar, robs, scal
            double precision, intent(in)  :: velocity(3)
            double precision, intent(out) :: f1234(4), lambda, q, ptotal
            logical, intent(in), optional :: debug

        logical :: dbg
        double precision :: mu, rin, rout
        double precision :: pem, pem_all
        dbg = .false.
        if (present(debug)) dbg = debug
        if (dbg) write(*,*) "[DBG] ray_init: alpha =", alpha, " beta =", beta
        if (dbg) write(*,*) "[DBG] ray_init: sinobs =", sinobs, " muobs =", muobs
        if (dbg) write(*,*) "[DBG] ray_init: a_spin =", a_spin, " robs =", robs
        if (dbg) write(*,*) "[DBG] ray_init: robs =", robs, " scal =", scal
        if (dbg) write(*,*) "[DBG] ray_init: velocity =", velocity(1), velocity(2), velocity(3)
        if (dbg) write(*,*) " "

        call lambdaq(alpha, beta, robs, sinobs, muobs, a_spin, lambdaBar, scal, velocity, f1234, lambda, q)

        if (dbg) write(*,*) "[DBG] ray_init: f1234", f1234
        if (dbg) write(*,*) "[DBG] ray_init: lambda", lambda
        if (dbg) write(*,*) "[DBG] ray_init: q", q
        if (dbg) write(*,*) " "

        mu = zero
        rin = (one + lambdaBar**two) + sqrt(max(0.0d0, (one + lambdaBar**two)**two - a_spin**two)) + 0.01d0
        rout = 30.0d0
        pem = Pemdisk(f1234, lambda, q, sinobs, muobs, a_spin, lambdaBar, robs, scal, mu, rout, rin)
        pem_all = Pemdisk_all(f1234, lambda, q, sinobs, muobs, a_spin, lambdaBar, robs, scal, mu, rout, rin)
        ! assign ptotal from computed pem (was previously left unset)
        ptotal = -1.0d0 ! pem 
        ! Override ptotal to allow long tracing (past turning point)
        ! ptotal = pem
        if (dbg) then
            write(*,'(A,1X,ES24.16)') "[DBG] pemd =", pem
            write(*,'(A,1X,ES24.16)') "[DBG] pemd_all =", pem_all
        endif

        ! ptotal = pem

        ! fallback if ptotal collapses to denormal
        if (abs(ptotal) < PTOL_MIN) then
            ptotal = -1.0d0   ! flag: use fixed dp stepping
        endif
    end subroutine ray_init

    subroutine trace_ray_core(sinobs, muobs, a_spin, lambdaBar, robs, scal, rhorizon, rthroat, &
                            f1234, lambda, q, ptotal, nsteps, &
                            x, y, z, tarr, sarr, nout, debug)
        implicit none
        integer, intent(in) :: nsteps
        double precision, intent(in) :: sinobs, muobs, a_spin, lambdaBar, robs, scal, rhorizon
        double precision, intent(in) :: rthroat
        double precision, intent(in) :: f1234(4), lambda, q, ptotal
        double precision, intent(out) :: x(nsteps), y(nsteps), z(nsteps)
        double precision, intent(out) :: tarr(nsteps), sarr(nsteps)
        integer, intent(out) :: nout
        logical, intent(in), optional :: debug

        logical :: dbg
        integer :: k
        double precision :: p, p_prev, ra, mua, phya, timea, sigmaa
        double precision :: rr_xy, dp_base, dp_dyn, frac
        ! variables for phi unwrapping
        double precision :: phya_raw_prev, phya_raw, phya_cont, dphi, phia_offset
        ! adaptive refinement stack and temporaries
        integer, parameter :: MAX_STACK = 512
        double precision :: pstack(MAX_STACK)
        integer :: sp
        double precision :: pcur, pmid
        double precision :: xcur, ycur, zcur, ra_local, mua_cur, phraw_local, phcont_local, time_local, sigma_local
        double precision :: dist, dist_thresh, rr_xy_cur
        ! bound orbit detection
        double precision, parameter :: MAX_ORBIT_PHI = 188.5d0  ! 30 full orbits (60*pi)
        integer, parameter :: MAX_BOUND_POINTS = 100000  ! limit points for bound orbits
        double precision :: phi_initial
        logical :: is_bound_orbit
        ! Track radius change to detect true bound orbits vs spiraling rays
        double precision :: r_min_recent
        integer :: points_since_r_check

        dbg = .false.
        if (present(debug)) dbg = debug

        ! nout = 0
        
        dp_base = 200.0d0 / dble(nsteps) * scal   ! fallback base step size
        p_prev = 0.0d0

        nout = 1
        call YNOGK(0.0d0, f1234, lambda, q, sinobs, muobs, a_spin, lambdaBar, robs, scal, &
                    ra, mua, phya, timea, sigmaa)
        ! initialize phi unwrap state
        phia_offset = 0.0d0
        phya_raw_prev = phya
        phya_cont = phya + phia_offset
        phi_initial = phya_cont   ! store initial phi for bound orbit detection
        is_bound_orbit = .false.
        r_min_recent = ra
        points_since_r_check = 0
        rr_xy = sqrt(ra*ra + a_spin*a_spin) * sqrt(max(0.0d0, one - mua*mua))
        x(nout)    =  rr_xy * cos(phya_cont)
        y(nout)    = -rr_xy * sin(phya_cont)
        z(nout)    =  ra * mua
        tarr(nout) =  timea
        sarr(nout) =  sigmaa
        
        if (dbg) then
            write(*,'(A,ES20.10)') "[DBG] rhorizon =", rhorizon
            write(*,'(A,ES20.10)') "[DBG] rthroat  =", rthroat
            write(*,'(A,ES20.10)') "[DBG] robs     =", robs
            ! call YNOGK(0.0d0, f1234, lambda, q, sinobs, muobs, a_spin, lambdaBar, robs, scal, ra, mua, phya, timea, sigmaa)
            write(*,'(A,5ES20.10)') "[DBG] YNOGK(p=0): ra mua phya time sig =", ra, mua, phya, timea, sigmaa
        end if


        if (dbg) write(*,*) "[trace_ray_core] Run Debug"
        do k = 2, nsteps
            if (ptotal > 0.0d0) then
                p = dble(k-1) * ptotal / dble(nsteps-1)
            else
                frac = dble(k-1) / dble(nsteps-1)
                dp_dyn = dp_base / (1.0d0 + 4.0d0*frac)  ! shrink dp as p grows
                p = p_prev + dp_dyn
            endif

            ! We'll refine the interval [p_prev, p] by bisection if the straight
            ! line between last saved point and the new point is too long. This
            ! avoids large jumps when the trajectory curvature is high.
            ! push the target p onto stack
            sp = 1
            pstack(sp) = p

            ! current last-saved coordinates and p_prev are the reference
            xcur = x(nout)
            ycur = y(nout)
            zcur = z(nout)

            do while (sp > 0)
                pcur = pstack(sp)
                sp = sp - 1
                call YNOGK(pcur, f1234, lambda, q, sinobs, muobs, a_spin, lambdaBar, robs, scal, &
                        ra_local, mua_cur, phraw_local, time_local, sigma_local)
                if (.not. ieee_is_finite(ra_local) .or. .not. ieee_is_finite(phraw_local) .or. &
                    .not. ieee_is_finite(time_local) .or. .not. ieee_is_finite(sigma_local)) then
                    sp = 0
                    exit
                endif

                ! unwrap phi consistently using previous raw phi
                dphi = phraw_local - phya_raw_prev
                if (dphi > pi) then
                    phia_offset = phia_offset - twopi
                else if (dphi < -pi) then
                    phia_offset = phia_offset + twopi
                end if
                phcont_local = phraw_local + phia_offset

                rr_xy_cur = sqrt(ra_local*ra_local + a_spin*a_spin) * sqrt(max(0.0d0, one - mua_cur*mua_cur))
                xcur = rr_xy_cur * cos(phcont_local)
                ycur = -rr_xy_cur * sin(phcont_local)
                zcur = ra_local * mua_cur

                ! adapt threshold: absolute minimum plus small fraction of radius
                dist = sqrt((xcur - x(nout))**2 + (ycur - y(nout))**2)
                dist_thresh = max(0.2d0, 0.01d0 * rr_xy_cur)

                if (dist > dist_thresh .and. sp+2 < MAX_STACK .and. dabs(pcur - p_prev) > 1e-12) then
                    ! subdivide: process midpoint first, then the current pcur
                    pmid = 0.5d0*(p_prev + pcur)
                    sp = sp + 1
                    pstack(sp) = pcur
                    sp = sp + 1
                    pstack(sp) = pmid
                else
                    ! accept point
                    phya_raw_prev = phraw_local
                    phya_cont = phcont_local
                    rr_xy = rr_xy_cur
                    if (ptotal <= 0.0d0) p_prev = pcur

                    nout = nout + 1
                    if (nout > nsteps) then
                        nout = nsteps
                        exit
                    endif
                    x(nout)    = xcur
                    y(nout)    = ycur
                    z(nout)    = zcur
                    tarr(nout) = time_local
                    sarr(nout) = sigma_local
                    
                    ! Track minimum radius to detect true bound orbits
                    points_since_r_check = points_since_r_check + 1
                    if (ra_local < r_min_recent - 0.01d0) then
                        r_min_recent = ra_local
                        points_since_r_check = 0
                    endif
                    
                    ! True bound orbit: many points without getting closer to throat
                    ! Only terminate if radius hasn't decreased significantly in many points
                    if (points_since_r_check > 8000 .and. ra_local > 1.05d0*rthroat) then
                        is_bound_orbit = .true.
                        sp = 0  ! clear stack to exit inner while loop
                    endif
                    ! Also limit by total orbit count if radius is stuck at large value
                    if (dabs(phya_cont - phi_initial) > MAX_ORBIT_PHI .and. ra_local > 1.2d0*rthroat) then
                        is_bound_orbit = .true.
                        sp = 0  ! clear stack to exit inner while loop
                    endif
                    if (nout > MAX_BOUND_POINTS) then
                        is_bound_orbit = .true.
                        sp = 0  ! clear stack to exit inner while loop
                    endif
                endif
            end do
            
            ! Exit outer loop if bound orbit was detected
            if (is_bound_orbit) then
                if (dbg) write(*,*) "[DBG] Bound orbit detected, terminating"
                exit
            endif

            ! stop conditions evaluated after points accepted (use last-evaluated ra_local)
            if (ra_local < 1.005d0*rthroat) exit   ! stop very close to throat
            if (ra_local >= robs*ESCAPE_FACTOR) exit

        end do

        if (nout == 0) then
            nout = 1
            x(1) = -100.0d0; y(1) = -100.0d0; z(1) = -100.0d0
            tarr(1) = 0.0d0
            sarr(1) = 0.0d0
        end if
    end subroutine trace_ray_core

    subroutine trace_single_ray(alpha, beta, sinobs, muobs, a_spin, lambdaBar, robs, scal, nsteps, &
                        x, y, z, tarr, sarr, nout, debug, out_file)
        implicit none
            integer, intent(in) :: nsteps
            double precision, intent(in) :: alpha, beta
            double precision, intent(in) :: sinobs, muobs, a_spin, lambdaBar, robs, scal
            logical, intent(in), optional :: debug
            character(len=*), intent(in), optional :: out_file

            double precision, intent(out) :: x(nsteps), y(nsteps), z(nsteps)
            double precision, intent(out) :: tarr(nsteps), sarr(nsteps)
            integer, intent(out) :: nout
            
            logical :: dbg
            integer :: i, unit, ios
            double precision :: velocity(3), rhorizon, rthroat
            double precision :: f1234(4), lambda, q, ptotal

        ! ---------- make output directory (requires Fortran 2008) ----------
        call execute_command_line("mkdir -p output")

        dbg = .false.
        if (present(debug)) dbg = debug

        ! one-time-per-call setup
        call camera_setup(sinobs, muobs, a_spin, lambdaBar, robs, velocity, rhorizon, rthroat)
        call ray_init(alpha, beta, sinobs, muobs, a_spin, lambdaBar, robs, scal, velocity, f1234, lambda, q, ptotal, dbg)

        call trace_ray_core(sinobs, muobs, a_spin, lambdaBar, robs, scal, rhorizon, rthroat, &
                            f1234, lambda, q, ptotal, nsteps, &
                            x, y, z, tarr, sarr, nout, debug)

        ! Smooth out numerical phi kinks from YNOGK Weierstrass function precision
        call smooth_phi_kinks(x, y, z, nout, a_spin)

        ! optional output
        if (present(out_file)) then
            unit = 99
            open(unit=unit, file=trim(out_file), status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write(*,*) "[ERR] Cannot open output file:", trim(out_file)
                stop 2
            endif

            write(unit,'(A)') "# Single ray output"
            write(unit,'(A,1X,ES26.16)') "# alpha =", alpha
            write(unit,'(A,1X,ES26.16)') "# beta  =", beta
            write(unit,'(A,1X,ES26.16)') "# labmdaBar=", lambdaBar
            write(unit,'(A,1X,ES26.16)') "# a_spin=", a_spin
            write(unit,'(A,1X,ES26.16)') "# robs  =", robs
            write(unit,'(A,1X,ES26.16)') "# scal  =", scal
            write(unit,'(A,1X,ES26.16)') "# sinobs=", sinobs
            write(unit,'(A,1X,ES26.16)') "# muobs =", muobs
            write(unit,'(A,1X,ES26.16)') "# ptotal=", ptotal
            write(unit,'(A,1X,ES26.16)') "# rthroat=", rthroat
            write(unit,'(A,1X,ES26.16)') "# rhorizon=", rhorizon
            write(unit,'(A)') "# Columns: x y z timea sigmaa"

            do i = 1, nout
                write(unit,'(5ES26.16)') x(i), y(i), z(i), tarr(i), sarr(i)
            end do

            close(unit)

            if (dbg) write(*,*) "[OK] Wrote:", trim(out_file)
        end if
    end subroutine trace_single_ray



    subroutine smooth_phi_kinks(x, y, z, n, a_spin)
        ! Detect and fix isolated single-point phi kinks from YNOGK numerical
        ! precision issues. Only phi has small discontinuities; r and mu are accurate.
        ! Detection criterion: a kink at point i causes dphi(i-1->i) and dphi(i->i+1)
        ! to have equal-and-opposite errors relative to the smooth trend, while the
        ! sum phi(i+1)-phi(i-1) remains approximately 2*dphi_expected.
        ! This distinguishes single-point errors from genuine trajectory oscillations.
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: a_spin
        double precision, intent(inout) :: x(n), y(n), z(n)

        integer :: i, nfixed
        double precision :: phi(n), rr_xy(n)
        double precision :: dphi_left, dphi_right, dphi_in_l, dphi_in_r
        double precision :: asymmetry, scale, phi_corrected
        ! Minimum absolute phi error (radians) to trigger correction (~0.006°)
        double precision, parameter :: MIN_PHI_ERR = 1.0d-4

        if (n < 7) return

        ! Extract rr_xy and phi from Cartesian coordinates
        do i = 1, n
            rr_xy(i) = sqrt(x(i)**2 + y(i)**2)
            if (rr_xy(i) > 1.0d-15) then
                phi(i) = atan2(-y(i), x(i))
            else
                phi(i) = 0.0d0
            end if
        end do

        ! Unwrap phi
        do i = 2, n
            do while (phi(i) - phi(i-1) > PI)
                phi(i) = phi(i) - twopi
            end do
            do while (phi(i) - phi(i-1) < -PI)
                phi(i) = phi(i) + twopi
            end do
        end do

        ! Detect and fix single-point kinks using 4-point stencil
        ! For point i, compare the two outer dphi segments (i-2->i-1, i+1->i+2)
        ! with the two inner segments (i-1->i, i->i+1).
        ! A single-point error at i causes:
        !   dphi_in_l = dphi_expected + epsilon
        !   dphi_in_r = dphi_expected - epsilon
        ! So the inner segments are anti-correlated while the outers are smooth.
        ! Key test: the asymmetry |dphi_in_l - dphi_in_r| is much larger than
        !           |dphi_left - dphi_right| (the outer segment difference).
        nfixed = 0
        do i = 3, n-2
            dphi_left  = phi(i-1) - phi(i-2)   ! outer left segment
            dphi_right = phi(i+2) - phi(i+1)   ! outer right segment
            dphi_in_l  = phi(i)   - phi(i-1)   ! inner left (through suspect point)
            dphi_in_r  = phi(i+1) - phi(i)     ! inner right (through suspect point)

            ! The asymmetry of inner segments
            asymmetry = abs(dphi_in_l - dphi_in_r)
            ! The scale: typical dphi magnitude from outer segments
            scale = (abs(dphi_left) + abs(dphi_right)) / 2.0d0

            ! Skip if scale is too small
            if (scale < 1.0d-15) cycle

            ! A kink means inner segments are much more different from each other
            ! than the outer segments are from each other.
            ! Also require that the phi error is above the absolute minimum threshold.
            ! The phi error is approximately asymmetry/2.
            if (asymmetry/2.0d0 > MIN_PHI_ERR .and. &
                asymmetry > 3.0d0 * abs(dphi_left - dphi_right)) then
                ! Correct phi(i) by linear interpolation between i-1 and i+1
                ! weighted by the outer segment dphi average
                phi_corrected = phi(i-1) + (dphi_left + dphi_right) / 2.0d0
                phi(i) = phi_corrected
                nfixed = nfixed + 1
            end if
        end do

        if (nfixed == 0) return

        ! Reconstruct x, y from corrected phi, preserving the ORIGINAL rr_xy
        do i = 1, n
            if (rr_xy(i) > 1.0d-15) then
                x(i) =  rr_xy(i) * cos(phi(i))
                y(i) = -rr_xy(i) * sin(phi(i))
            end if
        end do

    end subroutine smooth_phi_kinks


    subroutine raytrace_camera(sinobs, muobs, a_spin, lambdaBar, robs, scal, nx, ny, alpha_max, beta_max, debug)
        implicit none

        integer, intent(in) :: nx, ny
        double precision, intent(in) :: sinobs, muobs, a_spin, lambdaBar, robs, scal
        double precision, intent(in) :: alpha_max, beta_max
        logical, intent(in), optional :: debug

        integer :: ix, iy, i, unit, ios
        integer :: ray_id
        double precision :: alpha, beta, dalpha, dbeta
        

        ! ray integration output
        integer, parameter :: nsteps_cam = 10001
        double precision :: x(nsteps_cam), y(nsteps_cam), z(nsteps_cam)
        double precision :: tarr(nsteps_cam), sarr(nsteps_cam)
        integer :: nout

        ! cached observer setup
        double precision :: velocity(3), rhorizon, rthroat

        ! per-ray constants
        double precision :: f1234(4), lambda, q, ptotal

        character(len=256) :: fname
        logical :: dbg
        dbg = .false.
        if (present(debug)) dbg = debug

        ! ---------- make output directory (requires Fortran 2008) ----------
        call execute_command_line("mkdir -p output")

        ! ---------- setup once ----------
        call camera_setup(sinobs, muobs, a_spin, lambdaBar, robs, velocity, rhorizon, rthroat)

        dalpha = (2.0d0*alpha_max) / dble(nx-1)
        dbeta  = (2.0d0*beta_max ) / dble(ny-1)

        ray_id = 0

        do iy = 0, ny-1
            beta = -beta_max + dble(iy)*dbeta
            do ix = 0, nx-1
                alpha = -alpha_max + dble(ix)*dalpha
                ray_id = ray_id + 1

                ! init constants for this ray
                call ray_init(alpha, beta, sinobs, muobs, a_spin, lambdaBar, robs, scal, velocity, &
                            f1234, lambda, q, ptotal, dbg)

                ! trace (core returns x,y,z,timea,sigmaa)
                call trace_ray_core(sinobs, muobs, a_spin, lambdaBar, robs, scal, rhorizon, rthroat, &
                                    f1234, lambda, q, ptotal, nsteps_cam, &
                                    x, y, z, tarr, sarr, nout, dbg)
                
                ! filename: output/ray{n}.txt
                write(fname,'("output/ray",I0,".txt")') ray_id

                unit = 20
                open(unit=unit, file=trim(fname), status="replace", action="write", iostat=ios)
                if (ios /= 0) then
                    write(*,*) "[ERR] Cannot open file:", trim(fname)
                    stop 2
                endif

                ! header (similar to single-ray)
                write(unit,'(A)') "# Camera ray output"
                write(unit,'(A,1X,I0)')      "# ray_id =", ray_id
                write(unit,'(A,1X,I0,1X,I0)')"# iy ix  =", iy, ix
                write(unit,'(A,1X,ES26.16)') "# alpha  =", alpha
                write(unit,'(A,1X,ES26.16)') "# beta   =", beta
                write(unit,'(A,1X,ES26.16)') "# a_spin =", a_spin
                write(unit,'(A,1X,ES26.16)') "# robs   =", robs
                write(unit,'(A,1X,ES26.16)') "# scal   =", scal
                write(unit,'(A,1X,ES26.16)') "# sinobs =", sinobs
                write(unit,'(A,1X,ES26.16)') "# muobs  =", muobs
                write(unit,'(A,1X,ES26.16)') "# ptotal =", ptotal
                write(unit,'(A)') "# Columns: x y z timea sigmaa"

                do i = 1, nout
                    write(unit,'(5ES26.16)') x(i), y(i), z(i), tarr(i), sarr(i)
                end do

                close(unit)

            end do
        end do

        write(*,*) "[OK] Wrote ", ray_id, " ray files into output/"

    end subroutine raytrace_camera

    
end module mod6_raytrace_camera
