program trace_one_ray
    use mod6_raytrace_camera, only: trace_single_ray
    use constants
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: nsteps = 1000000
    double precision :: alpha, beta
    double precision :: sinobs, muobs
    double precision :: a_spin, lambdaBar, robs, scal, cobs
    integer :: ios

    double precision :: x(nsteps), y(nsteps), z(nsteps), tarr(nsteps), sarr(nsteps)
    integer :: nout

    ! integer :: i, unit

    read(unit=*, fmt=*, iostat=ios) alpha, beta, a_spin, lambdaBar, cobs, robs, scal
    if (ios /= 0) then
        write(*,*) "Usage:"
        write(*,*) "  echo 'alpha beta a_spin lambdaBar cobs_deg robs scal' | ./trace_one_ray"
        stop 1
    endif

    muobs  = cos(cobs * dtor)
    sinobs = sqrt(max(0.0d0, 1.0d0 - muobs*muobs))

    write(*,*) "[INFO] Tracing single ray with parameters:"
    write(*,*) "       alpha  =", alpha
    write(*,*) "       beta   =", beta
    write(*,*) "       a_spin =", a_spin
    write(*,*) "       lambdaBar =", lambdaBar
    write(*,*) "       cobs   =", cobs
    write(*,*) "       robs   =", robs
    write(*,*) "       scal   =", scal

    call trace_single_ray(alpha, beta, sinobs, muobs, a_spin, lambdaBar, robs, scal, nsteps, &
               x, y, z, tarr, sarr, nout, debug=.true., out_file="output/ray1.txt")
    ! call trace_single_ray(alpha, beta, sinobs, muobs, a_spin, robs, scal, nsteps, &
    !            x, y, z, tarr, sarr, nout, debug=.true.)
    write(*,*) "[OK] Ray traced. Points =", nout


end program trace_one_ray
