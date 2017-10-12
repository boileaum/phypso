subroutine timeloop(nmax, tmax, xm, wn)

    implicit none

    integer, intent(in) :: nmax
    double precision, intent(in) :: tmax
    double precision, dimension(0:nmax+1), intent(inout) :: xm
    double precision, dimension(0:nmax+1), intent(inout) :: wn

    double precision, dimension(0:nmax+1) :: wnp1
    double precision, parameter :: xmin = -1.d0
    double precision, parameter :: xmax = 2.d0

    integer :: i
    double precision :: t, cfl, dt, dx, w, vmax

    dx = (xmax - xmin)/nmax
    cfl = 0.8d0

    do i = 0, nmax+1
        xm(i) = xmin + (i - 0.5d0)*dx
        call sol_exact(xm(i), 0.d0, w)
        wn(i) = w
        wnp1(i) = wn(i)
    enddo
    t = 0.d0
    do while(t < tmax)
        vmax = maxval(abs(wn))
        dt = cfl*dx/vmax
        do i = 1, nmax
            call riemann(wn(i), wn(i+1), 0.d0, w)  ! Right flux
            wnp1(i) = wnp1(i) - dt/dx*w*w*0.5d0
            call riemann(wn(i-1), wn(i), 0.d0, w)  ! Left flux
            wnp1(i) = wnp1(i) + dt/dx*w*w*0.5d0
        enddo
        wn(:) = wnp1(:)  ! Update
        t = t + dt
        !write(*,*) 't =', t
    enddo

end subroutine


subroutine riemann(wl, wr, xi, w)

    implicit none
    double precision, intent(in) :: xi, wl, wr
    double precision, intent(out) :: w

    double precision :: sigma

    if (wl > wr) then
        ! shock case
        sigma = 0.5d0*(wl + wr)  ! shock velocity
        if (xi < sigma) w = wl
        if (xi >= sigma) w = wr
    else
        ! rarefaction case
        if (xi <= wl) w = wl
        if (xi >= wr) w = wr
        if (xi > wl .and. xi < wr) w = xi
    endif

end subroutine


subroutine sol_exact(x, t, w)

    implicit none
    double precision, intent(in) :: x, t
    double precision, intent(out) :: w

    if (t <= 1.d0) then
        if (x <= t) w = 1.d0
        if (x >= 1.d0) w = 0.d0
        if (x > t .and. x < 1) w = (1.d0 - x)/(1.d0 - t)
    else
        if ((x - 1.d0) <= 0.5d0*(t - 1.d0)) then
            w = 1.d0
        else
            w = 0.d0
        endif
    endif

end subroutine
