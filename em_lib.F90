module em_lib
  use props
  implicit none

  ! MPI Variables:
  integer :: offx, offy

  ! Constants
  real(8), parameter:: a = 0.99 / sqrt(2.0)

  integer, allocatable    :: rij(:,:), rij_em(:,:)
  real(8), allocatable    :: Hz(:,:), Ey(:,:), Ex(:,:), &
                             Jx(:,:), Jy(:,:), wp(:,:), eps(:,:), &
                             lbc(:,:,:), rbc(:,:,:), tbc(:,:,:), bbc(:,:,:), &
                             flx(:), frq(:), x(:), y(:), nu(:,:), &
                             temp(:), ne_avg(:), nte_avg(:), &
                             Em(:), Em_full(:), r_full(:), r_em(:)
  complex(8), allocatable :: a_ey(:), a_hz(:)

  integer :: nx, ny, nt, nf, fskip
  real(8) :: dx, dt, np, wp2

  public  :: flx, frq, em_init, em_run, Em
  private :: px, py, nx, ny, nt, nf, fskip, np, pi, eps0, mu0, c0, e, kb, me, a, Hz, Ey, Ex, &
             Jx, Jy, dx, dt, a_ey, a_hz, wp, eps, lbc, rbc, tbc, bbc, nu, em_step, x, y, rij, &
             get_nu, wp2, Em_full, r_em, r_full, rij_em

contains

! *** Initialization ***
  subroutine em_init(g, nf0, fskip0)
    type(grid), intent(in) :: g
    integer, intent(in) :: nf0, fskip0
    real(8) :: l = 15.8d-3 * 2d0, xtemp, r, min_val, max_val
    real(8), allocatable:: r_temp(:), r_uniq(:)
    integer :: i, j, k, i1=1, i2=1, j1=1

    nx = 2 * g%ny / px
    ny = g%ny / py

    offx = rx * nx
    offy = ry * ny

    allocate(Hz(nx+2, ny+2), Ex(nx+2, ny+2), Ey(nx+2, ny+2), &
             Jx(nx+2, ny+1), Jy(nx+1, ny+2), wp(nx+2,ny+2), eps(nx+2,ny+2), &
             lbc(3,ny+2,2), rbc(3,ny+2,2), tbc(nx+2,3,2), bbc(nx+2,3,2), &
             x(nx+2), y(ny+2), rij(nx+2, ny+2), nu(nx+2, ny+2), &
             ne_avg(g%ny), nte_avg(g%ny), temp(g%by), &
             Em(g%by), r_full(g%ny), r_temp(nx*ny), r_uniq(nx*ny), &
             rij_em(nx,ny))

    Em = 0
    Hz = 0
    Ex = 0
    Ey = 0
    Jx = 0
    Jy = 0
    lbc = 0
    rbc = 0
    tbc = 0
    bbc = 0

    wp = 0
    eps = 1
    nu = 0

    dx = l / 2d0 / float(py*ny+1)
    dt = dx / c0 / sqrt(2.0)
    np  = 1.0 / (13.49e9 * dt)
    nt  = 20 * int(np)
    wp2 = wp02 * dt**2 / t0**2

    nf = nf0
    fskip = fskip0
    if (nf > 0) then
      allocate(frq(nf), flx(nf), a_ey(nf), a_hz(nf))
      frq = 0
      flx = 0
    end if

    do i = 1, nx+2
      x(i) = dx * (offx + i - 1)
    end do

    do j = 1, ny+2
      y(j) = dx * (offy + j - 1)
    end do

    xtemp = x(nx+2)
    call MPI_Bcast(xtemp, 1, MPI_Real8, nproc-1, comm, ierr)
    x = x - xtemp / 2d0
    x = x / x0

    xtemp = y(ny+2)
    call MPI_Bcast(xtemp, 1, MPI_Real8, nproc-1, comm, ierr)
    y = y - xtemp / 2d0
    y = y / x0

    call MPI_AllGather(g%r, g%by, etype, r_full, g%by, etype, rxComm, ierr)

    rij = -1
    do j = 1, ny+2
      do i = 1, nx+2
        r = sqrt(x(i)**2 + y(j)**2)
        do k = 1, g%ny
          if ((r_full(k) < r) .and. (r < g%w)) then
            rij(i,j) = max(min(k,g%by), 2)
          else
            exit
          end if
        end do
      end do
    end do

    r_temp = 0
    do j = 1, ny
      do i = 1, nx
        r_temp(i+(j-1)*nx) = sqrt((dx/x0*(i-0.5))**2 + (dx/x0*(j-0.5))**2)
      end do
    end do

    min_val = minval(r_temp)-1
    max_val = maxval(r_temp)
    i = 0
    do while (min_val<g%w)
        i = i+1
        min_val = minval(r_temp, mask=r_temp>(min_val+0.1))
        r_uniq(i) = min_val
    end do
    allocate(r_em(i-1), source=r_uniq(1:i-1))
    allocate(Em_full(i-1))
    Em_full = 0
    deallocate(r_uniq, r_temp)

    rij_em = -1
    do j = 1, ny
      do i = 1, nx
        r = sqrt(x(i+1)**2 + y(j+1)**2)
        if (r < g%w) rij_em(i,j) = minloc((r_em-r)**2, dim=1)
      end do
    end do

    do j = 1, ny+2
      do i = 1, nx+2
        if (((abs(x(i)) - 15.8d-3 / x0 / 2d0) > 0d0) &
          .and. ((20.8d-3 / x0 / 2d0 - abs(x(i))) > 0d0) &
          .and. ((abs(y(j)) - 15.8d-3 / x0 / 3d0) < 0d0)) then
          eps(i,j) = -1d40
        end if
      end do
    end do

    call MPI_File_Open(comm, 'Output/Hz.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (myId == 0) call MPI_File_Delete('Output/Hz.dat', info, ierr);
        call MPI_File_Open(comm, 'Output/Hz.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)
  end subroutine

! *** Run ***
  subroutine em_run(g, ne, nte, ampl)
    type(grid), intent(in) :: g
    real(8), intent(in) :: ne(:,:), nte(:,:), ampl
    integer :: i, j, t, k
    real(8) :: r, r1, r0, Te, temp1, temp2
    complex(8) :: etemp, htemp

    ! average ne(y,r) to ne_avg(r)
    do j = 2, g%by+1
      temp(j-1) = sum(ne(2:g%bx+1,j)) / float(g%bx)
      call MPI_AllReduce(MPI_In_Place, temp(j-1), 1, etype, MPI_Sum, ryComm, ierr)
      temp(j-1) = temp(j-1) / float(px)
    end do
    call MPI_AllGather(temp, g%by, etype, ne_avg, g%by, etype, rxComm, ierr)

    do j = 2, g%by+1
      temp(j-1) = sum(nte(2:g%bx+1,j)) / float(g%bx)
      call MPI_AllReduce(MPI_In_Place, temp(j-1), 1, etype, MPI_Sum, ryComm, ierr)
      temp(j-1) = temp(j-1) / float(px)
    end do
    call MPI_AllGather(temp, g%by, etype, nte_avg, g%by, etype, rxComm, ierr)

    ! convert ne_avg(r) to wp(x,y)
    do j = 2, ny+1
      do i = 2, nx+1
        k = rij(i,j)
        if (k > 0) then
          r = sqrt(x(i)**2 + y(j)**2)
          temp1 = (ne_avg(k-1) +  (r - r_full(k-1)) * (ne_avg(k) - ne_avg(k-1)) &
                  / (r_full(k) - r_full(k-1)))
          temp2 = (nte_avg(k-1) +  (r - r_full(k-1)) * (nte_avg(k) - nte_avg(k-1)) &
                  / (r_full(k) - r_full(k-1)))
          wp(i,j) = wp2 * temp1
          Te = temp1 / temp2
          nu(i,j) = get_nu(Te)
        end if
      end do
    end do

    call comm_real(nx, ny, wp)
    call comm_real(nx, ny, eps)

    if (nf > 0) then
      a_ey = 0
      a_hz = 0
    end if

    do t = 1, nt
      call em_step(t, fskip, ampl)
      ! if (t > nt-200) call savedat('Output/Hz.dat', Ey)
      if (t > (nt-int(np/2d0))) then
        do j = 2, ny+1
          do i = 2, nx+1
            k = rij_em(i-1,j-1)
            if (k > 0) then
              Em_full(k) = Em_full(k) + 0.05 * (Ey(i,j)**2 - Em_full(k))
            end if
          end do
        end do
      end if
    end do

    if (nf > 0) then
      do k = 1, nf
        frq(k) = float(k + fskip)  / (nt * dt)
        call MPI_Reduce(a_ey(k), etemp, 1, MPI_Complex16, MPI_SUM, 0, comm, ierr)
        call MPI_Reduce(a_hz(k), htemp, 1, MPI_Complex16, MPI_SUM, 0, comm, ierr)
        flx(k) = 2 * RealPart(etemp * conjg(htemp)) / float(nt)
        write(*,33) frq(k)/1d9, flx(k)
      end do
    end if

    ! convert Em_full(r_em) to Em_r(r)
    do j = 1, size(Em_full)
      call MPI_AllReduce(MPI_In_Place, Em_full(j), 1, etype, MPI_Sum, comm, ierr)
    end do
    Em_full = Em_full / float(nproc)

    do j = 1, g%by
      k = minloc((g%r(j+1) - r_em)**2, dim=1) !replace this with stored array
      Em(j) = (Em_full(k) +  (g%r(j+1) - r_em(k)) * (Em_full(k+1) - Em_full(k)) &
              / (r_em(k+1) - r_em(k)))
    end do

    ! write(*,34) Em
    ! 34 format(10f6.3)
    
    33 format(2es11.3)
  end subroutine

! *** Timestepping Routine ***
  subroutine em_step(t, fskip, ampl)
    integer, intent(in) :: t, fskip
    real(8), intent(in) :: ampl
    integer :: i, j, k, j1, j2
    real(8) :: arg, md = 1
    complex :: ii = (0,1)
    ! .H0 |E0 .H1 |E1 .H2 |E2 ... |En-2 .Hn-2 |En-1 .Hn-1 |En .Hn+1
    ! H ranges 0 -> n+1, with n values and 2 boundary conditions
    ! E ranges from 0 -> n for n+1 values and no boundary conditions

    ! Ex and Ey Update (whole domain):
    do i = 1, nx+1
        Ey(i,:) = Ey(i,:) - 2. * a * (Hz(i+1,:) - Hz(i,:) - Jy(i,:)) &
                / (eps(i+1,:) + eps(i,:))
    end do
    do j = 1, ny+1
        Ex(:,j) = Ex(:,j) + 2. * a * (Hz(:,j+1) - Hz(:,j) + Jx(:,j)) &
                / (eps(:,j+1) + eps(:,j))
    end do

    ! Hz wavelet source:
    arg = pi**2 * (float(t) / np - md)**2
    j1 = 1
    j2 = ny+2
    if (ry == 0) j1 = 5
    if (ry == py-1) j2 = ny-3
    if (nf == 0) then
      Hz(4,j1:j2) = Hz(4,j1:j2) + ampl * sin(2d0 * pi * float(t) / np)
    else
      if (rx == 0) Hz(4,j1:j2) = Hz(4,j1:j2) + ampl * (1d0 - 2d0 * arg) * exp(-arg)
    end if

    ! Hz update (whole domain):
    do j = 2, ny+1
      do i = 2, nx+1
        Hz(i,j) = Hz(i,j) + a * (Ex(i,j) - Ex(i,j-1) &
                                - Ey(i,j) + Ey(i-1,j))
      end do
    end do

    ! ABC on Hz:
    ! - top -
    if ((ry == 0) .and. (.False.)) then
      Hz(:,1) = -1. / (a + 1./a + 2.) * ((a + 1./a - 2.) * (Hz(:,3) + tbc(:,1,2)) &
                +2. * (a - 1./a) * (tbc(:,1,1) + tbc(:,3,1) - Hz(:,2) - tbc(:,2,2)) &
                -4. * (a + 1./a) * tbc(:,2,1)) - tbc(:,3,2)
      tbc(:,:,2) = tbc(:,:,1)
      tbc(:,:,1) = Hz(:,1:3)
    end if
    ! - bottom -
    if ((ry == py-1) .and. (.False.)) then
      Hz(:,ny+2) = -1. / (a + 1./a + 2.) * ((a + 1./a - 2.) * (Hz(:,ny) + bbc(:,3,2)) &
                   +2. * (a - 1./a) * (bbc(:,3,1) + bbc(:,1,1) - Hz(:,ny+1) - bbc(:,2,2)) &
                   -4. * (a + 1./a) * bbc(:,2,1)) - bbc(:,1,2)
      bbc(:,:,2) = bbc(:,:,1)
      bbc(:,:,1) = Hz(:,ny:ny+2)
    end if
    ! - left -
    if (rx == 0) then
      Hz(1,:) = -1. / (a + 1./a + 2.) * ((a + 1./a - 2.) * (Hz(3,:) + lbc(1,:,2)) &
                +2. * (a - 1./a) * (lbc(1,:,1) + lbc(3,:,1) - Hz(2,:) - lbc(2,:,2)) &
                -4. * (a + 1./a) * lbc(2,:,1)) - lbc(3,:,2)
      lbc(:,:,2) = lbc(:,:,1)
      lbc(:,:,1) = Hz(1:3,:)
    end if
    ! - right -
    if (rx == px-1) then
      Hz(nx+2,:) = -1. / (a + 1./a + 2.) * ((a + 1./a - 2.) * (Hz(nx,:) + rbc(3,:,2)) &
                   +2. * (a - 1./a) * (rbc(3,:,1) + rbc(1,:,1) - Hz(nx+1,:) - rbc(2,:,2)) &
                   -4. * (a + 1./a) * rbc(2,:,1)) - rbc(1,:,2)
      rbc(:,:,2) = rbc(:,:,1)
      rbc(:,:,1) = Hz(nx:nx+2,:)
    end if

    ! Perfectly conducting (E = 0) boundaryies:
    ! x-dir
    !Hz(1,:) = Hz(1,:) - a * Ey(1,:)
    !Hz(nx,:) = Hz(nx,:) + a * Ey(nx-1,:)

    ! y-dir
    !Hz(:,1) = Hz(:,1) + a * Ex(:,1)
    !Hz(:,ny) = Hz(:,ny) - a * Ex(:,ny-1)

    ! Jx and Jy Update (whole domain):
    do i = 1, nx+1
      Jy(i,:) = (Jy(i,:) - 0.5 * a * (wp(i+1,:) + wp(i,:)) * Ey(i,:)) &
                * (1-nu(i,:))/(1+nu(i,:))
    end do
    do j = 1, ny+1
      Jx(:,j) = (Jx(:,j) - 0.5 * a * (wp(:,j+1) + wp(:,j)) * Ex(:,j)) &
                * (1-nu(:,j))/(1+nu(:,j))
    end do

    call comm_real(nx, ny, Hz)

    ! Update Fourier coefficients
    if (nf > 0) then
      do k = 1, nf
        do j = 2, ny+1
          a_ey(k) = a_ey(k) + 0.5 * (Ey(nx, j) + Ey(nx+1, j)) &
                      * exp( -ii * (k+fskip) * 2 * pi * t / float(nt))
          a_hz(k) = a_hz(k) + Hz(nx+1, j) &
                      * exp( -ii * (k+fskip) * 2 * pi * t / float(nt))
        end do
      end do
    end if
  end subroutine

  function get_nu(T)
    real(8):: get_nu
    real(8), intent(in):: T
    real(8):: x, &
              a = -32.275912575,      &
              b =   1.45173283977,    &
              c =   0.00936933121094, &
              d =   0.129397015353,   &
              f =  -0.0414865809044,  &
              g =  -0.0582934303409,  &
              h =   0.0309832277826,  &
              i =  -0.00542014733763,  &
              j =   0.000325615321708

    x = log(max(2.34d-1, min(1.57d2, T * ph0)))

    get_nu = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. + g*x**5. &
                 + h*x**6. + i*x**7 + j*x**8.) / x0**3 * ninf * dt
    return
  end function get_nu
end module
