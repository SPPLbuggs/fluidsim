module em_lib
  use props
  implicit none

  integer, allocatable :: rij(:,:)
  real(8), allocatable :: Hz(:,:),   Ey(:,:),  Ex(:,:),  Em(:), &
                          Jx(:,:),   Jy(:,:),  wp(:,:), nu(:,:), eps(:,:), &
                          tbc(:,:,:), rbc(:,:,:), lbc(:,:,:), bbc(:,:,:), &
                          x(:), y(:), frq(:), flx(:)
  complex(8), allocatable :: a_ey(:), a_hz(:)

  ! properties
  real(8), parameter:: a = 0.99 / sqrt(2.0)
  real(8) :: em_amp = 1d0, freq = 1d0, dt, np
  integer :: nx, ny, nf, nt, fskip

  ! variables
  public  :: Em, em_amp, freq
  private :: Hz, Ey, Ex, Jx, Jy, a, &
             tbc, rbc, lbc, bbc, &
             wp, nu, x, y, eps, dt, &
             nx, ny, rij, nf, nt, fskip

  ! subroutines and functions
  public  :: em_step, em_init
  private :: em_run

contains

  subroutine em_step(g, ne, nte)
    type(grid), intent(in) :: g
    real(8), intent(in) :: ne(:,:), nte(:,:)
    integer :: i, j, k, t
    complex(8) :: etemp, htemp

    Hz = 0
    Ex = 0
    Ey = 0
    Jx = 0
    Jy = 0

    do t = 1, nt
      call em_run(t, fskip)
      call savedat('Output/f2.dat', Hz)
    end do

    do k = 1, nf
      frq(k) = float(k + fskip)  / (nt * dt)
      call MPI_Reduce(a_ey(k), etemp, 1, MPI_Complex16, MPI_SUM, 0, comm, ierr)
      call MPI_Reduce(a_hz(k), htemp, 1, MPI_Complex16, MPI_SUM, 0, comm, ierr)
      flx(k) = 2 * RealPart(etemp * conjg(htemp)) / float(nt)
    end do
  end subroutine

  subroutine em_run(t, fskip)
    integer, intent(in) :: t, fskip
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
    j1 = 1
    j2 = ny+2
    if (ry == 0) j1 = 5
    if (ry == py-1) j2 = ny-3
    ! Hz(4,j1:j2) = Hz(4,j1:j2) + em_amp * sin(2d0 * pi * float(t) * freq)
    ! Hz wavelet source:
    arg = pi**2 * (float(t) / np - md)**2
    if (rx == 0) Hz(4,j1:j2) = Hz(4,j1:j2) + (1d0 - 2d0 * arg) * exp(-arg)

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

    do k = 1, nf
      do j = 2, ny+1
        a_ey(k) = a_ey(k) + 0.5 * (Ey(nx, j) + Ey(nx+1, j)) &
                    * exp( -ii * (k+fskip) * 2 * pi * t / float(nt))
        a_hz(k) = a_hz(k) + Hz(nx+1, j) &
                    * exp( -ii * (k+fskip) * 2 * pi * t / float(nt))
      end do
    end do
  end subroutine

! *** Initialization ***
  subroutine em_init(g)
    type(grid), intent(in) :: g
    integer :: i, j, offx, offy, k
    real(8) :: Ly = 15.8e-3 / x0, r, dx, xtemp

    nx = 2 * g%ny / px
    ny = g%ny / py

    offx = rx * nx
    offy = ry * ny

    allocate(Hz(nx+2, ny+2), Ex(nx+2, ny+2), Ey(nx+2, ny+2), Jx(nx+2, ny+1), &
             Jy(nx+1, ny+2), wp(nx+2,ny+2), eps(nx+2,ny+2), nu(nx+2, ny+2), &
             lbc(3,ny+2,2), rbc(3,ny+2,2), tbc(nx+2,3,2), bbc(nx+2,3,2), &
             x(nx+2), y(ny+2), rij(nx+2,ny+2), Em(g%by+2))
    allocate(a_ey(nf), a_hz(nf))

    Hz = 0
    Ex = 0
    Ey = 0
    Jx = 0
    Jy = 0
    eps = 1
    wp = 0
    nu = 0
    lbc = 0
    rbc = 0
    tbc = 0
    bbc = 0
    a_ey = 0
    a_hz = 0

    do i = 1, nx+2
      x(i) = 2d0 * Ly / float(nx+1) * (offx + i - 1)
    end do

    xtemp = x(1)
    call MPI_Bcast(xtemp, 1, MPI_Real8, 0, comm, ierr)
    x = x - xtemp

    xtemp = x(nx+2)
    call MPI_Bcast( xtemp, 1, MPI_Real8, nproc-1, comm, ierr)
    x = x / xtemp
    x = x * Ly * 2d0 - Ly

    do j = 1, ny+2
      y(j) = Ly / float(ny+1) * (offy + j - 1)
    end do

    xtemp = y(1)
    call MPI_Bcast(xtemp, 1, MPI_Real8, 0, comm, ierr)
    y = y - xtemp

    xtemp = y(ny+2)
    call MPI_Bcast( xtemp, 1, MPI_Real8, nproc-1, comm, ierr)
    y = y / xtemp
    y = y * Ly - Ly / 2d0

    rij = -1
    do j = 1, ny+2
      do i = 1, nx+2
        r = sqrt(x(i)**2 + y(j)**2)
        do k = 1, g%by+2
          if ((g%r(k) < r) .and. (r < g%w)) then
            rij(i,j) = k
          else
            exit
          end if
        end do
      end do
    end do

    ! do j = 1, ny+2
    !   do i = 1, nx+2
    !     if (((abs(x(i)) - 15.8d-3 / x0 / 2d0) > 0d0) &
    !       .and. ((24.8d-3 / x0 / 2d0 - abs(x(i))) > 0d0) &
    !       .and. ((abs(y(j)) - 15.8d-3 / x0 / 4d0) < 0d0)) then
    !       eps(i,j) = -1d40
    !     end if
    !   end do
    ! end do

    nf = 200
    fskip = 55

    allocate(frq(nf), flx(nf))

    dx  = abs(x(2) - x(1)) * x0
    dt  = dx / c0 / sqrt(2.0)
    np  = 1.0 / (14e9 * dt)
    nt  = 200 * int(np)
  end subroutine
end module
