module em_lib
  use props
  implicit none

  type :: em_grid
    integer :: nx, ny, bx, by, &
               offx, offy, &
               nf, nt, fskip
    real(8) :: dt, dx, w, np
    real(8), allocatable :: x(:), y(:), r(:,:)
  end type

  real(8), allocatable :: Hz(:,:), Chzh(:,:), Chze(:,:), &
                          HzLeft(:,:,:), HzRight(:,:,:), &
                          Ex(:,:), Cexh(:,:), Cexe(:,:), &
                          Ey(:,:), Ceyh(:,:), Ceye(:,:), &
                          ExTemp(:,:), EyTemp(:,:), Em(:), &
                          Jx(:,:), Cjxj(:,:), Cjxe(:,:), &
                          Jy(:,:), Cjyj(:,:), Cjye(:,:), &
                          eps(:,:), sig(:,:), &
                          Nw(:,:), Ng(:,:), &
                          frq(:), flx(:)
  complex(8), allocatable :: F_ey(:), F_hz(:)

  integer :: t
  real(8) :: coef0, coef1, coef2, c, fldAmp = 0
  real(8), parameter :: a = 1d0 / sqrt(2d0), imp0=1.0
  logical :: Finit = .True., ABCinit = .True.

  public  :: em_init, em_step, Em
  private :: updateH, updateE, updateF, updateJ, &
             setC, src, abc, getNu, &
             Hz, Chzh, Chze,  HzLeft, HzRight, &
             Ex, Cexh, Cexe, &
             Ey, Ceyh, Ceye, &
             Jx, Cjxj, Cjxe, &
             Jy, Cjyj, Cjye, &
             coef0, coef1, coef2, ABCinit, &
             c, eps, sig, a, Nw, Ng, &
             frq, flx, F_ey, F_hz, Finit

contains
  ! *** Initialization
  subroutine em_init(g, eg)
    type(grid), intent(in) :: g
    type(em_grid), intent(inout) :: eg
    integer :: i, j
    real(8) :: xtemp

    ! Setup domain
    eg%nx = 3 * g%ny
    eg%ny = g%ny

    eg%bx = eg%nx / px
    eg%by = eg%ny / py

    eg%offx = rx * eg%bx
    eg%offy = ry * eg%by

    eg%w = 2d0 * g%w
    eg%dx = eg%w / float(eg%ny + 1)
    eg%dt = x0 * eg%dx / c0 / sqrt(2d0)
    c = c0 * eg%dt / eg%dx

    eg%np = 1d0 / (13.622d9 * eg%dt)
    eg%nt = 200 * int(eg%np)
    eg%nf = 0
    eg%fskip = 175

    allocate(eg%x(eg%bx+2))
    do i = 1, eg%bx+2
      eg%x(i) = eg%dx * (eg%offx + i - 1)
    end do
    xtemp = eg%x(eg%bx+2)
    call MPI_Bcast(xtemp, 1, MPI_Real8, nproc-1, comm, ierr)
    eg%x = eg%x - xtemp / 2d0

    allocate(eg%y(eg%by+2))
    do j = 1, eg%by+2
      eg%y(j) = eg%dx * (eg%offy + j - 1)
    end do
    xtemp = eg%y(eg%by+2)
    call MPI_Bcast(xtemp, 1, MPI_Real8, nproc-1, comm, ierr)
    eg%y = eg%y - xtemp / 2d0

    allocate(eg%r(eg%bx+2, eg%by+2))
    do j = 1, eg%by+2
      do i = 1, eg%bx+2
        eg%r(i,j) = sqrt(eg%x(i)**2 + eg%y(j)**2)
      end do
    end do

    ! Setup field variables
    allocate(Hz(eg%bx+2, eg%by+2), Chzh(eg%bx+2, eg%by+2), Chze(eg%bx+2, eg%by+2), &
             HzLeft(3, 2, eg%by+2), HzRight(3, 2, eg%by+2), &
             Ex(eg%bx+2, eg%by+2), Cexh(eg%bx+2, eg%by+2), Cexe(eg%bx+2, eg%by+2), &
             Ey(eg%bx+2, eg%by+2), Ceyh(eg%bx+2, eg%by+2), Ceye(eg%bx+2, eg%by+2), &
             ExTemp(eg%bx+2, eg%by+2), EyTemp(eg%bx+2, eg%by+2), &
             Jx(eg%bx+2, eg%by+2), Cjxj(eg%bx+2, eg%by+2), Cjxe(eg%bx+2, eg%by+2), &
             Jy(eg%bx+2, eg%by+2), Cjyj(eg%bx+2, eg%by+2), Cjye(eg%bx+2, eg%by+2), &
             eps(eg%bx+2, eg%by+2), sig(eg%bx+2, eg%by+2), &
             Nw(eg%bx+2, eg%by+2), Ng(eg%bx+2, eg%by+2), Em(g%by+2))
    Hz = 0;  HzLeft = 0; HzRight = 0
    Ex = 0;  Ey = 0; Em = 0
    Jx = 0;  Jy = 0
    eps = 1; sig = 0
    Nw = 0;  Ng = 0

    do j = 1, eg%by+2
      do i = 1, eg%bx+2
        if (((abs(eg%x(i)) - 15.8d-3 / x0 / 2d0) > 0d0) &
          .and. ((23.8d-3 / x0 / 2d0 - abs(eg%x(i))) > 0d0) &
          .and. ((abs(eg%y(j)) - 15.8d-3 / x0 / 3d0) < 0d0)) then
          eps(i,j) = -1d40
        end if
      end do
    end do

    call comm_real(eg%bx, eg%by, eps)

    do j = 1, eg%by+2
      do i = 1, eg%bx+2
        Chzh(i,j) = 1d0
        Chze(i,j) = a / imp0
      end do
    end do

    call comm_real(eg%bx, eg%by, Chzh)
    call comm_real(eg%bx, eg%by, Chze)

    call MPI_File_Open(comm, 'Output/Hz.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (myId == 0) call MPI_File_Delete('Output/Hz.dat', info, ierr);
        call MPI_File_Open(comm, 'Output/Hz.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)
  end subroutine

  subroutine setC(eg, fp, gp)
    type(em_grid), intent(in) :: eg
    real(8), intent(in) :: fp, gp
    integer :: i, j

    Nw = 1d40
    Ng = 1d40
    do j = 1, eg%by+2
      do i = 1, eg%bx+2
        if (eg%r(i,j) < eg%w/2d0) then
          Nw(i,j) = 1d0 / cos(pi * eg%r(i,j) / eg%w)**2 / (fp * eg%dt)
          Ng(i,j) = 1d0 / cos(pi * eg%r(i,j) / eg%w)**2 / (gp * eg%dt)
        end if
      end do
    end do

    do j = 1, eg%by+1
      do i = 1, eg%bx+2
        Cjxj(i,j) = (1d0 - 1d0 / (Ng(i,j+1) + Ng(i,j))) &
               / (1d0 + 1d0 / (Ng(i,j+1) + Ng(i,j)))
        Cjxe(i,j) = 8d0 * pi**2 * a / imp0 / (Nw(i,j+1) + Nw(i,j))**2 &
               / (1d0 + 1d0 / (Ng(i,j+1) + Ng(i,j)))
      end do
    end do


    do j = 1, eg%by+2
      do i = 1, eg%bx+1
        Cjyj(i,j) = (1d0 - 1d0 / (Ng(i+1,j) + Ng(i,j))) &
               / (1d0 + 1d0 / (Ng(i+1,j) + Ng(i,j)))
        Cjye(i,j) = 8d0 * pi**2 * a / imp0 / (Nw(i+1,j) + Nw(i,j))**2 &
               / (1d0 + 1d0 / (Ng(i+1,j) + Ng(i,j)))
      end do
    end do

    do j = 1, eg%by+1
      do i = 1, eg%bx+2
        Cexh(i,j) = 2d0 * a * imp0 / (eps(i,j+1) + eps(i,j)) &
                    / (1d0 + (sig(i,j+1) + sig(i,j)) / (eps(i,j+1) + eps(i,j)) / 2d0 &
                    + a * imp0 * Cjxe(i,j) / (eps(i,j+1) + eps(i,j)))
        Cexe(i,j) = (1d0 - (sig(i,j+1) + sig(i,j)) / (eps(i,j+1) + eps(i,j)) / 2d0 &
                    - a * imp0 * Cjxe(i,j) / (eps(i,j+1) + eps(i,j))) &
                    / (1d0 + (sig(i,j+1) + sig(i,j)) / (eps(i,j+1) + eps(i,j)) / 2d0 &
                    + a * imp0 * Cjxe(i,j) / (eps(i,j+1) + eps(i,j)))
      end do
    end do

    do j = 1, eg%by+2
      do i = 1, eg%bx+1
        Ceyh(i,j) = 2d0 * a * imp0 / (eps(i+1,j) + eps(i,j)) &
                    / (1d0 + (sig(i+1,j) + sig(i,j)) / (eps(i+1,j) + eps(i,j)) / 2d0 &
                    + a * imp0 * Cjye(i,j) / (eps(i+1,j) + eps(i,j)))
        Ceye(i,j) = (1d0 - (sig(i+1,j) + sig(i,j)) / (eps(i+1,j) + eps(i,j)) / 2d0 &
                    - a * imp0 * Cjye(i,j) / 2d0) &
                    / (1d0 + (sig(i+1,j) + sig(i,j)) / (eps(i+1,j) + eps(i,j)) / 2d0 &
                    + a * imp0 * Cjye(i,j) / (eps(i+1,j) + eps(i,j)))
      end do
    end do
  end subroutine

  subroutine addCyl(eg, x0, y0, rad, epsCyl)
    type(em_grid), intent(in) :: eg
    real(8), intent(in) :: x0, y0, rad, epsCyl
    integer :: i, j
    real(8) :: r

    do j = 1, eg%by+2
      do i = 1, eg%bx+2
        r = sqrt((eg%x(i) - x0)**2 + (eg%y(j) - y0)**2)
        if (r .le. rad) eps(i,j) = epsCyl
      end do
    end do
  end subroutine

  subroutine addBlock(eg, x0, y0, dx, dy, epsBlk)
    type(em_grid), intent(in) :: eg
    real(8), intent(in) :: x0, y0, dx, dy, epsBlk
    integer :: i, j
    real(8) :: x, y

    do j = 1, eg%by+2
      do i = 1, eg%bx+2
        x = abs(eg%x(i) - x0)
        y = abs(eg%y(j) - y0)
        if ((x .le. dx) .and. (y .le. dy)) eps(i,j) = epsBlk
      end do
    end do
  end subroutine

  ! *** Timestepping Routine ***
  subroutine em_step(g, eg, ne, Te, ampIn, ampOut)
    type(grid), intent(in) :: g
    type(em_grid), intent(inout) :: eg
    real(8), intent(in) :: ne(:,:), Te(:,:), ampIn
    real(8), intent(inout) :: ampOut
    integer :: j, k, ttemp
    real(8) :: fp, gp, temp
    complex(8) :: etemp, htemp

    fp = sqrt(ne(g%bx/2, 2) * e**2 / (eps0 * x0**3 * me)) * eg%dt
    gp = getNu(Te(g%bx/2,2), eg%dt)

    call setC(eg, fp, gp)

    t = mod(t, 20*int(eg%np))
    do ttemp = 1, max(int(g%dt * t0 / eg%dt), 1000*int(eg%np))!int(eg%np)/2)
      t = t + 1
      call updateE(eg)
      call updateJ(eg)
      call src(eg, t, ampIn)
      call updateH(eg)
      call abc(eg)
      if (eg%nf > 0) call updateF(eg, t)
      ! if (t > eg%nt-200) call savedat('Output/Hz.dat', Ey)
      fldAmp = fldAmp + (Ey(eg%nx/2, eg%ny/2)**2 - fldAmp) / eg%np
    end do

    do j = 1, g%by+2
      Em(j) = fldAmp * cos(pi * g%r(j) / g%w)
    end do

    if (eg%nf > 0) then
      do k = 1, eg%nf
        frq(k) = float(k + eg%fskip)  / (eg%nt * eg%dt)
        call MPI_Reduce(F_ey(k), etemp, 1, MPI_Complex16, MPI_SUM, 0, comm, ierr)
        call MPI_Reduce(F_hz(k), htemp, 1, MPI_Complex16, MPI_SUM, 0, comm, ierr)
        flx(k) = 2 * RealPart(etemp * conjg(htemp)) / float(eg%nt)
        write(*,33) frq(k)/1d9, flx(k)
      end do
    end if

    temp = sum(-5d-1 * (Ey(eg%bx-2,3:eg%by) + Ey(eg%bx-1,3:eg%by)) * Hz(eg%bx-1,3:eg%by))
    ampOut = ampOut + 0.01 * (temp - ampOut)

    write(*,*) fldAmp

    fp = 10d9 * eg%dt
    gp = 1d9 * eg%dt

    call setC(eg, fp, gp)

    t = mod(t, 20*int(eg%np))
    do ttemp = 1, max(int(g%dt * t0 / eg%dt), 1000*int(eg%np))!int(eg%np)/2)
      t = t + 1
      call updateE(eg)
      call updateJ(eg)
      call src(eg, t, ampIn)
      call updateH(eg)
      call abc(eg)
      if (eg%nf > 0) call updateF(eg, t)
      ! if (t > eg%nt-200) call savedat('Output/Hz.dat', Ey)
      fldAmp = fldAmp + (Ey(eg%nx/2, eg%ny/2)**2 - fldAmp) / eg%np
    end do

    temp = sum(-5d-1 * (Ey(eg%bx-2,3:eg%by) + Ey(eg%bx-1,3:eg%by)) * Hz(eg%bx-1,3:eg%by))

    write(*,*) fldAmp
    stop



    33 format(f7.3, 10es11.3)
  end subroutine

  ! *** Update Hz Field ***
  subroutine updateH(eg)
    type(em_grid), intent(in) :: eg
    integer :: i, j

    do j = 2, eg%by+1
      do i = 2, eg%bx+1
        Hz(i,j) = Chzh(i,j) * Hz(i,j) &
                + Chze(i,j) * ((Ex(i,j) - Ex(i,j-1)) - (Ey(i,j) - Ey(i-1,j)))
      end do
    end do

    call comm_real(eg%bx, eg%by, Hz)
  end subroutine

  ! *** Update Ex and Ey Fields ***
  subroutine updateE(eg)
    type(em_grid), intent(in) :: eg
    integer :: i, j

    ExTemp = Ex
    EyTemp = Ey

    do j = 1, eg%by+1
      do i = 1, eg%bx+2
        Ex(i,j) = Cexe(i,j) * Ex(i,j) + Cexh(i,j) * (Hz(i,j+1) - Hz(i,j) &
                  - 5d-1 * (1d0 + Cjxj(i,j)) * Jx(i,j))
      end do
    end do

    do j = 1, eg%by+2
      do i = 1, eg%bx+1
        Ey(i,j) = Ceye(i,j) * Ey(i,j) - Ceyh(i,j) * (Hz(i+1,j) - Hz(i,j) &
                  + 5d-1 * (1d0 + Cjyj(i,j)) * Jy(i,j))
      end do
    end do
  end subroutine

  ! *** Update Jx and Jy Currents ***
  subroutine updateJ(eg)
    type(em_grid), intent(in) :: eg
    integer :: i, j

    do j = 1, eg%by+1
      do i = 1, eg%bx+2
        Jx(i,j) = Cjxj(i,j) * Jx(i,j) + Cjxe(i,j) * (Ex(i,j) + ExTemp(i,j))
      end do
    end do

    do j = 1, eg%by+2
      do i = 1, eg%bx+1
        Jy(i,j) = Cjyj(i,j) * Jy(i,j) + Cjye(i,j) * (Ey(i,j) + EyTemp(i,j))
      end do
    end do
  end subroutine

  ! *** Apply Sources ***
  subroutine src(eg, t, amp)
    type(em_grid), intent(in) :: eg
    integer, intent(in) :: t
    real(8), intent(in) :: amp
    integer :: j1, j2
    real(8) :: arg, md = 1d0

    ! Hz wavelet source:
    arg = pi**2 * (float(t) / eg%np - md)**2
    j1 = 1
    j2 = eg%by+2
    if (ry == 0) j1 = 5
    if (ry == py-1) j2 = eg%by-3
    if (eg%nf == 0) then
      Hz(4,j1:j2) = Hz(4,j1:j2) + amp * sin(2d0 * pi * float(t) / eg%np)
    else
      if (rx == 0) Hz(4,j1:j2) = Hz(4,j1:j2) + amp * (1d0 - 2d0 * arg) * exp(-arg)
    end if
  end subroutine

  ! *** Apply Absorbing Boundary Conditions ***
  subroutine abc(eg)
    type(em_grid), intent(in) :: eg
    real(8) :: temp1, temp2
    integer :: i, j

    if (ABCinit) then
      temp1 = sqrt(Cexh(1,1) * Chze(1,1))
      temp2 = 1d0 / temp1 + 2d0 + temp1

      coef0 = -(1d0 / temp1 - 2d0 + temp1) / temp2
      coef1 = -2d0 * (temp1 - 1d0 / temp1) / temp2
      coef2 = 4d0 * (temp1 + 1d0 / temp1) / temp2

      ABCinit = .False.
    end if

    ! ABC on left side
    if (rx == 0) then
      do j = 2, eg%by+1
        Hz(1, j) = coef0 * (Hz(3,j) + HzLeft(1,2,j)) &
                 + coef1 * (HzLeft(1,1,j) + HzLeft(3,1,j) &
                          - Hz(2,j) - HzLeft(2,2,j)) &
                 + coef2 * HzLeft(2,1,j) - HzLeft(3,2,j)
        ! Save old fields
        do i = 1, 3
          HzLeft(i,2,j) = HzLeft(i,1,j)
          HzLeft(i,1,j) = Hz(i,j)
        end do
      end do
    end if

    ! ABC on right side
    if (rx == px - 1) then
      do j = 2, eg%by+1
        Hz(eg%bx+2, j) = coef0 * (Hz(eg%bx,j) + HzRight(1,2,j)) &
                 + coef1 * (HzRight(1,1,j) + HzRight(3,1,j) &
                          - Hz(eg%bx+1,j) - HzRight(2,2,j)) &
                 + coef2 * HzRight(2,1,j) - HzRight(3,2,j)
        ! Save old fields
        do i = 1, 3
          HzRight(i,2,j) = HzRight(i,1,j)
          HzRight(i,1,j) = Hz(eg%bx+3-i,j)
        end do
      end do
    end if

    call comm_real(eg%bx, eg%by, Hz)
  end subroutine

  ! *** Update Fourier Components ***
  subroutine updateF(eg, t)
    type(em_grid), intent(in) :: eg
    integer, intent(in) :: t
    integer :: k, j
    complex :: ii = (0,1)

    if (Finit) then
      allocate(frq(eg%nf), flx(eg%nf), F_ey(eg%nf), F_hz(eg%nf))
      frq = 0
      flx = 0
      F_ey = 0
      F_hz = 0
      Finit = .False.
    end if

    do k = 1, eg%nf
      do j = 3, eg%by
        F_ey(k) = F_ey(k) + 5d-1 * (Ey(eg%bx-2, j) + Ey(eg%bx-1, j)) &
                * exp( -ii * (k+eg%fskip) * 2d0 * pi * t / float(eg%nt))
        F_hz(k) = F_hz(k) + Hz(eg%bx-1, j) &
                * exp( -ii * (k+eg%fskip) * 2d0 * pi * t / float(eg%nt))
      end do
    end do
  end subroutine

  function getNu(T,dt)
    real(8):: getNu
    real(8), intent(in):: T, dt
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

    getNu = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. + g*x**5. &
                + h*x**6. + i*x**7 + j*x**8.) / x0**3 * ninf * dt
    return
  end function getNu
end module
