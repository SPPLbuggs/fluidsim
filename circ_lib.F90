! *** External Circuit Module ***
module circ_lib
  use props
  use ptcl_props
  implicit none

  real(8) :: Vmax, Vd_pl, Vd_mi, Vd_or, Id, Vsrc, Cap, Res, kv(5) = 0

  ! variables
  public  :: Vd_pl, Vd_mi, Vd_or, Id, Res
  private :: Vmax, Vsrc, Cap

  ! subroutines
  public  :: circ_step, circ_init

contains

  ! *** External Circuit Timestepping ***
  subroutine circ_step(g, stage, ph, ni, ne, nt)
    type(grid), intent(in) :: g
    integer, intent(in) :: stage
    real(8), intent(inout) :: ph(:,:), ni(:,:), ne(:,:), nt(:,:)
    integer :: i, j
    real(8) :: a, Te, Ex, mue, ve, dr, flxi, flxe

    if (rf == 0) then
      ! ** DC **
      Vsrc = Vmax
    else if (rf == 1) then
      ! ** AC **
      Vsrc = Vmax * sin(pi * g%t)
    else if (rf == 2) then
      ! ** Pulse **
      if (g%t < 2.5d-1) then
        Vsrc = Vmax * sin(4 * pi * g%t)**3
      else
        Vsrc = 0d0
      end if
    else if (rf == 3) then
      ! ** Real Pulse **
      if (g%t < 2.3d-1) then
        Vsrc = Vmax * sin(4d0 * pi * g%t)**2
      else if (g%t < 3.35d-1) then
        Vsrc = Vmax * (sin(3.81d0 * pi * (g%t + 5.85d-2)) * 7.19d-1 + 2.81d-1)
      else if (g%t < 8.2d-1) then
        Vsrc = Vmax * (sin(3.81d0 * pi * (g%t + 1.2d-1) / 3.5d0)**3 * (-4.38d-1))
      else
        Vsrc = 0d0
      end if


      ! ** Taemin's Pulse **
      ! if (g%t < 20d-3) then
      !   Vsrc = Vmax * sin(50d0 * pi * g%t)**3
      ! else
      !   Vsrc = 0d0
      ! end if
    end if

    !Res = R0 * e / (ph0 * t0)

    i = 2
    Id = 0d0
    do j = 2, g%by+1
      if (g%ny > 1) then
        if (cyl) then
          dr = g%dy(j-1) * g%r(j-1) * 2d0 * pi
        else
          dr = g%dy(j-1) * g%w
        end if
      else
        dr = g%w**2 * pi
      end if

      if (g%type_x(i-1,j-1) == -2) then
        if (Vd_mi > ph(i,j)) then
            a = 1d0 ! electrons drift
        else
            a = 0d0 ! ions drift
        end if

        Ex = (Vd_mi - ph(i,j)) / g%dx(i)
        Te = get_Te(nt(i,j), ne(i,j))
        mue = get_mue(Te)
        ve = sqrt((16d0 * e * ph0 * Te) / (3d0 * pi * me)) * t0 / x0

        ! Flux at i + 1/2
        flxi = (1d0 - a) * mui * Ex * ni(i,j) &
               - 2.5d-1 * vi * ni(i,j)
        flxe = - a * mue * Ex * ne(i,j) &
               - 2.5d-1 * ve * ne(i,j) &
               - gam * flxi

        Id = Id + dr * (flxi - flxe)
      end if
    end do

    call MPI_Allreduce(MPI_In_Place, Id, 1, etype, MPI_Sum, comm, ierr)

    kv(stage) = -(Id - (Vsrc - Vd_mi) / Res) / Cap

    if (stage == 1) then
      Vd_mi = Vd_or + g%dt * kv(1) / 3d0

    else if (stage == 2) then
      Vd_mi = Vd_or + g%dt * ( kv(1) + kv(2)) / 6d0

    else if (stage == 3) then
      Vd_mi = Vd_or + g%dt * (kv(1) + kv(3) * 3d0) / 8d0

    else if (stage == 4) then
      Vd_mi = Vd_or + g%dt * (kv(1) &
              - kv(3) * 3d0 + kv(4) * 4d0 ) / 2d0

    else
      Vd_pl = Vd_or + g%dt * (kv(1) &
              + kv(4) * 4d0 + kv(5)) / 6d0
    end if

    if (rf .ne. -1) Vd_pl = Vsrc
    if (rf .ne. -1) Vd_mi = Vsrc
    i = 2
    do j = 2, g%by+1
      if (g%type_x(i-1,j-1) == -2) then
        if (stage == 5) then
          ph(i-1,j) = Vd_pl
        else
          ph(i-1,j) = Vd_mi
        end if
      end if
    end do
  end subroutine

  ! *** External Circuit Initialization ***
  subroutine circ_init(Vset)
    real(8), intent(in) :: Vset

    Vmax  = Vset
    Vsrc  = Vset
    Vd_pl = Vset
    Vd_mi = Vset
    Vd_or = Vset
    Cap = 1d-12 * ph0 / e
    !Res = R0 * e / (ph0 * t0)
  end subroutine
end module
