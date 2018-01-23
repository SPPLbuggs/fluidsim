! *** External Circuit Module ***
module circ_lib
  use props
  use ptcl_lib
  implicit none

  real(8) :: Vmax, Vd_pl, Id, Vsrc, Cap, Res

  ! variables
  public  :: Vd_pl, Id
  private :: Vmax, Vsrc, Cap, Res

  ! subroutines
  public  :: circ_step, circ_init

contains

  ! *** External Circuit Timestepping ***
  subroutine circ_step(g, rf, ph, R0)
    type(grid), intent(in) :: g
    logical, intent(in) :: rf
    real(8), intent(in) :: R0
    real(8), intent(inout) :: ph(:,:)
    integer :: i, j
    real(8) :: num, den, mue, Te, ve, a, temp, dr

    if (rf) then
      ! Vsrc = Vmax
      Vsrc = Vmax * cos(pi * g%t)
      ! if (g%t < 2.5d-1) then
      !   Vsrc = Vmax * sin(4 * pi * g%t)**3
      ! else
      !   Vsrc = 0
      ! end if
    end if

    Res = R0 * e / (ph0 * t0)

    i = 2
    num = 0
    den = 0
    Id = 0d0
    do j = 2, g%by+1
      if (g%type_x(i-1,j-1) == -2) then
        Te = get_Te(nt(i,j,1), ne(i,j,1))
        mue = get_mue(Te)
        ve = sqrt((16d0 * e * ph0 * Te) / (3d0 * pi * me)) * t0 / x0

        a = 0d0
        if (Vd_pl < ph(i,j)) a = 1d0

        temp = (a * (1d0 + gam) * mui * ni(i,j,1) &
               + (1d0 - a) * mue * ne(i,j,1)) / g%dx(i-1)

        if (g%ny > 1) then
          if (cyl) then
            dr = g%dy(j-1) * g%r(j-1) * 2d0 * pi
          else
            dr = g%dy(j-1) * g%w
          end if
        else
          dr = g%w**2 * pi
        end if

        den = den + g%dt / Cap * temp * dr

        num = num + g%dt / Cap * (ph(i,j) * temp &
              + (1d0 + gam) / 4d0 * vi * ni(i,j,1) &
              - 1d0 / 4d0 * ve * ne(i,j,1)) * dr

        Id = Id + ((Vd_pl - ph(i,j)) * temp &
             - (1 + gam) / 4.0 * vi * ni(i,j,1) &
             + 1.0 / 4.0 * ve * ne(i,j,1)) * dr
      end if
    end do

    call MPI_Allreduce(MPI_In_Place, num, 1, etype, MPI_Sum, comm, ierr)
    call MPI_Allreduce(MPI_In_Place, den, 1, etype, MPI_Sum, comm, ierr)

    num = num + Vd_pl + g%dt / (Cap * Res) * Vsrc
    den = den + 1d0 + g%dt / (Cap * Res)

    Vd_pl = num / den

    if (rf) Vd_pl = Vsrc
    do j = 2, g%by+1
      if (g%type_x(i-1,j-1) == -2) ph(i-1,j) = Vd_pl
    end do
  end subroutine

  ! *** External Circuit Initialization ***
  subroutine circ_init(Vset, R0)
    real(8), intent(in) :: Vset, R0

    Vmax  = Vset
    Vsrc  = Vset
    Vd_pl = Vset
    Cap = 1d-12 * ph0 / e
    Res = R0 * e / (ph0 * t0)
  end subroutine
end module
