! *** External Circuit Module ***
module circ_lib
  use props
  use ptcl_props
  implicit none

  real(8) :: Vmax, Vd_pl, Vd_mi, Vd_or, Ida, Idc, Vsrc, Cap, Res, kv(5) = 0

  ! variables
  public  :: Vd_pl, Vd_mi, Vd_or, Ida, Idc, Res
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
    real(8) :: Te(2), Ex, mue, De, dr, flxi, flxe

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
    end if

    ! i = 2
    ! Ida = 0d0
    ! do j = 2, g%by+1
    !   if ((g%type_x(i-1,j-1) == -2) .and. (g%r(j) < g%ew)) then
    !     if (g%ny > 1) then
    !       if (cyl) then
    !         dr = g%dy(j-1) * g%r(j-1) * 2d0 * pi
    !       else
    !         dr = g%dy(j-1) * g%w
    !       end if
    !     else
    !       dr = g%w**2 * pi
    !     end if
    !     if (Vd_mi > ph(i,j)) then
    !         a = 1d0 ! electrons drift
    !     else
    !         a = 0d0 ! ions drift
    !     end if
    !
    !     Ex = (Vd_mi - ph(i,j)) / g%dx(i)
    !     Te = get_Te(nt(i,j), ne(i,j))
    !     mue = get_mue(Te)
    !     ve = sqrt((16d0 * e * ph0 * Te) / (3d0 * pi * me)) * t0 / x0
    !
    !     ! Flux at i + 1/2
    !     flxi = (1d0 - a) * mui * Ex * ni(i,j) &
    !            - 2.5d-1 * vi * ni(i,j)
    !     flxe = - a * mue * Ex * ne(i,j) &
    !            - 2.5d-1 * ve * ne(i,j) &
    !            - gam * flxi
    !
    !     Ida = Ida + dr * (flxi - flxe)
    !   end if
    ! end do

    i = 3
    Ida = 0
    do j = 2, g%by+1
      if (g%type_x(i-2,j-1) == -2) then
        if (g%ny > 1) then
          if (cyl) then
            dr = g%dy(j-1) * g%r(j-1) * 2d0 * pi
          else
            dr = g%dy(j-1) * g%w
          end if
          if (g%r(j) < g%ew) dr = 0d0
        else
          dr = g%w**2 * pi
        end if

        Ex = -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)

        Te(1) = get_Te(nt(i-1,j), ne(i-1,j))
        Te(2) = get_Te(nt(i,j),   ne(i,j))

        mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        De = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))

        ! Flux at i - 1/2
        call get_flux(flxi, Ex, g%dx(i-1),  1, mui, Di, &
                      ni(i-1,j), ni(i,j))
        call get_flux(flxe, Ex, g%dx(i-1), -1, mue, De, &
                      ne(i-1,j), ne(i,j))

        Ida = Ida + dr * (flxi - flxe)
      end if
    end do

    i = g%bx
    Idc = 0
    do j = 2, g%by+1
      if (g%type_x(i,j-1) == 2) then
        if (g%ny > 1) then
          if (cyl) then
            dr = g%dy(j-1) * g%r(j-1) * 2d0 * pi
          else
            dr = g%dy(j-1) * g%w
          end if
          if (g%r(j) < g%ew) dr = 0d0
        else
          dr = g%w**2 * pi
        end if

        Ex = -(ph(i+1,j) - ph(i,j)) / g%dx(i)

        Te(1) = get_Te(nt(i,j), ne(i,j))
        Te(2) = get_Te(nt(i+1,j),   ne(i+1,j))

        mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        De = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))

        ! Flux at i - 1/2
        call get_flux(flxi, Ex, g%dx(i),  1, mui, Di, &
                      ni(i,j), ni(i+1,j))
        call get_flux(flxe, Ex, g%dx(i), -1, mue, De, &
                      ne(i,j), ne(i+1,j))

        Idc = Idc + dr * (flxi - flxe)
      end if
    end do

    call MPI_Allreduce(MPI_In_Place, Ida, 1, etype, MPI_Sum, comm, ierr)
    call MPI_Allreduce(MPI_In_Place, Idc, 1, etype, MPI_Sum, comm, ierr)

    kv(stage) = -(Idc - (Vsrc - Vd_mi) / Res) / Cap

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

    if (rf .ne. 1) Vd_pl = Vsrc
    if (rf .ne. 1) Vd_mi = Vsrc
    i = 2

    if (g%ny > 1) then
      do j = 2, g%by+1
        if (g%r(j) .le. g%ew) then
          if (g%type_x(i-1,j-1) == -2) then
            if (stage == 5) then
              ph(i-1,j) = Vd_pl / 2d0
            else
              ph(i-1,j) = Vd_mi / 2d0
            end if
          end if
          if (g%type_x(g%bx,j-1) == 2) then
            if (stage == 5) then
              ph(g%bx+2,j) = Vd_pl / (-2d0)
            else
              ph(g%bx+2,j) = Vd_mi / (-2d0)
            end if
          end if
        end if
      end do
    else
      if (g%type_x(1,1) == -2) ph(1,2) = Vsrc / 2d0
      if (g%type_x(g%bx,1) == 2) ph(g%bx+2,2) = Vsrc / (-2d0)
    end if


    ! do j = 1, g%by+2
    !   write(*,32) ph(:,j)
    !   32 format(12f7.2)
    ! end do
    ! read(*,*) i
  end subroutine

  ! *** External Circuit Initialization ***
  subroutine circ_init(Vset)
    real(8), intent(in) :: Vset

    Vmax  = Vset
    Vsrc  = Vset
    Vd_pl = Vset
    Vd_mi = Vset
    Vd_or = Vset
    Cap = 1d-11 * ph0 / e
    !Res = R0 * e / (ph0 * t0)
  end subroutine
end module
