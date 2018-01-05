! *** Dielectric Surface Charge Module ***
module sfc_lib
  use props
  use ptcl_lib
  implicit none

  real(8), allocatable :: sig(:)

  ! variables
  public  :: sig
  !private ::

  ! subroutines
  public  :: sfc_step, sfc_init
  !private ::

contains

  ! *** Surface Charge Timestepping ***
  subroutine sfc_step(g, ph)
    type(grid), intent(in) :: g
    real(8), intent(inout) :: ph(:,:)
    integer :: i, j
    real(8) :: a, Te, mue, ve, Ey, fluxi, fluxe

    if (ry == py-1) then
      j = g%by+1
      do i = 2, g%bx + 1
        if (ph(i,j-1) > ph(i,j)) then
          a = 1d0 ! ions drifting
        else
          a = 0d0 ! electrons drifting
        end if

        Te  = get_Te(nt(i,j,1), ne(i,j,1))
        mue = get_mue(Te)
        ve  = sqrt((16d0 * e * ph0 * Te) / (3d0 * pi * me)) * t0 / x0
        Ey  = -sig(i)!(ph(i,j-1) - ph(i,j)) / g%dy(j-1)

        fluxi = a * mui * Ey * ni(i,j,1) + 0.25 * vi * ni(i,j,1)
        fluxe = (a - 1) * mue * Ey * ne(i,j,1) + 0.25 * ve * ne(i,j,1)
        sig(i) = sig(i) + g%dt * (fluxi - fluxe)

        ! sig(i) = (sig(i) + g%dt / 4d0 * (vi * ni(i,j,1) + ve * ne(i,j,1))) &
        !          / (1d0 + g%dt * (a * mui * ni(i,j,1) - (a - 1d0) * mue * ne(i,j,1)))
      end do

      call MPI_Send(sig(g%bx+1), 1, etype, east, 9, comm, ierr)
      call MPI_Send(sig(2),      1, etype, west, 9, comm, ierr)
      call MPI_Recv(sig(g%bx+2), 1, etype, east, 9, comm, stat, ierr)
      call MPI_Recv(sig(1),      1, etype, west, 9, comm, stat, ierr)
    end if

    call MPI_Barrier(comm, ierr)
  end subroutine

  ! *** Surface Charge Initialization ***
  subroutine sfc_init(g)
    type(grid), intent(inout) :: g

    allocate(sig(g%bx+2))

    sig = 0d0

  end subroutine
end module
