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
  subroutine sfc_step(g)
    type(grid), intent(in) :: g
    integer :: i, j
    real(8) :: Te, De, mue, fluxi, fluxe

    if (ry == py-1) then
      j = g%by+1
      do i = 2, g%bx + 1
        Te = get_Te(nt(i,j,1),   ne(i,j,1))
        mue = get_mue(Te)
        De = get_De(Te)

        ! Flux at j + 1/2
        call get_flux(fluxi, -sig(i), g%dy(j),  1, mui, Di, &
                      ni(i,j,1), 0d0)
        call get_flux(fluxe, -sig(i), g%dy(j), -1, mue, De, &
                      ne(i,j,1), 0d0)

        sig(i) = sig(i) + g%dt * (fluxi - fluxe)
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
