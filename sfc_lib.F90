! *** Dielectric Surface Charge Module ***
module sfc_lib
  use props
  use ptcl_props
  implicit none

  real(8), allocatable :: sig_pl(:), sig_mi(:), sig_or(:), ks(:,:)

  ! variables
  public  :: sig_mi, sig_pl, sig_or
  !private ::

  ! subroutines
  public  :: sfc_step, sfc_init
  !private ::

contains

  ! *** Surface Charge Timestepping ***
  subroutine sfc_step(g, stage, ne, ni, nt)
    type(grid), intent(in) :: g
    integer, intent(in) :: stage
    real(8), intent(in) :: ne(:,:), ni(:,:), nt(:,:)
    integer :: i, j
    real(8) :: Te, De, mue, fluxi, fluxe

    j = g%by+1
    do i = 2, g%bx + 1
      if (g%type_y(i-1,j) == 3) then
        Te = get_Te(nt(i,j),   ne(i,j))
        mue = get_mue(Te)
        De = get_De(Te)

        ! Flux at j + 1/2
        call get_flux(fluxi, -sig_mi(i), g%dy(j),  1, mui, Di, &
                      ni(i,j), 0d0)
        call get_flux(fluxe, -sig_mi(i), g%dy(j), -1, mue, De, &
                      ne(i,j), 0d0)

        ks(i,stage) = fluxi - fluxe
      end if
    end do

    do i = 2, g%bx + 1
      if (g%type_y(i-1,j) == 3) then
        if (stage == 1) then
          sig_mi(i) = sig_or(i) + g%dt * ks(i,1) / 3d0

        else if (stage == 2) then
          sig_mi(i) = sig_or(i) + g%dt * ( ks(i,1) + ks(i,2)) / 6d0

        else if (stage == 3) then
          sig_mi(i) = sig_or(i) + g%dt * (ks(i,1) + ks(i,3) * 3d0) / 8d0

        else if (stage == 4) then
          sig_mi(i) = sig_or(i) + g%dt * (ks(i,1) &
                      - ks(i,3) * 3d0 + ks(i,4) * 4d0 ) / 2d0

        else
          sig_pl(i) = sig_or(i) + g%dt * (ks(i,1) &
                      + ks(i,4) * 4d0 + ks(i,5)) / 6d0
        end if
      end if
    end do

    if (stage == 5) then
      call MPI_Send(sig_pl(g%bx+1), 1, etype, east, 9, comm, ierr)
      call MPI_Send(sig_pl(2),      1, etype, west, 9, comm, ierr)
      call MPI_Recv(sig_pl(g%bx+2), 1, etype, east, 9, comm, stat, ierr)
      call MPI_Recv(sig_pl(1),      1, etype, west, 9, comm, stat, ierr)
    else
      call MPI_Send(sig_mi(g%bx+1), 1, etype, east, 9, comm, ierr)
      call MPI_Send(sig_mi(2),      1, etype, west, 9, comm, ierr)
      call MPI_Recv(sig_mi(g%bx+2), 1, etype, east, 9, comm, stat, ierr)
      call MPI_Recv(sig_mi(1),      1, etype, west, 9, comm, stat, ierr)
    end if

    call MPI_Barrier(comm, ierr)
  end subroutine

  ! *** Surface Charge Initialization ***
  subroutine sfc_init(g)
    type(grid), intent(inout) :: g

    allocate(sig_pl(g%bx+2), sig_mi(g%bx+2), sig_or(g%bx+2), ks(g%bx+2, 5))

    sig_pl = 0d0
    sig_mi = 0d0
    sig_or = 0d0

  end subroutine
end module
