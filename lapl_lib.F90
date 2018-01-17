module lapl_lib
#include <petsc/finclude/petscksp.h>
  use petscksp
  use props
  use ptcl_props
  implicit none

  Mat A
  Vec b
  KSP ksp
  PetscInt Istart, Iend, ii, jj

  integer, allocatable :: nn(:)
  real(8), allocatable :: ph(:,:), ph_mi(:,:), soln(:)
  logical :: assem = .True.

  public :: ph
  private :: ph_mi

contains

  subroutine lapl_init(g)
    type(grid), intent(in) :: g

    allocate(ph(g%bx+2, g%by+2), ph_mi(g%bx+2, g%by+2), &
             nn(g%bx), soln(g%bx))
    ph = 0
    ph_mi = 0

    ! Petsc Objects A and b
    call MatCreate(comm, A, ierr)
    call MatSetSizes(A, g%nloc, g%nloc, g%nglob, g%nglob, ierr)
    call MatSetUp(A, ierr)
    call MatSetType(A, MATAIJ, ierr)
    call MatSetFromOptions(A, ierr)
    call MatSeqAIJSetPreallocation(A, g%dof*5, petsc_null_integer, ierr)
    call MatMPIAIJSetPreallocation(A, g%dof*5, petsc_null_integer, &
                                   g%dof*2, petsc_null_integer, ierr)
    call MatSetOption(A, mat_ignore_zero_entries, petsc_true, ierr)

    ! Find parallel partitioning range
    call MatGetOwnershipRange(A, Istart, Iend, ierr)

    ! Create parallel vectors
    call VecCreateMPI(comm, g%nloc, g%nglob, b, ierr)
    call VecSetFromOptions(b, ierr)
    call VecSetOption(b, vec_ignore_negative_indices, petsc_true, ierr)

    ! Create Linear Solver
    call KSPCreate(comm, ksp, ierr)
    call KSPSetOperators(ksp, A, A, ierr)
    !call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr) !can be helpful
    call KSPSetTYpe(ksp, KSPIBCGS, ierr) !works well for poisson
    call KSPSetFromOptions(ksp, ierr)
  end subroutine

  subroutine lapl_solve(g, ne, ni, nt, sig, lerr)
    type(grid), intent(in) :: g
    integer, intent(out) :: lerr
    real(8), intent(in) :: ne(:,:), ni(:,:), nt(:,:), sig(:)
    integer :: i, j, cols(5), rows(1), conv
    real(8) :: b_temp(1), A_temp(5), sig0

    ! Assemble A and b
    do j = 2, g%by+1
        do i = 2, g%bx+1
          cols = -1
          rows = g%node(i,j,1)
          A_temp = 0d0
          b_temp = 0d0

          if (g%type_y(i-1,j-1) == 3) then
            sig0 = sig(i)
          else
            sig0 = 0d0
          end if

          call laplEqn(g, i, j, ne, ni, nt, sig0, b_temp(1))
          call laplJ(g, i, j, ne, ni, nt, sig0, A_temp, cols)

          ii = g%dof
          jj = g%dof*5
          call VecSetValues(b, ii, rows, -b_temp, Insert_Values, ierr)
          call MatSetValues(A, ii, rows, jj, cols, A_temp, &
                            Insert_Values, ierr)
      end do
    end do

    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)

    call MatAssemblyBegin(A, Mat_Final_Assembly, ierr)
    call MatAssemblyEnd(A, Mat_Final_Assembly, ierr)

    if (assem) then
      call MatSetOption(A, Mat_New_Nonzero_Locations, PETSc_False, ierr)
    else
      assem = .false.
    end if

    ! Solve system
    call KSPSetOperators(ksp, A, A, ierr)
    call KSPSolve(ksp, b, b, ierr)
    call KSPGetConvergedReason(ksp, conv, ierr)

    if (conv < 0) then
      if (myId == 0) write(*,*) 'PETSc KSP diverged. Rejecting time step.'
      lerr = 1
    else
      lerr = 0

      do j = 2, g%by+1
        nn = g%node(2:g%bx+1,j,1)
        call VecGetValues(b, g%bx, nn, soln, ierr)
        ph(2:g%bx+1,j) = ph(2:g%bx+1,j) + soln
      end do

      call comm_real(g%bx, g%by, ph)
      ph_mi = ph
    end if
  end subroutine

  subroutine laplEqn(g, i, j, ne, ni, nt, sig, b)
    type(grid), intent(in) :: g
    integer, intent(in)  :: i, j
    real(8), intent(in)  :: ne(:,:), ni(:,:), nt(:,:), sig
    real(8), intent(out) :: b
    real(8) :: dfdx = 0, dfdy = 0, &
               dflxe_x = 0, dflxe_y = 0, &
               dflxi_x = 0, dflxi_y = 0, &
               Ex(2), Ey(2), mue(2), De(2), Te(3), ve, a, &
               flxi(2), flxe(2)

    if (g%nx > 1) then
      ! X-dir fields:
      if (g%type_x(i-1,j-1) == -1) then
        Ex(1) = 0d0
        Ex(2) = -(ph(i+1,j) - ph(i,j)) / g%dx(i)
      else if (g%type_x(i-1,j-1) == 1) then
        Ex(1) = -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)
        Ex(2) = 0d0
      else
        Ex(1) = -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)
        Ex(2) = -(ph(i+1,j) - ph(i,j)) / g%dx(i)
      end if

      ! X-dir Fluxes:
      if (g%type_x(i-1,j-1) == 0) then
        ! rates and coefficients
        Te(1) = get_Te(nt(i-1,j), ne(i-1,j))
        Te(2) = get_Te(nt(i,j),   ne(i,j))
        Te(3) = get_Te(nt(i+1,j), ne(i+1,j))

        mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))

        De(1) = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))
        De(2) = 5d-1 * (get_De(Te(2)) + get_De(Te(3)))

        ! Flux at i - 1/2
        call get_flux(flxi(1), Ex(1), g%dx(i-1),  1, mui, Di, &
                      ni(i-1,j), ni(i,j))
        call get_flux(flxe(1), Ex(1), g%dx(i-1), -1, mue(1), De(1), &
                      ne(i-1,j), ne(i,j))

        ! Flux at i + 1/2
        call get_flux(flxi(2), Ex(2), g%dx(i),  1, mui, Di, &
                      ni(i,j), ni(i+1,j))
        call get_flux(flxe(2), Ex(2), g%dx(i), -1, mue(2), De(2), &
                      ne(i,j), ne(i+1,j))
      else if (g%type_x(i-1,j-1) < 0) then
        ! rates and coefficients
        Te(2) = get_Te(nt(i,j),   ne(i,j))
        Te(3) = get_Te(nt(i+1,j), ne(i+1,j))

        mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))
        De(2) = 5d-1 * (get_De(Te(2)) + get_De(Te(3)))

        ! Flux at i + 1/2
        call get_flux(flxi(2), Ex(2), g%dx(i),  1, mui, Di, &
                      ni(i,j), ni(i+1,j))
        call get_flux(flxe(2), Ex(2), g%dx(i), -1, mue(2), De(2), &
                      ne(i,j), ne(i+1,j))

        ! - electrode -
        if (g%type_x(i-1,j-1) == -2) then
          if (ph(i-1,j) > ph(i,j)) then
            a = 1d0 ! electrons drift
          else
            a = 0d0 ! ions drift
          end if

          mue(1) = get_mue(Te(2))
          ve = sqrt((16d0 * e * ph0 * Te(2)) / (3d0 * pi * me)) * t0 / x0

          ! Flux at i - 1/2
          flxi(1) = (1d0 - a) * mui * Ex(1) * ni(i,j) &
                    - 2.5d-1 * vi * ni(i,j)

          flxe(1) = - a * mue(1) * Ex(1) * ne(i,j) &
                    - 2.5d-1 * ve * ne(i,j) &
                    - gam * flxi(1)

        ! - vacuum -
        else if (g%type_x(i-1,j-1) == -1) then
          ! Flux at i - 1/2
          flxi(1) = 0d0
          flxe(1) = 0d0
        end if
      else if (g%type_x(i-1,j-1) > 0) then
        ! rates and coefficients
        Te(1) = get_Te(nt(i-1,j), ne(i-1,j))
        Te(2) = get_Te(nt(i,j),   ne(i,j))
        mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        De(1) = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))

        ! Flux at i - 1/2
        call get_flux(flxi(1), Ex(1), g%dx(i-1),  1, mui, Di, &
                      ni(i-1,j), ni(i,j))
        call get_flux(flxe(1), Ex(1), g%dx(i-1), -1, mue(1), De(1), &
                      ne(i-1,j), ne(i,j))

        ! - electrode -
        if (g%type_x(i-1,j-1) == 2) then
          if (ph(i+1,j) > ph(i,j)) then
              a = 1d0 ! electrons drift
          else
              a = 0d0 ! ions drift
          end if

          mue(2) = get_mue(Te(2))
          ve = sqrt((16d0 * e * ph0 * Te(2)) / (3d0 * pi * me)) * t0 / x0

          ! Flux at i + 1/2
          flxi(2) = (1d0 - a) * mui * Ex(2) * ni(i,j) &
                    + 2.5d-1 * vi * ni(i,j)

          flxe(2) = - a * mue(2) * Ex(2) * ne(i,j) &
                    + 2.5d-1 * ve * ne(i,j) &
                    - gam * flxi(2)

        ! - vacuum -
        else if (g%type_x(i-1,j-1) == 1) then
          ! Flux at i + 1/2
          flxi(2) = 0d0
          flxe(2) = 0d0
        end if
      end if

      dfdx = -(Ex(2) - Ex(1)) / g%dlx(i-1)
      dflxe_x = (flxe(2) - flxe(1)) / g%dlx(i-1)
      dflxi_x = (flxi(2) - flxi(1)) / g%dlx(i-1)
    end if

    if (g%ny > 1) then
      ! Y-dir fields:
      if (g%type_y(i-1,j-1) == -1) then
        Ey(1) = 0d0
        Ey(2) = -(ph(i,j+1) - ph(i,j)) / g%dy(j)
      else if (g%type_y(i-1,j-1) == 1) then
        Ey(1) = -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)
        Ey(2) = 0d0
      else if (g%type_y(i-1,j-1) == 3) then
        Ey(1) = -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)
        Ey(2) = -sig
      else
        Ey(1) = -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)
        Ey(2) = -(ph(i,j+1) - ph(i,j)) / g%dy(j)
      end if

      ! Y-dir Fluxes
      if (g%type_y(i-1,j-1) == 0) then
        ! rates and coefficients
        Te(1) = get_Te(nt(i,j-1), ne(i,j-1))
        Te(2) = get_Te(nt(i,j),   ne(i,j))
        Te(3) = get_Te(nt(i,j+1), ne(i,j+1))

        mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))
        De(1) = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))
        De(2) = 5d-1 * (get_De(Te(2)) + get_De(Te(3)))

        ! Flux at j - 1/2
        call get_flux(flxi(1), Ey(1), g%dy(j-1),  1, mui, Di, &
                      ni(i,j-1), ni(i,j))
        call get_flux(flxe(1), Ey(1), g%dy(j-1), -1, mue(1), De(1), &
                      ne(i,j-1), ne(i,j))

        ! Flux at j + 1/2
        call get_flux(flxi(2), Ey(2), g%dy(j),  1, mui, Di, &
                      ni(i,j), ni(i,j+1))
        call get_flux(flxe(2), Ey(2), g%dy(j), -1, mue(2), De(2), &
                      ne(i,j), ne(i,j+1))
      else if (g%type_y(i-1,j-1) < 0) then
        ! rates and coefficients
        Te(2) = get_Te(nt(i,j),   ne(i,j))
        Te(3) = get_Te(nt(i,j+1), ne(i,j+1))
        mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))
        De(2) = 5d-1 * (get_De(Te(2)) + get_De(Te(3)))

        ! Flux at j + 1/2
        call get_flux(flxi(2), Ey(2), g%dy(j),  1, mui, Di, &
                      ni(i,j), ni(i,j+1))
        call get_flux(flxe(2), Ey(2), g%dy(j), -1, mue(2), De(2), &
                      ne(i,j), ne(i,j+1))

        ! Flux at j - 1/2
        flxi(1) = 0d0
        flxe(1) = 0d0
      else if (g%type_y(i-1,j-1) > 0) then
        ! rates and coefficients
        Te(1) = get_Te(nt(i,j-1), ne(i,j-1))
        Te(2) = get_Te(nt(i,j),   ne(i,j))
        mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        De(1) = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))

        ! Flux at j - 1/2
        call get_flux(flxi(1), Ey(1), g%dy(j-1),  1, mui, Di, &
                      ni(i,j-1), ni(i,j))
        call get_flux(flxe(1), Ey(1), g%dy(j-1), -1, mue(1), De(1), &
                      ne(i,j-1), ne(i,j))

        ! Flux at j + 1/2
        if (g%type_y(i-1,j-1) == 3) then
          ! if (Ey(2) < 0d0) then
          !   a = 1d0 ! electrons drift
          ! else
          !   a = 0d0 ! ions drift
          ! end if
          !
          ! mue(2) = get_mue(Te(2))
          ! ve = sqrt((16d0 * e * ph0 * Te(2)) / (3d0 * pi * me)) * t0 / x0
          !
          ! ! Flux at j + 1/2
          ! flxi(2) = - (1d0 - a) * mui * Ey(2) * ni(i,j) &
          !           + 2.5d-1 * vi * ni(i,j)
          !
          ! flxe(2) = - a * mue(2) * Ey(2) * ne(i,j) &
          !           + 2.5d-1 * ve * ne(i,j)

          Te(2) = get_Te(nt(i,j),   ne(i,j))
          mue(2) = get_mue(Te(2))
          De(2) = get_De(Te(2))

          ! Flux at j + 1/2
          call get_flux(flxi(2), Ey(2), g%dy(j),  1, mui, Di, &
                        ni(i,j), 0d0)
          call get_flux(flxe(2), Ey(2), g%dy(j), -1, mue(2), De(2), &
                        ne(i,j), 0d0)
        else
          flxi(2) = 0d0
          flxe(2) = 0d0
        end if
      end if

      if (cyl) then
        dfdy = -(Ey(2) * (g%r(j+1) + g%r(j)) / 2d0   &
                - Ey(1) * (g%r(j) + g%r(j-1)) / 2d0) &
                / g%dly(j-1) / g%r(j)
        dflxi_y = (flxi(2) * (g%r(j+1) + g%r(j)) / 2d0   &
                  - flxi(1) * (g%r(j) + g%r(j-1)) / 2d0) &
                  / g%dly(j-1) / g%r(j)
        dflxe_y = (flxe(2) * (g%r(j+1) + g%r(j)) / 2d0   &
                  - flxe(1) * (g%r(j) + g%r(j-1)) / 2d0) &
                  / g%dly(j-1) / g%r(j)
      else
        dfdy = -(Ey(2) - Ey(1)) / g%dly(j-1)
        dflxi_y = (flxi(2) - flxi(1)) / g%dly(j-1)
        dflxe_y = (flxe(2) - flxe(1)) / g%dly(j-1)
      end if
    end if

    b = dfdx + dfdy + ni(i,j) - ne(i,j) &
        - g%dt * (dflxi_x + dflxi_y - dflxe_x - dflxe_y)
    if (isnan(b)) then
      write(*,*) 'Lapl b is NaN. Proc', myId, 'loc', i, j
    end if
  end subroutine

  subroutine laplJ(g, i, j, ne, ni, nt, sig, A, cols)
    type(grid), intent(in) :: g
    integer, intent(in)  :: i, j
    real(8), intent(in)  :: ne(:,:), ni(:,:), nt(:,:), sig
    real(8), intent(out) :: A(:)
    integer, intent(out) :: cols(:)
    integer :: k, ac
    real(8) E_pl, E_mi, dfxi_pl, dfxi_mi, dfxe_pl, dfxe_mi, &
            dfyi_pl, dfyi_mi, dfye_pl, dfye_mi, Te(2), mue, De

    k = 1
    if (g%type_x(i-1,j-1) .ge. 0) then
      Te(1) = get_Te(nt(i-1,j), ne(i-1,j))
      Te(2) = get_Te(nt(i,j), ne(i,j))

      mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
      De = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))

      E_pl = -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)
      E_mi = -(ph_mi(i,j) - ph_mi(i-1,j)) / g%dx(i-1)

      call dFdPh(dfxe_mi, E_pl, E_mi, -1, mue, De, g%dx(i-1), ne(i-1,j), ne(i,j))
      dfxe_mi = -mue / g%dx(i-1) * dfxe_mi
      call dFdPh(dfxi_mi, E_pl, E_mi,  1, mui, Di, g%dx(i-1), ni(i-1,j), ni(i,j))
      dfxi_mi =  mui / g%dx(i-1) * dfxi_mi

      cols(k) = g%node(i-1,j,1)
      A(k) = (1d0 / g%dx(i-1) + g%dt * (dfxi_mi - dfxe_mi)) / g%dlx(i-1)
      k = k + 1
    else if (g%type_x(i-1,j-1) .eq. -2) then
      if (ph(i-1,j) > ph(i,j)) then
        ac = 1 ! electrons streaming into wall
      else
        ac = 0 ! ions streaming into wall
      end if

      Te(2) = get_Te(nt(i,j), ne(i,j))
      mue = get_mue(Te(2))

      dfxi_mi = (1 - ac) * mui * ni(i,j) / g%dx(i-1)
      dfxe_mi = -ac * mue * ne(i,j) / g%dx(i-1) - gam * dfxi_mi
    else
      dfxi_mi = 0d0
      dfxe_mi = 0d0
    end if

    if (g%type_x(i-1,j-1) .le. 0) then
      Te(1) = get_Te(nt(i,j), ne(i,j))
      Te(2) = get_Te(nt(i+1,j), ne(i+1,j))

      mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
      De = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))

      E_pl = -(ph(i+1,j) - ph(i,j)) / g%dx(i)
      E_mi = -(ph_mi(i+1,j) - ph_mi(i,j)) / g%dx(i)

      call dFdPh(dfxe_pl, E_pl, E_mi, -1, mue, De, g%dx(i), ne(i,j), ne(i+1,j))
      dfxe_pl = mue / g%dx(i) * dfxe_pl
      call dFdPh(dfxi_pl, E_pl, E_mi,  1, mui, Di, g%dx(i), ni(i,j), ni(i+1,j))
      dfxi_pl = -mui / g%dx(i) * dfxi_pl

      cols(k) = g%node(i+1,j,1)
      A(k) = (1d0 / g%dx(i) - g%dt * (dfxi_pl - dfxe_pl)) / g%dlx(i-1)
      k = k + 1
    else if (g%type_x(i-1,j-1) .eq. 2) then
      if (ph(i+1,j) > ph(i,j)) then
        ac = 1 ! electrons streaming into wall
      else
        ac = 0 ! ions streaming into wall
      end if

      Te(2) = get_Te(nt(i,j), ne(i,j))
      mue = get_mue(Te(2))

      dfxi_pl = (ac - 1) * mui * ni(i,j) / g%dx(i)
      dfxe_pl = ac * mue * ne(i,j) / g%dx(i) - gam * dfxi_pl
    else
      dfxi_pl = 0d0
      dfxe_pl = 0d0
    end if

    if (g%ny > 1) then
      if (g%type_y(i-1,j-1) .ge. 0) then
        Te(1) = get_Te(nt(i,j-1), ne(i,j-1))
        Te(2) = get_Te(nt(i,j), ne(i,j))

        mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        De = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))

        E_pl = -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)
        E_mi = -(ph_mi(i,j) - ph_mi(i,j-1)) / g%dy(j-1)

        call dFdPh(dfye_mi, E_pl, E_mi, -1, mue, De, g%dy(j-1), ne(i,j-1), ne(i,j))
        dfye_mi = -mue / g%dy(j-1) * dfye_mi
        call dFdPh(dfyi_mi, E_pl, E_mi,  1, mui, Di, g%dy(j-1), ni(i,j-1), ni(i,j))
        dfyi_mi =  mui / g%dy(j-1) * dfyi_mi

        cols(k) = g%node(i,j-1,1)

        if (cyl) then
          A(k) = 5d-1 * (g%r(j) + g%r(j-1)) / g%dy(j-1) &
                 / g%dly(j-1) / g%r(j) &
                 + g%dt * (dfyi_mi - dfye_mi) * (g%r(j) + g%r(j-1)) / 2d0 &
                 / g%dly(j-1) / g%r(j)
        else
          A(k) = 1d0 / g%dy(j-1) / g%dly(j-1) &
                 + g%dt * (dfyi_mi - dfye_mi) / g%dly(j-1)
        end if

        k = k + 1
      else if (g%type_y(i-1,j-1) .eq. -2) then
        if (ph(i,j-1) > ph(i,j)) then
          ac = 1 ! electrons streaming into wall
        else
          ac = 0 ! ions streaming into wall
        end if

        Te(2) = get_Te(nt(i,j), ne(i,j))
        mue = get_mue(Te(2))

        dfyi_mi = (1 - ac) * mui * ni(i,j) / g%dy(j-1)
        dfye_mi = -ac * mue * ne(i,j) / g%dy(j-1) - gam * dfyi_mi
      else
        dfyi_mi = 0d0
        dfye_mi = 0d0
      end if

      if (g%type_y(i-1,j-1) .le. 0) then
        Te(1) = get_Te(nt(i,j), ne(i,j))
        Te(2) = get_Te(nt(i,j+1), ne(i,j+1))

        mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        De = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))

        E_pl = -(ph(i,j+1) - ph(i,j)) / g%dy(j)
        E_mi = -(ph_mi(i,j+1) - ph_mi(i,j)) / g%dy(j)

        call dFdPh(dfye_pl, E_pl, E_mi, -1, mue, De, g%dy(j), ne(i,j), ne(i,j+1))
        dfye_pl =  mue / g%dy(j) * dfye_pl
        call dFdPh(dfyi_pl, E_pl, E_mi,  1, mui, Di, g%dy(j), ni(i,j), ni(i,j+1))
        dfyi_pl = -mui / g%dy(j) * dfyi_pl

        cols(k) = g%node(i,j+1,1)

        if (cyl) then
          A(k) = 5d-1 * (g%r(j+1) + g%r(j)) / g%dy(j) &
                 / g%dly(j-1) / g%r(j) &
                 - g%dt * (dfyi_pl - dfye_pl) * (g%r(j+1) + g%r(j)) / 2d0 &
                 / g%dly(j-1) / g%r(j)
        else
          A(k) = 1d0 / g%dy(j) / g%dly(j-1) &
                 - g%dt * (dfyi_pl - dfye_pl) / g%dly(j-1)
        end if

        k = k + 1
      else if (g%type_y(i-1,j-1) .eq. 2) then
        if (ph(i,j+1) > ph(i,j)) then
          ac = 1 ! electrons streaming into wall
        else
          ac = 0 ! ions streaming into wall
        end if

        Te(2) = get_Te(nt(i,j), ne(i,j))
        mue = get_mue(Te(2))

        dfyi_pl = (ac - 1) * mui * ni(i,j) / g%dy(j)
        dfye_pl = ac * mue * ne(i,j) / g%dy(j) - gam * dfyi_pl
      else if (g%type_y(i-1,j-1) == 3) then
        if (-sig < 0d0) then
          ac = 1 ! electrons streaming into wall
        else
          ac = 0 ! ions streaming into wall
        end if

        Te(2) = get_Te(nt(i,j), ne(i,j))
        mue = get_mue(Te(2))

        dfyi_pl = 0d0!(ac - 1) * mui * ni(i,j) / g%dy(j)
        dfye_pl = 0d0!ac * mue * ne(i,j) / g%dy(j)
      else
        dfyi_pl = 0d0
        dfye_pl = 0d0
      end if
    end if

    cols(k) = g%node(i,j,1)

    A(k) = g%dt * (dfxi_pl - dfxi_mi - dfxe_pl + dfxe_mi) / g%dlx(i-1)
    if (g%type_x(i-1,j-1) .eq. -1) then
      A(k) = A(k) - 1d0 / g%dx(i-1) / g%dlx(i-1)
    else if (g%type_x(i-1,j-1) .eq. 1 ) then
      A(k) = A(k) - 1d0 / g%dx(i-1) / g%dlx(i-1)
    else
      A(k) = A(k) - 1d0 / g%dlx(i-1) * (1d0 / g%dx(i) + 1d0 / g%dx(i-1))
    end if

    if (g%ny > 1) then
      if (cyl) then
        A(k) = A(k) + g%dt * ((dfyi_pl - dfye_pl) * (g%r(j+1) + g%r(j)) / 2d0 &
               - (dfyi_mi - dfye_mi) * (g%r(j) + g%r(j-1)) / 2d0)             &
               / g%dly(j-1) / g%r(j)
      else
        A(k) = A(k) + g%dt * (dfyi_pl - dfyi_mi - dfye_pl + dfye_mi) &
               / g%dly(j-1)
      end if
      if (g%type_y(i-1,j-1) .eq. -1) then
        if (cyl) then
          A(k) = A(k) - 1d0 / g%dy(j-1) / g%dly(j-1) * (g%r(j+1) + g%r(j)) &
                 / 2d0 / g%r(j)
        else
          A(k) = A(k) - 1d0 / g%dy(j-1) / g%dly(j-1)
        end if
      else if (g%type_y(i-1,j-1) .ge. 1 ) then
        if (cyl) then
          A(k) = A(k) - 1d0 / g%dy(j-1) / g%dly(j-1) * (g%r(j) + g%r(j-1)) &
                 / 2d0 / g%r(j)
        else
          A(k) = A(k) - 1d0 / g%dy(j-1) / g%dly(j-1)
        end if
      else
        if (cyl) then
          A(k) = A(k) - 5d-1 * ((g%r(j+1) + g%r(j)) / g%dy(j) &
                 + (g%r(j) + g%r(j-1)) / g%dy(j-1)) / g%dly(j-1) / g%r(j)
        else
          A(k) = A(k) - 1d0 / g%dly(j-1) * (1d0 / g%dy(j) + 1d0 / g%dy(j-1))
        end if
      end if
    end if
  end subroutine

  subroutine dFdPh(df, E_pl, E_mi, q, mu, D, dx, n1, n2)
      real(8), intent(in)  :: E_pl, E_mi, mu, D, dx, n1, n2
      integer, intent(in)  :: q
      real(8), intent(out) :: df
      real(8) :: v_pl, v_mi, tol, a_pl, a_mi

      tol = 1d-12
      v_mi = q * mu * E_mi
      v_pl = q * mu * E_pl
      a_mi = v_mi * dx / D
      a_pl = v_pl * dx / D

      if (abs(v_mi) < tol) then
        df = 0d0
      else if (a_mi > 0) then
        df = (n1 - n2 * exp(-a_mi)) / (1d0 - exp(-a_mi)) &
             + v_mi * exp(-a_mi) * (n2 - n1) / (1d0 - exp(-a_mi))**2 * dx / D
      else
        df = (n2 - n1 * exp( a_mi)) / (1d0 - exp( a_mi)) &
             + v_mi * exp( a_mi) * (n2 - n1) / (1d0 - exp( a_mi))**2 * dx / D
      end if
      if (isnan(df)) then
        write(*,*) 'flx is NaN', v_mi, v_pl, a_mi, a_pl, n1, n2
        stop
      end if
  end subroutine
end module
