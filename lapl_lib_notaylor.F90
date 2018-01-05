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

  public  :: ph
  private :: ph_mi

contains

  subroutine lapl_init(g)
    type(grid), intent(in) :: g
    integer :: i, j, cols(5), rows(1)
    real(8) :: A_temp(5)

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

    do j = 2, g%by+1
        do i = 2, g%bx+1
          cols = -1
          rows = g%node(i,j,1)
          A_temp = 0d0

          call laplJ(g, i, j, A_temp, cols)

          ii = g%dof
          jj = g%dof*5
          call MatSetValues(A, ii, rows, jj, cols, A_temp, &
                            Insert_Values, ierr)
      end do
    end do

    call MatAssemblyBegin(A, Mat_Final_Assembly, ierr)
    call MatAssemblyEnd(A, Mat_Final_Assembly, ierr)

    call MatSetOption(A, Mat_New_Nonzero_Locations, PETSc_False, ierr)
    call KSPSetOperators(ksp, A, A, ierr)
  end subroutine

  subroutine lapl_solve(g, ne, ni, sig)
    type(grid), intent(in) :: g
    real(8), intent(in) :: ne(:,:), ni(:,:), sig(:)
    integer :: i, j, rows(1), conv
    real(8) :: b_temp(1), sig0

    ! Assemble A and b
    do j = 2, g%by+1
        do i = 2, g%bx+1
          rows = g%node(i,j,1)
          b_temp = 0d0

          if (g%type_y(i-1,j-1) == 3) then
            sig0 = sig(i)
          else
            sig0 = 0d0
          end if

          call laplEqn(g, i, j, ne(i,j), ni(i,j), sig0, b_temp(1))

          ii = g%dof
          jj = g%dof*5
          call VecSetValues(b, ii, rows, -b_temp, Insert_Values, ierr)
      end do
    end do

    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)

    ! Solve system
    call KSPSolve(ksp, b, b, ierr)
    call KSPGetConvergedReason(ksp, conv, ierr)
    if ((myId == 0) .and. (conv < 0)) then
      write(*,*) 'PETSc KSP diverged. Stopping'
      call MPI_Abort(comm, 1, ierr)
      stop
    end if

    ! Update Solution
    ! do j = 2, g%by+1
    !   do i = 2, g%bx+1
    !     nn = g%node(i,j,1)
    !     call VecGetValues(b, 1, nn, b_temp, ierr)
    !     ph(i,j) = ph(i,j) + b_temp(1)
    !   end do
    ! end do
    do j = 2, g%by+1
      nn = g%node(2:g%bx+1,j,1)
      call VecGetValues(b, g%bx, nn, soln, ierr)
      ph(2:g%bx+1,j) = ph(2:g%bx+1,j) + soln
    end do

    call comm_real(g%bx, g%by, ph)
    ph_mi = ph
  end subroutine

  subroutine laplEqn(g, i, j, ne, ni, sig, b)
    type(grid), intent(in) :: g
    integer, intent(in)  :: i, j
    real(8), intent(in)  :: ne, ni, sig
    real(8), intent(out) :: b
    real(8) :: dfdx = 0, dfdy = 0

    if (g%nx > 1) then
      ! Left symmetry boundary
      if (g%type_x(i-1,j-1) == -1) then
        dfdx = (ph(i+1,j) - ph(i,j)) / g%dx(i) / g%dlx(i-1)

      ! Right symmetry boundary
      else if (g%type_x(i-1,j-1) == 1) then
        dfdx = (ph(i-1,j) - ph(i,j)) / g%dx(i-1) / g%dlx(i-1)

      ! Domain center
      else
        dfdx = ((ph(i+1,j) - ph(i,j)) / g%dx(i) &
               -(ph(i,j) - ph(i-1,j)) / g%dx(i-1)) &
               / g%dlx(i-1)
     end if
    end if

    if (g%ny > 1) then
      ! Left symmetry boundary
      if (g%type_y(i-1,j-1) == -1) then
        if (cyl) then
          dfdy = (ph(i,j+1) - ph(i,j)) / g%dy(j) * (g%r(j+1) + g%r(j)) &
                 / 2d0 / g%dly(j-1) / g%r(j)
        else
          dfdy = (ph(i,j+1) - ph(i,j)) / g%dy(j) / g%dly(j-1)
        end if

      ! Right symmetry boundary
      else if (g%type_y(i-1,j-1) == 1) then
        if (cyl) then
          dfdy = (ph(i,j-1) - ph(i,j)) / g%dy(j-1) * (g%r(j) + g%r(j-1)) &
                 / 2d0 / g%dly(j-1) / g%r(j)
        else
          dfdy = (ph(i,j-1) - ph(i,j)) / g%dy(j-1) / g%dly(j-1)
        end if

      ! Right dielectric wall
      else if (g%type_y(i-1,j-1) == 3) then
        if (cyl) then
          dfdy = -(sig * (g%r(j+1) + g%r(j))     &
                 + (ph(i,j) - ph(i,j-1)) / g%dy(j-1) * (g%r(j) + g%r(j-1))) &
                 / 2d0 / g%dly(j-1) / g%r(j)
        else
          dfdy = -(sig &
                 + (ph(i,j) - ph(i,j-1)) / g%dy(j-1)) &
                 / g%dly(j-1)
        end if

      ! Domain center
      else
        if (cyl) then
          dfdy = ((ph(i,j+1) - ph(i,j)) / g%dy(j) * (g%r(j+1) + g%r(j))     &
                 - (ph(i,j) - ph(i,j-1)) / g%dy(j-1) * (g%r(j) + g%r(j-1))) &
                 / 2d0 / g%dly(j-1) / g%r(j)
        else
          dfdy = ((ph(i,j+1) - ph(i,j)) / g%dy(j) &
                 -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)) &
                 / g%dly(j-1)
        end if
     end if
    end if

    b = dfdx + dfdy + ni - ne
  end subroutine

  subroutine laplJ(g, i, j, A, cols)
    type(grid), intent(in) :: g
    integer, intent(in)  :: i, j
    real(8), intent(out) :: A(:)
    integer, intent(out) :: cols(:)
    integer :: k

    k = 1
    if (g%type_x(i-1,j-1) .ge. 0) then
      cols(k) = g%node(i-1,j,1)
      A(k) = 1d0 / g%dx(i-1) / g%dlx(i-1)
      k = k + 1
    end if

    if (g%type_x(i-1,j-1) .le. 0) then
      cols(k) = g%node(i+1,j,1)
      A(k) = 1d0 / g%dx(i) / g%dlx(i-1)
      k = k + 1
    end if

    if (g%ny > 1) then
      if (g%type_y(i-1,j-1) .ge. 0) then
        cols(k) = g%node(i,j-1,1)

        if (cyl) then
          A(k) = 5d-1 * (g%r(j) + g%r(j-1)) / g%dy(j-1) &
                 / g%dly(j-1) / g%r(j)
        else
          A(k) = 1d0 / g%dy(j-1) / g%dly(j-1)
        end if

        k = k + 1
      end if

      if (g%type_y(i-1,j-1) .le. 0) then
        cols(k) = g%node(i,j+1,1)

        if (cyl) then
          A(k) = 5d-1 * (g%r(j+1) + g%r(j)) / g%dy(j) &
                 / g%dly(j-1) / g%r(j)
        else
          A(k) = 1d0 / g%dy(j) / g%dly(j-1)
        end if

        k = k + 1
      end if
    end if

    cols(k) = g%node(i,j,1)

    A(k) = 0d0
    if (g%type_x(i-1,j-1) .eq. -1) then
      A(k) = A(k) - 1d0 / g%dx(i-1) / g%dlx(i-1)
    else if (g%type_x(i-1,j-1) .eq. 1 ) then
      A(k) = A(k) - 1d0 / g%dx(i-1) / g%dlx(i-1)
    else
      A(k) = A(k) - 1d0 / g%dlx(i-1) * (1d0 / g%dx(i) + 1d0 / g%dx(i-1))
    end if

    if (g%ny > 1) then
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
end module
