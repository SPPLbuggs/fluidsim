module props
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none

  real(8) :: t_m = 1e9

  ! mpi variables
  integer :: comm, rxComm, ryComm, myId, nproc, ierr, &
             rx, ry, px, py, north, south, east, west, &
             fh, etype, amode, info, stat(MPI_Status_Size), &
             core_array, glob_array

  integer(kind=MPI_Offset_Kind) :: dispx, dispy

  type :: grid
    integer :: nx, ny, bx, by, offx, offy, nglob, nloc, dof
    integer, allocatable :: node(:,:,:), type_x(:,:), type_y(:,:)
    real(8) :: l, w, ew, gap, dt, t
    real(8), allocatable :: dx(:), dlx(:), dy(:), dly(:), r(:)
  end type

  ! fundamental constants
  real(8), parameter:: pi     = 4d0 * atan(1d0),        &
                       eps0   = 8.85418782e-12,         &
                       mu0    = 4 * pi * 1e-7,          &
                       c0     = 1d0 / sqrt(eps0 * mu0), &
                       e      = 1.60217646e-19,         &
                       kb     = 1.3806503e-23,          &
                       me     = 9.10938188d-31

  ! non-dimensional parameters
  real(8), parameter:: x0   = 1e-4, &
                       ph0  = e / (eps0 * x0), &
                       t0   = 1e-6, &
                       wp02 = t0**2 * e**2 / (eps0 * x0**3 * me), &
                       w2   = (t0 * 14e9)**2

  ! case properties
  real(8), parameter:: Tg     = 300, & ! kelvin
                       p      = 2,   & ! torr
                       ninf   = p * 101325d0 / 760d0 / kb / Tg * x0**3, &
                       n_zero = 1e8 * x0**3
  real(8) :: n_init = 1e11 * x0**3

  integer :: rf = 0
  logical :: unif = .True., cyl = .True., rwall = .False.

contains

  ! *** Initialize Grid ***
  subroutine g_init(g, nx, ny, px, py, dof, l, w, ew, gap, path)
    type(grid), intent(inout) :: g
    real(8), intent(in) :: l, w, ew, gap
    integer, intent(in) :: nx, ny, px, py, dof
    character(*), intent(in) :: path
    integer :: i, j, d
    real(8) :: xtemp, ytemp
    real(8), allocatable :: x(:), y(:)

    g%nx = nx
    g%ny = ny

    if (g%ny == 1) cyl = .False.
    if (g%ny == 1) rwall = .False.

    ! Coordinates of processor (myId = ry * px + rx)
    rx = mod(myId, px)
    ry = myId / px

    ! Determine neighbor processors
    north = (ry - 1) * px + rx
    if (ry - 1 < 0) north = MPI_Proc_Null
    south = (ry + 1) * px + rx
    if (ry + 1 >= py) south = MPI_Proc_Null
    west = ry * px + rx - 1
    if (rx - 1 < 0) west = MPI_Proc_Null
    east = ry * px + rx + 1
    if (rx + 1 >= px) east = MPI_Proc_Null

    call MPI_Comm_split(comm, ry, myId, ryComm, ierr)
    call MPI_Comm_split(comm, rx, myId, rxComm, ierr)

    ! Domain decomposition
    g%bx = nx / px
    g%by = ny / py

    g%offx = rx * g%bx
    g%offy = ry * g%by

    g%l   = l
    g%w   = w
    g%ew  = ew
    g%gap = gap

    xtemp = 2.0 * 1.175 / float(g%nx+1)
    allocate(x(g%bx+2), g%dx(g%bx+1), g%dlx(g%bx) )

    if (unif) then
      do i = 1, g%bx+2
        x(i) = g%l / float(g%bx+1) * (g%offx + i - 1)
      end do
    else
      do i = 1, g%bx+2
        x(i) = tanh(-1.175 + xtemp * (g%offx + i - 1))
      end do
    end if

    xtemp = x(1)
    call MPI_Bcast(xtemp, 1, MPI_Real8, 0, comm, ierr)
    x = x - xtemp

    xtemp = x(g%bx+2)
    call MPI_Bcast( xtemp, 1, MPI_Real8, nproc-1, comm, ierr)
    x = x / xtemp
    x = x * g%l

    do i = 1, g%bx+1
      g%dx(i) = x(i+1) - x(i)
    end do

    do i = 1, g%bx
      g%dlx(i) = 0.5 * (g%dx(i+1) + g%dx(i))
    end do

    xtemp = minval(g%dx)
    ytemp = maxval(g%dx)
    call MPI_Allreduce(MPI_In_Place, xtemp, 1, MPI_REAL8, MPI_Min, comm, ierr)
    call MPI_Allreduce(MPI_In_Place, ytemp, 1, MPI_REAL8, MPI_Max, comm, ierr)
    if ((.not. unif) .and. (myID == 0)) &
      write(*,31) g%l / float(g%nx+1) / xtemp, ytemp / xtemp
    31 format('Non-uniform Grid, dx_u / dx_min:', f6.2, '   dx_max / dx_min:', f6.2)

    if (g%ny > 1) then
      allocate( y(g%by+2), g%dy(g%by+1), g%dly(g%by) )
      ytemp = 0.886 / float(g%ny+1)

      if (unif) then
        do j = 1, g%by+2
          y(j) = g%w / float(g%by+1) * (g%offy + j - 1)
        end do
      else
        do j = 1, g%by+2
          y(j) = tanh(-0.886 + ytemp * (g%offy + j - 1))
        end do
      end if

      ytemp = y(1)
      call MPI_Bcast( ytemp, 1, MPI_Real8, 0, comm, ierr)
      y = y - ytemp

      ytemp = y(g%by+2)
      call MPI_Bcast( ytemp, 1, MPI_Real8, nproc-1, comm, ierr)
      y = y / ytemp
      y = y * g%w

      do j = 1, g%by+1
        g%dy(j) = y(j+1) - y(j)
      end do

      do j = 1, g%by
        g%dly(j) = 0.5 * (g%dy(j+1) + g%dy(j))
      end do

      allocate(g%r(g%by+2))
      g%r = y

      xtemp = minval(g%dy)
      ytemp = maxval(g%dy)
      call MPI_Allreduce(MPI_In_Place, xtemp, 1, MPI_REAL8, MPI_Min, comm, ierr)
      call MPI_Allreduce(MPI_In_Place, ytemp, 1, MPI_REAL8, MPI_Max, comm, ierr)
      if( (.not. unif) .and. (myID == 0)) &
        write(*,32) g%w / float(g%ny+1) / xtemp, ytemp / xtemp
      32 format('                  dy_u / dy_min:', f6.2, '   dy_max / dy_min:', f6.2)
    else
      allocate(g%dy(1), y(1))
      g%dy = g%dx(1)
      y = 1
    end if

    ! Define node types
    allocate(g%type_x(g%bx,g%by), g%type_y(g%bx,g%by))
    g%type_x = 0
    g%type_y = 0

    do j = 1, g%by
      do i = 1, g%bx
        if ((ry == 0) .and. (j == 1)) then
          g%type_y(i,j) = -1
        end if

        if ((ry == py-1) .and. (j == g%by)) then
          if (rwall) then
            g%type_y(i,j) = 3
          else
            g%type_y(i,j) = 2
          end if
        end if

        if ((rx == 0) .and. (i == 1)) then
          if ((y(j) .le. g%ew) .or. &
              (y(j) .ge. g%gap)) then
            g%type_x(i,j) = -2
          else
            g%type_x(i,j) = -1
          end if
        end if

        if ((rx == px-1) .and. (i == g%bx)) then
          if ((y(j) .le. g%ew) .or. &
              (y(j) .ge. g%gap)) then
            g%type_x(i,j) = 2
          else
            g%type_x(i,j) = 1
          end if
        end if
      end do
    end do

    g%dof   = dof
    g%nloc  = g%bx * g%by * g%dof
    g%nglob = g%nx * g%ny * g%dof

    allocate( g%node(g%bx+2, g%by+2, g%dof) )
    g%node = 0
    do d = 1, g%dof
      do j = 2, g%by+1
        do i = 2, g%bx+1
          g%node(i,j,d) = g%nloc * myId + (d - 1) &
                          + ((i - 2) + (j - 2) * g%bx) * g%dof
        end do
      end do

      call comm_int(g%bx, g%by, g%node(:,:,d))
    end do

    ! MPI-IO Variables
    amode = MPI_Mode_WRonly + MPI_Mode_Create + MPI_Mode_EXCL
    etype = MPI_Real8
    stat  = MPI_Status_Ignore
    info  = MPI_Info_Null
    dispx  = rx*g%bx*8
    dispy  = ry*g%by*8

    ! Save mesh to disk
    call MPI_File_Open(comm, path//'meshx.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'meshx.dat', info, ierr);
      call MPI_File_Open(comm, path//'meshx.dat', amode,  info, fh, ierr)
    end if

    call MPI_File_Set_View(fh, dispx, etype, etype, 'native', info, ierr)
    if (ry == 0) call MPI_File_Write(fh, x(2:g%bx+1) * x0, g%bx, etype, stat, ierr)
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, path//'meshy.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'meshy.dat', info, ierr);
      call MPI_File_Open(comm, path//'meshy.dat', amode,  info, fh, ierr)
    end if

    if (g%ny > 1) then
      call MPI_File_Set_View(fh, dispy, etype, etype, 'native', info, ierr)
      if (rx == 0) call MPI_File_Write(fh, y(2:g%by+1) * x0, g%by, etype, stat, ierr)
      call MPI_File_Close(fh, ierr)
    end if

    call MPI_File_Open(comm, path//'time.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'time.dat', info, ierr);
      call MPI_File_Open(comm, path//'time.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, path//'dt.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'dt.dat', info, ierr);
      call MPI_File_Open(comm, path//'dt.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, path//'vd.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'vd.dat', info, ierr);
      call MPI_File_Open(comm, path//'vd.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, path//'ida.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'ida.dat', info, ierr);
      call MPI_File_Open(comm, path//'ida.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, path//'idc.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'idc.dat', info, ierr);
      call MPI_File_Open(comm, path//'idc.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, path//'f1.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'f1.dat', info, ierr);
      call MPI_File_Open(comm, path//'f1.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, path//'f2.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'f2.dat', info, ierr);
      call MPI_File_Open(comm, path//'f2.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, path//'f3.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'f3.dat', info, ierr);
      call MPI_File_Open(comm, path//'f3.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, path//'f4.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'f4.dat', info, ierr);
      call MPI_File_Open(comm, path//'f4.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, path//'f5.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      if (myId == 0) call MPI_File_Delete(path//'f5.dat', info, ierr);
      call MPI_File_Open(comm, path//'f5.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_Type_Create_Subarray(2, (/ g%bx+2, g%by+2 /), (/ g%bx, g%by /), &
      (/ 1, 1 /), MPI_Order_Fortran, etype, core_array, ierr)
    call MPI_Type_Commit(core_array, ierr)

    call MPI_Type_Create_Subarray(2, (/ g%nx, g%ny /), (/ g%bx, g%by /), &
      (/ rx * g%bx, ry * g%by /), MPI_Order_Fortran, etype, glob_array, ierr)
    call MPI_Type_Commit(glob_array, ierr)

    call MPI_Barrier(comm, ierr)
  end subroutine

  ! *** Communicate data across processors
  subroutine comm_real(bx, by, A)
    integer, intent(in) :: bx, by
    real(8), intent(inout) :: A(:,:)

    call MPI_Send(A(2:bx+1,2), bx, etype, north, 9, comm, ierr)
    call MPI_Send(A(2:bx+1,by+1), bx, etype, south, 9, comm, ierr)
    call MPI_Recv(A(2:bx+1,1), bx, etype, north, 9, comm, stat, ierr)
    call MPI_Recv(A(2:bx+1,by+2), bx, etype, south, 9, comm, stat, ierr)

    call MPI_Send(A(bx+1,2:by+1), by, etype, east, 9, comm, ierr)
    call MPI_Send(A(2,2:by+1), by, etype, west, 9, comm, ierr)
    call MPI_Recv(A(bx+2,2:by+1), by, etype, east, 9, comm, stat, ierr)
    call MPI_Recv(A(1,2:by+1), by, etype, west, 9, comm, stat, ierr)
  end subroutine

  ! *** Communicate data across processors
  subroutine comm_int(bx, by, A)
    integer, intent(in) :: bx, by
    integer, intent(inout) :: A(:,:)

    call MPI_Send(A(2:bx+1,2), bx, MPI_INT, north, 9, comm, ierr)
    call MPI_Send(A(2:bx+1,by+1), bx, MPI_INT, south, 9, comm, ierr)
    call MPI_Recv(A(2:bx+1,1), bx, MPI_INT, north, 9, comm, stat, ierr)
    call MPI_Recv(A(2:bx+1,by+2), bx, MPI_INT, south, 9, comm, stat, ierr)

    call MPI_Send(A(bx+1,2:by+1), by, MPI_INT, east, 9, comm, ierr)
    call MPI_Send(A(2,2:by+1), by, MPI_INT, west, 9, comm, ierr)
    call MPI_Recv(A(bx+2,2:by+1), by, MPI_INT, east, 9, comm, stat, ierr)
    call MPI_Recv(A(1,2:by+1), by, MPI_INT, west, 9, comm, stat, ierr)
  end subroutine

  ! *** Save Data ***
  subroutine savedat(path,dat)
    character(*), intent(in) :: path
    real(8), intent(in) :: dat(:,:)
    integer (kind = MPI_Offset_Kind) :: offset

    call MPI_File_Open(comm, path, MPI_MODE_RDWR + MPI_MODE_APPEND, info, fh, ierr)
    call MPI_File_Get_Position(fh, offset, ierr)
    call MPI_File_Set_View(fh, offset, etype, glob_array, 'native', info, ierr)
    call MPI_File_Write_All(fh, dat, 1, core_array, stat, ierr)
    call MPI_File_Close(fh, ierr)
  end subroutine
end module
