    module mdl
    use mpi
    implicit none

    ! MPI Variables:
    integer :: comm, my_id, nproc, ierr, px, py, rx, ry, offx, offy, &
               north, south, east, west, &
               fh, etype, amode, info, core_array, glob_array

    ! Constants
    real(8), parameter:: pi   = 4d0 * atan(1d0),        &
                         eps0 = 8.85418782e-12,         &
                         mu0  = pi * 4e-7,              &
                         c0   = 1d0 / sqrt(eps0 * mu0), &
                         e    = 1.60217646e-19,         &
                         kb   = 1.3806503e-23,          &
                         me   = 9.10938356e-31,         &
                         a    = 0.99 / sqrt(2.0)

    real(8), allocatable    :: Hz(:,:), Ey(:,:), Ex(:,:), &
                               Jx(:,:), Jy(:,:), wp(:,:), eps(:,:), &
                               lbc(:,:,:), rbc(:,:,:), tbc(:,:,:), bbc(:,:,:)
    complex(8), allocatable :: a_ey(:), a_hz(:)

    integer :: nx, ny, nt, nf
    real(8) :: dx, dt, nu, np
    logical :: plt_flds = .False.

    public  :: plt_flds
    private :: px, py, nx, ny, nt, nf, np, pi, eps0, mu0, c0, e, kb, me, a, Hz, Ey, Ex, &
               Jx, Jy, dx, dt, a_ey, a_hz, wp, eps, lbc, rbc, tbc, bbc, nu

    contains

! *** Initialization ***
    subroutine init(comm0, px0, py0, wp0, eps0, nt0, np0, l, nx0, ny0)
    integer, intent(in) :: comm0, px0, py0, nx0, ny0, nt0
    real, intent(in) :: l
    real(8), intent(in) :: np0, wp0(nx0,ny0), eps0(nx0,ny0)

    ! Inputs:
    !              px = # of processors in x-dir
    !              py = # of processors in y-dir
    !     wp2(nx, ny) = input array of square of plasma frequency
    !     eps(nx, ny) = input array of dielectric constant (metal is -inf)
    !     nt = # points in time
    !     np = # points in source wavelet
    !      l = length in x-dir in meters

    comm = comm0
    call MPI_Comm_rank(comm, my_id, ierr)
    call MPI_Comm_size(comm, nproc, ierr)

    px = px0
    py = py0
    nx = nx0
    ny = ny0
    nt = nt0
    np = np0

    ! Coordinates of processor (my_id = ry * px + rx)
    rx = mod(my_id, px)
    ry = my_id / px

    ! Determine neighbor processors
    north = (ry - 1) * px + rx
    if (ry - 1 < 0) north = MPI_Proc_Null
    south = (ry + 1) * px + rx
    if (ry + 1 >= py) south = MPI_Proc_Null
    west = ry * px + rx - 1
    if (rx - 1 < 0) west = MPI_Proc_Null
    east = ry * px + rx + 1
    if (rx + 1 >= px) east = MPI_Proc_Null

    ! Domain decomposition
    offx = rx * nx
    offy = ry * ny

    allocate(Hz(nx+2, ny+2), Ex(nx+2, ny+2), Ey(nx+2, ny+2), & !Ey(nx+1, ny+2)
             Jx(nx+2, ny+1), Jy(nx+1, ny+2), wp(nx+2,ny+2), eps(nx+2,ny+2), &
             lbc(3,ny+2,2), rbc(3,ny+2,2), tbc(nx+2,3,2), bbc(nx+2,3,2))

    Hz = 0
    Ex = 0
    Ey = 0
    Jx = 0
    Jy = 0
    lbc = 0
    rbc = 0
    tbc = 0
    bbc = 0

    wp = 0
    wp(2:nx+1, 2:ny+1)  = wp0**2
    eps = 1
    eps(2:nx+1, 2:ny+1) = eps0
    nu  = maxval(wp)/10.

    dx = l / float(nx-1)
    dt = dx / c0 / sqrt(2.0)

    ! MPI-IO Variables
    amode = MPI_Mode_WRonly + MPI_Mode_Create + MPI_Mode_EXCL
    etype = MPI_Real8
    info  = MPI_Info_Null

    call comm_real(nx, ny, wp)
    call comm_real(nx, ny, eps)

    call MPI_File_Open(comm, 'output/Hz.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/Hz.dat', info, ierr);
        call MPI_File_Open(comm, 'output/Hz.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, 'output/Ey.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/Ey.dat', info, ierr);
        call MPI_File_Open(comm, 'output/Ey.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, 'output/Ex.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/Ex.dat', info, ierr);
        call MPI_File_Open(comm, 'output/Ex.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)

    call MPI_Type_Create_Subarray(2, (/ nx+2, ny+2 /), (/ nx, ny /), &
        (/ 1, 1 /), MPI_Order_Fortran, etype, core_array, ierr)
    call MPI_Type_Commit(core_array, ierr)

    call MPI_Type_Create_Subarray(2, (/ nx * px, ny * py /), (/ nx, ny /), &
        (/ rx * nx, ry * ny /), MPI_Order_Fortran, etype, glob_array, ierr)
    call MPI_Type_Commit(glob_array, ierr)

    call MPI_Barrier(comm, ierr)
    end subroutine

! *** Update Wp and Eps ***
    subroutine set(wp0, eps0, nx0, ny0)
    integer, intent(in) :: nx0, ny0
    real(8), intent(in) :: wp0(nx0,ny0), eps0(nx0,ny0)

    ! Inputs:
    !     wp2(nx, ny) = input array of square of plasma frequency
    !     eps(nx, ny) = input array of dielectric constant (metal is -inf)

    wp(2:nx+1, 2:ny+1)  = wp0**2
    eps(2:nx+1, 2:ny+1) = eps0
    nu  = maxval(wp)/10.

    call comm_real(nx, ny, wp)
    call comm_real(nx, ny, eps)
    end subroutine

! *** Run ***
    subroutine run(frq, flx, fskip, nf0)
    integer, intent(in)  :: fskip, nf0
    real(8), intent(inout) :: frq(nf0), flx(nf0)
    integer :: t, k
    complex(8) :: etemp, htemp

    ! Inputs: frq(nf) = output array of frequency points
    !         flx(nf) = output array of fluxes
    !           fskip = # frequency points to skip

    nf = nf0
    allocate(a_ey(nf), a_hz(nf))
    a_ey = 0
    a_hz = 0

    Hz = 0
    Ex = 0
    Ey = 0
    Jx = 0
    Jy = 0

    do t = 1, nt
        call step(t, fskip)
        if ((mod(t,5) == 0) .and. (plt_flds)) call savedat('output/Hz.dat', Hz)
        if ((mod(t,5) == 0) .and. (plt_flds)) call savedat('output/Ey.dat', Ey)
        if ((mod(t,5) == 0) .and. (plt_flds)) call savedat('output/Ex.dat', Ex)
        if ((mod(t,1000) == 0) .and. (my_id == 0)) &
             write(*,1) t, int(ceiling(nt/1000.)*1000)
    end do
    1 format('        ', i0, ' out of ', i0, ' steps')

    do k = 1, nf
!         frq(k) = float(k + fskip)  / (nt * dt * 100.)
        frq(k) = float(k + fskip)  / (nt * dt)
        call MPI_Reduce(a_ey(k), etemp, 1, MPI_Complex16, MPI_SUM, 0, comm, ierr)
        call MPI_Reduce(a_hz(k), htemp, 1, MPI_Complex16, MPI_SUM, 0, comm, ierr)
        flx(k) = 2 * RealPart(etemp * conjg(htemp)) / float(nt)
    end do

    deallocate(a_ey, a_hz)
    end subroutine

! *** Timestepping Routine ***
    subroutine step(t, fskip)
    integer, intent(in) :: t, fskip
    integer :: i, j, k, j1, j2
    real(8) :: arg, md = 1
    complex :: ii = (0,1)
    ! .H0 |E0 .H1 |E1 .H2 |E2 ... |En-2 .Hn-2 |En-1 .Hn-1 |En .Hn+1
    ! H ranges 0 -> n+1, with n values and 2 boundary conditions
    ! E ranges from 0 -> n for n+1 values and no boundary conditions

    ! Ex and Ey Update (whole domain):
    do i = 1, nx+1
        Ey(i,:) = Ey(i,:) - 2. * a * (Hz(i+1,:) - Hz(i,:) - Jy(i,:)) &
                / (eps(i+1,:) + eps(i,:))
    end do
    do j = 1, ny+1
        Ex(:,j) = Ex(:,j) + 2. * a * (Hz(:,j+1) - Hz(:,j) + Jx(:,j)) &
                / (eps(:,j+1) + eps(:,j))
    end do

    ! Hz wavelet source:
    arg = pi**2 * (float(t) / np - md)**2
    j1 = 1
    j2 = ny+2
    if (ry == 0) j1 = 5
    if (ry == py-1) j2 = ny-3
    if (plt_flds) then
      Hz(4,j1:j2) = Hz(4,j1:j2) + sin(2d0 * pi * float(t) / np)
    else
      if (rx == 0) Hz(4,j1:j2) = Hz(4,j1:j2) + (1d0 - 2d0 * arg) * exp(-arg)
    end if

!     Hz(nx/2,ny/2) = Hz(nx/2,ny/2) + (1d0 - 2d0 * arg) * exp(-arg)
    ! Hz update (whole domain):
    do j = 2, ny+1
        do i = 2, nx+1
            Hz(i,j) = Hz(i,j) + a * (Ex(i,j) - Ex(i,j-1) &
                                    - Ey(i,j) + Ey(i-1,j))
        end do
    end do

    ! ABC on Hz:
    ! - top -
    if ((ry == 0) .and. (.True.)) then
        Hz(:,1) = -1. / (a + 1./a + 2.) * ((a + 1./a - 2.) * (Hz(:,3) + tbc(:,1,2)) &
                  +2. * (a - 1./a) * (tbc(:,1,1) + tbc(:,3,1) - Hz(:,2) - tbc(:,2,2)) &
                  -4. * (a + 1./a) * tbc(:,2,1)) - tbc(:,3,2)
        tbc(:,:,2) = tbc(:,:,1)
        tbc(:,:,1) = Hz(:,1:3)
    end if
    ! - bottom -
    if ((ry == py-1) .and. (.True.)) then
        Hz(:,ny+2) = -1. / (a + 1./a + 2.) * ((a + 1./a - 2.) * (Hz(:,ny) + bbc(:,3,2)) &
                     +2. * (a - 1./a) * (bbc(:,3,1) + bbc(:,1,1) - Hz(:,ny+1) - bbc(:,2,2)) &
                     -4. * (a + 1./a) * bbc(:,2,1)) - bbc(:,1,2)
        bbc(:,:,2) = bbc(:,:,1)
        bbc(:,:,1) = Hz(:,ny:ny+2)
    end if
    ! - left -
    if (rx == 0) then
        Hz(1,:) = -1. / (a + 1./a + 2.) * ((a + 1./a - 2.) * (Hz(3,:) + lbc(1,:,2)) &
                  +2. * (a - 1./a) * (lbc(1,:,1) + lbc(3,:,1) - Hz(2,:) - lbc(2,:,2)) &
                  -4. * (a + 1./a) * lbc(2,:,1)) - lbc(3,:,2)
        lbc(:,:,2) = lbc(:,:,1)
        lbc(:,:,1) = Hz(1:3,:)
    end if
    ! - right -
    if (rx == px-1) then
        Hz(nx+2,:) = -1. / (a + 1./a + 2.) * ((a + 1./a - 2.) * (Hz(nx,:) + rbc(3,:,2)) &
                     +2. * (a - 1./a) * (rbc(3,:,1) + rbc(1,:,1) - Hz(nx+1,:) - rbc(2,:,2)) &
                     -4. * (a + 1./a) * rbc(2,:,1)) - rbc(1,:,2)
        rbc(:,:,2) = rbc(:,:,1)
        rbc(:,:,1) = Hz(nx:nx+2,:)
    end if

    ! Perfectly conducting (E = 0) boundaryies:
    ! x-dir
    !Hz(1,:) = Hz(1,:) - a * Ey(1,:)
    !Hz(nx,:) = Hz(nx,:) + a * Ey(nx-1,:)

    ! y-dir
    !Hz(:,1) = Hz(:,1) + a * Ex(:,1)
    !Hz(:,ny) = Hz(:,ny) - a * Ex(:,ny-1)

    ! Jx and Jy Update (whole domain):
    do i = 1, nx+1
        Jy(i,:) = (Jy(i,:) - 0.5 * a * (wp(i+1,:) + wp(i,:)) * Ey(i,:)) &
                  * (1-nu)/(1+nu)
    end do
    do j = 1, ny+1
        Jx(:,j) = (Jx(:,j) - 0.5 * a * (wp(:,j+1) + wp(:,j)) * Ex(:,j)) &
                  * (1-nu)/(1+nu)
    end do

    call comm_real(nx, ny, Hz)

    ! Update Fourier coefficients
!     if ((rx == px-1) .and. (mod(t,100) == 0)) then
!       if (mod(t,100) == 0) then
        do k = 1, nf
            do j = 2, ny+1
                a_ey(k) = a_ey(k) + 0.5 * (Ey(nx, j) + Ey(nx+1, j)) &
                            * exp( -ii * (k+fskip) * 2 * pi * t / float(nt))
                a_hz(k) = a_hz(k) + Hz(nx+1, j) &
                            * exp( -ii * (k+fskip) * 2 * pi * t / float(nt))
            end do
        end do
!     end if
    end subroutine

! *** Communicate data across processors
  subroutine comm_real(bx, by, A)
    integer, intent(in) :: bx, by
    real(8), intent(inout) :: A(:,:)

    call MPI_Send(A(2:bx+1,2), bx, etype, north, 9, comm, ierr)
    call MPI_Send(A(2:bx+1,by+1), bx, etype, south, 9, comm, ierr)
    call MPI_Recv(A(2:bx+1,1), bx, etype, north, 9, comm, MPI_Status_Ignore, ierr)
    call MPI_Recv(A(2:bx+1,by+2), bx, etype, south, 9, comm, MPI_Status_Ignore, ierr)

    call MPI_Send(A(bx+1,2:by+1), by, etype, east, 9, comm, ierr)
    call MPI_Send(A(2,2:by+1), by, etype, west, 9, comm, ierr)
    call MPI_Recv(A(bx+2,2:by+1), by, etype, east, 9, comm, MPI_Status_Ignore, ierr)
    call MPI_Recv(A(1,2:by+1), by, etype, west, 9, comm, MPI_Status_Ignore, ierr)
  end subroutine

  ! *** Save Data ***
  subroutine savedat(path,dat)
    character(*), intent(in) :: path
    real(8), intent(in) :: dat(:,:)
    integer (kind = MPI_Offset_Kind) :: offset

    call MPI_File_Open(comm, path, MPI_MODE_RDWR + MPI_MODE_APPEND, info, fh, ierr)
    call MPI_File_Get_Position(fh, offset, ierr)
    call MPI_File_Set_View(fh, offset, etype, glob_array, 'native', info, ierr)
    call MPI_File_Write_All(fh, dat, 1, core_array, MPI_Status_Ignore, ierr)
    call MPI_File_Close(fh, ierr)
  end subroutine
end module
