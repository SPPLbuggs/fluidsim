program main
  use props
  use lapl_lib
  use ptcl_lib
  use circ_lib
  use sfc_lib
  implicit none

  type(grid) :: g
  integer :: ts = 0, ts1, ts2, nx, ny, dof, lerr, rIdx = 1
  real(8) :: l, w, ew, vl, dt, t_fin, t_pr, t_sv, t_sv0, &
             sim_start, time1, time2, t_rk = 0, t_ph = 0, t1, t2, Id_mi = 0, Vd_mi2 = 0
  character(80):: path

  ! Initialize PETSc and MPI
  call PetscInitialize(petsc_null_character, ierr)
  comm = PETSC_COMM_WORLD
  call MPI_Comm_rank(comm, myId, ierr)
  call MPI_Comm_size(comm, nproc, ierr)

  call cpu_time(sim_start)
  call cpu_time(time1)
  ts1 = 0

  ! Default properties
  nx = 100
  ny = 1
  px = 1
  py = 1
  dof = 1
  l  = 1e-2   / x0
  w  = 1.5e-2 / x0
  ew = 2e-2   / x0
  dt = 5e-6
  t_fin = 1000
  t_pr = 0d0
  t_sv = 1e-2
  t_sv0 = 1e-2
  vl = 500 / ph0
  res = 1e7

  ! Read input arguments
  call read_in
  if (rf) then
    rwall = .False.
    w = 1e-2 / x0
    ew = 1.25e-3 / x0
  end if

  ! rwall = .false.
  ! Initialize grid and arrays
  path = 'Output/'
  call g_init(g, nx, ny, px, py, dof, l, w, ew, trim(path))
  call lapl_init(g)
  call ptcl_init(g)
  call circ_init(vl)
  call sfc_init(g)

  g%t  = 0
  g%dt = dt

  ! vl = 350 / ph0
  ! ph(1,:) = vl

  do
    if (g%t >= t_fin) exit

    ! Solve ne system
    t_m = 1e9

    ! if (myID == 0) call cpu_time(t1)
    if (ts > 1) call ptcl_step(g, ph)
    ! if (myID == 0) call cpu_time(t2)
    t_rk = t_rk + t2 - t1

    ! if (.not. rf) then
    !   if (rIdx == 1) then
    !     res = 6.3e6
    !   else if (rIdx == 2) then
    !     res = 2e6
    !   else if (rIdx == 3) then
    !     res = 6.3e5
    !   else if (rIdx == 4) then
    !     res = 2e5
    !   else if (rIdx == 5) then
    !     res = 6.3e4
    !   else if (rIdx == 6) then
    !     res = 2e4
    !   else if (rIdx == 7) then
    !     res = 6.3e3
    !   else if (rIdx == 8) then
    !     res = 2e3
    !   else
    !     exit
    !   end if
    ! end if

    !call circ_step(g, rf, ph, res)
    !call circ_step(g, 1, ph, ni(:,:,2), ne(:,:,2), nt(:,:,2))

    ! Solve surface charge system
    ! if (rwall) call sfc_step(g)

    ! Solve ph system
    ! if (myID == 0) call cpu_time(t1)
    call lapl_solve(g, ne(:,:,1), ni(:,:,1), nt(:,:,1), sig_pl, lerr)
    ! if (myID == 0) call cpu_time(t2)
    t_ph = t_ph + t2 - t1

    ! Accept step
    if (lerr == 0) then
      ts = ts + 1
      g%t = g%t + g%dt

      ni(:,:,3) = ni(:,:,2)
      ne(:,:,3) = ne(:,:,2)
      nm(:,:,3) = nm(:,:,2)
      nt(:,:,3) = nt(:,:,2)

      ni(:,:,2) = ni(:,:,1)
      ne(:,:,2) = ne(:,:,1)
      nm(:,:,2) = nm(:,:,1)
      nt(:,:,2) = nt(:,:,1)

      Vd_or = Vd_mi
      Vd_mi = Vd_pl

      sig_or = sig_mi
      sig_mi = sig_pl

    ! Reject step
    else
      g%dt = g%dt / 2d0

      ni(:,:,1) = ni(:,:,3)
      ne(:,:,1) = ne(:,:,3)
      nm(:,:,1) = nm(:,:,3)
      nt(:,:,1) = nt(:,:,3)

      ni(:,:,2) = ni(:,:,3)
      ne(:,:,2) = ne(:,:,3)
      nm(:,:,2) = nm(:,:,3)
      nt(:,:,2) = nt(:,:,3)

      Vd_pl = Vd_or
      Vd_mi = Vd_or

      sig_pl = sig_or
      sig_mi = sig_or
    end if

    ! Print out some information
    if ((t_pr >= 1.5d0) .and. (myId == 0)) then
      call cpu_time(time2)
      ts2 = ts
      write(*,*)
      write(*,11) ts, g%t, (time2 - time1) / g%dt / float(ts2-ts1) / 60.
      write(*,12) g%dt, t_m
      ! write(*,13) t_rk / (time2 - time1), t_ph / (time2 - time1)
      write(*,14) Vd_pl * ph0, Id * e / t0, abs((Vd_pl - Vd_mi2) / Vd_pl)
      t_pr = 0
      time1 = time2
      ts1 = ts
      t_rk = 0
      t_ph = 0

    else if (mod(ts,50) == 0) then
      call cpu_time(time2)
      t_pr = time2 - time1
    end if

    ! Save data
    if (t_sv <= g%t) then
      call writeOut

      t_sv  = t_sv + t_sv0
      t_sv0 = t_sv0 * 1.01
      t_sv0 = min(t_sv0, 1d-1)
    end if

    if ((.not. rf) .and. (g%t > 2d0)) then
      if ((abs((Id - Id_mi) / Id) .le. (g%dt * 1d-3)) .and. &
          (abs((Vd_pl - Vd_mi2) / Vd_pl) .le. (g%dt * 1d-3))) rIdx = rIdx + 1
      Id_mi = Id
      Vd_mi2 = Vd_pl
    end if
    call MPI_Allreduce(MPI_In_Place, rIdx, 1, MPI_Int, MPI_Max, comm, ierr)
    if (rIdx > 1) exit
  end do

  if (myId == 0) then
      call cpu_time(time1)
      write(*,*)
      write(*,9) int(time1 - sim_start) / 3600, &
                 mod(int(time1 - sim_start)/60,60)
      write(*,*)
  end if

  call KSPDestroy(ksp,ierr)
  call VecDestroy(b,ierr)
  call MatDestroy(A,ierr)
  call PetscFinalize(ierr)

11 format('Timestep:', i7, '  Time:', es9.2, '  time/us:', f6.2, ' min')
12 format('   dT:  ', es9.2, '   tm:', es9.2)
13 format('  trk:  ', f5.2, '      tph:', f5.2)
14 format('   Vd:', f7.2, '       Id:', es9.2, '     diff:' es9.2)
9  format('Simulation finished in ', i0, ' hr ', i0, ' min')

contains

  subroutine read_in
    integer :: i, narg
    character(80) :: arg

    ! Check for -help
    narg = iargc()
    if (mod(narg,2) .ne. 0) then
      if (myId == 0) then
        write(*,*)
        write(*,*) 'Usage:   mpiexec -n <nproc> ./main <options>'
        write(*,*) 'Options: -nx <nx>, -ny <ny>, -px <px>, -py <py>'
        write(*,*) '         -l <l>, -w <w>, -t <t>, -unif <T/F>'
        write(*,*)
      end if
      call MPI_Finalize(ierr)
      stop
    end if

    ! Read input arguments
    do i = 1, narg/2
      call getarg(2 * (i - 1) + 1, arg)
      select case (arg)
        case ('-nx')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) nx
        case ('-ny')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) ny
        case ('-px')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) px
        case ('-py')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) py
        case ('-l')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) l
            l = l / x0
        case ('-w')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) w
          w = w / x0
        case ('-ew')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) ew
          ew = ew / x0
        case ('-t')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) t_fin
        case ('-dt')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) dt
        case ('-v')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) vl
          vl = vl / ph0
        case ('-r')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) res
        case ('-unif')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) unif
        case ('-cyl')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) cyl
        case ('-wall')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) rwall
        case ('-n0')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) n_init
          n_init = n_init * x0**3
        case ('-rf')
          call getarg(2 * (i - 1) + 2, arg)
          read(arg,*) rf
      end select
    end do

    if (px * py .ne. nproc) then
      if (myId == 0) then
        write(*,*) 'Error: px * py must equal nproc. Stop.'
        write(*,*) '       Enter ./main -help for usage.'
        call MPI_Abort(comm,2,ierr)
      end if
    end if
    if (mod(nx,px) .ne. 0) then
      if (myId == 0) then
        write(*,*) 'Error: px needs to divide nx. Stop.'
        write(*,*) '       Enter ./main -help for usage.'
        call MPI_Abort(comm,3,ierr)
      end if
    end if
    if (mod(ny,py) .ne. 0) then
      if (myId == 0) then
        write(*,*) 'Error: py needs to divide ny. Stop.'
        write(*,*) '       Enter ./main -help for usage.'
        call MPI_Abort(comm,4,ierr)
      end if
    end if

    if (.not. rf) then
      if (ny > 1) then
          write(path,41) int(res / 10**floor(log10(res))), floor(log10(res))
      else
          write(path,42) int(res / 10**floor(log10(res))), floor(log10(res))
      end if
    else
      if (ny > 1) then
          write(path,43) int(vl*ph0/10.0)*10
      else
          write(path,44) int(vl*ph0/10.0)*10
      end if
    end if

    res = res * e / (ph0 * t0)


    if (myId == 0) call system('mkdir '//trim(path))

    41 format('Output/2d_res_',i0,'e',i0,'/')
    42 format('Output/1d_res_',i0,'e',i0,'/')
    43 format('Output/2d_pulse_',i0,'V/')
    44 format('Output/1d_pulse_',i0,'V/')
  end subroutine

  subroutine writeOut
    call savedat(trim(path)//'f1.dat', ph * ph0)
    call savedat(trim(path)//'f2.dat', ne(:,:,1) / x0**3)
    call savedat(trim(path)//'f3.dat', ni(:,:,1) / x0**3)
    call savedat(trim(path)//'f4.dat', nt(:,:,1) / x0**3 * ph0 / 1.5)
    call savedat(trim(path)//'f5.dat', nm(:,:,1) / x0**3)

    call MPI_File_Open(comm, trim(path)//'time.dat', &
        MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
    if (myId == 0) call MPI_File_Write(fh, g%t, 1, etype, stat, ierr)
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, trim(path)//'vd.dat', &
        MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
    if (myId == 0) call MPI_File_Write(fh, Vd_pl*ph0, 1, etype, stat, ierr)
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, trim(path)//'id.dat', &
        MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
    if (myId == 0) call MPI_File_Write(fh, Id*e/t0, 1, etype, stat, ierr)
    call MPI_File_Close(fh, ierr)

    call MPI_File_Open(comm, trim(path)//'dt.dat', &
        MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
    if (myId == 0) call MPI_File_Write(fh, g%dt, 1, etype, stat, ierr)
    call MPI_File_Close(fh, ierr)
  end subroutine
end program
