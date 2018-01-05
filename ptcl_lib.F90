module ptcl_lib
  use ptcl_props
  implicit none

  real(8), allocatable :: ki(:,:,:), ni(:,:,:), &
                          ke(:,:,:), ne(:,:,:), &
                          km(:,:,:), nm(:,:,:), &
                          kt(:,:,:), nt(:,:,:)
  real(8) :: err_prev = 1

contains

  subroutine ptcl_init(g)
    type(grid), intent(inout) :: g
    allocate( ki(g%bx+2, g%by+2, 5), ni(g%bx+2, g%by+2, 3), &
              ke(g%bx+2, g%by+2, 5), ne(g%bx+2, g%by+2, 3), &
              kt(g%bx+2, g%by+2, 5), nt(g%bx+2, g%by+2, 3), &
              km(g%bx+2, g%by+2, 5), nm(g%bx+2, g%by+2, 3) )
    ki = 0
    ke = 0
    kt = 0
    km = 0
    ni = n_init
    ne = n_init
    nt = n_init / ph0 / 1d2
    nm = n_init
  end subroutine

  subroutine ptcl_step(g, ph, sig)
    type(grid), intent(inout) :: g
    real(8), intent(in) :: ph(:,:), sig(:)
    integer :: i, j, stage
    real(8) :: err_cur, scfac, err_ni, err_ne, err_nt, err_nm

    stage = 0
    do
      stage = stage + 1
      if (stage > 5) exit

      ! Update particle densities
      do j = 2, g%by+1
        do i = 2, g%bx+1
          call ptclEqn(g, i, j, stage, ph, sig)
        end do
      end do

      call rk_step(g, ni, ki, stage, g%dt, err_ni, 4, n_zero)
      call rk_step(g, ne, ke, stage, g%dt, err_ne, 4, n_zero)
      call rk_step(g, nt, kt, stage, g%dt, err_nt, 4, n_zero/ph0/1d2)
      call rk_step(g, nm, km, stage, g%dt, err_nm, 4, n_zero)

      if (stage == 5) then
        err_cur = (err_ni**2 + err_ne**2 + &
                   err_nm**2 + err_nt**2)**0.5
        scfac = 8d-1 * err_cur**(-7d-1 / 4d0) &
                * err_prev**( 4d-1 / 4d0)
        scfac = min(2.5d0, max(3d-1, scfac))
        g%dt = scfac * g%dt

        ! call MPI_Allreduce(MPI_In_Place, tau_m, 1, etype, &
        !                    MPI_Min, comm, ierr)
        ! if (g%ny > 1) then
        !   g%dt = min(max(g%dt, 1d-12), tau_m*sqrt(2d0))
        ! else
        !   g%dt = min(max(g%dt, 1d-12), tau_m*2d0)
        ! end if

        call MPI_Bcast(g%dt, 1, etype, 0, comm, ierr)

        if (g%dt <= 1.1d-12) then
          write(*,*) 'minimum timestep reached; finishing simulation'
          write(*,'(es10.2)') err_cur
          stop
        end if

        if (err_cur .le. 1d0) then
          err_prev = err_cur
        else
          stage = 0

          ni(:,:,1) = ni(:,:,3)
          ne(:,:,1) = ne(:,:,3)
          nm(:,:,1) = nm(:,:,3)
          nt(:,:,1) = nt(:,:,3)

          ni(:,:,2) = ni(:,:,3)
          ne(:,:,2) = ne(:,:,3)
          nm(:,:,2) = nm(:,:,3)
          nt(:,:,2) = nt(:,:,3)
        end if
      end if
    end do

    ni(:,:,3) = ni(:,:,2)
    ne(:,:,3) = ne(:,:,2)
    nm(:,:,3) = nm(:,:,2)
    nt(:,:,3) = nt(:,:,2)

    ni(:,:,2) = ni(:,:,1)
    ne(:,:,2) = ne(:,:,1)
    nm(:,:,2) = nm(:,:,1)
    nt(:,:,2) = nt(:,:,1)
  end subroutine

  subroutine ptclEqn(g, i, j, stage, ph, sig)
    type(grid), intent(in) :: g
    integer, intent(in)    :: i, j, stage
    real(8), intent(in)    :: ph(:,:), sig(:)
    real(8) :: a, Ex(2) = 0, Ey(2) = 0, Te(3), mue(2), mut(2), De(2), Dt(2), &
               dfluxi_dx = 0, dfluxi_dy = 0, fluxi_x(2) = 0, fluxi_y(2) = 0, &
               dfluxe_dx = 0, dfluxe_dy = 0, fluxe_x(2) = 0, fluxe_y(2) = 0, &
               dfluxt_dx = 0, dfluxt_dy = 0, fluxt_x(2) = 0, fluxt_y(2) = 0, &
               dfluxm_dx = 0, dfluxm_dy = 0, fluxm_x(2) = 0, fluxm_y(2) = 0, &
               term_sie, term_st1, term_st2, term_st3, term_sm, &
               k_ir, k_ex, k_sc, k_si, nu, ve

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
    ! - center -
    if (g%type_x(i-1,j-1) == 0) then
      ! rates and coefficients
      Te(1) = get_Te(nt(i-1,j,2), ne(i-1,j,2))
      Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))
      Te(3) = get_Te(nt(i+1,j,2), ne(i+1,j,2))

      mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
      mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))

      De(1) = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))
      De(2) = 5d-1 * (get_De(Te(2)) + get_De(Te(3)))

      mut(1) = 5d-1 * (get_mut(Te(1)) + get_mut(Te(2)))
      mut(2) = 5d-1 * (get_mut(Te(2)) + get_mut(Te(3)))

      Dt(1) = 5d-1 * (get_Dt(Te(1)) + get_Dt(Te(2)))
      Dt(2) = 5d-1 * (get_Dt(Te(2)) + get_Dt(Te(3)))

      ! Flux at i - 1/2
      call get_flux(fluxi_x(1), Ex(1), g%dx(i-1),  1, mui, Di, &
                    ni(i-1,j,2), ni(i,j,2))
      call get_flux(fluxe_x(1), Ex(1), g%dx(i-1), -1, mue(1), De(1), &
                    ne(i-1,j,2), ne(i,j,2))
      call get_flux(fluxt_x(1), Ex(1), g%dx(i-1), -1, mut(1), Dt(1), &
                    nt(i-1,j,2), nt(i,j,2))
      fluxm_x(1) = -Dm * (nm(i,j,2) - nm(i-1,j,2)) / g%dx(i-1)

      ! Flux at i + 1/2
      call get_flux(fluxi_x(2), Ex(2), g%dx(i),  1, mui, Di, &
                    ni(i,j,2), ni(i+1,j,2))
      call get_flux(fluxe_x(2), Ex(2), g%dx(i), -1, mue(2), De(2), &
                    ne(i,j,2), ne(i+1,j,2))
      call get_flux(fluxt_x(2), Ex(2), g%dx(i), -1, mut(2), Dt(2), &
                    nt(i,j,2), nt(i+1,j,2))
      fluxm_x(2) = -Dm * (nm(i+1,j,2) - nm(i,j,2)) / g%dx(i)

    ! - left -
    else if (g%type_x(i-1,j-1) < 0) then
      ! rates and coefficients
      Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))
      Te(3) = get_Te(nt(i+1,j,2), ne(i+1,j,2))

      mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))
      De(2) = 5d-1 * (get_De(Te(2)) + get_De(Te(3)))

      mut(2) = 5d-1 * (get_mut(Te(2)) + get_mut(Te(3)))
      Dt(2) = 5d-1 * (get_Dt(Te(2)) + get_Dt(Te(3)))

      ! Flux at i + 1/2
      call get_flux(fluxi_x(2), Ex(2), g%dx(i),  1, mui, Di, &
                    ni(i,j,2), ni(i+1,j,2))
      call get_flux(fluxe_x(2), Ex(2), g%dx(i), -1, mue(2), De(2), &
                    ne(i,j,2), ne(i+1,j,2))
      call get_flux(fluxt_x(2), Ex(2), g%dx(i), -1, mut(2), Dt(2), &
                    nt(i,j,2), nt(i+1,j,2))
      fluxm_x(2) = -Dm * (nm(i+1,j,2) - nm(i,j,2)) / g%dx(i)

      ! - electrode -
      if (g%type_x(i-1,j-1) == -2) then
        if (ph(i-1,j) > ph(i,j)) then
          a = 1d0 ! electrons drift
        else
          a = 0d0 ! ions drift
        end if

        mue(1) = get_mue(Te(2))
        mut(1) = get_mut(Te(2))
        ve = sqrt((16d0 * e * ph0 * Te(2)) / (3d0 * pi * me)) * t0 / x0

        ! Flux at i - 1/2
        fluxi_x(1) = (1d0 - a) * mui * Ex(1) * ni(i,j,2) &
                     - 2.5d-1 * vi * ni(i,j,2)

        fluxe_x(1) = - a * mue(1) * Ex(1) * ne(i,j,2) &
                     - 2.5d-1 * ve * ne(i,j,2) &
                     - gam * fluxi_x(1)

        fluxt_x(1) = - a * mut(1) * Ex(1) * nt(i,j,2) &
                     - 1d0/3d0 * ve * nt(i,j,2) &
                     - gam * Te(2) * fluxi_x(1)

        fluxm_x(1) = - 2.5d-1 * vi * nm(i,j,2)

      ! - vacuum -
      else if (g%type_x(i-1,j-1) == -1) then
        ! Flux at i - 1/2
        fluxi_x(1) = 0d0
        fluxe_x(1) = 0d0
        fluxt_x(1) = 0d0
        fluxm_x(1) = 0d0
      end if

    ! - right -
    else if (g%type_x(i-1,j-1) > 0) then
      ! rates and coefficients
      Te(1) = get_Te(nt(i-1,j,2), ne(i-1,j,2))
      Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))

      mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
      De(1) = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))

      mut(1) = 5d-1 * (get_mut(Te(1)) + get_mut(Te(2)))
      Dt(1) = 5d-1 * (get_Dt(Te(1)) + get_Dt(Te(2)))

      ! Flux at i - 1/2
      call get_flux(fluxi_x(1), Ex(1), g%dx(i-1),  1, mui, Di, &
                    ni(i-1,j,2), ni(i,j,2))
      call get_flux(fluxe_x(1), Ex(1), g%dx(i-1), -1, mue(1), De(1), &
                    ne(i-1,j,2), ne(i,j,2))
      call get_flux(fluxt_x(1), Ex(1), g%dx(i-1), -1, mut(1), Dt(1), &
                    nt(i-1,j,2), nt(i,j,2))
      fluxm_x(1) = -Dm * (nm(i,j,2) - nm(i-1,j,2)) / g%dx(i-1)

      ! - electrode -
      if (g%type_x(i-1,j-1) == 2) then
        if (ph(i+1,j) > ph(i,j)) then
            a = 1d0 ! electrons drift
        else
            a = 0d0 ! ions drift
        end if

        mue(2) = get_mue(Te(2))
        mut(2) = get_mut(Te(2))
        ve = sqrt((16d0 * e * ph0 * Te(2)) / (3d0 * pi * me)) * t0 / x0

        ! Flux at i + 1/2
        fluxi_x(2) = (1d0 - a) * mui * Ex(2) * ni(i,j,2) &
                     + 2.5d-1 * vi * ni(i,j,2)

        fluxe_x(2) = - a * mue(2) * Ex(2) * ne(i,j,2) &
                     + 2.5d-1 * ve * ne(i,j,2) &
                     - gam * fluxi_x(2)

        fluxt_x(2) = - a * mut(2) * Ex(2) * nt(i,j,2) &
                     + 1d0/3d0 * ve * nt(i,j,2) &
                     - gam * Te(2) * fluxi_x(2)

        fluxm_x(2) = 2.5d-1 * vi * nm(i,j,2)

      ! - vacuum -
      else if (g%type_x(i-1,j-1) == 1) then
        ! Flux at i + 1/2
        fluxi_x(2) = 0d0
        fluxe_x(2) = 0d0
        fluxt_x(2) = 0d0
        fluxm_x(2) = 0d0
      end if
    end if

    ! Y-dir Fluxes
    if (g%ny > 1) then
      if (g%type_y(i-1,j-1) == -1) then
        Ey(1) = 0d0
        Ey(2) = -(ph(i,j+1) - ph(i,j)) / g%dy(j)
      else if (g%type_y(i-1,j-1) == 1) then
        Ey(1) = -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)
        Ex(2) = 0d0
      else if (g%type_y(i-1,j-1) == 3) then
        Ey(1) = -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)
        Ey(2) = -sig(i)
      else
        Ey(1) = -(ph(i,j) - ph(i,j-1)) / g%dy(j-1)
        Ey(2) = -(ph(i,j+1) - ph(i,j)) / g%dy(j)
      end if

      ! - center -
      if (g%type_y(i-1,j-1) == 0) then
        ! rates and coefficients
        Te(1) = get_Te(nt(i,j-1,2), ne(i,j-1,2))
        Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))
        Te(3) = get_Te(nt(i,j+1,2), ne(i,j+1,2))

        mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))

        De(1) = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))
        De(2) = 5d-1 * (get_De(Te(2)) + get_De(Te(3)))

        mut(1) = 5d-1 * (get_mut(Te(1)) + get_mut(Te(2)))
        mut(2) = 5d-1 * (get_mut(Te(2)) + get_mut(Te(3)))

        Dt(1) = 5d-1 * (get_Dt(Te(1)) + get_Dt(Te(2)))
        Dt(2) = 5d-1 * (get_Dt(Te(2)) + get_Dt(Te(3)))

        ! Flux at j - 1/2
        call get_flux(fluxi_y(1), Ey(1), g%dy(j-1),  1, mui, Di, &
                      ni(i,j-1,2), ni(i,j,2))
        call get_flux(fluxe_y(1), Ey(1), g%dy(j-1), -1, mue(1), De(1), &
                      ne(i,j-1,2), ne(i,j,2))
        call get_flux(fluxt_y(1), Ey(1), g%dy(j-1), -1, mut(1), Dt(1), &
                      nt(i,j-1,2), nt(i,j,2))
        fluxm_y(1) = -Dm * (nm(i,j,2) - nm(i,j-1,2)) / g%dy(j-1)


        ! Flux at j + 1/2
        call get_flux(fluxi_y(2), Ey(2), g%dy(j),  1, mui, Di, &
                      ni(i,j,2), ni(i,j+1,2))
        call get_flux(fluxe_y(2), Ey(2), g%dy(j), -1, mue(2), De(2), &
                      ne(i,j,2), ne(i,j+1,2))
        call get_flux(fluxt_y(2), Ey(2), g%dy(j), -1, mut(2), Dt(2), &
                      nt(i,j,2), nt(i,j+1,2))
        fluxm_y(2) = -Dm * (nm(i,j+1,2) - nm(i,j,2)) / g%dy(j)

      ! - left -
      else if (g%type_y(i-1,j-1) < 0) then
        ! rates and coefficients
        Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))
        Te(3) = get_Te(nt(i,j+1,2), ne(i,j+1,2))

        mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))
        De(2) = 5d-1 * (get_De(Te(2)) + get_De(Te(3)))

        mut(2) = 5d-1 * (get_mut(Te(2)) + get_mut(Te(3)))
        Dt(2) = 5d-1 * (get_Dt(Te(2)) + get_Dt(Te(3)))

        ! Flux at j + 1/2
        call get_flux(fluxi_y(2), Ey(2), g%dy(j),  1, mui, Di, &
                      ni(i,j,2), ni(i,j+1,2))
        call get_flux(fluxe_y(2), Ey(2), g%dy(j), -1, mue(2), De(2), &
                      ne(i,j,2), ne(i,j+1,2))
        call get_flux(fluxt_y(2), Ey(2), g%dy(j), -1, mut(2), Dt(2), &
                      nt(i,j,2), nt(i,j+1,2))
        fluxm_y(2) = -Dm * (nm(i,j+1,2) - nm(i,j,2)) / g%dy(j)

        ! Flux at j - 1/2
        fluxi_y(1) = 0d0
        fluxe_y(1) = 0d0
        fluxt_y(1) = 0d0
        fluxm_y(1) = 0d0

      ! - right -
      else if (g%type_y(i-1,j-1) > 0) then
        ! rates and coefficients
        Te(1) = get_Te(nt(i,j-1,2), ne(i,j-1,2))
        Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))

        mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        De(1) = 5d-1 * (get_De(Te(1)) + get_De(Te(2)))

        mut(1) = 5d-1 * (get_mut(Te(1)) + get_mut(Te(2)))
        Dt(1) = 5d-1 * (get_Dt(Te(1)) + get_Dt(Te(2)))

        ! Flux at j - 1/2
        call get_flux(fluxi_y(1), Ey(1), g%dy(j-1),  1, mui, Di, &
                      ni(i,j-1,2), ni(i,j,2))
        call get_flux(fluxe_y(1), Ey(1), g%dy(j-1), -1, mue(1), De(1), &
                      ne(i,j-1,2), ne(i,j,2))
        call get_flux(fluxt_y(1), Ey(1), g%dy(j-1), -1, mut(1), Dt(1), &
                      nt(i,j-1,2), nt(i,j,2))
        fluxm_y(1) = -Dm * (nm(i,j,2) - nm(i,j-1,2)) / g%dy(j-1)

        ! Flux at j + 1/2
        if (g%type_y(i-1,j-1) == 3) then
          if (Ey(2) < 0d0) then
            a = 1d0 ! electrons drift
          else
            a = 0d0 ! ions drift
          end if

          mue(2) = get_mue(Te(2))
          mut(2) = get_mut(Te(2))
          ve = sqrt((16d0 * e * ph0 * Te(2)) / (3d0 * pi * me)) * t0 / x0

          ! Flux at j + 1/2
          fluxi_y(2) = - (1d0 - a) * mui * Ey(2) * ni(i,j,2) &
                       + 2.5d-1 * vi * ni(i,j,2)

          fluxe_y(2) = - a * mue(2) * Ey(2) * ne(i,j,2) &
                       + 2.5d-1 * ve * ne(i,j,2)

          fluxt_y(2) = - a * mut(2) * Ey(2) * nt(i,j,2) &
                       + 1d0/3d0 * ve * nt(i,j,2)

          fluxm_y(2) = 2.5d-1 * vi * nm(i,j,2)
        else
          fluxi_y(2) = 0d0
          fluxe_y(2) = 0d0
          fluxt_y(2) = 0d0
          fluxm_y(2) = 0d0
        end if
      end if
    end if

    ! X-dir flux gradients
    dfluxi_dx = (fluxi_x(2) - fluxi_x(1)) / g%dlx(i-1)
    dfluxe_dx = (fluxe_x(2) - fluxe_x(1)) / g%dlx(i-1)
    dfluxt_dx = (fluxt_x(2) - fluxt_x(1)) / g%dlx(i-1)
    dfluxm_dx = (fluxm_x(2) - fluxm_x(1)) / g%dlx(i-1)

    ! X-dir midpoint fluxes
    fluxe_x(1) = (fluxe_x(2) + fluxe_x(1)) / 2d0
    Ex(1) = (Ex(1) + Ex(2)) / 2d0

    if (g%ny > 1) then
      ! Y-dir flux gradients
      if (cyl) then
        dfluxi_dy = (fluxi_y(2) * (g%r(j+1) + g%r(j)) / 2d0   &
                     - fluxi_y(1) * (g%r(j) + g%r(j-1)) / 2d0) &
                     / g%dly(j-1) / g%r(j)
        dfluxe_dy = (fluxe_y(2) * (g%r(j+1) + g%r(j)) / 2d0   &
                     - fluxe_y(1) * (g%r(j) + g%r(j-1)) / 2d0) &
                     / g%dly(j-1) / g%r(j)
        dfluxt_dy = (fluxt_y(2) * (g%r(j+1) + g%r(j)) / 2d0   &
                     - fluxt_y(1) * (g%r(j) + g%r(j-1)) / 2d0) &
                     / g%dly(j-1) / g%r(j)
        dfluxm_dy = (fluxm_y(2) * (g%r(j+1) + g%r(j)) / 2d0   &
                     - fluxm_y(1) * (g%r(j) + g%r(j-1)) / 2d0) &
                     / g%dly(j-1) / g%r(j)
      else
        dfluxi_dy = (fluxi_y(2) - fluxi_y(1)) / g%dly(j-1)
        dfluxe_dy = (fluxe_y(2) - fluxe_y(1)) / g%dly(j-1)
        dfluxt_dy = (fluxt_y(2) - fluxt_y(1)) / g%dly(j-1)
        dfluxm_dy = (fluxm_y(2) - fluxm_y(1)) / g%dly(j-1)
      end if

      ! Y-dir midpoint fluxes
      fluxe_y(1) = (fluxe_y(2) + fluxe_y(1))/ 2d0
      Ey(1) = (Ey(1) + Ey(2)) / 2d0
    end if

    ! rates and coefficients
    k_ir = get_k_ir(Te(2))
    k_sc = get_k_sc(Te(2))
    k_si = get_k_si(Te(2))
    k_ex = get_k_ex(Te(2))
    nu   = get_nu(Te(2))

    ! evaluate source terms
    term_sie =   k_ir * ninf * ne(i,j,2) &
               - beta * ni(i,j,2) * ne(i,j,2) &
               + k_si * nm(i,j,2) * ne(i,j,2) &
               + k_mp * nm(i,j,2)**2

    term_sm =   k_ex * ninf * ne(i,j,2) &
              - k_si * nm(i,j,2) * ne(i,j,2) &
              - k_sc * nm(i,j,2) * ne(i,j,2) &
              - k_r  * nm(i,j,2) * ne(i,j,2) &
              - 2d0  * k_mp * nm(i,j,2)**2 &
              - k_2q * ninf * nm(i,j,2) &
              - k_3q * ninf**2 * nm(i,j,2)

    ! -e flux_e . E
    term_st1 = - (fluxe_x(1) * Ex(1) + fluxe_y(1) * Ey(1))

    ! -me/mg nu_e (Te - Tg)
    term_st2 = - nt(i,j,2) * nu * me / mi

    ! reactions
    term_st3 = - h_ir * k_ir * ninf * ne(i,j,2) &
               - h_si * k_si * nm(i,j,2) * ne(i,j,2) &
               - h_ex * k_ex * ninf * ne(i,j,2) &
               - h_sc * k_sc * nm(i,j,2) * ne(i,j,2)

    ! evaluate expression
    ki(i,j,stage) = -dfluxi_dx - dfluxi_dy + term_sie
    ke(i,j,stage) = -dfluxe_dx - dfluxe_dy + term_sie
    kt(i,j,stage) = -dfluxt_dx - dfluxt_dy + term_st1 + term_st2 + term_st3
    km(i,j,stage) = -dfluxm_dx - dfluxm_dy + term_sm

    if (stage == 5) then
        t_m = min(t_m, 1d0 / (mue(2) * ne(i,j,2) + mui * ni(i,j,2)))
    end if

    if (isnan(ki(i,j,stage))) then
        write(*,*) 'i,j,stage ', i,j,stage
        write(*,*) 'ki(i,j,stage)',ki(i,j,stage)
        call MPI_Finalize(ierr)
    end if
    if (isnan(ke(i,j,stage))) then
        write(*,*) 'i,j,stage ', i,j,stage
        write(*,*) 'ke(i,j,stage)',ke(i,j,stage)
        call MPI_Finalize(ierr)
    end if
    if (isnan(kt(i,j,stage))) then
        write(*,*) 'i,j,stage ', i,j,stage
        write(*,*) 'kt(i,j,stage)',kt(i,j,stage)
        call MPI_Finalize(ierr)
    end if
    if (isnan(km(i,j,stage))) then
        write(*,*) 'i,j,stage ', i,j,stage
        write(*,*) 'km(i,j,stage)',km(i,j,stage)
        call MPI_Finalize(ierr)
    end if
  end subroutine
end module
