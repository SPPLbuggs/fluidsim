module ptcl_props
    use props
    implicit none

    ! electron properties
    real(8), parameter:: me = 9.10938188d-31

    ! positive ion properties
    real(8), parameter:: mi  = 1.67262178d-27 * 39.948d0, &
                         mui = 1.45d3 / p * 1d-4 * ph0 * t0 / x0**2, &
                         Ti  = Tg, &
                         vi  = Sqrt( (8d0 * kb*Ti) / (pi * mi) ) * t0 / x0, &
                         Di  = mui * (kb*Ti / e) / ph0

    ! metastable argon properties
    real(8), parameter:: mm   = mi, &
                         Dm   = 2.42d18 * x0 / ninf / 1d2 * t0, &
                         k_r  = 2d-7 / 1d6 / x0**3 * t0, &
                         k_mp = 6.2d-10 / 1.0d6 / x0**3 * t0, &
                         k_2q = 3d-15 / 1d6 / x0**3 * t0, &
                         k_3q = 1.1d-31 / 1.0d12 / x0**6 * t0

    ! reactions
    real(8), parameter:: H_ir  =  15.80d0 / ph0, &
                         H_ex  =  11.56d0 / ph0, &
                         H_si  =   4.14d0 / ph0, &
                         H_sc  = -11.56d0 / ph0, &
                         beta  =  1d-13 * t0 / x0**3, &
                         gam   =   0.10d0

contains
! *** Runge-Kutta Adaptive Timestepping ***
  subroutine rk_step(g, n, k, stage, dt, nerr_n, order, n_min)
    type(grid), intent(in) :: g
    integer, intent(in) :: stage, order
    real(8), intent(in) :: dt, k(:,:,:), n_min
    real(8), intent(inout) :: n(:,:,:), nerr_n
    real(8) :: err_n(g%bx+2, g%by+2), abs_tol = 1d-2, rel_tol = 1d-2
    integer :: i,j

    if (order == 1) then
    ! euler scheme
      do j = 2, g%by+1
        do i = 2, g%bx+1
          n(i,j,1) = n(i,j,3) + k(i,j,1) * dt
        end do
      end do

      call comm_real(g%bx, g%by, n(:,:,1))

    else if (order == 2) then
    ! 2nd order
      if (stage == 1) then
        do j = 2, g%by+1
          do i = 2, g%bx+1
            n(i,j,1) = n(i,j,3) + k(i,j,1) * dt / 2d0
            n(i,j,2) = n(i,j,3) + k(i,j,1) * dt
          end do
        end do

      else
        do j = 2, g%by+1
          do i = 2, g%bx+1
            n(i,j,1) = n(i,j,1) + k(i,j,2) * dt / 2d0
          end do
        end do
      end if

      n(:,:,1:2) = max(n(:,:,1:2), n_min)

      call comm_real(g%bx, g%by, n(:,:,1))
      call comm_real(g%bx, g%by, n(:,:,2))

    else if (order == 4) then
      ! merson 4("5") adaptive time-stepping
      if (stage == 1) then
        do j = 2, g%by+1
          do i = 2, g%bx+1
            n(i,j,2) = n(i,j,3) + k(i,j,1) * dt / 3d0
          end do
        end do

        n(:,:,2) = max(n(:,:,2), n_min)
        call comm_real(g%bx, g%by, n(:,:,2))

      else if (stage == 2) then
        do j = 2, g%by+1
          do i = 2, g%bx+1
            n(i,j,2) = n(i,j,3) + dt * ( k(i,j,1) + k(i,j,2) ) / 6d0
           end do
        end do

        n(:,:,2) = max(n(:,:,2), n_min)

        call comm_real(g%bx, g%by, n(:,:,2))

      else if (stage == 3) then
        do j = 2, g%by+1
          do i = 2, g%bx+1
            n(i,j,2) = n(i,j,3) + dt * (k(i,j,1) + k(i,j,3) * 3d0) / 8d0
          end do
        end do

        n(:,:,2) = max(n(:,:,2), n_min)
        call comm_real(g%bx, g%by, n(:,:,2))

    else if (stage == 4) then
        do j = 2, g%by+1
          do i = 2, g%bx+1
            n(i,j,2) = n(i,j,3) + dt * (k(i,j,1) &
                       - k(i,j,3) * 3d0 + k(i,j,4) * 4d0 ) / 2d0
          end do
        end do

        n(:,:,2) = max(n(:,:,2), n_min)
        call comm_real(g%bx, g%by, n(:,:,2))
    else
        do j = 2, g%by+1
          do i = 2, g%bx+1
            n(i,j,1) = n(i,j,3) + dt * (k(i,j,1) &
                       + k(i,j,4) * 4d0 + k(i,j,5)) / 6d0
          end do
        end do

        n(:,:,1) = max(n(:,:,1), n_min)
        call comm_real(g%bx, g%by, n(:,:,1))

        err_n = 0

        err_n = abs(dt * (k(:,:,1) * 2d0 / 30d0 - k(:,:,3) * 3d0 / 10d0 &
                    + k(:,:,4) * 4d0 / 15d0 - k(:,:,5) / 30d0 ))

        ! nerr_n = maxval(err_n/(abs_tol+rel_tol*abs(n(:,:,3))))
        ! call MPI_Allreduce(MPI_In_Place, nerr_n, 1, etype, &
        !                    MPI_Max, comm, ierr)

        nerr_n = sum((err_n/(abs_tol+rel_tol*abs(n(:,:,3))))**2)
        call MPI_Allreduce(MPI_In_Place, nerr_n, 1, etype, &
                           MPI_Sum, comm, ierr)
        ! nerr_n = sqrt(nerr_n)
      end if
    end if
  end subroutine

  subroutine get_Flux(flux, E, dh, q, mu, D, n_left, n_right)
    real(8), intent(inout) :: flux
    integer, intent(in) :: q
    real(8), intent(in) :: E, dh, mu, D, n_left, n_right
    real(8) :: v, tol, arg

    tol = 1e-12
    v = q * mu * E
    arg = v * dh / D

    if (abs(q*E) < tol) then
      flux = D * (n_left - n_right) / dh

    ! Positive exponentials blow up,
    !  so rewrite always as negative exp.
    !  this is the same analytical expression
    else if (arg > 0) then
      flux = v * (n_left - n_right * exp(max(-arg,-100d0))) &
           / (1d0 - exp(max(-arg, -100d0)))
    else
      flux = v * (n_right - n_left * exp(max(arg, -100d0))) &
           / (1d0  - exp( max(arg,-100d0)))
    end if
  end subroutine

  function get_Te(nte,ne)
      real(8):: get_Te
      real(8), intent(in) :: nte, ne

      if ((nte >= 0d0) .and. (ne >= 0d0)) then
          get_Te = nte / ne
      else
          write(*,*) "Error, negative density. Stop."
          stop
          get_Te = 1d-8
      end if

      return
  end function

  function get_mue(T)
    real(8):: get_mue
    real(8), intent(in):: T
    real(8):: x, &
              a1 =  58.1133016145,    &
              b1 =  -0.984082217962,  &
              c1 =  -0.164770900119,  &
              d1 =  -0.142058042584,  &
              f1 =  -0.0637079234081, &
              g1 =   0.0436642742558, &
              a2 =  70.7846754548,    &
              b2 = -22.7558138237,    &
              c2 =  12.8366493242,    &
              d2 =  -3.57242763244,   &
              f2 =   0.486276623664,  &
              g2 =  -0.0259045697422

    x = log(max(2.34d-1, min(1.57d2, T * ph0)))

    if (x < log(5.119)) then
        get_mue = exp(a1 + b1*x + c1*x**2 + d1*x**3 + f1*x**4 &
                      + g1*x**5) * x0 / ninf * t0 * ph0
    else
        get_mue = exp(a2 + b2*x + c2*x**2 + d2*x**3 + f2*x**4 &
                      + g2*x**5) * x0 / ninf * t0 * ph0
    end if
    return
  end function get_mue

  function get_mut(T)
    real(8):: get_mut
    real(8), intent(in):: T
    real(8):: x, &
              a1 =  58.1354622038,    &
              b1 =  -1.29918391109,   &
              c1 =  -0.106559658509,  &
              d1 =  -0.0202265418394, &
              f1 =  -0.0286667333915, &
              g1 =   0.022457210094,  &
              a2 =  68.2627663433,    &
              b2 = -19.3524000831,    &
              c2 =  11.3251980583,    &
              d2 =  -3.24070386509,   &
              f2 =   0.449848156054,  &
              g2 =  -0.0242585762405

    x = log(max(2.34d-1, min(1.57d2, T * ph0)))

    if (x < log(5.119)) then
        get_mut = exp(a1 + b1*x + c1*x**2 + d1*x**3 + f1*x**4 &
                      + g1*x**5) * x0 / ninf * t0 * ph0
    else
        get_mut = exp(a2 + b2*x + c2*x**2 + d2*x**3 + f2*x**4 &
                      + g2*x**5) * x0 / ninf * t0 * ph0
    end if
    return
  end function get_mut

  function get_mui(E)
    real(8) :: get_mui
    real(8), intent(in) :: E
    real(8) :: x, &
               a =  0.476914771378,   &
               b = -0.227345106588,   &
               c =  0.0530705372631,  &
               d = -0.00521915548354, &
               f =  0.000176442914616

    x = log(min(1d5, max(E * ph0 / x0, 1d2)))

    get_mui = (a + b*x + c*x**2 + d*x**3 + f*x**4) * ph0 * t0 / x0**2
    return
  end function

  function get_De(T)
    real(8):: get_De
    real(8), intent(in):: T
    real(8):: x, &
              a = 57.7067018813,      &
              b = -0.0699892875381,   &
              c = -0.117645585949,    &
              d =  0.105390295278,    &
              f = -0.102862612604,    &
              g = -0.0469171521686,   &
              h =  0.0584908312121,   &
              i =  0.000578601715687, &
              j = -0.0122860884883,   &
              k =  0.00426793748856,  &
              l = -0.000590000082557, &
              m = 3.00706533201e-05

    x = log(max(2.34d-1, min(1.57d2, T * ph0)))

    get_De = exp(a + b*x + c*x**2 + d*x**3 + f*x**4 + g*x**5 &
                 + h*x**6 + i*x**7 + j*x**8 + k*x**9 &
                 + l*x**10 + m*x**11) * x0 / ninf * t0
    return
  end function get_De

  function get_Dt(T)
    real(8) :: get_Dt
    real(8), intent(in):: T
    real(8):: x, &
              a = 57.7214952841,     &
              b = -0.353731464348,   &
              c = -0.00505156731795, &
              d =  0.161173720584,   &
              f = -0.196610345869,   &
              g = -0.0719115643218,  &
              h =  0.11874085868,    &
              i = -0.00669712724784, &
              j = -0.0236445308504 , &
              k =  0.00917326671764, &
              l = -0.00135453096483, &
              m =  7.26379684461e-05

    x = log(max(2.34d-1, min(1.57d2, T * ph0)))

    get_Dt = exp(a + b*x + c*x**2 + d*x**3 + f*x**4 + g*x**5 &
                 + h*x**6 + i*x**7 + j*x**8 + k*x**9 &
                 + l*x**10 + m*x**11) * x0 / ninf * t0
    return
  end function get_Dt

  function get_k_ex2(T)
    real(8):: get_k_ex2
    real(8), intent(in):: T
    real(8):: x, &
              a1 =  -54.4513453969,  &
              b1 =   17.1472529739,  &
              c1 =   -1.05188982824, &
              d1 =  -10.5053010711,  &
              f1 =    7.51502551486, &
              g1 =   -1.44763942525, &
              a2 = -150.266369415,   &
              b2 =  158.67012941,    &
              c2 =  -85.5163470769,  &
              d2 =   23.0655081891,  &
              f2 =   -3.0908688844,  &
              g2 =    0.16403226902

    x = log(max(9.611d-1, min(1.57d2, T * ph0)))

    if (x < log(5.667)) then
        get_k_ex2 = exp(a1 + b1*x + c1*x**2 + d1*x**3 + f1*x**4 &
                      + g1*x**5) * t0 / x0**3
    else
        get_k_ex2 = exp(a2 + b2*x + c2*x**2 + d2*x**3 + f2*x**4 &
                      + g2*x**5) * t0 / x0**3
    end if
    return
  end function get_k_ex2

  function get_k_ir(T)
    real(8):: get_k_ir
    real(8), intent(in):: T
    real(8):: x, &
              a = -109.77,  &
              b =   81.803, &
              c =  -32.3,   &
              d =    5.729, &
              f =   -0.3803

    x = log(max(1d0, min(1.1d2, T * ph0)))

    get_k_ir = exp(a + b*x + c*x**2 + d*x**3 + f*x**4) * t0 / x0**3
    return
  end function get_k_ir

  function get_k_ir2(T)
    real(8):: get_k_ir2
    real(8), intent(in):: T
    real(8):: x, &
              a1 = -73.27080258,   &
              b1 = -30.1312514721, &
              c1 = 138.743842326,  &
              d1 = -161.378167125, &
              f1 = 81.9325486193,  &
              g1 = -14.9080865672, &
              a2 = -168.925304837, &
              b2 = 177.113452524,  &
              c2 = -92.3047166553, &
              d2 = 24.1814225482,  &
              f2 = -3.15463258104, &
              g2 = 0.16331188192

    x = log(max(1.37d0, min(1.57d2, T * ph0)))

    if (x < log(7.5)) then
        get_k_ir2 = exp(a1 + b1*x + c1*x**2 + d1*x**3 + f1*x**4 &
                      + g1*x**5) * t0 / x0**3
    else
        get_k_ir2 = exp(a2 + b2*x + c2*x**2 + d2*x**3 + f2*x**4 &
                      + g2*x**5) * t0 / x0**3
    end if
    return
  end function get_k_ir2

  function get_nu(T)
    real(8):: get_nu
    real(8), intent(in):: T
    real(8):: x, &
              a = -32.275912575,      &
              b =   1.45173283977,    &
              c =   0.00936933121094, &
              d =   0.129397015353,   &
              f =  -0.0414865809044,  &
              g =  -0.0582934303409,  &
              h =   0.0309832277826,  &
              i =  -0.00542014733763,  &
              j =   0.000325615321708

    x = log(max(2.34d-1, min(1.57d2, T * ph0)))

    get_nu = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. + g*x**5. &
                 + h*x**6. + i*x**7 + j*x**8.) / x0**3 * ninf * t0
    return
  end function get_nu

  function get_k_sc(T)
    real(8):: get_k_sc
    real(8), intent(in):: T
    real(8):: x, &
              a = -21.4827864151,   &
              b =   0.457356923276, &
              c =  -0.555439231606, &
              d =   1.27257798891,  &
              f =  -0.67840685073,  &
              g =   0.10591014464

    x = log(min(16d0, max(5d-1, T * ph0)))

    get_k_sc = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. &
                   + g*x**5. ) / 1.0d6 / x0**3 * t0
    return
  end function get_k_sc

  function get_k_si(T)
    real(8):: get_k_si
    real(8), intent(in):: T
    real(8):: x, &
              a = -43.1347385848, &
              b =  43.9905424566, &
              c = -28.1169537586, &
              d =   8.28853856817, &
              f =  -0.931626144207

    x = log(min(16d0, max(5d-1, T * ph0)))

    get_k_si = exp(a + b*x + c*x**2. + d*x**3. + f*x**4.) / 1.0d6 / x0**3 * t0
    return
  end function get_k_si

  ! old funcs
  function get_k_ex(T)
    real(8):: get_k_ex
    real(8), intent(in):: T
    real(8):: x,a,b,c,d,e,f,g,h,k
    a = -50.3785102239
    b =  19.0129183764
    c = -11.7950315424
    d =   7.41674013553
    e =  -3.84148086698
    f =   1.2962229976
    g =  -0.259359346989
    h =   0.0279182131315
    k =  -0.00124438710099
    x = log(T * ph0)
    if (x > log(2d2))  x = log(2d2)
    if (x < -5d-1) then
        get_k_ex = 0
    else
        get_k_ex = exp(a + b*x + c*x**2. + d*x**3. + e*x**4. &
                 + f*x**5. + g*x**6. + h*x**7. + k*x**8.) * t0 / x0**3
    end if
    return
  end function get_k_ex
end module
