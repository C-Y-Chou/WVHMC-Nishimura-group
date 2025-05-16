module solve_flow
   use param_mod
   use utils
   use model
contains
   subroutine dop54_step(f, y, h, res, er1, k1)
      implicit none
      real(dp), intent(in)::y(:), h
      real(dp), intent(out)::res(:), er1
      real(dp), intent(inout)::k1(:)
      real(dp), allocatable::k(:, :), y2(:)
      real(dp)::a(6, 6), c(6), b(7), d(6), sk
      integer::i, j
      interface
         function f(y) result(res)
            use, intrinsic::iso_fortran_env, only: real64
            integer, parameter::dp = real64
            real(dp), intent(in)::y(:)
            real(dp)::res(size(y))
         end function f
      end interface
      res = y
      a = reshape((/1.0_dp/5, 3.0_dp/40.0_dp, 44.0_dp/45, 19372.0_dp/6561, 9017.0_dp/3168, 35.0_dp/384, &
                    0.0_dp, 9.0_dp/40.0_dp, -56.0_dp/15, -25360.0_dp/2187, -355.0_dp/33, 0.0_dp, &
                    0.0_dp, 0.0_dp, 32.0_dp/9, 64448.0_dp/6561, 46732.0_dp/5247, 500.0_dp/1113, &
                    0.0_dp, 0.0_dp, 0.0_dp, -212.0_dp/729, 49.0_dp/176, 125.0_dp/192, &
                    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -5103.0_dp/18656, -2187.0_dp/6784, &
                    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 11.0_dp/84 &
                    /), shape(a))
      c = (/1.0_dp/5, 3.0_dp/10.0_dp, 4.0_dp/5, 8.0_dp/9, 1.0_dp, 1.0_dp/)
      b = (/5179.0_dp/57600.0_dp, 0.0_dp, 7571.0_dp/16695, 393.0_dp/640.0_dp, &
            -92097.0_dp/339200.0_dp, 187.0_dp/2100.0_dp, 1.0_dp/40/)
      d = (/-0.08536_dp, 0.088_dp, -0.0096_dp, 0.0052_dp, 0.00576_dp, -0.004_dp/)
      allocate (k(7, size(y)), y2(size(y)))
      if (abs(k1(1)) < 1d-14) then
         k(1, :) = f(res)
      else
         k(1, :) = k1
      end if
      do i = 1, 6
         k(i + 1, :) = 0
         do j = 1, i
            k(i + 1, :) = k(i + 1, :) + a(i, j)*k(j, :)
         end do
         k(i + 1, :) = f(res + h*k(i + 1, :))
      end do
      k1 = k(7, :)
      y2 = res
      do i = 1, 7
         y2 = y2 + b(i)*k(i, :)*h
      end do
      do i = 1, 6
         res = res + a(6, i)*k(i, :)*h
      end do
      er1 = 0
      do i = 1, size(y)
         sk = 1.0e-12_dp + 1.0e-10_dp*max(abs(res(i)), abs(y2(i)))
         er1 = er1 + ((res(i) - y2(i))/sk)**2/size(y)
      end do
      er1 = sqrt(er1)
      deallocate (k, y2)
   end subroutine dop54_step

   subroutine dop853_step(f, y, h, res, err, k1)
      implicit none

      !> Procedure interface for RHS
      interface
         function f(y) result(fy)
            import :: dp
            real(dp), intent(in) :: y(:)
            real(dp)             :: fy(size(y))
         end function f
      end interface

      !> Inputs
      real(dp), intent(in)          :: y(:), h
      real(dp), intent(inout)       :: k1(:)   !! Possibly precomputed f(y)
      !> Outputs
      real(dp), intent(out)         :: res(:) !! 8th-order solution
      real(dp), intent(out)         :: err     !! local error estimate

      !> Locals
      real(dp), allocatable         :: k(:, :)
      real(dp), allocatable         :: y2(:)
      real(dp)                      :: sk, err2
      integer                       :: n, i, j

      ! ----------------------------------------------------------------------
      !! DOP853 coefficients
      !!  c_i: the times (fractions of h) at which stages are taken
      !!  a_ij: the Butcher tableau for internal stages
      !!  b8_i: final combination for 8th-order solution
      !!  b5_i: combination for the 5th-order embedded solution (or use er_i)
      !!  er_i: direct difference coefficients (another common approach)
      !!
      !! We store them in arrays for clarity. Indices:
      !!   a(i, j) means the coefficient for stage i, summation of k_j
      !!   c(i)    time fraction for stage i
      !!   b8(i)   weight for stage i in final 8th-order solution
      !!   b5(i)   weight for stage i in the embedded 5th-order solution
      !!

      !> For 12 stages, define c(1) = 0.0 by convention, c(2) = 0.0526..., etc.
      !> We only need c(2..12). c(1) is conceptually 0.0 (the start).
      real(dp), parameter :: c(12) = &
                             (/0.0_dp, &
                               0.0526001519587677318785587544488_dp, &
                               0.0789002279381515978178381316732_dp, &
                               0.118350341907227396726757197510_dp, &
                               0.281649658092772603273242802490_dp, &
                               0.333333333333333333333333333333_dp, &
                               0.25_dp, &
                               0.307692307692307692307692307692_dp, &
                               0.651282051282051282051282051282_dp, &
                               0.6_dp, &
                               0.857142857142857142857142857142_dp, &
                               1.0_dp/)  !! some forms use c(11)=0.857..., c(12)=1.0

      !> a(i, j) for i=2..12, j=1..(i-1). Zeros otherwise.
      real(dp), parameter :: a(11, 11) = reshape((/ &
                                                 5.26001519587677318785587544488D-2, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, & !
                                                 1.97250569845378994544595329183D-2, &
                                                 5.91751709536136983633785987549D-2, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &!
                                                 2.95875854768068491816892993775D-2, &
                                                 0.0_dp, &
                                                 8.87627564304205475450678981324D-2, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &!
                                                 2.41365134159266685502369798665D-1, &
                                                 0.0_dp, &
                                                 -8.84549479328286085344864962717D-1, &
                                                 9.24834003261792003115737966543D-1, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &!
                                                 3.7037037037037037037037037037D-2, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 1.70828608729473871279604482173D-1, &
                                                 1.25467687566822425016691814123D-1, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &!
                                                 3.7109375D-2, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 1.70252211019544039314978060272D-1, &
                                                 6.02165389804559606850219397283D-2, &
                                                 -1.7578125D-2, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &!
                                                 3.70920001185047927108779319836D-2, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 1.70383925712239993810214054705D-1, &
                                                 1.07262030446373284651809199168D-1, &
                                                 -1.53194377486244017527936158236D-2, &
                                                 8.27378916381402288758473766002D-3, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &!
                                                 6.24110958716075717114429577812D-1, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 -3.36089262944694129406857109825D0, &
                                                 -8.68219346841726006818189891453D-1, &
                                                 2.75920996994467083049415600797D1, &
                                                 2.01540675504778934086186788979D1, &
                                                 -4.34898841810699588477366255144D1, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 0.0_dp, &!
                                                 4.77662536438264365890433908527D-1, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 -2.48811461997166764192642586468D0, &
                                                 -5.90290826836842996371446475743D-1, &
                                                 2.12300514481811942347288949897D1, &
                                                 1.52792336328824235832596922938D1, &
                                                 -3.32882109689848629194453265587D1, &
                                                 -2.03312017085086261358222928593D-2, &
                                                 0.0_dp, &
                                                 0.0_dp, &!
                                                 -9.3714243008598732571704021658D-1, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 5.18637242884406370830023853209D0, &
                                                 1.09143734899672957818500254654D0, &
                                                 -8.14978701074692612513997267357D0, &
                                                 -1.85200656599969598641566180701D1, &
                                                 2.27394870993505042818970056734D1, &
                                                 2.49360555267965238987089396762D0, &
                                                 -3.0467644718982195003823669022D0, &
                                                 0.0_dp, &
                                                 2.27331014751653820792359768449D0, &
                                                 0.0_dp, &
                                                 0.0_dp, &
                                                 -1.05344954667372501984066689879D1, &
                                                 -2.00087205822486249909675718444D0, &
                                                 -1.79589318631187989172765950534D1, &
                                                 2.79488845294199600508499808837D1, &
                                                 -2.85899827713502369474065508674D0, &
                                                 -8.87285693353062954433549289258D0, &
                                                 1.23605671757943030647266201528D1, &
                                                 6.43392746015763530355970484046D-1/), shape=(/11, 11/))

      !> b8: combination for the 8th-order solution
      real(dp), parameter :: b8(12) = (/5.42937341165687622380535766363D-2, &
                                        0.0_dp, &
                                        0.0_dp, &
                                        0.0_dp, &
                                        0.0_dp, &
                                        4.45031289275240888144113950566D0, &
                                        1.89151789931450038304281599044D0, &
                                        -5.8012039600105847814672114227D0, &
                                        3.1116436695781989440891606237D-1, &
                                        -1.52160949662516078556178806805D-1, &
                                        2.01365400804030348374776537501D-1, &
                                        4.47106157277725905176885569043D-2/)
      real(dp), parameter :: bhh(12) = (/ &
                             0.244094488188976377952755905512_dp, &
                             0.0_dp, &
                             0.0_dp, &
                             0.0_dp, &
                             0.0_dp, &
                             0.0_dp, &
                             0.0_dp, &
                             0.0_dp, &
                             0.733846688281611857341361741547_dp, &
                             0.0_dp, &
                             0.0_dp, &
                             0.0220588235294117647058823529412_dp/)

      !> The embedded 5th-order version can be formed by:
      !!   b5(i) = b8(i) + er_i   or from the original DOP853 paper
      !!   we can define them directly:
      real(dp), parameter :: er(12) = (/0.1312004499419488073250102996D-01, &
                                        0.0_dp, &
                                        0.0_dp, &
                                        0.0_dp, &
                                        0.0_dp, &
                                        -0.1225156446376204440720569753D+01, &
                                        -0.4957589496572501915214079952D+00, &
                                        0.1664377182454986536961530415D+01, &
                                        -0.3503288487499736816886487290D+00, &
                                        0.3341791187130174790297318841D+00, &
                                        0.8192320648511571246570742613D-01, &
                                        -0.2235530786388629525884427845D-01/)

      n = size(y)
      res = y
      allocate (k(12, n), y2(n))
      k(1, :) = f(res)
      do i = 1, 11
         y2 = 0
         do j = 1, i
            y2 = y2 + a(j, i)*k(j, :)
         end do
         k(i + 1, :) = f(res + h*y2)
      end do
      k1 = k(12, :)
      k(4, :) = 0
      do i = 1, 12
         k(4, :) = k(4, :) + b8(i)*k(i, :)
      end do
      k(5, :) = y + h*k(4, :)
      y2 = 0
      do i = 1, 12
         y2 = y2 + er(i)*k(i, :)
      end do
      err = 0.0_dp
      err2 = 0.0_dp
      do i = 1, n
         sk = at + rt*max(abs(y(i)), abs(k(5, i)))
         err = err + (y2(i)/sk)**2
         err2 = err2 + ((k(4, i) - DOT_PRODUCT(bhh, k(:, i)))/sk)**2
      end do
      ! err = abs(h)*err/sqrt(err + 0.01*err2)
      err = abs(h)*sqrt(err)
      res = k(5, :)
      deallocate (k, y2)
   end subroutine dop853_step

   subroutine odex_step(f, y, h, k, res, err)
      implicit none
      interface
         function f(y) result(dy)
            import :: dp
            real(dp), intent(in) :: y(:)
            real(dp) :: dy(size(y))
         end function f
      end interface

      integer, intent(inout) :: k
      real(dp), intent(inout) :: h
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: res(:), err

      integer :: i, j, l, n, ni
      real(dp) :: dt, scale, errsum, wk1, wk2, hk1, hk2
      real(dp), allocatable :: T(:, :, :), yprev(:), ycurr(:), ynext(:), fval(:)
      integer, allocatable :: Nsteps(:)

      n = size(y)
      allocate (T(k + 1, k + 1, n), Nsteps(k + 1), yprev(n), ycurr(n), ynext(n), fval(n))
      res = y

      ! Define substep counts n_i = 2*i
      i = 2
      Nsteps(1) = 2
      do
         Nsteps(i) = Nsteps(i - 1)*2
         if (i == k) exit
         Nsteps(i + 1) = Nsteps(i - 1)*3
         if (i + 1 == k) exit
         i = i + 2
      end do

      ! Loop over i = 1 to k to compute T(i,1)
      !$omp parallel do private(i, ni, dt, fval, yprev, ycurr, ynext) shared(T, Nsteps, h, y)
      do i = 1, k
         ni = Nsteps(i)
         dt = h/ni

         fval = f(y)
         yprev = y
         ycurr = y + dt*fval

         do l = 2, ni
            fval = f(ycurr)
            ynext = yprev + 2.0_dp*dt*fval
            yprev = ycurr
            ycurr = ynext
         end do
         fval = f(ycurr)
         T(i, 1, :) = 0.5_dp*(yprev + ycurr + dt*fval)
      end do
      !$omp end parallel do
      ! Richardson extrapolation to fill T(i,j)
      do j = 2, k
         do i = j, k
            T(i, j, :) = T(i, j - 1, :) + (T(i, j - 1, :) - T(i - 1, j - 1, :))/ &
                         ((real(Nsteps(i))/real(Nsteps(i - j + 1)))**2 - 1.0_dp)

            if (i == k - 1 .and. j == k - 1) then
               errsum = 0.0_dp
               do l = 1, n
                  scale = at + rt*max(abs(T(k - 2, k - 2, l)), abs(T(k - 2, k - 3, l)))
                  errsum = errsum + ((T(k - 2, k - 2, l) - T(k - 2, k - 3, l))/scale)**2
               end do
               err = sqrt(errsum/n)
               wk2 = calculate_wk(h, err, k - 2)
               errsum = 0.0_dp
               do l = 1, n
                  scale = at + rt*max(abs(T(k - 1, k - 1, l)), abs(T(k - 1, k - 2, l)))
                  errsum = errsum + ((T(k - 1, k - 1, l) - T(k - 1, k - 2, l))/scale)**2
               end do
               err = sqrt(errsum/n)
               wk1 = calculate_wk(h, err, k - 1)
               hk1 = calculate_hk(h, err, k - 1)
               if (err < 1.0_dp) then
                  res = T(k - 1, k - 1, :)
                  if (wk1 > 0.9*wk2) then
                     k = max(4, k - 1)
                     h = hk1
                  else
                     h = hk1*calculate_ak(k)/calculate_ak(k - 1)
                  end if
                  deallocate (T, Nsteps, yprev, ycurr, ynext, fval)
                  return
               elseif (err > (k*k + 1)**2) then
                  k = max(4, k - 1)
                  h = hk1
                  res = y
                  deallocate (T, Nsteps, yprev, ycurr, ynext, fval)
                  return
               end if
            end if
         end do
      end do
      ! Final result and error estimate
      errsum = 0.0_dp
      do i = 1, n
         scale = at + rt*max(abs(T(k, k, i)), abs(T(k, k - 1, i)))
         errsum = errsum + ((T(k, k, i) - T(k, k - 1, i))/scale)**2
      end do
      err = sqrt(errsum/n)
      hk2 = calculate_hk(h, err, k)
      wk2 = calculate_wk(h, err, k)
      if (err < 1) then
         res = T(k, k, :)
         if (wk1 < 0.9*wk2) then
            k = max(4, k - 1)
            h = hk1
         elseif (wk2 < 0.9*wk1) then
            k = k + 1
            h = hk2*calculate_ak(k + 1)/calculate_ak(k)
         else
            h = hk2
         end if
         deallocate (T, Nsteps, yprev, ycurr, ynext, fval)
         return
      end if
      do i = k + 1, k + 1
         ni = Nsteps(i)
         dt = h/ni

         fval = f(y)
         yprev = y
         ycurr = y + dt*fval

         do l = 2, ni
            fval = f(ycurr)
            ynext = yprev + 2.0_dp*dt*fval
            yprev = ycurr
            ycurr = ynext
         end do

         ! Final smoothing step
         fval = f(ycurr)
         T(i, 1, :) = 0.5_dp*(yprev + ycurr + dt*fval)
      end do
      do j = 2, k + 1
         do i = k + 1, k + 1
            T(i, j, :) = T(i, j - 1, :) + (T(i, j - 1, :) - T(i - 1, j - 1, :))/ &
                         ((real(Nsteps(i))/real(Nsteps(i - j + 1)))**2 - 1.0_dp)
         end do
      end do
      errsum = 0.0_dp
      do i = 1, n
         scale = at + rt*max(abs(T(k + 1, k + 1, i)), abs(T(k + 1, k, i)))
         errsum = errsum + ((T(k + 1, k + 1, i) - T(k + 1, k, i))/scale)**2
      end do
      err = sqrt(errsum/n)
      if (err < 1) then
         res = T(k + 1, k + 1, :)
         if (wk1 < 0.9*wk2) then
            k = max(4, k - 1)
            h = hk1
         elseif (wk1 < 0.9*wk2) then
            hk1 = calculate_hk(h, err, k + 1)
            k = k + 1
            h = hk1
         else
            h = hk2
         end if
         deallocate (T, Nsteps, yprev, ycurr, ynext, fval)
         return
      else
         res = y
         k = max(4, k - 1)
         h = hk1
         deallocate (T, Nsteps, yprev, ycurr, ynext, fval)
         return
      end if

      deallocate (T, Nsteps, yprev, ycurr, ynext, fval)
   end subroutine odex_step
   function calculate_wk(h, er1, k) result(wk)
      implicit none
      real(dp), intent(in)::h, er1
      integer, intent(in)::k
      integer :: kc
      real(dp) :: hk, wk, ak
      kc = max(1, k)
      hk = h*0.94_dp*(0.65_dp/er1)**(1.0_dp/(2.0_dp*kc - 1.0_dp))
      if (er1 == 0) hk = h
      ak = calculate_ak(kc)
      wk = ak/hk
   end function calculate_wk
   function calculate_hk(h, er1, k) result(hk)
      implicit none
      real(dp), intent(in)::h, er1
      integer, intent(in)::k
      integer :: kc
      real(dp) :: hk, wk, ak
      kc = max(1, k)
      hk = h*0.94_dp*(0.65_dp/er1)**(1.0_dp/(2.0_dp*kc - 1.0_dp))
      if (er1 == 0) hk = h
      ak = calculate_ak(kc)
      wk = ak/hk
   end function calculate_hk
   function calculate_ak(k) result(ak)
      implicit none
      integer, intent(in)::k
      integer :: kc, wk(k), i
      real(dp) ::ak
      kc = max(1, k)
      wk = 0
      i = 2
      wk(1) = 2
      do
         wk(i) = wk(i - 1)*2
         if (i == kc) exit
         wk(i + 1) = wk(i - 1)*3
         if (i + 1 == kc) exit
         i = i + 2
      end do
      ak = 1 + sum(wk)
   end function calculate_ak
   subroutine intode(f, y, t, res, error_flag)
      implicit none
      real(dp), intent(in)::y(:), t
      real(dp), intent(out)::res(:)
      logical, intent(out)::error_flag
      logical :: stiff, notlastloop
      real(dp)::h, tc, er1, h_min, d0, d1, temp, t1, t2, tnew
      real(dp), allocatable::yc(:), yf(:), k1(:), kt(:)
      complex(dp) :: hess(n_size - 1, n_size - 1), zc(n_size - 1)
      integer::ct, k
      real(dp) :: h_min_fp, h_min_tol, h_min_span
      real(dp), parameter :: c_fp = 16.0_dp          ! round-off safety
      real(dp), parameter :: c_tol = 0.01_dp          ! 1 % of atol/rtol
      real(dp), parameter :: c_span = 1.0e-12_dp       ! 1 × 10⁻¹² of span
      interface
         function f(y) result(res)
            use, intrinsic::iso_fortran_env, only: real64
            integer, parameter::dp = real64
            real(dp), intent(in)::y(:)
            real(dp)::res(size(y))
         end function f
      end interface
      if (t == 0) then
         res = y
         error_flag = .false.
         return
      end if
      ALLOCATE (yc(size(y)), yf(size(y)), k1(size(y)), kt((size(y))))
      h_min_fp = c_fp*epsilon(1.0_dp)*max(1.0_dp, abs(t))
      h_min_tol = c_tol*1.0_dp
      h_min_span = c_span*abs(t - 0.0_dp)

      h_min = max(h_min_fp, min(h_min_tol, h_min_span))
      tc = 0
      yc = y
      error_flag = .true.
      stiff = .false.
      k1 = 0
      kt = 0
      d0 = norm2(yc)
      h = t
      ct = 0
      k = 4
      notlastloop = .true.
      stiff = .false.
      do while (notlastloop)
      if (tc + h > t .or. tc + h == t) then
         notlastloop = .false.
         h = t - tc
      end if
      tnew = tc + h

      ! call real_to_complex(yc(1:2*size(zc)), zc)
      ! call hessian(zc, hess)
      ! hess = matmul(conjg(transpose(hess)), hess)
      ! call isStiff_Hermitian(hess, n_size - 1, h, stiff, 0.9_dp*6.394_dp)
      ! if (stiff) then
      !    if (ctsf < 3) then
      !       ctsf = ctsf + 1
      !       stiff = .false.
      !    else
      !       ctsf = 0
      !    end if
      ! end if
      if (stiff) then

      else
         call odex_step(f, yc, h, k, yf, er1)
         ! call dop853_step(f, yc, h, yf, er1, k1)
      end if
      if (isnan(norm2(yf))) then
         error_flag = .true.
         deallocate (yc, yf)
         return
      end if
      if (er1 < 1.0) then
         kt = k1
         tc = tnew
         yc = yf
         ct = ct + 1
      else
         notlastloop = .true.
      end if
      k1 = kt
      ! h = h*min(10.0_dp, max(0.1_dp, 0.9_dp*er1**(-1.0_dp/6.0_dp)))
      if (h < h_min .or. ct > 10000 .or. isnan(h)) then
         error_flag = .true.
         deallocate (yc, yf, k1, kt)
         return
      end if
      end do
      ! print *, ct, k
      res = yc
      error_flag = .false.
      deallocate (yc, yf, k1, kt)
   end subroutine intode

   subroutine flowz(x, z, error)
      real(dp), intent(in)::x(:)
      complex(dp), intent(inout)::z(:)
      logical, intent(out)::error
      integer::n
      real(dp), dimension(:), allocatable::y, yf
      real(dp)::t0, t1

      n = size(z)*2
      t0 = 0.0_dp
      t1 = x(1)
      error = .false.

      allocate (y(n), yf(n))

      z = x(2:)
      call complex_to_real(z, y(1:n))
      call intode(rhs_function, y, t1, yf, error)
      y = yf
      call real_to_complex(y(1:), z)
      deallocate (y, yf)

   contains
      function rhs_function(y) result(f)
         real(dp), dimension(:), intent(in)::y
         real(dp), dimension(size(y))::f
         complex(dp), dimension(size(z))::z_temp, ds_val

         call real_to_complex(y(1:), z_temp)
         call ds(z_temp, ds_val)
         call complex_to_real(conjg(ds_val), f(1:))

      end function rhs_function
   end subroutine flowz

   subroutine flowzr(x, z, error)
      real(dp), intent(in)::x(:)
      complex(dp), intent(inout)::z(:)
      logical, intent(out)::error
      integer::n
      real(dp), dimension(:), allocatable::y, yf
      real(dp)::t0, t1

      n = size(z)*2
      t0 = 0.0_dp
      t1 = x(1)
      error = .false.

      allocate (y(n), yf(n))

      call complex_to_real(z, y(1:n))
      call intode(rhs_function, y, t1, yf, error)
      y = yf
      call real_to_complex(y(1:), z)
      deallocate (y, yf)

   contains
      function rhs_function(y) result(f)
         real(dp), dimension(:), intent(in)::y
         real(dp), dimension(size(y))::f
         complex(dp), dimension(size(z))::z_temp, ds_val
         call real_to_complex(y(1:), z_temp)
         call ds(z_temp, ds_val)
         call complex_to_real(-conjg(ds_val), f(1:))

      end function rhs_function
   end subroutine flowzr

   subroutine flow(x, z, j, error)
      real(dp), intent(in)::x(:)
      complex(dp), intent(inout)::z(:)
      complex(dp), dimension(:, :), intent(inout)::j
      logical, intent(out)::error
      integer::n, m
      real(dp), dimension(:), allocatable::y, yf
      real(dp)::t0, t1

      n = size(z)*2
      m = size(j, 1)*size(j, 2)*2
      t0 = 0.0_dp
      t1 = x(1)
      error = .false.

      allocate (y(n + m), yf(n + m))

      z = x(2:)
      j = identity_matrix(n/2)
      call complex_to_real(z, y(1:n))
      call map_to_real(j, y(n + 1:))
      call intode(rhs_function, y, t1, yf, error)
      y = yf
      call real_to_complex(y(1:n), z)
      call map_to_complex(y(n + 1:), j)
      deallocate (y, yf)

   contains
      function rhs_function(y) result(f)
         real(dp), dimension(:), intent(in)::y
         real(dp)::f(size(y)), ti, tf
         complex(dp), dimension(size(z))::z_temp, ds_val
         complex(dp), dimension(size(j, 1), size(j, 2))::j_temp, djdt, j_temp2, j_temp3
         integer :: n2, i1, i2
         real(dp) :: t1, t2
         n2 = size(djdt, 1)
         call real_to_complex(y(1:n), z_temp)
         call map_to_complex(y(n + 1:), j_temp)
         call ds(z_temp, ds_val)
         call complex_to_real(conjg(ds_val), f(1:n))
         call hessian(z_temp, djdt)
         call zgemm('N', 'N', n2, n2, n2, (1.0d0, 0.0d0), djdt, n2, j_temp, n2, (0.0d0, 0.0d0), j_temp2, n2)
         djdt = conjg(j_temp2)
         call map_to_real(djdt, f(n + 1:))

      end function rhs_function
   end subroutine flow

   subroutine flowr(x, z, j, error)
      real(dp), intent(in)::x(:)
      complex(dp), intent(inout)::z(:)
      complex(dp), dimension(:, :), intent(inout)::j
      logical, intent(out)::error
      integer::n, m
      real(dp), dimension(:), allocatable::y, yf
      real(dp)::t0, t1

      n = size(z)*2
      m = size(j, 1)*size(j, 2)*2
      t0 = 0.0_dp
      t1 = x(1)
      error = .false.

      allocate (y(n + m), yf(n + m))

      call complex_to_real(z, y(1:n))
      call map_to_real(j, y(n + 1:))
      call intode(rhs_function, y, t1, yf, error)
      y = yf
      call real_to_complex(y(1:n), z)
      call map_to_complex(y(n + 1:), j)
      deallocate (y, yf)

   contains
      function rhs_function(y) result(f)
         real(dp), dimension(:), intent(in)::y
         real(dp)::f(size(y)), ti, tf
         complex(dp), dimension(size(z))::z_temp, ds_val
         complex(dp), dimension(size(j, 1), size(j, 2))::j_temp, djdt, j_temp2
         integer :: n2
         n2 = size(djdt, 1)
         call real_to_complex(y(1:n), z_temp)
         call map_to_complex(y(n + 1:), j_temp)
         call ds(z_temp, ds_val)
         call complex_to_real(-conjg(ds_val), f(1:n))
         call hessian(z_temp, djdt)
         call zgemm('N', 'N', n2, n2, n2, (1.0d0, 0.0d0), djdt, n2, j_temp, n2, (0.0d0, 0.0d0), j_temp2, n2)
         djdt = conjg(j_temp2)
         call map_to_real(-djdt, f(n + 1:))

      end function rhs_function
   end subroutine flowr
   subroutine flowv(x, uv, ju, jv, error)
      real(dp), intent(in)::x(:)
      real(dp), intent(in)::uv(:)
      real(dp), dimension(:), intent(out)::ju, jv
      logical, intent(out)::error
      complex(dp), dimension(:), allocatable::z
      integer::n
      real(dp), dimension(:), allocatable::y, yf
      real(dp)::t0, t1

      n = (size(x) - 1)*2
      t0 = 0.0_dp
      t1 = x(1)
      error = .false.

      allocate (y(3*n), z(n/2), yf(3*n))

      if (size(uv) /= n) then
         write (*, *) "Error:v must be 2n."
         return
      end if

      z = x(2:)
      ju = uv
      call real_vec(ju)
      jv = uv - ju
      call complex_to_real(z, y(1:n))
      y(n + 1:2*n) = ju
      y(2*n + 1:3*n) = jv
      call intode(rhs_function, y, t1, yf, error)
      y = yf
      call real_to_complex(y(1:n), z)
      ju = y(n + 1:2*n)
      jv = y(2*n + 1:3*n)
      deallocate (y, yf, z)

   contains
      function rhs_function(y) result(f)
         real(dp), dimension(:), intent(in)::y
         real(dp)::f(size(y))
         complex(dp), dimension(size(z))::z_temp, ds_val, ju_temp, jv_temp
         complex(dp), dimension(size(z), size(z))::h

         call real_to_complex(y(1:n), z_temp)
         call real_to_complex(y(n + 1:2*n), ju_temp)
         call real_to_complex(y(2*n + 1:3*n), jv_temp)

         call ds(z_temp, ds_val)
         call complex_to_real(conjg(ds_val), f(1:n))
         call hessian(z_temp, h)
         call complex_to_real(conjg(matmul(h, ju_temp)), f(n + 1:2*n))
         call complex_to_real(-conjg(matmul(h, jv_temp)), f(2*n + 1:3*n))

      end function rhs_function
   end subroutine flowv

   subroutine isStiff_Hermitian(J, n, h, stiff, Cstab)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: h
      complex(dp), intent(in) :: J(n, n)   ! Hermitian
      logical, intent(out) :: stiff
      real(dp), optional, intent(in) :: Cstab

      ! LAPACK constants
      character(len=1), parameter :: JOBZ = 'N', RANGE = 'A', UPLO = 'U'
      real(dp), parameter :: ABSTOL = 0.0_dp

      ! Work arrays
      complex(dp), allocatable :: Jcopy(:, :), work(:), Zdummy(:, :)
      real(dp), allocatable :: eig(:), rwork(:)
      integer, allocatable :: iwork(:), isuppz(:)
      integer :: lwork, lrwork, liwork, m, info
      real(dp) :: lambda_min, lambda_max, thresh

      ! Allocate and copy
      allocate (Jcopy(n, n))
      Jcopy = J

      ! Workspace query
      allocate (work(1), rwork(1), iwork(1), eig(2), Zdummy(1, 1), isuppz(4))
      call zheevr(JOBZ, RANGE, UPLO, n, Jcopy, n, 0.0_dp, 0.0_dp, 0, 0, &
                  ABSTOL, m, eig, Zdummy, 1, isuppz, work, -1, rwork, -1, &
                  iwork, -1, info)
      if (info /= 0) then
         stiff = .true.
         deallocate (work, rwork, iwork, eig, Zdummy, isuppz, Jcopy)
         return
      end if

      lwork = int(real(work(1), dp))
      lrwork = int(rwork(1))
      liwork = iwork(1)

      deallocate (work, rwork, iwork, eig, Zdummy, isuppz)

      ! Allocate real work arrays
      allocate (work(lwork), rwork(lrwork), iwork(liwork), eig(n), Zdummy(1, 1), isuppz(2*n))

      Jcopy = J
      call zheevr(JOBZ, RANGE, UPLO, n, Jcopy, n, 0.0_dp, 0.0_dp, 0, 0, &
                  ABSTOL, m, eig, Zdummy, 1, isuppz, work, lwork, rwork, &
                  lrwork, iwork, liwork, info)

      if (info /= 0 .or. m /= n) then
         stiff = .true.
      else
         lambda_min = abs(eig(1))
         lambda_max = abs(eig(n))
         thresh = merge(Cstab, 1.0_dp, present(Cstab))
         stiff = sqrt(lambda_max/lambda_min) > 1.0e4_dp
      end if

      ! Clean up
      deallocate (work, rwork, iwork, eig, Zdummy, isuppz, Jcopy)
   end subroutine isStiff_Hermitian

end module solve_flow
