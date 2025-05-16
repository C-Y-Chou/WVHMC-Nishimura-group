module quasi_newton_module
   use utils
   use param_mod
   use model
   use solve_flow
   use mpi

contains

   subroutine calculate_fqv2(xt, z, xi, fq, del_z, ierr, Jl, jac)
      implicit none
      real(dp), intent(in) :: xt(:), xi(:), del_z(:)
      complex(dp), intent(in) :: z(:), jac(:, :)
      real(dp), intent(out) :: fq(:), Jl(:)
      logical, intent(out) :: ierr
      complex(dp), allocatable :: ds_val(:), E0(:), Jlc(:), z_new(:)
      real(dp), allocatable  ::  aux(:), lambdai(:), z_real(:), f(:)
      ALLOCATE (ds_val(size(z)), E0(size(z)), aux(2*size(z)), Jlc(size(z)), z_new(size(z)))
      ALLOCATE (lambdai(2*size(z)), z_real(2*size(z)))
      Jlc = matmul(jac, cmplx(0, 1)*xi(1:size(z)) + xi(size(z) + 1:))
      call complex_to_real(Jlc, Jl)
      call real_to_complex(del_z, z_new)
      z_new = z + z_new + Jlc
      call flowzr(xt, z_new, ierr)
      if (ierr) then
         call flowz(xt, z_new, ierr)
         fq = 1.0E+10
         ierr = .true.
         deallocate (ds_val, E0, aux, Jlc, z_new, lambdai, z_real)
         return
      end if
      fq(1:size(z)) = aimag(z_new)
      fq(size(z) + 1:) = xi(size(z) + 1:)
      deallocate (ds_val, E0, aux, Jlc, z_new, lambdai, z_real)
   end subroutine calculate_fqv2

   subroutine calculate_fqm(xt, z, xi, fq, fqinv, del_z, ierr, Jl, jac)
      implicit none
      real(dp), intent(in) :: xt(:), xi(:), del_z(:)
      complex(dp), intent(in) :: z(:), jac(:, :)
      real(dp), intent(out) :: fq(:, :), Jl(:), fqinv(:, :)
      logical, intent(out) :: ierr
      complex(dp), allocatable :: ds_val(:), E0(:), Jlc(:), z_new(:), jprime(:, :)
      real(dp), allocatable  ::  aux(:), lambdai(:), z_real(:), f(:)
      integer :: i, j, work
      real(dp) :: diff, error
      integer :: info
      INTEGER, ALLOCATABLE :: ipiv(:)
      ALLOCATE (ds_val(size(z)), E0(size(z)), aux(2*size(z)), Jlc(size(z)), z_new(size(z)))
      ALLOCATE (lambdai(2*size(z)), z_real(2*size(z)), jprime(size(z), size(z)))
      Jlc = matmul(jac, cmplx(0, 1)*xi(1:size(z)) + xi(size(z) + 1:))
      call complex_to_real(Jlc, Jl)
      call real_to_complex(del_z, z_new)
      z_new = z + z_new + Jlc
      if (ierr) then
         fq = 1.0E+10
         ierr = .true.
         deallocate (ds_val, E0, aux, Jlc, z_new, lambdai, z_real, jprime)
         return
      end if
      fq = 0
      do i = size(z) + 1, 2*size(z)
         fq(i, i) = 1
      end do
      jprime = cmplx(0, 1)*jac
      call flowr(xt, z_new, jprime, ierr)
      fq(1:size(z), 1:size(z)) = aimag(jprime)
      jprime = jac
      call real_to_complex(del_z, z_new)
      z_new = z + z_new + Jlc
      call flowr(xt, z_new, jprime, ierr)
      fq(1:size(z), size(z) + 1:) = aimag(jprime)
      deallocate (ds_val, E0, aux, Jlc, z_new, lambdai, z_real, jprime)
      ! Compute the inverse of fq using LAPACK (fqinv = inv(fq))
      allocate (ipiv(size(fq, 1)))

      ! Copy fq to fqinv since LAPACK will overwrite the input
      fqinv = fq

      ! LU factorization of fqinv
      ! call dgetrf(size(fq, 1), size(fq, 2), fqinv, size(fq, 1), ipiv, info)
      ! if (info /= 0) then
      !    ierr = .true.
      !    fqinv = 1.0E+10
      !    deallocate (ipiv)
      !    return
      ! end if
      ! ALLOCATE (aux(1))
      ! Compute the inverse using the LU factorization
      ! call dgetri(size(fq, 1), fqinv, size(fq, 1), ipiv, aux, -1, info)
      ! if (info /= 0) then
      !    ierr = .true.
      !    fqinv = 1.0E+10
      ! end if
      ! work = aux(1)
      ! DEALLOCATE (aux)
      ! ALLOCATE (aux(work))
      ! call dgetri(size(fq, 1), fqinv, size(fq, 1), ipiv, aux, work, info)
      ! if (info /= 0) then
      !    ierr = .true.
      !    fqinv = 1.0E+10
      ! end if
      ! deallocate (ipiv, aux)

      deallocate (ipiv)
   end subroutine calculate_fqm

   subroutine quasi_newtonv2(f, tol, max_iter, xt, z, del_z, ierr, Jl, x_new, jac)
      implicit none

      ! Arguments
      integer, intent(in)          :: max_iter
      real(dp), intent(in)         :: tol
      logical, intent(out)         :: ierr
      real(dp), intent(in) :: xt(:), del_z(:)
      complex(dp), intent(in) :: z(:)
      real(dp), intent(out) :: Jl(:)
      complex(dp), intent(in) :: jac(:, :)
      real(dp) :: ti, tf

      interface
         subroutine f(xt, z, xi, fq, del_z, ierr, Jl, jac)
            use, intrinsic :: iso_fortran_env, only: real64
            integer, parameter :: dp = real64
            real(dp), intent(in) :: xt(:), xi(:), del_z(:)
            complex(dp), intent(in) :: z(:), jac(:, :)
            real(dp), intent(out) :: fq(:), Jl(:)
            logical, intent(out) :: ierr
         end subroutine f
      end interface

      ! Local variables
      integer :: n, iter, i
      real(dp) :: alpha
      real(dp), ALLOCATABLE :: x(:)
      real(dp), allocatable :: dx(:), Hm(:, :), d(:), y(:), grad(:), grad_new(:), fx_val(:), Bm(:, :)
      complex(dp), allocatable :: tempj(:, :)
      real(dp), INTENT(OUT):: x_new(:)
      complex(dp):: z_new(size(z)), z_new2(size(z))
      integer, allocatable :: ipiv(:)
      integer :: info, ct
      real(dp) :: z_new_real(size(z)*2), rho, delta, delta_max, tempf, fx_max, fx_his(30)
      real(dp) :: gamma, rhok, ldk, sgm, dta, etak, lk, phi, tk
      logical :: converged
      ti = MPI_Wtime()
      gamma = 0.5
      rhok = 0.5
      sgm = 0.5
      dta = 0.25
      n = 2*size(z)
      allocate (x(n), tempj(n/2, n/2), ipiv(n/2))
      allocate (Hm(n, n), Bm(n, n), dx(n), d(n), y(n), grad(n), grad_new(n), fx_val(n))

      call real_to_complex(-del_z, z_new)
      tempj = jac
      call zgesv(n/2, 1, tempj, n/2, ipiv, z_new, n/2, info)
      if (info /= 0) z_new = 0
      DEALLOCATE (ipiv)
      x(1:n/2) = AIMAG(z_new)
      x(n/2 + 1:) = real(z_new)
      call f(xt, z, x, fx_val, del_z, ierr, Jl, jac)
      phi = norm2(fx_val)
      ! Evaluate f at the initial guess.
      call calculate_fqm(xt, z, x, Bm, Hm, del_z, ierr, Jl, jac)
      if (ierr) then
         print *, "bad2"
         Bm = 0
         do i = 1, n
            Bm(i, i) = 1
         end do
      end if

      ! call numerical_jacobian(f, xt, z, x, fx_val, del_z, Jl, jac, Bm)
      iter = 0
      ierr = .false.
      !-----------------------------------------------------------------
      ! Main iteration loop
      !-----------------------------------------------------------------
      ct = 0
      do while (ct < max_iter)
         iter = iter + 1
         etak = 1/real(iter + 1, dp)**2
         call inv(fx_val, Bm, d, ierr)
         y = fx_val
         call f(xt, z, x + d, fx_val, del_z, ierr, Jl, jac)
         converged = (norm2(fx_val) < (gamma*norm2(y) - rhok*norm2(d)**2))
         if (ierr) converged = .false.
         ldk = 1
         lk = 0
         do while (.not. converged)
            lk = lk + 1
            ldk = dta**lk
            call f(xt, z, x + ldk*d, fx_val, del_z, ierr, Jl, jac)
            converged = (norm2(fx_val) < ((1.0_dp + etak)*phi - sgm*norm2(dta**lk*d)**2))
            if (ierr) converged = .false.
            if (lk > 20 .and. norm2(y) > 1d-8) then
               ierr = .true.
               deallocate (Hm, Bm, dx, d, y, grad, grad_new, tempj, x)
               return
            end if
         end do
         delta = norm2(ldk*d)
         if (lk > 5) then
            ct = ct + 1
         end if
         dx = ldk*d

         y = fx_val - y
         tk = ((1.0_dp + etak)*phi + 1.0_dp)*norm2(fx_val)/(norm2(fx_val) + 1.0_dp)
         phi = 0.7*tk + 0.3*norm2(fx_val)
         ! Update x and evaluate f(x).
         x = x + dx

         Bm = Bm + outer_product((y - matmul(Bm, dx)), dx)/dot_product(dx, dx)

         tf = MPI_Wtime()
         if (mod(ct, 10) == 0 .and. lk > 5) then
            ! print *, "restart"
            call calculate_fqm(xt, z, x, Bm, Hm, del_z, ierr, Jl, jac)
         end if
         ! print *, "Iteration:", iter, "time:", tf - ti, "||f(x)|| =", norm2(fx_val), "|d|", norm2(dx), "lk", lk
         if (norm2(fx_val) < tol) then
            call real_to_complex(del_z + Jl, z_new)
            z_new = z + z_new
            call flowzr(xt, z_new, ierr)
            x_new = xt
            x_new(2:) = real(z_new)
            call flowz(x_new, z_new2, ierr)
            call real_to_complex(del_z + Jl, z_new)
            z_new = z + z_new
            ! print*, sqrt(dot_product(z_new2-z_new,z_new2-z_new)),"real error"
            deallocate (Hm, Bm, dx, d, y, grad, grad_new, tempj, x)
            ierr = .false.
            return
         end if
      end do
      ierr = .true.
      deallocate (Hm, Bm, dx, d, y, grad, grad_new, tempj, x)
   end subroutine quasi_newtonv2
   subroutine dogleg(Fval, B, Hm, delta, d, ierr)
      implicit none

      !------------------------------------------------------------------
      ! Inputs
      !------------------------------------------------------------------
      real(dp), intent(in)  :: Fval(:)
      ! B is assumed to be size(N,N).
      real(dp), intent(in)  :: B(:, :), Hm(:, :)
      real(dp), intent(in)  :: delta    ! Trust-region radius

      !------------------------------------------------------------------
      ! Outputs
      !------------------------------------------------------------------
      real(dp), intent(out) :: d(:)     ! Computed dogleg step
      logical, intent(out)  :: ierr     ! Error flag

      !------------------------------------------------------------------
      ! Local variables
      !------------------------------------------------------------------
      integer            :: n, info
      integer, allocatable :: ipiv(:)
      real(dp), allocatable :: g(:), d_sd(:), d_newton(:), H(:, :), Hfact(:, :), p(:)
      real(dp)           :: alpha, norm_sd, norm_newton
      real(dp)           :: a_coef, b_coef, c_coef, disc, tau
      real(dp)           :: dot_val

      !------------------------------------------------------------------
      ! Initialization
      !------------------------------------------------------------------
      ierr = .false.
      n = size(Fval)
      allocate (g(n), d_sd(n), d_newton(n), p(n), ipiv(n))
      allocate (H(n, n), Hfact(n, n))

      !------------------------------------------------------------------
      ! 1. Build gradient and Hessian approximation
      !    If Fval = (F - B*d), then the negative gradient is:
      !         - grad(phi) = B^T * Fval   (up to a factor of 2).
      !    We'll store g as the gradient for convenience:
      !         g = - (B^T Fval).
      !    Then H ~ B^T B.
      !------------------------------------------------------------------
      g = matmul(transpose(B), Fval)    ! "gradient" ignoring factor of 2
      H = matmul(transpose(B), B)        ! approximate Hessian ignoring factor of 2

      !------------------------------------------------------------------
      ! 2. Steepest-descent (Cauchy) step
      !    d_sd = - ( g^T g / (g^T H g) ) * g
      !------------------------------------------------------------------
      dot_val = norm2(matmul(transpose(B), Fval))
      alpha = dot_val/norm2(matmul(B, matmul(transpose(B), Fval)))
      d_sd = -alpha*g
      norm_sd = sqrt(dot_product(d_sd, d_sd))
      !------------------------------------------------------------------
      ! 3. Full Newton step:
      !    Solve  H * d_newton = -g   with one DGESV call
      !------------------------------------------------------------------
      d_newton = -Fval             ! RHS now contains  -g
      Hfact = B             ! DGESV overwrites the matrix with its LU factors
      call dgesv(n, 1, Hfact, n, ipiv, d_newton, n, info)

      if (info /= 0) then       ! info>0 ⇒ singular; info<0 ⇒ illegal arg
         ierr = .true.
         return
      end if
      norm_newton = norm2(d_newton)
      !------------------------------------------------------------------
      ! 4. Dogleg logic with trust-region radius 'delta'
      !------------------------------------------------------------------
      if (norm_newton <= delta) then
         ! Full Newton step is within the trust region
         d = d_newton

      else if (norm_sd >= delta) then
         ! Even the steepest-descent step is too big; scale it down
         d = (delta/norm_sd)*d_sd

      else
         ! Cauchy step is inside delta, but Newton step is outside
         ! --> find a tau in [0,1] s.t. ||d_sd + tau*(d_newton - d_sd)|| = delta
         p = d_newton - d_sd
         a_coef = dot_product(p, p)
         b_coef = 2.0_dp*dot_product(d_sd, p)
         c_coef = dot_product(d_sd, d_sd) - delta**2

         disc = b_coef**2 - 4.0_dp*a_coef*c_coef
         if (disc < 0.0_dp) then
            ! Should be rare for a well-defined trust region problem
            ierr = .true.
            d = d_sd
         else
            tau = (-b_coef + sqrt(disc))/(2.0_dp*a_coef)
            d = d_sd + tau*p
         end if
      end if

      !------------------------------------------------------------------
      ! Cleanup
      !------------------------------------------------------------------
      deallocate (g, d_sd, d_newton, p, ipiv, H, Hfact)

   end subroutine dogleg
   subroutine inv(Fval, B, d, ierr)
      implicit none

      !------------------------------------------------------------------
      ! Inputs
      !------------------------------------------------------------------
      real(dp), intent(in)  :: Fval(:)
      ! B is assumed to be size(N,N).
      real(dp), intent(in)  :: B(:, :)

      !------------------------------------------------------------------
      ! Outputs
      !------------------------------------------------------------------
      real(dp), intent(out) :: d(:)     ! Computed dogleg step
      logical, intent(out)  :: ierr     ! Error flag

      !------------------------------------------------------------------
      ! Local variables
      !------------------------------------------------------------------
      integer            :: n, info
      integer, allocatable :: ipiv(:)
      real(dp), allocatable :: g(:), d_sd(:), d_newton(:), H(:, :), Hfact(:, :), p(:)
      real(dp)           :: alpha, norm_sd, norm_newton
      real(dp)           :: a_coef, b_coef, c_coef, disc, tau
      real(dp)           :: dot_val

      !------------------------------------------------------------------
      ! Initialization
      !------------------------------------------------------------------
      ierr = .false.
      n = size(Fval)
      allocate (g(n), d_sd(n), d_newton(n), p(n), ipiv(n))
      allocate (H(n, n), Hfact(n, n))

      d = -Fval             ! RHS now contains  -g
      Hfact = B             ! DGESV overwrites the matrix with its LU factors
      call dgesv(n, 1, Hfact, n, ipiv, d, n, info)

      if (info /= 0) then       ! info>0 ⇒ singular; info<0 ⇒ illegal arg
         ierr = .true.
         return
      end if

      deallocate (g, d_sd, d_newton, p, ipiv, H, Hfact)

   end subroutine inv

   subroutine quasi_newtonwv(f, tol, max_iter, xt, z, del_z, ierr, Jl, x_new, jac)
      implicit none

      ! Arguments
      integer, intent(in)          :: max_iter
      real(dp), intent(in)         :: tol
      logical, intent(out)         :: ierr
      real(dp), intent(in) :: xt(:), del_z(:)
      complex(dp), intent(in) :: z(:)
      real(dp), intent(out) :: Jl(:)
      complex(dp), intent(in) :: jac(:, :)
      real(dp) :: ti, tf

      interface
         subroutine f(xt, z, xi, fq, del_z, ierr, Jl, jac)
            use, intrinsic :: iso_fortran_env, only: real64
            integer, parameter :: dp = real64
            real(dp), intent(in) :: xt(:), xi(:), del_z(:)
            complex(dp), intent(in) :: z(:), jac(:, :)
            real(dp), intent(out) :: fq(:), Jl(:)
            logical, intent(out) :: ierr
         end subroutine f
      end interface

      ! Local variables
      integer :: n, iter, i
      real(dp) :: alpha
      real(dp), ALLOCATABLE :: x(:)
      real(dp), allocatable :: dx(:), Hm(:, :), d(:), y(:), grad(:), grad_new(:), fx_val(:), Bm(:, :)
      complex(dp), allocatable :: tempj(:, :)
      real(dp), INTENT(OUT):: x_new(:)
      complex(dp):: z_new(size(z))
      integer, allocatable :: ipiv(:)
      integer :: info, ct
      real(dp) :: z_new_real(size(z)*2), rho, delta, delta_max, tempf, fx_max, fx_his(30)
      real(dp) :: gamma, rhok, ldk, sgm, dta, etak, lk, phi, tk
      logical :: converged
      real(dp), allocatable :: Jnum(:, :)
      gamma = 0.5
      rhok = 0.5
      sgm = 0.5
      dta = 0.25
      n = 2*size(z)
      allocate (x(n + 1), tempj(n/2, n/2), ipiv(n/2))
      allocate (Hm(n + 1, n + 1), Bm(n + 1, n + 1), dx(n + 1), d(n + 1), y(n + 1), fx_val(n + 1))

      call real_to_complex(-del_z, z_new)
      tempj = jac
      call zgesv(n/2, 1, tempj, n/2, ipiv, z_new, n/2, info)
      if (info /= 0) z_new = 0
      DEALLOCATE (ipiv)
      x = 0
      x(1:n/2) = AIMAG(z_new)
      x(n/2 + 1:n) = real(z_new)
      call f(xt, z, x, fx_val, del_z, ierr, Jl, jac)
      phi = norm2(fx_val)

      call calculate_fqmwv(xt, z, x, Bm, Hm, del_z, ierr, Jl, jac)
      if (ierr) then
         print *, "bad2"
         Bm = 0
         do i = 1, n+1
            Bm(i, i) = 1
         end do
      end if

      ! call numerical_jacobian(f, xt, z, x, fx_val, del_z, Jl, jac, Bm)
      ct = 0
      iter = 0
      call CPU_TIME(ti)
      do while (iter < max_iter)
         iter = iter + 1
         etak = 1/real(iter + 1, dp)**2
         call inv(fx_val, Bm, d, ierr)
         y = fx_val
         call f(xt, z, x + d, fx_val, del_z, ierr, Jl, jac)
         converged = (norm2(fx_val) < (gamma*norm2(y) - rhok*norm2(d)**2))
         if (ierr) converged = .false.
         ldk = 1
         lk = 0
         do while (.not. converged)
            lk = lk + 1
            ldk = dta**lk
            call f(xt, z, x + ldk*d, fx_val, del_z, ierr, Jl, jac)
            converged = (norm2(fx_val) < ((1.0_dp + etak)*phi - sgm*norm2(dta**lk*d)**2))
            if (ierr) converged = .false.
            if (lk > 20 .and. norm2(y) > 1d-8) then
               ierr = .true.
               deallocate (Hm, Bm, dx, d, y, tempj, x)
               return
            end if
         end do
         if (lk > 5) then
            ct = ct + 1
         end if
         dx = ldk*d

         y = fx_val - y
         tk = ((1.0_dp + etak)*phi + 1.0_dp)*norm2(fx_val)/(norm2(fx_val) + 1.0_dp)
         phi = 0.7*tk + 0.3*norm2(fx_val)
         ! Update x and evaluate f(x).
         x = x + dx

         Bm = Bm + outer_product((y - matmul(Bm, dx)), dx)/dot_product(dx, dx)

         call CPU_TIME(tf)
         ! if (mod(ct, 10) == 0 .and. lk>5) then
         !    ! print *, "restart"
         !    call calculate_fqmwv(xt, z, x, Bm, Hm, del_z, ierr, Jl, jac)
         ! end if

         ! print *, "Iteration:", iter, "time:", tf - ti, "||f(x)|| =", norm2(fx_val), "|d|", norm2(dx), "lk", lk, "x0", xt(1) + x(n+1)
         if (norm2(fx_val) < tol) then
            call real_to_complex(del_z + Jl, z_new)
            z_new = z + z_new
            call flowzr(xt + x(n + 1), z_new, ierr)
            x_new = xt + x(n + 1)
            x_new(2:) = real(z_new)
            deallocate (Hm, Bm, dx, d, y, tempj, x)
            ierr = .false.
            return
         end if
      end do
      ierr = .true.
      deallocate (Hm, Bm, dx, d, y, tempj, x)
   end subroutine quasi_newtonwv

   subroutine calculate_fqwv(xt, z, xi, fq, del_z, ierr, Jl, jac)
      implicit none
      real(dp), intent(in) :: xt(:), xi(:), del_z(:)
      complex(dp), intent(in) :: z(:), jac(:, :)
      real(dp), intent(out) :: fq(:), Jl(:)
      logical, intent(out) :: ierr
      complex(dp), allocatable :: ds_val(:), E0(:), Jlc(:), z_new(:)
      real(dp), allocatable  ::  aux(:), lambdai(:), z_real(:), f(:)
      ALLOCATE (ds_val(size(z)), E0(size(z)), aux(2*size(z)), Jlc(size(z)), z_new(size(z)))
      ALLOCATE (lambdai(2*size(z)), z_real(2*size(z)))
      fq = 0
      Jlc = matmul(jac, cmplx(0, 1)*xi(1:size(z)) + xi(size(z) + 1:))
      call ds(z, ds_val)
      fq(2*size(z) + 1) = real(sum(ds_val*Jlc))
      call complex_to_real(Jlc, Jl)
      call real_to_complex(del_z, z_new)
      z_new = z + z_new + Jlc
      call flowzr(xt + xi(2*size(z) + 1), z_new, ierr)
      if (ierr) then
         fq = 1.0E+10
         ierr = .true.
         deallocate (ds_val, E0, aux, Jlc, z_new, lambdai, z_real)
         return
      end if
      fq(1:size(z)) = aimag(z_new)
      fq(size(z) + 1:2*size(z)) = xi(size(z) + 1:2*size(z))
      deallocate (ds_val, E0, aux, Jlc, z_new, lambdai, z_real)
   end subroutine calculate_fqwv

   subroutine calculate_fqmwv(xt, z, xi, fq, fqinv, del_z, ierr, Jl, jac)
      implicit none
      real(dp), intent(in) :: xt(:), xi(:), del_z(:)
      complex(dp), intent(in) :: z(:), jac(:, :)
      real(dp), intent(out) :: fq(:, :), Jl(:), fqinv(:, :)
      logical, intent(out) :: ierr
      complex(dp), allocatable :: ds_val(:), E0(:), Jlc(:), z_new(:), jprime(:, :)
      real(dp), allocatable  ::  aux(:), lambdai(:), z_real(:), f(:)
      integer :: i, j, work
      real(dp) :: diff, error
      integer :: info
      INTEGER, ALLOCATABLE :: ipiv(:)
      ALLOCATE (ds_val(size(z)), E0(size(z)), aux(2*size(z)), Jlc(size(z)), z_new(size(z)))
      ALLOCATE (lambdai(2*size(z)), z_real(2*size(z)), jprime(size(z), size(z)))
      fq = 0
      Jlc = matmul(jac, cmplx(0, 1)*xi(1:size(z)) + xi(size(z) + 1:))
      call complex_to_real(Jlc, Jl)
      call real_to_complex(del_z, z_new)
      z_new = z + z_new + Jlc
      call ds(z_new, ds_val)
      ds_val = matmul(jac, ds_val)
      fq(2*size(z) + 1, 1:size(z)) = aimag(ds_val)
      if (ierr) then
         fq = 1.0E+10
         ierr = .true.
         deallocate (ds_val, E0, aux, Jlc, z_new, lambdai, z_real, jprime)
         return
      end if

      do i = 1, 2*size(z)
         fq(i, i) = 1
      end do
      ! do i = size(z) + 1, 2*size(z)
      !    fq(i, i) = 1
      ! end do
      ! jprime = cmplx(0, 1)*jac
      ! call flowr(xt, z_new, jprime, ierr)
      call flowzr(xt + xi(2*size(z) + 1), z_new, ierr)
      call ds(z_new, ds_val)
      fq(1:size(z), 2*size(z) + 1) = -aimag(ds_val)
      ! fq(1:size(z), 1:size(z)) = aimag(jprime)
      ! jprime = jac
      ! call real_to_complex(del_z, z_new)
      ! z_new = z + z_new + Jlc
      ! call flowr(xt, z_new, jprime, ierr)
      ! fq(1:size(z), size(z) + 1:2*size(z)) = aimag(jprime)
      ! fqinv = fq
      deallocate (ds_val, E0, aux, Jlc, z_new, lambdai, z_real, jprime)
   end subroutine calculate_fqmwv
   subroutine numerical_jacobian(f, xt, z, x, fx_val, del_z, Jl, jac, J_num)
      implicit none
      integer, parameter :: dp = real64
      real(dp), intent(in) :: xt(:), x(:), del_z(:)
      complex(dp), intent(in) :: z(:), jac(:, :)
      real(dp), intent(in) :: fx_val(:)
      real(dp), intent(out) :: J_num(:, :)
      real(dp), intent(inout) :: Jl(:)
      logical :: ierr
      integer :: i, n
      real(dp) :: h
      real(dp), allocatable :: x_pert(:)
      real(dp), allocatable :: f_plus2h(:), f_plush(:), f_minush(:), f_minus2h(:)

      interface
         subroutine f(xt, z, xi, fq, del_z, ierr, Jl, jac)
            use, intrinsic :: iso_fortran_env, only: real64
            integer, parameter :: dp = real64
            real(dp), intent(in) :: xt(:), xi(:), del_z(:)
            complex(dp), intent(in) :: z(:), jac(:, :)
            real(dp), intent(out) :: fq(:), Jl(:)
            logical, intent(out) :: ierr
         end subroutine f
      end interface

      n = size(x)
      allocate (x_pert(n), f_plus2h(n + 1), f_plush(n + 1), f_minush(n + 1), f_minus2h(n + 1))
      J_num = 0.0_dp
      h = 1.0d-6

      do i = 1, n
         ! f(x + 2h)
         x_pert = x
         x_pert(i) = x(i) + 2*h
         call f(xt, z, x_pert, f_plus2h, del_z, ierr, Jl, jac)
         if (ierr) then
            print *, "Error in f(x+2h) during Jacobian estimation"
            return
         end if

         ! f(x + h)
         x_pert = x
         x_pert(i) = x(i) + h
         call f(xt, z, x_pert, f_plush, del_z, ierr, Jl, jac)
         if (ierr) then
            print *, "Error in f(x+h) during Jacobian estimation"
            return
         end if

         ! f(x - h)
         x_pert = x
         x_pert(i) = x(i) - h
         call f(xt, z, x_pert, f_minush, del_z, ierr, Jl, jac)
         if (ierr) then
            print *, "Error in f(x-h) during Jacobian estimation"
            return
         end if

         ! f(x - 2h)
         x_pert = x
         x_pert(i) = x(i) - 2*h
         call f(xt, z, x_pert, f_minus2h, del_z, ierr, Jl, jac)
         if (ierr) then
            print *, "Error in f(x-2h) during Jacobian estimation"
            return
         end if

         ! 5-point stencil
         J_num(:, i) = (-f_plus2h + 8*f_plush - 8*f_minush + f_minus2h)/(12*h)
      end do

      deallocate (x_pert, f_plus2h, f_plush, f_minush, f_minus2h)
   end subroutine numerical_jacobian

   subroutine quasi_newtons(f, tol, max_iter, xt, z, del_z, ierr, Jl, x_new, jac)
      implicit none

      ! Arguments
      integer, intent(in)          :: max_iter
      real(dp), intent(in)         :: tol
      logical, intent(out)         :: ierr
      real(dp), intent(in) :: xt(:), del_z(:)
      complex(dp), intent(in) :: z(:)
      real(dp), intent(out) :: Jl(:)
      complex(dp), intent(in) :: jac(:, :)
      real(dp) :: ti, tf

      interface
         subroutine f(xt, z, xi, fq, del_z, ierr, Jl, jac)
            use, intrinsic :: iso_fortran_env, only: real64
            integer, parameter :: dp = real64
            real(dp), intent(in) :: xt(:), xi(:), del_z(:)
            complex(dp), intent(in) :: z(:), jac(:, :)
            real(dp), intent(out) :: fq(:), Jl(:)
            logical, intent(out) :: ierr
         end subroutine f
      end interface

      ! Local variables
      integer :: n, iter, i
      real(dp) :: alpha
      real(dp), ALLOCATABLE :: x(:)
      real(dp), allocatable :: dx(:), Hm(:, :), d(:), y(:), grad(:), grad_new(:), fx_val(:), Bm(:, :)
      complex(dp), allocatable :: tempj(:, :)
      real(dp), INTENT(OUT):: x_new(:)
      complex(dp):: z_new(size(z))
      integer, allocatable :: ipiv(:)
      integer :: info, ct
      real(dp) :: z_new_real(size(z)*2), rho, delta, delta_max, tempf, fx_max, fx_his(30)
      real(dp) :: gamma, rhok, ldk, sgm, dta, etak, lk, phi, tk
      logical :: converged
      ti = MPI_Wtime()
      gamma = 0.5
      rhok = 0.5
      sgm = 0.5
      dta = 0.25
      n = 2*size(z)
      allocate (x(n), tempj(n/2, n/2), ipiv(n/2))
      allocate (Hm(n, n), Bm(n, n), dx(n), d(n), y(n), grad(n), grad_new(n), fx_val(n))

      x = 0
      fx_val = -del_z
      phi = norm2(fx_val)
      ! Evaluate f at the initial guess.
      call calculate_fqms(xt, z, x, Bm, Hm, del_z, ierr, Jl, jac)
      ! call numerical_jacobian(f, xt, z, x, fx_val, del_z, Jl, jac, Bm)
      iter = 0
      ierr = .false.
      !-----------------------------------------------------------------
      ! Main iteration loop
      !-----------------------------------------------------------------
      ct = 0
      do while (ct < max_iter)
         iter = iter + 1
         etak = 1/real(iter + 1, dp)**2
         call inv(fx_val, Bm, d, ierr)
         y = fx_val
         call f(xt, z, x + d, fx_val, del_z, ierr, Jl, jac)
         converged = (norm2(fx_val) < (gamma*norm2(y) - rhok*norm2(d)**2))
         if (ierr) converged = .false.
         ldk = 1
         lk = 0
         do while (.not. converged)
            lk = lk + 1
            ldk = dta**lk
            call f(xt, z, x + ldk*d, fx_val, del_z, ierr, Jl, jac)
            converged = (norm2(fx_val) < ((1.0_dp + etak)*phi - sgm*norm2(dta**lk*d)**2))
            if (ierr) converged = .false.
            if (lk > 20 .and. norm2(y) > 1d-8) then
               ierr = .true.
               deallocate (Hm, Bm, dx, d, y, grad, grad_new, tempj, x)
               return
            end if
         end do
         if (lk > 5) then
            ct = ct + 1
         end if
         dx = ldk*d

         y = fx_val - y
         tk = ((1.0_dp + etak)*phi + 1.0_dp)*norm2(fx_val)/(norm2(fx_val) + 1.0_dp)
         phi = 0.7*tk + 0.3*norm2(fx_val)
         ! Update x and evaluate f(x).
         x = x + dx

         Bm = Bm + outer_product((y - matmul(Bm, dx)), dx)/dot_product(dx, dx)

         tf = MPI_Wtime()
         ! if (mod(ct, 10) == 0 .and. lk > 5) then
         !    ! print *, "restart"
         !    call calculate_fqms(xt, z, x, Bm, Hm, del_z, ierr, Jl, jac)
         ! end if
         ! print *, "Iteration:", iter, "time:", tf - ti, "||f(x)|| =", norm2(fx_val), "|d|", norm2(dx), "lk", lk
         if (norm2(fx_val) < tol) then
            x_new = xt
            x_new(2:) = x(1:n/2)
            deallocate (Hm, Bm, dx, d, y, grad, grad_new, tempj, x)
            ierr = .false.
            return
         end if
      end do
      ierr = .true.
      deallocate (Hm, Bm, dx, d, y, grad, grad_new, tempj, x)
   end subroutine quasi_newtons
   subroutine calculate_fqms(xt, z, xi, fq, fqinv, del_z, ierr, Jl, jac)
      implicit none
      real(dp), intent(in) :: xt(:), xi(:), del_z(:)
      complex(dp), intent(in) :: z(:), jac(:, :)
      real(dp), intent(out) :: fq(:, :), Jl(:), fqinv(:, :)
      logical, intent(out) :: ierr
      real(dp), allocatable :: x_new(:)
      complex(dp), allocatable :: z_new(:), jprime(:, :)
      allocate (x_new(size(xt)), z_new(size(z)), jprime(size(jac, 1), size(jac, 2)))
      x_new = xt
      x_new(2:) = x_new(2:) + xi(1:size(z))
      call flow(x_new, z_new, jprime, ierr)
      fq(1:size(z), 1:size(z)) = real(jprime)
      fq(size(z) + 1:2*size(z), 1:size(z)) = aimag(jprime)
      fq(1:size(z), size(z) + 1:2*size(z)) = -aimag(jac)
      fq(size(z) + 1:2*size(z), size(z) + 1:2*size(z)) = real(jac)
      fqinv = fq
   end subroutine calculate_fqms

   subroutine calculate_fqvs(xt, z, xi, fq, del_z, ierr, Jl, jac)
      implicit none
      real(dp), intent(in) :: xt(:), xi(:), del_z(:)
      complex(dp), intent(in) :: z(:), jac(:, :)
      real(dp), intent(out) :: fq(:), Jl(:)
      logical, intent(out) :: ierr
      complex(dp), allocatable :: Jlc(:), z_new(:), z_new2(:)
      real(dp), allocatable :: x_new(:)
      ALLOCATE (Jlc(size(z)), z_new(size(z)), z_new2(size(z)), x_new(size(xt)))
      Jlc = matmul(jac, cmplx(0, 1)*xi(size(z) + 1:))
      call complex_to_real(Jlc, Jl)
      call real_to_complex(del_z, z_new)
      z_new = z + z_new - Jlc
      x_new = xt
      x_new(2:) = x_new(2:) + xi(1:size(z))
      call flowz(x_new, z_new2, ierr)
      z_new = z_new2 - z_new
      fq(1:size(z)) = real(z_new)
      fq(size(z) + 1:) = aimag(z_new)
      deallocate (Jlc, z_new, z_new2)
   end subroutine calculate_fqvs

end module quasi_newton_module
