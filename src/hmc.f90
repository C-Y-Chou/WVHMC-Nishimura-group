
module hmc
   use solve_flow
   use param_mod
   use utils
   use model
   use quasi_newton_module
   use mpi
   implicit none

contains
   subroutine rattle(state_x, state_z, step_size, num_steps, &
                     final_x, final_z, initial_hamiltonian, final_hamiltonian, jaci, jacf)
      implicit none

      !-------------------------------------------------------------------
      ! Input parameters
      !-------------------------------------------------------------------
      real(dp), intent(in)      :: state_x(:)
      complex(dp), intent(in)   :: state_z(:)
      real(dp), intent(in)      :: step_size
      integer, intent(in)       :: num_steps
      complex(dp), intent(in)   :: jaci(:, :)
      complex(dp), intent(out)   :: jacf(:, :)

      !-------------------------------------------------------------------
      ! Output parameters
      !-------------------------------------------------------------------
      real(dp), allocatable, intent(out)   :: final_x(:)
      complex(dp), allocatable, intent(out):: final_z(:)
      real(dp), intent(out)                :: initial_hamiltonian
      real(dp), intent(out)                :: final_hamiltonian

      !-------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------
      integer                :: step, state_size, wi
      real(dp)               :: integration_step_size, ti, tf
      logical                :: method_converged, has_error, reversed

      ! Momentum and force
      real(dp), allocatable  :: momentum(:), momentumuv(:), momentumu(:), momentumv(:), uv(:)
      real(dp), allocatable  :: dV(:), del_z(:)
      complex(dp), allocatable :: ds_val(:), E0(:)
      real(dp), allocatable  :: E0_real(:), E0_perp(:), E0_parl(:)

      ! Temporary states
      real(dp), allocatable  :: temp_x(:), z_real(:), z_new_real(:)
      complex(dp), allocatable :: temp_z(:), temp_jac(:, :)

      ! RATTLE-specific variables
      real(dp), allocatable  :: lambda(:)
      real(dp), allocatable  :: Jl(:), xiq(:)
      ! real(dp), parameter :: w_rattle(7) = (/ &
      !                        0.102799849391985_dp, &
      !                        -0.196061023297549E1_dp, &
      !                        0.193813913762276E1_dp, &
      !                        -0.158240635368243_dp, &
      !                        -0.144485223686048E1_dp, &
      !                        0.253693336566229_dp, &
      !                        0.914844246229740_dp/)
      ! real(dp), parameter :: w_rattle(3) = (/ &
      !                        -0.117767998417887E1_dp, &
      !                        0.235573213359357_dp, &
      !                        0.784513610477560_dp/)
      real(dp), parameter :: w_rattle(1) = (/1.0_dp/(2.0_dp - 2.0_dp**(1.0_dp/3.0_dp))/)

      !-------------------------------------------------------------------
      ! Initialization
      !-------------------------------------------------------------------
      has_error = .false.
      method_converged = .false.
      reversed = .false.
      state_size = size(state_z)

      allocate (final_x(size(state_x)), final_z(state_size))
      allocate (momentum(2*state_size))
      allocate (momentumuv(2*state_size), momentumu(2*state_size), momentumv(2*state_size))
      allocate (dV(2*state_size), del_z(2*state_size))
      allocate (ds_val(state_size), E0(state_size))
      allocate (E0_real(2*state_size), E0_perp(2*state_size), E0_parl(2*state_size))
      allocate (temp_x(size(state_x)), temp_z(state_size))
      allocate (z_real(2*state_size), z_new_real(2*state_size))
      allocate (lambda(2*state_size), Jl(state_size*2), xiq(1 + state_size))
      allocate (temp_jac(size(jaci, 1), size(jaci, 2)))
      allocate (uv(2*state_size))

      final_x = state_x
      final_z = state_z
      temp_jac = jaci

      ! Initial momentum and Hamiltonian
      call grand(momentum)
      if (istest) momentum = testmom
      if (wv) then
         call ds(state_z, ds_val)
         E0 = conjg(ds_val)
         call complex_to_real(E0, E0_real)
         call decompose2(state_x, E0_real, uv, E0_parl, E0_perp, temp_jac, has_error)
         call decompose2(state_x, momentum, momentumuv, momentumu, momentumv, temp_jac, has_error)
         momentum = momentumu + dot_product(momentumv, E0_perp)/norm2(E0_perp)**2*E0_perp
      else
         call decompose2(state_x, momentum, momentumuv, momentumu, momentumv, temp_jac, has_error)
         momentum = momentumu
      end if
      call calculate_hamiltonian(state_x, state_z, momentum, initial_hamiltonian)

      ! RATTLE loop
      do step = 1, num_steps
         if (eo) then
         do wi = 0, size(w_rattle) - 1
            temp_x = final_x
            temp_z = final_z
            integration_step_size = w_rattle(size(w_rattle) - wi)*step_size/real(num_steps, dp)
            call rattle_step(temp_x, temp_z, integration_step_size, final_x, final_z, &
                             temp_jac, jacf, momentum, method_converged)
            if (.not. method_converged) then
               final_hamiltonian = 0
               jacf = jaci
               call deallocate_all()
               return
            end if
            temp_jac = jacf
         end do
         temp_x = final_x
         temp_z = final_z
         integration_step_size = (1.0_dp - 2*(sum(w_rattle)))*step_size/real(num_steps, dp)
         call rattle_step(temp_x, temp_z, integration_step_size, final_x, final_z, &
                          temp_jac, jacf, momentum, method_converged)
         if (.not. method_converged) then
            final_hamiltonian = 0
            jacf = jaci
            call deallocate_all()
            return
         end if

         temp_jac = jacf

         do wi = 1, size(w_rattle)
            temp_x = final_x
            temp_z = final_z
            integration_step_size = w_rattle(wi)*step_size/real(num_steps, dp)
            call rattle_step(temp_x, temp_z, integration_step_size, final_x, final_z, &
                             temp_jac, jacf, momentum, method_converged)
            if (.not. method_converged) then
               final_hamiltonian = 0
               jacf = jaci
               call deallocate_all()
               return
            end if

            temp_jac = jacf
         end do
         else
         integration_step_size = step_size/real(num_steps, dp)
         temp_x = final_x
         temp_z = final_z
         call rattle_step(temp_x, temp_z, integration_step_size, final_x, final_z, temp_jac, jacf, momentum, method_converged)
         temp_jac = jacf
         end if

      end do
      if (final_x(2) == temp_x(2)) then
         final_hamiltonian = 0
         jacf = jaci
         call deallocate_all()
         return
      end if
      call ds(final_z, ds_val)
      E0 = conjg(ds_val)
      call complex_to_real(E0, E0_real)
      if (wv) then
         call decompose2(final_x, E0_real, uv, E0_parl, E0_perp, temp_jac, has_error)
         call decompose2(final_x, momentum, momentumuv, momentumu, momentumv, temp_jac, has_error)
         if (has_error) then
            method_converged = .false.
            call deallocate_all()
            return
         end if
         momentum = momentumu + dot_product(momentumv, E0_perp)/norm2(E0_perp)**2*E0_perp
      else
         call decompose2(final_x, momentum, momentumuv, momentumu, momentumv, temp_jac, has_error)
         if (has_error) then
            method_converged = .false.
            call deallocate_all()
            return
         end if
         momentum = momentumu
      end if
      !if (method_converged) then
      !momentumu = -momentum
      !call crv(final_x, final_z, integration_step_size, num_steps, &
      !temp_x, temp_z, temp_jac, momentumu)
      !if (norm2(abs(temp_z - state_z))>1d-9) then
      !print*,"denied by reversibility",norm2(abs(temp_z - state_z))
      !method_converged = .false.
      !endif
      !endif
      ! Final Hamiltonian
      if (method_converged) then
         call calculate_hamiltonian(final_x, final_z, momentum, final_hamiltonian)
      else
         final_hamiltonian = 0
         jacf = jaci
      end if

      call deallocate_all()
      return

   contains

      subroutine deallocate_all()
         if (allocated(momentum)) deallocate (momentum)
         if (allocated(momentumuv)) deallocate (momentumuv)
         if (allocated(momentumu)) deallocate (momentumu)
         if (allocated(momentumv)) deallocate (momentumv)
         if (allocated(dV)) deallocate (dV)
         if (allocated(del_z)) deallocate (del_z)
         if (allocated(ds_val)) deallocate (ds_val)
         if (allocated(E0)) deallocate (E0)
         if (allocated(E0_real)) deallocate (E0_real)
         if (allocated(E0_perp)) deallocate (E0_perp)
         if (allocated(E0_parl)) deallocate (E0_parl)
         if (allocated(temp_x)) deallocate (temp_x)
         if (allocated(temp_z)) deallocate (temp_z)
         if (allocated(z_real)) deallocate (z_real)
         if (allocated(z_new_real)) deallocate (z_new_real)
         if (allocated(lambda)) deallocate (lambda)
         if (allocated(Jl)) deallocate (Jl)
         if (allocated(xiq)) deallocate (xiq)
         if (allocated(temp_jac)) deallocate (temp_jac)
         if (allocated(uv)) deallocate (uv)
      end subroutine deallocate_all

   end subroutine rattle
   subroutine rattle_step(state_x, state_z, step_size, &
                          final_x, final_z, jaci, jacf, momentum, method_converged)
      implicit none

      !-------------------------------------------------------------------
      ! Input parameters
      !-------------------------------------------------------------------
      real(dp), intent(in)      :: state_x(:)
      complex(dp), intent(in)   :: state_z(:)
      real(dp), intent(in)      :: step_size
      real(dp), intent(inout)    :: momentum(:)
      complex(dp), intent(in)   :: jaci(:, :)
      complex(dp), intent(out)   :: jacf(:, :)

      !-------------------------------------------------------------------
      ! Output parameters
      !-------------------------------------------------------------------
      real(dp), allocatable, intent(out)   :: final_x(:)
      complex(dp), allocatable, intent(out):: final_z(:)

      !-------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------
      integer                :: step, state_size
      real(dp)               :: integration_step_size, ti, tf
      logical                :: method_converged, has_error, reversed

      ! Momentum and force
      real(dp), allocatable  :: momentumuv(:), momentumu(:), momentumv(:), uv(:)
      real(dp), allocatable  :: dV(:), del_z(:)
      complex(dp), allocatable :: ds_val(:), E0(:)
      real(dp), allocatable  :: E0_real(:), E0_perp(:), E0_parl(:)

      ! Temporary states
      real(dp), allocatable  :: temp_x(:), z_real(:), z_new_real(:)
      complex(dp), allocatable :: temp_z(:), temp_jac(:, :)

      ! RATTLE-specific variables
      real(dp), allocatable  :: lambda(:)
      real(dp), allocatable  :: Jl(:), xiq(:)

      !-------------------------------------------------------------------
      ! Initialization
      !-------------------------------------------------------------------
      has_error = .false.
      method_converged = .true.
      reversed = .false.
      state_size = size(state_z)
      integration_step_size = step_size

      allocate (final_x(size(state_x)), final_z(state_size))
      allocate (momentumuv(2*state_size), momentumu(2*state_size), momentumv(2*state_size))
      allocate (dV(2*state_size), del_z(2*state_size))
      allocate (ds_val(state_size), E0(state_size))
      allocate (E0_real(2*state_size), E0_perp(2*state_size), E0_parl(2*state_size))
      allocate (temp_x(size(state_x)), temp_z(state_size))
      allocate (z_real(2*state_size), z_new_real(2*state_size))
      allocate (lambda(2*state_size), Jl(state_size*2), xiq(1 + state_size))
      allocate (temp_jac(size(jaci, 1), size(jaci, 2)))
      allocate (uv(2*state_size))

      final_x = state_x
      final_z = state_z
      temp_jac = jaci

      call ds(state_z, ds_val)
      E0 = conjg(ds_val)
      call complex_to_real(E0, E0_real)
      if (wv) then
         call decompose2(state_x, E0_real, uv, E0_parl, E0_perp, temp_jac, has_error)
         call calculate_dV(state_x, E0_real, E0_perp, dV, has_error)
         del_z = integration_step_size*momentum - integration_step_size**2*dV
      else
         call calculate_dV(state_x, E0_real, E0_perp, dV, has_error)
         del_z = integration_step_size*momentum - integration_step_size**2*dV
      end if

      ! RATTLE loop
      do step = 1, 1
         temp_x = final_x
         temp_z = final_z
         lambda = 0
         ! call simp_newton(cttol, 100, temp_x, temp_z, del_z, has_error, Jl, final_z, temp_jac)
         ! if (has_error) then
         !    call quasi_newtonv2(calculate_fqv2, cttol, 1000, temp_x, temp_z, del_z, has_error, Jl, final_z, temp_jac)
         ! end if
         if (wv) then
            call simp_newtonwv(cttol, 100, temp_x, temp_z, del_z, has_error, Jl, final_x, temp_jac)
            ! has_error = .true.
            if (has_error) then
               ! call quasi_newtons(calculate_fqvs, cttol, 100, temp_x, temp_z, del_z, has_error, Jl, final_x, temp_jac)
               call quasi_newtonwv(calculate_fqwv, cttol, 100, temp_x, temp_z, del_z, has_error, Jl, final_x, temp_jac)
            end if

            if (final_x(1) < t0 - d0 .or. final_x(1) > t1 + d0) then
               ! print *, "denied", final_x(1)
               final_x = temp_x
               final_z = temp_z
               momentum = -momentum
               call deallocate_all()
               return
            end if
         else
            call simp_newton(cttol, 100, temp_x, temp_z, del_z, has_error, Jl, final_x, temp_jac)
            ! has_error = .true.
            if (has_error) then
               ! call quasi_newtons(calculate_fqvs, cttol, 100, temp_x, temp_z, del_z, has_error, Jl, final_x, temp_jac)
               call quasi_newtonv2(calculate_fqv2, cttol, 100, temp_x, temp_z, del_z, has_error, Jl, final_x, temp_jac)
            end if
         end if
         method_converged = .not. has_error
         if (has_error) then
            method_converged = .false.
            call deallocate_all()
            return
         end if
         call flow(final_x, final_z, temp_jac, has_error)
         if (has_error) then
            method_converged = .false.
            call deallocate_all()
            return
         end if
         call complex_to_real((final_z - temp_z)/integration_step_size, momentum)
         call ds(final_z, ds_val)
         E0 = conjg(ds_val)
         call complex_to_real(E0, E0_real)
         if (wv) then
            call decompose2(final_x, E0_real, uv, E0_parl, E0_perp, temp_jac, has_error)
            call calculate_dV(final_x, E0_real, E0_perp, dV, has_error)
            momentum = momentum - integration_step_size*dV
         else
            call calculate_dV(final_x, E0_real, E0_perp, dV, has_error)
            momentum = momentum - integration_step_size*dV
         end if
      end do
      jacf = temp_jac
      method_converged = .true.
      call deallocate_all()
   contains

      subroutine deallocate_all()
         if (allocated(momentumuv)) deallocate (momentumuv)
         if (allocated(momentumu)) deallocate (momentumu)
         if (allocated(momentumv)) deallocate (momentumv)
         if (allocated(dV)) deallocate (dV)
         if (allocated(del_z)) deallocate (del_z)
         if (allocated(ds_val)) deallocate (ds_val)
         if (allocated(E0)) deallocate (E0)
         if (allocated(E0_real)) deallocate (E0_real)
         if (allocated(E0_perp)) deallocate (E0_perp)
         if (allocated(E0_parl)) deallocate (E0_parl)
         if (allocated(temp_x)) deallocate (temp_x)
         if (allocated(temp_z)) deallocate (temp_z)
         if (allocated(z_real)) deallocate (z_real)
         if (allocated(z_new_real)) deallocate (z_new_real)
         if (allocated(lambda)) deallocate (lambda)
         if (allocated(Jl)) deallocate (Jl)
         if (allocated(xiq)) deallocate (xiq)
         if (allocated(temp_jac)) deallocate (temp_jac)
         if (allocated(uv)) deallocate (uv)
      end subroutine deallocate_all

   end subroutine rattle_step

   subroutine rattle2(state_x, state_z, step_size, num_steps, &
                      final_x, final_z, initial_hamiltonian, final_hamiltonian, jaci)
      implicit none

      !-------------------------------------------------------------------
      ! Input parameters
      !-------------------------------------------------------------------
      real(dp), intent(in)      :: state_x(:)
      complex(dp), intent(in)   :: state_z(:)
      real(dp), intent(in)      :: step_size
      integer, intent(in)       :: num_steps
      complex(dp), intent(inout)   :: jaci(:, :)

      !-------------------------------------------------------------------
      ! Output parameters
      !-------------------------------------------------------------------
      real(dp), allocatable, intent(out)   :: final_x(:)
      complex(dp), allocatable, intent(out):: final_z(:)
      real(dp), intent(out)                :: initial_hamiltonian
      real(dp), intent(out)                :: final_hamiltonian

      !-------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------
      integer                :: step, state_size
      real(dp)               :: integration_step_size
      logical                :: method_converged, has_error

      ! Momentum and force
      real(dp), allocatable  :: momentum(:), momentumuv(:), momentumu(:), momentumv(:)
      real(dp), allocatable  :: dV(:), del_z(:)
      complex(dp), allocatable :: ds_val(:), E0(:)
      real(dp), allocatable  :: E0_real(:), E0_perp(:)

      ! Temporary states
      real(dp), allocatable  :: temp_x(:), z_real(:), z_new_real(:)
      complex(dp), allocatable :: temp_z(:), temp_jac(:, :)

      ! RATTLE-specific variables
      real(dp), allocatable  :: lambda(:)
      real(dp), allocatable  :: Jl(:), xiq(:)

      !-------------------------------------------------------------------
      ! Initialization
      !-------------------------------------------------------------------
      has_error = .false.
      method_converged = .false.
      state_size = size(state_z)
      integration_step_size = step_size/real(num_steps, dp)

      allocate (final_x(size(state_x)), final_z(state_size))
      allocate (momentum(2*state_size))
      allocate (momentumuv(2*state_size), momentumu(2*state_size), momentumv(2*state_size))
      allocate (dV(2*state_size), del_z(2*state_size))
      allocate (ds_val(state_size), E0(state_size))
      allocate (E0_real(2*state_size), E0_perp(2*state_size))
      allocate (temp_x(size(state_x)), temp_z(state_size))
      allocate (z_real(2*state_size), z_new_real(2*state_size))
      allocate (lambda(2*state_size), Jl(state_size*2), xiq(1 + state_size))
      allocate (temp_jac(size(jaci, 1), size(jaci, 2)))

      final_x = state_x
      final_z = state_z
      temp_jac = jaci

      ! Initial momentum and Hamiltonian
      momentum = 0
      call calculate_hamiltonian(final_x, final_z, momentum, initial_hamiltonian)
      ! Initial force
      call ds(state_z, ds_val)
      E0 = conjg(ds_val)
      call complex_to_real(E0, E0_real)
      call calculate_dV(state_x, E0_real, E0_perp, dV, has_error)
      del_z = integration_step_size*momentum - integration_step_size**2*dV

      ! RATTLE loop
      do step = 1, num_steps
         temp_x = final_x
         temp_z = final_z
         lambda = 0

         ! call simp_newton(cttol, 100, temp_x, temp_z, del_z, has_error, Jl, final_z, temp_jac)
         ! if (has_error) then
         !    call quasi_newtonv2(calculate_fqv2, cttol, 1000, temp_x, temp_z, del_z, has_error, Jl, final_z, temp_jac)
         ! end if
         call quasi_newtonv2(calculate_fqv2, cttol, 1000, temp_x, temp_z, del_z, has_error, Jl, final_x, temp_jac)
         if (.not. has_error) method_converged = .true.
         if (has_error) then
            method_converged = .false.
            return
         end if
         lambda = Jl
         call flow(final_x, final_z, temp_jac, has_error)
         if (has_error) then
            method_converged = .false.
            return
         end if
      end do
      ! Final Hamiltonian
      if (method_converged) then
         momentum = 0
         call calculate_hamiltonian(final_x, final_z, momentum, final_hamiltonian)
      else
         final_hamiltonian = 0
      end if

      call deallocate_all()
      return

   contains

      subroutine deallocate_all()
         if (allocated(momentum)) deallocate (momentum)
         if (allocated(momentumuv)) deallocate (momentumuv)
         if (allocated(momentumu)) deallocate (momentumu)
         if (allocated(momentumv)) deallocate (momentumv)
         if (allocated(dV)) deallocate (dV)
         if (allocated(del_z)) deallocate (del_z)
         if (allocated(ds_val)) deallocate (ds_val)
         if (allocated(E0)) deallocate (E0)
         if (allocated(E0_real)) deallocate (E0_real)
         if (allocated(E0_perp)) deallocate (E0_perp)
         if (allocated(temp_x)) deallocate (temp_x)
         if (allocated(temp_z)) deallocate (temp_z)
         if (allocated(z_real)) deallocate (z_real)
         if (allocated(z_new_real)) deallocate (z_new_real)
         if (allocated(lambda)) deallocate (lambda)
         if (allocated(Jl)) deallocate (Jl)
         if (allocated(xiq)) deallocate (xiq)
         if (allocated(temp_jac)) deallocate (temp_jac)
      end subroutine deallocate_all

   end subroutine rattle2

   subroutine calculate_dV(x, E0_real, E0_perp, dV, error)
      implicit none

! Inputs
      real(dp), intent(in) :: x(:), E0_real(:), E0_perp(:)

! Outputs
      real(dp), dimension(:), intent(out) :: dV     ! Derivative of the potential (2n length)
      logical, intent(out) :: error  ! Error flag

! Locals
      integer :: n
      real(dp) :: dW

! Basic checks
      n = size(x) - 1  ! z is length n
      error = .false.
      if (size(dV) /= 2*n) then
         write (*, *) "Error(calculate_dV): dV must have length 2*size(z)."
         error = .true.
         return
      end if

! Derivative of weight function W(t)
      call W_t_derivative(x(1), dW)

! Calculate dV = E0_real + dW * E0_perp / sum(E0_perp^2)
      if (wv) then
         dV = (E0_real + dW*E0_perp/norm2(E0_perp)**2)/2.0_dp
      else
         dV = E0_real/2.0_dp
      end if

   end subroutine calculate_dV
!  subroutine decompose(xt, b, x, au, av, max_iter, tol, ierr)
!        implicit none
!        real(dp), intent(in)    :: xt(:)    ! Possibly geometry or other data
!        real(dp), intent(in)    :: b(:)    ! Right-hand side or initial residual
!        real(dp), intent(inout) :: x(:)   ! Approximate solution
!        real(dp), intent(out)   :: au(:), av(:)  ! Optional outputs from flowv
!        integer, intent(in)     :: max_iter       ! Maximum iterations
!        real(dp), intent(in)    :: tol            ! Convergence tolerance
!        logical, intent(out)    :: ierr          ! Error flag
!        real(dp):: res

!        integer :: i, j, k, outer, j_iter, m, n, iter
!        real(dp) :: beta_d, resid, temp, denom

!        ! Allocate arrays for the Krylov basis (V), the upper Hessenberg matrix (H),
!        ! the right-hand side of the least–squares problem (g), the rotation parameters,
!        ! and some work vectors.
!        real(dp), allocatable :: V(:, :), H(:, :), g(:), y(:), w(:)
!        real(dp), allocatable :: cs(:), sn(:)
!        real(dp), allocatable :: aux(:), a_x(:)
!        n = size(b)
!        m = 100
!        ierr = .true.
!        allocate (aux(n), a_x(n))
!        allocate (V(n, m + 1))
!        allocate (H(m + 1, m))
!        allocate (g(m + 1))
!        allocate (cs(m))
!        allocate (sn(m))
!        allocate (w(n))
!        allocate (y(m))

!        !--------------------------------------------------------------------
!        ! Compute the initial residual r = b - A*x and its norm beta_d.
!        !--------------------------------------------------------------------
!        x = 0
!        w = b
!        beta_d = sqrt(dot_product(w, w))
!        res = beta_d
!        if (beta_d < tol) then
!           iter = 0
!           deallocate (V, H, g, cs, sn, w, y, aux, a_x)
!           return
!        end if
!        iter = 0

!        !--------------------------------------------------------------------
!        ! Outer (restart) loop.
!        !--------------------------------------------------------------------
!        do outer = 1, max_iter
!           ! Set the first Krylov vector.
!           V(:, 1) = w/beta_d

!           ! Initialize the right–hand side of the least–squares problem.
!           g = 0.0d0
!           g(1) = beta_d

!           ! Clear H and the rotation parameters.
!           H = 0.0d0
!           cs = 0.0d0
!           sn = 0.0d0

!           !-----------------------------------------------------------------
!           ! Build an m–dimensional Krylov subspace.
!           !-----------------------------------------------------------------
!           do j = 1, m
!              j_iter = j         ! record the current subspace dimension
!              iter = iter + 1

!              ! w = A * V(:,j)
!              call flowv(xt, V(:, j), aux, w, ierr)
!              if (ierr) then
!                 print *, "error 1"
!                 return
!              end if
!              w = w + aux

!              ! Modified Gram–Schmidt: orthogonalize w against V(:,1:j)
!              do i = 1, j
!                 H(i, j) = dot_product(w, V(:, i))
!                 w = w - H(i, j)*V(:, i)
!              end do
!              H(j + 1, j) = sqrt(dot_product(w, w))

!              if (H(j + 1, j) == 0.0d0) then
!                 ! Happy breakdown: exact invariant subspace found.
!                 exit
!              else
!                 V(:, j + 1) = w/H(j + 1, j)
!              end if

!              !-----------------------------------------------------------------
!              ! Apply all previously computed Givens rotations to the new column H(:,j)
!              !-----------------------------------------------------------------
!              do i = 1, j - 1
!                 temp = cs(i)*H(i, j) + sn(i)*H(i + 1, j)
!                 H(i + 1, j) = -sn(i)*H(i, j) + cs(i)*H(i + 1, j)
!                 H(i, j) = temp
!              end do

!              !-----------------------------------------------------------------
!              ! Compute the j-th Givens rotation to eliminate H(j+1,j)
!              !-----------------------------------------------------------------
!              denom = sqrt(H(j, j)**2 + H(j + 1, j)**2)
!              if (denom == 0.0d0) then
!                 cs(j) = 1.0d0
!                 sn(j) = 0.0d0
!              else
!                 cs(j) = H(j, j)/denom
!                 sn(j) = H(j + 1, j)/denom
!              end if

!              ! Apply the rotation to H(j,j) and update the right–hand side g.
!              H(j, j) = cs(j)*H(j, j) + sn(j)*H(j + 1, j)
!              g(j + 1) = -sn(j)*g(j)
!              g(j) = cs(j)*g(j)

!              resid = abs(g(j + 1))
!              ! Check for convergence.
!              if (resid < tol) then
!                 ! Solve the upper triangular system H(1:j,1:j)*y = g(1:j)
!                 do k = j, 1, -1
!                    g(k) = g(k)/H(k, k)
!                    do i = 1, k - 1
!                       g(i) = g(i) - H(i, k)*g(k)
!                    end do
!                 end do
!                 y(1:j) = g(1:j)
!                 ! Update the approximate solution: x = x + V(:,1:j)*y
!                 x = x + matmul(V(:, 1:j), y(1:j))
!                 ierr = .false.
!                 call flowv(xt, x, au, av, ierr)
!                 if (ierr) then
!                    print *, "error 2"
!                    return
!                 end if
!                 res = resid
!                 deallocate (V, H, g, cs, sn, w, y, aux, a_x)
!                 return
!              end if
!           end do  ! End of inner loop.

!           !-----------------------------------------------------------------
!           ! If we reached here, we did not converge within the inner loop.
!           ! Solve the least–squares problem for the current subspace dimension j_iter.
!           !-----------------------------------------------------------------
!           do k = j_iter, 1, -1
!              g(k) = g(k)/H(k, k)
!              do i = 1, k - 1
!                 g(i) = g(i) - H(i, k)*g(k)
!              end do
!           end do
!           y(1:j_iter) = g(1:j_iter)
!           x = x + matmul(V(:, 1:j_iter), y(1:j_iter))

!           ! Compute the new residual.
!           call flowv(xt, x, aux, w, ierr)
!           if (ierr) then
!              print *, "error 3"
!              return
!           end if
!           w = b - (w + aux)
!           beta_d = sqrt(dot_product(w, w))
!           res = beta_d
!           if (beta_d < tol) then
!              ierr = .false.
!              call flowv(xt, x, au, av, ierr)
!              call flowv(xt, x, aux, w, ierr)
!              if (ierr) then
!                 print *, "error 4"
!                 return
!              end if
!              deallocate (V, H, g, cs, sn, w, y, aux, a_x)
!              return
!           end if
!           ! Prepare for the next restart by setting beta_d and continuing.
!        end do

!        ! If we reach here, max_iter was exceeded.
!        deallocate (V, H, g, cs, sn, w, y, aux, a_x)
!        return
!     end subroutine decompose
   subroutine decompose2(xt, b, x, au, av, jac, ierr)
      implicit none
      real(dp), intent(in)    :: xt(:)
      real(dp), intent(in)    :: b(:)
      real(dp), intent(inout) :: x(:)
      real(dp), intent(out)   :: au(:), av(:)
      complex(dp), intent(in) :: jac(:, :)
      logical, intent(out)    :: ierr

      ! Local variables
      real(dp), allocatable :: jacr(:, :)
      integer :: n, info
      integer, allocatable :: ipiv(:)
      allocate (jacr(2*size(jac, 1), 2*size(jac, 2)))
      ! Convert complex matrix to real matrix
      call map_to_real_mat(jac, jacr)
      ierr = .false.

      ! Size of system
      n = size(b)
      allocate (ipiv(n))
      x = b
      ! Solve the system jacr * x = b using LAPACK (DGETRF + DGETRS)
      call dgesv(n, 1, jacr, n, ipiv, x, n, info)
      if (info /= 0) then
         ierr = .true.
         return
      end if

      ! Assign outputs
      au = x
      call real_vec(au)
      av = x - au
      call map_to_real_mat(jac, jacr)
      au = matmul(jacr, au)
      av = b - au
      deallocate (jacr)
   end subroutine decompose2
   subroutine calculate_hamiltonian(x, z, p, h)
      implicit none

! Inputs
      real(dp), intent(in) :: x(:)
      complex(dp), intent(in) :: z(:)
      real(dp), intent(in) :: p(:)

! Outputs
      real(dp), intent(out) :: h

      complex(dp)  :: s
      real(dp)     :: w

      if (2*size(z) /= size(p)) then
         write (*, *) "Warning: z and p differ in length in calculate_hamiltonian."
      end if
      if (wv) then
         call calculate_action(z, s)
         call W_function(x(1), w)
         h = 0.5_dp*norm2(p)**2 + real(s, dp) + w
      else
         call calculate_action(z, s)
         h = 0.5_dp*norm2(p)**2 + real(s, dp)
      end if
   end subroutine calculate_hamiltonian

   subroutine simp_newton(tol, max_iter, xt, z, del_z, ierr, Jl, x_new, jac)
      implicit none

      ! Arguments
      integer, intent(in)          :: max_iter
      real(dp), intent(in)         :: tol
      logical, intent(out)         :: ierr
      real(dp), intent(in) :: xt(:), del_z(:)
      complex(dp), intent(in) :: z(:)
      real(dp), intent(out) :: x_new(:)
      real(dp), intent(out) :: Jl(:)
      complex(dp), intent(in) :: jac(:, :)

      ! locals
      real(dp), allocatable :: B(:)
      real(dp), allocatable :: xtu(:), u(:), dxi(:), au(:), av(:)
      real(dp)::ti, tf, B_old
      complex(dp), allocatable :: dummy(:), ld(:)
      integer :: n, iter, ct
      complex(dp) :: z_new(size(z))

      n = size(z)

      allocate (B(2*n), xtu(1 + n), ld(n), u(n), dxi(2*n), au(2*n), av(2*n), dummy(n))
      xtu = xt
      B = del_z
      B_old = norm2(B)
      u = 0
      ld = 0
      ! print *, norm2(B)
      iter = 0
      ct = 0
      ierr = .true.
      call cpu_time(ti)
      do while (iter < max_iter)
         iter = iter + 1
         call decompose2(xt, B, dxi, au, av, jac, ierr)
         call real_to_complex(dxi, dummy)
         u = u + real(dummy, dp)
         call real_to_complex(av, dummy)
         ld = ld + dummy

         xtu(2:) = xt(2:) + u
         call flowz(xtu, z_new, ierr)
         z_new = z - z_new - ld
         call complex_to_real(z_new, B)
         B = B + del_z
         call cpu_time(tf)
         ! print *, tf - ti, iter, norm2(B)
         if (norm2(B) < tol) then
            x_new = xtu
            return
         end if
         if (norm2(B) > B_old) then
            ct = ct + 1
         end if
         if (ct > 3) then
            ierr = .true.
            return
         end if
      end do
      ierr = .true.

   end subroutine simp_newton

   subroutine simp_newtonwv(tol, max_iter, xt, z, del_z, ierr, Jl, x_new, jac)
      implicit none

      ! Arguments
      integer, intent(in)          :: max_iter
      real(dp), intent(in)         :: tol
      logical, intent(out)         :: ierr
      real(dp), intent(in) :: xt(:), del_z(:)
      complex(dp), intent(in) :: z(:)
      real(dp), intent(out) :: x_new(:)
      real(dp), intent(out) :: Jl(:)
      complex(dp), intent(in) :: jac(:, :)

      ! locals
      real(dp), allocatable :: B(:)
      real(dp), allocatable :: xtu(:), u(:), dxi(:), au(:), av(:), E0_real(:), uv(:), E0_parl(:), E0_perp(:)
      real(dp)::ti, tf, B_old, cb
      complex(dp), allocatable :: dummy(:), ld(:), z_new(:)
      integer :: n, iter, ct

      n = size(z)

      allocate (B(2*n), xtu(1 + n), ld(n), u(n + 1), dxi(2*n), au(2*n), av(2*n), dummy(n))
      allocate (E0_real(2*n), uv(2*n), E0_parl(2*n), E0_perp(2*n), z_new(n))
      xtu = xt
      B = del_z
      B_old = norm2(B)
      u = 0
      ld = 0
      ! print *, norm2(B)
      iter = 0
      ct = 0
      ierr = .true.
      call ds(z, dummy)
      dummy=conjg(dummy)
      call complex_to_real(dummy, E0_real)
      call decompose2(xt, E0_real, uv, E0_parl, E0_perp, jac, ierr)
      call cpu_time(ti)
      do while (iter < max_iter)
         iter = iter + 1
         cb = sum(B*E0_perp)/sum(E0_perp*E0_perp)
         u(1) = u(1) + cb
         call decompose2(xt, B, dxi, au, av, jac, ierr)
         call real_to_complex(dxi - cb*uv, dummy)
         u(2:) = u(2:) + real(dummy, dp)
         call real_to_complex(av - cb*E0_perp, dummy)
         ld = ld + dummy
         xtu = xt + u
         call flowz(xtu, z_new, ierr)
         z_new = z - z_new - ld
         call complex_to_real(z_new, B)
         B = B + del_z
         call cpu_time(tf)
         ! print *, tf - ti, iter, norm2(B), xtu(1)
         if (norm2(B) < tol) then
            x_new = xtu
            return
         end if
         if (norm2(B) > B_old) then
            ct = ct + 1
         end if
         if (ct > 3) then
            ierr = .true.
            return
         end if
      end do
      ierr = .true.

   end subroutine simp_newtonwv
end module hmc
