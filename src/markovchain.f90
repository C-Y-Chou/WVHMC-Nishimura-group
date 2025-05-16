module markovchain
   use solve_flow
   use param_mod
   use model
   use hmc
   use mt95
   use utils
   use mpi
   implicit none

contains
   subroutine metropolis_step(x, z, j, total_step_size, num_steps, x_new, z_new, j_new, accept)
      ! Perform one step of the Metropolis procedure using the RATTLE integrator.
      implicit none

      ! Inputs
      real(dp), intent(in) :: x(:)                         ! Current state vector
      complex(dp), intent(in) :: z(:)                      ! Current state vector in the complex space
      complex(dp), dimension(:, :), intent(in) :: j        ! Current Jacobian matrix
      real(dp), intent(in) :: total_step_size              ! Total step size for RATTLE
      integer, intent(in) :: num_steps                         ! Number of RATTLE steps

      ! Outputs
      real(dp), allocatable, intent(out) :: x_new(:)       ! Updated state vector
      complex(dp), allocatable, intent(out) :: z_new(:)    ! Updated state vector in the complex space
      complex(dp), allocatable, dimension(:, :), intent(out) :: j_new  ! Updated Jacobian matrix
      logical, intent(out) :: accept                           ! Whether the proposed state was accepted

      ! Local variables
      real(dp) :: h_initial                                ! Initial Hamiltonian
      real(dp) :: h_final                                  ! Final Hamiltonian after RATTLE
      real(dp) :: accept_probability                       ! Metropolis acceptance probability
      real(dp) :: rand
      logical :: rattle_step_error                             ! Error flag for RATTLE

      ! Initialize outputs
      accept = .false.

      ! Allocate outputs
      allocate (x_new(size(x)))
      allocate (z_new(size(z)))
      allocate (j_new(size(j, 1), size(j, 2)))
      ! Perform the RATTLE integration to propose a new state
      call rattle(x, z, total_step_size, num_steps, x_new, z_new, h_initial, h_final, j, j_new)
      ! Compute Metropolis acceptance probability
      accept_probability = exp(-(h_final - h_initial))
      rand = grnd()

      ! Decide whether to accept the proposed state
      if (accept_probability >= 1.0_8 .or. rand <= accept_probability) then
         accept = .true.
      else
         accept = .false.
      end if
      if (h_final == 0) then
         accept = .false.
         print *, "not converged"
      end if
   end subroutine metropolis_step
   subroutine generate_markov_chain(chain_length, x_initial, total_step_size, num_steps, &
                                    x_history_file, z_history_file, phi_history_file)
      ! Generate a Markov chain and save the history of x, z, and j.
      implicit none

      ! Inputs
      integer, intent(in) :: chain_length
      real(dp), intent(in), allocatable :: x_initial(:)
      real(dp), intent(in) :: total_step_size
      integer, intent(in) :: num_steps
      character(len=*), intent(in) :: x_history_file, z_history_file, phi_history_file

      ! Local variables
      real(dp), allocatable :: x(:), x_new(:)
      complex(dp), allocatable :: z(:), z_new(:)
      complex(dp), allocatable :: j(:, :), j_new(:, :)
      real(dp) :: h_initial, h_proposed
      logical :: accept, flow_error_flag, project_error, error
      integer :: i, n, row, col, k, accepted_steps
      real(dp) :: acceptance_rate
      complex(dp) :: det_j, phi, s_i
      complex(dp), dimension(size(x_initial) - 1) :: ds_val, E0     ! Derivative and conjugate components
      real(dp), dimension(size(x_initial)*2 - 2) :: E0_real, E0_perp, uv, E0_parl
      real(dp), dimension(size(M_diag)) :: M_variance, z_real, M_mean
      real(dp), dimension(size(M_diag), n_warm) :: M_his
      complex(dp), dimension(size(M_diag)/2) :: M_temp
      ! >>> (1) 定義時間相關變數 <<<
      real(dp) :: time_start, time_current, time_elapsed, time_per_iter, time_left
      real(dp) :: progress

      ! Allocate and initialize variables
      n = size(x_initial)
      allocate (x(n), x_new(n))
      allocate (z(n - 1), z_new(n - 1))
      allocate (j(n - 1, n - 1), j_new(n - 1, n - 1))
      accepted_steps = 0

      ! Open files to save history
      open (unit=10, file=x_history_file, status='replace', access='stream', form='unformatted')
      open (unit=20, file=z_history_file, status='replace', access='stream', form='unformatted')
      open (unit=30, file=phi_history_file, status='replace', access='stream', form='unformatted')

      ! >>> (2) 取得初始時間 <<<
      time_start = MPI_Wtime()

      ! Perform the initial flow calculation
      x = x_initial
      call flow(x, z, j, flow_error_flag)
      if (flow_error_flag) then
         print *, "Error: Initial flow calculation failed."
         stop
      end if
      h_proposed = h_initial*2
      i = 0
      do while ((h_proposed < h_initial .or. i < 5) .and. i < n_warm)
         i = i + 1
         acceptance_rate = total_step_size/num_steps
         call rattle2(x, z, acceptance_rate, 1, x_new, z_new, h_initial, h_proposed, j)
         x = x_new
         if (h_proposed == 0) acceptance_rate = acceptance_rate/2
         call flow(x, z, j, flow_error_flag)
         print *, i
         print *, h_proposed
      end do
      call save_initial_condition(initial_file, x)

      if (wv) then
         call ds(z, ds_val)
         E0 = conjg(ds_val)
         call complex_to_real(E0, E0_real)
         call decompose2(x, E0_real, uv, E0_parl, E0_perp, j, error)
         call determinant(j, det_j, error)
         call calculate_action(z, s_i)
         s_i = aimag(s_i)
         phi = exp(cmplx(0.0, -1.0, dp)*s_i + cmplx(0.0, 1.0, dp)*aimag(det_j))/norm2(E0_perp)
      else
         call determinant(j, det_j, error)
         call calculate_action(z, s_i)
         s_i = aimag(s_i)
         phi = exp(cmplx(0.0, -1.0, dp)*s_i + cmplx(0.0, 1.0, dp)*aimag(det_j))
      end if

      ! Save the initial state
      write (10) x(1)
      write (20) z
      write (30) phi

      ! Markov chain generation
      do i = 2, chain_length

         ! Perform one step of the Metropolis procedure
         do k = 1, hmc_step
            call metropolis_step(x, z, j, total_step_size, num_steps, x_new, z_new, j_new, accept)

            if (accept) then
               ! Accept the new state
               x = x_new
               z = z_new
               j = j_new

               ! Update x_initial.dat
               call save_initial_condition(initial_file, x)
               accepted_steps = accepted_steps + 1
            end if

         end do

         if (wv) then
            call ds(z, ds_val)
            E0 = conjg(ds_val)
            call complex_to_real(E0, E0_real)
            call decompose2(x, E0_real, uv, E0_parl, E0_perp, j, error)
            call determinant(j, det_j, error)
            call calculate_action(z, s_i)
            s_i = aimag(s_i)
            phi = exp(cmplx(0.0, -1.0, dp)*s_i + cmplx(0.0, 1.0, dp)*aimag(det_j))/norm2(E0_perp)
         else
            call determinant(j, det_j, error)
            call calculate_action(z, s_i)
            s_i = aimag(s_i)
            phi = exp(cmplx(0.0, -1.0, dp)*s_i + cmplx(0.0, 1.0, dp)*aimag(det_j))
         end if
         ! Save the current state to the history files
         write (10) x(1)
         FLUSH (10)
         write (20) z
         FLUSH (20)
         write (30) phi
         FLUSH (30)

         if (mod(i, 10) == 0) then
            time_current = MPI_Wtime()
            time_elapsed = time_current - time_start
            time_per_iter = time_elapsed/real(i, dp)
            time_left = time_per_iter*real(chain_length - i, dp)
            progress = 100.0_dp*real(i, dp)/real(chain_length, dp)
            print *, "Progress: ", progress, "%  |  Elapsed:", time_elapsed, "sec  |  ETA:", time_left, "sec"
            print *, "Acceptance rate", real(accepted_steps)/(real(i)*hmc_step)
            print *, "current tolerance", cttol
            ckrv = .True.
         end if

      end do

      ! Close files
      close (10)
      close (20)
      close (30)
      print *, "Acceptance rate", real(accepted_steps)/(real(chain_length)*hmc_step)
   end subroutine generate_markov_chain

   !==============================================================
   ! Subroutine: execute_generate_markov_chain
   !==============================================================
   subroutine execute_generate_markov_chain()
      implicit none
      real(dp), allocatable :: x_initial(:)

      ! Read bω parameters (for completeness, though may not be used directly)
      print *, "Reading bw parameters from:", bw_file
      call read_bw_parameters(bw_file, bw, size(bw))
      print *, "Current bw:", bw

      ! Read the initial state vector
      print *, "Reading initial state from:", initial_file
      call read_initial_condition(initial_file, x_initial)
      x_initial(1) = T0
      ! x_initial(17:) = 0

      ! Generate the Markov chain
      print *, "Generating Markov chain..."
      call generate_markov_chain(chain_length_bw, x_initial, total_step_size, num_steps, &
                                 x_history_file, z_history_file, phi_history_file)
      print *, "Markov chain generation complete."
      print *, "  x-history =>", x_history_file
      print *, "  z-history =>", z_history_file
      print *, "  phi-history =>", phi_history_file

      ! Cleanup
      if (allocated(x_initial)) deallocate (x_initial)
   end subroutine execute_generate_markov_chain

   !==============================================================
   ! Function: update_bw
   !   Returns .true. if we should continue cycling,
   !           .false. if we should stop.
   !==============================================================
   logical function update_bw()
      implicit none

      integer :: i, j
      integer :: bw_size, n_segments
      real(dp), allocatable :: x_history(:)
      real(dp), allocatable :: bw_shifted(:)
      real(dp) :: min_bw, average_value, ratio
      real(dp), allocatable :: bin_edges(:), bw2(:)
      integer, allocatable :: bin_counts(:)

      ! Default return value is .true.; if condition fails, we set it to .false.
      update_bw = .true.
      allocate (x_history(chain_length_bw))
      ! Read current bω from file and x_history
      call read_bw_parameters(bw_file, bw, size(bw))
      call read_x_history(x_history_file, x_history)
      bw_size = size(bw)
      n_segments = bw_size - 1

      print *, "Initial bω:", bw
      print *, "Number of steps in x_history =", chain_length_bw

      ! Extract x(1) values

      ! Build bin edges for [T0, T1] across bw_size bins
      allocate (bin_edges(bw_size + 1))
      do i = 1, bw_size + 1
         bin_edges(i) = T0 + (i - 1)*(T1 - T0)/real(bw_size, dp)
      end do
      bin_edges(1) = bin_edges(1) - d0
      bin_edges(bw_size + 1) = bin_edges(bw_size + 1) + d0

      ! Initialize counts and accumulate
      allocate (bin_counts(bw_size), bw2(bw_size))
      bin_counts = 0

      do i = 1, chain_length_bw
         do j = 1, bw_size
            if (x_history(i) >= bin_edges(j) .and. x_history(i) < bin_edges(j + 1)) then
               bin_counts(j) = bin_counts(j) + 1
               exit
            end if
         end do
         if (x_history(i) == T1) then
            bin_counts(bw_size) = bin_counts(bw_size) + 1
         end if
      end do

      print *, "Bin edges:", bin_edges
      print *, "Bin counts:", bin_counts

      !------------------------------------------------------------
      !  Check condition:
      !    (1/n_segments)*Σ [ (h(l+1)-h(l)) / (h(l+1)-2*h(l)) ]^2  > delta^2 ?
      !  We skip or define ratio=0 if denominators are zero, etc.
      !------------------------------------------------------------
      if (bw_size < 2) then
         print *, "WARNING: Less than 2 bins => cannot do ratio check."
         ! We'll just keep going, no update or skip?
         ! update_bw remains .true.
      else
         average_value = 0.0_dp
         do i = 1, n_segments
            if (bin_counts(i) == 0 .and. bin_counts(i + 1) == 0) then
               ratio = 0.0_dp
            else if ((bin_counts(i + 1) - 2*bin_counts(i)) == 0) then
               ratio = 0.0_dp
            else
               ratio = real(bin_counts(i + 1) - bin_counts(i), dp)/ &
                       real(bin_counts(i + 1) + bin_counts(i), dp)*2.0_dp
            end if
            average_value = average_value + ratio**2
         end do
         average_value = average_value/real(n_segments, dp)

         print *, "Check: average(ratio^2) =", average_value, " vs delta^2 =", delta**2
      end if

      if (average_value <= delta**2) then
         print *, "Condition NOT satisfied => skipping bw update."
         update_bw = .false.
         deallocate (x_history, bin_edges, bin_counts)
         return
      end if
      !------------------------------------------------------------
      ! Update bω if condition is satisfied
      !------------------------------------------------------------

      do i = 1, bw_size
         call W_function(T0 + (i - 0.5_dp)*((T1 - T0)/real(size(bw), dp)), bw2(i))
      end do
      bw = bw2
      do i = 1, bw_size
         bw(i) = bw(i) + 0.1*log(real(bin_counts(i), dp) + 1.0_dp)
      end do
      ! Shift so min(bw) = 0
      min_bw = minval(bw)
      allocate (bw_shifted(bw_size))
      do i = 1, bw_size
         bw_shifted(i) = bw(i) - min_bw
      end do
      bw = bw_shifted

      print *, "Shifted bω:", bw
      call save_bw(bw_file, bw, bw_size)
      print *, "Updated bω saved to:", bw_file

      ! Cleanup
      deallocate (x_history, bin_edges, bin_counts, bw_shifted, bw2)
   end function update_bw
end module markovchain
