!======================================================
! File: execute_generate_markov_chain.f90
!======================================================
program generate_markov_chain_pro
   use param_mod        ! Get chain_length, total_step_size, num_steps,
   !    initial_file, bw_file, x_history_file, ...
   use model            ! Provides: real(dp), T0, T1, and bw(:)
   use hmc
   use utils
   use markovchain
   use mpi
   implicit none
   ! Variables
   real(dp), allocatable :: x_initial(:)  ! Initial state vector
   integer :: seed, ierr

   ! Initialize the RNG with the current time as the seed
   seed = getseed()
   call sgrnd(seed)
   call read_parameters()
   call MPI_Init(ierr)

   ! Read the bw parameters from a file
   if (wv) then
      print *, "Reading bw parameters from file:", bw_file
      call read_bw_parameters(bw_file, bw, size(bw))
      print *, "bw parameters:", bw
   else
      bw = 0
   end if

   ! Read the initial state vector from a file
   print *, "Reading initial state from file:", initial_file
   call read_initial_condition(initial_file, x_initial)
   print *, size(x_initial)
   !  if (x_initial(1) <= 0.0_dp) then
   !     print *, "Error: x_initial(1) must be positive. Exiting."
   !     stop
   !  end if
   ! Execute the Markov chain generation using the integrated parameters
   print *, "Generating Markov chain..."
   call generate_markov_chain(chain_length, x_initial, total_step_size, num_steps, &
                              x_history_file, z_history_file, phi_history_file)

   print *, "Markov chain generation complete."
   print *, "History of x saved to", x_history_file
   print *, "History of z saved to", z_history_file
   print *, "History of phi saved to", phi_history_file

   ! Deallocate x_initial
   if (allocated(x_initial)) deallocate (x_initial)
   call MPI_Finalize(ierr)
   print *, char(7)
end program generate_markov_chain_pro
