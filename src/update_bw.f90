!======================  combined_code.f90  ======================

program update_bw_pro
   ! This main program will do repeated cycles of:
   !   1) Generate the Markov chain (execute_generate_markov_chain)
   !   2) Update bω based on the generated chain (update_bw)
   ! until either update_bw signals we should stop or max_cycles is reached.

   use param_mod         ! Shared parameters and filenames
   use solve_flow
   use model
   use hmc
   use utils
   use markovchain
   use mpi
   implicit none

   logical :: keep_going
   integer :: seed, ierr, cycle

   seed = getseed()
   call sgrnd(seed)
   call read_parameters()
   call MPI_Init(ierr)

   ! Optionally, do any setup code here (e.g., printing a header).
   print *, "======================================="
   print *, " Starting cyclical Markov-chain update "
   print *, "======================================="
   cycle = 0
   keep_going = .true.

   do while (keep_going)
      cycle = cycle + 1
      print *, " "
      print *, "================================="
      print *, "   CYCLE = ", cycle
      print *, "================================="

      ! 1) Generate the Markov chain
      call execute_generate_markov_chain()

      ! 2) Update bω; if update_bw() returns .false., we break the loop
      keep_going = update_bw()

   end do

   print *, "Update signaled stop => bω update condition not satisfied."
   print *, char(7)
   call MPI_Finalize(ierr)
end program update_bw_pro
