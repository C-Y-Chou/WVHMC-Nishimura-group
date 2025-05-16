program test_rattle_hamiltonian_conservation
   use param_mod
   use mt95
   use hmc
   use utils
   use mpi
   implicit none

   ! Parameters
   integer, parameter :: max_steps = 10  ! Maximum number of substeps for testing
   integer :: num_step  ! Number of substeps
   real(kind=8), dimension(:), allocatable :: x_state, x_state_prime
   complex(kind=8), dimension(:), allocatable :: z_state, z_state_prime
   complex(kind=8), dimension(:, :), allocatable :: j_matrix, j_matrix_prime
   real(kind=8) :: h_initial, h_final  ! Initial and final Hamiltonians
   real(kind=8) :: ti, tf
   real(kind=8), dimension(max_steps) :: hamiltonian_difference  ! Difference in Hamiltonians
   logical :: rattle_error, flow_error_flag
   integer :: seed
   complex(dp) :: det
   integer :: i
   integer :: ierr
   call MPI_Init(ierr)


   ! Get the current system time as an integer seed
   ti = MPI_Wtime()
   ! Initialize the RNG with the current time as the seed
   seed = getseed()
   call sgrnd(seed)
   call read_parameters()
   ! Allocate variables
    allocate(x_state(n_size -1  + 1), z_state(n_size -1 ), z_state_prime(n_size -1 ), &
    j_matrix(n_size -1 , n_size -1 ), j_matrix_prime(n_size -1 , n_size -1 ), x_state_prime(n_size -1  + 1))
   call grand(testmom)
   ! Randomize x_state and ensure x_state(1) is positive
   call read_initial_condition(initial_file, x_state)
   call read_bw_parameters(bw_file, bw, size(bw))
   print *, "Current bw:", bw
   istest = .True.
   ! Determine z_state and j_matrix using the flow subroutine
   call flow(x_state, z_state, j_matrix, flow_error_flag)
   if (flow_error_flag) then
      print *, "Error: Flow subroutine failed."
   end if
   ! Loop over the number of substeps to test Hamiltonian conservation
   do num_step = 1, max_steps
      ckrv = .true.
      ! Call the RATTLE subroutine
      call rattle(x_state, z_state, total_step_size, num_step, x_state_prime, &
      z_state_prime, h_initial, h_final, j_matrix, j_matrix_prime)
      print*,x_state_prime(1)
      if (h_final==0) then
         print *, "Error: RATTLE subroutine failed for num_step =", num_step
         hamiltonian_difference(num_step) = -1.0_8  ! Mark as invalid
         cycle
      end if

      ! Compute the difference in Hamiltonians
      hamiltonian_difference(num_step) = abs(h_final - h_initial)
      tf = MPI_Wtime()
      print *, "ETA :", (tf - ti)*(1.0/(real((1 + num_step)*num_step)/(1 + max_steps)/max_steps) - 1), "sec"," ", "del_h = "&
       ,abs(h_final - h_initial)
   end do
   print *, "Time elapsed:", (tf - ti), "s"
   print *, "order" , log(hamiltonian_difference(1)/hamiltonian_difference(max_steps))/log(real(max_steps))

   call save_and_plot_hamiltonian_loglog(max_steps, hamiltonian_difference)
   call MPI_Finalize(ierr)
   ! Deallocate variables
   deallocate (x_state, z_state, z_state_prime, j_matrix, j_matrix_prime, x_state_prime)
end program test_rattle_hamiltonian_conservation

! Subroutine to save Hamiltonian difference data and generate a log-log plot
!---------------------------------------------------------------------
!  save_and_plot_hamiltonian_loglog — head-less version
!
!  * Writes data to  hamitonian_conservation.dat
!  * Writes a small gnuplot script that produces
!      hamiltonian_conservation_loglog.png
!  * Runs gnuplot non-interactively and reports any failure.
!---------------------------------------------------------------------
subroutine save_and_plot_hamiltonian_loglog(max_entries, h_data)
   use iso_fortran_env, only : output_unit
   implicit none
   integer,   intent(in) :: max_entries
   real(8),   intent(in) :: h_data(max_entries)

   integer :: i, istat

   ! ---------- 1. dump the data -------------------------------------
   open(unit=20, file="hamiltonian_conservation.dat", status="replace", action="write")
   write(20, '(A)') "# num_step  Hamiltonian_Difference"
   do i = 1, max_entries
      write(20, '(I6,1X,E22.15)') i, h_data(i)
   end do
   close(20)

   ! ---------- 2. prepare a *head-less* gnuplot script --------------
open(unit=30, file="hamiltonian_conservation_loglog.gp", status="replace", action="write")
write(30,'(A)') "set terminal pngcairo size 1200,800 enhanced font 'Helvetica,14'"
write(30,'(A)') "set output 'hamiltonian_conservation_loglog.png'"
   write(30,'(A)')  "set logscale x"
   write(30,'(A)')  "set logscale y"
   write(30,'(A)')  "set xlabel 'Number of Steps'"
   write(30,'(A)')  "set ylabel 'Hamiltonian |ΔH|'"
   write(30,'(A)')  "set title 'RATTLE Hamiltonian Conservation'"
   write(30,'(A)')  "set grid"
   write(30,'(A)')  "plot 'hamiltonian_conservation.dat' using 1:2 with linespoints lw 2 title 'ΔH'"
   write(30,'(A)')  "unset output"
   close(30)

   ! ---------- 3. run gnuplot (non-interactive) ----------------------
   print *, "Generating Hamiltonian conservation log-log PNG..."
   call execute_command_line("gnuplot hamiltonian_conservation_loglog.gp", exitstat=istat)

   if (istat /= 0) then
      write(output_unit,'(A,I0)') "Error: gnuplot exited with status ", istat
   else
      write(output_unit,'(A)') "Plot written to hamiltonian_conservation_loglog.png"
   end if
end subroutine save_and_plot_hamiltonian_loglog