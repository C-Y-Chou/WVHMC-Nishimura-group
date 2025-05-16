module param_mod
   use utils   ! Ensure utils module is properly implemented for any auxiliary functions
   implicit none

   ! Constants
   real(dp) :: T0, T1               ! Start and End time
   integer :: bwn
   real(dp), ALLOCATABLE :: bw(:)                 ! Coefficients for W(t)
   integer  :: order                  ! Order of W(t) polynomial approximation
   real(dp) :: d0, gw, cw
   logical  :: istest, wv, ckrv, tra2, eo             ! Test mode flag
   real(dp), allocatable :: testmom(:)

   ! Markov Chain Parameters
   integer  :: chain_length_bw
   integer  :: chain_length
   integer  :: hmc_step
   real(dp) :: total_step_size
   integer  :: num_steps
   real(dp) :: cttol
   real(dp) :: at, rt

   ! Update_bw parameters
   integer  :: n_size, n_warm
   real(dp) :: delta

   integer  :: DN, nh
   real(dp) :: gm, gm2
   complex(dp) :: alpha, beta
   real(dp), allocatable :: M_diag(:), Minv_diag(:)

   ! File names
   character(len=256) :: initial_file
   character(len=256) :: bw_file
   character(len=256) :: x_history_file
   character(len=256) :: z_history_file
   character(len=256) :: phi_history_file

contains

   subroutine read_parameters()
      implicit none
      integer :: ios
      character(len=256) :: filename
      filename = "../data/parameters.dat"  ! Adjust based on actual file location

      open (unit=10, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print *, "Error: Unable to open parameter file: ", trim(filename)
         stop
      end if

      ! Read values from file
      read (10, *) T0, T1
      read (10, *) order
      read (10, *) bwn
      read (10, *) d0, gw, cw
      read (10, *) istest, wv, tra2, eo
      read (10, *) chain_length_bw, chain_length, hmc_step
      read (10, *) total_step_size, num_steps
      read (10, *) n_size, delta
      read (10, *) n_warm
      read (10, *) DN, alpha, beta, gm, gm2
      read (10, *) at, rt, cttol
      read (10, '(A)') initial_file
      read (10, '(A)') bw_file
      read (10, '(A)') x_history_file
      read (10, '(A)') z_history_file
      read (10, '(A)') phi_history_file

      close (10)
      gm = gm*sqrt(real(DN/3, dp))
      ALLOCATE (bw(bwn))
      chain_length_bw = bwn*chain_length_bw
      allocate (M_diag(2*(n_size - 1)), Minv_diag(2*(n_size - 1)))
      M_diag = 1
      Minv_diag = 1
      nh = 1
      allocate (testmom((n_size - 1)*2))
   end subroutine read_parameters
end module param_mod
