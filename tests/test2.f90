program test_hessian_and_ds_comparison
   use model
   use param_mod
   implicit none

   complex(dp), dimension(:), ALLOCATABLE :: z, ds_result, numerical_ds
   complex(dp), dimension(:), ALLOCATABLE :: s_p2, s_p1, s_m1, s_m2
   complex(dp), dimension(:, :), ALLOCATABLE :: hessian_result, numerical_hessian
   complex(dp), dimension(:), ALLOCATABLE :: ds_p2, ds_p1, ds_m1, ds_m2
   complex(dp) :: s_val, epsilon_t
   real(dp) :: norm_ds_diff, norm_hessian_diff
   real(dp), dimension(:), ALLOCATABLE :: rand_real, rand_imag
   integer :: i, j
   complex(dp) :: z_orig
   real(dp) :: h_opt, ti, tf
    integer :: seed_size
    integer, allocatable :: seed(:)

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    call system_clock(count=i)
    seed = [(i + 37 * j, j = 1, seed_size)]  ! Better variation using array constructor
    call random_seed(put=seed)

   call read_parameters()
   ALLOCATE (z(n_size - 1), ds_result(n_size - 1), numerical_ds(n_size - 1))
   ALLOCATE (hessian_result(n_size - 1, n_size - 1), numerical_hessian(n_size - 1, n_size - 1))
   ALLOCATE (rand_real(n_size - 1), rand_imag(n_size - 1))
   ALLOCATE (s_p2(n_size - 1), s_p1(n_size - 1), s_m1(n_size - 1), s_m2(n_size - 1))
   ALLOCATE (ds_p2(n_size - 1), ds_p1(n_size - 1), ds_m1(n_size - 1), ds_m2(n_size - 1))

   ! Generate random values for initialization
   call random_number(rand_real) ! Real part
   call random_number(rand_imag) ! Imaginary part

   ! Scale and combine into a complex array
   z = [(cmplx(rand_real(i), rand_imag(i), dp), i=1, n_size - 1)]

   ! Compute the optimal epsilon_t based on the machine precision.
   ! For the five-point stencil, optimal h ~ (machine epsilon)^(1/5).
   h_opt = epsilon(1.0_dp)**(0.2_dp)
   ! We use epsilon_t as a complex number (perturbation in the imaginary direction)
   epsilon_t = cmplx(0.0_dp, h_opt, dp)

   ! --- Step 1: Compare Numerical and Analytical ds using five–point approximation ---
   ! Compute the analytical derivative:
   call CPU_TIME(ti)
   call ds(z, ds_result)
   call CPU_TIME(tf)
   print*,tf-ti
   call CPU_TIME(ti)

   numerical_ds = (0.0_dp, 0.0_dp)
   norm_ds_diff = 0.0_dp

   do i = 1, n_size - 1
      ! Save the original value of z(i)
      z_orig = z(i)

      ! Compute S(z) at four perturbed points using the five–point stencil:
      z(i) = z_orig + 2*epsilon_t
      call calculate_action(z, s_val)
      s_p2(i) = s_val

      z(i) = z_orig + epsilon_t
      call calculate_action(z, s_val)
      s_p1(i) = s_val

      z(i) = z_orig - epsilon_t
      call calculate_action(z, s_val)
      s_m1(i) = s_val

      z(i) = z_orig - 2*epsilon_t
      call calculate_action(z, s_val)
      s_m2(i) = s_val

      ! Reset the value of z(i)
      z(i) = z_orig

      ! Use the five–point formula for the first derivative:
      numerical_ds(i) = (-s_p2(i) + 8*s_p1(i) - 8*s_m1(i) + s_m2(i))/(12.0_dp*epsilon_t)
      norm_ds_diff = norm_ds_diff + abs(ds_result(i) - numerical_ds(i))**2
   end do

   norm_ds_diff = sqrt(norm_ds_diff)

   print *, "Component   ds_result           numerical_ds         Difference"
   do i = 1, n_size - 1
      print *, i, ds_result(i), numerical_ds(i), abs(ds_result(i) - numerical_ds(i))
   end do
   
   ! --- Step 2: Compare Numerical and Analytical Hessian using five–point approximation ---
   ! Compute the analytical Hessian:
   call hessian(z, hessian_result)
   numerical_hessian = (0.0_dp, 0.0_dp)
   norm_hessian_diff = 0.0_dp

   do i = 1, n_size - 1
      ! Save the original value for the i-th component of z
      z_orig = z(i)

      ! Perturb only the i-th component and compute the ds vector at four offset values.
      z(i) = z_orig + 2*epsilon_t
      call ds(z, ds_p2)
      z(i) = z_orig + epsilon_t
      call ds(z, ds_p1)
      z(i) = z_orig - epsilon_t
      call ds(z, ds_m1)
      z(i) = z_orig - 2*epsilon_t
      call ds(z, ds_m2)
      z(i) = z_orig  ! Reset

      do j = 1, n_size - 1
         ! Approximate the second derivative of S with respect to z(i) and z(j)
         numerical_hessian(i, j) = (-ds_p2(j) + 8*ds_p1(j) - 8*ds_m1(j) + ds_m2(j))/(12.0_dp*epsilon_t)
         norm_hessian_diff = norm_hessian_diff + abs(hessian_result(i, j) - numerical_hessian(i, j))**2
      end do
   end do

   norm_hessian_diff = sqrt(norm_hessian_diff)

   print *, "Component   hessian_result      numerical_hessian     Difference"
   do i = 1, n_size - 1
      do j = 1, n_size - 1
         print *, i, j, hessian_result(i, j), numerical_hessian(i, j), &
            abs(hessian_result(i, j) - numerical_hessian(i, j))
      end do
   end do

   print *, "Norm of ds difference: ", norm_ds_diff
   print *, "Norm of Hessian difference: ", norm_hessian_diff

end program test_hessian_and_ds_comparison
