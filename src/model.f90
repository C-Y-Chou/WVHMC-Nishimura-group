module model
   use param_mod
   use utils
   use mt95

   implicit none

contains

   ! ======================== RANDOM NUMBER GENERATION ==========================
   ! Generates Gaussian random numbers using the Box-Muller method
   subroutine grand(gaus_rand)
      real(dp), dimension(:), intent(inout) :: gaus_rand
      integer :: n, i

      n = size(gaus_rand)
      do i = 1, n
         gaus_rand(i) = gaussrnd()
      end do
   end subroutine grand

   ! ======================== ACTION AND DERIVATIVES ==========================
   ! Computes the derivative of the action S(z) with respect to z
   subroutine calculate_action(z, s)
      ! Input/Output parameters
      complex(dp), dimension(:), intent(in) :: z
      complex(dp), intent(out) :: s

      ! Local variables
      complex(dp) :: detm
      complex(dp), allocatable :: shifted_z(:)
      integer :: i, j, k, l

      ! Allocate memory for shifted_z and initialize it with z
      allocate (shifted_z(size(z)))
      shifted_z = z*gm2

      ! Initialize action
      s = 0.0_dp

      ! Compute the action
      do i = 1, 3
         s = s + cmplx(0.0_dp, -1.0_dp, dp)*(shifted_z(i)**2)*2*gm
      end do

      do i = 4, DN
         s = s + cmplx(0.0_dp, -1.0_dp, dp)*(-shifted_z(i)**2)*2*gm
      end do

      do i = 1, 3
         do j = 1, DN/3 - 1
         do k = 1, 3
            if (k == i) cycle
            s = s + cmplx(0.0_dp, 1.0_dp, dp)* &
                (shifted_z(i)**2*shifted_z(3*j + k)**2 - &
                 shifted_z(i)*shifted_z(k)*shifted_z(3*j + i)*shifted_z(3*j + k))*8
         end do
         end do
      end do

      do j = 1, DN/3 - 1
         do l = j, DN/3 - 1
         do i = 1, 3
            do k = 1, 3
               if (k == i) cycle
               s = s + cmplx(0.0_dp, -1.0_dp, dp)* &
                   (shifted_z(3*l + i)**2*shifted_z(3*j + k)**2 - &
                    shifted_z(3*l + i)*shifted_z(3*l + k)*shifted_z(3*j + i)*shifted_z(3*j + k))*8
            end do
         end do
         end do
      end do

      do i = 1, DN/3 - 1
         s = s + alpha*shifted_z(DN + i)* &
             (shifted_z(1)*shifted_z(3*i + 1) + &
              shifted_z(2)*shifted_z(3*i + 2) + &
              shifted_z(3)*shifted_z(3*i + 3))*2
      end do

      do i = 1, DN/3 - 1
         s = s + beta* &
             ((shifted_z(1)*shifted_z(3*i + 1) + &
               shifted_z(2)*shifted_z(3*i + 2) + &
               shifted_z(3)*shifted_z(3*i + 3))*2)**2
      end do

      call m(shifted_z, detm)
      s = s - detm

      ! Deallocate memory
      deallocate (shifted_z)

   end subroutine calculate_action

   subroutine m(z, s)
      complex(kind=dp), dimension(:), intent(in)::z
      complex(kind=dp), intent(out)::s
      complex(kind=dp), dimension(DN/3 - 1, DN/3 - 1)::fpm
      integer::i, j, k
      logical :: error_det
      fpm = 0
      do k = 1, 3
         do i = 1, DN/3 - 1
            fpm(i, i) = fpm(i, i) + 2*z(k)**2
            do j = 1, DN/3 - 1
               fpm(i, j) = fpm(i, j) + 2*z(3*i + k)*z(3*j + k)
            end do
         end do
      end do
      call determinant(fpm, s, error_det)

      if (error_det) then
         print *, "Determinant error"
      end if
      ! s = log(s)
   end subroutine

   subroutine ds(z, s)
      ! Input and Output parameters
      complex(dp), dimension(:), intent(in) :: z
      complex(dp), dimension(:), intent(out) :: s

      ! Local variables
      complex(dp) :: ddetm(DN)
      complex(dp), allocatable :: shifted_z(:)
      integer :: i, j, k, l

      ! Allocate and shift z
      allocate (shifted_z(size(z)))
      shifted_z = z*gm2

      ! Initialize output array
      s = 0.0_dp

      ! First computation loop
      do i = 1, 3
         s(i) = s(i) + cmplx(0.0_dp, -1.0_dp, dp)*(shifted_z(i)*2)*2*gm
      end do
      do i = 4, DN
         s(i) = s(i) + cmplx(0.0_dp, -1.0_dp, dp)*(-shifted_z(i)*2)*2*gm
      end do

      ! Nested loops for interaction terms
      do i = 1, 3
         do j = 1, DN/3 - 1
         do k = 1, 3
            if (k == i) cycle
            s(i) = s(i) + cmplx(0.0_dp, 1.0_dp, dp)* &
                   (shifted_z(i)*2*shifted_z(3*j + k)**2 - shifted_z(k)*shifted_z(3*j + i)*shifted_z(3*j + k))*8
            s(3*j + k) = s(3*j + k) + cmplx(0.0_dp, 1.0_dp, dp)* &
                         (shifted_z(i)**2*shifted_z(3*j + k)*2 - shifted_z(i)*shifted_z(k)*shifted_z(3*j + i))*8
            s(k) = s(k) + cmplx(0.0_dp, 1.0_dp, dp)* &
                   (-shifted_z(i)*shifted_z(3*j + i)*shifted_z(3*j + k))*8
            s(3*j + i) = s(3*j + i) + cmplx(0.0_dp, 1.0_dp, dp)* &
                         (-shifted_z(i)*shifted_z(k)*shifted_z(3*j + k))*8
         end do
         end do
      end do

      do j = 1, DN/3 - 1
         do l = j, DN/3 - 1
         do i = 1, 3
            do k = 1, 3
               if (k == i) cycle
               s(3*l + i) = s(3*l + i) + cmplx(0.0_dp, -1.0_dp, dp)* &
                            (shifted_z(3*l + i)*2*shifted_z(3*j + k)**2 &
                             - shifted_z(3*l + k)*shifted_z(3*j + i)*shifted_z(3*j + k))*8
               s(3*j + k) = s(3*j + k) + cmplx(0.0_dp, -1.0_dp, dp)* &
                            (shifted_z(3*l + i)**2*shifted_z(3*j + k)*2 &
                             - shifted_z(3*l + i)*shifted_z(3*l + k)*shifted_z(3*j + i))*8
               s(3*l + k) = s(3*l + k) + cmplx(0.0_dp, -1.0_dp, dp)* &
                            (-shifted_z(3*l + i)*shifted_z(3*j + i)*shifted_z(3*j + k))*8
               s(3*j + i) = s(3*j + i) + cmplx(0.0_dp, -1.0_dp, dp)* &
                            (-shifted_z(3*l + i)*shifted_z(3*l + k)*shifted_z(3*j + k))*8
            end do
         end do
         end do
      end do

      ! Additional computation loop
      do i = 1, 3
         do j = 1, DN/3 - 1
            s(i) = s(i) + alpha*shifted_z(DN + j)*shifted_z(3*j + i)*2
            s(3*j + i) = s(3*j + i) + alpha*shifted_z(DN + j)*shifted_z(i)*2
         end do
      end do

      do i = 1, DN/3 - 1
         s(DN + i) = s(DN + i) + alpha* &
                     (shifted_z(1)*shifted_z(3*i + 1) + shifted_z(2)*shifted_z(3*i + 2) + shifted_z(3)*shifted_z(3*i + 3))*2
      end do

      do i = 1, DN/3 - 1
      do j = 1, 3
         s(j) = s(j) + beta* &
                ((shifted_z(1)*shifted_z(3*i + 1) + &
                  shifted_z(2)*shifted_z(3*i + 2) + &
                  shifted_z(3)*shifted_z(3*i + 3))*2)*2* &
                (shifted_z(3*i + j)*2)
         s(3*i + j) = s(3*i + j) + beta* &
                      ((shifted_z(1)*shifted_z(3*i + 1) + &
                        shifted_z(2)*shifted_z(3*i + 2) + &
                        shifted_z(3)*shifted_z(3*i + 3))*2)*2* &
                      (shifted_z(j)*2)
      end do
      end do

      ! Call dm to compute determinant contribution
      call dm(shifted_z, ddetm)

      ! Update s with determinant result
      s(1:DN) = s(1:DN) - ddetm
      s = s*gm2
      ! Deallocate shifted_z
      deallocate (shifted_z)
   end subroutine ds

   subroutine dm(z, s)

      ! Input parameters
      complex(dp), dimension(:), intent(in) :: z

      ! Output parameters
      complex(dp), dimension(DN), intent(out) :: s

      ! Local variables
      complex(dp), dimension(DN/3 - 1, DN/3 - 1) :: fpm, fpm2
      complex(dp), dimension(DN/3 - 1, DN/3 - 1, DN) :: dfpm
      integer :: i, j, k, ipiv(DN/3 - 1), info
      complex(dp), allocatable :: work(:)

      ! Allocate work array based on problem size
      allocate (work((DN/3 - 1)*DN/3 - 1))

      ! Initialization
      fpm = 0.0_dp
      fpm2 = 0.0_dp
      dfpm = 0.0_dp

      ! Loop over k, updating fpm
      do k = 1, 3
         do i = 1, DN/3 - 1
            fpm(i, i) = fpm(i, i) + 2.0_dp*z(k)**2
            do j = 1, DN/3 - 1
               fpm(i, j) = fpm(i, j) + 2.0_dp*z(3*i + k)*z(3*j + k)
            end do
         end do
      end do

      ! Store original fpm for reuse
      fpm2 = fpm

      ! LU factorization and inversion of fpm
      call zgetrf(DN/3 - 1, DN/3 - 1, fpm, DN/3 - 1, ipiv, info)
      if (info /= 0) print *, "Error in LU factorization: ", info
      call zgetri(DN/3 - 1, fpm, DN/3 - 1, ipiv, work, size(work), info)
      if (info /= 0) print *, "Error in matrix inversion: ", info

      ! Compute derivatives of fpm
      do k = 1, 3
         do i = 1, DN/3 - 1
            dfpm(i, i, k) = dfpm(i, i, k) + 4.0_dp*z(k)
            do j = 1, DN/3 - 1
               dfpm(i, j, 3*i + k) = dfpm(i, j, 3*i + k) + 2.0_dp*z(3*j + k)
               dfpm(j, i, 3*i + k) = dfpm(j, i, 3*i + k) + 2.0_dp*z(3*j + k)
            end do
         end do
      end do

      ! Compute output s
      s = 0.0_dp
      fpm2 = 0.0_dp
      do i = 1, DN
         fpm2 = matmul(fpm, dfpm(:, :, i))
         do j = 1, DN/3 - 1
            s(i) = s(i) + fpm2(j, j)
         end do
      end do

      ! Deallocate work array
      deallocate (work)

   end subroutine dm

   subroutine hessian(z, h)
      ! Input and Output parameters
      complex(dp), dimension(:), intent(in) :: z
      complex(dp), dimension(:, :), intent(out) :: h

      ! Local variables
      integer :: i, j, k, l
      complex(dp), allocatable :: shifted_z(:)
      complex(dp) :: dddetm(DN, DN)

      ! Allocate and shift z
      allocate (shifted_z(size(z)))
      shifted_z = z*gm2

      ! Initialize Hessian matrix
      h = 0.0_dp

      ! Diagonal terms for the first loop
      do i = 1, 3
         h(i, i) = h(i, i) + cmplx(0.0_dp, -1.0_dp, dp)*(2)*2*gm
      end do

      do i = 4, DN
         h(i, i) = h(i, i) + cmplx(0.0_dp, -1.0_dp, dp)*(-2)*2*gm
      end do

      ! Interaction terms
      do i = 1, 3
         do j = 1, DN/3 - 1
         do k = 1, 3
            if (k == i) cycle
            h(i, i) = h(i, i) + cmplx(0.0_dp, 1.0_dp, dp)*(2*shifted_z(3*j + k)**2)*8
            h(i, 3*j + k) = h(i, 3*j + k) + cmplx(0.0_dp, 1.0_dp, dp)* &
                            (shifted_z(i)*2*shifted_z(3*j + k)*2 - shifted_z(k)*shifted_z(3*j + i))*8
            h(i, k) = h(i, k) + cmplx(0.0_dp, 1.0_dp, dp)* &
                      (-shifted_z(3*j + i)*shifted_z(3*j + k))*8
            h(i, 3*j + i) = h(i, 3*j + i) + cmplx(0.0_dp, 1.0_dp, dp)* &
                            (-shifted_z(k)*shifted_z(3*j + k))*8
            h(3*j + k, 3*j + k) = h(3*j + k, 3*j + k) + cmplx(0.0_dp, 1.0_dp, dp)*(shifted_z(i)**2*2)*8
            h(3*j + k, i) = h(3*j + k, i) + cmplx(0.0_dp, 1.0_dp, dp)* &
                            (shifted_z(i)*2*shifted_z(3*j + k)*2 - shifted_z(k)*shifted_z(3*j + i))*8
            h(3*j + k, k) = h(3*j + k, k) + cmplx(0.0_dp, 1.0_dp, dp)* &
                            (-shifted_z(i)*shifted_z(3*j + i))*8
            h(3*j + k, 3*j + i) = h(3*j + k, 3*j + i) + cmplx(0.0_dp, 1.0_dp, dp)* &
                                  (-shifted_z(i)*shifted_z(k))*8
            h(k, i) = h(k, i) + cmplx(0.0_dp, 1.0_dp, dp)* &
                      (-shifted_z(3*j + i)*shifted_z(3*j + k))*8
            h(k, 3*j + i) = h(k, 3*j + i) + cmplx(0.0_dp, 1.0_dp, dp)* &
                            (-shifted_z(i)*shifted_z(3*j + k))*8
            h(k, 3*j + k) = h(k, 3*j + k) + cmplx(0.0_dp, 1.0_dp, dp)* &
                            (-shifted_z(i)*shifted_z(3*j + i))*8
            h(3*j + i, i) = h(3*j + i, i) + cmplx(0.0_dp, 1.0_dp, dp)* &
                            (-shifted_z(k)*shifted_z(3*j + k))*8
            h(3*j + i, k) = h(3*j + i, k) + cmplx(0.0_dp, 1.0_dp, dp)* &
                            (-shifted_z(i)*shifted_z(3*j + k))*8
            h(3*j + i, 3*j + k) = h(3*j + i, 3*j + k) + cmplx(0.0_dp, 1.0_dp, dp)* &
                                  (-shifted_z(i)*shifted_z(k))*8
         end do
         end do
      end do

      do j = 1, DN/3 - 1
         do l = j, DN/3 - 1
         do i = 1, 3
            do k = 1, 3
               if (k == i) cycle
               h(3*l + i, 3*l + i) = h(3*l + i, 3*l + i) + cmplx(0.0_dp, -1.0_dp, dp)*(2*shifted_z(3*j + k)**2)*8
               h(3*l + i, 3*j + k) = h(3*l + i, 3*j + k) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (shifted_z(3*l + i)*2*shifted_z(3*j + k)*2 - shifted_z(3*l + k)*shifted_z(3*j + i))*8
               h(3*l + i, 3*j + i) = h(3*l + i, 3*j + i) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (-shifted_z(3*l + k)*shifted_z(3*j + k))*8
               h(3*l + i, 3*l + k) = h(3*l + i, 3*l + k) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (-shifted_z(3*j + i)*shifted_z(3*j + k))*8
               h(3*j + k, 3*j + k) = h(3*j + k, 3*j + k) + cmplx(0.0_dp, -1.0_dp, dp)*(shifted_z(3*l + i)**2*2)*8
               h(3*j + k, 3*l + i) = h(3*j + k, 3*l + i) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (shifted_z(3*l + i)*2*shifted_z(3*j + k)*2 - shifted_z(3*l + k)*shifted_z(3*j + i))*8
               h(3*j + k, 3*l + k) = h(3*j + k, 3*l + k) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (-shifted_z(3*l + i)*shifted_z(3*j + i))*8
               h(3*j + k, 3*j + i) = h(3*j + k, 3*j + i) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (-shifted_z(3*l + i)*shifted_z(3*l + k))*8
               h(3*l + k, 3*l + i) = h(3*l + k, 3*l + i) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (-shifted_z(3*j + i)*shifted_z(3*j + k))*8
               h(3*l + k, 3*j + i) = h(3*l + k, 3*j + i) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (-shifted_z(3*l + i)*shifted_z(3*j + k))*8
               h(3*l + k, 3*j + k) = h(3*l + k, 3*j + k) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (-shifted_z(3*l + i)*shifted_z(3*j + i))*8
               h(3*j + i, 3*l + i) = h(3*j + i, 3*l + i) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (-shifted_z(3*l + k)*shifted_z(3*j + k))*8
               h(3*j + i, 3*l + k) = h(3*j + i, 3*l + k) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (-shifted_z(3*l + i)*shifted_z(3*j + k))*8
               h(3*j + i, 3*j + k) = h(3*j + i, 3*j + k) + cmplx(0.0_dp, -1.0_dp, dp)* &
                                     (-shifted_z(3*l + i)*shifted_z(3*l + k))*8
            end do
         end do
         end do
      end do

      do i = 1, 3
         do j = 1, DN/3 - 1
            h(i, DN + j) = h(i, DN + j) + alpha*(shifted_z(3*j + i))*2
            h(DN + j, i) = h(i, DN + j)
            h(i, 3*j + i) = h(i, 3*j + i) + alpha*shifted_z(DN + j)*(1)*2
            h(3*j + i, i) = h(i, 3*j + i)
            h(i + 3*j, DN + j) = h(i + 3*j, DN + j) + alpha*(shifted_z(i))*2
            h(DN + j, i + 3*j) = h(i + 3*j, DN + j)
         end do
      end do

      do i = 1, DN/3 - 1
      do j = 1, 3
         h(j, j) = h(j, j) + beta* &
                   ((shifted_z(3*i + j))*2)*2* &
                   ((shifted_z(3*i + j))*2)
         h(3*i + j, 3*i + j) = h(3*i + j, 3*i + j) + beta* &
                               ((shifted_z(j))*2)*2* &
                               ((shifted_z(j))*2)
         h(j, 3*i + j) = h(j, 3*i + j) + beta* &
                         ((shifted_z(j)*shifted_z(3*i + j))*2)*2* &
                         (2)
         h(3*i + j, j) = h(j, 3*i + j)
         do k = 1, 3
            h(j, 3*i + j) = h(j, 3*i + j) + beta* &
                            ((shifted_z(k)*shifted_z(3*i + k))*2)*2* &
                            (2)
            h(3*i + j, j) = h(j, 3*i + j)
            if (j == k) cycle
            h(j, 3*i + k) = h(j, 3*i + k) + beta* &
                            (shifted_z(k)*2)* &
                            (shifted_z(3*i + j)*2)
            h(3*i + k, j) = h(j, 3*i + k)
            h(j, k) = h(j, k) + beta* &
                      (shifted_z(3*i + k)*2)* &
                      (shifted_z(3*i + j)*2)
            h(k, j) = h(j, k)
            h(3*i + j, k) = h(3*i + j, k) + beta* &
                            (shifted_z(3*i + k)*2)* &
                            (shifted_z(j)*2)
            h(k, 3*i + j) = h(3*i + j, k)
            h(3*i + j, 3*i + k) = h(3*i + j, 3*i + k) + beta* &
                                  (shifted_z(k)*2)* &
                                  (shifted_z(j)*2)
            h(3*i + k, 3*i + j) = h(3*i + j, 3*i + k)
         end do
      end do
      end do
      ! Call ddm to compute second derivative contributions
      call ddm(shifted_z, dddetm)
      ! Subtract determinant contributions
      h(1:DN, 1:DN) = h(1:DN, 1:DN) - dddetm
      h = h*gm2**2
      ! Deallocate shifted_z
      deallocate (shifted_z)

   end subroutine hessian

   subroutine ddm(z, s)
      ! Input parameters
      complex(dp), dimension(:), intent(in) :: z

      ! Output parameters
      complex(dp), dimension(DN, DN), intent(out) :: s

      ! Local variables
      complex(dp), dimension(DN/3 - 1, DN/3 - 1) :: fpm, fpm2
      complex(dp), dimension(DN/3 - 1, DN/3 - 1, DN) :: dfpm
      complex(dp), dimension(DN/3 - 1, DN/3 - 1, DN, DN) :: ddfpm
      integer :: i, j, k, ipiv(DN/3 - 1), info
      complex(dp), allocatable :: work(:)

      ! Allocate work array based on problem size
      allocate (work((DN/3 - 1)*DN/3 - 1))

      ! Initialization
      s = 0.0_dp
      fpm = 0.0_dp
      fpm2 = 0.0_dp
      dfpm = 0.0_dp
      ddfpm = 0.0_dp

      ! Compute fpm
      do k = 1, 3
         do i = 1, DN/3 - 1
            fpm(i, i) = fpm(i, i) + 2.0_dp*z(k)**2
            do j = 1, DN/3 - 1
               fpm(i, j) = fpm(i, j) + 2.0_dp*z(3*i + k)*z(3*j + k)
            end do
         end do
      end do

      ! Store original fpm for reuse
      fpm2 = fpm

      ! LU factorization and inversion of fpm
      call zgetrf(DN/3 - 1, DN/3 - 1, fpm, DN/3 - 1, ipiv, info)
      if (info /= 0) print *, "Error in LU factorization: ", info
      call zgetri(DN/3 - 1, fpm, DN/3 - 1, ipiv, work, size(work), info)
      if (info /= 0) print *, "Error in matrix inversion: ", info

      ! Compute derivatives of fpm
      do k = 1, 3
         do i = 1, DN/3 - 1
            dfpm(i, i, k) = dfpm(i, i, k) + 4.0_dp*z(k)
            do j = 1, DN/3 - 1
               dfpm(i, j, 3*i + k) = dfpm(i, j, 3*i + k) + 2.0_dp*z(3*j + k)
               dfpm(j, i, 3*i + k) = dfpm(j, i, 3*i + k) + 2.0_dp*z(3*j + k)
            end do
         end do
      end do

      ! Compute second derivatives of fpm
      do k = 1, 3
         do i = 1, DN/3 - 1
            ddfpm(i, i, k, k) = ddfpm(i, i, k, k) + 4.0_dp
            do j = 1, DN/3 - 1
               ddfpm(i, j, 3*i + k, 3*j + k) = ddfpm(i, j, 3*i + k, 3*j + k) + 2.0_dp
               ddfpm(j, i, 3*i + k, 3*j + k) = ddfpm(j, i, 3*i + k, 3*j + k) + 2.0_dp
            end do
         end do
      end do

      ! Compute output s
      do i = 1, DN
         do j = 1, DN
            fpm2 = matmul(matmul(matmul(fpm, dfpm(:, :, j)), fpm), dfpm(:, :, i))
            fpm2 = -fpm2 + matmul(fpm, ddfpm(:, :, i, j))
            do k = 1, DN/3 - 1
               s(i, j) = s(i, j) + fpm2(k, k)
            end do
         end do
      end do

      ! Deallocate work array
      deallocate (work)

   end subroutine ddm

   subroutine W_function(t, W)
      ! Inputs
      real(dp), intent(in)  :: t
      ! Output
      real(dp), intent(out) :: W

      ! Locals
      real(dp), allocatable :: a(:), coeffs(:)
      real(dp) :: h
      integer :: i, n

      n = size(bw)
      h = (T1 - T0)/real(n, dp)

      ! Create x-values in [T0, T1]
      allocate (a(n))
      do i = 1, n
         a(i) = T0 + (i - 0.5_dp)*h
      end do
      ! We fit a polynomial of order (n-2)
      allocate (coeffs(order + 1))
      call least_squares_fit(a, bw, order, coeffs)
      ! Evaluate the fitted polynomial at t
      W = 0.0_dp
      do i = 0, order
         W = W + coeffs(i + 1)*t**i
      end do

      deallocate (a, coeffs)

      ! if (t > t1) then
      !    W = -gw*(t - t0) + cw*(exp((t - t1)**2.0_dp/2.0_dp/d0**2.0_dp) - 1.0_dp)
      ! else if (t < t0) then
      !    W = -gw*(t - t0) + cw*(exp((t - t0)**2.0_dp/2.0_dp/d0**2.0_dp) - 1.0_dp)
      ! else
      !    W = -gw*(t - t0)
      ! end if
   end subroutine W_function

   ! -------------------------------------------------------------
   subroutine least_squares_fit(x, y, poly_order, coeffs)
      !
      !  x(:):         the independent variable (length = n)
      !  y(:):         the data to be fit (length = n)
      !  poly_order:   polynomial order to fit (e.g., n-2)
      !  coeffs(:):    output array for the fitted coefficients

      real(dp), intent(in)  :: x(:), y(:)
      integer, intent(in)  :: poly_order
      real(dp), intent(out) :: coeffs(:)

      ! Locals
      integer :: n, m, lda, ldb, info
      real(dp), allocatable :: design_matrix(:, :), rhs(:)
      real(dp), allocatable :: work(:)
      integer :: lwork, i, j

      n = size(x)               ! number of data points
      m = poly_order + 1        ! number of unknowns = (order+1)

      ! 1) Allocate matrix A (of shape n x m) for design matrix,
      !    and vector RHS (of length n) which will hold y data
      allocate (design_matrix(n, m))
      allocate (rhs(n))

      ! Fill design matrix with x(i)**(j-1), j=1..m
      do i = 1, n
         do j = 1, m
            design_matrix(i, j) = x(i)**(j - 1)
         end do
      end do

      ! Set up the right-hand side
      rhs = y

      lda = n   ! leading dimension of design_matrix
      ldb = n   ! leading dimension of rhs (since it’s n x 1)

      ! 2) Query DGELS for optimal workspace size
      lwork = -1
      allocate (work(1))   ! temporary
      call dgels('N', n, m, 1, design_matrix, lda, rhs, ldb, work, lwork, info)
      if (info /= 0) then
         print *, "DGELS workspace query failed with INFO = ", info
         stop
      end if

      ! The first entry of WORK now holds the optimal size
      lwork = int(work(1))
      deallocate (work)
      allocate (work(lwork))

      ! 3) Call DGELS to perform the least-squares fit.
      !    On exit, 'rhs' will contain the solution in its first 'm' entries.
      call dgels('N', n, m, 1, design_matrix, lda, rhs, ldb, work, lwork, info)
      if (info /= 0) then
         print *, "DGELS failed with INFO = ", info
         stop
      end if

      ! 4) Extract polynomial coefficients from 'rhs'.
      !    The solution is in rhs(1:m). The rest of rhs can hold residuals.
      coeffs = rhs(1:m)

      ! Cleanup
      deallocate (design_matrix, rhs, work)
   end subroutine least_squares_fit

   subroutine W_t_derivative(t, W)
      ! Inputs
      real(dp), intent(in)  :: t
      ! Output
      real(dp), intent(out) :: W

      ! Locals
      real(dp), allocatable :: a(:), coeffs(:)
      real(dp) :: h
      integer :: i, n

      n = size(bw)
      h = (T1 - T0)/real(n, dp)

      ! Create x-values in [T0, T1]
      allocate (a(n))
      do i = 1, n
         a(i) = T0 + (i - 0.5_dp)*h
      end do
      ! We fit a polynomial of order (n-2)
      allocate (coeffs(order + 1))
      call least_squares_fit(a, bw, order, coeffs)
      ! Evaluate the fitted polynomial at t
      W = 0.0_dp
      do i = 0, order
         W = W + coeffs(i + 1)*t**(i - 1)*i
      end do

      deallocate (a, coeffs)
      ! if (t > t1) then
      !    W = -gw + cw*exp((t - t1)**2.0_dp/2.0_dp/d0**2.0_dp)*(t - t1)/d0**2.0_dp
      ! else if (t < t0) then
      !    W = -gw + cw*exp((t - t0)**2.0_dp/2.0_dp/d0**2.0_dp)*(t - t0)/d0**2.0_dp
      ! else
      !    W = -gw
      ! end if
   end subroutine W_t_derivative
end module model
