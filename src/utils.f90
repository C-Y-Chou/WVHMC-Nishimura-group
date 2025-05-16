module utils
   use, intrinsic :: iso_fortran_env, only: real64
   integer, parameter :: dp = real64
contains

   !-----------------------------------------------------------------------
   !> Generate an identity matrix of size (size x size).
   !-----------------------------------------------------------------------
   function identity_matrix(size) result(identity)
      integer, intent(in) :: size
      real(dp), dimension(size, size) :: identity
      integer :: i

      identity = 0.0_dp
      do i = 1, size
         identity(i, i) = 1.0_dp
      end do
   end function identity_matrix
   !-----------------------------------------------------------------------
   !> Convert an (n x n) complex matrix cmat to a (2n x 2n) real matrix rmat.
   !  If cmat(i,j) = a + i*b, then:
   !    rmat(2i-1, 2j-1) = a
   !    rmat(2i-1, 2j  ) = -b
   !    rmat(2i,   2j-1) =  b
   !    rmat(2i,   2j  ) =  a
   !-----------------------------------------------------------------------
   subroutine map_to_real_mat(cmat, rmat)
      implicit none
      ! Inputs
      complex(dp), intent(in)  :: cmat(:, :)
      ! Outputs
      real(dp), intent(out) :: rmat(:, :)

      ! Locals
      integer :: n, i, j

      n = size(cmat, 1)
      if (size(cmat, 2) /= n) then
         write (*, *) "Error(map_to_real_mat): cmat is not square."
         return
      end if

      if (size(rmat, 1) /= 2*n .or. size(rmat, 2) /= 2*n) then
         write (*, *) "Error(map_to_real_mat): rmat must be (2n x 2n)."
         return
      end if

      rmat = 0.0_dp
      do i = 1, n
         do j = 1, n
            rmat(2*i - 1, 2*j - 1) = real(cmat(i, j), dp)
            rmat(2*i - 1, 2*j) = -aimag(cmat(i, j))
            rmat(2*i, 2*j - 1) = aimag(cmat(i, j))
            rmat(2*i, 2*j) = real(cmat(i, j), dp)
         end do
      end do
   end subroutine map_to_real_mat

   !-----------------------------------------------------------------------
   !> Convert a (2n x 2n) real matrix rmat (from map_to_real_mat) back to
   !  an (n x n) complex matrix cmat.
   !-----------------------------------------------------------------------
   subroutine map_to_complex_mat(rmat, cmat)
      implicit none
      ! Inputs
      real(dp), intent(in)  :: rmat(:, :)
      ! Outputs
      complex(dp), intent(out) :: cmat(:, :)

      ! Locals
      integer :: n, i, j

      n = size(cmat, 1)
      if (size(cmat, 2) /= n) then
         write (*, *) "Error(map_to_complex_mat): cmat is not square."
         return
      end if
      if (size(rmat, 1) /= 2*n .or. size(rmat, 2) /= 2*n) then
         write (*, *) "Error(map_to_complex_mat): rmat must be (2n x 2n)."
         return
      end if

      cmat = cmplx(0.0_dp, 0.0_dp, dp)
      do i = 1, n
         do j = 1, n
            cmat(i, j) = cmplx(rmat(2*i - 1, 2*j - 1), &
                               rmat(2*i, 2*j - 1), dp)
         end do
      end do
   end subroutine map_to_complex_mat

   !-----------------------------------------------------------------------
   !> Convert a complex vector c of length n to a real vector r of length 2n,
   !  in interleaved form: [Re(c(1)), Im(c(1)), Re(c(2)), Im(c(2)), ...]
   !-----------------------------------------------------------------------
   subroutine complex_to_real(c, r)
      implicit none
      ! Inputs
      complex(dp), intent(in)  :: c(:)
      ! Outputs
      real(dp), intent(out) :: r(:)

      ! Locals
      integer :: i, n

      n = size(c)
      if (size(r) /= 2*n) then
         write (*, *) "Error(complex_to_real): r must have length 2*n."
         return
      end if

      do i = 1, n
         r(2*i - 1) = real(c(i), dp)
         r(2*i) = aimag(c(i))
      end do
   end subroutine complex_to_real

   !-----------------------------------------------------------------------
   !> Convert a real vector r of length 2n (interleaved real/im parts) to
   !  a complex vector c of length n.
   !-----------------------------------------------------------------------
   subroutine real_to_complex(r, c)
      implicit none
      ! Inputs
      real(dp), intent(in)  :: r(:)
      ! Outputs
      complex(dp), intent(out) :: c(:)

      ! Locals
      integer :: i, n

      n = size(c)
      if (size(r) /= 2*n) then
         write (*, *) "Error(real_to_complex): r must have length 2*n."
         return
      end if

      do i = 1, n
         c(i) = cmplx(r(2*i - 1), r(2*i), dp)
      end do
   end subroutine real_to_complex

   !-----------------------------------------------------------------------
   !> Flatten an (n x n) complex matrix mat into a real vector vec of
   !  length 2*n*n, in row-major order:
   !  [Re(mat(1,1)), Im(mat(1,1)), Re(mat(1,2)), Im(mat(1,2)), ...]
   !-----------------------------------------------------------------------
   subroutine map_to_real(mat, vec)
      implicit none
      ! Inputs
      complex(dp), intent(in)  :: mat(:, :)
      ! Outputs
      real(dp), intent(out) :: vec(:)

      ! Locals
      integer :: i, j, n

      n = size(mat, 1)
      if (size(mat, 2) /= n) then
         write (*, *) "Error(map_to_real): mat is not square."
         return
      end if
      if (size(vec) /= 2*n*n) then
         write (*, *) "Error(map_to_real): vec must have length=2*n*n."
         return
      end if

      do i = 1, n
         do j = 1, n
            vec(2*((i - 1)*n + j) - 1) = real(mat(i, j), dp)
            vec(2*((i - 1)*n + j)) = aimag(mat(i, j))
         end do
      end do
   end subroutine map_to_real

   !-----------------------------------------------------------------------
   !> Zero out the imaginary parts in a (2n x 2n) real-block matrix
   !  produced by map_to_real_mat, effectively forcing it to be purely real.
   !-----------------------------------------------------------------------
   subroutine real_mat(mat)
      implicit none
      real(dp), intent(inout) :: mat(:, :)

      integer :: i, j, n

      n = size(mat, 1)/2
      if (size(mat, 2) /= 2*n) then
         write (*, *) "Error(real_mat): mat must be (2n x 2n)."
         return
      end if

      do i = 1, n
         do j = 1, n
            mat(2*i - 1, 2*j) = 0.0_dp
            mat(2*i, 2*j - 1) = 0.0_dp
         end do
      end do
   end subroutine real_mat

   !-----------------------------------------------------------------------
   !> Extract the imaginary part from a (2n x 2n) real-block matrix
   !  (originally from map_to_real_mat) and store it in the "real" slots.
   !  Replaces imaginary slots with 0. This is a trick to handle e.g.
   !  matrix( a + i b ) => ( b ) in the real part, zero out old a.
   !-----------------------------------------------------------------------
   subroutine im_mat(mat)
      implicit none
      real(dp), intent(inout) :: mat(:, :)

      integer :: i, j, n

      n = size(mat, 1)/2
      if (size(mat, 2) /= 2*n) then
         write (*, *) "Error(im_mat): mat must be (2n x 2n)."
         return
      end if

      do i = 1, n
         do j = 1, n
            ! Transfer the imaginary part to real slot
            mat(2*i - 1, 2*j - 1) = mat(2*i - 1, 2*j)
            mat(2*i, 2*j) = mat(2*i, 2*j - 1)
            ! Zero out the old imaginary positions
            mat(2*i - 1, 2*j) = 0.0_dp
            mat(2*i, 2*j - 1) = 0.0_dp
         end do
      end do
   end subroutine im_mat

   !-----------------------------------------------------------------------
   !> Force a real vector (2n) representing interleaved complex data
   !  [a1 b1 a2 b2 ...] to have zero imaginary parts => [a1 0 a2 0 ...].
   !-----------------------------------------------------------------------
   subroutine real_vec(vec)
      implicit none
      real(dp), intent(inout) :: vec(:)

      integer :: i, n

      n = size(vec)/2
      do i = 1, n
         vec(2*i) = 0.0_dp
      end do
   end subroutine real_vec

   !-----------------------------------------------------------------------
   !> Force a real vector (2n) [a1 b1 a2 b2 ...] to keep only imaginary
   !  parts => [b1 0 b2 0 ...].
   !-----------------------------------------------------------------------
   subroutine im_vec(vec)
      implicit none
      real(dp), intent(inout) :: vec(:)

      integer :: i, n

      n = size(vec)/2
      do i = 1, n
         vec(2*i - 1) = vec(2*i)
         vec(2*i) = 0.0_dp
      end do
   end subroutine im_vec

   subroutine conjg_vec(vec)
      implicit none
      real(dp), intent(inout) :: vec(:)

      integer :: i, n

      n = size(vec)/2
      do i = 1, n
         vec(2*i) = -vec(2*i)
      end do
   end subroutine conjg_vec

   subroutine timesi_vec(vec)
      implicit none
      real(dp), intent(inout) :: vec(:)
      real(dp) :: temp

      integer :: i, n

      n = size(vec)/2
      do i = 1, n
         temp = vec(2*i - 1)
         vec(2*i - 1) = -vec(2*i)
         vec(2*i) = temp
      end do
   end subroutine timesi_vec
   !-----------------------------------------------------------------------
   !> Reconstruct an (n x n) complex matrix mat from a real vector vec
   !  of length 2*n*n, which was presumably produced by map_to_real().
   !-----------------------------------------------------------------------
   subroutine map_to_complex(vec, mat)
      implicit none
      real(dp), intent(in)  :: vec(:)
      complex(dp), intent(out) :: mat(:, :)

      integer :: i, j, n

      n = size(mat, 1)
      if (size(mat, 2) /= n) then
         write (*, *) "Error(map_to_complex): mat is not square."
         return
      end if
      if (size(vec) /= 2*n*n) then
         write (*, *) "Error(map_to_complex): vec must have length 2*n*n."
         return
      end if

      do i = 1, n
         do j = 1, n
            mat(i, j) = cmplx(vec(2*((i - 1)*n + j) - 1), &
                              vec(2*((i - 1)*n + j)), dp)
         end do
      end do
   end subroutine map_to_complex

   !-----------------------------------------------------------------------
   !> Factorial function using recursion. Factorial(0)=1, etc.
   !-----------------------------------------------------------------------
   recursive function factorial(n) result(fact)
      integer, intent(in) :: n
      real(dp)            :: fact

      if (n < 0) then
         print *, "Error: Factorial is undefined for negative n."
         fact = 0.0_dp
      else if (n == 0) then
         fact = 1.0_dp
      else
         fact = real(n, dp)*factorial(n - 1)
      end if
   end function factorial

   !-----------------------------------------------------------------------
   !> Reads a 1D real array 'x' from file:
   !   1st line: integer n => size of x
   !   Next n lines: each line has a real number
   !-----------------------------------------------------------------------
   subroutine read_initial_condition(filename, x)
      character(len=*), intent(in) :: filename
      real(dp), allocatable, intent(out) :: x(:)
      integer :: n, i, ios

      open (unit=60, file=filename, status="old", action="read", iostat=ios)
      if (ios /= 0) then
         write (*, *) "Error(read_initial_condition): cannot open", filename
         return
      end if

      read (60, *, iostat=ios) n
      if (ios /= 0) then
         write (*, *) "Error: reading size of x."
         close (60)
         return
      end if

      allocate (x(n))
      do i = 1, n
         read (60, *, iostat=ios) x(i)
         if (ios /= 0) then
            write (*, *) "Error: reading x(", i, ") from file."
            deallocate (x)
            close (60)
            return
         end if
      end do
      close (60)
   end subroutine read_initial_condition

   !-----------------------------------------------------------------------
   !> Saves a 1D real array 'x' to a file:
   !   1st line: integer => size(x)
   !   Next lines: each line has 1 real number (E15.8 format).
   !-----------------------------------------------------------------------
   subroutine save_initial_condition(filename, x)
      character(len=*), intent(in) :: filename
      real(dp), dimension(:), intent(in) :: x
      integer :: i, ios

      open (unit=40, file=filename, status="replace", action="write", iostat=ios)
      if (ios /= 0) then
         write (*, *) "Error(save_initial_condition): cannot open", filename
         return
      end if

      write (40, *) size(x)
      do i = 1, size(x)
         write (40, '(E25.17)') x(i)
      end do
      close (40)
   end subroutine save_initial_condition

   !-----------------------------------------------------------------------
   !> Reads 8 real numbers (bw parameters) from a file.
   !-----------------------------------------------------------------------
   subroutine read_bw_parameters(filename, bw, n)
      character(len=*), intent(in)  :: filename
      integer, intent(in) :: n
      real(dp), dimension(n), intent(inout) :: bw
      integer :: i, ios

      open (unit=50, file=filename, status="old", action="read", iostat=ios)
      if (ios /= 0) then
         write (*, *) "Error(read_bw_parameters): cannot open", filename
         return
      end if

      do i = 1, n
         read (50, *, iostat=ios) bw(i)
         if (ios /= 0) then
            write (*, *) "Error(read_bw_parameters): reading param", i
            close (50)
            return
         end if
      end do
      close (50)
   end subroutine read_bw_parameters

   subroutine read_x_history(filename, x_history)
      implicit none

      character(len=*), intent(in) :: filename
      real(real64), intent(out) :: x_history(:)

      integer :: ios, unit, max_size
      real(real64), allocatable :: buffer(:)
      real(real64) :: temp
      integer :: i

      ! Open the file in stream, unformatted mode (same as writing)
      unit = 23
      inquire (iolength=max_size) temp  ! Get record size
      open (unit=unit, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
      if (ios /= 0) then
         write (*, *) "Error(read_x_history): cannot open", filename
         return
      end if

      read (unit) x_history

      close (unit)
   end subroutine read_x_history

   !-----------------------------------------------------------------------
   !> Save an array bw(:) of length n to a file, one per line in E15.8.
   !-----------------------------------------------------------------------
   subroutine save_bw(filename, bw, n)
      character(len=*), intent(in) :: filename
      real(dp), intent(in)        :: bw(:)
      integer, intent(in)         :: n
      integer :: i, ios

      open (unit=30, file=filename, status="replace", action="write", iostat=ios)
      if (ios /= 0) then
         write (*, *) "Error(save_bw): cannot open", filename
         return
      end if

      do i = 1, n
         write (30, '(E25.17)') bw(i)
      end do
      close (30)
   end subroutine save_bw

   !-----------------------------------------------------------------------
   !> Compute determinant of a complex (n x n) matrix via LU (zgetrf).
   !  det = product(diagonal) * (-1)^(# of row interchanges).
   !-----------------------------------------------------------------------
   subroutine determinant(matrix, det, error_flag)
      implicit none

      ! Inputs
      complex(dp), intent(in) :: matrix(:, :)

      ! Outputs
      complex(dp), intent(out) :: det
      logical, intent(out) :: error_flag

      ! Locals
      integer :: n, lda, info, i
      integer, allocatable :: ipiv(:)
      complex(dp), allocatable :: lu_matrix(:, :)
      integer :: sign
      error_flag = .false.
      sign = 1

      n = size(matrix, 1)
      if (n /= size(matrix, 2)) then
         write (*, *) "Error(determinant): matrix must be square."
         error_flag = .true.
         return
      end if
      lda = n

      allocate (ipiv(n), lu_matrix(n, n))
      lu_matrix = matrix

      ! Perform LU decomposition
      call zgetrf(n, n, lu_matrix, lda, ipiv, info)
      if (info /= 0) then
         write (*, *) "Error(determinant): LU decomposition failed, info=", info
         error_flag = .true.
         deallocate (ipiv, lu_matrix)
         return
      end if

      det = cmplx(0.0_dp, 0.0_dp, dp)
      do i = 1, n
         det = det + log(lu_matrix(i, i))
         if (ipiv(i) /= i) sign = -sign
      end do
      det = det + log(sign*cmplx(1, 0, dp))

      deallocate (ipiv, lu_matrix)
   end subroutine determinant

   subroutine shift(z)
      ! Input/Output parameters
      complex(dp), dimension(:), intent(inout) :: z

      ! Local variables
      real(dp), allocatable :: a(:)
      integer :: i, n
      character(len=50) :: filename

      ! Define the input file path
      filename = "../data/shift.dat"

      ! Determine the size of z and allocate a
      n = size(z)
      allocate (a(n))

      ! Read the shift values from the file
      open (unit=11, file=filename, status='old', action='read')
      read (11, *) a
      close (11)

      ! Apply the shift to each z(i)
      do i = 1, n
         z(i) = z(i) + a(i)*cmplx(0.0_dp, 1.0_dp, kind=dp)
      end do

      ! Deallocate the array
      deallocate (a)

   end subroutine shift
   subroutine compute_preconditioner(H_ref, P0)
      implicit none
      complex(dp), intent(in)  :: H_ref(:, :)
      complex(dp), intent(out) :: P0(:, :)
      complex(dp), allocatable :: U(:, :)
      complex(dp), allocatable :: A_copy(:, :)
      real(dp) :: epsilon
      integer :: info, lwork, n, i
      complex(dp), allocatable :: work(:)
      real(dp), allocatable    :: s(:)       ! singular values
      real(dp), allocatable    :: rwork(:)
      complex(dp), allocatable :: diag_inv_sqrt(:, :)
      complex(dp), allocatable :: VT(:, :)

      ! Regularization parameter to avoid numerical issues
      epsilon = 1e-15_dp
      n = size(H_ref, 1)

      allocate (s(n))
      allocate (U(n, n))
      allocate (VT(n, n))
      allocate (diag_inv_sqrt(n, n))
      allocate (A_copy(n, n))

      ! Make a copy of H_ref as ZGESVD will overwrite its input.
      A_copy = H_ref

      ! Add a small value to the diagonal to regularize A_copy
      do i = 1, n
         A_copy(i, i) = A_copy(i, i) + epsilon
      end do

      ! Workspace query for ZGESVD
      allocate (work(1))
      allocate (rwork(5*n))  ! Minimum size required for ZGESVD
      call zgesvd('A', 'A', n, n, A_copy, n, s, U, n, VT, n, work, -1, rwork, info)
      if (info /= 0) then
         print *, "Error in LAPACK ZGESVD workspace query, info =", info
      end if
      lwork = int(real(work(1)))
      deallocate (work)
      allocate (work(lwork))

      ! Perform the singular value decomposition: A_copy = U * diag(s) * VT
      call zgesvd('A', 'A', n, n, A_copy, n, s, U, n, VT, n, work, lwork, rwork, info)
      if (info /= 0) then
         print *, "Error in LAPACK ZGESVD, info =", info
      end if

      ! Build the diagonal matrix containing the inverse square roots of the singular values.
      diag_inv_sqrt = (0.0_dp, 0.0_dp)
      do i = 1, n
         diag_inv_sqrt(i, i) = 1.0_dp/s(i)
      end do

      ! Form the preconditioner P0 = U * diag_inv_sqrt * U^H.
      ! (For a symmetric/Hermitian matrix, the left singular vectors U can serve in place of the eigenvectors.)
      P0 = matmul(transpose(conjg(VT)), matmul(diag_inv_sqrt, VT))

      deallocate (s, U, VT, diag_inv_sqrt, work, rwork, A_copy)
   end subroutine compute_preconditioner

   !-----------------------------------------------------------------
   ! Function: outer_product
   ! Purpose : Compute the outer product of vectors a and b,
   !           i.e. result(i,j) = a(i)*b(j)
   !-----------------------------------------------------------------
   function outer_product(a, b) result(mat)
      implicit none
      real(dp), intent(in) :: a(:), b(:)
      real(dp), allocatable :: mat(:, :)
      integer :: n, m, i, j

      n = size(a)
      m = size(b)
      allocate (mat(n, m))

      do i = 1, n
         do j = 1, m
            mat(i, j) = a(i)*b(j)
         end do
      end do
   end function outer_product
end module utils
