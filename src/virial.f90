program virial
   use param_mod
   use utils
   use model
   implicit none

   integer :: Dz, i, j, maxW
   integer :: iotatus, z_rows, phiows
   real(dp) :: maxJKError
   complex(dp) :: ratio
   complex(dp), allocatable :: Ob(:), numData(:), denData(:)
   real(dp), allocatable :: walues(:), jkRatios(:, :)
   real(dp) :: ratioCIRealLowBB, ratioCIRealHighBB
   real(dp) :: ratioCIImagLowBB, ratioCIImagHighBB
   ! For the block bootstrap
   integer :: nBoot, B, maxLag
   complex(dp) :: ratioMeanBB
   real(dp) :: ratioErrReBB, ratioErrImBB
   real(dp) :: tauInt  ! integrated autocorrelation time
   complex(dp), allocatable :: z_all(:, :), phi_all(:)
   integer :: nows, iostat

   ! 1) Read user parameters and data
   call read_parameters()
   Dz = n_size - 1
   call read_z_history(z_history_file, z_all, nows, iostat)
   call read_phi_history(phi_history_file, phi_all, nows, iostat)
   allocate (Ob(size(z_all, 2)))
   allocate (numData(size(phi_all)))
   allocate (denData(size(phi_all)))
   ! 2) Compute observable array O(z) for each row
   !$omp parallel do default(shared) private(j)
   do j = 1, size(z_all, 2)
      Ob(j) = O(z_all(:, j))
   end do
   !$omp end parallel do
   ! 3) Build the ratio data arrays
   !    numData(i) = Ob(i) * phi(i)
   !    denData(i) = phi(i)
   numData = Ob*phi_all
   denData = phi_all

   maxLag = 1  ! a guess. can be smaller or bigger depending on N
   tauInt = 1
   do while (tauInt < size(numData))
      call autocorr_time_complex(numData/sum(denData), tauInt, maxLag)
      if (tauInt < 0) then
         exit
      elseif (maxLag < tauInt*9) then
         maxLag = 10*tauInt
      else
         exit
      end if
   end do

   print *, "Estimated autocorrelation time =", tauInt

   maxW = tauInt*2
   ratio = sum(numData)/sum(denData)

   print *, "Real Ratio:", real(ratio), "Imag Ratio:", aimag(ratio)

   !=====================================================
   ! 4) BLOCK BOOTSTRAP with automatic block-size selection
   !=====================================================

   ! 4a) Estimate the autocorrelation time from (for example) the numerator
   !     data's real part. Alternatively, you could do the full complex or the
   !     imaginary part, or the entire ratio time series if you store it.
   !
   !     We define a large 'maxLag' so that the correlation is likely negligible
   !     beyond that lag.  You can experiment with the choice.

   ! 4b) Pick block size B ~ a few times tauInt. Must be integer, bounded by array size
   B = int(2.0_dp*tauInt)  ! e.g. 2*tauInt
   if (B < 1) B = 1
   if (B > size(numData)) B = size(numData)

   print *, "Chosen block size B =", B

   ! 4c) Perform block bootstrap
   nBoot = 10000
   call ratioAndBlockBootstrap(numData, denData, nBoot, B, &
                               ratioMeanBB, ratioErrReBB, ratioErrImBB, &
                               ratioCIRealLowBB, ratioCIRealHighBB, &
                               ratioCIImagLowBB, ratioCIImagHighBB)
   print *, "Block-Bootstrap ratio (real)  =", real(ratioMeanBB)
   print *, "Block-Bootstrap ratio (imag)  =", aimag(ratioMeanBB)
   print *, "Block-Bootstrap error (real)  =", ratioErrReBB
   print *, "Block-Bootstrap error (imag)  =", ratioErrImBB
   print *, "68% CI for Re(ratio) = [", ratioCIRealLowBB, ",", ratioCIRealHighBB, "]"
   print *, "68% CI for Im(ratio) = [", ratioCIImagLowBB, ",", ratioCIImagHighBB, "]"
   print *, "point needed:", NINT(1/abs(sum(denData)/size(denData))**2*((B + 1)/2))
   print *, "current point:", size(denData)
   call save_virial("../output/virial.dat", Ob)
   call execute_command_line("gnuplot ../output/plot.gp")

contains
   function O(z_row) result(res)
      complex(dp), intent(in) :: z_row(:)
      complex(dp) :: res
      if (tra2) then
         call calculate_a2(z_row, res)
      else
         call calculate_virial(z_row, res)
      end if
   end function

   subroutine calculate_a2(z, s)
      ! Input/Output parameters
      complex(dp), dimension(:), intent(in) :: z
      complex(dp), intent(out) :: s

      ! Local variables
      complex(dp) :: detm
      complex(dp), allocatable :: shifted_z(:)
      integer :: i, k, l

      ! Allocate memory for shifted_z and initialize it with z
      allocate (shifted_z(size(z)))
      shifted_z = z*gm2

      ! Initialize action
      s = 0.0_dp
      ! Compute the action
      do i = 1, 3
         s = s + (shifted_z(i)**2)
      end do

      do i = 4, DN
         s = s + (-shifted_z(i)**2)
      end do
      s = -s/sqrt(5.0)*2

      ! Deallocate memory
      deallocate (shifted_z)
   end subroutine calculate_a2

   subroutine calculate_virial(z, s)
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
         s = s + cmplx(0.0_dp, -1.0_dp, dp)*(shifted_z(i)**2)*2*gm*2
      end do

      do i = 4, DN
         s = s + cmplx(0.0_dp, -1.0_dp, dp)*(-shifted_z(i)**2)*2*gm*2
      end do

      do i = 1, 3
         do j = 1, DN/3 - 1
         do k = 1, 3
            if (k == i) cycle
            s = s + cmplx(0.0_dp, 1.0_dp, dp)* &
                (shifted_z(i)**2*shifted_z(3*j + k)**2 - &
                 shifted_z(i)*shifted_z(k)*shifted_z(3*j + i)*shifted_z(3*j + k))*8*4
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
                    shifted_z(3*l + i)*shifted_z(3*l + k)*shifted_z(3*j + i)*shifted_z(3*j + k))*8*4
            end do
         end do
         end do
      end do

      do i = 1, DN/3 - 1
         s = s + alpha*shifted_z(DN + i)* &
             (shifted_z(1)*shifted_z(3*i + 1) + &
              shifted_z(2)*shifted_z(3*i + 2) + &
              shifted_z(3)*shifted_z(3*i + 3))*2*2
      end do

      do i = 1, DN/3 - 1
         s = s + beta* &
             ((shifted_z(1)*shifted_z(3*i + 1) + &
               shifted_z(2)*shifted_z(3*i + 2) + &
               shifted_z(3)*shifted_z(3*i + 3))*2)**2*4
      end do
      s = s/(DN + (DN/3 - 1)*2) - 1
      deallocate (shifted_z)

   end subroutine calculate_virial

   subroutine read_z_history(filename, z_all, nows, iotatus)
      use iso_fortran_env, only: dp => real64
      implicit none
      character(len=*), intent(in) :: filename
      complex(dp), allocatable, intent(out) :: z_all(:, :)
      integer, intent(out) :: nows, iotatus

      integer :: unit, filesize_bytes
      integer, parameter :: element_size = 2*8  ! 16 bytes per complex(dp)
      integer :: n_z
      complex(dp), allocatable :: flat_z(:)

      n_z = n_size - 1
      unit = 20

      ! Open with correct mode
      open (unit=unit, file=filename, access='stream', form='unformatted', status='old', iostat=iotatus)
      if (iotatus /= 0) then
         print *, "Error opening z file: ", filename
         return
      end if

      ! Get file size in bytes
      inquire (unit=unit, size=filesize_bytes)
      if (mod(filesize_bytes, n_z*element_size) /= 0) then
         print *, "File size not divisible by z vector size"
         iotatus = 1
         close (unit)
         return
      end if

      nows = filesize_bytes/(n_z*element_size)
      allocate (flat_z(n_z*nows), stat=iotatus)
      if (iotatus /= 0) then
         print *, "Allocation failed for flat_z"
         close (unit)
         return
      end if

      ! Read entire file
      read (unit) flat_z
      close (unit)

      ! Reshape into z_all(n_z, nows)
      allocate (z_all(n_z, nows), stat=iotatus)
      if (iotatus /= 0) then
         print *, "Allocation failed for z_all"
         return
      end if

      z_all = reshape(flat_z, [n_z, nows])
   end subroutine
   subroutine read_phi_history(filename, phi_all, nows, iotatus)
      implicit none
      character(len=*), intent(in) :: filename
      complex(dp), allocatable, intent(out) :: phi_all(:)
      integer, intent(out) :: nows, iotatus

      integer :: unit, filesize_bytes
      integer, parameter :: element_size = 2*8  ! complex(dp)

      unit = 21
      open (unit=unit, file=filename, access='stream', form='unformatted', status='old', iostat=iotatus)
      if (iotatus /= 0) return

      inquire (unit=unit, size=filesize_bytes)
      nows = filesize_bytes/element_size

      allocate (phi_all(nows), stat=iotatus)
      if (iotatus /= 0) return

      read (unit) phi_all
      close (unit)
   end subroutine
   subroutine ratioAndError(numData, denData, maxW, walues, jkRatios)
      complex(dp), intent(in) :: numData(:), denData(:)
      integer, intent(in) :: maxW
      real(dp), intent(out) :: walues(:), jkRatios(:)
      complex(dp) :: fullRatio
      integer :: w, i, n_blocks, blockize
      complex(dp), allocatable :: blockNum(:), blockDen(:)

      do w = 1, maxW
         blockize = w
         if (blockize < 1) exit

         allocate (blockNum(size(numData)/w), blockDen(size(numData)/w))
         n_blocks = size(numData)/w

         do i = 1, n_blocks
            blockNum(i) = sum(numData((i - 1)*blockize + 1:i*blockize))
            blockDen(i) = sum(denData((i - 1)*blockize + 1:i*blockize))
         end do
         blockNum = sum(blockNum) - blockNum
         blockDen = sum(blockDen) - blockDen
         fullRatio = sum(blockNum/blockDen)/size(blockNum)

         jkRatios(w) = sqrt(sum(real(blockNum/blockDen - fullRatio)**2)/n_blocks*(n_blocks - 1))
         walues(w) = w

         deallocate (blockNum, blockDen)
      end do
   end subroutine

   subroutine ratioAndErrorim(numData, denData, maxW, walues, jkRatios)
      complex(dp), intent(in) :: numData(:), denData(:)
      integer, intent(in) :: maxW
      real(dp), intent(out) :: walues(:), jkRatios(:)
      complex(dp) :: fullRatio
      integer :: w, i, n_blocks, blockize
      complex(dp), allocatable :: blockNum(:), blockDen(:)

      do w = 1, maxW
         blockize = w
         if (blockize < 1) exit

         allocate (blockNum(size(numData)/w), blockDen(size(numData)/w))
         n_blocks = size(numData)/w

         do i = 1, n_blocks
            blockNum(i) = sum(numData((i - 1)*blockize + 1:i*blockize))
            blockDen(i) = sum(denData((i - 1)*blockize + 1:i*blockize))
         end do
         blockNum = sum(blockNum) - blockNum
         blockDen = sum(blockDen) - blockDen
         fullRatio = sum(blockNum/blockDen)/size(blockNum)

         jkRatios(w) = sqrt(sum(aimag(blockNum/blockDen - fullRatio)**2)/n_blocks*(n_blocks - 1))
         walues(w) = w

         deallocate (blockNum, blockDen)
      end do
   end subroutine
   subroutine save_virial(file_name, walues)
      character(len=*), intent(in) :: file_name
      complex(dp), intent(in) :: walues(:)
      integer :: i

      open (30, file=file_name, status='replace')
      do i = 1, size(walues)
         write (30, *) real(walues(i)), aimag(walues(i))
      end do
      close (30)
   end subroutine
   subroutine ratioAndBlockBootstrap(numData, denData, nBoot, blockSize, ratioMean, ratioErrReal, ratioErrImag, &
                                     ratioCIRealLow, ratioCIRealHigh, ratioCIImagLow, ratioCIImagHigh)
      implicit none

      integer, intent(in)                :: nBoot, blockSize
      complex(dp), intent(in)            :: numData(:), denData(:)
      complex(dp), intent(out)           :: ratioMean
      real(dp), intent(out)              :: ratioErrReal, ratioErrImag
      real(dp), intent(out) :: ratioCIRealLow, ratioCIRealHigh
      real(dp), intent(out) :: ratioCIImagLow, ratioCIImagHigh
      integer                            :: i, b, s, nData, nBlocks
      integer                            :: startIndex
      real(dp)                           :: rRand
      complex(dp)                        :: sumNum, sumDen, tmpRatio
      complex(dp), allocatable           :: ratioVals(:)

      integer :: idxLow, idxHigh
      real(dp)                           :: rMean, iMean, rVar, iVar, rPart, iPart

      real(dp), allocatable :: realVals(:), imagVals(:)
      nData = size(numData)
      nBlocks = nData/blockSize   ! floor division

      if (nBlocks < 1) then
         print *, "ERROR: blockSize is too large compared to nData!"
         stop
      end if

      allocate (ratioVals(nBoot))

      ! ======================================================
      !  Loop over nBoot bootstrap resamples
      ! ======================================================
      do b = 1, nBoot

         sumNum = (0.0_dp, 0.0_dp)
         sumDen = (0.0_dp, 0.0_dp)

         ! Construct the entire "resampled" time series
         ! by picking nBlocks blocks, each of length blockSize
         do i = 1, nBlocks
            call random_number(rRand)
            ! pick a random start index in [1 .. N-blockSize+1]
            startIndex = 1 + int(rRand*real(nData - blockSize + 1, dp))

            ! accumulate the block
            do s = 0, blockSize - 1
               sumNum = sumNum + numData(startIndex + s)
               sumDen = sumDen + denData(startIndex + s)
            end do
         end do

         ! Compute ratio for this bootstrap sample
         if ((real(sumDen) /= 0.0_dp) .or. (aimag(sumDen) /= 0.0_dp)) then
            tmpRatio = sumNum/sumDen
         else
            tmpRatio = (0.0_dp, 0.0_dp)
         end if

         ratioVals(b) = tmpRatio
      end do

      ! ======================================================
      !  Compute mean ratio
      ! ======================================================
      ratioMean = (0.0_dp, 0.0_dp)
      do b = 1, nBoot
         ratioMean = ratioMean + ratioVals(b)
      end do
      ratioMean = ratioMean/real(nBoot, dp)

      ! ======================================================
      !  Compute stdev of Re(ratio) and Im(ratio)
      ! ======================================================
      rMean = real(ratioMean)
      iMean = aimag(ratioMean)
      rVar = 0.0_dp
      iVar = 0.0_dp

      do b = 1, nBoot
         rPart = real(ratioVals(b)) - rMean
         iPart = aimag(ratioVals(b)) - iMean
         rVar = rVar + rPart*rPart
         iVar = iVar + iPart*iPart
      end do

      if (nBoot > 1) then
         rVar = rVar/real(nBoot - 1, dp)
         iVar = iVar/real(nBoot - 1, dp)
      else
         rVar = 0.0_dp
         iVar = 0.0_dp
      end if

      ratioErrReal = sqrt(rVar)
      ratioErrImag = sqrt(iVar)

      allocate (realVals(nBoot))
      allocate (imagVals(nBoot))

! Extract real and imaginary parts
      do b = 1, nBoot
         realVals(b) = real(ratioVals(b))
         imagVals(b) = aimag(ratioVals(b))
      end do

! Sort the values
      call mergesort(realVals, 1, nBoot)
      call mergesort(imagVals, 1, nBoot)

      ! idxLow = max(1, int(0.025_dp*real(nBoot)))
      ! idxHigh = min(nBoot, int(0.975_dp*real(nBoot)))
      idxLow = int(0.16_dp*real(nBoot))    ! 16% lower bound
      idxHigh = int(0.84_dp*real(nBoot))    ! 84% upper bound
! Extract percentile bounds
      ratioCIRealLow = realVals(idxLow)
      ratioCIRealHigh = realVals(idxHigh)
      ratioCIImagLow = imagVals(idxLow)
      ratioCIImagHigh = imagVals(idxHigh)

      deallocate (realVals, imagVals)
      deallocate (ratioVals)
   end subroutine ratioAndBlockBootstrap
   subroutine autocorr_time_complex(Ob, tauInt, maxLag)
      implicit none

      complex(dp), intent(in) :: Ob(:)
      integer, intent(in) :: maxLag
      real(dp), intent(out):: tauInt

      integer :: N, t, k
      complex(dp) :: meanOb, tmp
      complex(dp), allocatable :: ObPrime(:)
      real(dp), allocatable   :: C(:)
      real(dp) :: norm0, normk

      ! 1) Basic setup
      N = size(Ob)
      if (N < 2 .or. maxLag < 1) then
         tauInt = 1.0_dp
         return
      end if

      ! 2) Compute the mean
      meanOb = (0.0_dp, 0.0_dp)
      do t = 1, N
         meanOb = meanOb + Ob(t)
      end do
      meanOb = meanOb/real(N, dp)

      ! 3) Subtract the mean -> ObPrime
      allocate (ObPrime(N))
      do t = 1, N
         ObPrime(t) = Ob(t) - meanOb
      end do

      ! 4) Allocate correlation array C(0..maxLag)
      allocate (C(0:maxLag))
      C(0) = 0.0_dp

      ! 5) c(0) = average of Re[ conj(ObPrime(t))*ObPrime(t) ], i.e. |Ob'(t)|^2
      do t = 1, N
         ! conj(ObPrime(t))*ObPrime(t) is |ObPrime(t)|^2, but let's keep it explicit:
         tmp = conjg(ObPrime(t))*ObPrime(t)
         C(0) = C(0) + real(tmp)
      end do
      C(0) = C(0)/real(N, dp)

      ! 6) c(k) for k=1..maxLag
      do k = 1, maxLag
         C(k) = 0.0_dp
         if (k < N) then
            do t = 1, N - k
               tmp = conjg(ObPrime(t))*ObPrime(t + k)
               C(k) = C(k) + real(tmp)
            end do
            C(k) = C(k)/real(N - k, dp)
         else
            ! If k >= N, no valid pairs, set c(k)=0
            C(k) = 0.0_dp
         end if
      end do

      ! 7) Sum up the normalized correlation for the integrated autocorrelation time
      !    tauInt = 1 + 2 * sum_{k=1..maxLag} [C(k)/C(0)].
      norm0 = C(0)
      if (norm0 <= 0.0_dp) then
         ! Degenerate case: all zero or something unphysical
         tauInt = 1.0_dp
      else
         tauInt = 1.0_dp
         do k = 1, maxLag
            normk = C(k)/norm0
            tauInt = tauInt + 2.0_dp*normk
         end do
      end if
      ! 8) Cleanup
      deallocate (C, ObPrime)
   end subroutine autocorr_time_complex
   recursive subroutine mergesort(x, left, right)
      ! ----------------------------------------------------------------
      ! Sorts x(left:right) in ascending order using MergeSort.
      ! stable, robust, O(n log n) worst-case
      ! ----------------------------------------------------------------
      implicit none
      integer, intent(in)       :: left, right
      real(dp), intent(inout)   :: x(:)
      integer                   :: mid

      if (left >= right) return   ! 0 or 1 element => already sorted

      mid = (left + right)/2
      call mergesort(x, left, mid)
      call mergesort(x, mid + 1, right)
      call merge_sub(x, left, mid, right)
   end subroutine mergesort

   recursive subroutine merge_sub(x, left, mid, right)
      ! ----------------------------------------------------------------
      ! Merges the two sorted segments:
      !   x(left:mid) and x(mid+1:right)
      ! into a single sorted segment x(left:right).
      ! ----------------------------------------------------------------
      implicit none
      integer, intent(in)      :: left, mid, right
      real(dp), intent(inout)  :: x(:)
      integer                  :: i, j, k, nL, nR
      real(dp), allocatable :: L(:), R(:)

      ! Lengths of the left and right halves
      nL = mid - left + 1
      nR = right - mid

      ! Allocate temporary arrays
      allocate (L(nL), R(nR))

      ! Copy data into L and R
      do i = 1, nL
         L(i) = x(left + i - 1)
      end do
      do j = 1, nR
         R(j) = x(mid + j)
      end do

      ! Merge L and R back into x(left:right)
      i = 1
      j = 1
      do k = left, right
         if (i > nL) then
            ! Left array exhausted, copy from right
            x(k) = R(j)
            j = j + 1
         else if (j > nR) then
            ! Right array exhausted, copy from left
            x(k) = L(i)
            i = i + 1
         else
            ! Compare L(i) and R(j)
            if (L(i) <= R(j)) then
               x(k) = L(i)
               i = i + 1
            else
               x(k) = R(j)
               j = j + 1
            end if
         end if
      end do

      deallocate (L, R)
   end subroutine merge_sub
end program virial
