! Subroutine sort2 from Numerical Recipes in Fortran 90
! Slightly modified to include passing-in of array sizes
!
! dwhipp - 04/08

      SUBROUTINE sort2(arr,n1,slave,n2)

      USE nrutil, ONLY : assert_eq

      IMPLICIT NONE

      ! Variable declaration
      integer, parameter :: sp = selected_real_kind(6, 37)

      INTEGER(KIND=sp) :: ndum,n1,n2
      INTEGER(KIND=sp),DIMENSION(n1) :: index
      REAL(KIND=sp),DIMENSION(n1),INTENT(INOUT) :: arr
      REAL(KIND=sp),DIMENSION(n2),INTENT(INOUT) :: slave
      ndum=assert_eq(size(arr),size(slave),'sort2')
      call indexx(arr,n1,index)
      arr=arr(index)
      slave=slave(index)
      END SUBROUTINE sort2