! Subroutine kptwo.f90 - 2 Sample Kuiper test
! Modified from the 2 sample K-S test subroutine in Numerical Recipes in Fortran
!
! Input: 2 cumulative PDFs (data1, data2) and their size (n1,n2)
! Output: d - Kuiper's statistic (measure of max +/- distance between CDF
!   curves)
!   prob - Significance level of the given value of d
!   h - Whether or not the null hypothesis is rejected at significance level
!   alpha
!
! dwhipp - 04/08


      SUBROUTINE kptwo(data1,n1,data2,n2,nsamp,d,prob,h)
      IMPLICIT NONE

! Variable declaration
      REAL*4,INTENT(OUT) :: d,prob
      INTEGER :: n1,n2,n3,n4,h,nsamp,i
      REAL*4 :: en1,en2,en,probkp,alpha,prob2,ns
      REAL*4,DIMENSION(n1),INTENT(IN) :: data1
      REAL*4,DIMENSION(n2),INTENT(IN) :: data2
      REAL*4,DIMENSION(size(data1)+size(data2)) :: dat,org,orgout,orgout2

! Variable initialization/definition
      alpha=0.05

! Fill dat and org arrays
      ns=real(nsamp)
      en1=n1
      en2=n2
      dat(1:n1)=data1
      dat(n1+1:n1+n2)=data2
      org(1:n1)=0.0
      org(n1+1:n1+n2)=1.0
      n3=n1+n2
      n4=n1+n2

! Sort dat and org arrays; calculate array cumulative sums
      call sort2(dat,n3,org,n4)
      call cumsum(org,n4,orgout)
      call cumsum(1.-org,n4,orgout2)

!      open(1111,file='test_cdf.dat',status='unknown')
!      do i=1,size(orgout)
!        write(1111,*) i,' ',orgout(i)/en2
!        !write(1111,*) orgout(i)/en2,' ',i
!      enddo
!      close(1111)
!
!      open(1111,file='test_cdf2.dat',status='unknown')
!      do i=1,size(orgout2)
!        write(1111,*) i,' ',orgout2(i)/en1
!        !write(1111,*) orgout2(i)/en1,' ',i
!      enddo
!      close(1111)

! Determine Kuiper's statistic, significance level and if the null hypothesis
! is rejected at level alpha
      !d=maxval(abs((orgout/en2)-(orgout2/en1)))
      d=maxval((orgout/en2)-(orgout2/en1))+maxval((orgout2/en1)-(orgout/en2))
      !write (*,*) 'd: ',d
      
!       do i=1,n3
!         print *,dat(i)
!       enddo
      !en=sqrt(en1*en2/(en1+en2))
      en=sqrt(ns**2/(2*ns))
      prob=probkp((en+0.155+0.24/en)*d)
      !write (*,*) 'prob: ',prob
!       print *,''
!       print *,'d: ',d
!       print *,'prob: ',prob
      h=0
      if (alpha.ge.prob) h=1
      END SUBROUTINE kptwo

! Subroutine cumsum calculates the cumulative sum of an array of dimensions
! [1 x n]. Modified from the version in Numerical Recipes in Fortran 90. This
! code simply adds the previous array values to the current and stores them
! in a new array of the same size.
!
! Input: Array arr of size n
! Output: ans - Cumulative sum array of arr
!
! dwhipp - 04/08

      SUBROUTINE cumsum(arr,n,ans)
      REAL*4,DIMENSION(n),INTENT(IN) :: arr
      REAL*4,DIMENSION(n),INTENT(OUT) :: ans
      INTEGER :: n,j
      REAL*4 :: sd
      if (n == 0) return
      sd=0.
      ans(1)=arr(1)+sd
      do j=2,n
        ans(j)=ans(j-1)+arr(j)
      end do
      END SUBROUTINE cumsum

! Function probkp calculates the approximate significance level of the given
! value of d. Modified from the version in Numerical Recipes in Fortran 90.
!
! dwhipp - 04/08

      FUNCTION probkp(alam)
      IMPLICIT NONE
      REAL*4, INTENT(IN) :: alam
      REAL*4 :: probkp
      REAL*4,PARAMETER :: EPS1=0.001,EPS2=1.0e-8
      INTEGER,PARAMETER :: NITER=100
      INTEGER :: j
      REAL*4 :: a2,fac,term,kterm,termbf
      a2=-2.0*alam*alam
      fac=2.0
      probkp=0.0
      termbf=0.0
      do j=1,NITER
        kterm=(4*j*j*alam*alam)-1
        term=fac*kterm*exp(a2*j*j)
        probkp=probkp+term
        if (abs(term) <= EPS1*termbf .or. abs(term) <= EPS2*probkp) RETURN
        termbf=abs(term)
      end do
      probkp=1.0
      END FUNCTION probkp