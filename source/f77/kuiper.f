      SUBROUTINE kuiper(data1,n1,data2,n2,d,prob)
! Kuiper test code modified from the 2 sample K-S test from Numerical Recipes
! in Fortran
!
! dwhipp - 04/08

      INTEGER n1,n2
      REAL d,prob,data1(n1),data2(n2)
CU    USES probkuiper,sort
      INTEGER j1,j2
      REAL d1,d2,dt,en1,en2,en,fn1,fn2,probkuiper
      call sort(n1,data1)
      call sort(n2,data2)
      en1=n1
      en2=n2
      j1=1
      j2=1
      fn1=0.
      fn2=0.
      d=0.
1     if(j1.le.n1.and.j2.le.n2)then
        d1=data1(j1)
        d2=data2(j2)
        if(d1.le.d2)then
          fn1=j1/en1
          j1=j1+1
        endif
        if(d2.le.d1)then
          fn2=j2/en2
          j2=j2+1
        endif
        dt=max(fn1 - fn2) + max(fn2 - fn1);
        if(dt.gt.d)d=dt
      goto 1
      endif
      en=sqrt(en1*en2/(en1+en2))
      prob=probkuiper((en+0.155+0.24/en)*d)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
