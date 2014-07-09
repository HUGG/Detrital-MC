      FUNCTION probkuiper(alam)
! Function used in Kuiper test code. Modified from probks.f from Numerical
! Recipes in Fortran
!
! dwhipp - 04/08

      REAL probkuiper,alam,EPS1,EPS2
      PARAMETER (EPS1=0.001, EPS2=1.e-8)
      INTEGER j
      REAL a2,fac,term,termbf
      a2=-2.*alam**2
      fac=2.
      probkuiper=0.
      termbf=0.
      do 11 j=1,100
        term=fac*exp(a2*j**2)
        probkuiper=probkuiper+term
        if (abs(term).le.EPS1*termbf .or. abs(term).le.EPS2*probkuiper) return
        fac=-fac
        termbf=abs(term)
11    continue
      probkuiper=1.
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
