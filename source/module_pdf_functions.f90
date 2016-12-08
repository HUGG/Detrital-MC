      module pdf_functions

      contains
      subroutine get_pdf_size(age,ageu,lc,num,pdfmin,pdfmax,dx,calc_pdf_range)

      USE definitions

      IMPLICIT NONE

      ! Passed in/out variable declaration
      integer(kind=sp) :: lc,num
      real(kind=sp),dimension(:) :: age(lc),ageu(lc)
      real(kind=sp) :: pdfmin,pdfmax,dx
      logical :: calc_pdf_range

      if (calc_pdf_range) then
        !pdfmin=minval(age)-2*maxval(ageu)                                      ! Find range of ages + uncertainties
        pdfmin=minval(age)-10*maxval(ageu)                                      ! Find range of ages + uncertainties
        pdfmax=maxval(age)+10*maxval(ageu)
      endif

      num=nint((pdfmax-pdfmin)/dx)                                              ! Find number of values in PDF arrays

      return
      end subroutine get_pdf_size

      subroutine make_age_pdf(age,ageu,alpha,eratesc,lc,num,n,pdf,pdfmin,dx,   &
                              pdfvsc,pi,cnt,fullpdf)

      USE definitions

      IMPLICIT NONE

      ! Passed in/out variable declaration
      integer(kind=sp) :: lc,num,cnt,fullpdf
      integer(kind=sp), dimension(:) :: eratesc(lc)
      real(kind=sp) :: pdfmin,pdfmax,dx,pdfvsc,pi,alpha
      real(kind=sp),dimension(:) :: age(lc),ageu(lc),n(num+1),pdf(num+1)

      ! Internal subroutine variables
      integer(kind=sp) :: hm,i,j,k,agecnt
      real(kind=sp) :: sum,amin
      real(kind=sp),dimension(:),allocatable :: psum,p

      allocate(psum(num+1),p(num+1))

      do i=1,num+1
        n(i)=pdfmin+real(i-1)*dx                                                  ! Fill age range array
      enddo
      psum=0.
      agecnt=0

      do i=1,lc                                                                 ! Loop over number of sample ages
        if (fullpdf == 1) then
          do j=1,eratesc(i)                                                     ! Loop over scaled erosion rate for given sample
            do k=1,num+1                                                        ! Loop over age range in PDF
              p(k)=(1./(alpha*ageu(i)*sqrt(2.*pi)))*exp(-0.5*((n(k)-age(i))/&   ! Fill probability array
                   (alpha*ageu(i)))**2.)
              psum(k)=psum(k)+p(k)                                              ! Fill sum array to check area under array curve
            enddo
            agecnt=agecnt+1
          enddo
        else
          do k=1,num+1                                                          ! Loop over age range in PDF
            p(k)=(1./(alpha*ageu(i)*sqrt(2.*pi)))*exp(-0.5*((n(k)-age(i))/     &! Fill probability array
                 (alpha*ageu(i)))**2.)
            psum(k)=psum(k)+p(k)                                                ! Fill sum array to check area under array curve
          enddo
          agecnt=agecnt+1
        endif
      enddo

      sum=0.
      do i=1,num+1                                                              ! Loop over age range in PDF
        pdf(i)=(psum(i)/real(agecnt))                                           ! Scale PDF array to normalize area under PDF curve
        sum=sum+pdf(i)*dx                                                       ! Calculate area under curve
      enddo

      ! Generate data PDF vector
      cnt=0
      do i=1,num+1
        hm=nint(pdfvsc*pdf(i))                                                  ! Set number of occurances of given age at current probability
        do j=1,hm
          cnt=cnt+1                                                             ! Count total number of ages in PDF vector for allocation below
        enddo
      enddo

      deallocate(psum,p)

      return
      end subroutine make_age_pdf

      subroutine make_age_pdfv(num,pdfvsc,n,pdf,pdfv,cnt)

      USE definitions

      IMPLICIT NONE

      ! Passed in/out variable declaration
      integer(kind=sp) num,cnt
      real(kind=sp) pdfvsc
      real(kind=sp),dimension(:) :: n(num+1),pdf(num+1),pdfv(cnt)

      ! Internal subroutine variables
      integer(kind=sp) hm,i,j

      cnt=0
      pdfv=0.
      do i=1,num+1
        hm=nint(pdfvsc*pdf(i))
        do j=1,hm
          cnt=cnt+1
          pdfv(cnt)=n(i)                                                        ! Fill PDF vector scaling number of ages by the probability they occur at given age
        enddo
      enddo

      return
      end subroutine make_age_pdfv

      subroutine make_age_cdf(pdf,num,dx,cdf)

      USE definitions

      IMPLICIT NONE

      ! Passed in/out variable declaration
      integer(kind=sp) num
      real(kind=sp) :: dx
      real(kind=sp),intent(in) :: pdf(num+1)
      real(kind=sp),intent(out) :: cdf(num+1)

      ! Internal subroutine variables
      integer(kind=sp) i

      cdf(1)=pdf(1)
      do i=2,num+1
        ! Simple integration using trapezoidal rule
        cdf(i)=cdf(i-1)+(pdf(i)+pdf(i-1))/2.0*dx
      enddo

      return
      end subroutine make_age_cdf

      subroutine make_age_ecdf(age,eratesc,eratesum,n,lc,num,ecdf)

      USE definitions

      IMPLICIT NONE

      ! Passed in/out variable declaration
      integer(kind=sp),intent(in) :: lc,num,eratesc(lc),eratesum
      real(kind=sp),intent(in) :: n(num+1)
      real(kind=sp),intent(inout) :: age(lc)
      real(kind=sp),intent(out) :: ecdf(num+1)

      ! Internal subroutine variables
      integer(kind=sp) i,j,agei,agecnt
      real(kind=sp),allocatable :: erateages(:)

      ! Create age array scaled by erosion rates
      allocate(erateages(eratesum))
      agecnt=0
      do i=1,lc
        do j=1,eratesc(i)
          agecnt=agecnt+1
          erateages(agecnt)=age(i)
        enddo
      enddo

      ! Sort age array (?)
      call insertion_sort(erateages)

      ! Initialize counter variable for position in incoming raw age array
      agei=1

      ! Loop over all ages in PDF age range and check to see if the age in the
      ! raw age array is less than that in the PDF age range. If so, increment
      ! age counter and increase value in ecdf.
      do i=1,num+1
        do while (erateages(agei) <= n(i) .and. agei <= eratesum)
          agei=agei+1
        enddo
        ecdf(i)=real(agei-1)/real(eratesum)
      enddo

      deallocate(erateages)

      return
      end subroutine make_age_ecdf

      function median(a, found)

        USE definitions

        IMPLICIT NONE

        real(kind=sp), dimension(:), intent(in) :: a
        ! the optional found argument can be used to check
        ! if the function returned a valid value; we need this
        ! just if we suspect our "vector" can be "empty"
        logical, optional, intent(out) :: found
        real(kind=sp) :: median

        integer(kind=sp) :: l
        real(kind=sp), dimension(size(a,1)) :: ac

        if ( size(a,1) < 1 ) then
          if ( present(found) ) found = .false.
        else
          ac = a
          ! this is not an intrinsic
          call insertion_sort(ac)

          l = size(a,1)
          if ( mod(l, 2) == 0 ) then
            median = (ac(l/2+1) + ac(l/2))/2.0
          else
            median = ac(l/2+1)
          endif

          if ( present(found) ) found = .true.
        endif

      end function median

      subroutine insertion_sort(a)

        USE definitions

        IMPLICIT NONE

        real(kind=sp), intent(in out), dimension(:) :: a
        real(kind=sp) :: temp
        integer(kind=sp) :: i, j

        do i = 2, size(a)
          j = i - 1
          temp = a(i)
          do while (j>=1 .and. a(j)>temp)
            a(j+1) = a(j)
            j = j - 1
          enddo
          a(j+1) = temp
        enddo
      end subroutine insertion_sort

      function kuiper(alpha,d,nsamp)

        USE definitions

        IMPLICIT NONE

        ! Values passed in/returned
        real(kind=sp) :: alpha,d
        integer(kind=sp) :: nsamp,kuiper
        ! Internal values
        real(kind=sp) :: prob,ns,probkp,en

        ns=real(nsamp)
        en=sqrt(ns**2/(2*ns))
        ! probkp is a function from kptwo.f90
        prob=probkp((en+0.155+0.24/en)*d)
        !write (*,*) 'probnew: ',prob
        kuiper=0
        if (alpha.ge.prob) kuiper=1
        return
      end function kuiper

      end module pdf_functions
