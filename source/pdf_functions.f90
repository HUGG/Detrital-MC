      module pdf_functions

      contains
      subroutine get_pdf_size(age,ageu,lc,num,pdfmin,pdfmax,dx,calc_pdf_range)

      IMPLICIT NONE

      ! Variable declaration
      integer, parameter :: sp = selected_real_kind(6, 37)

      ! Passed in/out variable declaration
      integer(kind=sp) :: lc,num
      real(kind=sp) :: age(lc),ageu(lc)
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
                              pdfvsc,pi,cnt)

      IMPLICIT NONE

      integer, parameter :: sp = selected_real_kind(6, 37)

      ! Passed in/out variable declaration
      integer(kind=sp) :: lc,num,cnt
      integer(kind=sp) :: eratesc(lc)
      real(kind=sp) :: pdfmin,dx,pdfvsc,pi,alpha
      real(kind=sp) :: age(lc),ageu(lc),n(num+1),pdf(num+1)

      ! Internal subroutine variables
      integer(kind=sp) :: hm,i,j,k,agecnt
      real(kind=sp) ::  sum
      real(kind=sp),allocatable :: psum(:),p(:)
      
      allocate(psum(num+1),p(num+1))

      do i=1,num+1
        n(i)=pdfmin+real(i-1)*dx                                                  ! Fill age range array
      enddo
      psum=0.
      agecnt=0
      do i=1,lc
        do j=1,eratesc(i)
          do k=1,num+1
            p(k)=(1./(alpha*ageu(i)*sqrt(2.*pi)))*exp(-0.5*((n(k)-age(i))/&     ! Fill probability array
                 (alpha*ageu(i)))**2.)
            psum(k)=psum(k)+p(k)                                                ! Fill sum array to check area under array curve
          enddo
          agecnt=agecnt+1
        enddo
      enddo

      sum=0.
      do i=1,num+1
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

      IMPLICIT NONE

      integer, parameter :: sp = selected_real_kind(6, 37)

      ! Passed in/out variable declaration
      integer(kind=sp) :: num,cnt
      real(kind=sp) ::  pdfvsc
      real(kind=sp) :: n(num+1),pdf(num+1),pdfv(cnt)

      ! Internal subroutine variables
      integer(kind=sp) :: hm,i,j

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
      end module pdf_functions