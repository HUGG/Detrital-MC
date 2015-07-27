      module pdf_functions

      contains
      subroutine get_pdf_size(age,ageu,lc,num,pdfmin,pdfmax,dx,calc_pdf_range)

      implicit none

      ! Passed in/out variable declaration
      integer lc,num
      real*4,dimension(:) :: age(lc),ageu(lc)
      real*4 pdfmin,pdfmax,dx
      logical :: calc_pdf_range
      
      if (calc_pdf_range) then
        !pdfmin=minval(age)-2*maxval(ageu)                                      ! Find range of ages + uncertainties
        pdfmin=minval(age)-10*maxval(ageu)                                      ! Find range of ages + uncertainties
        pdfmax=maxval(age)+10*maxval(ageu)
      endif
      
      num=nint((pdfmax-pdfmin)/dx)                                              ! Find number of values in PDF arrays

!        write (*,*) 'dx in pdf functions: ',dx
!        write (*,*) 'num in pdf functions: ',num
!        write (*,*) 'age: ',age
!        write (*,*) 'ageu: ',ageu
!        write (*,*) 'lc: ',lc
!        write (*,*) 'pdfmin: ',pdfmin
!        write (*,*) 'pdfmax: ',pdfmax

      return
      end subroutine get_pdf_size
      
      subroutine make_age_pdf(age,ageu,alpha,eratesc,lc,num,n,pdf,pdfmin,pdfmax,dx,    &
                              pdfvsc,pi,cnt)

      implicit none

      ! Passed in/out variable declaration
      integer lc,num,cnt
      integer, dimension(:) :: eratesc(lc)
      real*4  pdfmin,pdfmax,dx,pdfvsc,pi,alpha
      real*4,dimension(:) :: age(lc),ageu(lc),n(num+1),pdf(num+1)

      ! Internal subroutine variables
      integer hm,i,j,k,agecnt
      real*4  sum,amin
      real*4,dimension(:),allocatable :: psum,p
      
      allocate(psum(num+1),p(num+1))

      !amin=minval(age)-2*maxval(ageu)                                          ! Find range of ages + uncertainties
      !amin=minval(age)-10*maxval(ageu)                                         ! Find range of ages + uncertainties
      !amin=pdfmin

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

      implicit none

      ! Passed in/out variable declaration
      integer num,cnt
      real*4  pdfvsc
      real*4,dimension(:) :: n(num+1),pdf(num+1),pdfv(cnt)

      ! Internal subroutine variables
      integer hm,i,j

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