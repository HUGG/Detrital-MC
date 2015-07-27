      module pdf_functions

      contains
      subroutine get_pdf_size(age,ageu,lc,num,pdfmin,pdfmax,dx,calc_pdf_range)

      implicit none

      ! Passed in/out variable declaration
      integer :: lc,num
      real*4,dimension(:) :: age(lc),ageu(lc)
      real*4 :: pdfmin,pdfmax,dx
      logical :: calc_pdf_range
      
      if (calc_pdf_range) then
        !pdfmin=minval(age)-2*maxval(ageu)                                      ! Find range of ages + uncertainties
        pdfmin=minval(age)-10*maxval(ageu)                                      ! Find range of ages + uncertainties
        pdfmax=maxval(age)+10*maxval(ageu)
      endif
      
      num=nint((pdfmax-pdfmin)/dx)                                              ! Find number of values in PDF arrays

      return
      end subroutine get_pdf_size
      
      subroutine make_age_pdf(age,ageu,alpha,eratesc,lc,num,n,pdf,pdfmin,pdfmax,dx,    &
                              pdfvsc,pi,cnt)

      implicit none

      ! Passed in/out variable declaration
      integer :: lc,num,cnt
      integer, dimension(:) :: eratesc(lc)
      real*4  :: pdfmin,pdfmax,dx,pdfvsc,pi,alpha
      real*4,dimension(:) :: age(lc),ageu(lc),n(num+1),pdf(num+1)

      ! Internal subroutine variables
      integer :: hm,i,j,k,agecnt
      real*4  :: sum,amin
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
        !pdf(i)=(psum(i)/real(lc))*dx                                           ! Scale PDF array to normalize area under PDF curve
        !pdf(i)=(psum(i)/real(lc))                                              ! Scale PDF array to normalize area under PDF curve
        pdf(i)=(psum(i)/real(agecnt))                                           ! Scale PDF array to normalize area under PDF curve
        sum=sum+pdf(i)*dx                                                       ! Calculate area under curve
      enddo
      !write (*,*) 'sum: ',sum

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

      subroutine make_age_cdf(pdf,num,dx,cdf)

      implicit none
      
      ! Passed in/out variable declaration
      integer num
      real*4 :: dx
      real*4,intent(in) :: pdf(num+1)
      real*4,intent(out) :: cdf(num+1)

      ! Internal subroutine variables
      integer i

      cdf(1)=pdf(1)
      do i=2,num+1
        ! Simple integration using trapezoidal rule
        cdf(i)=cdf(i-1)+(pdf(i)+pdf(i-1))/2.0*dx
      enddo
      
      return
      end subroutine make_age_cdf

      subroutine make_age_ecdf(age,eratesc,eratesum,n,lc,num,ecdf)

      implicit none

      ! Passed in/out variable declaration
      integer,intent(in) :: lc,num,eratesc(lc),eratesum
      real*4,intent(in) :: n(num+1)
      real*4,intent(inout) :: age(lc)
      real*4,intent(out) :: ecdf(num+1)

      ! Internal subroutine variables
      integer i,j,agei,agecnt
      real*4,allocatable :: erateages(:)

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
        real, dimension(:), intent(in) :: a
        ! the optional found argument can be used to check
        ! if the function returned a valid value; we need this
        ! just if we suspect our "vector" can be "empty"
        logical, optional, intent(out) :: found
        real :: median
 
        integer :: l
        real, dimension(size(a,1)) :: ac
 
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
        real, intent(in out), dimension(:) :: a
        real :: temp
        integer :: i, j
 
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
        implicit none
        ! Values passed in/returned
        real*4 :: alpha,d
        integer :: nsamp,kuiper
        ! Internal values
        real*4 :: prob,ns,probkp
        integer :: en,h

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