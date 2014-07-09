subroutine make_age_pdf(age,ageu,lc,pdfmin,pdfmax,dx,pdfvsc,pi,num,n,pdf,pdfv,&
                       cnt)

      implicit none

      ! Passed in/out variable declaration
      integer lc,i,j,k,num,cnt
      real*4  pdfmin,pdfmax,dx,pdfvsc,pi
      real*4,dimension(:) :: age(lc),ageu(lc)
      real*4,dimension(:),allocatable,intent(inout) :: n,pdf,pdfv

      ! Internal subroutine variables
      integer hm
      real*4  amin,amax,sum
      real*4,dimension(:),allocatable :: psum,p

      print *,'Entered make_age_pdfs'

      !amin=minval(age)-2*maxval(ageu)                                          ! Find range of ages + uncertainties
      amin=minval(age)-10*maxval(ageu)                                          ! Find range of ages + uncertainties
      !amin=pdfmin
      amax=maxval(age)+10*maxval(ageu)
      !amax=pdfmax
      num=nint((amax-amin)/dx)                                                  ! Find number of values in PDF arrays
      print *,'amax: ',amax
      print *,'amin: ',amin
      print *,'dx: ',dx
      print *,'num: ',num
      print *,'nint((amax-amin)/dx): ',nint((amax-amin)/dx)
      !do j=1,lc                                                                ! Scale uncertainties using optimal scaling factor of 0.6
      !  ageu(j)=ageu(j)*0.6                                                    ! See Brandon, M., Probability Density Plot for Fission-Track Grain-Age Samples,
      !enddo                                                                    ! Radiation Measurements, Vol. 26, No. 5, pp. 663-676, 1996
      print *,'Before allocate'
      allocate(n(num+1))
      print *,'After allocate1'
      allocate(psum(num+1))
      print *,'After allocate2'
      allocate(pdf(num+1))
      print *,'After allocate3'
      allocate(p(num+1))
      print *,'After allocate4'
      allocate(n(num+1),psum(num+1),pdf(num+1),p(num+1))                        ! Allocate data PDF arrays
      print *,'After allocate'
      do j=1,num+1
        n(j)=amin+real(j-1)*dx                                                  ! Fill age range array
      enddo
      psum=0.
      do j=1,lc
        do k=1,num+1
          p(k)=(1./(ageu(j)*sqrt(2.*pi)))*exp(-0.5*((n(k)-age(j))/&             ! Fill probability array
                (ageu(j)))**2.)
          psum(k)=psum(k)+p(k)                                                  ! Fill sum array to check area under array curve
        enddo
      enddo

      sum=0.
      do j=1,num+1
        !pdf(j)=(psum(j)/real(lc))*dx                                           ! Scale PDF array to normalize area under PDF curve
        pdf(j)=(psum(j)/real(lc))                                               ! Scale PDF array to normalize area under PDF curve
        sum=sum+pdf(j)                                                          ! Calculate area under curve
      enddo

      ! Generate data PDF vector
      cnt=0
      do j=1,num+1
        hm=nint(pdfvsc*pdf(j))                                                  ! Set number of occurances of given age at current probability
        do k=1,hm
          cnt=cnt+1                                                             ! Count total number of ages in PDF vector for allocation below
        enddo
      enddo
      allocate(pdfv(cnt))                                                       ! Allocate PDF vector
      cnt=0
      pdfv=0.
      do j=1,num+1
        hm=nint(pdfvsc*pdf(j))
        do k=1,hm
          cnt=cnt+1
          pdfv(cnt)=n(j)                                                        ! Fill PDF vector scaling number of ages by the probability they occur at given age
        enddo
      enddo
      
      deallocate(psum)

return
end