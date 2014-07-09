! detrital_mc.f90
!
! This fortran code uses a monte carlo simulation method to randomly grab n
! samples from a distribution of model predicted cooling ages, assign a
! designated uncertainty to the samples, generate a PDF of the model predicted
! cooling ages, compare that PDF to an observed age PDF using a Kuiper test,
! and finally record the results of that test. This process is repeated a
! large number of times (~10000).
!
! dwhipp - 04/08
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program detrital_mc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none
! Variable declaration
      real*4,dimension(:),allocatable :: oage,oageu,on,opsum,op,opdf,opdfv
      real*4,dimension(:),allocatable :: page,pageu,perate,pagesc,pageusc
      real*4,dimension(:),allocatable :: pagemc,pageumc,pnmc,ppsummc,ppmc,ppdfmc
      real*4,dimension(:),allocatable :: ppdfvmc,kpct,lsage,lsageu,lsagesc
      real*4,dimension(:),allocatable :: lsageusc,pagesct,pageusct,lserate
      real*4,dimension(:),allocatable :: pagetot,pageutot,pn,ppsum,ppdf,pp,ppdfv
      real*4,dimension(:),allocatable :: pagemc2,pageumc2,pnmc2,ppsummc2,ppmc2
      real*4,dimension(:),allocatable :: ppdfmc2,ppdfvmc2
      real*4,dimension(:,:),allocatable :: pmc
      integer,dimension(:),allocatable :: peratesc,kuiper_res,lseratesc,eratetot
      !integer,dimension(38) :: numsamp
      !integer,dimension(21) :: numsamp
      integer,dimension(3) :: numsamp
      real*4 :: d,prob,pagemu,pagemed,pagesd,pdfvsc,mc_iterf,jf,pageus,pageup
      real*4 :: dx,osum,agenow,psummc,pi,peratemin,peratescl,psum,d1,d2,d3,d4,d5
      real*4 :: d6,d7,psummc2
      real*8 :: randflt
      integer :: olc,onum,oamin,oamax,h,i,j,k,l,mc_iter,basnum,pamin,pamax,pnum
      integer :: plc,plcsc,cnt,paminmc,pamaxmc,pnummc,cnt2,cnt3,hm,hm2,cnt4,m
      integer :: lsc,cnt5,cnt6,lscsc,lctot,pdfnum,pdfmin,pdfmax,mcsamp,cnt7,nsc
      integer :: paminmc2,pamaxmc2,pnummc2,cnt8
      integer*8 :: eratesum,eratechk,rint
      logical lsero,datacomp,mcboth
      character :: buffer*8,obasin*12,pbasin*12,jc*3,hc*5,dump*80,mcschar*5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! User-defined variables
      datacomp=.false.                                                          ! Compare model output to data?
      lsero=.true.                                                              ! Simulate landslide erosion rate scaling?
      mcboth=.false.                                                            ! Use MC model to compare 2 age PDFs?
      !numsamp=0
      !numsamp=(/1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,& ! Number of samples to use if only simulating model-predicted ages
      !         105,110,115,120,125,130,135,140,145,150,200,250,300,350,400,450,&
      !         500/)
      !numsamp=(/1,5,10,15,20,25,30,40,50,60,70,80,90,100,125,150,200,250,300,&
      !         400,500/)
      numsamp=(/10,100,500/)
      basnum=16                                                                 ! Number of basins to analyze
      !mc_iter=10000                                                             ! Number of iterations in the Monte Carlo simulation
      mc_iter=1                                                             ! Number of iterations in the Monte Carlo simulation
      pageus=5.                                                                 ! Percent age uncertainty if not comparing to data
      dx=0.01                                                                   ! Specify x spacing for data PDF generation
      pdfmin=0.                                                                 ! Minimum age for PDF calculation
      pdfmax=15.                                                                ! Maximum age for PDF calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

! Initialize random number generator
      call init_random_seed()
      !call random_seed()

! Variable initialization/declaration
      pi=atan(1.)*4.                                                            ! Define pi
      pdfvsc=50.                                                                ! Approximate number of values in scaled PDFs

! Allocate kuiper test results summary array
      allocate(kpct(basnum))

! Read in basin info
      !open(10,file='basin_summary_info.txt',status='old')
      open(10,file='new_basin_summary_info.txt',status='old')

      do i=1,basnum                                                             ! Loop through number of subbasins to analyze
        plc=0                                                                   ! Reset predicted age line counter
        lsc=0																	! Reset landslide age line counter
        read(10,*) obasin,pbasin,olc,pagemu,pagemed,pagesd                      ! Read input basin filenames, number of samples, mean 1 sigma and associated s.d.
        
        !if (i.eq.5 .or. i.eq.15 .or. i.eq.16) then
        !if (i.eq.5) then
        if (i.eq.12) then
          print *,'Processing basin ',i,' of ',basnum

          open(21,file='kmc_percent_pass_summary_'//obasin,status='unknown')      ! Open basin summary results file
          write(21,'(a27)') '            Percent passing'                         ! Write header
          write(21,'(a23)') 'Basin       Kuiper test'

          if (datacomp) then                                                    ! If comparing to data, read in observed cooling ages and generate PDF
            print *,'Data comparison requested, reading observed ages...'
            open (11,file='observed_det_ages/'//obasin,status='old')
            allocate(oage(olc),oageu(olc))                                      ! Allocate observed age and age uncertainty arrays

            do j=1,olc
              read(11,*) oage(j),oageu(j)                                       ! Fill observed sample arrays
            enddo
            close(11)
          endif

          print *,'Reading predicted ages...'
          open (12,file='detrital_ages/det_ages_'//pbasin,status='old')
    3     read (12,*,end=4) dump
          plc=plc+1                                                             ! Find number of lines in model age files
          goto 3
    4     rewind (12)
          allocate(page(plc),pageu(plc),perate(plc),peratesc(plc))              ! Allocate model age/uncertainty and erate arrays
          do j=1,plc
            read(12,*) d1,d2,d3,d4,d5,d6,d7,page(j)                             ! Read in model predicted ages
            if (datacomp) then                                                  ! Assign predicted age uncertainties based on dataset, if doing a data comparison
              !pageup=pagemu                                                    ! Assign uncertainty (using mean measured value for dataset)
              pageup=pagemed                                                    ! Assign uncertainty (using median measured value for dataset)
              !pageu(j)=page(i)*0.022                                           ! Include standard deviation of errors (pagesd is the uncertainty s.d.)
              pageu(j)=page(j)*(pageup/100.)                                    ! Assign uncertainty to predicted ages
            else
              pageu(j)=page(j)*(pageus/100.)                                    ! Assign predicted age uncertainties based on user-specified value above
            endif
          enddo
          close(12)

          print *,'Reading predicted erosion rates...'
          open (13,file='detrital_ages/erates_'//pbasin,status='old')
          do j=1,plc
            read(13,*) d1,d2,d3,perate(j),dump                                  ! Read in erosion rates for scaling age distribution
            !peratesc(j)=nint(perate(j)*100.)                                   ! Generate scaling factors by multiplying erosion rates by 100 and converting to integers
          enddo
          close(13)
          peratemin=minval(perate)                                              ! Determine minimum erosion rate in model domain
          peratescl=1./peratemin                                                ! Set scaling value to ensure at least 1 occurance of min rate ages
          peratesc=nint(perate*peratescl)                                       ! Generate scaling factors by multiplying erosion rates by peratescl and converting to integers

          if (lsero) then                                                       ! If simulating landsliding, then read in the landslide ages and erosion rates
            print *,'Landslide simulation requested, reading landslide ',&
                    'erosion rates...'
            open(14,file='detrital_ages/ls_ages_'//pbasin,status='old')
    5       read(14,*,end=6) dump
            lsc=lsc+1
            goto 5
    6       rewind (14)
            allocate(lsage(lsc),lsageu(lsc),lserate(lsc))
            do j=1,lsc
              read(14,*) lsage(j),lserate(j)
              if (datacomp) then                                                ! If doing a data comparison, then use the data uncertainty specified above
                lsageu(j)=lsage(j)*(pageup/100.)
              else                                                              ! Otherwise, use the user-specified value above
                lsageu(j)=lsage(j)*(pageus/100.)
              endif
            enddo
            close(14)
            allocate(lseratesc(lsc))
            lseratesc=nint(lserate*peratescl)                                   ! Scale landslide erosion rates
          endif


! Generate data PDF
          if (datacomp) then                                                    ! Generate data PDF if doing a data comparison
            print *,'Generating observed age PDFs...'
            !oamin=nint(minval(oage)-2*maxval(oageu))                           ! Find range of ages + uncertainties
            oamin=nint(minval(oage)-10*maxval(oageu))                           ! Find range of ages + uncertainties
            !oamin=pdfmin
            oamax=nint(maxval(oage)+10*maxval(oageu))
            !oamax=pdfmax
            onum=(oamax-oamin)/dx                                               ! Find number of values in PDF arrays
            !do j=1,olc                                                         ! Scale uncertainties using optimal scaling factor of 0.6
            !  oageu(j)=oageu(j)*0.6                                            ! See Brandon, M., Probability Density Plot for Fission-Track Grain-Age Samples,
            !enddo                                                              ! Radiation Measurements, Vol. 26, No. 5, pp. 663-676, 1996
            allocate(on(onum+1),opsum(onum+1),opdf(onum+1),op(onum+1))          ! Allocate data PDF arrays
            do j=1,onum+1
              on(j)=oamin+(j-1)*dx                                              ! Fill age range array
            enddo
            opsum=0.
            do j=1,olc
              do k=1,onum+1
                op(k)=(1./(oageu(j)*sqrt(2.*pi)))*exp(-0.5*((on(k)-oage(j))/&   ! Fill probability array
                      (oageu(j)))**2.)
                opsum(k)=opsum(k)+op(k)                                         ! Fill sum array to check area under array curve
              enddo
            enddo

            osum=0.
            do j=1,onum+1
              !opdf(j)=(opsum(j)/olc)*dx                                        ! Scale PDF array to normalize area under PDF curve
              opdf(j)=(opsum(j)/olc)                                            ! Scale PDF array to normalize area under PDF curve
              osum=osum+opdf(j)                                                 ! Calculate area under curve
            enddo

            ! Generate data PDF vector
            cnt2=0
            do j=1,onum+1
              hm=nint(pdfvsc*opdf(j))                                           ! Set number of occurances of given age at current probability
              do k=1,hm
                cnt2=cnt2+1                                                     ! Count total number of ages in PDF vector for allocation below
              enddo
            enddo
            allocate(opdfv(cnt2))                                               ! Allocate PDF vector
            cnt2=0
            opdfv=0.
            do j=1,onum+1
              hm=nint(pdfvsc*opdf(j))
              do k=1,hm
                cnt2=cnt2+1
                opdfv(cnt2)=on(j)                                               ! Fill PDF vector scaling number of ages by the probability they occur at given age
              enddo
            enddo
          endif


! Generate combined kinematic/landslide age/frequency distributions
          lctot=plc+lsc
          eratesum=0
          allocate(pagetot(lctot),pageutot(lctot),eratetot(lctot))
          do j=1,plc
            pagetot(j)=page(j)
            pageutot(j)=pageu(j)
            eratetot(j)=peratesc(j)
            eratesum=eratesum+eratetot(j)
          enddo
          do j=1,lsc
            pagetot(plc+j)=lsage(j)
            pageutot(plc+j)=lsageu(j)
            eratetot(plc+j)=lseratesc(j)
            eratesum=eratesum+eratetot(j)
          enddo


! Generate full predicted age PDF (if not comparing to data)
!
! NOTE: The full predicted age PDF does not include landslides (can uncomment
! the commented lines to do that).
          if (.not.datacomp) then                                               ! Generate predicted age PDF if not doing a data comparison
            if (.not.mcboth) then
              print *,'Generating full predicted age PDFs...'
              !pamin=nint(minval(page)-2*maxval(pageu))                         ! Find range of ages + uncertainties
              pamin=nint(minval(page)-10*maxval(pageu))                         ! Find range of ages + uncertainties
              !pamin=pdfmin
              pamax=nint(maxval(page)+10*maxval(pageu))
              !pamax=pdfmax
              pnum=(pamax-pamin)/dx                                             ! Find number of values in PDF arrays
              !do j=1,plc                                                       ! Scale uncertainties using optimal scaling factor of 0.6
              !  pageu(j)=pageu(j)*0.6                                          ! See Brandon, M., Probability Density Plot for Fission-Track Grain-Age Samples,
              !enddo                                                            ! Radiation Measurements, Vol. 26, No. 5, pp. 663-676, 1996
              allocate(pn(pnum+1),ppsum(pnum+1),ppdf(pnum+1),pp(pnum+1))        ! Allocate predicted age PDF arrays
              do j=1,pnum+1
                pn(j)=pamin+(j-1)*dx                                            ! Fill age range array
              enddo
              ppsum=0.
              !do j=1,lctot
              do j=1,plc
                !do k=1,eratetot(j)+1                                           ! Oops!
                do k=1,eratetot(j)
                  do l=1,pnum+1
                    pp(l)=(1./(pageu(j)*sqrt(2.*pi)))*exp(-0.5*((pn(l)-&        ! Fill probability array
                          page(j))/(pageu(j)))**2.)
                    ppsum(l)=ppsum(l)+pp(l)                                     ! Fill sum array to check area under array curve
                  enddo
                enddo
              enddo

              psum=0.
              do j=1,pnum+1
                !ppdf(j)=(ppsum(j)/lctot)                                        ! Scale PDF array to normalize area under PDF curve
                ppdf(j)=(ppsum(j)/plc)                                          ! Scale PDF array to normalize area under PDF curve
                psum=psum+ppdf(j)                                               ! Calculate area under curve
              enddo

              ! Generate data PDF vector
              cnt7=0
              do j=1,pnum+1
                hm=nint(pdfvsc*ppdf(j))                                         ! Set number of occurances of given age at current probability
                do k=1,hm
                  cnt7=cnt7+1                                                   ! Count total number of ages in PDF vector for allocation below
                enddo
              enddo
              allocate(ppdfv(cnt7))                                             ! Allocate PDF vector
              cnt7=0
              ppdfv=0.
              do j=1,pnum+1
                hm=nint(pdfvsc*ppdf(j))
                do k=1,hm
                  cnt7=cnt7+1
                  ppdfv(cnt7)=pn(j)                                             ! Fill PDF vector scaling number of ages by the probability they occur at given age
                enddo
              enddo
            endif
          endif


!!!
! Start monte carlo runs for testing model/data fits with Kuiper test
!!!!
          print *,'Starting Monte Carlo PDF generation...'
          allocate(kuiper_res(mc_iter))                                         ! Allocate kuiper test array
          !allocate(pmc(int(pdfvsc),2*mc_iter))
          if (datacomp) then
            nsc=1
          else
            nsc=size(numsamp)
          endif
          do m=1,nsc                                                            ! Loop through number of desired samples
            !pmc=0.
            cnt4=0
            if (datacomp) then
              mcsamp=olc
            else
              mcsamp=numsamp(m)
            endif
            print *,'Running Monte Carlo simulation for ',mcsamp,'samples'
            allocate(pagemc(mcsamp),pageumc(mcsamp))
            if (mcboth) allocate(pagemc2(mcsamp),pageumc2(mcsamp))
            write (mcschar,'(i5)') mcsamp
            mcschar=adjustl(mcschar)
            do j=1,mc_iter                                                      ! Model age distribution; This loop runs mc_iter times (usually ~10000)
              !call init_random_seed()
              jf=real(j)

! Randomly grab mcsamp (n) grains from model distribution
              do k=1,mcsamp
                call random_number(randflt)                                     ! Generate random number [0,1)
                rint=int8(randflt*(eratesum))+1                                 ! Get random integer value within range of size of scaled age dist.
                eratechk=rint
                cnt=0
                do while (eratechk.gt.0)
                  cnt=cnt+1
                  eratechk=eratechk-eratetot(cnt)
                enddo
                pagemc(k)=pagetot(cnt)                                          ! Add random age to monte carlo age array
                pageumc(k)=pageutot(cnt)                                        ! Add associated uncertainty to monte carlo uncertainty array
              enddo

! Generate MC model PDF using mcsamp grains
              !paminmc=nint(minval(pagemc)-2*maxval(pageumc))                   ! Find range of ages + uncertainties
              paminmc=nint(minval(pagemc)-10*maxval(pageumc))                   ! Find range of ages + uncertainties
              !paminmc=pdfmin
              pamaxmc=nint(maxval(pagemc)+10*maxval(pageumc))
              !pamaxmc=pdfmax
              pnummc=(pamaxmc-paminmc)/dx                                       ! Find number of values in MC model PDF arrays
              !do k=1,olc                                                       ! Scale uncertainties using optimal scaling factor of 0.6
              !  pageumc(k)=pageumc(k)*0.6                                      ! See Brandon, M., Probability Density Plot for Fission-Track Grain-Age Samples,
              !enddo                                                            ! Radiation Measurements, Vol. 26, No. 5, pp. 663-676, 1996
              allocate(pnmc(pnummc+1),ppsummc(pnummc+1),ppdfmc(pnummc+1))
              allocate(ppmc(pnummc+1))
              do k=1,pnummc+1
                pnmc(k)=paminmc+(k-1)*dx                                        ! Fill age range array
              enddo
              ppsummc=0.
              do k=1,mcsamp
                do l=1,pnummc+1
                  ppmc(l)=(1./(pageumc(k)*sqrt(2.*pi)))*exp(-0.5*((pnmc(l)-&    ! Fill probability array
                          pagemc(k))/(pageumc(k)))**2.)
                  ppsummc(l)=ppsummc(l)+ppmc(l)                                 ! Fill sum array to check area under array curve
                enddo
              enddo

              psummc=0.
              do k=1,pnummc+1
                ppdfmc(k)=(ppsummc(k)/mcsamp)                                   ! Scale PDF array to normalize area under PDF curve
                psummc=psummc+ppdfmc(k)                                         ! Calculate area under curve
              enddo

! Generate MC model PDF vector
              cnt3=0
              do k=1,pnummc+1
                hm2=nint(pdfvsc*ppdfmc(k))                                      ! Set number of occurances of given age at current probability
                do l=1,hm2
                  cnt3=cnt3+1                                                   ! Count total number of ages in PDF vector for allocation below
                enddo
              enddo
              allocate(ppdfvmc(cnt3))                                           ! Allocate PDF vector
              cnt3=0
              ppdfvmc=0.
              do k=1,pnummc+1
                hm2=nint(pdfvsc*ppdfmc(k))
                do l=1,hm2
                  cnt3=cnt3+1
                  ppdfvmc(cnt3)=pnmc(k)                                         ! Fill PDF vector scaling number of ages by the probability they occur at given age
                enddo
              enddo

! If a comparison of 2 MC PDFs is desired, generate second PDF
!
! NOTE: This does not include the landslides in the second PDF!
! Uncomment the lines below if that is desired
              if (mcboth) then
                ! Randomly grab mcsamp (n) grains from model distribution
                do k=1,mcsamp
                  call random_number(randflt)                                   ! Generate random number [0,1)
                  !rint=int8(randflt*(eratesum))+1                               ! Get random integer value within range of size of scaled age dist.
                  rint=int8(randflt*(sum(peratesc)))+1                          ! Get random integer value within range of size of scaled age dist.
                  eratechk=rint
                  cnt=0
                  do while (eratechk.gt.0)
                    cnt=cnt+1
                    eratechk=eratechk-eratetot(cnt)
                  enddo
                  pagemc2(k)=pagetot(cnt)                                       ! Add random age to monte carlo age array
                  pageumc2(k)=pageutot(cnt)                                     ! Add associated uncertainty to monte carlo uncertainty array
                enddo

                ! Generate MC model PDF using mcsamp grains
                !paminmc2=nint(minval(pagemc2)-2*maxval(pageumc2))              ! Find range of ages + uncertainties
                paminmc2=nint(minval(pagemc2)-10*maxval(pageumc2))              ! Find range of ages + uncertainties
                !paminmc=pdfmin
                pamaxmc2=nint(maxval(pagemc2)+10*maxval(pageumc2))
                !pamaxmc=pdfmax
                pnummc2=(pamaxmc2-paminmc2)/dx                                  ! Find number of values in MC model PDF arrays
                !do k=1,olc                                                     ! Scale uncertainties using optimal scaling factor of 0.6
                !  pageumc2(k)=pageumc2(k)*0.6                                  ! See Brandon, M., Probability Density Plot for Fission-Track Grain-Age Samples,
                !enddo                                                          ! Radiation Measurements, Vol. 26, No. 5, pp. 663-676, 1996
                allocate(pnmc2(pnummc2+1),ppsummc2(pnummc2+1))
                allocate(ppdfmc2(pnummc2+1),ppmc2(pnummc2+1))
                do k=1,pnummc2+1
                  pnmc2(k)=paminmc2+(k-1)*dx                                    ! Fill age range array
                enddo
                ppsummc2=0.
                do k=1,mcsamp
                  do l=1,pnummc2+1
                    ppmc2(l)=(1./(pageumc2(k)*sqrt(2.*pi)))*exp(-0.5*&          ! Fill probability array
                    ((pnmc2(l)-pagemc2(k))/(pageumc2(k)))**2.)
                    ppsummc2(l)=ppsummc2(l)+ppmc2(l)                            ! Fill sum array to check area under array curve
                  enddo
                enddo

                psummc2=0.
                do k=1,pnummc2+1
                  ppdfmc2(k)=(ppsummc2(k)/mcsamp)                               ! Scale PDF array to normalize area under PDF curve
                  psummc2=psummc2+ppdfmc2(k)                                    ! Calculate area under curve
                enddo

                ! Generate MC model PDF vector
                cnt8=0
                do k=1,pnummc2+1
                  hm2=nint(pdfvsc*ppdfmc2(k))                                   ! Set number of occurances of given age at current probability
                  do l=1,hm2
                    cnt8=cnt8+1                                                 ! Count total number of ages in PDF vector for allocation below
                  enddo
                enddo
                allocate(ppdfvmc2(cnt8))                                        ! Allocate PDF vector
                cnt8=0
                ppdfvmc2=0.
                do k=1,pnummc2+1
                  hm2=nint(pdfvsc*ppdfmc2(k))
                  do l=1,hm2
                    cnt8=cnt8+1
                    ppdfvmc2(cnt8)=pnmc2(k)                                     ! Fill PDF vector scaling number of ages by the probability they occur at given age
                  enddo
                enddo
              endif


! Run Kuiper test to get misfit between data and model
              if (datacomp) then
                call kptwo(opdfv,cnt2,ppdfvmc,cnt3,mcsamp,d,prob,h)
              else
                if (mcboth) then
                  call kptwo(ppdfvmc,cnt3,ppdfvmc2,cnt8,mcsamp,d,prob,h)
                else
                  call kptwo(ppdfv,cnt7,ppdfvmc,cnt3,mcsamp,d,prob,h)
                endif
              endif
              kuiper_res(j)=h                                                   ! Store kuiper test result (0=pass;1=fail) for this iteration in kuiper results array
              if (h.eq.0) cnt4=cnt4+1                                           ! Increment counter for number of models that pass Kuiper test

! Write out data PDF
              if (j.eq.1 .and. m.eq.1) then
                if (datacomp) then
                  open(22,file='age_pdf_output/data_age_PDF_'//obasin,&
                      status='unknown')
                  do k=1,onum+1
                    write(22,*) on(k),opdf(k)
                  enddo
                  close(22)
                else
                  if (.not.mcboth) then
                    open(23,file='age_pdf_output/full_predicted_age_PDF_'&
                        //obasin,status='unknown')
                    do k=1,pnum+1
                      write(23,*) pn(k),ppdf(k)
                    enddo
                    close(23)
                  endif
                endif
              endif

! Write out 100 PDFs that pass the Kuiper test
              !if (h.eq.0 .and. cnt4.le.100) then
              if (cnt4.le.1) then
                write(hc,'(i5)') j
                if (j.lt.10) hc(1:4)='0000'
                if (j.lt.100) hc(1:3)='000'
                if (j.lt.1000) hc(1:2)='00'
                if (j.lt.10000) hc(1:1)='0'
                open(24,file='age_pdf_output/pass_mc_age_PDF_'//hc//'_'&
                    //trim(mcschar)//'_samples_'//obasin,status='unknown')
                do k=1,pnummc+1
                  write(24,*) pnmc(k),ppdfmc(k)
                enddo
                close(24)
                if (mcboth) then
                  open(27,file='age_pdf_output/pass_mc2_age_PDF_'//hc//'_'&
                      //trim(mcschar)//'_samples_'//obasin,status='unknown')
                  do k=1,pnummc2+1
                    write(27,*) pnmc2(k),ppdfmc2(k)
                  enddo
                  close(27)
                endif
              endif


! Write out first 500 monte carlo PDFs
            !if (j.le.500) then
            !  write(jc,'(i3)') j
            !  if (j.lt.10) jc(1:2)='00'
            !  if (j.lt.100) jc(1:1)='0'
            !  open(25,file='age_pdf_output/mc_age_PDF_'//jc//'_'//obasin,&
            !      status='unknown')
            !  do k=1,pnummc+1
            !    write(25,*) pnmc(k),ppdfmc(k)
            !  enddo
            !  close(25)
            !endif

! Deallocate arrays
              deallocate(pnmc,ppsummc,ppdfmc,ppmc,ppdfvmc)                      ! Deallocate arrays reallocated during monte carlo sim
              if (mcboth) then
                deallocate(pnmc2,ppsummc2,ppdfmc2,ppmc2,ppdfvmc2)               ! Deallocate arrays reallocated during monte carlo sim
              endif
            enddo

            mc_iterf=real(mc_iter)
            kpct(i)=(1-(sum(kuiper_res)/mc_iterf))*100.                         ! Store percent of models that passed kuiper test for given basin

! Write output files
            open(20,file='kuiper_mc_results_'//trim(mcschar)//'_samples_'&
                 //obasin,status='unknown')
            do j=1,mc_iter
              write(20,*) kuiper_res(j)                                         ! Write out individual subset kuiper test result values to file for this basin
            enddo
            close(20)

            write(21,'(a12,a5,f10.2)') obasin,mcschar,kpct(i)                    ! Write out the summary percent of models that passed the Kuiper test

            !open(26,file='age_pdf_output/all_PDF_results_'//obasin,status='unknown')
            !write(26,*) pmc                                                    ! Write out individual subset kuiper test result values to file for this basin
            !close(26)

            ! Deallocate arrays
            deallocate(pagemc,pageumc)
            if (mcboth) deallocate(pagemc2,pageumc2)

! End of main loop
          enddo

! Deallocate arrays
          if (datacomp) deallocate(oage,oageu,on,opdf,op,opsum,opdfv)
          if (.not.mcboth) deallocate(pn,ppsum,ppdf,pp,ppdfv)
          deallocate(page,pageu)
          deallocate(perate,peratesc)
          if (lsero) deallocate(lsage,lsageu,lserate,lseratesc)
          deallocate(pagetot,pageutot,eratetot)
          deallocate(kuiper_res)

! Close open files
          close(21)
        endif
      enddo

! Exit
      end
