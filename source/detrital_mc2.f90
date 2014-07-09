! detrital_mc2.f90
!
! This fortran code uses a monte carlo simulation method to randomly grab n
! samples from a distribution of model predicted cooling ages, assign a
! designated uncertainty to the samples, generate a PDF of the model predicted
! cooling ages, compare that PDF to an observed age PDF using a Kuiper test,
! and finally record the results of that test. This process is repeated a
! large number of times (~10000).
!
! dwhipp - 10.13
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program detrital_mc2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use definitions
      use pdf_functions

      implicit none

      type (detrital_params) params

      ! Variable declaration
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: qp = selected_real_kind(33, 4931)

      real(kind=sp), allocatable :: oage(:),oageu(:),on(:),opsum(:),op(:),opdf(:),opdfv(:)
      real(kind=sp), allocatable :: page(:),pageu(:),perate(:),pagesc(:),pageusc(:)
      real(kind=sp), allocatable :: pagemc(:),pageumc(:),pnmc(:),ppsummc(:),ppmc(:),ppdfmc(:)
      real(kind=sp), allocatable :: ppdfvmc(:),kpct(:),lsage(:),lsageu(:),lsagesc(:)
      real(kind=sp), allocatable :: lsageusc(:),pagesct(:),pageusct(:),lserate(:)
      real(kind=sp), allocatable :: pagetot(:),pageutot(:),pn(:),ppsum(:),ppdf(:),pp(:),ppdfv(:)
      real(kind=sp), allocatable :: pagemc2(:),pageumc2(:),pnmc2(:),ppsummc2(:),ppmc2(:)
      real(kind=sp), allocatable :: ppdfmc2,ppdfvmc2(:)
      real(kind=sp), allocatable :: pmc(:)
      integer,dimension(:),allocatable :: peratesc,kuiper_res,lseratesc
      !integer*8, allocatable :: lseratesc(:)
      real*4 :: d,prob,pagemu,pagemed,pagesd,mc_iterf,jf,pageup
      real*4 :: osum,agenow,psummc,pi,peratemin,peratescl,d1,d2,d3,d4,d5
      real*4 :: d6,d7,psummc2
      real*8 :: randflt
      integer :: olc,onum,h,i,j,k,pamin,pamax,pnum
      integer :: plc,plcsc,cnt,paminmc,pamaxmc,pnummc,cnt2,cnt3,hm,hm2,cnt4,m
      integer :: lsc,cnt5,cnt6,lscsc,lctot,pdfnum,mcsamp,cnt7,nsc,ocnt,pcnt
      integer :: paminmc2,pamaxmc2,pnummc2,cnt8,curbasin,pcntmc
      integer :: eratesum,peratesum,lseratesum,rint,eratechk
      !integer*8 :: lseratesum,rint,eratechk
      logical datacomp
      !logical lsppdf,lspdf_out
      character(len=5)  :: hc,jc,mcschar
      character(len=8)  :: buffer
      character(len=10) :: onumc,pnumc,pnummcc
      character(len=12) :: obasin,pbasin
      character(len=80) :: dump

      mcboth=.false.

! Initialize random number generator
      call init_random_seed()

! Variable initialization/declaration
      pi=atan(1.)*4.

! Read in the input file
      call read_input_file(params)

! Allocate kuiper test results summary array
      allocate(kpct(params%num_basins))

! Read in basin info
      open(10,file=trim(params%basin_summary_file),status='old')

      do i=1,params%num_basins                                                             ! Loop through number of subbasins to analyze
        plc=0                                                                   ! Reset predicted age line counter
        lsc=0                                                                   ! Reset landslide age line counter
        read(10,*) obasin,pbasin,olc,pagemu,pagemed,pagesd                      ! Read input basin filenames, number of samples, mean 1 sigma and associated s.d.
        
        !if (i.eq.1 .or. i.eq.3 .or. i.eq.10 .or. i.eq.13) then
        if (i.eq.12) then
        !if (i.eq.curbasin) then
        !if (i.gt.0) then
          write (*,'(a,i3,a,i3)') 'Processing basin ',i,' of ',params%num_basins

          if (params%datamcpdfs .or. params%ppdfmcpdfs) then
            open(21,file='kmc_percent_pass_summary_'//trim(obasin)//'.dat',&    ! Open basin summary results file
                 status='unknown')
            write(21,'(a27)') '            Percent passing'                     ! Write header
            write(21,'(a23)') 'Basin       Kuiper test'
          endif

          if (params%datapdf .or. params%datappdf) then                                         ! If using data, read in observed cooling ages and generate PDF
            write (*,'(a)') 'Data comparison requested, reading observed ages...'
            open (11,file='observed_det_ages/'//trim(obasin)//'.dat',&
                  status='old')
            allocate(oage(olc),oageu(olc))                                      ! Allocate observed age and age uncertainty arrays
            do j=1,olc
              read(11,*) oage(j),oageu(j)                                       ! Fill observed sample arrays
              if (j==olc) write (*,'(a,i3,a,a)') 'Read ',olc,' ages for basin ',trim(obasin)
            enddo
            close(11)
          endif

          write (*,'(a)') 'Reading predicted ages...'
          !open (12,file='detrital_ages/det_ages_'//trim(pbasin)//'.dat',&
          !      status='old')
          open (12,file='detrital_ages/ba_basins/det_ages_'//trim(pbasin)//'.dat',&
                status='old')
          !open (12,file='detrital_ages/nd_basins/det_ages_'//trim(pbasin)//'.dat',&
          !      status='old')
    3     read (12,*,end=4) dump
          plc=plc+1                                                             ! Find number of lines in model age files
          goto 3
    4     rewind (12)
          allocate(page(plc),pageu(plc),perate(plc),peratesc(plc))              ! Allocate model age/uncertainty and erate arrays
          do j=1,plc
            read(12,*) d1,d2,d3,d4,d5,d6,d7,page(j)                             ! Read in model predicted ages
            if (params%datapdf) then                                                   ! Assign predicted age uncertainties based on dataset, if doing a data comparison
              !pageup=pagemu                                                    ! Assign uncertainty (using mean measured value for dataset)
              pageup=pagemed                                                    ! Assign uncertainty (using median measured value for dataset)
              !pageu(j)=page(i)*0.022                                           ! Include standard deviation of errors (pagesd is the uncertainty s.d.)
              pageu(j)=page(j)*(pageup/100.)                                    ! Assign uncertainty to predicted ages
            else
              pageu(j)=page(j)*(params%pdf_pct_uncert/100.)                                    ! Assign predicted age uncertainties based on user-specified value above
            endif
            if (j==plc) write (*,'(a,i6,a,a)') 'Read ',plc,' predicted ages for basin ',trim(pbasin)
          enddo
          close(12)

          write (*,'(a)') 'Reading predicted erosion rates...'
          open (13,file='detrital_ages/erates_'//trim(pbasin)//'.dat',&
                status='old')
          !open (13,file='detrital_ages/nd_basins/erates_'//trim(pbasin)//'.dat',&
          !      status='old')
          do j=1,plc
            read(13,*) d1,d2,d3,perate(j),dump                                  ! Read in erosion rates for scaling age distribution
            !peratesc(j)=nint(perate(j)*100.)                                   ! Generate scaling factors by multiplying erosion rates by 100 and converting to integers
            if (j==plc) write (*,'(a,i6,a,a)') 'Read ',plc,' predicted erates for basin ',trim(pbasin)
          enddo
          close(13)
          peratemin=minval(perate)                                              ! Determine minimum erosion rate in model domain
          peratescl=1./peratemin                                                ! Set scaling value to ensure at least 1 occurance of min rate ages
          peratesc=nint(perate*peratescl)                                       ! Generate scaling factors by multiplying erosion rates by peratescl and converting to integers
          peratesum=sum(peratesc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DAVE: SHOULD THIS BE MODIFIED TO UPDATE THE LS ERATE FILE FOR MULTIPLE SETS?
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--- LANDSLIDE AGES/RATES ARE NOW ONLY READ DOWN BELOW FOR THE MONTE CARLO STUFF
!          if (lsppdf) then                                        ! If simulating landsliding, then read in the landslide ages and erosion rates
!            write (*,*) 'Landslide simulation requested, reading landslide ',&
!                    'erosion rates...'
!            open(14,file='detrital_ages/ls_ages_'//trim(pbasin)//'_'&
!                 //trim(simyr)//'yrs_iter1.dat',status='old')
!    5       read(14,*,end=6) dump
!            lsc=lsc+1
!            goto 5
!    6       rewind (14)
!            if (lsc.eq.0) then
!              write (*,*) '*** Warning: No landslide ages found in given catchment ***'
!              write (*,*) 'Assigning junk age of',lsagejunk,'and erosion rate of',&
!                      lseratejunk
!              lsc=lsc+1
!              allocate(lsage(lsc),lsageu(lsc),lserate(lsc))
!              do j=1,lsc
!                lsage(j)=lsagejunk
!                lsageu(j)=lsagejunk*(pageus/100)
!                lserate(j)=lseratejunk
!              enddo
!            else
!              allocate(lsage(lsc),lsageu(lsc),lserate(lsc))
!              do j=1,lsc
!                read(14,*) lsage(j),lserate(j)
!                if (datappdf .or. datamcpdfs) then                              ! If doing a data comparison, then use the data uncertainty specified above
!                  lsageu(j)=lsage(j)*(pageup/100.)
!                else                                                            ! Otherwise, use the user-specified value above
!                  lsageu(j)=lsage(j)*(pageus/100.)
!                endif
!              enddo
!            endif
!            close(14)
!            allocate(lseratesc(lsc))
!            lseratesc=nint(lserate*peratescl)                                   ! Scale landslide erosion rates
!            lseratesum=sum(lseratesc)
!          endif

! Generate data PDF
          if (params%datapdf) then                                                     ! Generate data PDF if doing a data comparison
            write (*,'(a)') 'Generating observed age PDF...'
            call get_pdf_size(oage,oageu,olc,onum,params%pdfmin,params%pdfmax, &
                              params%dx,params%calc_pdf_range)
            allocate(on(onum+1),opdf(onum+1))                                   ! Allocate data PDF arrays
            call make_age_pdf(oage,oageu,olc,onum,on,opdf,params%pdfmin,       &
                              params%pdfmax,params%dx,params%pdfscl,pi,ocnt)
            allocate(opdfv(ocnt))                                               ! Allocate PDF vector
            call make_age_pdfv(onum,params%pdfscl,on,opdf,opdfv,ocnt)
          endif

! Generate full predicted PDF
          if (params%fullppdf) then
            write (*,'(a)') 'Generating full predicted age PDF...'
            ! I've commented this out for not as I don't think it makes any
            ! sense to generate full landslide age PDFs (since the landslides
            ! are random). If uncommented, it needs to be fixed to work
            ! properly.
            ! dwhipp - 06/10

            !if (lsero) then
            !  call get_pdf_size(lsage,lsageu,lsc,pnum,pdfmin,pdfmax,dx)
            !  allocate(pn(pnum+1),ppdf(pnum+1))                                   ! Allocate data PDF arrays
            !  call make_age_pdf(lsage,lsageu,lsc,pnum,pn,ppdf,pdfmin,pdfmax,dx,&
            !                     pdfvsc,pi,pcnt)
            !  allocate(ppdfv(pcnt))                                               ! Allocate PDF vector
            !  call make_age_pdfv(pnum,pdfvsc,pn,ppdf,ppdfv,pcnt)
            !else
              call get_pdf_size(page,pageu,plc,pnum,params%pdfmin,             &
                                params%pdfmax,params%dx,params%calc_pdf_range)
              allocate(pn(pnum+1),ppdf(pnum+1))                                   ! Allocate data PDF arrays
              call make_age_pdf(page,pageu,plc,pnum,pn,ppdf,params%pdfmin,     &
                                params%pdfmax,params%dx,params%pdfscl,pi,pcnt)
              allocate(ppdfv(pcnt))                                               ! Allocate PDF vector
              call make_age_pdfv(pnum,params%pdfscl,pn,ppdf,ppdfv,pcnt)
            !endif
          endif


!!!
! Start monte carlo runs for testing model/data fits with Kuiper test
!!!!
          if (params%mcpdfs .or. params%datamcpdfs .or. params%ppdfmcpdfs) then
            write (*,'(a)') 'Starting Monte Carlo PDF generation...'
            if (params%datamcpdfs .or. params%ppdfmcpdfs) allocate(kuiper_res(params%mc_iter))       ! Allocate kuiper test array
            !allocate(pmc(int(pdfvsc),2*mc_iter))
            if (params%datamcpdfs) then
              nsc=1
            else
              nsc=size(params%numsamp)
            endif
            do m=1,nsc                                                            ! Loop through number of desired samples
              !pmc=0.
              cnt4=0
! COMMENTED OUT FOR TESTING
              if (params%datamcpdfs) then
                mcsamp=olc
              else
                mcsamp=params%numsamp(m)
              endif
              write (*,'(a,i7,a)') 'Running Monte Carlo simulation for ',mcsamp,' samples'
              allocate(pagemc(mcsamp),pageumc(mcsamp))
              if (mcboth) allocate(pagemc2(mcsamp),pageumc2(mcsamp))
              write (mcschar,'(i5)') mcsamp
              mcschar=adjustl(mcschar)
              do j=1,params%mc_iter                                             ! Model age distribution; This loop runs mc_iter times (usually ~10000)
                !call init_random_seed()
                jf=real(j)
                ! Read landslide ages here now!
                if (params%lsero) then                                                 ! If simulating landsliding, then read in the landslide ages and erosion rates
!                  write (*,*) 'Landslide simulation requested, reading landslide ',&
!                          'erosion rates...'
                  lsc=0                                                           ! Reset landslide age line counter
                  write(jc,'(i5)') j                                              ! Write iteration number to character
                  jc=adjustl(jc)
                  !open(14,file='detrital_ages/ls_ages_'//trim(pbasin)//'_'&
                  !     //trim(simyr)//'yrs_iter'//trim(jc)//'.dat',status='old')
                  !open(14,file='ls_ages/'//trim(simyr)//'yrs/ls_ages_'//trim(pbasin)//'_'&
                  !    //trim(simyr)//'yrs_iter'//trim(jc)//'.dat',status='old')
                  !open(14,file='ls_ages/'//trim(simyr)//'yrs/ls_ages_ndf1_'&
                  !    //trim(simyr)//'yrs_iter'//trim(jc)//'.dat',status='old')
                  open(14,file='ls_ages/ls_ages_'//trim(pbasin)//'_'&
                       //trim(params%simyr)//'yrs_iter'//trim(jc)//'.dat',status='old')
    7             read(14,*,end=8) dump
                  lsc=lsc+1
                  goto 7
    8             rewind (14)
                  if (lsc.eq.0) then
                    !write (*,*) '*** Warning: No landslide ages found in given catchment ***'
                    !write (*,*) 'Assigning junk age of',lsagejunk,'and erosion rate of',&
                    !        lseratejunk
                    lsc=lsc+1
                    allocate(lsage(lsc),lsageu(lsc),lserate(lsc))
                    do k=1,lsc
                      lsage(k)=params%lsagejunk
                      lsageu(k)=params%lsagejunk*(params%pdf_pct_uncert/100)
                      lserate(k)=params%lseratejunk
                    enddo
                  else
                    allocate(lsage(lsc),lsageu(lsc),lserate(lsc))
                    do k=1,lsc
                      read(14,*) lsage(k),lserate(k)
                      if (params%datappdf .or. params%datamcpdfs) then                        ! If doing a data comparison, then use the data uncertainty specified above
                        lsageu(k)=lsage(k)*(pageup/100.)
                      else                                                      ! Otherwise, use the user-specified value above
                        lsageu(k)=lsage(k)*(params%pdf_pct_uncert/100.)
                      endif
                    enddo
                  endif
                  close(14)
                  allocate(lseratesc(lsc))
                  lseratesc=nint(lserate*peratescl)                                   ! Scale landslide erosion rates
                  lseratesum=sum(lseratesc,lsc)
                endif

! Randomly grab mcsamp (n) grains from model distribution

                do k=1,mcsamp
                  call random_number(randflt)                                     ! Generate random number [0,1)
                  if (params%lsero) then
                    rint=int8(randflt*(lseratesum))+1                             ! Get random integer value within range of size of ls age dist.
                    eratechk=rint
                    cnt=0
                    do while (eratechk.gt.0)
                      cnt=cnt+1
                      eratechk=eratechk-lseratesc(cnt)
                    enddo
                    pagemc(k)=lsage(cnt)
                    pageumc(k)=lsageu(cnt)                  
                  else
                    rint=int8(randflt*(peratesum))+1                              ! Get random integer value within range of size of predicted age dist.
                    eratechk=rint
                    cnt=0
                    do while (eratechk.gt.0)
                      cnt=cnt+1
                      eratechk=eratechk-peratesc(cnt)
                    enddo
                    pagemc(k)=page(cnt)
                    pageumc(k)=pageu(cnt)                  
                  endif
                  !rint=int8(randflt*(eratesum))+1                                 ! Get random integer value within range of size of scaled age dist.
                  !eratechk=rint
                  !cnt=0
                  !do while (eratechk.gt.0)
                    !cnt=cnt+1
                    !eratechk=eratechk-eratetot(cnt)
                  !enddo
                  !pagemc(k)=pagetot(cnt)                                          ! Add random age to monte carlo age array
                  !pageumc(k)=pageutot(cnt)                                        ! Add associated uncertainty to monte carlo uncertainty array
                enddo
                call get_pdf_size(pagemc,pageumc,mcsamp,pnummc,params%pdfmin,  &
                                  params%pdfmax,params%dx,params%calc_pdf_range)
                allocate(pnmc(pnummc+1),ppdfmc(pnummc+1))                       ! Allocate data PDF arrays
                call make_age_pdf(pagemc,pageumc,mcsamp,pnummc,pnmc,ppdfmc,&
                                  params%pdfmin,params%pdfmax,params%dx,params%pdfscl,&
                                  pi,pcntmc)
                allocate(ppdfvmc(pcntmc))                                       ! Allocate PDF vector
                call make_age_pdfv(pnummc,params%pdfscl,pnmc,ppdfmc,ppdfvmc,   &
                                   pcntmc)

                !if (scale_ls_pdfs) then
                !
                !endif

! Run Kuiper test to get misfit between data and model
                if (params%datappdf) call kptwo(opdfv,ocnt,ppdfv,pcnt,olc,d,prob,h)
                if (params%datamcpdfs) call kptwo(opdfv,ocnt,ppdfvmc,pcntmc,mcsamp,d,&
                                prob,h)
                if (params%ppdfmcpdfs) call kptwo(ppdfv,pcnt,ppdfvmc,pcntmc,mcsamp,d,&
                                prob,h)
!                if (mcboth) then
!                  call kptwo(ppdfvmc,cnt3,ppdfvmc2,cnt8,mcsamp,d,prob,h)
!                else
!                  call kptwo(ppdfv,pcnt,ppdfvmc,cnt3,mcsamp,d,prob,h)
!                endif
!              endif
                if (params%datamcpdfs .or. params%ppdfmcpdfs) then
                  kuiper_res(j)=h                                               ! Store kuiper test result (0=pass;1=fail) for this iteration in kuiper results array
                  if (h.eq.0) cnt4=cnt4+1                                       ! Increment counter for number of models that pass Kuiper test
                endif

! Write out select PDFs
                !if (h.eq.0 .and. cnt4.le.100) then                               ! First 100 that pass the Kuiper test
                !if (h.eq.0 .and. cnt4.le.1) then                                 ! First one that passes the Kuiper test
                if (params%mcpdfs_out) then
                  if (j.le.params%num_mc_out) then                              ! First num_mc_out PDFs
                    write(hc,'(i5)') j
                    if (j.lt.10) hc(1:4)='0000'
                    if (j.lt.100) hc(1:3)='000'
                    if (j.lt.1000) hc(1:2)='00'
                    if (j.lt.10000) hc(1:1)='0'
                    open(24,file='age_pdf_output/mc_age_PDF_'//hc//'_'&
                         //trim(mcschar)//'_samples_'//trim(obasin)//'.dat',&
                         status='unknown')
                    if (params%tec_header) then
                      write(pnummcc,'(i10)') pnummc+1
                      pnummcc=adjustl(pnummcc)
                      write(24,'(a37)') 'TITLE="Monte Carlo predicted age PDF"'
                      write(24,'(a34)') 'VARIABLES="Age [Ma]" "Probability"'
                      write(24,'(a80)') 'ZONE I='//trim(pnummcc)//&
                             ' DATAPACKING=POINT T="Monte Carlo predicted PDF '&
                             //trim(obasin)//'"'
                    endif
                    do k=1,pnummc+1
                      write(24,*) pnmc(k),ppdfmc(k)
                    enddo
                    close(24)
!                    if (mcboth) then
!                      open(27,file='age_pdf_output/pass_mc2_age_PDF_'//hc//'_'&
!                          //trim(mcschar)//'_samples_'//trim(obasin)//'.dat',&
!                          status='unknown')
!                      do k=1,pnummc2+1
!                        write(27,*) pnmc2(k),ppdfmc2(k)
!                      enddo
!                      close(27)
!                    endif
                  endif
                endif

! Deallocate arrays
                deallocate(pnmc,ppdfmc,ppdfvmc)                                 ! Deallocate arrays reallocated during monte carlo sim
!                if (mcboth) then
!                  deallocate(pnmc2,ppsummc2,ppdfmc2,ppmc2,ppdfvmc2)               ! Deallocate arrays reallocated during monte carlo sim
!                endif
                if (params%lsero) deallocate(lsage,lsageu,lserate,lseratesc)
              enddo

              mc_iterf=real(params%mc_iter)
              if (params%datamcpdfs .or. params%ppdfmcpdfs) then
                kpct(i)=(1-(sum(kuiper_res)/mc_iterf))*100.                     ! Store percent of models that passed kuiper test for given basin

                ! Write output files
                open(20,file='kuiper_mc_results_'//trim(mcschar)//'_samples_'&
                     //trim(obasin)//'.dat',status='unknown')
                do j=1,params%mc_iter
                  write(20,*) kuiper_res(j)                                         ! Write out individual subset kuiper test result values to file for this basin
                enddo
                close(20)

                write(21,'(a12,a5,f10.2)') obasin,mcschar,kpct(i)                    ! Write out the summary percent of models that passed the Kuiper test

                !open(26,file='age_pdf_output/all_PDF_results_'//trim(obasin)//&
                !     '.dat',status='unknown')
                !write(26,*) pmc                                                    ! Write out individual subset kuiper test result values to file for this basin
                !close(26)
              endif

              ! Deallocate arrays
              deallocate(pagemc,pageumc)
              !if (mcboth) deallocate(pagemc2,pageumc2)

! End of main loop
            enddo
          endif

! Write out PDFs
              ! Write out data PDF
              if (params%opdf_out) then
                open(22,file='age_pdf_output/data_age_PDF_'//trim(obasin)//&
                    '.dat',status='unknown')
                if (params%tec_header) then
                  write(onumc,'(i10)') onum+1
                  onumc=adjustl(onumc)
                  write(22,'(a20)') 'TITLE="Data age PDF"'
                  write(22,'(a34)') 'VARIABLES="Age [Ma]" "Probability"'
                  write(22,'(a60)') 'ZONE I='//trim(onumc)//&
                           ' DATAPACKING=POINT T="Data PDF '//trim(obasin)//'"'
                endif
                do k=1,onum+1
                  write(22,*) on(k),opdf(k)
                enddo
                close(22)
              endif

              ! Write out full predicted PDF
              if (params%ppdf_out) then
                open(23,file='age_pdf_output/full_predicted_age_PDF_'&
                     //trim(obasin)//'.dat',status='unknown')
                if (params%tec_header) then
                  write(pnumc,'(i10)') pnum+1
                  pnumc=adjustl(pnumc)
                  write(23,'(a25)') 'TITLE="Predicted age PDF"'
                  write(23,'(a34)') 'VARIABLES="Age [Ma]" "Probability"'
                  write(23,'(a65)') 'ZONE I='//trim(pnumc)//&
                  ' DATAPACKING=POINT T="Predicted PDF '//trim(obasin)//'"'
                endif
                do k=1,pnum+1
                  write(23,*) pn(k),ppdf(k)
                enddo
                close(23)
              endif

! Deallocate arrays
          if (params%datapdf) deallocate(oage,oageu,on,opdf,opdfv)
          !if (datapdf) deallocate(oage,oageu)
          if (params%fullppdf) deallocate(page,pageu,pn,ppdf,ppdfv)
          deallocate(perate,peratesc)
          !if (lsero) deallocate(lsage,lsageu,lserate,lseratesc)
          !deallocate(pagetot,pageutot,eratetot)
          if (params%datamcpdfs .or. params%ppdfmcpdfs) deallocate(kuiper_res)

! Close open files
          close(21)
        endif
      enddo

! Sign off
write (*,'(a)') 'Execution complete'

! Exit
      end
