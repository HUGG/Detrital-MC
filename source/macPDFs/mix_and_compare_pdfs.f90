!
! mix_and_compare_pdfs.f90
!
! This program reads in PDFs output from detrital_ls and mixes the sediment in
! the tectonically driven exhumation PDF with that of the Monte Carlo landslide
! pdfs over a range of mixing ratios. The resulting new PDFs are compared and
! output is written to a summary file.
!
! dwhipp 04.13

program mix_and_compare_pdfs

    use mod_constants
    use pdf_functions

    implicit none

    !--- Variable declaration ---
    integer :: i,j,k,l,m,cnt,num,hm,onum,ocnt,pnum,pcnt,num_tr,num_mc_out,old
    integer :: mcsamp,pnummc,h,mix_bins,header_len,olc,pcntmc
    integer,allocatable :: kuiper_res(:,:),kcnt(:)
    real(kind=sp) :: pdfvsc,d,prob,num_mc_outf,lsero_frac,bgero_frac
    real(kind=sp),allocatable :: on(:),opdf(:),opdfv(:)
    real(kind=sp),allocatable :: pn(:),ppdf(:),ppdfv(:)
    real(kind=sp),allocatable :: pnmc(:),ppdfmc(:),ppdfvmc(:)
    real(kind=sp),allocatable :: kmc_all(:,:),kpct(:),ppdfmcmix(:)
    logical :: mcecho
    character(len=5)  :: hc,mcschar,totalptsc
    character(len=10) :: simyr(10)=(/'1.0','2.5','5.0','10.0','25.0','50.0','100.0','250.0','500.0','1000.0'/)
    character(len=12) :: obasin,pbasin
    character(len=80) :: dump

    ! Initialize variables
    obasin='ba4'
    pdfvsc=50.
    header_len=3
    num_tr=10
    num_mc_out=1000
    olc=111
    mix_bins=11

    onum=0
    pnum=0

    pbasin=obasin

    ! Say hello
    write (*,'(a)') '--------------------------------------------------------------------------------'
    write (*,'(a)') '------------------- Starting Mix and Compare PDFs (macPDFs) --------------------'
    write (*,'(a)') '--------------------------------------------------------------------------------'
    write (*,'(a)') ''
    
    ! Read in data PDF
    write (*,'(2a)') 'Reading data PDF for basin ',trim(obasin)
    open (10,file='age_pdf_output/bg_exhum/data_age_PDF_'//trim(obasin)//      &
                  '.dat',status='old')
70  read (10,*,end=71) dump
    onum=onum+1                                                                   ! Find number of lines in model age files
    goto 70
71  rewind (10)
    ! Read header
    do i=1,header_len
      read (10,*) dump
    enddo
    onum=onum-header_len
    allocate(on(onum),opdf(onum))                           ! Allocate 
    do i=1,onum
      read(10,*) on(i),opdf(i)
      if (i==onum) write (*,'(a,i6,a,a)') 'Read ',onum,' data PDF values for basin ',trim(obasin)
    enddo
    close(10)

    ! Format data PDF for Kuiper test
    ! Generate data PDF vector
    ocnt=0
    do i=1,onum
      hm=nint(pdfvsc*opdf(i))                                                   ! Set number of occurances of given age at current probability
      do j=1,hm
        ocnt=ocnt+1                                                             ! Count total number of ages in PDF vector for allocation below
      enddo
    enddo
    allocate(opdfv(ocnt))
    call make_age_pdfv(onum-1,pdfvsc,on,opdf,opdfv,ocnt)

    ! Read in full predicted PDF
    write (*,'(2a)') 'Reading full predicted PDF for basin ',trim(pbasin)
    open (11,file='age_pdf_output/bg_exhum/full_predicted_age_PDF_'            &
                  //trim(pbasin)//'.dat',status='old')
72  read (11,*,end=73) dump
    pnum=pnum+1                                                                   ! Find number of lines in model age files
    goto 72
73  rewind (11)
    ! Read header
    do i=1,header_len
      read (11,*) dump
    enddo
    pnum=pnum-header_len
    allocate(pn(pnum),ppdf(pnum))                           ! Allocate 
    do i=1,pnum
      read(11,*) pn(i),ppdf(i)
      if (i==pnum) write (*,'(a,i6,a,a)') 'Read ',pnum,' predicted PDF values for basin ',trim(pbasin)
    enddo
    close(11)

    ! Format full predicted PDF for Kuiper test
    ! Generate predicted PDF vector
    pcnt=0
    do i=1,pnum
      hm=nint(pdfvsc*ppdf(i))                                                   ! Set number of occurances of given age at current probability
      do j=1,hm
        pcnt=pcnt+1                                                             ! Count total number of ages in PDF vector for allocation below
      enddo
    enddo
    allocate(ppdfv(pcnt))
    call make_age_pdfv(pnum-1,pdfvsc,pn,ppdf,ppdfv,pcnt)

    mcsamp=olc
    write (mcschar,'(i5)') mcsamp
    mcschar=adjustl(mcschar)

    allocate(kuiper_res(mix_bins,num_mc_out),kcnt(mix_bins))
    allocate(kpct(mix_bins),kmc_all(mix_bins,num_tr))

    ! Loop over different residence times
    do i=1,num_tr
      kuiper_res=0
      kcnt=0
      mcecho=.false.

      ! Loop over all all MC ls pdfs
      write (*,'(3a)') 'Reading and comparing MC PDFs for ',trim(simyr(i)),' year residence time'
      do j=1,num_mc_out
        pnummc=0
        write(hc,'(i5)') j
        if (j.lt.10) hc(1:4)='0000'
        if (j.lt.100) hc(1:3)='000'
        if (j.lt.1000) hc(1:2)='00'
        if (j.lt.10000) hc(1:1)='0'
        ! Read in ls PDF
        !write (*,'(a)') 'Reading full predicted PDF for basin ',trim(pbasin)
        open (12,file='age_pdf_output/ls_'//trim(simyr(i))//'yr/mc_age_PDF_'//hc//&
                      '_'//trim(mcschar)//'_samples_'//trim(obasin)//'.dat',   &
                      status='unknown')
74      read (12,*,end=75) dump
        pnummc=pnummc+1                                                         ! Find number of lines in model age files
        goto 74
75      rewind (12)
        ! Read header
        do k=1,header_len
          read (12,*) dump
        enddo
        pnummc=pnummc-header_len
        allocate(pnmc(pnummc),ppdfmc(pnummc),ppdfmcmix(pnummc))                 ! Allocate 
        do k=1,pnummc
          read(12,*) pnmc(k),ppdfmc(k)
          if (.not. mcecho) then
            if (k==pnummc) then
              write (*,'(a,i6,a,a)') 'Read ',pnummc,' predicted MC PDF values for basin ',trim(obasin)
              mcecho=.true.
            endif
          endif
        enddo
        close(12)

        lsero_frac=0.0
        do k=1,mix_bins
          bgero_frac=1.0-lsero_frac

          ! Mix PDFs
          do l=1,pnummc
            ppdfmcmix(l)=ppdfmc(l)*lsero_frac+ppdf(l)*bgero_frac
          enddo

          ! Format full predicted PDF for Kuiper test
          ! Generate predicted PDF vector
          pcntmc=0
          do l=1,pnummc
            hm=nint(pdfvsc*ppdfmcmix(l))                                        ! Set number of occurances of given age at current probability
            do m=1,hm
              pcntmc=pcntmc+1                                                   ! Count total number of ages in PDF vector for allocation below
            enddo
          enddo
          allocate(ppdfvmc(pcntmc))
          call make_age_pdfv(pnummc-1,pdfvsc,pnmc,ppdfmcmix,ppdfvmc,pcntmc)

          !write (*,*) 'age_pdf_output/ls_'//trim(simyr(i))//'yr/mc_age_PDF_'//hc//&
          !            '_'//trim(mcschar)//'_samples_'//trim(obasin)//'.dat'


          !write(*,*) 'opdfv(1:10): ',opdfv(1:10)
          !write(*,*) 'ppdfvmc(1:10): ',ppdfvmc(1:10)


          call kptwo(opdfv,ocnt,ppdfvmc,pcntmc,mcsamp,d,prob,h)
          deallocate(ppdfvmc)

          kuiper_res(k,j)=h                                                     ! Store kuiper test result (0=pass;1=fail) for this iteration in kuiper results array


          !write (*,*) 'd: ',d
          !write (*,*) 'prob: ',prob
          !write (*,*) 'h: ',h



          if (h.eq.0) kcnt(k)=kcnt(k)+1                                         ! Increment counter for number of models that pass Kuiper test

          lsero_frac=lsero_frac+(1./real(mix_bins-1))
        enddo
        deallocate(pnmc,ppdfmc,ppdfmcmix)                                       ! Deallocate 
      enddo

      num_mc_outf=real(num_mc_out)
      do j=1,mix_bins
        kpct(j)=(1.-(sum(kuiper_res(j,:))/num_mc_outf))*100.                     ! Store percent of models that passed kuiper test for given basin
        ! Add new kpct to 2d array
        kmc_all(j,i)=kpct(j)


        write (*,*) 'sum(kuiper_res(j,:)): ',sum(kuiper_res(j,:))


      enddo

      !! Write output files
      !open(20,file='kuiper_mc_results_'//trim(mcschar)//'_samples_'&
      !             //trim(obasin)//'.dat',status='unknown')
      !do j=1,mc_iter
      !  write(20,*) kuiper_res(j)                                               ! Write out individual subset kuiper test result values to file for this basin
      !enddo
      !close(20)
    enddo

    ! write out 2d kpct array
    write (totalptsc,'(i5)') mix_bins*num_tr
    open (20,file='kmc_mix_summary_'//trim(obasin)//'.dat',status='unknown')    ! Open basin summary results file
    write(20,'(a26)') 'TITLE="PDF mixing results"'
    write(20,'(a76)') 'VARIABLES="Residence time [a]" "LS fraction" "BG fraction" "Percent passing"'
    write(20,'(a65)') 'ZONE I='//trim(totalptsc)//&
                      ' DATAPACKING=POINT T="PDF mixing for basin '//trim(obasin)//'"'
    do i=1,num_tr
      lsero_frac=0.0
      do j=1,mix_bins
        bgero_frac=1.0-lsero_frac
        write(20,'(a,3f12.4)') trim(simyr(i)),lsero_frac,bgero_frac,kmc_all(j,i)
        lsero_frac=lsero_frac+(1./real(mix_bins-1))
      enddo
    enddo
    close(20)

    ! Sign off
    write (*,'(a)') ''
    write (*,'(a)') '--------------------------------------------------------------------------------'
    write (*,'(a)') '-------------------------- macPDF execution complete ---------------------------'
    write (*,'(a)') '------------------------------ Have a nice day! --------------------------------'
    write (*,'(a)') '--------------------------------------------------------------------------------'

    ! Exit
    end