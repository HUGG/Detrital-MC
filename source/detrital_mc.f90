! detrital_mc.f90
!
! This fortran code uses a monte carlo simulation method to randomly grab n
! samples from a distribution of model predicted cooling ages, assign a
! designated uncertainty to the samples, generate a PDF of the model predicted
! cooling ages, compare that PDF to an observed age PDF using a Kuiper test,
! and finally record the results of that test. This process is repeated a
! large number of times (~10000).
!
! This is the main program file for version 3.1 of Detrital MC.
!
! dwhipp - 06.16
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program detrital_mc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE definitions
      USE pdf_functions

      IMPLICIT NONE

      type (detrital_params) :: params
      type (basin_information) :: basin_info(100)

      real(kind=sp), allocatable :: oage(:),oageu(:),on(:)
      real(kind=sp), allocatable :: opdf(:),opdfv(:)
      real(kind=sp), allocatable :: page(:),pageu(:),perate(:)
      real(kind=sp), allocatable :: pagemc(:),pageumc(:),pnmc(:)
      real(kind=sp), allocatable :: ppdfmc(:),ppdfvmc(:),peratemc(:)
      real(kind=sp), allocatable :: lsage(:),lsageu(:)
      real(kind=sp), allocatable :: lserate(:)
      real(kind=sp), allocatable :: kpct(:),ppdf(:),ppdfv(:),pn(:)
      real(kind=sp), allocatable :: ocdf(:),pcdf(:),pcdfmc(:)
      real(kind=sp), allocatable :: oecdf(:),pecdf(:),pecdfmc(:)
      real(kind=sp), allocatable :: oageupct(:)
      integer(kind=sp), allocatable :: peratesc(:),kuiper_res(:),lseratesc(:)
      integer(kind=sp), allocatable :: oeratesc(:),peratescmc(:)
      !integer*8, allocatable :: lseratesc(:)
      real(kind=sp) :: d,prob,pagemu,pagemed,pagesd,mc_iterf,jf
      real(kind=sp) :: pi,peratemin,peratescl,dumpr,lseratescl
      real(kind=sp) :: pctpagemu,pctpagemed,pctpagesd
      real(kind=sp) :: dump1,dump2,dump3,dump4,dump5,dump6,dump7,dump8,dump9
      real(kind=sp) :: dump10,dump11,dump12,dump13,dump14,dump15,dump16,dump17
      real(kind=sp) :: dump18,dump19,dump20,dump21,dump22,dump23,dump24,dump25
      real(kind=sp) :: dump26,dump27,dump28,dump29,dump30,dump31,dump32
      real(kind=sp) :: randflt,eps,peratesum
      integer(kind=sp) :: onum,h,i,j,k,pnum,olc,oeratesum,peratescsummc
      integer(kind=sp) :: plc,cnt,pnummc,m,peratescsum
      integer(kind=sp) :: lsc,mcsamp,nsc,ocnt,pcnt,pcntmc
      !integer :: lseratesum,rint,eratechk
      !integer*8 :: peratesum,rint,eratechk
      integer(kind=dp) :: lseratesum,rint,eratechk
      logical, allocatable :: peratemask(:)
      logical datacomp,fileexist
      !logical lsppdf,lspdf_out
      character(len=5)  :: hc,jc,mcschar
      character(len=10) :: onumc,pnumc,pnummcc
      character(len=80) :: dump
!      character(len=12), allocatable :: obasin(:),pbasin(:)

      ! This stuff will be needed for auto-sizing the PDF arrays using the data
      ! uncertainties. Leaving it here for now...
      !real(kind=sp) :: probcut,oagemin,oagemax,pagemin,pagemax,pagemcmin
      !real(kind=sp) :: pagemcmax
      !probcut=0.005

      !mcboth=.false.
      eps=1.e-8

! Initialize random number generator
      call init_random_seed()

! Variable initialization/declaration
      pi=atan(1.)*4.

! Write program starting info
      write (*,'(a)') '#------------------------------------------------------------------------------#'
      write (*,'(a)') 'Detrital Monte Carlo PDF creator started'
      write (*,'(a)') 'Version 3.1 (dwhipp - June 2016)'

! Read in the input file
      call read_input_file(params,basin_info)

! Allocate kuiper test results summary array
      allocate(kpct(params%num_basins))

      do i=1,params%num_basins                                                  ! Loop through number of subbasins to analyze
        plc=0                                                                   ! Reset predicted age line counter
        lsc=0                                                                   ! Reset landslide age line counter

        write (*,*) ''
        write (*,'(a,i3,a,i3,a,a,a)') 'Processing basin ',i,' of ',params%num_basins,' (',trim(basin_info(i)%obasin_name),')'

        if (params%datamcpdfs .or. params%ppdfmcpdfs) then
          open(21,file='kmc_percent_pass_summary_'//trim(basin_info(i)%obasin_name)//'.dat',&      ! Open basin summary results file
               status='unknown')
          write(21,'(a27)') '            Percent passing'                       ! Write header
          write(21,'(a23)') 'Basin       Kuiper test'
        endif

        if (params%datapdf .or. params%datappdf) then                           ! If using data, read in observed cooling ages and generate PDF
          write (*,'(a)') 'Data comparison requested, reading observed ages...'
          open (11,file='data/observed_ages/'//trim(basin_info(i)%obasin_name)//'.dat',status='old')
          read(11,*) olc
          allocate(oage(olc),oageu(olc),oageupct(olc),oeratesc(olc))            ! Allocate observed age and age uncertainty arrays
          do j=1,olc
            read(11,*) oage(j),oageu(j)                                         ! Fill observed sample arrays
            oageupct(j)=oageu(j)/oage(j)*100.                                   ! Store percent uncertainties
            if (j==olc) write (*,'(a,i3,a,a)') 'Read ',olc,' ages for basin ',trim(basin_info(i)%obasin_name)
          enddo

          ! Fill erate array with dummy values
          oeratesc=1
          oeratesum=sum(oeratesc)

          ! Calculate percentage age uncertainties
          pagemu=sum(oageu(:))/real(olc)                                        ! Mean 1-sigma uncertainty in observed ages
          pctpagemu=sum(oageupct(:))/real(olc)                                  ! Mean percent uncertainty in observed ages
          pagemed=median(oageu(:))                                              ! Median 1-sigma uncertainty in observed ages
          pctpagemed=median(oageupct(:))                                        ! Median percent uncertainty in observed ages
          pagesd=sqrt(sum((oageu(:)-pagemu)**2)/real(olc))                      ! Standard deviation in 1-sigma uncertainty in observed ages
          pctpagesd=sqrt(sum((oageupct(:)-pctpagemu)**2)/real(olc))             ! Standard deviation in percent uncertainty in observed ages
          deallocate(oageupct)
          if (params%datapdf) then
            if (params%obs_uncert_type == 1) then
              write (*,'(a,f6.2,a)') 'Applying mean uncertainty (',pctpagemu,'%) to predicted ages'
            elseif (params%obs_uncert_type == 2) then
              write (*,'(a,f6.2,a)') 'Applying median uncertainty (',pctpagemed,'%) to predicted ages'
            elseif (params%obs_uncert_type == 3) then
              write (*,'(a,f6.2,a)') 'Applying standard deviation in uncertainty (',pctpagesd,'%) to predicted ages'
            elseif (params%obs_uncert_type == 4) then
              write (*,'(a,f6.2,a)') 'Applying user-specified uncertainty (',params%pdf_pct_uncert,'%) to predicted ages'
            endif
          else
            write (*,'(a,f6.2,a)') 'Applying user-specified uncertainty (',params%pdf_pct_uncert,'%) to predicted ages'
          endif
          close(11)
        endif

        write (*,'(a)') 'Reading predicted ages...'
        if (basin_info(i)%page_ftype==1) then
          open (12,file='data/predicted_ages/'//trim(basin_info(i)%pbasin_name)//'/Comparison.txt',&
                status='old')
          read(12,*) plc
          allocate(page(plc),pageu(plc),perate(plc),peratesc(plc))              ! Allocate model age/uncertainty and erate arrays
          do j=1,plc
            read(12,*) dump1,dump2,dump3,dump4,dump5,dump6,dump7,dump8,dump9,  &
                       dump10,dump11,dump12,dump13,dump14,dump15,dump
            if (basin_info(i)%page_sys == 1) page(j)=dump7                      ! Store predicted AHe ages
            if (basin_info(i)%page_sys == 2) page(j)=dump9                      ! Store predicted AFT ages
            if (basin_info(i)%page_sys == 3) page(j)=dump11                     ! Store predicted ZHe ages
            if (basin_info(i)%page_sys == 4) page(j)=dump13                     ! Store predicted ZFT ages
            if (basin_info(i)%page_sys == 5) page(j)=dump15                     ! Store predicted MAr ages
            perate(j)=dump5                                                     ! Store predicted erosion rate

            if (params%datapdf) then                                            ! Assign predicted age uncertainties based on dataset, if doing a data comparison
              if (params%obs_uncert_type == 1) then                             ! Assign uncertainty using mean percent uncertainty in observed ages
                pageu(j) = page(j)*(pctpagemu/100.)
              elseif (params%obs_uncert_type == 2) then                         ! Assign uncertainty using median percent uncertainty in observed ages
                pageu(j) = page(j)*(pctpagemed/100.)
              elseif (params%obs_uncert_type == 3) then                         ! Assign uncertainty using standard deviation in percent uncertainty in observed ages
                pageu(j) = page(j)*(pctpagesd/100.)
              elseif (params%obs_uncert_type == 4) then                         ! Assign user-specified percent uncertainty
                pageu(j)=page(j)*(params%pdf_pct_uncert/100.)
              endif
            else
              pageu(j)=page(j)*(params%pdf_pct_uncert/100.)                     ! Assign predicted age uncertainties based on user-specified value above
            endif
          enddo
        else
          open (12,file='data/predicted_ages/'//trim(basin_info(i)%pbasin_name)//'.csv',&
                status='old')
          read(12,*) plc
          allocate(page(plc),pageu(plc),perate(plc),peratesc(plc))              ! Allocate model age/uncertainty and erate arrays
          do j=1,plc
            read(12,*) dump1,dump2,dump3,dump4,dump5,dump6,dump7,dump8,dump9,  &
                       dump10,dump11,dump12,dump13,dump14,dump15,dump16,dump17,&
                       dump18,dump19,dump20,dump21,dump22,dump23,dump24,dump25,&
                       dump26,dump27,dump28,dump29,dump30,dump31,dump32
            if (basin_info(i)%page_col == 6) page(j)=dump6                      ! Store predicted AHe age
            if (basin_info(i)%page_col == 7) page(j)=dump7                      ! Store predicted ZHe age
            if (basin_info(i)%page_col == 8) page(j)=dump8                      ! Store predicted AFT age
            if (basin_info(i)%page_col == 9) page(j)=dump9                      ! Store predicted ZFT age
            if (basin_info(i)%page_col == 12) page(j)=dump12                    ! Store predicted MAr age
            if (basin_info(i)%perate_col == 5) perate(j)=dump5                  ! Store predicted erosion rate
            if (basin_info(i)%perate_col == 16) then                            ! Scale ages using bedrock geology
              if (dump16 > eps .and. dump16 < 7.0) perate(j)=basin_info(i)%geol_scale_factor(int(dump16))
              if (basin_info(i)%scale_by_uplift_velo) perate(j)=perate(j)*dump5
            endif
            if (basin_info(i)%perate_col == 17) then                            ! Scale ages using glaciers
              if (dump17 > eps) then
                perate(j)=basin_info(i)%geol_scale_factor(1)
              else
                perate(j)=basin_info(i)%geol_scale_factor(2)
              endif
              if (basin_info(i)%scale_by_uplift_velo) perate(j)=perate(j)*dump5
            endif
            if (basin_info(i)%perate_col == 18) then                            ! Scale ages using moraines
              if (dump18 > eps) then
                perate(j)=basin_info(i)%geol_scale_factor(1)
              else
                perate(j)=basin_info(i)%geol_scale_factor(2)
              endif
              if (basin_info(i)%scale_by_uplift_velo) perate(j)=perate(j)*dump5
            endif
            if (basin_info(i)%perate_col == 19) then                            ! Scale ages using rock glaciers
              if (dump19 > eps) then
                perate(j)=basin_info(i)%geol_scale_factor(1)
              else
                perate(j)=basin_info(i)%geol_scale_factor(2)
              endif
              if (basin_info(i)%scale_by_uplift_velo) perate(j)=perate(j)*dump5
            endif
            if (basin_info(i)%perate_col == 20) perate(j)=dump16                ! Scale ages using Ksn
            if (basin_info(i)%perate_col == 21) perate(j)=dump17                ! Scale ages using Ksn_t045
            if (basin_info(i)%perate_col == 22) perate(j)=dump18                ! Scale ages using Ksn_t2
            if (basin_info(i)%perate_col == 23) perate(j)=dump19                ! Scale ages using Ksn_t3
            if (basin_info(i)%perate_col == 31) perate(j)=dump27                ! Scale ages using ssp_t2b31
            if (basin_info(i)%perate_col == 32) perate(j)=dump28                ! Scale ages using ssp_t3b42
            if (basin_info(i)%perate_col == 99) then
              if (dump17 > eps) then
                perate(j)=basin_info(i)%geol_scale_factor(1)
              elseif (dump18 > eps) then
                perate(j)=basin_info(i)%geol_scale_factor(2)
              elseif (dump19 > eps) then
                perate(j)=basin_info(i)%geol_scale_factor(3)
              else
                perate(j)=basin_info(i)%geol_scale_factor(4)
              endif
              if (basin_info(i)%scale_by_uplift_velo) perate(j)=perate(j)*dump5
            endif
            if (perate(j) < eps) perate(j)=0.0                                  ! Ensure no erosion rate scalars are negative

            if (params%datapdf) then                                            ! Assign predicted age uncertainties based on dataset, if doing a data comparison
              if (params%obs_uncert_type == 1) then                             ! Assign uncertainty using mean percent uncertainty in observed ages
                pageu(j) = page(j)*(pctpagemu/100.)
              elseif (params%obs_uncert_type == 2) then                         ! Assign uncertainty using median percent uncertainty in observed ages
                pageu(j) = page(j)*(pctpagemed/100.)
              elseif (params%obs_uncert_type == 3) then                         ! Assign uncertainty using standard deviation in percent uncertainty in observed ages
                pageu(j) = page(j)*(pctpagesd/100.)
              elseif (params%obs_uncert_type == 4) then                         ! Assign user-specified percent uncertainty
                pageu(j)=page(j)*(params%pdf_pct_uncert/100.)
              endif
            else
              pageu(j)=page(j)*(params%pdf_pct_uncert/100.)                     ! Assign predicted age uncertainties based on user-specified value above
            endif
          enddo

          if (params%scale_erates) then
            if (params%scaletype == 1) then                                     ! Normalized
              perate=perate/maxval(perate)
            elseif (params%scaletype == 2) then                                 ! Power law exponent = 2.0
              perate=perate**2.0
              perate=perate/maxval(perate)
            elseif (params%scaletype == 3) then                                 ! Power law exponent = 5.0
              perate=perate**5.0
              perate=perate/maxval(perate)
            elseif (params%scaletype == 4) then
              write(*,*) 'Not yet implemented.'
            elseif (params%scaletype == 5) then
              write(*,*) 'Not yet implemented.'
            elseif (params%scaletype == 6) then
              write(*,*) 'Not yet implemented.'
            elseif (params%scaletype == 7) then
              write(*,*) 'Not yet implemented.'
            endif
          endif
        endif
        write (*,'(a,i6,a,a)') 'Read ',plc,' predicted ages/erosion rates for basin ',trim(basin_info(i)%pbasin_name)
        close(12)

        if (params%dist_size > 0) then                                          ! Scale age prevalence to fit desired number of ages
          peratesum = sum(perate)                                               ! Sum predicted erosion rates
          peratescl = real(params%dist_size) / peratesum                        ! Calculate age distribution scaling factor
        else
          peratemin=minval(perate)                                              ! Determine minimum erosion rate in model domain
          if (peratemin < eps) then
            allocate(peratemask(plc))
            peratemask=(perate > eps)
            peratemin=minval(perate,mask=peratemask)
            deallocate(peratemask)
          endif
          peratescl=1./peratemin                                                ! Set scaling value to ensure at least 1 occurance of min rate ages
        endif
        peratesc=nint(perate*peratescl)                                         ! Generate scaling factors by multiplying erosion rates by peratescl and converting to integers
        peratescsum=sum(peratesc)

! Generate data PDF
        if (params%datapdf) then                                                ! Generate data PDF if doing a data comparison
          write (*,'(a)') 'Generating observed age PDF...'
          call get_pdf_size(oage,oageu,olc,onum,params%pdfmin,params%pdfmax, &
                            params%dx,params%calc_pdf_range)
          allocate(on(onum+1),opdf(onum+1))                                     ! Allocate data PDF arrays
          ! I think the data PDF should still use alpha=1.0, so I've hard-coded that in
          call make_age_pdf(oage,oageu,1.0,oeratesc,olc,onum,on,opdf,          &
                            params%pdfmin,params%dx,params%pdfscl,pi,ocnt,0)
          ! Calculate cumulative distributions, if requested
          if (params%ocdf_out) then
            if (params%ecdfs) then
              allocate(oecdf(onum+1))
              call make_age_ecdf(oage,oeratesc,oeratesum,on,olc,onum,oecdf)
            else
              allocate(ocdf(onum+1))
              call make_age_cdf(opdf,onum,params%dx,ocdf)
            endif
          endif
          allocate(opdfv(ocnt))                                                 ! Allocate PDF vector
          call make_age_pdfv(onum,params%pdfscl,on,opdf,opdfv,ocnt)
        endif

! Generate full predicted PDF
        if (params%fullppdf) then
          write (*,'(a)') 'Generating full predicted age PDF...'
          call get_pdf_size(page,pageu,plc,pnum,params%pdfmin,                 &
                            params%pdfmax,params%dx,params%calc_pdf_range)
          allocate(pn(pnum+1),ppdf(pnum+1))                                     ! Allocate data PDF arrays
          if (params%alphain < 0.d0) then
            params%alpha = (4.0/(3.0*plc))**0.2
          else
            params%alpha = params%alphain
          endif
          call make_age_pdf(page,pageu,params%alpha,peratesc,plc,pnum,pn,      &
                            ppdf,params%pdfmin,params%dx,params%pdfscl,pi,     &
                            pcnt,1)
          if (params%pcdf_out) then
            if (params%ecdfs) then
              allocate(pecdf(pnum+1))
              call make_age_ecdf(page,peratesc,peratescsum,pn,plc,pnum,pecdf)
            else
              allocate(pcdf(pnum+1))
              call make_age_cdf(ppdf,pnum,params%dx,pcdf)
            endif
          endif
          allocate(ppdfv(pcnt))                                                 ! Allocate PDF vector
          call make_age_pdfv(pnum,params%pdfscl,pn,ppdf,ppdfv,pcnt)
        endif

!!!
! Start monte carlo runs for testing model/data fits with Kuiper test
!!!!
        if (params%mcpdfs .or. params%datamcpdfs .or. params%ppdfmcpdfs) then
          write (*,'(a)') 'Starting Monte Carlo PDF generation...'
          if (params%datamcpdfs .or. params%ppdfmcpdfs) allocate(kuiper_res(params%mc_iter))       ! Allocate kuiper test array
          if (params%datamcpdfs) then
            nsc=1
          else
            nsc=size(params%numsamp)
          endif
          do m=1,nsc                                                            ! Loop through number of desired samples
            if (params%datamcpdfs) then
              mcsamp=olc
            else
              mcsamp=params%numsamp(m)
            endif
            ! Calculate optimal alpha value, if input value was negative
            if (params%alphain < 0.d0) then
              params%alpha = (4.0/(3.0*mcsamp))**0.2
            else
              params%alpha = params%alphain
            endif
            write (*,'(a,i7,a)') 'Running Monte Carlo simulation for ',mcsamp,' samples'
            allocate(pagemc(mcsamp),pageumc(mcsamp),peratemc(mcsamp),peratescmc(mcsamp))
            !if (mcboth) allocate(pagemc2(mcsamp),pageumc2(mcsamp))
            write (mcschar,'(i5)') mcsamp
            mcschar=adjustl(mcschar)
            do j=1,params%mc_iter                                               ! Model age distribution; This loop runs mc_iter times (usually ~10000)
              if (mod(j,params%mc_iter/10) == 0) then
                write (*,'(f5.1,a,i5,a)') 100*real(j)/real(params%mc_iter),'% (',j,') Monte Carlo simulations complete...'
              endif
              !call init_random_seed()
              jf=real(j)
              ! Read landslide ages here now!
              if (params%lsero) then                                                 ! If simulating landsliding, then read in the landslide ages and erosion rates
                lsc=0                                                           ! Reset landslide age line counter
                write(jc,'(i5)') j                                              ! Write iteration number to character
                jc=adjustl(jc)
                !
                ! THIS NEEDS TO BE ADJUSTED TO WORK FOR ASCII OR BINARY INPUT FILES!!!!
                !
                open(14,file='data/ls_ages/ls_ages_'//trim(basin_info(i)%pbasin_name)//'_'&
                     //trim(adjustl(params%simyr))//'yrs_iter'//trim(jc)//'.dat',status='old')
  7             read(14,*,end=8) dump
                lsc=lsc+1
                goto 7
  8             rewind (14)
                if (lsc.eq.0) then
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
                    if (params%datappdf .or. params%datamcpdfs) then            ! Assign predicted age uncertainties based on dataset, if doing a data comparison
                      if (params%obs_uncert_type == 1) then                     ! Assign uncertainty using mean percent uncertainty in observed ages
                        lsageu(k) = lsage(k)*(pctpagemu/100.)
                      elseif (params%obs_uncert_type == 2) then                 ! Assign uncertainty using median percent uncertainty in observed ages
                        lsageu(k) = lsage(k)*(pctpagemed/100.)
                      elseif (params%obs_uncert_type == 3) then                 ! Assign uncertainty using standard deviation in percent uncertainty in observed ages
                        lsageu(k) = lsage(k)*(pctpagesd/100.)
                      endif
                    else
                      lsageu(k)=lsage(k)*(params%pdf_pct_uncert/100.)           ! Assign predicted age uncertainties based on user-specified value above
                    endif
                  enddo
                endif
                close(14)
                allocate(lseratesc(lsc))


                ! SHOULD ALL OF THIS BE IN TERMS OF PERATE, RATHER THAN LSERATE???
                if (params%dist_size > 0) then                                  ! Scale age prevalence to fit desired number of ages
                  lseratesum = sum(lserate)                                     ! Sum predicted landsliding erosion rates
                  peratescl = real(params%dist_size) / lseratesum               ! Re-calculate age distribution scaling factor
                endif
                lseratesc=nint(lserate*peratescl)                               ! Scale landslide erosion rates
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
                  peratemc(k)=lserate(cnt)
                else
                  rint=int8(randflt*(peratescsum))+1                            ! Get random integer value within range of size of predicted age dist.
                  eratechk=rint
                  cnt=0
                  do while (eratechk.gt.0)
                    cnt=cnt+1
                    eratechk=eratechk-peratesc(cnt)
                  enddo
                  pagemc(k)=page(cnt)
                  pageumc(k)=pageu(cnt)
                  peratemc(k)=perate(cnt)
                endif
              enddo
              



              ! SHOULD THIS BE DIFFERENT FOR LS VERSUS NO LS CASES???
              ! if (params%dist_size > 0) then
              !   peratesummc = sum(lserate
              ! else
              !   peratescmc=nint(peratemc*peratescl)                             ! Scale erosion rates
              !   peratescsummc=sum(peratescmc,mcsamp)
              ! endif

              peratescmc=nint(peratemc*peratescl)                             ! Scale erosion rates
              peratescsummc=sum(peratescmc,mcsamp)



              call get_pdf_size(pagemc,pageumc,mcsamp,pnummc,params%pdfmin,  &
                                params%pdfmax,params%dx,params%calc_pdf_range)
              allocate(pnmc(pnummc+1),ppdfmc(pnummc+1))                         ! Allocate data PDF arrays
              call make_age_pdf(pagemc,pageumc,params%alpha,peratescmc,mcsamp, &
                                pnummc,pnmc,ppdfmc,params%pdfmin,params%dx,    &
                                params%pdfscl,pi,pcntmc,0)
              if (params%mccdfs_out) then
                if (params%ecdfs) then
                  allocate(pecdfmc(pnummc+1))
                  call make_age_ecdf(pagemc,peratescmc,peratescsummc,pnmc,mcsamp,pnummc,pecdfmc)
                else
                  allocate(pcdfmc(pnummc+1))
                  call make_age_cdf(ppdfmc,pnummc,params%dx,pcdfmc)
                endif
              endif
              allocate(ppdfvmc(pcntmc))                                       ! Allocate PDF vector
              call make_age_pdfv(pnummc,params%pdfscl,pnmc,ppdfmc,ppdfvmc,   &
                                 pcntmc)

! Run Kuiper test to get misfit between data and model
              if (params%kuipernew) then
                ! Calculate h for comparison of observed and full predicted PDFs
                if (params%datappdf) then
                  if (params%ecdfs) then
                    d = maxval(oecdf-pecdf)+maxval(pecdf-oecdf)
                  else
                    d = maxval(ocdf-pcdf)+maxval(pcdf-ocdf)
                  endif
                  h = kuiper(params%kalpha,d,olc)
                endif
                ! Calculate h for comparison of observed and MC predicted PDFs
                if (params%datamcpdfs) then
                  if (params%ecdfs) then
                    d = maxval(oecdf-pecdfmc)+maxval(pecdfmc-oecdf)
                    !write(*,*) 'd_ecdf: ',d
                  else
                    d = maxval(ocdf-pcdfmc)+maxval(pcdfmc-ocdf)
                    !write(*,*) 'd_cdf: ',d
                  endif
                  h = kuiper(params%kalpha,d,mcsamp)
                endif
                ! Calculate h for comparison of full predicted and MC PDFs
                if (params%ppdfmcpdfs) then
                  if (params%ecdfs) then
                    d = maxval(pecdf-pecdfmc)+maxval(pecdfmc-pecdf)
                  else
                    d = maxval(pcdf-pcdfmc)+maxval(pcdfmc-pcdf)
                  endif
                  h = kuiper(params%kalpha,d,mcsamp)
                endif
              else
                ! Calculate h for comparison of observed and full predicted PDFs
                if (params%datappdf) call kptwo(opdfv,ocnt,ppdfv,pcnt,olc,d,   &
                                                prob,h)
                ! Calculate h for comparison of observed and MC predicted PDFs
                if (params%datamcpdfs) call kptwo(opdfv,ocnt,ppdfvmc,pcntmc,   &
                                                  mcsamp,d,prob,h)
                ! Calculate h for comparison of full predicted and MC PDFs
                if (params%ppdfmcpdfs) call kptwo(ppdfv,pcnt,ppdfvmc,pcntmc,   &
                                                  mcsamp,d,prob,h)
              endif

              if (params%datamcpdfs .or. params%ppdfmcpdfs) then
                kuiper_res(j)=h                                                 ! Store kuiper test result (0=pass;1=fail) for this iteration in kuiper results array
              endif

! Write out select PDFs
              if (params%mcpdfs_out .or. params%mccdfs_out) then
                if (j.le.params%num_mc_out) then                                ! First num_mc_out PDFs
                  write(hc,'(i5)') j
                  if (j.lt.10) hc(1:4)='0000'
                  if (j.lt.100) hc(1:3)='000'
                  if (j.lt.1000) hc(1:2)='00'
                  if (j.lt.10000) hc(1:1)='0'
                  ! Write out Monte Carlo PDFs
                  if (params%mcpdfs_out) then
                    if (params%veusz_output) then
                      if (j == 1) then
                        open(24,file='age_pdf_output/mc_age_PDFs_'//trim(mcschar)//&
                             '_samples_'//trim(basin_info(i)%obasin_name)//'_veusz.csv',&
                             status='unknown',action='write')
                        write(24,'(a20)') 'Age [Ma],Probability'
                        do k=1,pnummc+1
                          write(24,'(e13.6,a1,e13.6)') pnmc(k),',',ppdfmc(k)
                        enddo
                        close(24)
                      else
                        open(24,file='age_pdf_output/mc_age_PDFs_'//trim(mcschar)//&
                             '_samples_'//trim(basin_info(i)%obasin_name)//'_veusz.csv',&
                             status='old',position='append',action='write')
                        write(24,'(a7)') 'nan,nan'
                        do k=1,pnummc+1
                          write(24,'(e13.6,a1,e13.6)') pnmc(k),',',ppdfmc(k)
                        enddo
                        close(24)
                      endif
                    else
                      open(24,file='age_pdf_output/mc_age_PDF_'//hc//'_'&
                           //trim(mcschar)//'_samples_'//trim(basin_info(i)%obasin_name)//'.dat',&
                           status='unknown')
                      if (params%tec_header) then
                        write(pnummcc,'(i10)') pnummc+1
                        pnummcc=adjustl(pnummcc)
                        write(24,'(a38)') 'TITLE="Monte Carlo predicted age PDFs"'
                        write(24,'(a34)') 'VARIABLES="Age [Ma]" "Probability"'
                        write(24,'(a80)') 'ZONE I='//trim(pnummcc)//&
                               ' DATAPACKING=POINT T="Monte Carlo predicted PDFs '&
                               //trim(basin_info(i)%obasin_name)//'"'
                      endif
                      do k=1,pnummc+1
                        write(24,'(e13.6,e13.6)') pnmc(k),ppdfmc(k)
                      enddo
                      close(24)
                    endif
                  endif
                  ! Write out Monte Carlo CDFs or ECDFs
                  if (params%mccdfs_out) then
                    ! Write out Monte Carlo empirical cumulative distributions
                    if (params%ecdfs) then
                      if (params%veusz_output) then
                        if (j == 1) then
                          open(27,file='age_pdf_output/mc_age_ECDFs_'//trim(mcschar)//&
                               '_samples_'//trim(basin_info(i)%obasin_name)//'_veusz.csv',&
                               status='unknown',action='write')
                          write(27,'(a31)') 'Age [Ma],Cumulative probability'
                          do k=1,pnummc+1
                            write(27,'(e13.6,a1,e13.6)') pnmc(k),',',pecdfmc(k)
                          enddo
                          close(27)
                        else
                          open(27,file='age_pdf_output/mc_age_ECDFs_'//trim(mcschar)//&
                               '_samples_'//trim(basin_info(i)%obasin_name)//'_veusz.csv',&
                               status='old',position='append',action='write')
                          write(27,'(a7)') 'nan,nan'
                          do k=1,pnummc+1
                            write(27,'(e13.6,a1,e13.6)') pnmc(k),',',pecdfmc(k)
                          enddo
                          close(27)
                        endif
                      else
                        open(27,file='age_pdf_output/mc_age_ECDF_'//hc//'_'&
                             //trim(mcschar)//'_samples_'//trim(basin_info(i)%obasin_name)//'.dat',&
                             status='unknown')
                        if (params%tec_header) then
                          write(pnummcc,'(i10)') pnummc+1
                          pnummcc=adjustl(pnummcc)
                          write(27,'(a39)') 'TITLE="Monte Carlo predicted age ECDFs"'
                          write(27,'(a45)') 'VARIABLES="Age [Ma]" "Cumulative probability"'
                          write(27,'(a80)') 'ZONE I='//trim(pnummcc)//&
                                 ' DATAPACKING=POINT T="Monte Carlo predicted ECDFs '&
                                 //trim(basin_info(i)%obasin_name)//'"'
                        endif
                        do k=1,pnummc+1
                          write(27,'(e13.6,e13.6)') pnmc(k),pecdfmc(k)
                        enddo
                        close(27)
                      endif
                    ! Write out Monte Carlo cumulative density functions
                    else
                      if (params%veusz_output) then
                        if (j == 1) then
                          open(27,file='age_pdf_output/mc_age_CDFs_'//trim(mcschar)//&
                               '_samples_'//trim(basin_info(i)%obasin_name)//'_veusz.csv',&
                               status='unknown',action='write')
                          write(27,'(a31)') 'Age [Ma],Cumulative probability'
                          do k=1,pnummc+1
                            write(27,'(e13.6,a1,e13.6)') pnmc(k),',',pcdfmc(k)
                          enddo
                          close(27)
                        else
                          open(27,file='age_pdf_output/mc_age_CDFs_'//trim(mcschar)//&
                               '_samples_'//trim(basin_info(i)%obasin_name)//'_veusz.csv',&
                               status='old',position='append',action='write')
                          write(27,'(a7)') 'nan,nan'
                          do k=1,pnummc+1
                            write(27,'(e13.6,a1,e13.6)') pnmc(k),',',pcdfmc(k)
                          enddo
                          close(27)
                        endif
                      else
                        open(27,file='age_pdf_output/mc_age_CDF_'//hc//'_'&
                             //trim(mcschar)//'_samples_'//trim(basin_info(i)%obasin_name)//'.dat',&
                             status='unknown')
                        if (params%tec_header) then
                          write(pnummcc,'(i10)') pnummc+1
                          pnummcc=adjustl(pnummcc)
                          write(27,'(a38)') 'TITLE="Monte Carlo predicted age CDFs"'
                          write(27,'(a45)') 'VARIABLES="Age [Ma]" "Cumulative probability"'
                          write(27,'(a80)') 'ZONE I='//trim(pnummcc)//&
                                 ' DATAPACKING=POINT T="Monte Carlo predicted CDFs '&
                                 //trim(basin_info(i)%obasin_name)//'"'
                        endif
                        do k=1,pnummc+1
                          write(27,'(e13.6,e13.6)') pnmc(k),pcdfmc(k)
                        enddo
                        close(27)
                      endif
                    endif
                  endif
                endif
              endif

! Deallocate arrays
              deallocate(pnmc,ppdfmc,ppdfvmc)                                   ! Deallocate arrays reallocated during monte carlo sim
              if (params%mccdfs_out) then
                if (params%ecdfs) then
                  deallocate(pecdfmc)
                else
                  deallocate(pcdfmc)
                endif
              endif

              if (params%lsero) deallocate(lsage,lsageu,lserate,lseratesc)
            enddo

            write(*,'(a)') 'Done.'
            mc_iterf=real(params%mc_iter)
            if (params%datamcpdfs .or. params%ppdfmcpdfs) then
              kpct(i)=(1-(sum(kuiper_res)/mc_iterf))*100.                       ! Store percent of models that passed kuiper test for given basin
              write (*,'(a,f5.1,a)') 'Predicted Monte Carlo PDFs that passed the Kuiper test: ',kpct(i),'%'

              ! Write output files
              open(20,file='kuiper_mc_results_'//trim(mcschar)//'_samples_'&
                   //trim(basin_info(i)%obasin_name)//'.dat',status='unknown')
              do j=1,params%mc_iter
                write(20,*) kuiper_res(j)                                       ! Write out individual subset kuiper test result values to file for this basin
              enddo
              close(20)
              write(21,'(a12,a5,f10.2)') basin_info(i)%obasin_name,mcschar,kpct(i) ! Write out the summary percent of models that passed the Kuiper test
            endif

            ! Deallocate arrays
            deallocate(pagemc,pageumc,peratemc,peratescmc)
            !if (mcboth) deallocate(pagemc2,pageumc2)

! End of main loop
          enddo
        endif

! Write out PDFs/CDFs/ECDFs
            ! Write out data PDF
            if (params%opdf_out) then
              if (params%veusz_output) then
                open(22,file='age_pdf_output/data_age_PDF_'//trim(basin_info(i)%obasin_name)//&
                    '_veusz.csv',status='unknown',action='write')
                write(22,'(a20)') 'Age [Ma],Probability'
                do k=1,onum+1
                  write(22,'(e13.6,a1,e13.6)') on(k),',',opdf(k)
                enddo
                close(22)
              else
                open(22,file='age_pdf_output/data_age_PDF_'//trim(basin_info(i)%obasin_name)//&
                    '.dat',status='unknown')
                if (params%tec_header) then
                  write(onumc,'(i10)') onum+1
                  onumc=adjustl(onumc)
                  write(22,'(a20)') 'TITLE="Data age PDF"'
                  write(22,'(a34)') 'VARIABLES="Age [Ma]" "Probability"'
                  write(22,'(a60)') 'ZONE I='//trim(onumc)//&
                           ' DATAPACKING=POINT T="Data PDF '//trim(basin_info(i)%obasin_name)//'"'
                endif
                do k=1,onum+1
                  write(22,'(e13.6,e13.6)') on(k),opdf(k)
                enddo
                close(22)
              endif
            endif

            ! Write out data CDF or ECDF
            if (params%ocdf_out) then
              ! Write out data empirical cumulative distribution function
              if (params%ecdfs) then
                if (params%veusz_output) then
                  open(25,file='age_pdf_output/data_age_ECDF_'//trim(basin_info(i)%obasin_name)//&
                      '_veusz.csv',status='unknown',action='write')
                  write(25,'(a31)') 'Age [Ma],Cumulative probability'
                  do k=1,onum+1
                    write(25,'(e13.6,a1,e13.6)') on(k),',',oecdf(k)
                  enddo
                  close(25)
                else
                  open(25,file='age_pdf_output/data_age_ECDF_'//trim(basin_info(i)%obasin_name)//&
                      '.dat',status='unknown')
                  if (params%tec_header) then
                    write(onumc,'(i10)') onum+1
                    onumc=adjustl(onumc)
                    write(25,'(a21)') 'TITLE="Data age ECDF"'
                    write(25,'(a45)') 'VARIABLES="Age [Ma]" "Cumulative probability"'
                    write(25,'(a60)') 'ZONE I='//trim(onumc)//&
                             ' DATAPACKING=POINT T="Data ECDF '//trim(basin_info(i)%obasin_name)//'"'
                  endif
                  do k=1,onum+1
                    write(25,'(e13.6,e13.6)') on(k),oecdf(k)
                  enddo
                  close(25)
                endif
              ! Write out data cumulative density function
              else
                if (params%veusz_output) then
                  open(25,file='age_pdf_output/data_age_CDF_'//trim(basin_info(i)%obasin_name)//&
                      '_veusz.csv',status='unknown',action='write')
                  write(25,'(a31)') 'Age [Ma],Cumulative probabilty'
                  do k=1,onum+1
                    write(25,'(e13.6,a1,e13.6)') on(k),',',ocdf(k)
                  enddo
                  close(25)
                else
                  open(25,file='age_pdf_output/data_age_CDF_'//trim(basin_info(i)%obasin_name)//&
                      '.dat',status='unknown')
                  if (params%tec_header) then
                    write(onumc,'(i10)') onum+1
                    onumc=adjustl(onumc)
                    write(25,'(a20)') 'TITLE="Data age CDF"'
                    write(25,'(a45)') 'VARIABLES="Age [Ma]" "Cumulative probability"'
                    write(25,'(a60)') 'ZONE I='//trim(onumc)//&
                             ' DATAPACKING=POINT T="Data CDF '//trim(basin_info(i)%obasin_name)//'"'
                  endif
                  do k=1,onum+1
                    write(25,'(e13.6,e13.6)') on(k),ocdf(k)
                  enddo
                  close(25)
                endif
              endif
            endif

            ! Write out full predicted PDF
            if (params%ppdf_out) then
              if (params%veusz_output) then
                open(23,file='age_pdf_output/full_predicted_age_PDF_'&
                     //trim(basin_info(i)%obasin_name)//'_veusz.csv',          &
                     status='unknown',action='write')
                write(23,'(a20)') 'Age [Ma],Probability'
                do k=1,pnum+1
                  write(23,'(e13.6,a1,e13.6)') pn(k),',',ppdf(k)
                enddo
                close(23)
              else
                open(23,file='age_pdf_output/full_predicted_age_PDF_'&
                     //trim(basin_info(i)%obasin_name)//'.dat',status='unknown')
                if (params%tec_header) then
                  write(pnumc,'(i10)') pnum+1
                  pnumc=adjustl(pnumc)
                  write(23,'(a25)') 'TITLE="Predicted age PDF"'
                  write(23,'(a34)') 'VARIABLES="Age [Ma]" "Probability"'
                  write(23,'(a65)') 'ZONE I='//trim(pnumc)//&
                  ' DATAPACKING=POINT T="Predicted PDF '//trim(basin_info(i)%obasin_name)//'"'
                endif
                do k=1,pnum+1
                  write(23,'(e13.6,e13.6)') pn(k),ppdf(k)
                enddo
                close(23)
              endif
            endif

            ! Write out full predicted CDF or ECDF
            if (params%pcdf_out) then
              ! Write out data empirical cumulative distribution function
              if (params%ecdfs) then
                if (params%veusz_output) then
                  open(26,file='age_pdf_output/full_predicted_age_ECDF_'&
                       //trim(basin_info(i)%obasin_name)//'_veusz.csv',        &
                       status='unknown',action='write')
                  write(26,'(a31)') 'Age [Ma],Cumulative probabilty'
                  do k=1,pnum+1
                    write(26,'(e13.6,e13.6)') pn(k),pecdf(k)
                  enddo
                  close(26)
                else
                  open(26,file='age_pdf_output/full_predicted_age_ECDF_'&
                       //trim(basin_info(i)%obasin_name)//'.dat',status='unknown')
                  if (params%tec_header) then
                    write(pnumc,'(i10)') pnum+1
                    pnumc=adjustl(pnumc)
                    write(26,'(a26)') 'TITLE="Predicted age ECDF"'
                    write(26,'(a45)') 'VARIABLES="Age [Ma]" "Cumulative probability"'
                    write(26,'(a65)') 'ZONE I='//trim(pnumc)//&
                    ' DATAPACKING=POINT T="Predicted ECDF '//trim(basin_info(i)%obasin_name)//'"'
                  endif
                  do k=1,pnum+1
                    write(26,'(e13.6,e13.6)') pn(k),pecdf(k)
                  enddo
                  close(26)
                endif
              ! Write out data cumulative density function
              else
                if (params%veusz_output) then
                  open(26,file='age_pdf_output/full_predicted_age_CDF_'&
                       //trim(basin_info(i)%obasin_name)//'_veusz.csv',        &
                       status='unknown',action='write')
                  write(26,'(a31)') 'Age [Ma],Cumulative probability'
                  do k=1,pnum+1
                    write(26,'(e13.6,e13.6)') pn(k),pcdf(k)
                  enddo
                  close(26)
                else
                  open(26,file='age_pdf_output/full_predicted_age_CDF_'&
                       //trim(basin_info(i)%obasin_name)//'.dat',status='unknown')
                  if (params%tec_header) then
                    write(pnumc,'(i10)') pnum+1
                    pnumc=adjustl(pnumc)
                    write(26,'(a25)') 'TITLE="Predicted age CDF"'
                    write(26,'(a45)') 'VARIABLES="Age [Ma]" "Cumulative probability"'
                    write(26,'(a65)') 'ZONE I='//trim(pnumc)//&
                    ' DATAPACKING=POINT T="Predicted CDF '//trim(basin_info(i)%obasin_name)//'"'
                  endif
                  do k=1,pnum+1
                    write(26,'(e13.6,e13.6)') pn(k),pcdf(k)
                  enddo
                  close(26)
                endif
              endif
            endif

! Deallocate arrays
        if (params%datapdf) deallocate(oage,oageu,on,opdf,opdfv,oeratesc)
        if (params%ocdf_out) then
          if (params%ecdfs) then
            deallocate(oecdf)
          else
            deallocate(ocdf)
          endif
        endif
        if (params%fullppdf) deallocate(pn,ppdf,ppdfv)
        deallocate(page,pageu,perate,peratesc)
        if (params%datamcpdfs .or. params%ppdfmcpdfs) deallocate(kuiper_res)

! Close open files
        close(21)
      enddo

      ! Sign off
      write (*,*) ''
      write (*,'(a)') '#------------------------------------------------------------------------------#'
      write (*,'(a)') 'Execution complete'
      write (*,'(a)') '#------------------------------------------------------------------------------#'

! Exit
      end
