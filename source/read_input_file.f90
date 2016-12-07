subroutine read_input_file(params,basin_info)

use definitions

implicit none

! This subroutine reads the detrital_mc v2.0 input file (det_mc_input.txt)

type (detrital_params) :: params
type (basin_information) :: basin_info(100)

character(len=80) :: line
integer :: datapdf_in,fullppdf_in,mcpdfs_in,ecdfs_in
integer :: datappdf_in,datamcpdfs_in,ppdfmcpdfs_in
integer :: opdf_out_in,ppdf_out_in,mcpdfs_out_in
integer :: ocdf_out_in,pcdf_out_in,mccdfs_out_in
integer :: lsero_in,tec_header_in,calc_pdf_range_in
integer :: kuipernew_in,scale_erates_in
integer :: i,j,k,io,geol_units,veusz_output_in
logical :: fileexist,echo_vals
real(kind=sp) :: simyr_in

echo_vals=.true.

inquire(file='input/det_mc_input.txt',exist=fileexist)
if (.not.fileexist) then
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a)') 'Error: Cannot find input file.'
  write (*,'(a)') '       Is the input file (det_mc_input.txt) in the input directory?'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  stop
endif

open (unit=100,file='input/det_mc_input.txt',status='old')
open (unit=101,status='scratch')

write (*,'(a)') '#------------------------------------------------------------------------------#'
write (*,*) ''
write (*,'(a)', advance='no') 'Reading input file... '

do
  read(unit=100,fmt='(a)',iostat=io) line
  line=trim(line)

  if (io > 0) then
    write (*,'(a)') '#------------------------------------------------------------------------------#'
    write (*,'(a,i3,a)') 'Error: Problem with input file (iostat=',io,').'
    write (*,'(a)') ''
    write (*,'(a)') 'Program exited with an error'
    write (*,'(a)') '#------------------------------------------------------------------------------#'
    close(unit=100)
    close(unit=101)
    stop
  elseif (io == 0) then
    if (line(1:1) /= '$' .and. line(1:1) /= ' ') then
      if (scan(line,'$') /= 0) then
        do i=scan(line,'$'),1024
          line(i:i)=' '
        enddo
      endif
      k=1
      do j=1,len(line)
        if (line(j:j) == ' ' .or. line(j:j) == ',') then
          if (j /= k) write (unit=101,fmt='(a)') line(k:j-1)
          k=j+1
        endif
      enddo
    endif
  else
    ! End of file reached, break of of loop
    write (*,'(a)') 'Done.'
    exit
  endif
enddo

close (unit=100)
rewind(unit=101)

![int] params%num_basins is the number of basins to consider in the detrital age
!analysis
read (unit=101,fmt=*) params%num_basins
if (echo_vals) write (*,*) 'params%num_basins: ',params%num_basins
if (params%num_basins == 0) then
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a)') 'Error: Zero basins to analyze (num_basins = 0), nothing to do.'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
elseif (params%num_basins < 0) then
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a)') 'Error: Negative number of basins to analyze (num_basins < 0), nothing to do.'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
elseif (params%num_basins > 100) then
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a)') 'Error: Too many basins to analyze (num_basins > 100).'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
else
  do i=1,params%num_basins
    !allocate (basin_info(params%num_basins))
    ![char] basin_info(:)%obasin_name is the name of the observed age file
    read (unit=101,fmt=*) basin_info(i)%obasin_name
    if (echo_vals) write (*,*) 'basin_info(',i,')%obasin_name: ',trim(basin_info(i)%obasin_name)
    inquire(file='data/observed_ages/'//trim(basin_info(i)%obasin_name)//'.dat',exist=fileexist)
    if (.not.fileexist) then
      write (*,'(a)') '#------------------------------------------------------------------------------#'
      write (*,'(a)') 'Error: Cannot find Pecube observed age file.'
      write (*,'(a)') '       Is the age file '//trim(basin_info(i)%obasin_name)//'.dat in the data/observed_ages/ directory?'
      write (*,'(a)') ''
      write (*,'(a)') 'Program exited with an error'
      write (*,'(a)') '#------------------------------------------------------------------------------#'
      close(unit=101)
      stop
    endif

    ![int] basin_info(:)%page_ftype is the number corresponding to the input predicted
    !age file format
    !1 = Pecube Comparison.txt, 2 = Generic CSV file
    read (unit=101,fmt=*) basin_info(i)%page_ftype
    if (echo_vals) write (*,*) 'basin_info(',i,')%page_ftype: ',basin_info(i)%page_ftype
    if (basin_info(i)%page_ftype < 1 .or. basin_info(i)%page_ftype > 2) then
      write (*,'(a)') '#------------------------------------------------------------------------------#'
      write (*,'(a,i1,a,i3)') 'Error: Unsupported predicted age file format (',basin_info(i)%page_ftype,') for basin ',i,'.'
      write (*,'(a)') '       The listed predicted age file type must be an integer value between 1 and 2'
      write (*,'(a)') ''
      write (*,'(a)') 'Program exited with an error'
      write (*,'(a)') '#------------------------------------------------------------------------------#'
      close(unit=101)
      stop
    endif

    if (basin_info(i)%page_ftype == 1) then
      ![char] basin_info(:)%pbasin_name is the name of the predicted age output
      !directory
      read (unit=101,fmt=*) basin_info(i)%pbasin_name
      if (echo_vals) write (*,*) 'basin_info(',i,')%pbasin_name: ',trim(basin_info(i)%pbasin_name)
      inquire(file='data/predicted_ages/'//trim(basin_info(i)%pbasin_name)//'/Comparison.txt',exist=fileexist)
      if (.not.fileexist) then
        write (*,'(a)') '#------------------------------------------------------------------------------#'
        write (*,'(a)') 'Error: Cannot find Pecube predicted age file.'
        write (*,'(a)') '       Does the data/predicted_ages/'//trim(basin_info(i)%pbasin_name)//'/ directory exist?'
        write (*,'(a)') '       Is the age file (Comparison.txt) in the data/predicted_ages/'&
                                //trim(basin_info(i)%pbasin_name)//' directory?'
        write (*,'(a)') ''
        write (*,'(a)') 'Program exited with an error'
        write (*,'(a)') '#------------------------------------------------------------------------------#'
        close(unit=101)
        stop
      endif

      ![int] basin_info(:)%page_sys is the number corresponding to the desired predicted
      !thermochronometer age system to be compared to the data
      !1 = AHe, 2 = AFT, 3 = ZHe, 4 = ZFT, 5 = MAr
      read (unit=101,fmt=*) basin_info(i)%page_sys
      if (echo_vals) write (*,*) 'basin_info(',i,')%page_sys: ',basin_info(i)%page_sys
      if (basin_info(i)%page_sys < 1 .or. basin_info(i)%page_sys > 5) then
        write (*,'(a)') '#------------------------------------------------------------------------------#'
        write (*,'(a,i1,a,i3)') 'Error: Unsupported thermochronometer system (',basin_info(i)%page_sys,') for basin ',i,'.'
        write (*,'(a)') '       The listed thermochronometer system must be an integer value between 1 and 5'
        write (*,'(a)') ''
        write (*,'(a)') 'Program exited with an error'
        write (*,'(a)') '#------------------------------------------------------------------------------#'
        close(unit=101)
        stop
      endif
    else
      ![char] basin_info(:)%pbasin_name is the name of the predicted age input
      !CSV file
      read (unit=101,fmt=*) basin_info(i)%pbasin_name
      if (echo_vals) write (*,*) 'basin_info(',i,')%pbasin_name: ',trim(basin_info(i)%pbasin_name)
      inquire(file='data/predicted_ages/'//trim(basin_info(i)%pbasin_name)//'.csv',exist=fileexist)
      if (.not.fileexist) then
        write (*,'(a)') '#------------------------------------------------------------------------------#'
        write (*,'(a)') 'Error: Cannot find predicted age file.'
        write (*,'(a)') '       Does the file data/predicted_ages/'//trim(basin_info(i)%pbasin_name)//'.csv exist?'
        write (*,'(a)') ''
        write (*,'(a)') 'Program exited with an error'
        write (*,'(a)') '#------------------------------------------------------------------------------#'
        close(unit=101)
        stop
      endif

      ![int] basin_info(:)%page_col is the number corresponding to the column in
      !the CSV file containing the predicted ages
      read (unit=101,fmt=*) basin_info(i)%page_col
      if (echo_vals) write (*,*) 'basin_info(',i,')%page_col: ',basin_info(i)%page_col
      if (basin_info(i)%page_col < 1 .or. basin_info(i)%page_col > 1000) then
        write (*,'(a)') '#------------------------------------------------------------------------------#'
        write (*,'(a,i1,a,i3)') 'Error: Unsupported CSV file age column value (',basin_info(i)%page_col,') for basin ',i,'.'
        write (*,'(a)') '       The listed value must be a positive integer value less than 1000'
        write (*,'(a)') ''
        write (*,'(a)') 'Program exited with an error'
        write (*,'(a)') '#------------------------------------------------------------------------------#'
        close(unit=101)
        stop
      endif

      ![int] basin_info(:)%perate_col is the number corresponding to the column
      !in the CSV file containing the predicted erosion rates
      read (unit=101,fmt=*) basin_info(i)%perate_col
      if (echo_vals) write (*,*) 'basin_info(',i,')%perate_col: ',basin_info(i)%perate_col
      if (basin_info(i)%perate_col < 1 .or. basin_info(i)%perate_col > 1000) then
        write (*,'(a)') '#------------------------------------------------------------------------------#'
        write (*,'(a,i1,a,i3)') 'Error: Unsupported CSV file erate column value (',basin_info(i)%perate_col,') for basin ',i,'.'
        write (*,'(a)') '       The listed value must be a positive integer value less than 1000'
        write (*,'(a)') ''
        write (*,'(a)') 'Program exited with an error'
        write (*,'(a)') '#------------------------------------------------------------------------------#'
        close(unit=101)
        stop
      endif

      if (basin_info(i)%perate_col==16 .or. basin_info(i)%perate_col==17 .or.  &
          basin_info(i)%perate_col==18 .or. basin_info(i)%perate_col==19 .or.  &
          basin_info(i)%perate_col==98 .or. basin_info(i)%perate_col==99) then
        allocate(basin_info(i)%geol_scale_factor(10))
        basin_info(i)%geol_scale_factor=0.0
        if (basin_info(i)%perate_col==16) geol_units = 6
        if (basin_info(i)%perate_col==17) geol_units = 2
        if (basin_info(i)%perate_col==18) geol_units = 2
        if (basin_info(i)%perate_col==19) geol_units = 2
        if (basin_info(i)%perate_col==98) geol_units = 10
        if (basin_info(i)%perate_col==99) geol_units = 4
        do j=1,geol_units
          read (unit=101,fmt=*) basin_info(i)%geol_scale_factor(j)
          if (echo_vals) write (*,*) 'basin_info(',i,')%geol_scale_factor(',j,'): ',basin_info(i)%geol_scale_factor(j)
          if (basin_info(i)%geol_scale_factor(j) < 0) then
            write (*,'(a)') '#------------------------------------------------------------------------------#'
            write (*,'(a,i1,a,i3)') 'Error: Unsupported value for erosion rate scaling (',basin_info(i)%geol_scale_factor(j),') for basin ',i,'.'
            write (*,'(a)') '       The listed value must be positive or zero'
            write (*,'(a)') ''
            write (*,'(a)') 'Program exited with an error'
            write (*,'(a)') '#------------------------------------------------------------------------------#'
            close(unit=101)
            stop
          endif
        enddo
        read (unit=101,fmt=*) basin_info(i)%uplift_velo_scaling
        if (echo_vals) write (*,*) 'basin_info(',i,')%uplift_velo_scaling: ',basin_info(i)%uplift_velo_scaling
        if (basin_info(i)%uplift_velo_scaling < 0 .or. basin_info(i)%uplift_velo_scaling > 3) then
          write (*,'(a)') '#------------------------------------------------------------------------------#'
          write (*,'(a,i1)') 'Error: Bad value for type of uplift scaling to apply: ',basin_info(i)%uplift_velo_scaling
          write (*,'(a)') '       Value must be either "0", "1", "2", or "3"'
          write (*,'(a)') ''
          write (*,'(a)') 'Program exited with an error'
          write (*,'(a)') '#------------------------------------------------------------------------------#'
          close(unit=101)
          stop
        endif
      endif
    endif
  enddo
endif

!![char] params%basin_summary_file is the name of file containing the summary
!!info for the desired set of basins. It should be no more than 80 characters
!!long
!read (unit=101,fmt=*) params%basin_summary_file
!if (echo_vals) write (*,*) 'params%basin_summary_file: ',trim(params%basin_summary_file)

!![char] params%model_name is the name of output directory containing the Pecube
!!Comparison.txt file. It should be no more than 5 characters long.
!read (unit=101,fmt=*) params%model_name
!if (echo_vals) write (*,*) 'params%model_name: ',trim(params%model_name)

!inquire(file='data/'//trim(params%model_name)//'/Comparison.txt',exist=fileexist)
!if (.not.fileexist) then
!  write (*,'(a)') '#------------------------------------------------------------------------------#'
!  write (*,'(a)') 'Error: Cannot find Pecube predicted age file.'
!  write (*,'(a)') '       Does the data/'//trim(params%model_name)//'/ directory exist?'
!  write (*,'(a)') '       Is the age file (Comparison.txt) in the data/'//trim(params%model_name)//' directory?'
!  write (*,'(a)') ''
!  write (*,'(a)') 'Program exited with an error'
!  write (*,'(a)') '#------------------------------------------------------------------------------#'
!  stop
!endif

!![int] params%page_sys is the number corresponding to the desired predicted
!!thermochronometer age system to be compared to the data
!!1 = AHe, 2 = AFT, 3 = ZHe, 4 = ZFT, 5 = MAr
!read (unit=101,fmt=*) params%page_sys
!if (echo_vals) write (*,*) 'params%page_sys: ',params%page_sys

![int] params%nss is the number of sample sizes to consider
![int] params%numsamp() is the array of params%nss sample sizes
read (unit=101,fmt=*) params%nss
if (echo_vals) write (*,*) 'params%nss: ',params%nss
if (params%nss > 0) then
  allocate (params%numsamp(params%nss))
  do i=1,params%nss
    read (unit=101,fmt=*) params%numsamp(i)
    if (echo_vals) write (*,*) 'params%numsamp(i): ',params%numsamp(i)
  enddo
else
  read (unit=101,fmt=*)
endif

![int] datapdf_in is the flag for whether or not observed age PDFs should be
!generated. This value is stored as [bool] params%datapdf
read (unit=101,fmt=*) datapdf_in
if (echo_vals) write (*,*) 'datapdf_in: ',datapdf_in
if (datapdf_in == 0) then
  params%datapdf = .false.
elseif (datapdf_in == 1) then
  params%datapdf = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to generate observed age PDFs: ',datapdf_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] fullppdf_in is the flag for whether or not full predicted age PDFs should
!be generated. This value is stored as [bool] params%fullppdf
read (unit=101,fmt=*) fullppdf_in
if (echo_vals) write (*,*) 'fullppdf_in: ',fullppdf_in
if (fullppdf_in == 0) then
  params%fullppdf = .false.
elseif (fullppdf_in == 1) then
  params%fullppdf = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to generate full predicted age PDFs: ',fullppdf_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] mcpdfs_in is the flag for whether or not Monte Carlo predicted age PDFs
!should be generated. This value is stored as [bool] params%mcpdfs
read (unit=101,fmt=*) mcpdfs_in
if (echo_vals) write (*,*) 'mcpdfs_in: ',mcpdfs_in
if (mcpdfs_in == 0) then
  params%mcpdfs = .false.
elseif (mcpdfs_in == 1) then
  params%mcpdfs = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to generate Monte Carlo predicted age PDFs: ',mcpdfs_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] datappdf_in is the flag for whether or not to compare the observed age
!PDFs to the full predicted age PDFs. This value is stored as [bool]
!params%datappdf
read (unit=101,fmt=*) datappdf_in
if (echo_vals) write (*,*) 'datappdf_in: ',datappdf_in
if (datappdf_in == 0) then
  params%datappdf = .false.
elseif (datappdf_in == 1) then
  params%datappdf = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a)')   'Error: Bad value for flag to compare observed age to full predicted age PDFs:'
  write (*,'(a,i1)') '       ',datappdf_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] datamcpdfs_in is the flag for whether or not to compare the observed age
!PDFs to the Monte Carlo predicted age PDFs. This value is stored as [bool]
!params%datamcpdfs
read (unit=101,fmt=*) datamcpdfs_in
if (echo_vals) write (*,*) 'datamcpdfs_in: ',datamcpdfs_in
if (datamcpdfs_in == 0) then
  params%datamcpdfs = .false.
elseif (datamcpdfs_in == 1) then
  params%datamcpdfs = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a)') 'Error: Bad value for flag to compare observed age to Monte Carlo predicted age'
  write (*,'(a,i1)') '       PDFs: ',datamcpdfs_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] ppdfmcpdfs_in is the flag for whether or not to compare the full
!predicted age PDFs to the Monte Carlo predicted age PDFs. This value is stored
!as [bool] params%ppdfmcpdfs
read (unit=101,fmt=*) ppdfmcpdfs_in
if (echo_vals) write (*,*) 'ppdfmcpdfs_in: ',ppdfmcpdfs_in
if (ppdfmcpdfs_in == 0) then
  params%ppdfmcpdfs = .false.
elseif (ppdfmcpdfs_in == 1) then
  params%ppdfmcpdfs = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a)') 'Error: Bad value for flag to compare Monte Carlo predicted age to full'
  write (*,'(a,i1)') '       predicted age PDFs: ',ppdfmcpdfs_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] ecdfs_in is the flag for whether or not to use empirical cumulative
!distribution functions (ECDFs) rather that CSPDFs for the comparison between
!the predicted and observed age distributions. This value is stored as [bool]
!params%ecdfs
read (unit=101,fmt=*) ecdfs_in
if (echo_vals) write (*,*) 'ecdfs_in: ',ecdfs_in
if (ecdfs_in == 0) then
  params%ecdfs = .false.
elseif (ecdfs_in == 1) then
  params%ecdfs = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a)') 'Error: Bad value for flag to use ECDFs rather than CSPDFs for the comparison'
  write (*,'(a,i1)') 'of the observed and predicted age distributions: ',ecdfs_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] kuipernew_in is the flag for whether or not to use the new formulation of
!Kuiper's test, which is required for comparisons using ECDFs. This value is
!stored as [bool] params%kuipernew
read (unit=101,fmt=*) kuipernew_in
if (echo_vals) write (*,*) 'kuipernew_in: ',kuipernew_in
if (kuipernew_in == 0) then
  params%kuipernew = .false.
elseif (kuipernew_in == 1) then
  params%kuipernew = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to use new formulation of Kuipers test: ',kuipernew_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![flt] params%kalpha is the significance level for Kuiper's test
read (unit=101,fmt=*) params%kalpha
if (echo_vals) write (*,*) 'params%kalpha: ',params%kalpha

![int] lsero_in is the flag for whether or not to simulate landslide erosion
!This value is stored as [bool] params%lsero
read (unit=101,fmt=*) lsero_in
if (echo_vals) write (*,*) 'lsero_in: ',lsero_in
if (lsero_in == 0) then
  params%lsero = .false.
elseif (lsero_in == 1) then
  params%lsero = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to simulate landslide erosion: ',lsero_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![flt] params%lsagejunk is the junk age to be used when no landslides occur in
!the current catchment
read (unit=101,fmt=*) params%lsagejunk
if (echo_vals) write (*,*) 'params%lsagejunk: ',params%lsagejunk

![flt] params%lseratejunk is the junk erosion rate to be used when no landslides
!occur in the current catchment
read (unit=101,fmt=*) params%lseratejunk
if (echo_vals) write (*,*) 'params%lseratejunk: ',params%lseratejunk

![flt] simyr_in is the sediment residence time in the catchment in years. This
!value is converted to a character value and stored as [char] params%simyr
read (unit=101,fmt=*) simyr_in
if (echo_vals) write (*,*) 'simyr_in: ',simyr_in
write (unit=params%simyr,fmt='(f8.4)') simyr_in

![int] params%lsfiletype is the type of file used for the landslide input files
read (unit=101,fmt=*) params%lsfiletype
if (echo_vals) write (*,*) 'params%lsfiletype: ',params%lsfiletype
if (params%scaletype > 2 .or. params%scaletype < 0) then
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for landslide file type: ',params%lsfiletype
  write (*,'(a)') '       Value must be "1" or "2"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] opdf_out_in is the flag for whether or not observed age PDFs should be
!written to file(s). This value is stored as [bool] params%opdf_out
read (unit=101,fmt=*) opdf_out_in
if (echo_vals) write (*,*) 'opdf_out_in: ',opdf_out_in
if (opdf_out_in == 0) then
  params%opdf_out = .false.
elseif (opdf_out_in == 1) then
  params%opdf_out = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to output observed age PDFs: ',opdf_out_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] ppdf_out_in is the flag for whether or not full predicted age PDFs should
!be written to file(s). This value is stored as [bool] params%ppdf_out
read (unit=101,fmt=*) ppdf_out_in
if (echo_vals) write (*,*) 'ppdf_out_in: ',ppdf_out_in
if (ppdf_out_in == 0) then
  params%ppdf_out = .false.
elseif (ppdf_out_in == 1) then
  params%ppdf_out = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to output full predicted age PDFs: ',ppdf_out_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] mcpdfs_out_in is the flag for whether or not Monte Carlo predicted age
!PDFs should be written to file(s). This value is stored as [bool]
!params%mcpdfs_out
read (unit=101,fmt=*) mcpdfs_out_in
if (echo_vals) write (*,*) 'mcpdfs_out_in: ',mcpdfs_out_in
if (mcpdfs_out_in == 0) then
  params%mcpdfs_out = .false.
elseif (mcpdfs_out_in == 1) then
  params%mcpdfs_out = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to output Monte Carlo predicted age PDFs: ',mcpdfs_out_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] params%num_mc_out is the number of Monte Carlo predicted age PDFs to
!write to file
read (unit=101,fmt=*) params%num_mc_out
if (echo_vals) write (*,*) 'params%num_mc_out: ',params%num_mc_out

![int] ocdf_out_in is the flag for whether or not observed age CDFs/ECDFs should
!be written to file(s). This value is stored as [bool] params%ocdf_out
read (unit=101,fmt=*) ocdf_out_in
if (echo_vals) write (*,*) 'ocdf_out_in: ',ocdf_out_in
if (ocdf_out_in == 0) then
  params%ocdf_out = .false.
elseif (ocdf_out_in == 1) then
  params%ocdf_out = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to output observed age CDFs/ECDFs: ',ocdf_out_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] pcdf_out_in is the flag for whether or not full predicted age CDFs/ECDFs
!should be written to file(s). This value is stored as [bool] params%pcdf_out
read (unit=101,fmt=*) pcdf_out_in
if (echo_vals) write (*,*) 'pcdf_out_in: ',pcdf_out_in
if (pcdf_out_in == 0) then
  params%pcdf_out = .false.
elseif (pcdf_out_in == 1) then
  params%pcdf_out = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to output full predicted age CDFs/ECDFs: ',pcdf_out_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] mccdfs_out_in is the flag for whether or not Monte Carlo predicted age
!CDFs/ECDFs should be written to file(s). This value is stored as [bool]
!params%mccdfs_out
read (unit=101,fmt=*) mccdfs_out_in
if (echo_vals) write (*,*) 'mccdfs_out_in: ',mccdfs_out_in
if (mccdfs_out_in == 0) then
  params%mccdfs_out = .false.
elseif (mccdfs_out_in == 1) then
  params%mccdfs_out = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to output Monte Carlo predicted age CDFs/ECDFs: ',mccdfs_out_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] tec_header_in is the flag for whether or not to include a Tecplot header
!in the PDF output files. This value is stored as [bool] params%tec_header
read (unit=101,fmt=*) tec_header_in
if (echo_vals) write (*,*) 'tec_header_in: ',tec_header_in
if (tec_header_in == 0) then
  params%tec_header = .false.
elseif (tec_header_in == 1) then
  params%tec_header = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to include Tecplot headers in the output files: ',tec_header_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] veusz_output_in is the flag for whether or not to write output PDF/CDF/
!ECDF files in a file format for use with the Veusz plotting software
read (unit=101,fmt=*) veusz_output_in
if (echo_vals) write (*,*) 'veusz_output_in: ',veusz_output_in
if (veusz_output_in == 0) then
  params%veusz_output = .false.
elseif (veusz_output_in == 1) then
  params%veusz_output = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to write output in Veusz format: ',veusz_output_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] params%mc_iter is the number of Monte Carlo iterations to simulate
read (unit=101,fmt=*) params%mc_iter
if (echo_vals) write (*,*) 'params%mc_iter: ',params%mc_iter

![flt] params%dx is the age interval for PDF value calculation
read (unit=101,fmt=*) params%dx
if (echo_vals) write (*,*) 'params%dx: ',params%dx

![int] calc_pdf_range_in is the flag for whether or not the age range for the
!PDFs should be calculated based on the observed age range and uncertainties.
!This value is stored as [bool] params%calc_pdf_range
read (unit=101,fmt=*) calc_pdf_range_in
if (echo_vals) write (*,*) 'calc_pdf_range_in: ',calc_pdf_range_in
if (calc_pdf_range_in == 0) then
  params%calc_pdf_range = .false.
elseif (calc_pdf_range_in == 1) then
  params%calc_pdf_range = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to calculate the PDF age range: ',calc_pdf_range_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![flt] params%pdfmin is the minimum age for the PDF calculations
read (unit=101,fmt=*) params%pdfmin
if (echo_vals) write (*,*) 'params%pdfmin: ',params%pdfmin

![flt] params%pdfmax is the maximum age for the PDF calculations
read (unit=101,fmt=*) params%pdfmax
if (echo_vals) write (*,*) 'params%pdfmax: ',params%pdfmax

![int] params%obs_uncert_type is the type of age uncertainty to be used when
!comparing to data (1 = mean observed uncertainty, 2 = median observed
!uncertainty)
read (unit=101,fmt=*) params%obs_uncert_type
if (echo_vals) write (*,*) 'params%obs_uncert_type: ',params%obs_uncert_type
if (params%obs_uncert_type > 4 .or. params%obs_uncert_type < 0) then
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for PDF uncertainty type: ',params%obs_uncert_type
  write (*,'(a)') '       Value must be either "1", "2", "3", or "4"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![flt] params%pdf_pct_uncert is percent age uncertainty applied to the predicted
!ages when not making observed age PDFs
read (unit=101,fmt=*) params%pdf_pct_uncert
if (echo_vals) write (*,*) 'params%pdf_pct_uncert: ',params%pdf_pct_uncert

![flt] params%pdfscl is the approximate number of values in scaled age PDFs
read (unit=101,fmt=*) params%pdfscl
if (echo_vals) write (*,*) 'params%pdfscl: ',params%pdfscl

![flt] params%alpha is the standard deviation scaling factor for making the PDFs
read (unit=101,fmt=*) params%alphain
if (echo_vals) write (*,*) 'params%alphain: ',params%alphain

![int] scale_erates_in is the flag for whether or not the input erosion rates
!should be scaled
!This value is stored as [bool] params%scale_erates
read (unit=101,fmt=*) scale_erates_in
if (echo_vals) write (*,*) 'scale_erates_in: ',scale_erates_in
if (scale_erates_in == 0) then
  params%scale_erates = .false.
elseif (scale_erates_in == 1) then
  params%scale_erates = .true.
else
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for flag to scale erosion rates: ',scale_erates_in
  write (*,'(a)') '       Value must be either "0" or "1"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] params%scaletype is the type of scaling used when scaling the input
!erosion rates
read (unit=101,fmt=*) params%scaletype
if (echo_vals) write (*,*) 'params%scaletype: ',params%scaletype
if (params%scaletype > 3 .or. params%scaletype < 0) then
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a,i1)') 'Error: Bad value for erosion rate scaling type: ',params%scaletype
  write (*,'(a)') '       Value must be "1", "2" or "3"'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

![int] params%dist_size is the desired size of age distribution arrays
read (unit=101,fmt=*) params%dist_size
if (echo_vals) write (*,*) 'params%dist_size: ',params%dist_size

!Check for compatibility of input file options
if(.not.params%kuipernew .and. params%ecdfs) then
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a)') 'Error: Incompatible options for CDF comparison and Kuipers test'
  write (*,'(a)') ''
  write (*,'(a)') '       Comparison of ECDFs requires the new version of Kuipers test'
  write (*,'(a)') '       Either use the option for CDF comparison (d in section 4 of input file)'
  write (*,'(a)') '       or use the new version of Kuipers test (e in section 4 of input file)'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

if (params%kuipernew .and. params%calc_pdf_range) then
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  write (*,'(a)') 'Error: Incompatible options for PDF age range and Kuipers test'
  write (*,'(a)') ''
  write (*,'(a)') '       The new version of Kuipers test currently requires use of a fixed age'
  write (*,'(a)') '       range. Disable option (c) in section 7 of the input file and specify an'
  write (*,'(a)') '       age range using options (d) and (e) in section 7.'
  write (*,'(a)') ''
  write (*,'(a)') 'Program exited with an error'
  write (*,'(a)') '#------------------------------------------------------------------------------#'
  close(unit=101)
  stop
endif

close(unit=101)

end
