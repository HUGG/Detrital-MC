subroutine read_input_file(params)

use definitions

implicit none

! This subroutine reads the detrital_mc2 input file (dmc2_input.txt)

type (detrital_params) params

character(len=1024) :: line
integer :: datapdf_in,fullppdf_in,mcpdfs_in
integer :: datappdf_in,datamcpdfs_in,ppdfmcpdfs_in
integer :: opdf_out_in,ppdf_out_in,mcpdfs_out_in
integer :: lsero_in,tec_header_in,calc_pdf_range_in
integer :: i,j,k,io
logical :: fileexist
real(kind=sp) :: simyr_in

inquire(file='input/dmc2_input.txt',exist=fileexist)
if (.not.fileexist) then
  write(*,'(a)') 'Error: Cannot find input file.'
  write(*,'(a)') '       Is the input file (dmc2_input.txt) in the input directory?'
  stop
endif

open (unit=100,file='input/dmc2_input.txt',status='old')
open (unit=101,status='scratch')

write (*,'(a)') 'Reading input file...'

do
  read(unit=100,*,iostat=io) line
  if (io > 0) then
    write(*,'(a)') 'Error: Problem with input file.'
    stop
  elseif (io = 0) then
    if (line(1:1) /= '$' .and. line(1:1) /= ' ') then
      if (scan(line,'$') /= 0) then
        do i=scan(line,'$'),1024
          line(i:i)=' '
        enddo
      endif
      k=1
      do j=1,1024
        if (line(j:j) == ' ' .or. line(j:j) == ',') then
          if (j /= k) write (unit=101,'(a)') line(k:j-1)
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

![char] params%basin_summary_file is the name of file containing the summary
!info for the desired set of basins. It should be no more than 80 characters
!long
read (unit=101,*) params%basin_summary_file

![int] params%num_basins is the number of basins to consider in the detrital age
!analysis
read (unit=101,*) params%num_basins

![int] params%nss is the number of sample sizes to consider
![int] params%numsamp() is the array of params%nss sample sizes
read (unit=101,*) params%nss
if (params%nss > 0) then
  allocate (params%numsamp(params%nss))
  do i=1,params%nss
    read (unit=101,*) params%numsamp(i)
  enddo
else
  read (unit=101,*)
endif
    
![int] datapdf_in is the flag for whether or not observed age PDFs should be
!generated. This value is stored as [bool] params%datapdf
read (unit=101,*) datapdf_in
if (datapdf_in == 0) then
  params%datapdf = .false.
elseif (datapdf_in == 1) then
  params%datapdf = .true.
else
  write (*,'(a,i)') 'Error: Bad value for flag to generate observed age PDFs: ',datapdf_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![int] fullppdf_in is the flag for whether or not full predicted age PDFs should
!be generated. This value is stored as [bool] params%fullppdf
read (unit=101,*) fullppdf_in
if (fullppdf_in == 0) then
  params%fullppdf = .false.
elseif (fullppdf_in == 1) then
  params%fullppdf = .true.
else
  write (*,'(a,i)') 'Error: Bad value for flag to generate full predicted age PDFs: ',fullppdf_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![int] mcpdfs_in is the flag for whether or not Monte Carlo predicted age PDFs
!should be generated. This value is stored as [bool] params%mcpdfs
read (unit=101,*) mcpdfs_in
if (mcpdfs_in == 0) then
  params%mcpdfs = .false.
elseif (mcpdfs_in == 1) then
  params%mcpdfs = .true.
else
  write (*,'(a,i)') 'Error: Bad value for flag to generate Monte Carlo predicted age PDFs: ',mcpdfs_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![int] datappdf_in is the flag for whether or not to compare the observed age
!PDFs to the full predicted age PDFs. This value is stored as [bool]
!params%datappdf
read (unit=101,*) datappdf_in
if (datappdf_in == 0) then
  params%datappdf = .false.
elseif (datappdf_in == 1) then
  params%datappdf = .true.
else
  write (*,'(a)') 'Error: Bad value for flag to compare observed age to full predicted age PDFs: '
  write (*,'(a,i)') '       ',datappdf_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![int] datamcpdfs_in is the flag for whether or not to compare the observed age
!PDFs to the Monte Carlo predicted age PDFs. This value is stored as [bool]
!params%datamcpdfs
read (unit=101,*) datamcpdfs_in
if (datamcpdfs_in == 0) then
  params%datamcpdfs = .false.
elseif (datamcpdfs_in == 1) then
  params%datamcpdfs = .true.
else
  write (*,'(a)') 'Error: Bad value for flag to compare observed age to Monte Carlo predicted age'
  write (*,'(a,i)') '       PDFs: ',datamcpdfs_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![int] ppdfmcpdfs_in is the flag for whether or not to compare the full
!predicted age PDFs to the Monte Carlo predicted age PDFs. This value is stored
!as [bool] params%ppdfmcpdfs
read (unit=101,*) ppdfmcpdfs_in
if (ppdfmcpdfs_in == 0) then
  params%ppdfmcpdfs = .false.
elseif (ppdfmcpdfs_in == 1) then
  params%ppdfmcpdfs = .true.
else
  write (*,'(a)') 'Error: Bad value for flag to compare Monte Carlo predicted age to full'
  write (*,'(a,i)') '       predicted age PDFs: ',ppdfmcpdfs_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![int] lsero_in is the flag for whether or not to simulate landslide erosion
!This value is stored as [bool] params%lsero
read (unit=101,*) lsero_in
if (lsero_in == 0) then
  params%lsero = .false.
elseif (lsero_in == 1) then
  params%lsero = .true.
else
  write (*,'(a)') 'Error: Bad value for flag to simulate landslide erosion: ',lsero_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![flt] params%lsagejunk is the junk age to be used when no landslides occur in
!the current catchment
read (unit=101,*) params%lsagejunk

![flt] params%lseratejunk is the junk erosion rate to be used when no landslides
!occur in the current catchment
read (unit=101,*) params%lseratejunk

![flt] simyr_in is the sediment residence time in the catchment in years. This
!value is converted to a character value and stored as [char] params%simyr
read (unit=101,*) simyr_in
write (unit=params%simyr,fmt='(f8.4)') simyr_in

![int] opdf_out_in is the flag for whether or not observed age PDFs should be
!written to file(s). This value is stored as [bool] params%opdf_out
read (unit=101,*) opdf_out_in
if (opdf_out_in == 0) then
  params%opdf_out = .false.
elseif (opdf_out_in == 1) then
  params%opdf_out = .true.
else
  write (*,'(a,i)') 'Error: Bad value for flag to output observed age PDFs: ',opdf_out_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![int] ppdf_out_in is the flag for whether or not full predicted age PDFs should
!be written to file(s). This value is stored as [bool] params%ppdf_out
read (unit=101,*) ppdf_out_in
if (ppdf_out_in == 0) then
  params%ppdf_out = .false.
elseif (ppdf_out_in == 1) then
  params%ppdf_out = .true.
else
  write (*,'(a,i)') 'Error: Bad value for flag to output full predicted age PDFs: ',ppdf_out_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![int] mcpdfs_out_in is the flag for whether or not Monte Carlo predicted age
!PDFs should be written to file(s). This value is stored as [bool]
!params%mcpdfs_out
read (unit=101,*) mcpdfs_out_in
if (mcpdfs_out_in == 0) then
  params%mcpdfs_out = .false.
elseif (mcpdfs_out_in == 1) then
  params%mcpdfs_out = .true.
else
  write (*,'(a,i)') 'Error: Bad value for flag to output Monte Carlo predicted age PDFs: ',mcpdfs_out_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![int] params%num_mc_out is the number of Monte Carlo predicted age PDFs to
!write to file
read (unit=101,*) params%num_mc_out

![int] tec_header_in is the flag for whether or not to include a Tecplot header
!in the PDF output files. This value is stored as [bool] params%tec_header
read (unit=101,*) tec_header_in
if (tec_header_in == 0) then
  params%tec_header = .false.
elseif (tec_header_in == 1) then
  params%tec_header = .true.
else
  write (*,'(a,i)') 'Error: Bad value for flag to include Tecplot headers in the output files: ',tec_header_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![int] params%mc_iter is the number of Monte Carlo iterations to simulate
read (unit=101,*) params%mc_iter

![flt] params%dx is the age interval for PDF value calculation
read (unit=101,*) params%dx

![int] calc_pdf_range_in is the flag for whether or not the age range for the
!PDFs should be calculated based on the observed age range and uncertainties.
!This value is stored as [bool] params%calc_pdf_range
read (unit=101,*) calc_pdf_range_in
if (calc_pdf_range_in == 0) then
  params%calc_pdf_range = .false.
elseif (calc_pdf_range_in == 1) then
  params%calc_pdf_range = .true.
else
  write (*,'(a,i)') 'Error: Bad value for flag to calculate the PDF age range: ',calc_pdf_range_in
  write (*,'(a)')   '       Value must be either "0" or "1"'
  stop
endif

![flt] params%pdfmin is the minimum age for the PDF calculations
read (unit=101,*) params%pdfmin

![flt] params%pdfmax is the maximum age for the PDF calculations
read (unit=101,*) params%pdfmax

![flt] params%pdf_pct_uncert is percent age uncertainty applied to the predicted
!ages when not making observed age PDFs
read (unit=101,*) params%pdf_pct_uncert

![flt] params%pdfscl is the approximate number of values in scaled age PDFs
read (unit=101,*) params%pdf_pct_uncert

close(unit=101)

end