module definitions

! This module contains parameter values read from the input file for
! detrital_mc2 (dmc2_input.txt)

use nrutil

  type detrital_params
    character(len=8)  :: simyr
    character(len=80) :: basin_summary_file
    integer,dimension(:),pointer :: numsamp
    integer :: num_basins,nss,num_mc_out,mc_iter
    logical :: datapdf,fullppdf,mcpdfs,datappdf,datamcpdfs,ppdfmcpdfs
    logical :: opdf_out,ppdf_out,mcpdfs_out,lsero,tec_header,calc_pdf_range
    real(kind=SP) :: lsagejunk,lseratejunk,dx,pdfmin,pdfmax,pdf_pct_uncert
    real(kind=SP) :: pdfscl
  end type detrital_params

end module definitions