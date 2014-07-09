module definitions

! This module contains parameter values read from the input file for
! detrital_mc v2.0 (det_mc_input.txt)

  ! Variable declaration
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)

  type detrital_params
    character(len=8)  :: simyr
    integer,dimension(:),pointer :: numsamp,basin_numbers
    integer :: num_basins,nss,num_mc_out,mc_iter,obs_uncert_type
    logical :: datapdf,fullppdf,mcpdfs,datappdf,datamcpdfs,ppdfmcpdfs
    logical :: opdf_out,ppdf_out,mcpdfs_out,lsero,tec_header,calc_pdf_range
    real(kind=sp) :: lsagejunk,lseratejunk,dx,pdfmin,pdfmax,pdf_pct_uncert
    real(kind=sp) :: pdfscl
  end type detrital_params

  type basin_information
    character(len=80) :: obasin_name
    character(len=5)  :: pbasin_name
    integer :: page_sys
  end type basin_information

end module definitions
