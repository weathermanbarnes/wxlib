! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- configuration variables for all modules
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module config
  use kind
  !
  logical :: grid_cyclic_ew = .false.
  !
  ! Front / Line / Axis detection settings
  real(kind=nr) :: frint_thres =-11.7e-11_nr, &
          &        frspd_thres = 1.5_nr, &
          &          vor_thres = 0.7e-4_nr, &
          &          div_thres =-0.25e-4_nr, &
          &          def_thres = 0.7e-4_nr, &
          &       jetint_thres = 2.5e-4_nr, & 
          &          searchrad = 1.5_nr, & 
          &        smooth_coef = 0.25_nr, &
          &             minlen = 1.0e6
  integer(kind=ni) :: nsmooth = 2_ni 
end module
