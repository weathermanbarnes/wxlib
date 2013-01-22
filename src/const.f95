! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- some constants
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module consts
  use kind
  !
  character(len=48), parameter :: version  = 'dynlib ####VERREVI####'
  character(len=48), parameter :: verdate  = '####VERDATE####'
  character(len=48), parameter :: compiler = '####FCOMPIL####'
  character(len=48), parameter :: comptime = '####FCOMTIM####'
  character(len=48), parameter :: comphost = '####FCOMHOS####'
  !
  ! Some thermodynamic constants
  real(kind=nr), parameter :: cp = 1004.0_nr     ! specific heat constant of dry air for constant pressure
  real(kind=nr), parameter :: Rl = 287.04_nr     ! Gas constant for dry air
  !
  ! Other constants
  real(kind=nr), parameter :: p0 = 100000.0_nr   ! Reference pressure for T <-> Theta conversion
  real(kind=nr), parameter :: omega = 7.292123516990375e-05 ! angular velocity of the earth's rotation
  !
  ! Mathematical constants
  real(kind=nr), parameter :: pi = 3.14159265359_nr ! should be clear
  !
end module
