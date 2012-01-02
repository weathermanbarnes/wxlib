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
end module
