!> example of how the glimmer memory management routines are used
!! \author Magnus Hagdorn

! to start with the preprocessor file containing the memory management
! macros needs to be loaded
#include "glimmer_memory.inc"

program handlememory
  ! load the module containing the logging routines
  use glimmer_log, only : glimmer_allocErr, glimmer_deallocErr
  implicit none

  ! an example allocatable array
  real, dimension(:), allocatable :: a
  ! an example pointer
  real, dimension(:,:), pointer :: p

  ! need to define the flag which records the status
  integer merr

  ! allocate a 1D array of size 10
  GLIMMER_ALLOC1D(a,10)
  ! allocate a 2D array of shape (/10,2/)
  GLIMMER_ALLOC2D(p,10,2)

  ! and deallocate them again
  GLIMMER_DEALLOC(a)
  GLIMMER_DEALLOC(p)
end program handlememory
