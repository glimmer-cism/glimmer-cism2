! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_vertcoord.f90 - part of the Glimmer-CISM ice model + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004-9 Glimmer-CISM contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! The Glimmer-CISM maintainer is:
!
! Ian Rutt
! School of the Environment and Society
! Swansea University
! Singleton Park
! Swansea
! SA2 8PP
! UK
!
! email: <i.c.rutt@swansea.ac.uk> or <ian.rutt@physics.org>
!
! Glimmer-CISM is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glimmer_memory.inc"

!> module for handling vertical sigma coordinate system
!!
!! \author Ian Rutt
!! \date September 2009
module glimmer_vertcoord

  use glimmer_global, only: dp,sizek

  implicit none

  !> derived type describing sigma coordinate system
  type vertCoord_type
     integer                           :: upn  =  0        !< the number of sigma levels
     real(dp),dimension(:),    GC_DYNARRAY_ATTRIB :: sigma !< the sigma levels
     real(dp),dimension(:),    GC_DYNARRAY_ATTRIB :: dupa  !< factors used for FD computations on non-regular grid
     real(dp),dimension(:),    GC_DYNARRAY_ATTRIB :: dupb  !< factors used for FD computations on non-regular grid
     real(dp),dimension(:),    GC_DYNARRAY_ATTRIB :: dupc  !< factors used for FD computations on non-regular grid
     real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: dups  !< factors used for FD computations on non-regular grid
     real(dp)                          :: dupn =  0.0      !< thickness of bottom sigma level
  end type vertCoord_type

  !> generic interface for initialising a vertical coordsystem
  interface initVertCoord
     module procedure initVertCoord_upn, initVertCoord_file, initVertCoord_array, initVertCoord_vertCoord
  end interface

  private :: initFDFactors
  
contains

  !> allocate memory for vertCoord_type
  subroutine vertCoord_allocate(params,upn)
    use glimmer_log, only : glimmer_allocErr
    implicit none
    type(vertCoord_type)             :: params !< the vertical coordinate type
    integer(kind=sizek), intent(in)  :: upn    !< the number of sigma levels

    ! local variables
    integer merr

    params%upn = upn
    GLIMMER_ALLOC1D(params%sigma,params%upn)
    GLIMMER_ALLOC1D(params%dupa,params%upn)
    GLIMMER_ALLOC1D(params%dupb,params%upn)
    GLIMMER_ALLOC1D(params%dupc,params%upn)
    GLIMMER_ALLOC2D(params%dups,params%upn,3)
  end subroutine vertCoord_allocate

  !> deallocate memory of vertCoord_type
  subroutine vertCoord_destroy(params)
    use glimmer_log, only : glimmer_deallocErr
    implicit none
    type(vertCoord_type) :: params
    
    ! local variables
    integer merr

    GLIMMER_DEALLOC(params%sigma)
    GLIMMER_DEALLOC(params%dupa)
    GLIMMER_DEALLOC(params%dupb)
    GLIMMER_DEALLOC(params%dupc)
    GLIMMER_DEALLOC(params%dups)
  end subroutine vertCoord_destroy

  !> print to verticalCoord_type to log
  subroutine vertCoord_print(params)
    use glimmer_log, only : write_log,msg_length
    implicit none
    
    type(vertCoord_type) :: params !< the vertical coordinate type    

    integer :: i
    character(len=msg_length) :: msg

    call write_log('Sigma levels:')
    call write_log('-------------')
    msg = ''
    do i=1,params%upn
       write(msg,'(a,x,f5.2)') trim(msg),params%sigma(i)
    end do
    call write_log(trim(msg))
    call write_log('')
  end subroutine vertCoord_print

  !> compute sigma levels and initialise vertCoord_type
  subroutine initVertCoord_upn(params,upn)
    implicit none

    type(vertCoord_type) :: params !< the vertical coordinate type
    integer, intent(in)  :: upn    !< the number of sigma coordinates

    integer :: up

    call vertCoord_allocate(params, upn)

    do up=1,upn
       params%sigma(up) = calc_sigma(real(up-1,kind=dp)/real(upn-1,kind=dp),2.d0)
    end do

    call initFDFactors(params)

  contains
    !> compute sigma level
    function calc_sigma(x,n)
      implicit none
      real(kind=dp), intent(in) :: x
      real(kind=dp), intent(in) :: n
      
      real(kind=dp) :: calc_sigma
      
      calc_sigma = (1-(x+1)**(-n))/(1-2**(-n))
    end function calc_sigma
  end subroutine initVertCoord_upn
    
  !> initialise vertical coordinate type from file
  subroutine initVertCoord_file(params,sigfile,upn)
    use glimmer_log, only : write_log, GM_FATAL
    implicit none
    type(vertCoord_type)  :: params         !< the vertical coordinate type
    character(len=*), intent(in) :: sigfile !< name of file containing sigma coordinates
    integer, intent(in) :: upn              !< the number of sigma coordinates

    integer :: up
    logical :: there
    integer ios
    integer, parameter :: unit = 999
    integer, parameter :: maxSig = 10000
    real(kind=dp), dimension(maxSig) :: tempSigma

    inquire (exist=there,file=sigfile)
    if (.not.there) then
       call write_log('Sigma levels file: '//trim(sigfile)//' does not exist',GM_FATAL)
    end if
    call write_log('Reading sigma file: '//sigfile)
    open(unit,file=sigfile)
    up = 1
    do while (up.le.maxSig)
       read(unit,*,err=10,iostat=ios) tempSigma(up)
       if (ios.ne.0) then
          exit
       end if
       up = up + 1
    end do
    close(unit)    
    
    if (up .ne. upn) then
       call write_log('Wrong number of sigma coordinates read from file',GM_FATAL)
    end if

    call initVertCoord_array(params,tempSigma(1:up))

    return
10  call write_log('something wrong with sigma coord file',GM_FATAL)
  end subroutine initVertCoord_file

  !> initialise vertical coordinate type from array of sigma levels
  subroutine initVertCoord_array(params,sigma)
    implicit none
    type(vertCoord_type)  :: params !< the vertical coordinate type
    real(dp),dimension(:) :: sigma  !< array containing vertical sigma levels

    integer :: up

    call vertCoord_allocate(params,size(sigma))

    params%sigma(:) = sigma(:)

    call initFDFactors(params)

  end subroutine initVertCoord_array

  !> initialise vertical coordinate type from another vertical coord type
  subroutine initVertCoord_vertCoord(params,vertCoord)
    implicit none
    type(vertCoord_type)  :: params    !< the vertical coordinate type
    type(vertCoord_type)  :: vertCoord !< fully populated vertical coordinate type to be copied

    call vertCoord_allocate(params,vertCoord%upn)

    params%sigma(:) = vertCoord%sigma(:)

    call initFDFactors(params)

  end subroutine initVertCoord_vertCoord

  !> compute finite difference factors
  subroutine initFDFactors(params)
    implicit none
    type(vertCoord_type) :: params !< the vertical coordinate type

    integer up

    params%dupa = (/ &
         0.0d0,      &
         0.0d0,      &
         ((params%sigma(up)   - params%sigma(up-1)) / &
         ((params%sigma(up-2) - params%sigma(up-1)) * (params%sigma(up-2) - params%sigma(up))), &
         up=3,params%upn)  &
         /)

    params%dupb = (/ &
         0.0d0,      &
         0.0d0,      &
         ((params%sigma(up)   - params%sigma(up-2)) / &
         ((params%sigma(up-1) - params%sigma(up-2)) * (params%sigma(up-1) - params%sigma(up))), &
         up=3,params%upn)  &
         /)

    params%dupc = (/ &
         (params%sigma(2) - params%sigma(1)) / 2.0d0, &
         ((params%sigma(up+1) - params%sigma(up-1)) / 2.0d0, up=2,params%upn-1), &
         (params%sigma(params%upn) - params%sigma(params%upn-1)) / 2.0d0  &
         /)

    params%dups = 0.0d0
    do up = 2, params%upn-1
       params%dups(up,1) = 1.d0/((params%sigma(up+1) - params%sigma(up-1)) * (params%sigma(up)   - params%sigma(up-1)))
       params%dups(up,2) = 1.d0/((params%sigma(up+1) - params%sigma(up-1)) * (params%sigma(up+1) - params%sigma(up)))
       params%dups(up,3) = 1.d0/( params%sigma(up+1) - params%sigma(up-1))
    end do

    params%dupn = params%sigma(params%upn) - params%sigma(params%upn-1)
  end subroutine initFDFactors

end module glimmer_vertcoord
