! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_vertcoord.f90 - part of the Glimmer-CISM ice model + 
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

module glide_vertcoord

  use glimmer_global, only: dp

  implicit none

  type vertCoord
     integer                           :: upn  =  0
     real(dp),dimension(:),    pointer :: dupa => null()
     real(dp),dimension(:),    pointer :: dupb => null()
     real(dp),dimension(:),    pointer :: dupc => null()
     real(dp),dimension(:,:),  pointer :: dups => null()
     real(dp)                          :: dupn =  0.0
  end type vertCoord

contains

  subroutine initVertCoord(params,sigma)

    implicit none

    type(vertCoord),      intent(out) :: params
    real(dp),dimension(:),intent(in)  :: sigma

    integer :: up

    ! Find size of arrays
    params%upn = size(sigma)

    ! Deallocate arrays
    if (associated(params%dupa)) deallocate(params%dupa)
    if (associated(params%dupb)) deallocate(params%dupb)
    if (associated(params%dupc)) deallocate(params%dupc)
    if (associated(params%dups)) deallocate(params%dups)

    ! Allocate arrays
    allocate(params%dupa(params%upn))
    allocate(params%dupb(params%upn))
    allocate(params%dupc(params%upn))
    allocate(params%dups(params%upn,3))

    params%dupa = (/ &
         0.0d0,      &
         0.0d0,      &
         ((sigma(up)   - sigma(up-1)) / &
         ((sigma(up-2) - sigma(up-1)) * (sigma(up-2) - sigma(up))), &
         up=3,params%upn)  &
         /)

    params%dupb = (/ &
         0.0d0,      &
         0.0d0,      &
         ((sigma(up)   - sigma(up-2)) / &
         ((sigma(up-1) - sigma(up-2)) * (sigma(up-1) - sigma(up))), &
         up=3,params%upn)  &
         /)

    params%dupc = (/ &
         (sigma(2) - sigma(1)) / 2.0d0, &
         ((sigma(up+1) - sigma(up-1)) / 2.0d0, up=2,params%upn-1), &
         (sigma(params%upn) - sigma(params%upn-1)) / 2.0d0  &
         /)

    params%dups = 0.0d0
    do up = 2, params%upn-1
       params%dups(up,1) = 1.d0/((sigma(up+1) - sigma(up-1)) * (sigma(up)   - sigma(up-1)))
       params%dups(up,2) = 1.d0/((sigma(up+1) - sigma(up-1)) * (sigma(up+1) - sigma(up)))
       params%dups(up,3) = 1.d0/( sigma(up+1) - sigma(up-1))
    end do

    params%dupn = sigma(params%upn) - sigma(params%upn-1)

  end subroutine initVertCoord

end module glide_vertcoord
