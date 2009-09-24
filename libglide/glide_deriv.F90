! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_deriv.f90 - part of the Glimmer-CISM ice model     + 
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

module glide_deriv

  implicit none

contains

  subroutine timeders(thckwk,ipvr,opvr,mask,time,which)

    !*FD Calculates the time-derivative of a field. This subroutine is used by 
    !*FD the temperature solver only.

    use glimmer_global, only : dp, sp
    use paramets, only : conv
    use glide_types, only: glide_thckwk

    implicit none 

    type(glide_thckwk) :: thckwk    !*FD Derived-type containing work data
    real(dp), intent(out), dimension(:,:) :: opvr  !*FD Input field
    real(dp), intent(in),  dimension(:,:) :: ipvr  !*FD Output (derivative) field
    real(sp), intent(in)                  :: time  !*FD current time
    integer,  intent(in),  dimension(:,:) :: mask  !*FD mask for calculation
    integer,  intent(in)                  :: which !*FD selector for stored field

    real(sp) :: factor

    factor = (time - thckwk%oldtime)
    if (factor .eq.0) then
       opvr = 0.0d0
    else
       factor = 1./factor
       where (mask /= 0)
          opvr = conv * (ipvr - thckwk%olds(:,:,which)) * factor
       elsewhere
          opvr = 0.0d0
       end where
    end if

    thckwk%olds(:,:,which) = ipvr

    if (which == thckwk%nwhich) then
      thckwk%oldtime = time
    end if

  end subroutine timeders

end module glide_deriv
