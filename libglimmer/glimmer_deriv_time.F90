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

module glimmer_deriv_time

  use glimmer_global, only: dp, sp

  implicit none

  !> work array for time derivatives
  type timeders_type
     real(sp) :: oldtime = 0.0
     real(dp),dimension(:,:,:),GC_DYNARRAY_ATTRIB :: olds
     integer  :: nwhich  = 2
  end type timeders_type

contains

  !> initialise the time derivative structure
  subroutine timeders_init(params,ewn,nsn)

    type(timeders_type),intent(inout) :: params !< Derived-type containing work data
    integer,               intent(in) :: ewn    !< the number of points along the x axis
    integer,               intent(in) :: nsn    !< the number of points along the y axis

    allocate(params%olds(ewn,nsn,params%nwhich))
    params%olds = 0.0d0

  end subroutine timeders_init

  !-----------------------------------------------------------------------------

  !> Calculates the time-derivative of a field. 
  subroutine timeders(params,ipvr,opvr,time,which)

    use glimmer_global, only : dp, sp
    use glimmer_paramets, only : conv

    implicit none 

    type(timeders_type) :: params                  !< Derived-type containing work data
    real(dp), intent(out), dimension(:,:) :: opvr  !< Input field
    real(dp), intent(in),  dimension(:,:) :: ipvr  !< Output (derivative) field
    real(sp), intent(in)                  :: time  !< current time
    integer,  intent(in)                  :: which !< selector for stored field

    real(sp) :: factor

    factor = (time - params%oldtime)
    if (factor .eq.0) then
       opvr = 0.0d0
    else
       factor = 1./factor
       opvr = conv * (ipvr - params%olds(:,:,which)) * factor
    end if

    params%olds(:,:,which) = ipvr

    if (which == params%nwhich) then
      params%oldtime = time
    end if

  end subroutine timeders

  !-----------------------------------------------------------------------------

  !> cleanup time-derivative work data
  subroutine timeders_final(params)

     type(timeders_type) :: params  !< Derived-type containing work data

     deallocate(params%olds)

  end subroutine timeders_final

end module glimmer_deriv_time
