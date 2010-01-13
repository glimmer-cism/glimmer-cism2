! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  isostasy_types.f90 - part of the GLIMMER ice model       + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
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
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> module defining derived types used for isostasy calculations
!!
!! \author Magnus Hagdorn
!! \date 2006

module isostasy_types

  use glimmer_global, only : dp
  use glimmer_coordinates, only : horizCoord_type

  !> Holds data used by isostatic adjustment calculations for an elastic lithosphere
  type isostasy_elastic
     real(kind=dp) :: d = 0.24e25                           !< flexural rigidity
     real(kind=dp) :: lr                                    !< radius of relative stiffness
     real(kind=dp) :: a                                     !< radius of disk
     real(kind=dp) :: c1,c2,cd3,cd4                         !< coefficients
     real(kind=dp), dimension(:,:), GC_DYNARRAY_ATTRIB :: w !< matrix operator for lithosphere deformation
     integer :: wsize                                       !< size of operator (0:rbel_wsize, 0:rbel_wsize), operator is axis symmetric
  end type isostasy_elastic

  !> contains isostasy configuration
  type isos_type
     logical :: do_isos = .False.    !< set to .True. if isostatic adjustment should be handled
     !> method for calculating equilibrium bedrock depression:
     !!  - <tt>0</tt> local lithosphere, equilibrium bedrock depression is found using Archimedes' principle
     !!  - <tt>1</tt> elastic lithosphere, flexural rigidity is taken into account
     integer :: lithosphere = 0      
     !> method for approximating the mantle
     !!  - <tt>0</tt> fluid mantle, isostatic adjustment happens instantaneously
     !!  - <tt>1</tt> relaxing mantle, mantle is approximated by a half-space
     integer :: asthenosphere = 0
     real :: relaxed_tau = 4000. !< characteristic time constant of relaxing mantle
     real :: period = 500. !< lithosphere update period
     real :: next_calc !< when to upate lithosphere
     logical :: new_load=.false. !< set to true if there is a new surface load
     type(isostasy_elastic) :: rbel !< structure holding elastic lithosphere setup

     real(dp),dimension(:,:), GC_DYNARRAY_ATTRIB :: relx !< The elevation of the relaxed topography, by <tt>thck0</tt>.
     real(dp),dimension(:,:), GC_DYNARRAY_ATTRIB :: load !< the load imposed on lithosphere
     real(dp),dimension(:,:), GC_DYNARRAY_ATTRIB :: load_factors !< temporary used for load calculation

     type(horizCoord_type):: hCoord      !< the horizontal coordinate system
  end type isos_type  

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_ISOSTASY_TYPES
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_ISOSTASY_TYPES
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_ISOSTASY_TYPES
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_ISOSTASY_TYPES
!MH!#endif

  !> allocate data for isostasy calculations
   subroutine isos_allocate(isos, coords)
    use glimmer_coordinates, only : coordinates_type
    use glimmer_horizcoord, only : horizCoord_allocate
    implicit none
    type(isos_type) :: isos           !< structure holding isostasy configuration
    type(coordinates_type), intent(in)  :: coords  !< the glide coordinate systems

    isos%hCoord = coords%ice_grid

    call horizCoord_allocate(isos%hCoord,isos%relx)
    call horizCoord_allocate(isos%hCoord,isos%load)
    call horizCoord_allocate(isos%hCoord,isos%load_factors)
  end subroutine isos_allocate

  !> deallocate data for isostasy calculations
  subroutine isos_deallocate(isos)
    implicit none
    type(isos_type) :: isos                !< structure holding isostasy configuration
    
    deallocate(isos%relx)
    deallocate(isos%load)
    deallocate(isos%load_factors)
  end subroutine isos_deallocate
end module isostasy_types
