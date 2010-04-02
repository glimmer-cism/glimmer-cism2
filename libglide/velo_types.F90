! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  velo_types.f90 - part of the GLIMMER ice model           + 
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
#include <glimmer_memory.inc>

!> module defining derived types for velocity work arrays
!!
!! \author Magnus Hagdorn
!! \date February 2010


module velo_types

  use glimmer_global, only : dp
  use glimmer_coordinates

  !> holds variables for velocity computations
  type velo_type
     type(horizCoord_type) :: velo_grid  !< coordinate system of the velocity grid
     type(horizCoord_type) :: ice_grid   !< coordinate system of the ice thickness grid
     type(vertCoord_type)   :: sigma_grid !< the sigma coordinate system
     
     real(dp),dimension(:),  GC_DYNARRAY_ATTRIB :: depth   
     real(dp),dimension(:),  GC_DYNARRAY_ATTRIB :: dups    
     real(dp),dimension(:),  GC_DYNARRAY_ATTRIB :: dupsw   
     real(dp),dimension(:),  GC_DYNARRAY_ATTRIB :: depthw  
     real(dp),dimension(:),  GC_DYNARRAY_ATTRIB :: suvel   
     real(dp),dimension(:),  GC_DYNARRAY_ATTRIB :: svvel   
     real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: dintflwa
     real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: dthckdtm !< Temporal derivative of thickness.
     real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: dusrfdtm !< Temporal derivative of upper surface elevation.     
     real(dp),dimension(4) :: c    = 0.0
     real(dp) :: watwd  = 3.0d0
     real(dp) :: watct  = 10.0d0
     real(dp) :: trc0   = 0.0
     real(dp) :: trcmin = 0.0d0
     real(dp) :: marine = 1.0d0
     real(dp) :: trcmax = 10.0d0
     real(dp) :: btrac_const = 0.0d0
     real(dp) :: btrac_slope = 0.0d0
     real(dp) :: btrac_max = 0.d0     
  end type velo_type
  
contains

  !> allocate data for velocity computations
  subroutine velo_allocate(velo,coords)
    use glimmer_log, only : glimmer_allocErr
    use glimmer_horizcoord, only : horizCoord_allocate
    use glimmer_vertcoord, only : initVertCoord
    implicit none
    type(velo_type) :: velo                      !< the derived type holding the velocity grid
    type(coordinates_type), intent(in) :: coords !< derived type holding the model coordinate systems

    ! local variables
    integer merr

    ! copy coordinate systems
    velo%velo_grid = coords%velo_grid
    velo%ice_grid = coords%ice_grid
    call initVertCoord(velo%sigma_grid,coords%sigma_grid)

    GLIMMER_ALLOC1D(velo%depth ,velo%sigma_grid%upn)
    GLIMMER_ALLOC1D(velo%dups ,velo%sigma_grid%upn)
    GLIMMER_ALLOC1D(velo%dupsw ,velo%sigma_grid%upn)
    GLIMMER_ALLOC1D(velo%depthw ,velo%sigma_grid%upn)
    GLIMMER_ALLOC1D(velo%suvel ,velo%sigma_grid%upn)
    GLIMMER_ALLOC1D(velo%svvel ,velo%sigma_grid%upn)

    call horizCoord_allocate(velo%velo_grid, velo%dintflwa)
    call horizCoord_allocate(velo%ice_grid, velo%dthckdtm)
    call horizCoord_allocate(velo%ice_grid, velo%dusrfdtm)
  end subroutine velo_allocate

  !> deallocate data for velocity computations
  subroutine velo_deallocate(velo)
    use glimmer_log, only : glimmer_deallocErr
    implicit none
    type(velo_type) :: velo                      !< the derived type holding the velocity grid

    ! local variables
    integer merr

    GLIMMER_DEALLOC(velo%depth)
    GLIMMER_DEALLOC(velo%dups)
    GLIMMER_DEALLOC(velo%dupsw)
    GLIMMER_DEALLOC(velo%depthw)
    GLIMMER_DEALLOC(velo%suvel)
    GLIMMER_DEALLOC(velo%svvel)
    GLIMMER_DEALLOC(velo%dintflwa)
    GLIMMER_DEALLOC(velo%dthckdtm)
    GLIMMER_DEALLOC(velo%dusrfdtm)

  end subroutine velo_deallocate

end module velo_types
