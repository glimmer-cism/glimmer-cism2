! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_coordinates.f90 - part of the GLIMMER ice model  + 
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

!> module defining glimmer collection of coordinate systems
!!
!! \author Magnus Hagdorn
!! \date November 2009
module glimmer_coordinates
  
  use glimmer_vertcoord, only  : vertCoord_type
  use glimmer_horizcoord, only : horizCoord_type
  use glimmer_map_types, only : glimmap_proj

  !> define the various coordinate systems used by glide
  type coordinates_type
     type(horizCoord_type) :: ice_grid   !< coordinate system of the ice grid
     type(horizCoord_type) :: velo_grid  !< coordinate system of the velocity grid
     type(vertCoord_type)   :: sigma_grid !< the sigma coordinate system
     type(glimmap_proj) :: projection     !< the geographic projection
  end type coordinates_type
end module glimmer_coordinates
