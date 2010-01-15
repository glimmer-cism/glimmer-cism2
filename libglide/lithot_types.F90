! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  lithot_types.F90 - part of the GLIMMER ice model         + 
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

!> module defining derived types for the bedrock temperature computations
!!
!! \author Magnus Hagdorn
!! \date January 2010

module lithot_types

  use glimmer_global, only : dp
  use glimmer_coordinates, only : horizCoord_type

  !> holds variables for temperature calculations in the lithosphere
  type lithot_type

     real(dp) :: geot   = -5.0d-2  !< parameter for constant geothermal heat flux W m^{-2}

     logical :: do_lithot = .False.    !< set to .True. if geothermal heat flux should be calculated

     real(dp),dimension(:,:,:),GC_DYNARRAY_ATTRIB :: temp    !< Three-dimensional temperature field.
     logical, dimension(:,:), GC_DYNARRAY_ATTRIB :: mask     !< whether the point has been ice covered at some time

     integer :: num_dim = 1                                  !< either 1 or 3 for 1D/3D calculations

     ! The sparse matrix and linearised arrays
     integer :: all_bar_top
     real(dp), dimension(:), GC_DYNARRAY_ATTRIB :: rhs
     real(dp), dimension(:), GC_DYNARRAY_ATTRIB :: answer
     real(dp), dimension(:), GC_DYNARRAY_ATTRIB :: supd,diag,subd

     real(dp), dimension(:), GC_DYNARRAY_ATTRIB :: deltaz     !< array holding grid spacing in z
     real(dp), dimension(:,:), GC_DYNARRAY_ATTRIB :: zfactors !< array holding factors for finite differences of vertical diffu
     real(dp) :: xfactor,yfactor                              !< factors for finite differences of horizontal diffu


     real :: surft = 2.         !< surface temperature, used for calculating initial temperature distribution
     real :: mart  = 2.         !< sea floor temperature 
     integer :: nlayer = 20     !< number of layers in lithosphere
     real :: rock_base = -5000. !< depth below sea-level at which geothermal heat gradient is applied
     
     integer :: numt = 0        !< number time steps for spinning up GTHF calculations

     real(dp) :: rho_r = 3300.0d0 !< The density of lithosphere (kg m$^{-3}$)
     real(dp) :: shc_r = 1000.0d0 !< specific heat capcity of lithosphere (J kg$^{-1}$ K$^{-1}$)
     real(dp) :: con_r = 3.3d0    !< thermal conductivity of lithosphere (W m$^{-1}$ K$^{-1}$)

     real(dp) :: diffu = 0. !< diffusion coefficient

     type(horizCoord_type):: hCoord      !< the horizontal coordinate system     

  end type lithot_type
  
contains
  
  !> allocate data for bedrock temperature computations
  subroutine lithot_allocate(litho,coords)
    use glimmer_coordinates, only : coordinates_type
    use glimmer_horizcoord, only : horizCoord_allocate
    implicit none
    type(lithot_type) :: litho           !< structure holding bedrock temperature configuration 
    type(coordinates_type), intent(in)  :: coords  !< the glide coordinate systems

    litho%hCoord = coords%ice_grid

    allocate(litho%temp(1:coords%ice_grid%size(1),1:coords%ice_grid%size(2),litho%nlayer)); litho%temp = 0.0
    call horizCoord_allocate(coords%ice_grid, litho%mask)

  end subroutine lithot_allocate

  !> deallocate data for bedrock temperature computations
  subroutine lithot_deallocate(litho)
    implicit none
    type(lithot_type) :: litho           !< structure holding bedrock temperature configuration 

    deallocate(litho%temp)
    deallocate(litho%mask)
  end subroutine lithot_deallocate

end module lithot_types
