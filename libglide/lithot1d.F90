! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_lithot1d.f90 - part of the GLIMMER ice model       + 
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

!> module for 1D temperature calculations in the upper lithosphere
!!
!! \author Magnus Hagdorn
!! \date 2006

module lithot1d

  use lithot_types, only : lithot_type

contains
  !> initialise 1D geothermal heat flux
  subroutine init_lithot1d(litho)
    implicit none
    type(lithot_type) :: litho            !< structure holding bedrock temperature configuration    

    ! allocate memory for 1D code
    allocate(litho%rhs(litho%nlayer))
    allocate(litho%subd(litho%nlayer))
    allocate(litho%diag(litho%nlayer))
    allocate(litho%supd(litho%nlayer))
    
    ! setup coefficient matrix
    litho%subd(:) =    - litho%zfactors(1,:)
    litho%diag(:) = 1. + litho%zfactors(2,:)
    litho%supd(:) =    - litho%zfactors(3,:)
    ! and the boundary conditions
    ! top face
    ! simply match air temperature where no ice and basal temperature where ice
    litho%subd(1) = 0.
    litho%diag(1) = 1.
    litho%supd(1) = 0.
    ! bottom face
    ! keep constant
    litho%subd(litho%nlayer) = 0.
    litho%diag(litho%nlayer) = 1.
    litho%supd(litho%nlayer) = 0.
  end subroutine init_lithot1d

  !> compute 1D geothermal heat flux
  subroutine calc_lithot1d(litho,thkmask,base_temp,air_temp)
    use glimmer_utils
    use glimmer_mask
    implicit none
    type(lithot_type) :: litho            !< structure holding bedrock temperature configuration    
    integer, dimension(:,:), intent(in) :: thkmask !< surface type mask
    real(dp), dimension(0:,0:), intent(in) :: base_temp !< temperature at ice base
    real(sp), dimension(:,:), intent(in) :: air_temp  !< air temperature

    integer i,j,k

    ! loop over grid
    do j=1,litho%hCoord%size(2)
       do i=1,litho%hCoord%size(1)
          ! calculate RHS for upper BC
          if (is_ground(thkmask(i,j)) .and. .not. is_thin(thkmask(i,j)) ) then
             litho%rhs(1) = base_temp(i,j) ! ice basal temperature
             litho%mask(i,j) = .true.
          else
             if (litho%mask(i,j)) then
                if (is_ocean(thkmask(i,j))) then
                   litho%rhs(1) = litho%mart
                else if (is_land(thkmask(i,j))) then
                   litho%rhs(1) = air_temp(i,j) ! air temperature outside ice sheet
                end if
             end if
          end if

          if (litho%mask(i,j)) then
             ! calculate RHS for rest
             do k=2,litho%nlayer-1
                litho%rhs(k) = - litho%subd(k)*litho%temp(i,j,k-1) &
                     + (2.-litho%diag(k))*litho%temp(i,j,k) &
                     - litho%supd(k)*litho%temp(i,j,k+1)
             end do
             litho%rhs(litho%nlayer) = litho%temp(i,j,litho%nlayer)

             ! solve tri-diagonal matrix eqn
             call tridiag(litho%subd(1:), &
                  litho%diag(:), &
                  litho%supd(:litho%nlayer), &
                  litho%temp(i,j,:) ,                 &
                  litho%rhs(:))
          end if
       end do
    end do
  end subroutine calc_lithot1d

  !> deallocate memory for 1D geothermal heat flux computations
  subroutine finalise_lithot1d(litho)
    implicit none
    type(lithot_type) :: litho            !< structure holding bedrock temperature configuration    

    deallocate(litho%rhs)
    deallocate(litho%subd)
    deallocate(litho%diag)
    deallocate(litho%supd)
  end subroutine finalise_lithot1d

end module lithot1d
