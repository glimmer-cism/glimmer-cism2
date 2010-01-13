! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  isostasy.f90 - part of the GLIMMER ice model             + 
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

!> calculate isostatic adjustment due to changing surface loads
!!
!! \author Magnus Hagdorn
!! \date 2006
module isostasy
  
  use glimmer_global, only : dp,sp
  use isostasy_setup
  use isostasy_types
  use isostasy_el

  private :: relaxing_mantle
  
contains
  !> initialise isostasy calculations
  subroutine init_isostasy(isos,tstart)
    use physcon,  only: scyr
    use glimmer_paramets, only: tim0
    implicit none
    type(isos_type) :: isos !< structure holding isostasy configuration
    real(sp), intent(in) :: tstart !< time when first isostasy computation should be done

    if (isos%lithosphere .eq. 1) then
       call init_elastic(isos%rbel,isos%hCoord%delta(1))
    end if
    isos%next_calc = tstart

    ! scale tau
    isos%relaxed_tau = isos%relaxed_tau * scyr / tim0

  end subroutine init_isostasy
  
  !> calculate surface load factors due to water and ice distribution
  subroutine isos_icewaterload(isos,topg,thck,eus)
    use physcon
    implicit none
    type(isos_type) :: isos !< structure holding isostasy configuration
    real(dp), dimension(:,:), intent(in) :: topg !< bedrock topography
    real(dp), dimension(:,:), intent(in) :: thck !< ice thicknesses
    real(sp), intent(in) :: eus                  !< the current eustatic sealevel

    real(kind=dp) :: ice_mass, water_depth, water_mass
    integer :: ew,ns
  
     do ns=1,isos%hCoord%size(1)
       do ew=1,isos%hCoord%size(2)
          ice_mass = rhoi * thck(ew,ns)
          if (topg(ew,ns)-eus.lt.0) then             ! check if we are below sea level
             water_depth = eus - topg(ew,ns)
             water_mass = rhoo * water_depth
             ! Just the water load due to changes in sea-level
             isos%load_factors(ew,ns) = rhoo* eus/rhom
             ! Check if ice is not floating
             if ( ice_mass .gt. water_mass ) then
                isos%load_factors(ew,ns) = isos%load_factors(ew,ns) + (ice_mass - water_mass)/rhom
             end if
          else                                       ! bedrock is above sea level
             isos%load_factors(ew,ns) = ice_mass/rhom
          end if
       end do
    end do
    isos%new_load = .true.
  end subroutine isos_icewaterload

  !> calculate isostatic adjustment due to changing surface loads
  subroutine isos_isostasy(isos,topg,dt)
    implicit none
    type(isos_type) :: isos !< structure holding isostasy configuration
    real(dp), dimension(:,:), intent(inout) :: topg !< bedrock topography
    real(dp), intent(in) :: dt !< current time step

    ! update load if necessary
    if (isos%new_load) then
       call isos_lithosphere(isos)
       ! update bed rock with (non-viscous) fluid mantle
       if (isos%asthenosphere .eq. 0) then
          topg = isos%relx - isos%load
       end if
       isos%new_load = .false.
    end if
    ! update bed rock with relaxing mantle
    if (isos%asthenosphere .eq. 1) then
       call relaxing_mantle(isos,topg,dt)
    end if
  end subroutine isos_isostasy

  !> compute isostatic load
  subroutine isos_lithosphere(isos)
    implicit none
    type(isos_type) :: isos !< structure holding isostasy configuration

    if (isos%lithosphere .eq. 0) then
       ! local lithosphere
       isos%load = isos%load_factors
    else if (isos%lithosphere .eq. 1) then
       call calc_elastic(isos%rbel,isos%load,isos%load_factors)
    end if
  end subroutine isos_lithosphere

  !> Calculate the relaxed topography, assuming the isostatic depression
  !! is the equilibrium state for the current topography.
  subroutine isos_relaxed(isos,topg,thck,eus)
    implicit none
    type(isos_type) :: isos !< structure holding isostasy configuration    
    real(dp), dimension(:,:), intent(in) :: topg !< bedrock topography
    real(dp), dimension(:,:), intent(in) :: thck !< ice thicknesses
    real(sp), intent(in) :: eus                  !< the current eustatic sealevel

    ! Calculate the load
    call isos_icewaterload(isos,topg,thck,eus)
    ! Apply lithosphere model
    call isos_lithosphere(isos)
    ! Add to present topography to get relaxed topography
    isos%relx = topg + isos%load

  end subroutine isos_relaxed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> approximate mantle with a relaxing half-space: dh/dt=-1/tau*(w-h)
  subroutine relaxing_mantle(isos,topg,dt)
    implicit none
    type(isos_type) :: isos !< structure holding isostasy configuration
    real(dp), dimension(:,:), intent(inout) :: topg !< bedrock topography
    real(dp), intent(in) :: dt !< current time step
    
    integer :: ew,ns
    real(kind=dp) :: ft1, ft2

    ft1 = exp(-dt/isos%relaxed_tau)
    ft2 = 1. - ft1
    
    do ns=1,isos%hCoord%size(1)
       do ew=1,isos%hCoord%size(2)
          topg(ew,ns) = ft2*(isos%relx(ew,ns)-isos%load(ew,ns)) + ft1*topg(ew,ns)
       end do
    end do
  end subroutine relaxing_mantle

end module isostasy
