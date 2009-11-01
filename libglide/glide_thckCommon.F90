! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_thckCommon.f90 - part of the GLIMMER ice model     + 
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

module glide_thckCommon

  use glimmer_global, only: dp
  use physcon, only: rhoi, grav, gn
  use glimmer_paramets, only: vis0, thk0, vel0, len0

  implicit none

  integer, private, parameter :: p2 = gn-1
  integer, private, parameter :: p3 = 2*gn+1
  integer, private, parameter :: p4 = gn+2
  real(dp),private, parameter :: c = -2.0d0*vis0*(rhoi*grav)**gn*thk0**p3/(8.0d0*vel0*len0**gn)

contains

  subroutine velo_integrate_flwa(dups,depth,dintflwa,stagthck,flwa)
    
    !*FD this routine calculates the part of the vertically averaged velocity 
    !*FD field which solely depends on the temperature

    use glimmer_utils, only : hsum4, vertintg
    use glimmer_global, only: dp
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    real(dp),dimension(:),    intent(in)    :: dups
    real(dp),dimension(:),    intent(in)    :: depth
    real(dp),dimension(:,:),  intent(out)   :: dintflwa
    real(dp),dimension(:,:),  intent(in)    :: stagthck       !*FD ice thickness on staggered grid
    real(dp),dimension(:,:,:),intent(in)    :: flwa           !*FD ice flow factor

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------
    real(dp),dimension(size(flwa,1)) :: hrzflwa, intflwa 
    integer :: ew,ns,up,ewn,nsn,upn

    upn=size(flwa,1) ; ewn=size(flwa,2) ; nsn=size(flwa,3)

    do ns = 1,nsn-1
       do ew = 1,ewn-1
          if (stagthck(ew,ns) /= 0.0d0) then
             
             hrzflwa = hsum4(flwa(:,ew:ew+1,ns:ns+1))  
             intflwa(upn) = 0.0d0

             do up = upn-1, 1, -1
                intflwa(up) = intflwa(up+1) + depth(up) * (hrzflwa(up)+hrzflwa(up+1))
             end do

             dintflwa(ew,ns) = c * vertintg(dups,intflwa)

          else 

             dintflwa(ew,ns) = 0.0d0

          end if
       end do
    end do
  end subroutine velo_integrate_flwa

  !*****************************************************************************

  subroutine velo_calc_diffu(dintflwa,stagthck,dusrfdew,dusrfdns,diffu)

    !*FD calculate diffusivities
    use glimmer_global, only: dp

    implicit none
    
    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    real(dp),dimension(:,:),  intent(in)    :: dintflwa
    real(dp),dimension(:,:),  intent(in)    :: stagthck
    real(dp),dimension(:,:),  intent(in)    :: dusrfdew
    real(dp),dimension(:,:),  intent(in)    :: dusrfdns
    real(dp),dimension(:,:),  intent(out)   :: diffu

    where (stagthck .ne. 0.)
       diffu = dintflwa * stagthck**p4 * sqrt(dusrfdew**2 + dusrfdns**2)**p2 
    elsewhere
       diffu = 0.0d0
    end where

  end subroutine velo_calc_diffu

  !*****************************************************************************

  subroutine velo_calc_velo(dintflwa,depth,stagthck,dusrfdew,dusrfdns,flwa,diffu,ubas,vbas,uvel,vvel,uflx,vflx)

    !*FD calculate 3D horizontal velocity field and 2D flux field from diffusivity
    use glimmer_utils, only : hsum4
    use glimmer_global, only: dp
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    real(dp),dimension(:,:),  intent(in)    :: dintflwa
    real(dp),dimension(:),    intent(in)    :: depth
    real(dp),dimension(:,:),  intent(in)    :: stagthck
    real(dp),dimension(:,:),  intent(in)    :: dusrfdew
    real(dp),dimension(:,:),  intent(in)    :: dusrfdns
    real(dp),dimension(:,:,:),intent(in)    :: flwa
    real(dp),dimension(:,:),  intent(in)    :: diffu
    real(dp),dimension(:,:),  intent(in)    :: ubas
    real(dp),dimension(:,:),  intent(in)    :: vbas
    real(dp),dimension(:,:,:),intent(out)   :: uvel
    real(dp),dimension(:,:,:),intent(out)   :: vvel
    real(dp),dimension(:,:),  intent(out)   :: uflx
    real(dp),dimension(:,:),  intent(out)   :: vflx
    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------
    real(dp),dimension(size(flwa,1)) :: hrzflwa
    real(dp) :: factor
    real(dp),dimension(3)           :: const
    integer :: ew,ns,up,ewn,nsn,upn

    upn=size(flwa,1) ; ewn=size(stagthck,1) ; nsn=size(stagthck,2)
    
    do ns = 1,nsn
       do ew = 1,ewn
          if (stagthck(ew,ns) /= 0.0d0) then

             vflx(ew,ns) = diffu(ew,ns) * dusrfdns(ew,ns) + vbas(ew,ns) * stagthck(ew,ns)
             uflx(ew,ns) = diffu(ew,ns) * dusrfdew(ew,ns) + ubas(ew,ns) * stagthck(ew,ns)

             uvel(upn,ew,ns) = ubas(ew,ns)
             vvel(upn,ew,ns) = vbas(ew,ns)

             hrzflwa = hsum4(flwa(:,ew:ew+1,ns:ns+1))  

             factor = dintflwa(ew,ns)*stagthck(ew,ns)
             if (factor /= 0.0d0) then
                const(2) = c * diffu(ew,ns) / factor
                const(3) = const(2) * dusrfdns(ew,ns)  
                const(2) = const(2) * dusrfdew(ew,ns) 
             else
                const(2:3) = 0.0d0
             end if

             do up = upn-1, 1, -1
                const(1) = depth(up) * (hrzflwa(up)+hrzflwa(up+1))
                uvel(up,ew,ns) = uvel(up+1,ew,ns) + const(1) * const(2)
                vvel(up,ew,ns) = vvel(up+1,ew,ns) + const(1) * const(3) 
             end do

          else 

             uvel(:,ew,ns) = 0.0d0
             vvel(:,ew,ns) = 0.0d0
             uflx(ew,ns) = 0.0d0
             vflx(ew,ns) = 0.0d0 

          end if
       end do
    end do
  end subroutine velo_calc_velo

  subroutine glide_calclsrf(thck,topg,eus,lsrf)

    !*FD Calculates the elevation of the lower surface of the ice, 
    !*FD by considering whether it is floating or not.

    use glimmer_global, only : dp
    use physcon, only : rhoi, rhoo

    implicit none

    real(dp), intent(in),  dimension(:,:) :: thck !*FD Ice thickness
    real(dp), intent(in),  dimension(:,:) :: topg !*FD Bedrock topography elevation
    real, intent(in)                      :: eus  !*FD global sea level
    real(dp), intent(out), dimension(:,:) :: lsrf !*FD Lower ice surface elevation

    real(dp), parameter :: con = - rhoi / rhoo

    where (topg-eus < con * thck)
      lsrf = con * thck
    elsewhere
      lsrf = topg
    end where
  end subroutine glide_calclsrf


end module glide_thckCommon
