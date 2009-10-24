! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_temp.f90 - part of the GLIMMER ice model         + 
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

module glide_temp

contains

  subroutine calcTemp_asSurfTemp(temp,artm)

    !*FD Sets ice column temperature to surface temperature

    use glimmer_log, only: write_log, GM_FATAL
    use glimmer_global, only: dp,sp

    implicit none
   
    real(dp),dimension(:,0:,0:),intent(out) :: temp !< Ice temperature (upn,0:ewn+1,0:nsn+1)
    real(sp),dimension(:,:),    intent(in)  :: artm !< Surface air temperature (ewn,nsn)

    integer :: ns,ew,nsn,ewn

    ewn = size(artm,1)
    nsn = size(artm,2)

#ifdef DEBUG
    ! Check array sizes match
    if (ewn/=ubound(temp,2)-1.or.nsn/=ubound(temp,3)-1) then
       call write_log('Array shape mismatch in calcTemp_asSurfTemp', &
            GM_FATAL,__FILE__,__LINE__)
    end if
#endif

    do ns = 1,nsn
       do ew = 1,ewn
          temp(:,ew,ns) = dmin1(0.0d0,dble(artm(ew,ns)))
       end do
    end do

  end subroutine calcTemp_asSurfTemp

  !------------------------------------------------------------------------------------

  subroutine calcTemp_VerticalProfile(temp,artm,sigma,thck,bwat)

    !*FD Sets ice column temperature to surface temperature at top
    !*FD and zero at bottom, correcting for pressure melting point

    use glimmer_log, only: write_log, GM_FATAL
    use glimmer_global, only: dp,sp
    use glimmer_pmpt, only: corrpmpt

    implicit none

    real(dp),dimension(:,0:,0:),intent(out) :: temp  !< Ice temperature (upn,0:ewn+1,0:nsn+1)
    real(sp),dimension(:,:),    intent(in)  :: artm  !< Surface air temperature (ewn,nsn)
    real(dp),dimension(:),      intent(in)  :: sigma !< Sigma levels (upn)
    real(dp),dimension(:,:),    intent(in)  :: thck  !< Ice thickness (ewn,nsn)
    real(dp),dimension(:,:),    intent(in)  :: bwat  !< Basal water depth (ewn,nsn)

    integer :: ns,ew,nsn,ewn

    ewn = size(artm,1)
    nsn = size(artm,2)

#ifdef DEBUG
    ! Check array sizes match
    if ( ewn/=ubound(temp,2)-1.or. &
         nsn/=ubound(temp,3)-1.or. &
         size(sigma)/=size(temp,1).or. &
         shape(artm)/=shape(thck).or. &
         shape(artm)/=shape(bwat)) then
       call write_log('Array shape mismatch in calcTemp_VerticalProfile', &
            GM_FATAL,__FILE__,__LINE__)
    end if
#endif

    do ns = 1,nsn
       do ew = 1,ewn
          temp(:,ew,ns) = dmin1(0.0d0,dble(artm(ew,ns))) * (1.0d0 - sigma)
          call corrpmpt(temp(:,ew,ns),thck(ew,ns),bwat(ew,ns),sigma)
       end do
    end do

  end subroutine calcTemp_VerticalProfile

end module glide_temp

