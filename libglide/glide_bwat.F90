! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_bwat.f90 - part of the GLIMMER ice model           + 
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

module glide_bwat

  implicit none

contains

    !-------------------------------------------------------------------

  subroutine calcbwat(which,bmlt,bwat,thck,topg,btem, &
       floater,thklim,dt,ewn,nsn,dew,dns,periodic_ew, &
       bwat_smooth,hydtim)

    use glimmer_global,   only: dp 
    use glimmer_paramets, only: thk0, tim0, len0
    use glimmer_utils,    only: stagvarb
    use physcon,          only: rhoi, rhow, grav, scyr
    use glimmer_pmpt,     only: calcpmptb

    implicit none

    integer,                  intent(in)    :: which
    real(dp), dimension(:,:), intent(in)    :: bmlt
    real(dp), dimension(:,:), intent(inout) :: bwat
    real(dp), dimension(:,:), intent(in)    :: thck
    real(dp), dimension(:,:), intent(in)    :: topg
    real(dp), dimension(:,:), intent(in)    :: btem
    logical,  dimension(:,:), intent(in)    :: floater
    real(dp),                 intent(in)    :: thklim
    real(dp),                 intent(in)    :: dt
    integer,                  intent(in)    :: ewn
    integer,                  intent(in)    :: nsn
    real(dp),                 intent(in)    :: dew
    real(dp),                 intent(in)    :: dns
    integer,                  intent(in)    :: periodic_ew
    real(dp),                 intent(in)    :: bwat_smooth
    real(dp),                 intent(in)    :: hydtim

    real(dp), dimension(2), parameter :: &
         blim = (/ 0.00001 / thk0, 0.001 / thk0 /)

    real(dp) :: dwphidew, dwphidns, dwphi, pmpt, bave

    integer :: t_wat,ns,ew

    real(dp),dimension(ewn,nsn) :: smth
    real(dp),dimension(ewn,nsn) :: wphi
    real(dp),dimension(ewn,nsn) :: bwatu
    real(dp),dimension(ewn,nsn) :: bwatv
    real(dp),dimension(ewn,nsn) :: fluxew
    real(dp),dimension(ewn,nsn) :: fluxns
    real(dp),dimension(ewn,nsn) :: bint

    integer  :: nwat
    real(dp) :: watvel
    real(dp) :: dt_wat
    real(dp) :: hydtim_local
    real(dp) :: estimate
    real(dp),dimension(8) :: tempwk_c

    ! Initialisation

    select case(which)
    case(0)

       hydtim_local = tim0  / (hydtim * scyr)
       estimate     = 0.2d0 /  hydtim_local

       call find_dt_wat(dt,estimate,dt_wat,nwat) 

       tempwk_c = (/ &
            dt_wat, &
            1.0d0 - 0.5d0 * dt_wat * hydtim_local, &
            1.0d0 + 0.5d0 * dt_wat * hydtim_local, &
            0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /) 

    case(1)

       watvel = hydtim_local * tim0 / (scyr * len0)
       estimate = (0.2d0 * watvel) / min(dew,dns)

       call find_dt_wat(dt,estimate,dt_wat,nwat) 

       tempwk_c = (/ &
            rhow * grav, &
            rhoi * grav, &
            2.0d0  * dew, &
            2.0d0  * dns, &
            0.25d0 * dt_wat / dew, &
            0.25d0 * dt_wat / dns, &
            0.5d0  * dt_wat / dew, &
            0.5d0  * dt_wat / dns /)
    end select


    select case (which)
    case(0)

       do t_wat = 1,nwat
          do ns = 1,nsn
             do ew = 1,ewn

                if (thklim < thck(ew,ns) .and. .not. floater(ew,ns)) then
                   bwat(ew,ns) = (tempwk_c(1) * bmlt(ew,ns) + tempwk_c(2) * bwat(ew,ns)) / &
                        tempwk_c(3)
                   if (blim(1) > bwat(ew,ns)) then
                      bwat(ew,ns) = 0.0d0
                   end if
                else
                   bwat(ew,ns) = 0.0d0
                end if

             end do
          end do
       end do

       smth = 0.
       do ns = 2,nsn-1
          do ew = 2,ewn-1
             call smooth_bwat(ew-1,ew,ew+1,ns-1,ns,ns+1)
          end do
       end do
       ! apply periodic BC
       if (periodic_ew.eq.1) then
          do ns = 2,nsn-1
             call smooth_bwat(ewn-1,1,2,ns-1,ns,ns+1)
             call smooth_bwat(ewn-1,ewn,2,ns-1,ns,ns+1)
          end do
       end if

       bwat(1:ewn,1:nsn) = smth(1:ewn,1:nsn)

    case(1)
       ! apply periodic BC
       if (periodic_ew.eq.1) then
          write(*,*) 'Warning, periodic BC are not implement for this case yet'
       end if
       ! ** add any melt_water

       bwat = max(0.0d0,bwat + dt * bmlt)

       wphi = 0.
       bwatu = 0.
       bwatv = 0.
       fluxew = 0.
       fluxns = 0.
       bint = 0.


       ! ** split time evolution into steps to avoid CFL problems

       do t_wat = 1,nwat

          ! ** find potential surface using paterson p112, eq 4
          ! ** if no ice then set to sea level or land surface potential
          ! ** if frozen then set high 

          do ns = 1,nsn
             do ew = 1,ewn
                if (thklim < thck(ew,ns) .and. .not. floater(ew,ns)) then
                   call calcpmptb(pmpt,thck(ew,ns))
                   if (btem(ew,ns) == pmpt) then
                      wphi(ew,ns) = tempwk_c(1) * (topg(ew,ns) + bwat(ew,ns)) + tempwk_c(2) * thck(ew,ns)
                   else
                      wphi(ew,ns) = tempwk_c(1) * (topg(ew,ns) + thck(ew,ns))
                   end if
                else 
                   wphi(ew,ns) = max(tempwk_c(1) * topg(ew,ns),0.0d0)
                end if
             end do
          end do

          ! ** determine x,y components of water velocity assuming
          ! ** contstant velocity magnitude and using potential
          ! ** to determine direction

          do ns = 2,nsn-1
             do ew = 2,ewn-1
                if (thck(ew,ns) > thklim) then

                   dwphidew = (wphi(ew+1,ns) - wphi(ew-1,ns)) / tempwk_c(3)       
                   dwphidns = (wphi(ew,ns+1) - wphi(ew,ns-1)) / tempwk_c(4)  

                   dwphi = - watvel / sqrt(dwphidew**2 + dwphidns**2)

                   bwatu(ew,ns) = dwphi * dwphidew  
                   bwatv(ew,ns) = dwphi * dwphidns  

                else
                   bwatu(ew,ns) = 0.0d0
                   bwatv(ew,ns) = 0.0d0
                end if
             end do
          end do

          ! ** use two-step law wendroff to solve dW/dt = -dF/dx - dF/dy

          ! ** 1. find fluxes F=uW

          fluxew = bwat * bwatu
          fluxns = bwat * bwatv

          ! ** 2. do 1st LW step on staggered grid for dt/2

          do ns = 1,nsn-1
             do ew = 1,ewn-1

                bave = 0.25 * sum(bwat(ew:ew+1,ns:ns+1))

                if (bave > 0.0d0) then

                   bint(ew,ns) = bave - &
                        tempwk_c(5) * (sum(fluxew(ew+1,ns:ns+1)) - sum(fluxew(ew,ns:ns+1))) - &
                        tempwk_c(6) * (sum(fluxns(ew:ew+1,ns+1)) - sum(fluxns(ew:ew+1,ns)))

                else
                   bint(ew,ns) = 0.0d0
                end if
             end do
          end do

          ! ** 3. find fluxes F=uW on staggered grid griven new Ws

          fluxew(1:ewn-1,1:nsn-1) = bint * 0.25 * &
               (bwatu(1:ewn-1,1:nsn-1) + &
               bwatu(2:ewn,1:nsn-1) + &
               bwatu(1:ewn-1,2:nsn) + &
               bwatu(2:ewn,2:nsn))
          fluxns(1:ewn-1,1:nsn-1) = bint * 0.25 * &
               (bwatv(1:ewn-1,1:nsn-1) + &
               bwatv(2:ewn,1:nsn-1) + &
               bwatv(1:ewn-1,2:nsn) + &
               bwatv(2:ewn,2:nsn))

          ! ** 4. finally do 2nd LW step to get back on to main grid

          do ns = 2,nsn-1
             do ew = 2,ewn-1
                if (bwat(ew,ns) > 0.0d0) then

                   bwat(ew,ns) = bwat(ew,ns) - &
                        tempwk_c(7) * (sum(fluxew(ew,ns-1:ns)) - sum(fluxew(ew-1,ns-1:ns))) - &
                        tempwk_c(8) * (sum(fluxns(ew-1:ew,ns)) - sum(fluxns(ew-1:ew,ns-1)))

                else
                   bwat(ew,ns) = 0.0d0
                end if
             end do
          end do
       end do

       where (blim(1) > bwat) 
          bwat = 0.0d0
       end where

    case default

       bwat = 0.0d0

    end select

  contains

    subroutine smooth_bwat(ewm,ew,ewp,nsm,ns,nsp)
      ! smoothing basal water distrib
      implicit none
      integer, intent(in) :: ewm,ew,ewp,nsm,ns,nsp
      if (blim(2) < bwat(ew,ns)) then
         smth(ew,ns) = bwat(ew,ns) + bwat_smooth * &
              (bwat(ewm,ns) + bwat(ewp,ns) + bwat(ew,nsm) + bwat(ew,nsp) - 4.0d0 * bwat(ew,ns))
      else 
         smth(ew,ns) = bwat(ew,ns)
      end if   
    end subroutine smooth_bwat

  end subroutine calcbwat
  
  !-------------------------------------------------------------------

  subroutine find_dt_wat(dttem,estimate,dt_wat,nwat)
    
    use glimmer_global, only: dp

    implicit none
    
    real(dp), intent(out) :: dt_wat
    integer, intent(out) :: nwat
    real(dp), intent(in) :: dttem, estimate
    
    nwat = int(dttem/estimate) + 1
    dt_wat = dttem / nwat

  end subroutine find_dt_wat


end module glide_bwat
