! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_thckADI.f90 - part of the GLIMMER ice model        + 
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

module glide_thckADI

  implicit none

contains

  subroutine thckADI_init(params,ewn,nsn,dew,dns,thklim)

    use glide_types, only: thckADI_type
    use glimmer_global, only: dp

    implicit none

    type(thckADI_type),intent(out) :: params
    integer,           intent(in)  :: ewn
    integer,           intent(in)  :: nsn
    real(dp),          intent(in)  :: dew
    real(dp),          intent(in)  :: dns
    real(dp),          intent(in)  :: thklim

    params%ewn = ewn
    params%nsn = nsn
    params%dew = dew
    params%dns = dns
    params%thklim = thklim

  end subroutine thckADI_init

  !--------------------------------------------------------------------

  subroutine stagleapthck(params,model,newtemps,thck,acab,lsrf,topg,btrc,ubas,vbas,dt,eus)
    
    !*FD this subroutine solves the ice sheet thickness equation using the ADI scheme
    !*FD diffusivities are updated for each half time step

    use glide_setup, only: glide_calclsrf
    use glide_velo, only: velo_calc_velo, velo_integrate_flwa, velo_calc_diffu
    use glimmer_utils
    use glide_types, only: glide_global_type, thckADI_type
    use glimmer_global, only: dp, sp
    use physcon, only: rhoi, grav

    implicit none

    ! subroutine arguments
    type(thckADI_type),intent(inout) :: params
    type(glide_global_type) :: model
    logical,                                      intent(in)    :: newtemps !< true when we should recalculate Glen's A
    real(dp),dimension(params%ewn,params%nsn),    intent(inout) :: thck     !< Ice thickness
    real(sp),dimension(params%ewn,params%nsn),    intent(in)    :: acab     !< Mass balance
    real(dp),dimension(params%ewn,params%nsn),    intent(inout) :: lsrf     !< Lower surface elevation
    real(dp),dimension(params%ewn,params%nsn),    intent(in)    :: topg     !< basal topography
    real(dp),dimension(params%ewn-1,params%nsn-1),intent(in)    :: btrc     !< basal traction coefficient 
    real(dp),dimension(params%ewn-1,params%nsn-1),intent(out)   :: ubas     !< basal x velocity
    real(dp),dimension(params%ewn-1,params%nsn-1),intent(out)   :: vbas     !< basal y velocity
    real(dp),                                     intent(in)    :: dt       !< Timestep
    real(sp),                                     intent(in)    :: eus      !< Sea level

    ! local variables
    integer ew,ns, n
    real(dp), parameter :: rhograv = - rhoi * grav

    real(dp), dimension(max(params%ewn,params%nsn)) :: alpha
    real(dp), dimension(max(params%ewn,params%nsn)) :: beta
    real(dp), dimension(max(params%ewn,params%nsn)) :: gamma
    real(dp), dimension(max(params%ewn,params%nsn)) :: delta
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: fslip
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: stagthck

    if (all(thck==0.0)) then

       thck = dmax1(0.0d0,thck + acab * dt)
#ifdef DEBUG       
       print *, "* thck empty - net accumulation added"
#endif
    else

       ! Calculate staggered thickness
       call stagvarb(thck,stagthck,params%ewn,params%nsn)

       ! First part of basal velocity calculation
       fslip =  rhograv * btrc

       if (newtemps) then
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(     &
               model%velowk,            &
               stagthck,                &
               model%temper%flwa)
       end if

       ! Second part of velocity calculation
       where (params%thklim < stagthck)
          ubas = fslip * stagthck**2  
       elsewhere
          ubas = 0.0d0
       end where

       ! calculate diffusivity
       call velo_calc_diffu(         &
            model%velowk,            &
            stagthck, &
            model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns, &
            model%velocity%diffu)

       model%velocity%total_diffu(:,:) = model%velocity%diffu(:,:) + ubas(:,:)

       ! first ADI step, solve thickness equation along rows j
       n = params%ewn
       do ns=2,params%nsn-1
          call adi_tri ( alpha,                 &
                         beta,                  &
                         gamma,                 &
                         delta,                 &
                         thck(:,ns),            &
                         lsrf(:,ns),            &
                         acab(:,ns)-real(model%options%basal_mbal)*real(model%temper%bmlt(:,ns),sp),           &
                         model%velocity%vflx(:,ns),          &
                         model%velocity%vflx(:,ns-1),        &
                         model%velocity%total_diffu(:,ns),   &
                         model%velocity%total_diffu(:,ns-1), &
                         dt,                                 &
                         params%dew,                         &
                         params%dns )

          call tridiag(alpha(1:n),    &
                       beta(1:n),     &
                       gamma(1:n),    &
                       model%thckwk%oldthck(:,ns), &
                       delta(1:n))
       end do

       model%thckwk%oldthck(:,:) = max(model%thckwk%oldthck(:,:), 0.d0)

       ! second ADI step, solve thickness equation along columns i
       n = params%nsn
       do ew=2,params%ewn-1
          call adi_tri ( alpha,                 &
                         beta,                  &
                         gamma,                 &
                         delta,                 &
                         model%thckwk%oldthck(ew,:),         &
                         lsrf(ew, :),                        &
                         acab(ew, :)-real(model%options%basal_mbal)*real(model%temper%bmlt(ew, :),sp),          &
                         model%velocity%uflx(ew,:),          &
                         model%velocity%uflx(ew-1,:),        &
                         model%velocity%total_diffu(ew,:),   &
                         model%velocity%total_diffu(ew-1,:), &
                         dt,                                 &
                         params%dns,                         &
                         params%dew )

          call tridiag(alpha(1:n),    &
                       beta(1:n),     &
                       gamma(1:n),    &
                       thck(ew,:),    &
                       delta(1:n))
       end do

       thck(:,:) = max(thck(:,:), 0.d0)

       ! Apply boundary conditions
       thck(1,:) = 0.0
       thck(params%ewn,:) = 0.0
       thck(:,1) = 0.0
       thck(:,params%nsn) = 0.0

       ! Final part of basal velocity calculation
       where (params%thklim < stagthck)
          vbas = ubas *  model%geomderv%dusrfdns / stagthck
          ubas = ubas *  model%geomderv%dusrfdew / stagthck
       elsewhere
          ubas = 0.0d0
          vbas = 0.0d0
       end where

       call velo_calc_velo(model%velowk, &
            stagthck,                    &
            model%geomderv%dusrfdew,     &
            model%geomderv%dusrfdns,     &
            model%temper%flwa,           &
            model%velocity%diffu,        &
            ubas,                        &
            vbas,                        &
            model%velocity%uvel,         &
            model%velocity%vvel,         &
            model%velocity%uflx,         &
            model%velocity%vflx)
    end if

    !------------------------------------------------------------
    ! calculate upper and lower surface
    !------------------------------------------------------------
    call glide_calclsrf(thck, topg, eus, lsrf)
    model%geometry%usrf = max(0.d0,thck + lsrf)

  end subroutine stagleapthck

!---------------------------------------------------------------------------------

  subroutine adi_tri(a,b,c,d,thk,tpg,mb,flx_p,flx_m,dif_p,dif_m,dt,ds1, ds2)

    !*FD construct tri-diagonal matrix system for a column/row

    use glimmer_global, only : dp, sp

    implicit none
    
    real(dp), dimension(:), intent(out) :: a !*FD alpha (subdiagonal)
    real(dp), dimension(:), intent(out) :: b !*FD alpha (diagonal)
    real(dp), dimension(:), intent(out) :: c !*FD alpha (superdiagonal)
    real(dp), dimension(:), intent(out) :: d !*FD right-hand side
    
    real(dp), dimension(:), intent(in) :: thk   !*FD ice thickness
    real(dp), dimension(:), intent(in) :: tpg   !*FD lower surface of ice
    real(sp), dimension(:), intent(in) :: mb    !*FD mass balance
    real(dp), dimension(:), intent(in) :: flx_p !*FD flux +1/2
    real(dp), dimension(:), intent(in) :: flx_m !*FD flux -1/2
    real(dp), dimension(:), intent(in) :: dif_p !*FD diffusivity +1/2
    real(dp), dimension(:), intent(in) :: dif_m !*FD diffusivity -1/2
    
    real(dp), intent(in) :: dt !*FD time step
    real(dp), intent(in) :: ds1, ds2 !*FD spatial steps inline and transversal

    ! local variables
    real(dp) :: f1, f2, f3
    integer :: i,n
    
    n = size(thk)

    f1 = dt/(4*ds1*ds1)
    f2 = dt/(4*ds2)
    f3 = dt/2.

    a(:) = 0.
    b(:) = 0.
    c(:) = 0.
    d(:) = 0.

    a(1) = 0.
    do i=2,n
       a(i) = f1*(dif_m(i-1)+dif_p(i-1))
    end do
    do i=1,n-1
       c(i) = f1*(dif_m(i)+dif_p(i))
    end do
    c(n) = 0.
    b(:) = -(a(:)+c(:))

    ! calculate RHS
    do i=2,n-1
       d(i) = thk(i) - &
            f2 * (flx_p(i-1) + flx_p(i) - flx_m(i-1) - flx_m(i)) + &
            f3 * mb(i) - &
            a(i)*tpg(i-1) - b(i)*tpg(i) - c(i)*tpg(i+1)
    end do

    b(:) = 1.+b(:)

  end subroutine adi_tri



end module glide_thckADI
