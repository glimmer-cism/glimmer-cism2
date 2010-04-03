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

  use glimmer_global, only: dp

  implicit none

  type thckADI_type
     integer  :: ewn
     integer  :: nsn
     integer  :: upn
     real(dp) :: dew
     real(dp) :: dns
     real(dp) :: thklim
     real(dp) :: alpha
     integer  :: basal_mbal_flag
     real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: dintflwa  
     !< Vertically-integrated value of Glen's A
     real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: uflx  !< u flux 
     real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: vflx  !< v flux
     real(dp),dimension(:),  GC_DYNARRAY_ATTRIB :: velo_dups
     real(dp),dimension(:),  GC_DYNARRAY_ATTRIB :: depth
     real(dp) :: fc2_3, fc2_4 !< convenience constants 
     integer :: periodic_ew !< Set to indicate periodic BCs
  end type thckADI_type

contains

  subroutine thckADI_init(params,ewn,nsn,upn,sigma,dew,dns,thklim,alpha,basal_mbal_flag,periodic_ew)

    use glimmer_global, only: dp
    use physcon,        only: gn

    implicit none

    type(thckADI_type),intent(out) :: params
    integer,           intent(in)  :: ewn
    integer,           intent(in)  :: nsn
    integer,           intent(in)  :: upn
    real(dp),dimension(:),intent(in) :: sigma
    real(dp),          intent(in)  :: dew
    real(dp),          intent(in)  :: dns
    real(dp),          intent(in)  :: thklim
    real(dp),          intent(in)  :: alpha
    integer,           intent(in)  :: basal_mbal_flag
    integer,           intent(in)  :: periodic_ew

    integer :: up

    params%ewn = ewn
    params%nsn = nsn
    params%upn = upn
    params%dew = dew
    params%dns = dns
    params%thklim = thklim
    params%alpha = alpha
    params%basal_mbal_flag = basal_mbal_flag
    params%periodic_ew = periodic_ew

    if (GC_DYNARRAY_CHECK(params%dintflwa)) deallocate(params%dintflwa)
    if (GC_DYNARRAY_CHECK(params%uflx))     deallocate(params%uflx)
    if (GC_DYNARRAY_CHECK(params%vflx))     deallocate(params%vflx)
    if (GC_DYNARRAY_CHECK(params%velo_dups)) deallocate(params%velo_dups)
    if (GC_DYNARRAY_CHECK(params%depth))    deallocate(params%depth)

    allocate(params%dintflwa(ewn-1,nsn-1))
    allocate(params%uflx(ewn-1,nsn-1))
    allocate(params%vflx(ewn-1,nsn-1))
    allocate(params%velo_dups(upn))
    allocate(params%depth(upn))

    params%velo_dups = (/ (sigma(up+1) - sigma(up), up=1,upn-1),0.0d0 /)
    params%depth = (/ (((sigma(up+1)+sigma(up))/2.0d0)**gn *(sigma(up+1)-sigma(up)),up=1,upn-1),0.0d0 /)

    params%fc2_3 = (1.0d0 - alpha) / alpha
    params%fc2_4 =  1.0d0 / alpha

  end subroutine thckADI_init

  !--------------------------------------------------------------------

  subroutine thckADI_tstep(params,thck,acab,lsrf,usrf,topg,btrc,ubas,vbas,bmlt,flwa,uvel,vvel,diffu,dt,eus,newtemps)
    
    !*FD this subroutine solves the ice sheet thickness equation using the ADI scheme
    !*FD diffusivities are updated for each half time step

    use glide_thckCommon, only: velo_calc_velo, velo_integrate_flwa, velo_calc_diffu, glide_calclsrf
    use glimmer_utils,    only: stagvarb, tridiag
    use glimmer_global,   only: dp, sp
    use physcon,          only: rhoi, grav
    use glimmer_deriv,    only: df_field_2d_staggered

    implicit none

    ! subroutine arguments
    type(thckADI_type),                           intent(inout) :: params   !< Model parameter type
    real(dp),dimension(params%ewn,params%nsn),    intent(inout) :: thck     !< Ice thickness
    real(sp),dimension(params%ewn,params%nsn),    intent(in)    :: acab     !< Mass balance
    real(dp),dimension(params%ewn,params%nsn),    intent(inout) :: lsrf     !< Lower surface elevation
    real(dp),dimension(params%ewn,params%nsn),    intent(inout) :: usrf     !< Upper surface elevation
    real(dp),dimension(params%ewn,params%nsn),    intent(in)    :: topg     !< basal topography
    real(dp),dimension(params%ewn-1,params%nsn-1),intent(in)    :: btrc     !< basal traction coefficient 
    real(dp),dimension(params%ewn-1,params%nsn-1),intent(out)   :: ubas     !< basal x velocity
    real(dp),dimension(params%ewn-1,params%nsn-1),intent(out)   :: vbas     !< basal y velocity
    real(dp),dimension(params%ewn,params%nsn),    intent(in)    :: bmlt     !< basal melt rate
    real(dp),dimension(params%upn,params%ewn,params%nsn),    intent(in)  :: flwa !< Glen's A
    real(dp),dimension(params%upn,params%ewn-1,params%nsn-1),intent(out) :: uvel !< x-velocity
    real(dp),dimension(params%upn,params%ewn-1,params%nsn-1),intent(out) :: vvel !< y-velocity
    real(dp),dimension(params%ewn-1,params%nsn-1),intent(out)   :: diffu    !< diffusivity
    real(dp),                                     intent(in)    :: dt       !< Timestep
    real(sp),                                     intent(in)    :: eus      !< Sea level
    logical,                                      intent(in)    :: newtemps !< true when we should recalculate Glen's A

    ! local variables
    integer ew,ns, n
    real(dp), parameter :: rhograv = - rhoi * grav

    ! Local arrays
    real(dp), dimension(max(params%ewn,params%nsn)) :: alpha
    real(dp), dimension(max(params%ewn,params%nsn)) :: beta
    real(dp), dimension(max(params%ewn,params%nsn)) :: gamma
    real(dp), dimension(max(params%ewn,params%nsn)) :: delta
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: fslip
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: stagthck
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: total_diffu
    real(dp), dimension(params%ewn,  params%nsn  )  :: oldthck
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: dusrfdew
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: dusrfdns

    if (all(thck==0.0)) then

       thck = dmax1(0.0d0,thck + acab * dt)
#ifdef DEBUG       
       print *, "* thck empty - net accumulation added"
#endif
    else

       ! Calculate staggered thickness
       call stagvarb(thck,stagthck,params%ewn,params%nsn)

       ! Calculate spatial derivatives
       call df_field_2d_staggered(usrf,params%dew,params%dns,dusrfdew,dusrfdns,.false.,.false.)

       ! First part of basal velocity calculation
       fslip =  rhograv * btrc

       if (newtemps) then
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(     &
               params%velo_dups,        &
               params%depth,            &
               params%dintflwa,         &
               stagthck,                &
               flwa)
       end if

       ! Second part of velocity calculation
       where (params%thklim < stagthck)
          ubas = fslip * stagthck**2  
       elsewhere
          ubas = 0.0d0
       end where

       ! calculate diffusivity
       call velo_calc_diffu(         &
            params%dintflwa,   &
            stagthck, &
            dusrfdew, &
            dusrfdns, &
            diffu)

       total_diffu(:,:) = diffu(:,:) + ubas(:,:)

       ! first ADI step, solve thickness equation along rows j
       n = params%ewn
       do ns=2,params%nsn-1
          call adi_tri ( alpha,                 &
                         beta,                  &
                         gamma,                 &
                         delta,                 &
                         thck(:,ns),            &
                         lsrf(:,ns),            &
                         acab(:,ns)-real(params%basal_mbal_flag)*real(bmlt(:,ns),sp),           &
                         params%vflx(:,ns),          &
                         params%vflx(:,ns-1),        &
                         total_diffu(:,ns),                  &
                         total_diffu(:,ns-1),                &
                         dt,                                 &
                         params%dew,                         &
                         params%dns )

          call tridiag(alpha(1:n),    &
                       beta(1:n),     &
                       gamma(1:n),    &
                       oldthck(:,ns), &
                       delta(1:n))
       end do

       oldthck(:,:) = max(oldthck(:,:), 0.d0)

       ! second ADI step, solve thickness equation along columns i
       n = params%nsn
       do ew=2,params%ewn-1
          call adi_tri ( alpha,                 &
                         beta,                  &
                         gamma,                 &
                         delta,                 &
                         oldthck(ew,:),         &
                         lsrf(ew, :),                        &
                         acab(ew, :)-real(params%basal_mbal_flag)*real(bmlt(ew, :),sp),          &
                         params%uflx(ew,:),          &
                         params%uflx(ew-1,:),        &
                         total_diffu(ew,:),                  &
                         total_diffu(ew-1,:),                &
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
          vbas = ubas *  dusrfdns / stagthck
          ubas = ubas *  dusrfdew / stagthck
       elsewhere
          ubas = 0.0d0
          vbas = 0.0d0
       end where

       call velo_calc_velo( &
            params%dintflwa,             &
            params%depth,                &
            stagthck,                    &
            dusrfdew,                    &
            dusrfdns,                    &
            flwa,                        &
            diffu,                       &
            ubas,                        &
            vbas,                        &
            uvel,                        &
            vvel,                        &
            params%uflx,                 &
            params%vflx)
    end if

    !------------------------------------------------------------
    ! calculate upper and lower surface
    !------------------------------------------------------------
    call glide_calclsrf(thck, topg, eus, lsrf)
    usrf = max(0.d0,thck + lsrf)

  end subroutine thckADI_tstep

!---------------------------------------------------------------------------------

  subroutine get_uflx(params,uflx)

    use glimmer_global, only: dp
    implicit none

    type(thckADI_type),                           intent(in)  :: params
    real(dp),dimension(params%ewn-1,params%nsn-1),intent(out) :: uflx !< flux in x (note that these are 
    uflx = params%uflx

  end subroutine get_uflx

!---------------------------------------------------------------------------------

  subroutine get_vflx(params,vflx)
    
    use glimmer_global, only: dp
    implicit none

    type(thckADI_type),                           intent(in)  :: params
    real(dp),dimension(params%ewn-1,params%nsn-1),intent(out) :: vflx !< flux in y  for information only)
    vflx = params%vflx

  end subroutine get_vflx

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
