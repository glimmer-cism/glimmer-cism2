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

! some macros used to disable parts of the temperature equation
! vertical diffusion
#ifdef NO_VERTICAL_DIFFUSION
#define VERT_DIFF 0.
#else
#define VERT_DIFF 1.
#endif
! horizontal advection
#ifdef NO_HORIZONTAL_ADVECTION
#define HORIZ_ADV 0.
#else
#define HORIZ_ADV 1.
#endif
! vertical advection
#ifdef NO_VERICAL_ADVECTION
#define VERT_ADV 0.
#else
#define VERT_ADV 1.
#endif
! strain heating
#ifdef NO_STRAIN_HEAT
#define STRAIN_HEAT 0.
#else
#define STRAIN_HEAT 1.
#endif

module glide_temp

  use glide_types

  private :: find_dt_wat

contains

  subroutine init_temp(model)
    !*FD initialise temperature module
    use physcon, only : rhoi, shci, coni, scyr, grav, gn, lhci, rhow
    use glimmer_paramets, only : tim0, thk0, acc0, len0, vis0, vel0
    use glide_vertcoord, only: initVertCoord
    use glimmer_global, only : dp 
    use glimmer_log
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD Ice model parameters.

    integer up
    real(dp) :: estimate

    if (VERT_DIFF.eq.0.) call write_log('Vertical diffusion is switched off')
    if (HORIZ_ADV.eq.0.) call write_log('Horizontal advection is switched off')
    if (VERT_ADV.eq.0.) call write_log('Vertical advection is switched off')
    if (STRAIN_HEAT.eq.0.) call write_log('Strain heating is switched off')

    ! horizontal advection stuff
    allocate(model%tempwk%hadv_u(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%hadv_v(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%initadvt(model%general%upn,model%general%ewn,model%general%nsn))

    allocate(model%tempwk%inittemp(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%dissip(model%general%upn,model%general%ewn,model%general%nsn))

    allocate(model%tempwk%smth(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%wphi(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bwatu(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bwatv(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%fluxew(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%fluxns(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bint(model%general%ewn-1,model%general%nsn-1))

    call initVertCoord(model%zCoord,model%numerics%sigma)

    select case(model%options%whichbwat)
       case(0)
          model%paramets%hydtim = tim0 / (model%paramets%hydtim * scyr)
          estimate = 0.2d0 / model%paramets%hydtim
          call find_dt_wat(model%numerics%dttem,estimate,model%tempwk%dt_wat,model%tempwk%nwat) 
          
          model%tempwk%c = (/ model%tempwk%dt_wat, 1.0d0 - 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, &
               1.0d0 + 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /) 
       case(1)
          model%tempwk%watvel = model%paramets%hydtim * tim0 / (scyr * len0)
          estimate = (0.2d0 * model%tempwk%watvel) / min(model%numerics%dew,model%numerics%dns)
          call find_dt_wat(model%numerics%dttem,estimate,model%tempwk%dt_wat,model%tempwk%nwat) 
          
          print *, model%numerics%dttem*tim0/scyr, model%tempwk%dt_wat*tim0/scyr, model%tempwk%nwat

          model%tempwk%c = (/ rhow * grav, rhoi * grav, 2.0d0 * model%numerics%dew, 2.0d0 * model%numerics%dns, &
               0.25d0 * model%tempwk%dt_wat / model%numerics%dew, 0.25d0 * model%tempwk%dt_wat / model%numerics%dns, &
               0.5d0 * model%tempwk%dt_wat / model%numerics%dew, 0.5d0 * model%tempwk%dt_wat / model%numerics%dns /)
          
       end select
  end subroutine init_temp

  !------------------------------------------------------------------------------------

  subroutine calcTemp_asSurfTemp(model)

    !*FD Sets ice column temperature to surface temperature

    type(glide_global_type),intent(inout) :: model       !*FD Ice model parameters.

    integer :: ns,ew

    do ns = 1,model%general%nsn
       do ew = 1,model%general%ewn
          model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(model%climate%artm(ew,ns)))
       end do
    end do

  end subroutine calcTemp_asSurfTemp

  !------------------------------------------------------------------------------------

  subroutine calcTemp_VerticalProfile(model)

    !*FD Sets ice column temperature to surface temperature at top
    !*FD and zero at bottom, correcting for pressure melting point

    type(glide_global_type),intent(inout) :: model       !*FD Ice model parameters.

    integer :: ns,ew

    do ns = 1,model%general%nsn
       do ew = 1,model%general%ewn
          model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(model%climate%artm(ew,ns))) * (1.0d0 - model%numerics%sigma)
          call corrpmpt(model%temper%temp(:,ew,ns),model%geometry%thck(ew,ns),model%temper%bwat(ew,ns),&
               model%numerics%sigma,model%general%upn)
       end do
    end do

  end subroutine calcTemp_VerticalProfile

  !------------------------------------------------------------------------------------

  subroutine calcTemp_FullSolution(model,dt)

    !*FD Calculates the ice temperature - full solution

    use glimmer_utils,    only: hsum4,tridiag
    use glimmer_global,   only: dp
    use glimmer_paramets, only: thk0, acc0, tim0, len0, vis0, scyr
    use glide_thck,       only: stagvarb
    use glide_mask,       only: is_float, is_thin
    use physcon,          only: rhoi, shci, coni, grav, gn, lhci

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_global_type),intent(inout) :: model       !*FD Ice model parameters.
    real(dp),               intent(in)    :: dt          !*FD Timestep (years)

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp),dimension(size(model%numerics%sigma)) :: subd, diag, supd, rhsd
    real(dp),dimension(size(model%numerics%sigma)) :: prevtemp, iteradvt, diagadvt
    real(dp) :: tempresid

    integer :: iter
    integer :: ew,ns

    real(dp),parameter :: tempthres = 0.001d0, floatlim = 10.0d0 / thk0
    integer, parameter :: mxit = 100
    integer, parameter :: ewbc = 1, nsbc = 1 
    real(dp),parameter :: wmax = 5.0d0 * tim0 / (scyr * thk0)

    real(dp), dimension(size(model%numerics%sigma)) :: weff

    real(dp) :: advconst1,advconst2
    real(dp) :: cons1,cons2,cons3,cons4
    real(dp),dimension(size(model%numerics%sigma)) :: c1
    real(dp) :: f3
    real(dp) :: slidef1, slidef2

    !------------------------------------------------------------------------------------
    ! Set up various parameters
    !------------------------------------------------------------------------------------

    advconst1 = HORIZ_ADV * dt / (16.0d0 * model%numerics%dew)
    advconst2 = HORIZ_ADV * dt / (16.0d0 * model%numerics%dns)

    cons1 = 2.0d0 * tim0 * dt * coni / (2.0d0 * rhoi * shci * thk0**2)
    cons2 = dt / 2.0d0
    cons3 = VERT_DIFF * 2.0d0 * tim0 * dt / (thk0 * rhoi * shci)
    cons4 = VERT_ADV  * tim0 * acc0 * dt / coni

    c1 = STRAIN_HEAT * (model%numerics%sigma * rhoi * grav * thk0**2 / len0)**(gn+1) * &
         2.0d0 * vis0 * dt * tim0 / (16.0d0 * rhoi * shci)

    f3 = tim0 * thk0 * rhoi * shci /  (thk0 * tim0 * dt * lhci * rhoi)

    ! sliding contribution to basal heat flux

    slidef1 = VERT_DIFF * grav * thk0 * dt / shci                      ! vert diffusion
    slidef2 = VERT_ADV  * rhoi * grav * acc0 * thk0 * thk0 * dt / coni ! vert advection

    !------------------------------------------------------------------------------------
    ! ewbc/nsbc set the type of boundary condition applied at the end of
    ! the domain. a value of 0 implies zero gradient.
    !------------------------------------------------------------------------------------

    model%tempwk%inittemp = 0.0d0
    model%tempwk%initadvt = 0.0d0

    ! Calculate dissipative term --------------------------------------------------------

    call finddisp(model%tempwk%dissip, &
         model%geometry%thck,          &
         model%geomderv%stagthck,      &
         model%geomderv%dusrfdew,      &
         model%geomderv%dusrfdns,      &
         model%temper%flwa,            &
         model%general%ewn,            &
         model%general%nsn,            &
         c1,                           &
         model%numerics%thklim)

    ! translate velo field --------------------------------------------------------------

    do ns = 2,model%general%nsn-1
       do ew = 2,model%general%ewn-1
          model%tempwk%hadv_u(:,ew,ns) = advconst1 * hsum4(model%velocity%uvel(:,ew-1:ew,ns-1:ns))
          model%tempwk%hadv_v(:,ew,ns) = advconst2 * hsum4(model%velocity%vvel(:,ew-1:ew,ns-1:ns))
       end do
    end do

    ! Calculate initial upwinding terms -------------------------------------------------

    call hadvall(model%tempwk%initadvt, &
         model%temper%temp,             &
         model%geometry%thck,           &
         model%numerics%thklim,         &
         model%tempwk%hadv_u,           &
         model%tempwk%hadv_v,           &
         model%general%ewn,             &
         model%general%nsn,             &
         model%general%upn)

    ! Iterative temperature solution ----------------------------------------------------

    iter = 0
    tempresid = abs(tempthres*2.0)  ! To make sure the loop is executed at least once
   
    do while (tempresid.gt.tempthres .and. iter.le.mxit)
       tempresid = 0.0d0

       do ns = 2,model%general%nsn-1
          do ew = 2,model%general%ewn-1
             if(model%geometry%thck(ew,ns)>model%numerics%thklim) then

                ! Calculate effective vertical velocity
                weff = model%velocity%wvel(:,ew,ns) - model%velocity%wgrd(:,ew,ns)

                ! Set effective vertical velocity to zero if it exceeds a threshold
                if (maxval(abs(weff)) > wmax) then
                   weff = 0.0d0
                end if

                ! Calculate upwinded advection term
                call hadvpnt(iteradvt,                       &
                     diagadvt,                               &
                     model%temper%temp(:,ew-2:ew+2,ns),      &
                     model%temper%temp(:,ew,ns-2:ns+2),      &
                     model%tempwk%hadv_u(:,ew,ns),           &
                     model%tempwk%hadv_v(:,ew,ns),           &
                     model%general%upn)

                call findvtri(model%zCoord,      &
                     model%geometry%thck(ew,ns), &
                     subd,                       &
                     diag,                       &
                     supd,                       &
                     diagadvt,                   &
                     weff,                       &
                     is_float(model%geometry%thkmask(ew,ns)), &
                     model%general%upn,          &
                     cons1,                      &
                     cons2)

                if (iter==0) then
                   call findvtri_init(model%tempwk, &
                        model%tempwk%inittemp,      &
                        model%temper%bheatflx,      &
                        model%geomderv%dusrfdew,    &
                        model%geomderv%dusrfdns,    &
                        model%zCoord,               &
                        ew,                         &
                        ns,                         &
                        subd,                       &
                        diag,                       &
                        supd,                       &
                        weff,                       &
                        model%velocity%ubas,        &
                        model%velocity%vbas,        &
                        model%temper%temp(:,ew,ns), &
                        model%geometry%thck(ew,ns), &
                        is_float(model%geometry%thkmask(ew,ns)), &
                        model%general%upn,          &
                        cons3,                      &
                        cons4,                      &
                        slidef1,                    &
                        slidef2)
                end if

                call findvtri_rhs(model%tempwk%inittemp(:,ew,ns), &
                     model%climate%artm(ew,ns),                   &
                     iteradvt,                                    &
                     rhsd,                                        &
                     is_float(model%geometry%thkmask(ew,ns)),     &
                     model%general%upn)

                prevtemp = model%temper%temp(:,ew,ns)

                call tridiag(subd(1:model%general%upn),            &
                     diag(1:model%general%upn),                    &
                     supd(1:model%general%upn),                    &
                     model%temper%temp(1:model%general%upn,ew,ns), &
                     rhsd(1:model%general%upn))

                call corrpmpt(model%temper%temp(:,ew,ns),    &
                     model%geometry%thck(ew,ns),             &
                     model%temper%bwat(ew,ns),               &
                     model%numerics%sigma,                   &
                     model%general%upn)

                tempresid = max(tempresid,maxval(abs(model%temper%temp(:,ew,ns)-prevtemp(:))))
             endif
          end do
       end do

       iter = iter + 1
    end do

    model%temper%niter = max(model%temper%niter, iter )
       
    ! set temperature of thin ice to the air temperature and set ice free nodes to zero

    do ns = 1,model%general%nsn
       do ew = 1,model%general%ewn
          if (is_thin(model%geometry%thkmask(ew,ns))) then
             model%temper%temp(:,ew,ns) = min(0.0d0,dble(model%climate%artm(ew,ns)))
          else if (model%geometry%thkmask(ew,ns)<0) then
             model%temper%temp(:,ew,ns) = 0.0d0
          end if
       end do
    end do

    ! apply periodic ew BC ----------------------------------------------------------

    if (model%options%periodic_ew.eq.1) then
       model%temper%temp(:,0,:) = model%temper%temp(:,model%general%ewn-2,:)
       model%temper%temp(:,1,:) = model%temper%temp(:,model%general%ewn-1,:)
       model%temper%temp(:,model%general%ewn,:) = model%temper%temp(:,2,:)
       model%temper%temp(:,model%general%ewn+1,:) = model%temper%temp(:,3,:)
    end if

    ! Calculate basal melt rate --------------------------------------------------

    call calcbmlt(model%tempwk,           &
         model%zCoord,                    &
         model%temper%temp,               &
         model%geometry%thck,             &
         model%geomderv%stagthck,         &
         model%geomderv%dusrfdew,         &
         model%geomderv%dusrfdns,         &
         model%velocity%ubas,             &
         model%velocity%vbas,             &
         model%temper%bheatflx,           &
         model%temper%bmlt,               &
         is_float(model%geometry%thkmask),&
         model%general%upn,               &
         model%general%nsn,               &
         model%general%ewn,               &
         model%numerics%thklim,           &
         model%numerics%sigma,            &
         model%options%periodic_ew,       &
         f3)

    ! Calculate basal water depth ------------------------------------------------

    call calcbwat(model,                           &
         model%options%whichbwat,                  &
         model%temper%bmlt,                        &
         model%temper%bwat,                        &
         model%geometry%thck,                      &
         model%geometry%topg,                      &
         model%temper%temp(model%general%upn,:,:), &
         is_float(model%geometry%thkmask))

    ! Transform basal temperature and pressure melting point onto velocity grid -

    call stagvarb(model%temper%temp(model%general%upn,1:model%general%ewn,1:model%general%nsn), &
         model%temper%stagbtemp ,  &
         model%general%  ewn,      &
         model%general%  nsn)
       
    call calcbpmp(model,model%geometry%thck,model%temper%bpmp)

    call stagvarb(model%temper%bpmp, &
         model%temper%stagbpmp ,     &
         model%general%ewn,          &
         model%general%nsn)

    ! Output some information ----------------------------------------------------
#ifdef DEBUG
    print *, "* temp ", model%numerics%time, iter, model%temper%niter, &
         real(model%temper%temp(model%general%upn,model%general%ewn/2+1,model%general%nsn/2+1))
#endif

  end subroutine calcTemp_FullSolution

  !-------------------------------------------------------------------------

  subroutine hadvpnt(iteradvt,diagadvt,tempx,tempy,u,v,upn)

    use glimmer_global, only : dp

    implicit none

    integer,                    intent(in)  :: upn        ! Number of points in vertical
    real(dp), dimension(upn),   intent(out) :: iteradvt
    real(dp), dimension(upn),   intent(out) :: diagadvt
    real(dp), dimension(upn,5), intent(in)  :: tempx
    real(dp), dimension(upn,5), intent(in)  :: tempy
    real(dp), dimension(upn),   intent(in)  :: u
    real(dp), dimension(upn),   intent(in)  :: v

    iteradvt = 0.0d0
    diagadvt = 0.0d0

    if (u(1) > 0.0d0) then
       iteradvt = u * (- 4.0d0*tempx(:,2) + tempx(:,1))
       diagadvt = u * 3.0d0
    else if (u(1) < 0.0d0) then
       iteradvt = u * (4.0d0*tempx(:,4) - tempx(:,5))
       diagadvt = - u * 3.0d0
    end if

    if (v(1) > 0.0d0) then
       iteradvt = iteradvt + v * (- 4.0d0*tempy(:,2) + tempy(:,1))
       diagadvt = diagadvt + v * 3.0d0
    else if (v(1) < 0.0d0) then
       iteradvt = iteradvt + v * (4.0d0*tempy(:,4) - tempy(:,5))
       diagadvt = diagadvt - v * 3.0d0
    end if

  end subroutine hadvpnt

  !-------------------------------------------------------------------------

  subroutine fohadvpnt(advconst1,advconst2,iteradvt,diagadvt,tempx,tempy,uvel,vvel)

    use glimmer_global, only: dp
    use glimmer_utils,  only: hsum

    implicit none

    real(dp), intent(in) :: advconst1
    real(dp), intent(in) :: advconst2
    real(dp), dimension(:,:,:), intent(in) :: uvel, vvel
    real(dp), dimension(:,:), intent(in) :: tempx, tempy
    real(dp), dimension(:), intent(out) :: iteradvt, diagadvt

    real(dp), dimension(size(iteradvt)) :: u, v

    iteradvt = 0.0d0
    diagadvt = 0.0d0

    u = advconst1 * hsum(uvel(:,:,:))
    v = advconst2 * hsum(vvel(:,:,:))

    if (u(1) > 0.0d0) then
       iteradvt = - u * 2.0d0 * tempx(:,1)
       diagadvt = 2.0d0 * u 
    else if (u(1) < 0.0d0) then
       iteradvt = u * 2.0d0 * tempx(:,3)
       diagadvt = - 2.0d0 * u 
    end if

    if (v(1) > 0.0d0) then
       iteradvt = iteradvt - v * 2.0d0 * tempy(:,1) 
       diagadvt = diagadvt + 2.0d0 * v 
    else if (v(1) < 0.0d0) then
       iteradvt = iteradvt + v * 2.0d0 * tempy(:,3)
       diagadvt = diagadvt - 2.0d0 * v 
    end if

  end subroutine fohadvpnt

  !-------------------------------------------------------------------------

  subroutine hadvall(initadvt,temp,thck,thklim,hadv_u,hadv_v,ewn,nsn,upn)

    use glimmer_global, only : dp 

    implicit none

    integer,                                     intent(in)  :: ewn
    integer,                                     intent(in)  :: nsn
    integer,                                     intent(in)  :: upn
    real(dp), dimension(upn,ewn,nsn),            intent(out) :: initadvt
    real(dp), dimension(upn,0:ewn+1,0:nsn+1),    intent(in)  :: temp
    real(dp), dimension(ewn,nsn),                intent(in)  :: thck
    real(dp),                                    intent(in)  :: thklim
    real(dp), dimension(upn,ewn,nsn),            intent(in)  :: hadv_u
    real(dp), dimension(upn,ewn,nsn),            intent(in)  :: hadv_v

    real(dp), dimension(upn) :: diagadvt

    integer :: ew,ns

    initadvt = 0.0d0

    do ns = 2,nsn-1
       do ew = 2,ewn-1
          if (thck(ew,ns) > thklim) then

             call hadvpnt(initadvt(:,ew,ns),      &
                  diagadvt,                       &
                  temp(:,ew-2:ew+2,ns),           &
                  temp(:,ew,ns-2:ns+2),           &
                  hadv_u(:,ew,ns),                &
                  hadv_v(:,ew,ns),                &
                  upn)
          end if
       end do
    end do

  end subroutine hadvall

  !-------------------------------------------------------------------------

  subroutine findvtri(zCoord,thck,subd,diag,supd,diagadvt,weff,float,upn,cons1,cons2)

    use glide_vertcoord, only: vertCoord
    use glimmer_global,  only: dp

    implicit none

    type(vertCoord),        intent(in)  :: zCoord
    real(dp),               intent(in)  :: thck
    real(dp), dimension(:), intent(in)  :: weff
    real(dp), dimension(:), intent(in)  :: diagadvt
    real(dp), dimension(:), intent(out) :: subd
    real(dp), dimension(:), intent(out) :: diag
    real(dp), dimension(:), intent(out) :: supd
    logical,                intent(in)  :: float
    integer,                intent(in)  :: upn
    real(dp),               intent(in)  :: cons1
    real(dp),               intent(in)  :: cons2

    real(dp) :: diff_factor,adv_factor

    diff_factor = VERT_DIFF * cons1 / thck**2
    adv_factor  = VERT_ADV  * cons2 / thck
    
    subd(2:upn-1) =   adv_factor  * weff(2:upn-1) * zCoord%dups(2:upn-1,3)
    supd(2:upn-1) = - subd(2:upn-1) - diff_factor * zCoord%dups(2:upn-1,2)
    subd(2:upn-1) =   subd(2:upn-1) - diff_factor * zCoord%dups(2:upn-1,1)

    diag(2:upn-1) = 1.0d0 - subd(2:upn-1) - supd(2:upn-1) + diagadvt(2:upn-1)

    ! Upper surface: hold temperature constant

    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0

    ! now do the basal boundary
    ! for grounded ice, a heat flux is applied
    ! for floating ice, temperature held constant

    if (float) then
       supd(upn) = 0.0d0
       subd(upn) = 0.0d0
       diag(upn) = 1.0d0
    else 
       supd(upn) = 0.0d0 
       subd(upn) = -0.5*diff_factor/(zCoord%dupn**2)
       diag(upn) = 1.0d0 - subd(upn) + diagadvt(upn)
    end if

  end subroutine findvtri

  !-------------------------------------------------------------------------

  subroutine findvtri_init(tempwk,inittemp,bheatflx,dusrfdew,dusrfdns, &
       zCoord,ew,ns,subd,diag,supd,weff,ubas,vbas,temp,thck,float,upn, &
       cons3,cons4,slidef1,slidef2)

    !*FD called during first iteration to set inittemp

    use glimmer_global, only : dp

    implicit none

    type(glide_tempwk),       intent(in)  :: tempwk
    real(dp),dimension(:,:,:),intent(out) :: inittemp
    real(dp),dimension(:,:),  intent(in)  :: bheatflx
    real(dp),dimension(:,:),  intent(in)  :: dusrfdew
    real(dp),dimension(:,:),  intent(in)  :: dusrfdns
    type(vertCoord),          intent(in)  :: zCoord
    integer,                  intent(in)  :: ew
    integer,                  intent(in)  :: ns
    real(dp),dimension(:),    intent(in)  :: temp
    real(dp),dimension(:),    intent(in)  :: diag
    real(dp),dimension(:),    intent(in)  :: subd
    real(dp),dimension(:),    intent(in)  :: supd
    real(dp),dimension(:),    intent(in)  :: weff
    real(dp),dimension(:,:),  intent(in)  :: ubas
    real(dp),dimension(:,:),  intent(in)  :: vbas
    real(dp),                 intent(in)  :: thck
    logical,                  intent(in)  :: float    
    integer,                  intent(in)  :: upn
    real(dp),                 intent(in)  :: cons3
    real(dp),                 intent(in)  :: cons4
    real(dp),                 intent(in)  :: slidef1
    real(dp),                 intent(in)  :: slidef2

    ! local variables
    real(dp) :: slterm
    integer  :: ewp,nsp
    integer  :: slide_count

    ! Main body points
    inittemp(2:upn-1,ew,ns) = temp(2:upn-1) * (2.0d0 - diag(2:upn-1)) &
         - temp(1:upn-2) * subd(2:upn-1)         &
         - temp(3:upn)   * supd(2:upn-1)         & 
         - tempwk%initadvt(2:upn-1,ew,ns)        &
         + tempwk%dissip(2:upn-1,ew,ns)
    
    ! Basal boundary points
    if (float) then
       inittemp(upn,ew,ns) = pmpt(thck)
    else 
       ! sliding contribution to basal heat flux
       slterm = 0.
       slide_count = 0
       ! only include sliding contrib if temperature node is surrounded by sliding velo nodes
       do nsp = ns-1,ns
          do ewp = ew-1,ew
             if (abs(ubas(ewp,nsp)).gt.0.000001 .or. abs(vbas(ewp,nsp)).gt.0.000001) then
                slide_count = slide_count + 1
                slterm = slterm + (dusrfdew(ewp,nsp)*ubas(ewp,nsp) + dusrfdns(ewp,nsp)*vbas(ewp,nsp))
             end if
          end do
       end do

       if (slide_count.ge.4) then
          slterm = 0.25*slterm
       else
          slterm = 0.
       end if

       inittemp(upn,ew,ns) = temp(upn) * (2.0d0 - diag(upn)) &
            - temp(upn-1) * subd(upn)                        &
            - 0.5 * cons3 * bheatflx(ew,ns) / (thck * zCoord%dupn) &  ! geothermal heat flux (diff)
            - slidef1 * slterm / zCoord%dupn                 &        ! sliding heat flux    (diff)
            - cons4 * bheatflx(ew,ns) * weff(upn)            &        ! geothermal heat flux (adv)
            - slidef2 * thck * slterm * weff(upn)            &        ! sliding heat flux    (adv)
            - tempwk%initadvt(upn,ew,ns)                     &
            + tempwk%dissip(upn,ew,ns)

    end if

  end subroutine findvtri_init

  !-------------------------------------------------------------------------

  subroutine findvtri_rhs(inittemp,artm,iteradvt,rhsd,float,upn)

    !*FD RHS of temperature tri-diag system for a single column

    use glimmer_global, only : dp, sp 

    implicit none

    real(dp),dimension(:),intent(in)  :: inittemp
    real(sp),             intent(in)  :: artm 
    real(dp),dimension(:),intent(in)  :: iteradvt
    real(dp),dimension(:),intent(out) :: rhsd
    logical,              intent(in)  :: float
    integer,              intent(in)  :: upn

    ! upper boundary condition

    rhsd(1) = artm

    if (float) then
       rhsd(upn) = inittemp(upn)    
    else
       rhsd(upn) = inittemp(upn) - iteradvt(upn)
    end if

    rhsd(2:upn-1) = inittemp(2:upn-1) - iteradvt(2:upn-1)

  end subroutine findvtri_rhs

  !-----------------------------------------------------------------------

  subroutine finddisp(dissip,thck,stagthck,dusrfdew,dusrfdns,flwa,ewn,nsn,c1,thklim)

    use glimmer_global, only : dp
    use physcon, only : gn

    implicit none

    real(dp), dimension(:,:,:), intent(out) :: dissip
    real(dp), dimension(:,:),   intent(in)  :: thck
    real(dp), dimension(:,:),   intent(in)  :: stagthck
    real(dp), dimension(:,:),   intent(in)  :: dusrfdew
    real(dp), dimension(:,:),   intent(in)  :: dusrfdns
    real(dp), dimension(:,:,:), intent(in)  :: flwa
    integer,                    intent(in)  :: ewn,nsn
    real(dp), dimension(:),     intent(in)  :: c1
    real(dp),                   intent(in)  :: thklim

    integer, parameter :: p1 = gn + 1  
    integer :: ew,ns

    real(dp) :: c2

    ! find dissipation term at H-pts by averaging quantities from u-pts

    dissip = 0.0d0
    
    do ns = 2, nsn-1
       do ew = 2, ewn-1
          if (thck(ew,ns) > thklim) then
             
             c2 = (0.25*sum(stagthck(ew-1:ew,ns-1:ns)) * dsqrt((0.25*sum(dusrfdew(ew-1:ew,ns-1:ns)))**2 &
                  + (0.25*sum(dusrfdns(ew-1:ew,ns-1:ns)))**2))**p1
             
             dissip(:,ew,ns) = c2 * c1 * ( &
                  flwa(:,ew-1,ns-1) + flwa(:,ew-1,ns+1) + flwa(:,ew+1,ns+1) + flwa(:,ew+1,ns-1) + &
                  2*(flwa(:,ew-1,ns)+flwa(:,ew+1,ns)+flwa(:,ew,ns-1)+flwa(:,ew,ns+1)) + &
                  4*flwa(:,ew,ns))             
          end if
       end do
    end do

  end subroutine finddisp

  !-----------------------------------------------------------------------------------

  real(dp) function pmpt(thck)

    !*FD Wrapper function for pressure-melting-point calculation

    real(dp),intent(in) :: thck
    real(dp) :: tmp

    call calcpmptb(tmp,thck)
    pmpt=tmp

  end function pmpt

  !-----------------------------------------------------------------------------------

  subroutine calcbpmp(model,thck,bpmp)

    ! Calculate the pressure melting point at the base of the ice sheet

    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(in)  :: thck
    real(dp), dimension(:,:), intent(out) :: bpmp

    integer :: ew,ns

    bpmp = 0.0

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1
          call calcpmptb(bpmp(ew,ns),thck(ew,ns))
       end do
    end do

  end subroutine calcbpmp

  !-----------------------------------------------------------------------------------

  subroutine calcbmlt(tempwk,zCoord,temp,thck,stagthck,dusrfdew,dusrfdns,ubas,vbas,bheatflx, &
       bmlt,floater,upn,nsn,ewn,thklim,sigma,periodic_ew,f3)

    use glimmer_global,   only: dp
    use glimmer_paramets, only: thk0, tim0, vel0, len0
    use physcon,          only: coni, lhci, rhoi, grav

    implicit none 

    type(glide_tempwk),           intent(in)  :: tempwk
    type(vertCoord),              intent(in)  :: zCoord
    real(dp), dimension(:,0:,0:), intent(in)  :: temp
    real(dp), dimension(:,:),     intent(in)  :: thck
    real(dp), dimension(:,:),     intent(in)  :: stagthck
    real(dp), dimension(:,:),     intent(in)  :: dusrfdew
    real(dp), dimension(:,:),     intent(in)  :: dusrfdns
    real(dp), dimension(:,:),     intent(in)  :: ubas
    real(dp), dimension(:,:),     intent(in)  :: vbas
    real(dp), dimension(:,:),     intent(in)  :: bheatflx
    real(dp), dimension(:,:),     intent(out) :: bmlt
    logical,  dimension(:,:),     intent(in)  :: floater
    integer,                      intent(in)  :: upn
    integer,                      intent(in)  :: nsn
    integer,                      intent(in)  :: ewn
    real(dp),                     intent(in)  :: thklim
    real(dp), dimension(:),       intent(in)  :: sigma
    integer,                      intent(in)  :: periodic_ew
    real(dp),                     intent(in)  :: f3
    

    real(dp), dimension(upn) :: pmptemp
    real(dp) :: slterm, newmlt

    integer :: ewp, nsp,up,ew,ns

    real(dp),parameter :: f1 = tim0 * coni / (thk0**2 * lhci * rhoi)
    real(dp),parameter :: f2 = tim0 / (thk0 * lhci * rhoi)
    real(dp),parameter :: f4 = tim0 * thk0**2 * vel0 * grav * rhoi / (4.0d0 * thk0 * len0 * rhoi * lhci)

    do ns = 2, nsn-1
       do ew = 2, ewn-1
          if (thck(ew,ns) > thklim .and. .not. floater(ew,ns)) then

             call calcpmpt(pmptemp,thck(ew,ns),sigma)

             if (abs(temp(upn,ew,ns)-pmptemp(upn)) .lt. 0.001) then

                slterm = 0.0d0

                do nsp = ns-1,ns
                   do ewp = ew-1,ew
                      slterm = slterm - stagthck(ewp,nsp) * &
                           (dusrfdew(ewp,nsp) * ubas(ewp,nsp) + dusrfdns(ewp,nsp) * vbas(ewp,nsp))
                   end do
                end do

                bmlt(ew,ns) = 0.0d0
                newmlt = tempwk%f(4) * slterm - tempwk%f(2)*bheatflx(ew,ns) + tempwk%f(3) * &
                     zCoord%dupc(upn) * thck(ew,ns) * tempwk%dissip(upn,ew,ns)

                up = upn - 1

                do while (abs(temp(up,ew,ns)-pmptemp(up)) .lt. 0.001 .and. up .ge. 3)
                   bmlt(ew,ns) = bmlt(ew,ns) + newmlt
                   newmlt = f3 * zCoord%dupc(up) * thck(ew,ns) * tempwk%dissip(up,ew,ns)
                   up = up - 1
                end do

                up = up + 1

                if (up == upn) then
                   bmlt(ew,ns) = newmlt - &
                        tempwk%f(1) * ( (temp(up-2,ew,ns) - pmptemp(up-2)) * zCoord%dupa(up) &
                        + (temp(up-1,ew,ns) - pmptemp(up-1)) * zCoord%dupb(up) ) / thck(ew,ns) 
                else
                   bmlt(ew,ns) = bmlt(ew,ns) + max(0.0d0, newmlt - &
                        tempwk%f(1) * ( (temp(up-2,ew,ns) - pmptemp(up-2)) * zCoord%dupa(up) &
                        + (temp(up-1,ew,ns) - pmptemp(up-1)) * zCoord%dupb(up) ) / thck(ew,ns)) 
                end if

             else

                bmlt(ew,ns) = 0.0d0

             end if

          else

             bmlt(ew,ns) = 0.0d0

          end if
       end do
    end do

    ! apply periodic BC
    if (periodic_ew.eq.1) then
       do ns = 2,nsn-1
          bmlt(1,ns) = bmlt(ewn-1,ns)
          bmlt(ewn,ns) = bmlt(2,ns)
       end do
    end if
  end subroutine calcbmlt

  !-------------------------------------------------------------------

  subroutine calcbwat(model,which,bmlt,bwat,thck,topg,btem,floater)

    use glimmer_global, only : dp 
    use glimmer_paramets, only : thk0
    use glide_thck
    implicit none

    type(glide_global_type) :: model
    integer, intent(in) :: which

    real(dp), dimension(:,:), intent(inout) :: bwat
    real(dp), dimension(:,:), intent(in) :: bmlt, thck, topg, btem
    logical, dimension(:,:), intent(in) :: floater

    real(dp), dimension(2), parameter :: &
         blim = (/ 0.00001 / thk0, 0.001 / thk0 /)

    real(dp) :: dwphidew, dwphidns, dwphi, pmpt, bave

    integer :: t_wat,ns,ew

    select case (which)
    case(0)

       do t_wat = 1, model%tempwk%nwat
          do ns = 1,model%general%nsn
             do ew = 1,model%general%ewn

                if (model%numerics%thklim < thck(ew,ns) .and. .not. floater(ew,ns)) then
                   bwat(ew,ns) = (model%tempwk%c(1) * bmlt(ew,ns) + model%tempwk%c(2) * bwat(ew,ns)) / &
                        model%tempwk%c(3)
                   if (blim(1) > bwat(ew,ns)) then
                      bwat(ew,ns) = 0.0d0
                   end if
                else
                   bwat(ew,ns) = 0.0d0
                end if

             end do
          end do
       end do

       model%tempwk%smth = 0.
       do ns = 2,model%general%nsn-1
          do ew = 2,model%general%ewn-1
             call smooth_bwat(ew-1,ew,ew+1,ns-1,ns,ns+1)
          end do
       end do
       ! apply periodic BC
       if (model%options%periodic_ew.eq.1) then
          do ns = 2,model%general%nsn-1
             call smooth_bwat(model%general%ewn-1,1,2,ns-1,ns,ns+1)
             call smooth_bwat(model%general%ewn-1,model%general%ewn,2,ns-1,ns,ns+1)
          end do
       end if

       bwat(1:model%general%ewn,1:model%general%nsn) = model%tempwk%smth(1:model%general%ewn,1:model%general%nsn)

    case(1)
       ! apply periodic BC
       if (model%options%periodic_ew.eq.1) then
          write(*,*) 'Warning, periodic BC are not implement for this case yet'
       end if
       ! ** add any melt_water

       bwat = max(0.0d0,bwat + model%numerics%dttem * bmlt)

       model%tempwk%wphi = 0.
       model%tempwk%bwatu = 0.
       model%tempwk%bwatv = 0.
       model%tempwk%fluxew = 0.
       model%tempwk%fluxns = 0.
       model%tempwk%bint = 0.


       ! ** split time evolution into steps to avoid CFL problems

       do t_wat = 1,model%tempwk%nwat

          ! ** find potential surface using paterson p112, eq 4
          ! ** if no ice then set to sea level or land surface potential
          ! ** if frozen then set high 

          do ns = 1,model%general%nsn
             do ew = 1,model%general%ewn
                if (model%numerics%thklim < thck(ew,ns) .and. .not. floater(ew,ns)) then
                   call calcpmptb(pmpt,thck(ew,ns))
                   if (btem(ew,ns) == pmpt) then
                      model%tempwk%wphi(ew,ns) = model%tempwk%c(1) * (topg(ew,ns) + bwat(ew,ns)) + model%tempwk%c(2) * thck(ew,ns)
                   else
                      model%tempwk%wphi(ew,ns) = model%tempwk%c(1) * (topg(ew,ns) + thck(ew,ns))
                   end if
                else 
                   model%tempwk%wphi(ew,ns) = max(model%tempwk%c(1) * topg(ew,ns),0.0d0)
                end if
             end do
          end do

          ! ** determine x,y components of water velocity assuming
          ! ** contstant velocity magnitude and using potential
          ! ** to determine direction

          do ns = 2,model%general%nsn-1
             do ew = 2,model%general%ewn-1
                if (thck(ew,ns) > model%numerics%thklim) then

                   dwphidew = (model%tempwk%wphi(ew+1,ns) - model%tempwk%wphi(ew-1,ns)) / model%tempwk%c(3)       
                   dwphidns = (model%tempwk%wphi(ew,ns+1) - model%tempwk%wphi(ew,ns-1)) / model%tempwk%c(4)  

                   dwphi = - model%tempwk%watvel / sqrt(dwphidew**2 + dwphidns**2)

                   model%tempwk%bwatu(ew,ns) = dwphi * dwphidew  
                   model%tempwk%bwatv(ew,ns) = dwphi * dwphidns  

                else
                   model%tempwk%bwatu(ew,ns) = 0.0d0
                   model%tempwk%bwatv(ew,ns) = 0.0d0
                end if
             end do
          end do

          ! ** use two-step law wendroff to solve dW/dt = -dF/dx - dF/dy

          ! ** 1. find fluxes F=uW

          model%tempwk%fluxew = bwat * model%tempwk%bwatu
          model%tempwk%fluxns = bwat * model%tempwk%bwatv

          ! ** 2. do 1st LW step on staggered grid for dt/2

          do ns = 1,model%general%nsn-1
             do ew = 1,model%general%ewn-1

                bave = 0.25 * sum(bwat(ew:ew+1,ns:ns+1))

                if (bave > 0.0d0) then

                   model%tempwk%bint(ew,ns) = bave - &
                        model%tempwk%c(5) * (sum(model%tempwk%fluxew(ew+1,ns:ns+1)) - sum(model%tempwk%fluxew(ew,ns:ns+1))) - &
                        model%tempwk%c(6) * (sum(model%tempwk%fluxns(ew:ew+1,ns+1)) - sum(model%tempwk%fluxns(ew:ew+1,ns)))

                else
                   model%tempwk%bint(ew,ns) = 0.0d0
                end if
             end do
          end do

          ! ** 3. find fluxes F=uW on staggered grid griven new Ws

          model%tempwk%fluxew(1:model%general%ewn-1,1:model%general%nsn-1) = model%tempwk%bint * 0.25 * &
               (model%tempwk%bwatu(1:model%general%ewn-1,1:model%general%nsn-1) + &
               model%tempwk%bwatu(2:model%general%ewn,1:model%general%nsn-1) + &
               model%tempwk%bwatu(1:model%general%ewn-1,2:model%general%nsn) + &
               model%tempwk%bwatu(2:model%general%ewn,2:model%general%nsn))
          model%tempwk%fluxns(1:model%general%ewn-1,1:model%general%nsn-1) = model%tempwk%bint * 0.25 * &
               (model%tempwk%bwatv(1:model%general%ewn-1,1:model%general%nsn-1) + &
               model%tempwk%bwatv(2:model%general%ewn,1:model%general%nsn-1) + &
               model%tempwk%bwatv(1:model%general%ewn-1,2:model%general%nsn) + &
               model%tempwk%bwatv(2:model%general%ewn,2:model%general%nsn))

          ! ** 4. finally do 2nd LW step to get back on to main grid

          do ns = 2,model%general%nsn-1
             do ew = 2,model%general%ewn-1
                if (bwat(ew,ns) > 0.0d0) then

                   bwat(ew,ns) = bwat(ew,ns) - &
                        model%tempwk%c(7) * (sum(model%tempwk%fluxew(ew,ns-1:ns)) - sum(model%tempwk%fluxew(ew-1,ns-1:ns))) - &
                        model%tempwk%c(8) * (sum(model%tempwk%fluxns(ew-1:ew,ns)) - sum(model%tempwk%fluxns(ew-1:ew,ns-1)))

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

    ! How to call the flow router.
    ! call advectflow(bwat,phi,bmlt,model%geometry%mask)

    ! now also calculate basal water in velocity coord system
    call stagvarb(model%temper%bwat, &
         model%temper%stagbwat ,&
         model%general%  ewn, &
         model%general%  nsn)

  contains
    subroutine smooth_bwat(ewm,ew,ewp,nsm,ns,nsp)
      ! smoothing basal water distrib
      implicit none
      integer, intent(in) :: ewm,ew,ewp,nsm,ns,nsp
      if (blim(2) < bwat(ew,ns)) then
         model%tempwk%smth(ew,ns) = bwat(ew,ns) + model%paramets%bwat_smooth * &
              (bwat(ewm,ns) + bwat(ewp,ns) + bwat(ew,nsm) + bwat(ew,nsp) - 4.0d0 * bwat(ew,ns))
      else 
         model%tempwk%smth(ew,ns) = bwat(ew,ns)
      end if   
    end subroutine smooth_bwat
  end subroutine calcbwat
  
  subroutine find_dt_wat(dttem,estimate,dt_wat,nwat)
    
    implicit none
    
    real(dp), intent(out) :: dt_wat
    integer, intent(out) :: nwat
    real(dp), intent(in) :: dttem, estimate
    
    nwat = int(dttem/estimate) + 1
    dt_wat = dttem / nwat

  end subroutine find_dt_wat


  !-------------------------------------------------------------------

  subroutine corrpmpt(temp,thck,bwat,sigma,upn)

    use glimmer_global, only : dp

    implicit none 

    real(dp), dimension(:), intent(inout) :: temp
    real(dp), intent(in) :: thck, bwat
    integer,intent(in) :: upn
    real(dp),dimension(:),intent(in) :: sigma

    real(dp), dimension(:) :: pmptemp(size(temp))

    ! corrects a temperature column for melting point effects
    ! 1. if temperature at any point in column is above pressure melting point then 
    ! set temperature to pressure melting point 
    ! 2. if bed is wet set basal temperature to pressure melting point 

    call calcpmpt(pmptemp,thck,sigma)

    temp = dmin1(temp,pmptemp)

    if (bwat > 0.0d0) temp(upn) = pmptemp(upn)

  end subroutine corrpmpt

  !-------------------------------------------------------------------

  subroutine calcpmpt(pmptemp,thck,sigma)

    !*FD Returns the pressure melting point of water (degC)

    use glimmer_global, only : dp !, upn
    use physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    implicit none 

    real(dp), dimension(:), intent(out) :: pmptemp
    real(dp), intent(in) :: thck
    real(dp),intent(in),dimension(:) :: sigma

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp = fact * thck * sigma

  end subroutine calcpmpt

  !-------------------------------------------------------------------

  subroutine calcpmptb(pmptemp,thck)

    use glimmer_global, only : dp
    use physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    implicit none 

    real(dp), intent(out) :: pmptemp
    real(dp), intent(in) :: thck

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp = fact * thck 

  end subroutine calcpmptb

  !-------------------------------------------------------------------

  subroutine swchpnt(a,b,c,d,e)

    implicit none 

    integer, intent(inout) :: a, b, c 
    integer, intent(in) :: d, e

    if (a == d) then
       a = e
       b = d
       c = -1
    else
       a = d
       b = e
       c = 1
    end if

  end subroutine swchpnt

  !-------------------------------------------------------------------

  subroutine swapbndt(bc,a,b,c,d)

    use glimmer_global, only : dp

    implicit none

    real(dp), intent(out), dimension(:,:) :: a, c
    real(dp), intent(in), dimension(:,:) :: b, d
    integer, intent(in) :: bc

    if (bc == 0) then
       a = b
       c = d
    end if

  end subroutine swapbndt

end module glide_temp

