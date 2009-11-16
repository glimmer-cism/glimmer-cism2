! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_tempFullSoln.f90 - part of Glimmer-CISM            + 
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

module glide_tempFullSoln

  use glimmer_vertcoord, only: vertCoord_type
  use glimmer_coordinates, only : coordsystem_type
  use glimmer_global,    only: dp

  implicit none

  !> Holds parameters for the full temperature solution. Each grid
  !! which the model is to be used on requires an instance of this
  !! variable.
  type type_tempFullSoln
     private
     type(coordsystem_type) :: hCoord    !< the horizontal coordinate system
     type(vertCoord_type) :: zCoord      !< Vertical coordinate information
     real(dp)             :: thklim      !< Thickness threshold for calculation
     integer              :: periodic_ew !< Set to indicate periodic BCs
     integer              :: niter = 0   !< Maximum number of iterations
  end type type_tempFullSoln

  private
  public :: type_tempFullSoln, init_tempFullSoln, tstep_tempFullSoln, get_niter

contains

  !> Initialises parameters for full temperature solution
  subroutine init_tempFullSoln(params,hCoord,zCoord,thklim,periodic_ew)

    use glimmer_vertcoord, only: initVertCoord,initVertCoord
    use glimmer_log,     only: write_log, GM_FATAL

    implicit none

    type(type_tempFullSoln),intent(out) :: params  !< Temperature model parameters
    type(coordsystem_type), intent(in)  :: hCoord  !< the horizontal coordinate system
    type(vertCoord_type)                :: zCoord  !< the sigma coordinate system
    real(dp),               intent(in)  :: thklim  !< Thickness threshold for calculation
    integer,                intent(in)  :: periodic_ew !< Set to indicate periodic BCs

    if (VERT_DIFF.eq.0.)   call write_log('Vertical diffusion is switched off')
    if (HORIZ_ADV.eq.0.)   call write_log('Horizontal advection is switched off')
    if (VERT_ADV.eq.0.)    call write_log('Vertical advection is switched off')
    if (STRAIN_HEAT.eq.0.) call write_log('Strain heating is switched off')

    params%hCoord = hCoord
    call initVertCoord(params%zCoord,zCoord)

    params%thklim = thklim

    if (periodic_ew==0.or.periodic_ew==1) then
       params%periodic_ew = periodic_ew
    else
       call write_log('Unsupported value of periodic_ew in init_tempFullSoln',GM_FATAL)
    end if

  end subroutine init_tempFullSoln

  !------------------------------------------------------------------------------------

  !> Calculates the ice temperature - full solution
  subroutine tstep_tempFullSoln(params,temp,artm,thck,usrf,thkmask, &
       topg,uvel,vvel,ubas,vbas,wvel,wgrd,flwa,bheatflx,bwat,bmlt,dt)

    use glimmer_utils,    only: hsum4,tridiag,stagvarb
    use glimmer_deriv,    only: df_field_2d_staggered
    use glimmer_global,   only: dp,sp
    use glimmer_paramets, only: thk0, acc0, tim0, len0, vis0, scyr
    use glimmer_pmpt,     only: corrpmpt
    use glimmer_mask,     only: is_float, is_thin
    use physcon,          only: rhoi, shci, coni, grav, gn, lhci, rhow

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(type_tempFullSoln),      intent(inout) :: params   !< temperature model parameters
    real(dp), dimension(:,0:,0:), intent(inout) :: temp     !< Ice temperature (upn,0:ewn+1,0:nsn+1)
    real(sp), dimension(:,:),     intent(in)    :: artm     !< Surface air temperature (ewn,nsn)
    real(dp), dimension(:,:),     intent(in)    :: thck     !< Ice thickness (ewn,nsn)
    real(dp), dimension(:,:),     intent(in)    :: usrf     !< Upper surface elevation (ewn,nsn)
    integer,  dimension(:,:),     intent(in)    :: thkmask  !< Mask for ice thickness (ewn,nsn)
    real(dp), dimension(:,:),     intent(in)    :: topg     !< basal topography (ewn,nsn)
    real(dp), dimension(:,:,:),   intent(in)    :: uvel     !< x-velocity (ewn-1,nsn-1)
    real(dp), dimension(:,:,:),   intent(in)    :: vvel     !< y-velocity (ewn-1,nsn-1)
    real(dp), dimension(:,:),     intent(in)    :: ubas     !< basal x-velocity (ewn-1,nsn-1)
    real(dp), dimension(:,:),     intent(in)    :: vbas     !< basal y-velocity (ewn-1,nsn-1)
    real(dp), dimension(:,:,:),   intent(in)    :: wvel     !< Vertical velocity (upn,ewn,nsn)
    real(dp), dimension(:,:,:),   intent(in)    :: wgrd     !< Vertical grid velocity (upn,ewn,nsn)
    real(dp), dimension(:,:,:),   intent(in)    :: flwa     !< Glen's A (upn,ewn,nsn)
    real(dp), dimension(:,:),     intent(in)    :: bheatflx !< Basal heat flux (ewn,nsn)
    real(dp), dimension(:,:),     intent(in)    :: bwat     !< Basal water depth
    real(dp), dimension(:,:),     intent(out)   :: bmlt     !< Basal melt rate
    real(dp),                     intent(in)    :: dt       !< Timestep (years)

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp),dimension(params%zCoord%upn) :: subd, diag, supd, rhsd
    real(dp),dimension(params%zCoord%upn) :: prevtemp, iteradvt, diagadvt
    real(dp) :: tempresid

    integer :: iter
    integer :: ew,ns

    real(dp),parameter :: tempthres = 0.001d0, floatlim = 10.0d0 / thk0
    integer, parameter :: mxit = 100
    integer, parameter :: ewbc = 1, nsbc = 1 
    real(dp),parameter :: wmax = 5.0d0 * tim0 / (scyr * thk0)

    real(dp), dimension(params%zCoord%upn) :: weff

    real(dp) :: advconst1,advconst2
    real(dp) :: cons1,cons2,cons3,cons4
    real(dp),dimension(params%zCoord%upn) :: c1
    real(dp) :: f3
    real(dp) :: slidef1, slidef2

    real(dp),dimension(params%zCoord%upn,size(thck,1),size(thck,2)) :: initadvt
    real(dp),dimension(params%zCoord%upn,size(thck,1),size(thck,2)) :: inittemp
    real(dp),dimension(params%zCoord%upn,size(thck,1),size(thck,2)) :: dissip
    real(dp),dimension(params%zCoord%upn,size(thck,1),size(thck,2)) :: hadv_u
    real(dp),dimension(params%zCoord%upn,size(thck,1),size(thck,2)) :: hadv_v
    real(dp),dimension(size(thck,1)-1,size(thck,2)-1) :: stagthck
    real(dp),dimension(size(thck,1)-1,size(thck,2)-1) :: dusrfdew
    real(dp),dimension(size(thck,1)-1,size(thck,2)-1) :: dusrfdns

    !------------------------------------------------------------------------------------
    ! Set up various parameters
    !------------------------------------------------------------------------------------

    advconst1 = HORIZ_ADV * dt / (16.0d0 * params%hCoord%delta(1))
    advconst2 = HORIZ_ADV * dt / (16.0d0 * params%hCoord%delta(2))

    cons1 = 2.0d0 * tim0 * dt * coni / (2.0d0 * rhoi * shci * thk0**2)
    cons2 = dt / 2.0d0
    cons3 = VERT_DIFF * 2.0d0 * tim0 * dt / (thk0 * rhoi * shci)
    cons4 = VERT_ADV  * tim0 * acc0 * dt / coni

    c1 = STRAIN_HEAT * (params%zCoord%sigma * rhoi * grav * thk0**2 / len0)**(gn+1) * &
         2.0d0 * vis0 * dt * tim0 / (16.0d0 * rhoi * shci)

    f3 = tim0 * thk0 * rhoi * shci /  (thk0 * tim0 * dt * lhci * rhoi)

    ! sliding contribution to basal heat flux

    slidef1 = VERT_DIFF * grav * thk0 * dt / shci                      ! vert diffusion
    slidef2 = VERT_ADV  * rhoi * grav * acc0 * thk0 * thk0 * dt / coni ! vert advection

    !------------------------------------------------------------------------------------
    ! ewbc/nsbc set the type of boundary condition applied at the end of
    ! the domain. a value of 0 implies zero gradient.
    !------------------------------------------------------------------------------------

    inittemp = 0.0d0
    initadvt = 0.0d0

    ! Calculate staggered variables and spatial derivatives -----------------------------

    call stagvarb(thck,stagthck,params%hCoord%size(1),params%hCoord%size(2))
    call df_field_2d_staggered(usrf,params%hCoord%delta(1),params%hCoord%delta(2),dusrfdew,dusrfdns,.false.,.false.)

    ! Calculate dissipative term --------------------------------------------------------

    call finddisp(dissip,              &
         thck,                         &
         stagthck,                     &
         dusrfdew,                     &
         dusrfdns,                     &
         flwa,                         &
         params%hCoord%size(1),                   &
         params%hCoord%size(2),                   &
         c1,                           &
         params%thklim)

    ! translate velo field --------------------------------------------------------------

    do ns = 2,params%hCoord%size(2)-1
       do ew = 2,params%hCoord%size(1)-1
          hadv_u(:,ew,ns) = advconst1 * hsum4(uvel(:,ew-1:ew,ns-1:ns))
          hadv_v(:,ew,ns) = advconst2 * hsum4(vvel(:,ew-1:ew,ns-1:ns))
       end do
    end do

    ! Calculate initial upwinding terms -------------------------------------------------

    call hadvall(initadvt,              &
         temp,                          &
         thck,                          &
         params%thklim,                 &
         hadv_u,                        &
         hadv_v,                        &
         params%hCoord%size(1),                    &
         params%hCoord%size(2),                    &
         params%zCoord%upn)

    ! Iterative temperature solution ----------------------------------------------------

    iter = 0
    tempresid = abs(tempthres*2.0)  ! To make sure the loop is executed at least once
   
    do while (tempresid.gt.tempthres .and. iter.le.mxit)
       tempresid = 0.0d0

       do ns = 2,params%hCoord%size(2)-1
          do ew = 2,params%hCoord%size(1)-1
             if(thck(ew,ns) > params%thklim) then

                ! Calculate effective vertical velocity
                weff = wvel(:,ew,ns) - wgrd(:,ew,ns)

                ! Set effective vertical velocity to zero if it exceeds a threshold
                if (maxval(abs(weff)) > wmax) then
                   weff = 0.0d0
                end if

                ! Calculate upwinded advection term
                call hadvpnt(iteradvt,                       &
                     diagadvt,                               &
                     temp(:,ew-2:ew+2,ns),                   &
                     temp(:,ew,ns-2:ns+2),                   &
                     hadv_u(:,ew,ns),                        &
                     hadv_v(:,ew,ns),                        &
                     params%zCoord%upn)

                call findvtri(params%zCoord,     &
                     thck(ew,ns),                &
                     subd,                       &
                     diag,                       &
                     supd,                       &
                     diagadvt,                   &
                     weff,                       &
                     is_float(thkmask(ew,ns)),   &
                     cons1,                      &
                     cons2)

                if (iter==0) then
                   call findvtri_init(initadvt,     &
                        dissip,                     &
                        inittemp,                   &
                        bheatflx,                   &
                        dusrfdew,                   &
                        dusrfdns,                   &
                        params%zCoord,              &
                        ew,                         &
                        ns,                         &
                        subd,                       &
                        diag,                       &
                        supd,                       &
                        weff,                       &
                        ubas,                       &
                        vbas,                       &
                        temp(:,ew,ns),              &
                        thck(ew,ns),                &
                        is_float(thkmask(ew,ns)),   &
                        cons3,                      &
                        cons4,                      &
                        slidef1,                    &
                        slidef2)
                end if

                call findvtri_rhs(params%zCoord,    &
                     inittemp(:,ew,ns),             &
                     artm(ew,ns),                   &
                     iteradvt,                      &
                     rhsd,                          &
                     is_float(thkmask(ew,ns)))

                prevtemp = temp(:,ew,ns)

                call tridiag(subd(1:params%zCoord%upn),            &
                     diag(1:params%zCoord%upn),                    &
                     supd(1:params%zCoord%upn),                    &
                     temp(1:params%zCoord%upn,ew,ns),              &
                     rhsd(1:params%zCoord%upn))

                call corrpmpt(temp(:,ew,ns),         &
                     thck(ew,ns),                    &
                     bwat(ew,ns),                    &
                     params%zCoord%sigma)

                tempresid = max(tempresid,maxval(abs(temp(:,ew,ns)-prevtemp(:))))
             endif
          end do
       end do

       iter = iter + 1
    end do

    params%niter = max(params%niter, iter )
       
    ! set temperature of thin ice to the air temperature and set ice free nodes to zero

    do ns = 1,params%hCoord%size(2)
       do ew = 1,params%hCoord%size(1)
          if (is_thin(thkmask(ew,ns))) then
             temp(:,ew,ns) = min(0.0d0,dble(artm(ew,ns)))
          else if (thkmask(ew,ns)<0) then
             temp(:,ew,ns) = 0.0d0
          end if
       end do
    end do

    ! apply periodic ew BC ----------------------------------------------------------

    if (params%periodic_ew.eq.1) then
       temp(:,0           ,:) = temp(:,params%hCoord%size(1)-2,:)
       temp(:,1           ,:) = temp(:,params%hCoord%size(1)-1,:)
       temp(:,params%hCoord%size(1)  ,:) = temp(:,2           ,:)
       temp(:,params%hCoord%size(1)+1,:) = temp(:,3           ,:)
    end if

    ! Calculate basal melt rate --------------------------------------------------

    call calcbmlt(dissip,                 &
         params%zCoord,                   &
         temp,                            &
         thck,                            &
         stagthck,                        &
         dusrfdew,                        &
         dusrfdns,                        &
         ubas,                            &
         vbas,                            &
         bheatflx,                        &
         bmlt,                            &
         is_float(thkmask),               &
         params%hCoord%size(2),                      &
         params%hCoord%size(1),                      &
         params%thklim,                   &
         params%periodic_ew,              &
         f3)

    ! Output some information ----------------------------------------
#ifdef DEBUG
    print *, "* temp ",  iter, params%niter, &
         real(temp(params%zCoord%upn,params%hCoord%size(1)/2+1,params%hCoord%size(2)/2+1))
#endif

  end subroutine tstep_tempFullSoln

  !-------------------------------------------------------------------

  integer function get_niter(params)

    implicit none

    type(type_tempFullSoln), intent(in) :: params

    get_niter = params%niter

  end function get_niter

  !-------------------------------------------------------------------

  subroutine hadvpnt(iteradvt,diagadvt,tempx,tempy,u,v,upn)

    use glimmer_global, only: dp

    implicit none

    integer,                    intent(in)  :: upn        !< Number of points in vertical
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

    real(dp),                   intent(in)  :: advconst1
    real(dp),                   intent(in)  :: advconst2
    real(dp), dimension(:),     intent(out) :: iteradvt
    real(dp), dimension(:),     intent(out) :: diagadvt
    real(dp), dimension(:,:),   intent(in)  :: tempx
    real(dp), dimension(:,:),   intent(in)  :: tempy
    real(dp), dimension(:,:,:), intent(in)  :: uvel
    real(dp), dimension(:,:,:), intent(in)  :: vvel

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

    use glimmer_global, only: dp 

    implicit none

    integer,                                 intent(in)  :: ewn
    integer,                                 intent(in)  :: nsn
    integer,                                 intent(in)  :: upn
    real(dp), dimension(upn,ewn,nsn),        intent(out) :: initadvt
    real(dp), dimension(upn,0:ewn+1,0:nsn+1),intent(in)  :: temp
    real(dp), dimension(ewn,nsn),            intent(in)  :: thck
    real(dp),                                intent(in)  :: thklim
    real(dp), dimension(upn,ewn,nsn),        intent(in)  :: hadv_u
    real(dp), dimension(upn,ewn,nsn),        intent(in)  :: hadv_v

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

  subroutine findvtri(zCoord,thck,subd,diag,supd,diagadvt,weff,float,cons1,cons2)

    use glimmer_vertcoord, only: vertCoord_type
    use glimmer_global,    only: dp

    implicit none

    type(vertCoord_type),   intent(in)  :: zCoord
    real(dp),               intent(in)  :: thck
    real(dp), dimension(:), intent(in)  :: weff
    real(dp), dimension(:), intent(in)  :: diagadvt
    real(dp), dimension(:), intent(out) :: subd
    real(dp), dimension(:), intent(out) :: diag
    real(dp), dimension(:), intent(out) :: supd
    logical,                intent(in)  :: float
    real(dp),               intent(in)  :: cons1
    real(dp),               intent(in)  :: cons2

    real(dp) :: diff_factor,adv_factor

    diff_factor = VERT_DIFF * cons1 / thck**2
    adv_factor  = VERT_ADV  * cons2 / thck
    
    subd(2:zCoord%upn-1) =   adv_factor  * weff(2:zCoord%upn-1) * zCoord%dups(2:zCoord%upn-1,3)
    supd(2:zCoord%upn-1) = - subd(2:zCoord%upn-1) - diff_factor * zCoord%dups(2:zCoord%upn-1,2)
    subd(2:zCoord%upn-1) =   subd(2:zCoord%upn-1) - diff_factor * zCoord%dups(2:zCoord%upn-1,1)

    diag(2:zCoord%upn-1) = 1.0d0 - subd(2:zCoord%upn-1) - supd(2:zCoord%upn-1) + diagadvt(2:zCoord%upn-1)

    ! Upper surface: hold temperature constant

    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0

    ! now do the basal boundary
    ! for grounded ice, a heat flux is applied
    ! for floating ice, temperature held constant

    if (float) then
       supd(zCoord%upn) = 0.0d0
       subd(zCoord%upn) = 0.0d0
       diag(zCoord%upn) = 1.0d0
    else 
       supd(zCoord%upn) = 0.0d0 
       subd(zCoord%upn) = -0.5*diff_factor/(zCoord%dupn**2)
       diag(zCoord%upn) = 1.0d0 - subd(zCoord%upn) + diagadvt(zCoord%upn)
    end if

  end subroutine findvtri

  !-------------------------------------------------------------------------

  subroutine findvtri_init(initadvt,dissip,inittemp,bheatflx,dusrfdew,dusrfdns, &
       zCoord,ew,ns,subd,diag,supd,weff,ubas,vbas,temp,thck,float, &
       cons3,cons4,slidef1,slidef2)

    !*FD called during first iteration to set inittemp

    use glimmer_global, only: dp
    use glimmer_pmpt,   only: pmpt

    implicit none

    real(dp),dimension(:,:,:),intent(in)  :: initadvt
    real(dp),dimension(:,:,:),intent(in)  :: dissip
    real(dp),dimension(:,:,:),intent(out) :: inittemp
    real(dp),dimension(:,:),  intent(in)  :: bheatflx
    real(dp),dimension(:,:),  intent(in)  :: dusrfdew
    real(dp),dimension(:,:),  intent(in)  :: dusrfdns
    type(vertCoord_type),     intent(in)  :: zCoord
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
    real(dp),                 intent(in)  :: cons3
    real(dp),                 intent(in)  :: cons4
    real(dp),                 intent(in)  :: slidef1
    real(dp),                 intent(in)  :: slidef2

    ! local variables
    real(dp) :: slterm
    integer  :: ewp,nsp
    integer  :: slide_count

    ! Main body points
    inittemp(2:zCoord%upn-1,ew,ns) = temp(2:zCoord%upn-1) * (2.0d0 - diag(2:zCoord%upn-1)) &
         - temp(1:zCoord%upn-2) * subd(2:zCoord%upn-1)         &
         - temp(3:zCoord%upn)   * supd(2:zCoord%upn-1)         & 
         - initadvt(2:zCoord%upn-1,ew,ns)               &
         + dissip(2:zCoord%upn-1,ew,ns)
    
    ! Basal boundary points
    if (float) then
       inittemp(zCoord%upn,ew,ns) = pmpt(thck)
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

       inittemp(zCoord%upn,ew,ns) = temp(zCoord%upn) * (2.0d0 - diag(zCoord%upn)) &
            - temp(zCoord%upn-1) * subd(zCoord%upn)                        &
            - 0.5 * cons3 * bheatflx(ew,ns) / (thck * zCoord%dupn) &  ! geothermal heat flux (diff)
            - slidef1 * slterm / zCoord%dupn                 &        ! sliding heat flux    (diff)
            - cons4 * bheatflx(ew,ns) * weff(zCoord%upn)            &        ! geothermal heat flux (adv)
            - slidef2 * thck * slterm * weff(zCoord%upn)            &        ! sliding heat flux    (adv)
            - initadvt(zCoord%upn,ew,ns)                            &
            + dissip(zCoord%upn,ew,ns)

    end if

  end subroutine findvtri_init

  !-------------------------------------------------------------------------

  subroutine findvtri_rhs(zCoord,inittemp,artm,iteradvt,rhsd,float)

    !*FD RHS of temperature tri-diag system for a single column

    use glimmer_global, only: dp, sp 

    implicit none

    type(vertCoord_type), intent(in)  :: zCoord
    real(dp),dimension(:),intent(in)  :: inittemp
    real(sp),             intent(in)  :: artm 
    real(dp),dimension(:),intent(in)  :: iteradvt
    real(dp),dimension(:),intent(out) :: rhsd
    logical,              intent(in)  :: float

    ! upper boundary condition

    rhsd(1) = artm

    if (float) then
       rhsd(zCoord%upn) = inittemp(zCoord%upn)    
    else
       rhsd(zCoord%upn) = inittemp(zCoord%upn) - iteradvt(zCoord%upn)
    end if

    rhsd(2:zCoord%upn-1) = inittemp(2:zCoord%upn-1) - iteradvt(2:zCoord%upn-1)

  end subroutine findvtri_rhs

  !-----------------------------------------------------------------------

  subroutine finddisp(dissip,thck,stagthck,dusrfdew,dusrfdns,flwa,ewn,nsn,c1,thklim)

    use glimmer_global, only: dp
    use physcon,        only: gn

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

  subroutine calcbmlt(dissip,zCoord,temp,thck,stagthck,dusrfdew,dusrfdns,ubas,vbas,bheatflx, &
       bmlt,floater,nsn,ewn,thklim,periodic_ew,f3)

    use glimmer_global,   only: dp
    use glimmer_paramets, only: thk0, tim0, vel0, len0
    use physcon,          only: coni, lhci, rhoi, grav
    use glimmer_pmpt,     only: calcpmpt

    implicit none 

    real(dp), dimension(:,:,:),   intent(in)  :: dissip
    type(vertCoord_type),         intent(in)  :: zCoord
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
    integer,                      intent(in)  :: nsn
    integer,                      intent(in)  :: ewn
    real(dp),                     intent(in)  :: thklim
    integer,                      intent(in)  :: periodic_ew
    real(dp),                     intent(in)  :: f3
    

    real(dp), dimension(zCoord%upn) :: pmptemp
    real(dp) :: slterm, newmlt

    integer :: ewp, nsp,up,ew,ns

    real(dp),parameter :: f1 = tim0 * coni / (thk0**2 * lhci * rhoi)
    real(dp),parameter :: f2 = tim0 / (thk0 * lhci * rhoi)
    real(dp),parameter :: f4 = tim0 * thk0**2 * vel0 * grav * rhoi / (4.0d0 * thk0 * len0 * rhoi * lhci)

    do ns = 2, nsn-1
       do ew = 2, ewn-1
          if (thck(ew,ns) > thklim .and. .not. floater(ew,ns)) then

             call calcpmpt(pmptemp,thck(ew,ns),zCoord%sigma)

             if (abs(temp(zCoord%upn,ew,ns)-pmptemp(zCoord%upn)) .lt. 0.001) then

                slterm = 0.0d0

                do nsp = ns-1,ns
                   do ewp = ew-1,ew
                      slterm = slterm - stagthck(ewp,nsp) * &
                           (dusrfdew(ewp,nsp) * ubas(ewp,nsp) + dusrfdns(ewp,nsp) * vbas(ewp,nsp))
                   end do
                end do

                bmlt(ew,ns) = 0.0d0
                newmlt = f4 * slterm - f2 * bheatflx(ew,ns) &
                     + f3 * zCoord%dupc(zCoord%upn) * thck(ew,ns) * dissip(zCoord%upn,ew,ns)

                up = zCoord%upn - 1

                do while (abs(temp(up,ew,ns)-pmptemp(up)) .lt. 0.001 .and. up .ge. 3)
                   bmlt(ew,ns) = bmlt(ew,ns) + newmlt
                   newmlt = f3 * zCoord%dupc(up) * thck(ew,ns) * dissip(up,ew,ns)
                   up = up - 1
                end do

                up = up + 1

                if (up == zCoord%upn) then
                   bmlt(ew,ns) = newmlt - &
                        f1 * ( (temp(up-2,ew,ns) - pmptemp(up-2)) * zCoord%dupa(up) &
                        + (temp(up-1,ew,ns) - pmptemp(up-1)) * zCoord%dupb(up) ) / thck(ew,ns) 
                else
                   bmlt(ew,ns) = bmlt(ew,ns) + max(0.0d0, newmlt - &
                        f1 * ( (temp(up-2,ew,ns) - pmptemp(up-2)) * zCoord%dupa(up) &
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

end module glide_tempFullSoln
