! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_thck.f90 - part of the GLIMMER ice model         + 
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

module glide_thck

  private
  public :: thck_nonlin_evolve

contains

  subroutine thck_nonlin_evolve(params,thck,usrf,lsrf,topg,acab,bmlt,btrc,flwa, &
       uvel,vvel,diffu,ubas,vbas,eus,newtemps,linear,dt)

    !*FD this subroutine solves the ice thickness equation by doing an outer, 
    !*FD non-linear iteration to update the diffusivities and in inner, linear
    !*FD iteration to calculate the new ice thickness distrib

    use glimmer_global, only : dp,sp
    use glide_thckADI, only: thckADI_type
    use glimmer_utils, only: stagvarb
    use glimmer_deriv, only: df_field_2d_staggered 
    use glide_thckCommon, only: velo_calc_velo, velo_integrate_flwa, velo_calc_diffu, glide_calclsrf
    use physcon,        only: rhoi, grav

    implicit none
    ! subroutine arguments
    type(thckADI_type),     intent(inout) :: params
    real(dp),dimension(:,:),intent(inout) :: thck     !< ice thickness
    real(dp),dimension(:,:),intent(inout) :: usrf     !< Surface elevation
    real(dp),dimension(:,:),intent(inout) :: lsrf     !< Lower elevation
    real(dp),dimension(:,:),intent(in)    :: topg     !< Bedrock topography
    real(sp),dimension(:,:),intent(in)    :: acab     !< surface mass balance
    real(dp),dimension(:,:),intent(in)    :: bmlt     !< Basal melt rate
    real(dp),dimension(:,:),intent(in)    :: btrc     !< basal traction coefficient 
    real(dp),dimension(:,:,:),intent(in)  :: flwa     !< Glen's A
    real(dp),dimension(:,:,:),intent(out) :: uvel     !< x-velocity
    real(dp),dimension(:,:,:),intent(out) :: vvel     !< y-velocity
    real(dp),dimension(:,:),intent(out)   :: diffu    !< diffusivity
    real(dp),dimension(:,:),intent(out)   :: ubas     !< basal x velocity
    real(dp),dimension(:,:),intent(out)   :: vbas     !< basal y velocity
    real(sp),               intent(in)    :: eus      !< Sea level
    logical,                intent(in)    :: newtemps !< true when we should recalculate Glen's A
    logical,                intent(in)    :: linear   !< true if we do a linear thickness calc
    real(dp),               intent(in)    :: dt       !< Timestep

    ! local variables
    integer, parameter :: pmax=50                       !*FD maximum Picard iterations
    real(kind=dp), parameter :: tol=1.0d-6
    real(kind=dp) :: residual
    integer p
    logical first_p
    logical empty
    real(dp),dimension(:),allocatable :: rhs
    integer, dimension(size(thck,1),size(thck,2)) :: mask
    real(dp),dimension(size(thck,1)-1,size(thck,2)-1)  :: fslip
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: stagthck
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: dusrfdew
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: dusrfdns
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: dthckdew
    real(dp), dimension(params%ewn-1,params%nsn-1)  :: dthckdns
    real(dp), dimension(params%ewn,  params%nsn  )  :: oldthck
    real(dp), dimension(params%ewn,  params%nsn  )  :: oldthck2
    integer :: totpts
    real(dp), parameter :: rhograv = - rhoi * grav
    real(dp) :: fc2_1, fc2_5

    ! Timestep-dependent constants
    fc2_1 = params%alpha * dt / (2.0d0 * params%dew * params%dew)
    fc2_5 = params%alpha * dt / (2.0d0 * params%dns * params%dns)
    
    ! Calculate staggered thickness
    call stagvarb(thck,stagthck,params%ewn,params%nsn)

    ! Calculate spatial derivatives
    call df_field_2d_staggered(usrf,params%dew,params%dns,dusrfdew,dusrfdns,.false.,.false.)
    call df_field_2d_staggered(thck,params%dew,params%dns,dthckdew,dthckdns,.false.,.false.)

    ! Masking for matrix solution
    call glide_maskthck(thck,acab,mask,totpts,empty)
    allocate(rhs(totpts))

    if (empty) then

       thck = dmax1(0.0d0,thck + acab * dt)
#ifdef DEBUG
       print *, "* thck empty - net accumulation added"
#endif
    else

       ! calculate basal velos
       if (newtemps) then

          fslip = rhograv * btrc

         ! calculate Glen's A if necessary
          call velo_integrate_flwa(     &
               params%velo_dups,        &
               params%depth,            &
               params%dintflwa,         &
               stagthck,                &
               flwa)
       end if

       first_p = .true.
       oldthck = thck
       ! do Picard iteration
       do p=1,pmax
          if (.not.linear) then
          oldthck2 = thck

          call stagvarb(thck,stagthck,params%ewn,params%nsn)
          call df_field_2d_staggered(usrf,params%dew,params%dns,dusrfdew,dusrfdns,.false., .false.)
          call df_field_2d_staggered(thck,params%dew,params%dns,dthckdew,dthckdns,.false., .false.)
          end if

          where (params%thklim < stagthck)
             ubas = fslip * stagthck**2  
          elsewhere
             ubas = 0.0d0
          end where

          ! calculate diffusivity
          call velo_calc_diffu(         &
               params%dintflwa,         &
               stagthck,                &
               dusrfdew,                &
               dusrfdns,                &
               diffu)

          ! get new thicknesses
          call thck_evolve(           &
               params,                &
               rhs,                   &
               mask,                  &
               diffu,                 &
               ubas,                  &
               lsrf,                  &
               acab,                  &
               bmlt,                  &
               thck,                  &
               topg,                  &
               usrf,                  &
               eus,                   &
               totpts,                &
               first_p,               &
               oldthck,               &
               thck,                  &
               params%periodic_ew,    &
               params%basal_mbal_flag, &
               fc2_1,                 &
               fc2_5,                 &
               dt)

          first_p = .false.
          residual = maxval(abs(thck-oldthck2))
          if ((residual.le.tol).or.linear) then
             exit
          end if
          
       end do

       ! calculate horizontal velocity field

      where (params%thklim < stagthck)
        vbas = ubas *  dusrfdns / stagthck
        ubas = ubas *  dusrfdew / stagthck
      elsewhere
        ubas = 0.0d0
        vbas = 0.0d0
      end where


       call velo_calc_velo(          &
            params%dintflwa,         &
            params%depth,            &
            stagthck,                &
            dusrfdew,                &
            dusrfdns,                &
            flwa,                    &
            diffu,                   &
            ubas,                    &
            vbas,                    &
            uvel,                    &
            vvel,                    &
            params%uflx,             &
            params%vflx)
    end if
 
    deallocate(rhs)

  end subroutine thck_nonlin_evolve

!---------------------------------------------------------------------------------

  subroutine thck_evolve(params,rhsd,mask,diffu,ubas,lsrf,acab,bmlt,thck,topg,usrf,eus,&
       totpts,calc_rhs,old_thck,new_thck,periodic_ew,basal_mbal,fc2_1,fc2_5,dt)

    !*FD set up sparse matrix and solve matrix equation to find new ice thickness distribution
    !*FD this routine does not override the old thickness distribution

    use glide_thckCommon, only: glide_calclsrf
    use glide_thckADI, only: thckADI_type
    use glimmer_global, only : dp, sp
    use glimmer_slap, only: slapMatrix_type, slapMatrix_init, slapMatrix_insertElement, &
         slapSolve, slapMatrix_dealloc

    implicit none

    ! subroutine arguments -------------------------------------------------------------

    type(thckADI_type)      :: params
    
    logical,                   intent(in)    :: calc_rhs    !< set to true when rhs should be calculated 
                                                            !! i.e. when doing lin solution or first picard iteration
    integer, dimension(:,:),   intent(in)    :: mask        !< Index mask for matrix mapping
    real(dp),dimension(:,:),   intent(in)    :: diffu       !< Diffusivity
    real(dp),dimension(:,:),   intent(in)    :: ubas        !< Basal x-velocity
    real(dp),dimension(:,:),   intent(inout) :: lsrf        !< Lower surface elevation
    real(sp),dimension(:,:),   intent(in)    :: acab        !< Surface mass balance
    real(dp),dimension(:,:),   intent(in)    :: bmlt        !< Basal melt rate
    real(dp),dimension(:,:),   intent(in)    :: thck        !< Ice thickness
    real(dp),dimension(:,:),   intent(in)    :: topg        !< Basal topography
    real(dp),dimension(:,:),   intent(inout) :: usrf        !< Upper surface elevation
    real(sp),                  intent(in)    :: eus         !< Eustatic sea level
    integer,                   intent(in)    :: totpts      !< Number of non-zero points in mask
    real(dp),dimension(:,:),   intent(in)    :: old_thck    !< contains ice thicknesses from previous time step
    real(dp),dimension(:,:),   intent(inout) :: new_thck    !< on entry contains first guess for new ice thicknesses
                                                            !< on exit contains ice thicknesses of new time step
    real(dp),dimension(totpts),intent(inout) :: rhsd
    integer,                   intent(in)    :: periodic_ew !< Periodic boundary conditions flag
    integer,                   intent(in)    :: basal_mbal  !< Include basal melt in mass-balance
    real(dp),                  intent(in)    :: fc2_1       !< A timestep-dependent constant
    real(dp),                  intent(in)    :: fc2_5       !< A timestep-dependent constant
    real(dp),                  intent(in)    :: dt          !< Timestep
    
    ! local variables ------------------------------------------------------------------

    real(dp),dimension(5)      :: sumd 
    real(dp),dimension(totpts) :: answ
    real(dp) :: err
    integer  :: linit
    integer  :: ew,ns

    type(slapMatrix_type) :: matrix

    ! Initialise sparse matrix object
    call slapMatrix_init(matrix,totpts,params%ewn*params%nsn*5)

    answ = 0.0

    ! Boundary Conditions ---------------------------------------------------------------
    ! lower and upper BC
    do ew = 1,params%ewn
       ns=1
       if (mask(ew,ns) /= 0) then
          call slapMatrix_insertElement(matrix,1.0d0,mask(ew,ns),mask(ew,ns))
          if (calc_rhs) then
             rhsd(mask(ew,ns)) = old_thck(ew,ns) 
          end if
          answ(mask(ew,ns)) = new_thck(ew,ns)
       end if
       ns=params%nsn
       if (mask(ew,ns) /= 0) then
          call slapMatrix_insertElement(matrix,1.0d0,mask(ew,ns),mask(ew,ns))
          if (calc_rhs) then
             rhsd(mask(ew,ns)) = old_thck(ew,ns) 
          end if
          answ(mask(ew,ns)) = new_thck(ew,ns)
       end if
    end do

    !left and right BC
    if (periodic_ew.eq.1) then
       do ns=2,params%nsn-1
          ew = 1
          if (mask(ew,ns) /= 0) then
             call findsums(diffu, &
                  ubas, &
                  sumd, &
                  fc2_1, &
                  fc2_5, &
                  params%ewn-2,params%ewn-1,ns-1,ns)
             call generate_row(matrix,sumd,mask, &
                  lsrf, &
                  acab, &
                  bmlt, &
                  old_thck, &
                  new_thck, &
                  rhsd, &
                  answ, &
                  calc_rhs, &
                  dt, &
                  params%fc2_3, &
                  params%fc2_4, &
                  params%ewn-2,ew,ew+1,ns-1,ns,ns+1, &
                  basal_mbal)
          end if
          ew=params%ewn
          if (mask(ew,ns) /= 0) then
             call findsums(diffu, &
                  ubas, &
                  sumd, &
                  fc2_1, &
                  fc2_5, &
                  1,2,ns-1,ns)
             call generate_row(matrix,sumd,mask, &
                  lsrf, &
                  acab, &
                  bmlt, &
                  old_thck, &
                  new_thck, &
                  rhsd, &
                  answ, &
                  calc_rhs, &
                  dt, &
                  params%fc2_3, &
                  params%fc2_4, &
                  ew-1,ew,3,ns-1,ns,ns+1, &
                  basal_mbal)
          end if
       end do
    else
       do ns=2,params%nsn-1
          ew=1
          if (mask(ew,ns) /= 0) then
             call slapMatrix_insertElement(matrix,1.0d0,mask(ew,ns),mask(ew,ns))
             if (calc_rhs) then
                rhsd(mask(ew,ns)) = old_thck(ew,ns) 
             end if
             answ(mask(ew,ns)) = new_thck(ew,ns)
          end if
          ew=params%ewn
          if (mask(ew,ns) /= 0) then
             call slapMatrix_insertElement(matrix,1.0d0,mask(ew,ns),mask(ew,ns))
             if (calc_rhs) then
                rhsd(mask(ew,ns)) = old_thck(ew,ns) 
             end if
             answ(mask(ew,ns)) = new_thck(ew,ns)
          end if
       end do
    end if

    ! ice body -------------------------------------------------------------------------

    do ns = 2,params%nsn-1
       do ew = 2,params%ewn-1

          if (mask(ew,ns) /= 0) then
                
             call findsums(diffu, &
                  ubas, &
                  sumd, &
                  fc2_1, &
                  fc2_5, &
                  ew-1,ew,ns-1,ns)
             call generate_row(matrix,sumd,mask, &
                  lsrf, &
                  acab, &
                  bmlt, &
                  old_thck, &
                  new_thck, &
                  rhsd, &
                  answ, &
                  calc_rhs, &
                  dt, &
                  params%fc2_3, &
                  params%fc2_4, &
                  ew-1,ew,ew+1,ns-1,ns,ns+1, &
                  basal_mbal)

          end if
       end do
    end do

    ! Solve the system using SLAP
    call slapSolve(matrix,rhsd,answ,linit,err)   

    ! Rejig the solution onto a 2D array
    do ns = 1,params%nsn
       do ew = 1,params%ewn 

          if (mask(ew,ns) /= 0) then
             new_thck(ew,ns) = answ(mask(ew,ns))
          end if

       end do
    end do

    new_thck = max(0.0d0, new_thck)

#ifdef DEBUG
    print *, "* thck ", linit, totpts,real(thk0*new_thck(params%ewn/2+1,params%nsn/2+1)), &
         real(vel0*maxval(abs(ubas)))
#endif

    ! calculate upper and lower surface
    call glide_calclsrf(thck, topg, eus, lsrf)
    usrf = max(0.d0,thck + lsrf)

    ! Deallocate matrix storage
    call slapMatrix_dealloc(matrix)

  end subroutine thck_evolve

!-------------------------------------------------------------------------

  subroutine generate_row(matrix,sumd,mask,lsrf,acab,bmlt,old_thck,new_thck,rhsd,answ, &
       calc_rhs,dt,fc2_3,fc2_4,ewm,ew,ewp,nsm,ns,nsp,basal_mbal)
    ! calculate row of sparse matrix equation
    use glimmer_global, only: dp,sp
    use glimmer_slap, only: slapMatrix_type, slapMatrix_insertElement
    implicit none
    type(slapMatrix_type),  intent(inout) :: matrix
    real(dp),dimension(5),  intent(in)    :: sumd
    integer, dimension(:,:),intent(in)    :: mask        !< Index mask for matrix mapping
    real(dp),dimension(:,:),intent(in)    :: lsrf
    real(sp),dimension(:,:),intent(in)    :: acab
    real(dp),dimension(:,:),intent(in)    :: bmlt
    real(dp),dimension(:,:),intent(in)    :: old_thck !< contains ice thicknesses from previous time step
    real(dp),dimension(:,:),intent(inout) :: new_thck !< on entry contains first guess for new ice thicknesses
    real(dp),dimension(:),  intent(inout) :: rhsd
    real(dp),dimension(:),  intent(out)   :: answ
    logical,                intent(in)    :: calc_rhs
    real(dp),               intent(in)    :: dt
    real(dp),               intent(in)    :: fc2_3
    real(dp),               intent(in)    :: fc2_4
    integer, intent(in) :: ewm,ew,ewp  ! ew index to left, central, right node
    integer, intent(in) :: nsm,ns,nsp  ! ns index to lower, central, upper node
    integer, intent(in) :: basal_mbal  !< Set ==1 if basal mass balance is considered

    ! fill sparse matrix
    call slapMatrix_insertElement(matrix,sumd(1),mask(ewm,ns),mask(ew,ns))       ! point (ew-1,ns)
    call slapMatrix_insertElement(matrix,sumd(2),mask(ewp,ns),mask(ew,ns))       ! point (ew+1,ns)
    call slapMatrix_insertElement(matrix,sumd(3),mask(ew,nsm),mask(ew,ns))       ! point (ew,ns-1)
    call slapMatrix_insertElement(matrix,sumd(4),mask(ew,nsp),mask(ew,ns))       ! point (ew,ns+1)
    call slapMatrix_insertElement(matrix,1.0d0 + sumd(5),mask(ew,ns),mask(ew,ns))! point (ew,ns)

    ! calculate RHS
    if (calc_rhs) then
       rhsd(mask(ew,ns)) =                    &
            old_thck(ew,ns) * (1.0d0 - fc2_3 * sumd(5))     &
            - fc2_3 * (old_thck(ewm,ns) * sumd(1)             &
                     + old_thck(ewp,ns) * sumd(2)             &
                     + old_thck(ew,nsm) * sumd(3)             &
                     + old_thck(ew,nsp) * sumd(4))            &
            - fc2_4 * (lsrf(ew,ns)  * sumd(5)  &
                     + lsrf(ewm,ns) * sumd(1)  &
                     + lsrf(ewp,ns) * sumd(2)  &
                     + lsrf(ew,nsm) * sumd(3)  &
                     + lsrf(ew,nsp) * sumd(4)) &
            + acab(ew,ns) * dt
       if(basal_mbal==1) then
          rhsd(mask(ew,ns)) = rhsd(mask(ew,ns)) - bmlt(ew,ns) * dt ! basal melt is +ve for mass loss
       end if
    end if

    answ(mask(ew,ns)) = new_thck(ew,ns)

  end subroutine generate_row

!-------------------------------------------------------------------------

  subroutine findsums(diffu,ubas,sumd,fc2_1,fc2_5,ewm,ew,nsm,ns)
  
    use glimmer_global, only: dp,sp

    implicit none

    real(dp),dimension(:,:),intent(in)  :: diffu
    real(dp),dimension(:,:),intent(in)  :: ubas
    real(dp),dimension(5),  intent(out) :: sumd 
    real(dp),               intent(in)  :: fc2_1
    real(dp),               intent(in)  :: fc2_5
    integer,                intent(in)  :: ewm,ew  ! ew index to left, right
    integer,                intent(in)  :: nsm,ns  ! ns index to lower, upper

    ! calculate sparse matrix elements
    sumd(1) = fc2_1 * ((diffu(ewm,nsm) + diffu(ewm,ns)) + (ubas (ewm,nsm) + ubas (ewm,ns)))
    sumd(2) = fc2_1 * ((diffu(ew,nsm)  + diffu(ew,ns))  + (ubas (ew,nsm)  + ubas (ew,ns)))
    sumd(3) = fc2_5 * ((diffu(ewm,nsm) + diffu(ew,nsm)) + (ubas (ewm,nsm) + ubas (ew,nsm)))
    sumd(4) = fc2_5 * ((diffu(ewm,ns)  + diffu(ew,ns))  + (ubas (ewm,ns)  + ubas (ew,ns)))
    sumd(5) = - (sumd(1) + sumd(2) + sumd(3) + sumd(4))

  end subroutine findsums

!-------------------------------------------------------------------------

  subroutine glide_maskthck(crita,critb,pointno,totpts,empty)
    
    !*FD Calculates the contents of the mask array.

    use glimmer_global, only : dp, sp 

    implicit none

    !-------------------------------------------------------------------------
    ! Subroutine arguments
    !-------------------------------------------------------------------------

    real(dp),dimension(:,:),intent(in)  :: crita      !*FD Ice thickness
    real(sp),dimension(:,:),intent(in)  :: critb      !*FD Mass balance
    integer, dimension(:,:),intent(out) :: pointno    !*FD Output mask
    integer,                intent(out) :: totpts     !*FD Total number of points
    logical,                intent(out) :: empty      !*FD Set if no mask points set.

    !-------------------------------------------------------------------------
    ! Internal variables
    !-------------------------------------------------------------------------

    integer :: covtot 
    integer :: ew,ns,ewn,nsn

    !-------------------------------------------------------------------------

    ewn=size(crita,1) ; nsn=size(crita,2)

    pointno = 0
    covtot  = 0 

    !-------------------------------------------------------------------------

    empty = .true.

    do ns = 1,nsn
      do ew = 1,ewn
        if ( thckcrit(crita(max(1,ew-1):min(ewn,ew+1),max(1,ns-1):min(nsn,ns+1)),critb(ew,ns)) ) then

          covtot = covtot + 1
          pointno(ew,ns) = covtot 
          if (empty) empty  = .false.

        end if
      end do
    end do
  
    totpts = covtot
 
  end subroutine glide_maskthck
  
  logical function thckcrit(ca,cb)

    use glimmer_global, only: sp,dp

    implicit none

    real(dp),dimension(:,:),intent(in) :: ca 
    real(sp),               intent(in) :: cb

    ! If the thickness in the region under consideration
    ! or the mass balance is positive, thckcrit is .true.

    if ( any((ca(:,:) > 0.0d0)) .or. cb > 0.0 ) then
       thckcrit = .true.
    else
       thckcrit = .false.
    end if

  end function thckcrit


end module glide_thck

