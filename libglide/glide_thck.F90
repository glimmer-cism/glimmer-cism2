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

  use glide_types

  private
  public :: init_thck, thck_nonlin_evolve, thck_lin_evolve

#ifdef DEBUG_PICARD
  ! debugging Picard iteration
  integer, private, parameter :: picard_unit=101
  real, private, parameter    :: picard_interval=500.
  integer, private            :: picard_max=0
#endif

contains

  subroutine init_thck(model)
    !*FD initialise work data for ice thickness evolution
    use glimmer_log
    implicit none
    type(glide_global_type) :: model

    
    model%pcgdwk%fc2 = (/ model%numerics%alpha * model%numerics%dt / (2.0d0 * model%numerics%dew * model%numerics%dew), &
         model%numerics%dt, &
         (1.0d0-model%numerics%alpha) / model%numerics%alpha, &
         1.0d0 / model%numerics%alpha, model%numerics%alpha * model%numerics%dt / &
         (2.0d0 * model%numerics%dns * model%numerics%dns), &
         0.0d0 /) 

#ifdef DEBUG_PICARD
    call write_log('Logging Picard iterations')
    open(picard_unit,name='picard_info.data',status='unknown')
    write(picard_unit,*) '#time    max_iter'
#endif

  end subroutine init_thck

!---------------------------------------------------------------------------------

  subroutine thck_lin_evolve(model,newtemps)

    !*FD this subroutine solves the linearised ice thickness equation by computing the
    !*FD diffusivity from quantities of the previous time step

    use glide_velo
    use glide_thckCommon, only: velo_calc_velo, velo_integrate_flwa, velo_calc_diffu, glide_calclsrf    
    implicit none
    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A

    logical :: empty
    real(dp),dimension(:),allocatable :: rhs

#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_mask1)
#endif
    call glide_maskthck( &
         model%geometry% thck,      &
         model%climate%  acab,      &
         model%geometry% mask,      &
         model%geometry% totpts,    &
         empty)
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_mask1)
#endif

    allocate(rhs(model%geometry%totpts))

    if (empty) then

       model%geometry%thck = dmax1(0.0d0,model%geometry%thck + model%climate%acab * model%pcgdwk%fc2(2))
#ifdef DEBUG       
       print *, "* thck empty - net accumulation added", model%numerics%time
#endif
    else

       ! calculate basal velos
       if (newtemps) then
          call slipvelo(model%velowk,   &
               1,                       &
               model%numerics%thklim,   &
               model%geomderv%stagthck, &
               model%geomderv%dusrfdew, &
               model%geomderv%dusrfdns, &
               model%velocity% btrc,    &
               model%velocity% ubas,    &
               model%velocity% vbas)
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(     &
               model%velowk%dups,       &
               model%velowk%depth,      &
               model%velowk%dintflwa,   &
               model%geomderv%stagthck, &
               model%temper%flwa)
       end if
       call slipvelo(model%velowk,   &
            2,                       &
            model%numerics%thklim,   &
            model%geomderv%stagthck, &
            model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns, &
            model%velocity% btrc,    &
            model%velocity% ubas,    &
            model%velocity% vbas)

       ! calculate diffusivity
       call velo_calc_diffu(model%velowk%dintflwa,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%velocity%diffu)

       ! get new thicknesses
       call thck_evolve(model,rhs,model%geometry%mask,model%geometry%totpts, &
            .true.,model%geometry%thck,model%geometry%thck,model%options%periodic_ew)

       ! calculate horizontal velocity field
       call slipvelo(model%velowk,   &
            3,                       &
            model%numerics%thklim,   &
            model%geomderv%stagthck, &
            model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns, &
            model%velocity% btrc,    &
            model%velocity% ubas,    &
            model%velocity% vbas)
       call velo_calc_velo(model%velowk%dintflwa,model%velowk%depth,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%temper%flwa,model%velocity%diffu,model%velocity%ubas, &
            model%velocity%vbas,model%velocity%uvel,model%velocity%vvel,model%velocity%uflx,model%velocity%vflx)
    end if

    deallocate(rhs)

  end subroutine thck_lin_evolve

!---------------------------------------------------------------------------------

  subroutine thck_nonlin_evolve(model,newtemps,linear)

    !*FD this subroutine solves the ice thickness equation by doing an outer, 
    !*FD non-linear iteration to update the diffusivities and in inner, linear
    !*FD iteration to calculate the new ice thickness distrib

    use glimmer_global, only : dp
    use glide_velo
    use glide_setup
    use glimmer_utils, only: stagvarb
    use glimmer_deriv, only: df_field_2d_staggered 
    use glide_thckCommon, only: velo_calc_velo, velo_integrate_flwa, velo_calc_diffu, glide_calclsrf

    implicit none
    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A
    logical, intent(in) :: linear

    ! local variables
    integer, parameter :: pmax=50                       !*FD maximum Picard iterations
    real(kind=dp), parameter :: tol=1.0d-6
    real(kind=dp) :: residual
    integer p
    logical first_p
    logical empty
    real(dp),dimension(:),allocatable :: rhs

#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_mask1)
#endif
    call glide_maskthck( &
         model%geometry% thck,      &
         model%climate%  acab,      &
         model%geometry% mask,      &
         model%geometry% totpts,    &
         empty)
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_mask1)
#endif

    allocate(rhs(model%geometry%totpts))

    if (empty) then

       model%geometry%thck = dmax1(0.0d0,model%geometry%thck + model%climate%acab * model%pcgdwk%fc2(2))
#ifdef DEBUG
       print *, "* thck empty - net accumulation added", model%numerics%time
#endif
    else

       ! calculate basal velos
       if (newtemps) then
          call slipvelo(model%velowk,   &
               1,                       &
               model%numerics%thklim,   &
               model%geomderv%stagthck, &
               model%geomderv%dusrfdew, &
               model%geomderv%dusrfdns, &
               model%velocity% btrc,    &
               model%velocity% ubas,    &
               model%velocity% vbas)
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(     &
               model%velowk%dups,       &
               model%velowk%depth,      &
               model%velowk%dintflwa,   &
               model%geomderv%stagthck, &
               model%temper%flwa)
       end if

       first_p = .true.
       model%thckwk%oldthck = model%geometry%thck
       ! do Picard iteration
       do p=1,pmax
          model%thckwk%oldthck2 = model%geometry%thck

          call stagvarb(model%geometry% thck, &
               model%geomderv% stagthck,&
               model%general%  ewn, &
               model%general%  nsn)

          call df_field_2d_staggered(model%geometry%usrf, &
               model%numerics%dew, model%numerics%dns, &
               model%geomderv%dusrfdew, & 
               model%geomderv%dusrfdns, &
               .false., .false.)

          call df_field_2d_staggered(model%geometry%thck, &
               model%numerics%dew, model%numerics%dns, &
               model%geomderv%dthckdew, & 
               model%geomderv%dthckdns, &
               .false., .false.)

          call slipvelo(model%velowk,   &
               2,                       &
               model%numerics%thklim,   &
               model%geomderv%stagthck, &
               model%geomderv%dusrfdew, &
               model%geomderv%dusrfdns, &
               model%velocity% btrc,    &
               model%velocity% ubas,    &
               model%velocity% vbas)

          ! calculate diffusivity
          call velo_calc_diffu(model%velowk%dintflwa,model%geomderv%stagthck,model%geomderv%dusrfdew, &
               model%geomderv%dusrfdns,model%velocity%diffu)

          ! get new thicknesses
          call thck_evolve(model,rhs,model%geometry%mask,model%geometry%totpts, &
               first_p,model%thckwk%oldthck,model%geometry%thck,model%options%periodic_ew)

          first_p = .false.
          residual = maxval(abs(model%geometry%thck-model%thckwk%oldthck2))
          if ((residual.le.tol).or.linear) then
             exit
          end if
          
       end do
#ifdef DEBUG_PICARD
       picard_max=max(picard_max,p)
       if (model%numerics%tinc > mod(model%numerics%time,picard_interval)) then
          write(picard_unit,*) model%numerics%time,p
          picard_max = 0
       end if
#endif

       ! calculate horizontal velocity field
       call slipvelo(model%velowk,   &
            3,                       &
            model%numerics%thklim,   &
            model%geomderv%stagthck, &
            model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns, &
            model%velocity% btrc,    &
            model%velocity% ubas,    &
            model%velocity% vbas)
       call velo_calc_velo(model%velowk%dintflwa,model%velowk%depth,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%temper%flwa,model%velocity%diffu,model%velocity%ubas, &
            model%velocity%vbas,model%velocity%uvel,model%velocity%vvel,model%velocity%uflx,model%velocity%vflx)
    end if
 
    deallocate(rhs)

  end subroutine thck_nonlin_evolve

!---------------------------------------------------------------------------------

  subroutine thck_evolve(model,rhsd,mask,totpts,calc_rhs,old_thck,new_thck,periodic_ew)

    !*FD set up sparse matrix and solve matrix equation to find new ice thickness distribution
    !*FD this routine does not override the old thickness distribution

    use glide_thckCommon, only: glide_calclsrf
    use glimmer_global, only : dp
    use glimmer_slap, only: slapMatrix_type, slapMatrix_init, slapMatrix_insertElement, slapSolve

    implicit none

    ! subroutine arguments -------------------------------------------------------------

    type(glide_global_type) :: model
    
    logical,intent(in) :: calc_rhs                      !< set to true when rhs should be calculated 
                                                        !! i.e. when doing lin solution or first picard iteration
    integer,  intent(in), dimension(:,:) :: mask        !< Index mask for matrix mapping
    integer,  intent(in)                 :: totpts      !< Number of non-zero points in mask
    real(dp), intent(in), dimension(:,:) :: old_thck    !< contains ice thicknesses from previous time step
    real(dp), intent(inout), dimension(:,:) :: new_thck !< on entry contains first guess for new ice thicknesses
                                                        !< on exit contains ice thicknesses of new time step
    real(dp),intent(inout),dimension(totpts) :: rhsd
    integer,  intent(in)                 :: periodic_ew !< Periodic boundary conditions flag

    ! local variables ------------------------------------------------------------------

    real(dp), dimension(5) :: sumd 
    real(dp) :: err
    integer :: linit
    integer :: ew,ns
    integer :: ewn,nsn
    real(dp),dimension(totpts) :: answ

    type(slapMatrix_type) :: matrix

    ewn = size(old_thck,1)
    nsn = size(old_thck,2)

    ! Initialise sparse matrix object
    call slapMatrix_init(matrix,totpts,ewn*nsn*5)

    answ = 0.0

    ! Boundary Conditions ---------------------------------------------------------------
    ! lower and upper BC
    do ew = 1,ewn
       ns=1
       if (mask(ew,ns) /= 0) then
          call slapMatrix_insertElement(matrix,1.0d0,mask(ew,ns),mask(ew,ns))
          if (calc_rhs) then
             rhsd(mask(ew,ns)) = old_thck(ew,ns) 
          end if
          answ(mask(ew,ns)) = new_thck(ew,ns)
       end if
       ns=nsn
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
       do ns=2,nsn-1
          ew = 1
          if (mask(ew,ns) /= 0) then
             call findsums(model%velocity%diffu, &
                  model%velocity%ubas, &
                  sumd, &
                  model%pcgdwk%fc2(1), &
                  model%pcgdwk%fc2(5), &
                  ewn-2,ewn-1,ns-1,ns)
             call generate_row(ewn-2,ew,ew+1,ns-1,ns,ns+1)
          end if
          ew=ewn
          if (mask(ew,ns) /= 0) then
             call findsums(model%velocity%diffu, &
                  model%velocity%ubas, &
                  sumd, &
                  model%pcgdwk%fc2(1), &
                  model%pcgdwk%fc2(5), &
                  1,2,ns-1,ns)
             call generate_row(ew-1,ew,3,ns-1,ns,ns+1)
          end if
       end do
    else
       do ns=2,nsn-1
          ew=1
          if (mask(ew,ns) /= 0) then
             call slapMatrix_insertElement(matrix,1.0d0,mask(ew,ns),mask(ew,ns))
             if (calc_rhs) then
                rhsd(mask(ew,ns)) = old_thck(ew,ns) 
             end if
             answ(mask(ew,ns)) = new_thck(ew,ns)
          end if
          ew=ewn
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

    do ns = 2,nsn-1
       do ew = 2,ewn-1

          if (mask(ew,ns) /= 0) then
                
             call findsums(model%velocity%diffu, &
                  model%velocity%ubas, &
                  sumd, &
                  model%pcgdwk%fc2(1), &
                  model%pcgdwk%fc2(5), &
                  ew-1,ew,ns-1,ns)
             call generate_row(ew-1,ew,ew+1,ns-1,ns,ns+1)

          end if
       end do
    end do

    ! Solve the system using SLAP
    call slapSolve(matrix,rhsd,answ,linit,err)   

    ! Rejig the solution onto a 2D array
    do ns = 1,nsn
       do ew = 1,ewn 

          if (mask(ew,ns) /= 0) then
             new_thck(ew,ns) = answ(mask(ew,ns))
          end if

       end do
    end do

    new_thck = max(0.0d0, new_thck)

#ifdef DEBUG
    print *, "* thck ", model%numerics%time, linit, totpts, &
         real(thk0*new_thck(ewn/2+1,nsn/2+1)), &
         real(vel0*maxval(abs(model%velocity%ubas))), real(vel0*maxval(abs(model%velocity%vbas))) 
#endif

    ! calculate upper and lower surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

  contains

    subroutine generate_row(ewm,ew,ewp,nsm,ns,nsp)
      ! calculate row of sparse matrix equation
      implicit none
      integer, intent(in) :: ewm,ew,ewp  ! ew index to left, central, right node
      integer, intent(in) :: nsm,ns,nsp  ! ns index to lower, central, upper node

      ! fill sparse matrix
      call slapMatrix_insertElement(matrix,sumd(1),mask(ewm,ns),mask(ew,ns))       ! point (ew-1,ns)
      call slapMatrix_insertElement(matrix,sumd(2),mask(ewp,ns),mask(ew,ns))       ! point (ew+1,ns)
      call slapMatrix_insertElement(matrix,sumd(3),mask(ew,nsm),mask(ew,ns))       ! point (ew,ns-1)
      call slapMatrix_insertElement(matrix,sumd(4),mask(ew,nsp),mask(ew,ns))       ! point (ew,ns+1)
      call slapMatrix_insertElement(matrix,1.0d0 + sumd(5),mask(ew,ns),mask(ew,ns))! point (ew,ns)

      ! calculate RHS
      if (calc_rhs) then
         rhsd(mask(ew,ns)) =                    &
              old_thck(ew,ns) * (1.0d0 - model%pcgdwk%fc2(3) * sumd(5))     &
            - model%pcgdwk%fc2(3) * (old_thck(ewm,ns) * sumd(1)             &
                                   + old_thck(ewp,ns) * sumd(2)             &
                                   + old_thck(ew,nsm) * sumd(3)             &
                                   + old_thck(ew,nsp) * sumd(4))            &
            - model%pcgdwk%fc2(4) * (model%geometry%lsrf(ew,ns)  * sumd(5)  &
                                   + model%geometry%lsrf(ewm,ns) * sumd(1)  &
                                   + model%geometry%lsrf(ewp,ns) * sumd(2)  &
                                   + model%geometry%lsrf(ew,nsm) * sumd(3)  &
                                   + model%geometry%lsrf(ew,nsp) * sumd(4)) &
            + model%climate%acab(ew,ns) * model%pcgdwk%fc2(2)
         if(model%options%basal_mbal==1) then
            rhsd(mask(ew,ns)) =                    &
                 rhsd(mask(ew,ns))                 &
                 - model%temper%bmlt(ew,ns) * model%pcgdwk%fc2(2) ! basal melt is +ve for mass loss
         end if
      end if

      answ(mask(ew,ns)) = new_thck(ew,ns)      

    end subroutine generate_row

  end subroutine thck_evolve

!-------------------------------------------------------------------------

  subroutine findsums(diffu,ubas,sumd,fc2_1,fc2_5,ewm,ew,nsm,ns)
  
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

