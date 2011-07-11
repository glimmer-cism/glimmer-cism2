! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_temp_utils.f90 - part of the GLIMMER ice model     +
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!
! Copyright (C) 2010 GLIMMER contributors - see COPYRIGHT file 
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
! GLIMMER is hosted on NeSCForge:
!
! http://forge.nesc.ac.uk/projects/glimmer/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!whl - Created this module from glide_temp.F90, Dec. 2010
! This module includes useful subroutines that previously were part
! of glide_temp.F90.  I put them in a separate module so they can
! be called from an alternative temperature driver, glissade_temp.F90,
! without duplicating code.

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide_temp_utils

  use glide_types

contains

  !-----------------------------------------------------------------------

  subroutine finddisp(model,thck,whichdisp,efvs,stagthck,dusrfdew,dusrfdns,flwa)

    use glimmer_global, only : dp
    use glimmer_physcon, only : gn

    implicit none

    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(in) :: thck, stagthck, dusrfdew, dusrfdns
    real(dp), dimension(:,:,:), intent(in) :: flwa, efvs
    integer, intent(in) :: whichdisp

    integer, parameter :: p1 = gn + 1  
    integer :: ew, ns
    integer :: iew, ins    !*sfp* for HO and SSA dissip. calc.

    real(dp) :: c2 

    ! two methods of doing this. 
    ! 1. find dissipation at u-pts and then average
    ! 2. find dissipation at H-pts by averaging quantities from u-pts
    ! 2. works best for eismint divide (symmetry) but 1 likely to be better for full expts

    model%tempwk%dissip(:,:,:) = 0.0d0

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1
          if (thck(ew,ns) > model%numerics%thklim) then
             
             c2 = (0.25*sum(stagthck(ew-1:ew,ns-1:ns)) * dsqrt((0.25*sum(dusrfdew(ew-1:ew,ns-1:ns)))**2 &
                  + (0.25*sum(dusrfdns(ew-1:ew,ns-1:ns)))**2))**p1
             
             model%tempwk%dissip(:,ew,ns) = c2 * model%tempwk%c1(:) * ( &
                  flwa(:,ew-1,ns-1) + flwa(:,ew-1,ns+1) + flwa(:,ew+1,ns+1) + flwa(:,ew+1,ns-1) + &
                  2*(flwa(:,ew-1,ns)+flwa(:,ew+1,ns)+flwa(:,ew,ns-1)+flwa(:,ew,ns+1)) + &
                  4*flwa(:,ew,ns))             
          end if
       end do
    end do

  end subroutine finddisp

  !-----------------------------------------------------------------------------------

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
    !    set temperature to pressure melting point 
    ! 2. if bed is wet set basal temperature to pressure melting point 

    call calcpmpt(pmptemp,thck,sigma)

    temp = dmin1(temp,pmptemp)

    if (bwat > 0.0d0) temp(upn) = pmptemp(upn)

  end subroutine corrpmpt

  !-------------------------------------------------------------------

  subroutine calcpmpt(pmptemp,thck,sigma)

    use glimmer_global, only : dp !, upn
    use glimmer_physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    implicit none 

    real(dp), dimension(:), intent(out) :: pmptemp
    real(dp), intent(in) :: thck
    real(dp),intent(in),dimension(:) :: sigma

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp(:) = fact * thck * sigma(:)

  end subroutine calcpmpt

  !-------------------------------------------------------------------

  subroutine calcpmptb(pmptemp,thck)

    use glimmer_global, only : dp
    use glimmer_physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    implicit none 

    real(dp), intent(out) :: pmptemp
    real(dp), intent(in) :: thck

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp = fact * thck 

  end subroutine calcpmptb

  !-------------------------------------------------------------------

  subroutine calcflwa(sigma, thklim, flwa, temp, thck, flow_factor, default_flwa_arg, flag)

    !*FD Calculates Glenn's $A$ over the three-dimensional domain,
    !*FD using one of three possible methods.

    use glimmer_physcon
    use glimmer_paramets, only : thk0, vis0

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

!      Note: dummy argument sigma can be either model%numerics%sigma or model%numerics%stagsigma.
!      The flwa, temp, and sigma arrays must have the same vertical dimension
!       (1:upn on an unstaggered vertical grid, or 1:upn-1 on a staggered vertical grid).
!
    real(dp),dimension(:),      intent(in)    :: sigma     !*FD Vertical coordinate
    real(dp),                   intent(in)    :: thklim    !*FD thickness threshold
    real(dp),dimension(:,:,:),  intent(out)   :: flwa      !*FD The calculated values of $A$
!whl - Changed this so that the temp array is assumed to start with horizontal index 1 instead of 0.
!!    real(dp),dimension(:,0:,0:),intent(in)    :: temp      !*FD The 3D temperature field
    real(dp),dimension(:,:,:),intent(in)      :: temp      !*FD The 3D temperature field
    real(dp),dimension(:,:),    intent(in)    :: thck      !*FD The ice thickness
    real(dp)                                  :: flow_factor !*FD Fudge factor in arrhenius relationship
    real(dp),                   intent(in)    :: default_flwa_arg !*FD Glen's A to use in isothermal case 
    integer,                    intent(in)    :: flag      !*FD Flag to select the method
                                                           !*FD of calculation:
    !*FD \begin{description}
    !*FD \item[0] {\em Paterson and Budd} relationship.
    !*FD \item[1] {\em Paterson and Budd} relationship, with temperature set to
    !*FD -5$^{\circ}$C.
    !*FD \item[2] Set constant, {\em but not sure how this works at the moment\ldots}
    !*FD \end{description}

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp), parameter :: fact = grav * rhoi * pmlt * thk0
    real(dp), parameter :: contemp = -5.0d0
    real(dp) :: default_flwa
    real(dp),dimension(4) :: arrfact
    integer :: ew,ns,up,ewn,nsn,upn

    real(dp), dimension(size(sigma)) :: tempcor

    !------------------------------------------------------------------------------------ 
   
    default_flwa = flow_factor * default_flwa_arg / (vis0*scyr) 

    upn=size(flwa,1) ; ewn=size(flwa,2) ; nsn=size(flwa,3)

    !write(*,*)"Default flwa = ",default_flwa

    arrfact = (/ flow_factor * arrmlh / vis0, &   ! Value of a when T* is above -263K
                 flow_factor * arrmll / vis0, &   ! Value of a when T* is below -263K
                 -actenh / gascon,        &       ! Value of -Q/R when T* is above -263K
                 -actenl / gascon/)               ! Value of -Q/R when T* is below -263K
 
    select case(flag)
    case(0)

      ! This is the Paterson and Budd relationship

      do ns = 1,nsn
        do ew = 1,ewn
          if (thck(ew,ns) > thklim) then
            

            ! Calculate the corrected temperature

            tempcor(:) = min(0.0d0, temp(:,ew,ns) + thck(ew,ns) * fact * sigma(:))
            tempcor(:) = max(-50.0d0, tempcor(:))

            ! Calculate Glenn's A

            call patebudd(tempcor(:), flwa(:,ew,ns), arrfact) 
          else
            flwa(:,ew,ns) = default_flwa
          end if
        end do
      end do

    case(1)

      ! This is the Paterson and Budd relationship, but with the temperature held constant
      ! at -5 deg C

      do ns = 1,nsn
        do ew = 1,ewn
          if (thck(ew,ns) > thklim) then

            ! Calculate Glenn's A with a fixed temperature.

            call patebudd((/(contemp, up=1,upn)/),flwa(:,ew,ns),arrfact) 
          else
            flwa(:,ew,ns) = default_flwa
          end if
        end do
      end do

    case default 

      flwa(:,:,:) = default_flwa
  
    end select

  end subroutine calcflwa 

!------------------------------------------------------------------------------------------

  subroutine patebudd(tempcor,calcga,fact)

    !*FD Calculates the value of Glenn's $A$ for the temperature values in a one-dimensional
    !*FD array. The input array is usually a vertical temperature profile. The equation used
    !*FD is from \emph{Paterson and Budd} [1982]:
    !*FD \[
    !*FD A(T^{*})=a \exp \left(\frac{-Q}{RT^{*}}\right)
    !*FD \]
    !*FD This is equation 9 in {\em Payne and Dongelmans}. $a$ is a constant of proportionality,
    !*FD $Q$ is the activation energy for for ice creep, and $R$ is the universal gas constant.
    !*FD The pressure-corrected temperature, $T^{*}$ is given by:
    !*FD \[
    !*FD T^{*}=T-T_{\mathrm{pmp}}+T_0
    !*FD \] 
    !*FD \[
    !*FD T_{\mathrm{pmp}}=T_0-\sigma \rho g H \Phi
    !*FD \]
    !*FD $T$ is the ice temperature, $T_{\mathrm{pmp}}$ is the pressure melting point 
    !*FD temperature, $T_0$ is the triple point of water, $\rho$ is the ice density, and 
    !*FD $\Phi$ is the (constant) rate of change of melting point temperature with pressure.

    use glimmer_physcon, only : trpt

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    real(dp),dimension(:), intent(in)    :: tempcor  !*FD Input temperature profile. This is 
                                                     !*FD {\em not} $T^{*}$, as it has $T_0$
                                                     !*FD added to it later on; rather it is
                                                     !*FD $T-T_{\mathrm{pmp}}$.
    real(dp),dimension(:), intent(out)   :: calcga   !*FD The output values of Glenn's $A$.
    real(dp),dimension(4), intent(in)    :: fact     !*FD Constants for the calculation. These
                                                     !*FD are set when the velo module is initialised

    !------------------------------------------------------------------------------------

    ! Actual calculation is done here - constants depend on temperature -----------------

    where (tempcor >= -10.0d0)         
      calcga = fact(1) * exp(fact(3) / (tempcor + trpt))
    elsewhere
      calcga = fact(2) * exp(fact(4) / (tempcor + trpt))
    end where

  end subroutine patebudd

!------------------------------------------------------------------------------------------

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

end module glide_temp_utils
