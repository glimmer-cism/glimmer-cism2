! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_glenflow.f90 - part of the Glimmer-CISM ice model  + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004-9 Glimmer-CISM contributors - see COPYRIGHT file 
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
! The Glimmer-CISM maintainer is:
!
! Ian Rutt
! School of the Environment and Society
! Swansea University
! Singleton Park
! Swansea
! SA2 8PP
! UK
!
! email: <i.c.rutt@swansea.ac.uk> or <ian.rutt@physics.org>
!
! Glimmer-CISM is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide_glenflow

  use glimmer_global, only: dp

  implicit none

  type glenflow_params
     real(dp) :: fact1  ! Value of a when T* is above -263K
     real(dp) :: fact2  ! Value of a when T* is below -263K
     real(dp) :: fact3  ! Value of -Q/R when T* is above -263K
     real(dp) :: fact4  ! Value of -Q/R when T* is below -263K
     real(dp) :: fiddle
  end type glenflow_params

  private
  public :: glenflow_params, glenflow_init, calcflwa

contains

  subroutine glenflow_init(params,fiddle)

    use physcon,  only: arrmlh, arrmll, actenh, actenl, gascon
    use paramets, only: vis0

    implicit none

    type(glenflow_params),intent(out) :: params
    real(dp),             intent(in)  :: fiddle

    params%fact1  = fiddle * arrmlh / vis0   ! Value of a when T* is above -263K
    params%fact2  = fiddle * arrmll / vis0   ! Value of a when T* is below -263K
    params%fact3  = -actenh / gascon         ! Value of -Q/R when T* is above -263K
    params%fact4  = -actenl / gascon         ! Value of -Q/R when T* is below -263K
    params%fiddle = fiddle

  end subroutine glenflow_init

  !------------------------------------------------------------------------------------------

  subroutine calcflwa(params,flwa,temp,thck,flag,thklim,sigma)

    !*FD Calculates Glenn's $A$ over the three-dimensional domain,
    !*FD using one of three possible methods.
    !*FD \textbf{I'm unsure how this ties in with the documentation, since}
    !*FD \texttt{fiddle}\ \textbf{is set to 3.0. This needs checking} 

    use physcon, only : pmlt, rhoi, grav
    use glimmer_global, only: dp
    use paramets, only : thk0

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glenflow_params),      intent(in)    :: params    !*FD Derived type containing
                                                           !*FD parameters for this module
    real(dp),dimension(:,:,:),  intent(out)   :: flwa      !*FD The calculated values of $A$
    real(dp),dimension(:,0:,0:),intent(in)    :: temp      !*FD The 3D temperature field
    real(dp),dimension(:,:),    intent(in)    :: thck      !*FD The ice thickness
    integer,                    intent(in)    :: flag      !*FD Flag to select the method
                                                           !*FD of calculation:
    !*FD \begin{description}
    !*FD \item[0] {\em Paterson and Budd} relationship.
    !*FD \item[1] {\em Paterson and Budd} relationship, with temperature set to
    !*FD -5$^{\circ}$C.
    !*FD \item[2] Set constant, {\em but not sure how this works at the moment\ldots}
    !*FD \end{description}

    real(dp),                   intent(in)    :: thklim    !*FD Thickness limit for calculation
    real(dp),dimension(:),      intent(in)    :: sigma     !*FD Sigma vertical levels

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp), parameter :: fact = grav * rhoi * pmlt * thk0
    real(dp), parameter :: contemp = -5.0d0  
    real(dp), dimension(size(sigma)) :: tempcor

    integer :: ew,ns,up,ewn,nsn,upn

    !------------------------------------------------------------------------------------
    
    upn=size(flwa,1) ; ewn=size(flwa,2) ; nsn=size(flwa,3)

    !------------------------------------------------------------------------------------

    select case(flag)
    case(0)

      ! This is the Paterson and Budd relationship

      do ns = 1,nsn
        do ew = 1,ewn
          if (thck(ew,ns) > thklim) then
            ! Calculate the corrected temperature

            tempcor = min(0.0d0, temp(:,ew,ns) + thck(ew,ns) * fact * sigma)
            tempcor = max(-50.0d0, tempcor)

            ! Calculate Glenn's A
            call patebudd(tempcor,flwa(:,ew,ns),params)
          else
            flwa(:,ew,ns) = params%fiddle
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

            call patebudd((/(contemp, up=1,upn)/),flwa(:,ew,ns),params) 
          else
            flwa(:,ew,ns) = params%fiddle
          end if
        end do
      end do

    case default 

      ! Set A equal to the value of fiddle. According to the documentation, this
      ! option means A=10^-16 yr^-1 Pa^-n, but I'm not sure how this squares with
      ! the value of fiddle, which is currently set to three.

      flwa = params%fiddle
  
    end select

  end subroutine calcflwa

  !------------------------------------------------------------------------------------------

  subroutine patebudd(tempcor,calcga,params)

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

    use physcon, only : trpt
    use glimmer_global, only: dp

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    real(dp),dimension(:), intent(in)    :: tempcor  !*FD Input temperature profile. This is 
                                                     !*FD {\em not} $T^{*}$, as it has $T_0$
                                                     !*FD added to it later on; rather it is
                                                     !*FD $T-T_{\mathrm{pmp}}$.
    real(dp),dimension(:), intent(out)   :: calcga   !*FD The output values of Glenn's $A$.
    type(glenflow_params)                :: params   !*FD Constants for the calculation. 

    !------------------------------------------------------------------------------------
    ! Actual calculation is done here - constants depend on temperature -----------------

    where (tempcor >= -10.0d0)         
      calcga = params%fact1 * exp(params%fact3 / (tempcor + trpt))
    elsewhere
      calcga = params%fact2 * exp(params%fact4 / (tempcor + trpt))
    end where

  end subroutine patebudd

end module glide_glenflow
