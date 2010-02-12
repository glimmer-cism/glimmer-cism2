! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_lithot.f90 - part of the GLIMMER ice model         + 
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

!> module for temperature calculations in the upper lithosphere
!!
!! \author Magnus Hagdorn
!! \date 2006

module lithot

  use glimmer_global, only : dp,sp
  use lithot_types
  use lithot_setup

contains  
  !> initialise geothermal heat flux computations
  subroutine init_lithot(litho,dt,hotstart)
    use glimmer_paramets, only: tim0
    use glimmer_log
    use lithot1d
    implicit none
    type(lithot_type) :: litho            !< structure holding bedrock temperature configuration    
    real(dp), intent(in) :: dt !< current time step
    logical, intent(in) :: hotstart !< set to .True. if model is hotstarted

    ! local variables
    integer k
    real(kind=dp) :: factor

    ! allocate memory for common arrays
    allocate(litho%deltaz(litho%nlayer)); litho%deltaz = 0.0
    allocate(litho%zfactors(3,litho%nlayer)); litho%zfactors = 0.0    

    ! set up vertical grid
    do k=1,litho%nlayer
       litho%deltaz(k) = (1-calc_sigma(real(litho%nlayer-k)/real(litho%nlayer-1),2.)) &
            *litho%rock_base
    end do

    ! calculate diffusion coefficient
    litho%diffu = litho%con_r/(litho%rho_r*litho%shc_r)

    ! set up factors for vertical finite differences
    do k=2,litho%nlayer-1
       litho%zfactors(1,k) =  litho%diffu*tim0*dt / &
            ((litho%deltaz(k)-litho%deltaz(k-1)) * (litho%deltaz(k+1)-litho%deltaz(k-1)))
       litho%zfactors(2,k) = litho%diffu*tim0*dt / &
            ((litho%deltaz(k+1)-litho%deltaz(k)) * (litho%deltaz(k)-litho%deltaz(k-1)))
       litho%zfactors(3,k) = litho%diffu*tim0*dt / &
            ((litho%deltaz(k+1)-litho%deltaz(k)) * (litho%deltaz(k+1)-litho%deltaz(k-1)))
    end do
    k = litho%nlayer
    litho%zfactors(:,k) = 0.5*litho%diffu*tim0*dt / &
         (litho%deltaz(k)-litho%deltaz(k-1))**2

    if (.not.hotstart) then
       ! set initial temp distribution to thermal gradient
       factor = litho%geot/litho%con_r
       do k=1,litho%nlayer
          litho%temp(:,:,k) = litho%surft+litho%deltaz(k)*factor
       end do
    end if


    if (litho%num_dim.eq.1) then
       call init_lithot1d(litho)
    else
       call write_log('Wrong number of dimensions.',GM_FATAL,__FILE__,__LINE__)
    end if

    contains
      function calc_sigma(x,n)
        implicit none
        real :: calc_sigma,x,n
        
        calc_sigma = (1-(x+1)**(-n))/(1-2**(-n))
      end function calc_sigma

  end subroutine init_lithot    

  !> spinup geothermal heat flux computations
  subroutine spinup_lithot(litho,hotstart,thkmask,base_temp,air_temp,bheatflx)
    use glimmer_log
    use glimmer_mask
    implicit none
    type(lithot_type) :: litho            !< structure holding bedrock temperature configuration    
    logical, intent(in) :: hotstart !< set to .True. if model is hotstarted
    integer, dimension(:,:), intent(in) :: thkmask !< surface type mask
    real(dp), dimension(:,:), intent(in) :: base_temp !< temperature at ice base
    real(sp), dimension(:,:), intent(in) :: air_temp  !< air temperature    
    real(dp),dimension(:,:), intent(out) :: bheatflx !< field holding basal heat flux
    integer t

    if (.not.hotstart .and. litho%numt .gt. 0) then
       call write_log('Spinning up GTHF calculations',type=GM_INFO)

       do t=1,litho%numt
          call calc_lithot(litho,thkmask,base_temp,air_temp,bheatflx)
       end do

    end if
  end subroutine spinup_lithot

  !> compute temperature in the lithosphere
  subroutine calc_lithot(litho,thkmask,base_temp,air_temp,bheatflx)
    use glimmer_log
    use lithot1d
    implicit none
    type(lithot_type) :: litho            !< structure holding bedrock temperature configuration    
    integer, dimension(:,:), intent(in) :: thkmask !< surface type mask
    real(dp), dimension(:,:), intent(in) :: base_temp !< temperature at ice base
    real(sp), dimension(:,:), intent(in) :: air_temp  !< air temperature    
    real(dp),dimension(:,:), intent(out) :: bheatflx !< field holding basal heat flux

    if (litho%num_dim.eq.1) then
       call calc_lithot1d(litho,thkmask,base_temp,air_temp)
    else
       call write_log('Wrong number of dimensions.',GM_FATAL,__FILE__,__LINE__)
    end if
      
    call calc_geoth(litho,bheatflx)

  end subroutine calc_lithot

  !> calculate geothermal heat flux
  subroutine calc_geoth(litho,bheatflx)
    implicit none
    type(lithot_type) :: litho            !< structure holding bedrock temperature configuration    
    real(dp),dimension(:,:), intent(out) :: bheatflx !< field holding basal heat flux

    real(dp) factor

    factor = litho%con_r/(litho%deltaz(2)-litho%deltaz(1))
    bheatflx(:,:) = factor*(litho%temp(:,:,2)-litho%temp(:,:,1))
  end subroutine calc_geoth

  !> clean up lithot module
  subroutine finalise_lithot(litho)
    use lithot1d
    use glimmer_log
    implicit none
    type(lithot_type) :: litho            !< structure holding bedrock temperature configuration    

    deallocate(litho%deltaz)
    deallocate(litho%zfactors)

    if (litho%num_dim.eq.1) then
       call finalise_lithot1d(litho)
    else
       call write_log('Wrong number of dimensions.',GM_FATAL,__FILE__,__LINE__)
    end if
  end subroutine finalise_lithot

end module lithot
