! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_damage.f90 - part of the Glimmer-CISM ice model  + 
! +  glide_damage.f90 - part of the GLIMMER ice model         +
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010
! Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
! This file is part of Glimmer-CISM.
!
! Glimmer-CISM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or (at
! your option) any later version.
!
! Glimmer-CISM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
!
! Glimmer-CISM is hosted on BerliOS.de:
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! New file for damage mechanics stuff

module glide_damage

  use glide_types

  private

  public :: glide_damage_init, update_damage

contains

  subroutine glide_damage_init(model)

    implicit none
    type(glide_global_type), intent(inout) :: model  !*FD Ice model parameters.    

    integer, parameter :: which_scenario = 1  ! 0 => No scenario. Initialize all damage to zero.
                                              ! 1 => One point of damage 
                                              ! 2 => One column of damage
                                              ! 3 => One east-west strip of damage
                                              ! 4 => One north-south strip of damage
                                              ! 5 => A square/rectangle of damage
                                              ! 6 => Some initial distribution of damage

    ! Set all damage, fracture  to zero. 
    model%damage%sclr_damage(:,:,:) = 0.0d0
    model%damage%fractured = 0

    ! Test Case 1: Square of damage
    if (which_scenario==1) then
      model%damage%sclr_damage(:,15:20,25:27) = 0.1d0

    ! Test Case 2: A column of damage
    elseif (which_scenario==2) then
      model%damage%sclr_damage(:,50,50) = 0.5d0

    ! Test Case 3: One strip of east-west damage
    elseif (which_scenario==3) then
      model%damage%sclr_damage(2,:,50) = 0.5d0

    ! Test Case 4: One strip of north-south damage
    elseif (which_scenario==4) then
      model%damage%sclr_damage(2,50,:) = 0.5d0

    ! Test Case 5: A square of damage
    elseif (which_scenario==5) then
      model%damage%sclr_damage(2, 20:25, 20   ) = 0.5d0
      model%damage%sclr_damage(2, 20:25, 25   ) = 0.5d0
      model%damage%sclr_damage(2, 20   , 20:25) = 0.5d0
      model%damage%sclr_damage(2, 25   , 20:25) = 0.5d0

    ! Test Case 6: Some initial distribution of damage
    elseif (which_scenario==6) then
      print *, 'No intial damage distribution set'
      model%damage%sclr_damage = 0.0d0 ! To do
    end if

  end subroutine glide_damage_init


  subroutine update_damage(model)

    !use glimmer_physcon, only: rhoi, grav
    use glimmer_paramets, only: thk0,tau0
    use glimmer_global, only: dp

    implicit none

    !--------------------------------------------------
    ! Subroutine arguments
    !--------------------------------------------------
    type(glide_global_type), intent(inout) :: model

    !--------------------------------------------------
    ! Internal variables
    !--------------------------------------------------
    real(kind = dp), dimension(:,:,:), allocatable :: damageSource
    integer :: ewn, nsn, upn
    integer :: ew, ns, up

    !--------------------------------------------------
    ! Parameters
    !--------------------------------------------------
    logical, parameter :: use_healing = .FALSE.

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn-1

    allocate(damageSource(upn-1,ewn,nsn))
    damageSource = 0.0d0

    ! Step 01: Calculate dynamic function of damage at each point
    do ew=1, ewn
      do ns = 1, nsn
        ! If there is ice...
        if (model%geometry%mask(ew,ns)> 0) then
          do up = 1, upn

            damageSource(up,ew,ns) = damage_source(model%stress%tau%scalar(up,ew,ns)*tau0, &
                                                   model%stress%tau%xx(up,ew,ns)*tau0, &
                                                   model%stress%tau%yy(up,ew,ns)*tau0, &
                                                   model%stress%tau%xy(up,ew,ns)*tau0, &
                                                   model%stress%tau%yz(up,ew,ns)*tau0, &
                                                   model%stress%tau%xz(up,ew,ns)*tau0, &
                                                   model%damage%sclr_damage(up,ew,ns), &
                                                   model%geometry%thck(ew,ns)*thk0, &
                                                   up, & 
                                                   use_healing)

            model%damage%sclr_damage(up,ew,ns) = model%damage%sclr_damage(up,ew,ns) + damageSource(up,ew,ns)            

            if(model%damage%sclr_damage(up,ew,ns) >= 0.56d0) then
              model%damage%sclr_damage(up,ew,ns) = 1.0d0
              model%damage%fractured = 1
            end if
           
            if (model%damage%sclr_damage(up,ew,ns) < 0.0d0) then
              print *, 'Warning: Damage < 0.0 at up,ew,ns', up,ew,ns, '. Setting damage to 0.0\n'
              model%damage%sclr_damage = 0.0d0
            end if

          end do !upn
        else
          model%damage%sclr_damage(:,ew,ns) = 0.0d0 ! Set damage to zero if no ice
        end if
      end do !nsn
    end do !ewn

    deallocate(damageSource)
  end subroutine update_damage


  function damage_source(tau, tau_xx, tau_yy, tau_xy, tau_yz, tau_xz, damage,iceThickness, up, healing)

    !-------------------------------------------------------------------------------------+
    ! To Do: 1. Calculate parameters: alpha,beta,B,k1,r,stress_threshold                  +  
    !        2. Investigate Fracture Zone Characteristic Size (FZ_CharacSize) values      +
    !-------------------------------------------------------------------------------------+

    use glimmer_physcon,  only: rhoi, grav, pi
    use glimmer_paramets, only: len0, thk0, tau0
    use glimmer_global,   only: dp

    implicit none

    !--------------------------------------------------
    ! Subroutine arguments
    !--------------------------------------------------
    real (kind = dp), intent(in) :: tau, tau_xx, tau_yy, tau_xy, tau_yz, tau_xz, damage, iceThickness
    integer, intent(in) :: up
    logical, intent(in) :: healing

    real (kind = dp) :: damage_source

    !--------------------------------------------------
    ! Internal parameters
    !--------------------------------------------------
    real(dp), parameter :: alpha                = 0.21d0             !
    real(dp), parameter :: beta                 = 0.63d0             !
    real(dp), parameter :: B_parameter          = 1.7d0 * 10.0d-9    !
    real(dp), parameter :: healing_parameter    = 0.4d0              ! lambda_h in Pralong05
    real(dp), parameter :: k1_parameter         = 3.75d0 * 10.0d-3   ! 
    real(dp), parameter :: r                    = 0.43d0             ! 
    real(dp), parameter :: stressThreshold_Ref  = 3.3d5              ! Reference Stress Threshold
    real(dp), parameter :: FPZ_size             = 10.0d0             ! Fracture Process Zone size = L_f
    real(dp), parameter :: FZ_CharacSize        = 1000.0d0             ! FIX ME... Characteristic size of Fracture Zone = L_g
    real(dp), parameter :: weibullParam         = 8.0d0              ! Weibull parameter
    real(dp), parameter :: dimSimilarity        = 3.0d0              ! Dimensional similarity parameter
    real(dp), parameter :: refLength            = 10.d0               ! Reference Length

    !--------------------------------------------------
    ! Internal variables
    !--------------------------------------------------
    real(dp) :: max_principal_stress
    real(dp) :: invariant_I, invariant_II_dev
    real(dp) :: sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz   !effective stresses
    real(dp) :: tau_zz, hydrostatic_pressure	                             !zz-comp deviatoric stress tensor
    real(dp) :: psi, k_parameter, source, phi, coeff                         !k_parameter: damage source magnification
    real(dp) :: stressThreshold, characSize
    integer :: ewn, nsn, upn

    real(dp) :: stress_tensor(3,3)   ! Effective Stress Tensor  
    real(dp) :: principal_stress(3)  ! LAPACK - returns eigenvalues of stress_tensor to this
    real(dp) :: work(6)              ! LAPACK - solver vars
    integer  :: lwork, lwmax, info   ! LAPACK - solver vars

    !--------------------------------------------------
    ! Calculations
    !--------------------------------------------------

    ! Hydrostatic Pressure
    !                  P = (-1/3)*(sigma_xx + sigma_yy + sigma_zz)  pressure
    !           sigma_zz = tau_zz - P                               deviatoric stress
    !    (d sigma_zz)/dz = (d tau_zz)/dz - dP/dz = rho_ice * grav   hydrostatic approx
    !                  P = tau_zz + rho_ice*grav*(s-z)              integrate to surface
    !
    !              sigma = up = (s-z)/iceThickness                  (s-z) calc
    !              (s-z) = iceThickness*up 

    tau_zz = -1.0d0 * (tau_xx + tau_yy)
    hydrostatic_pressure = rhoi*grav*(iceThickness*up) + tau_zz

    ! Effective Stresses
    sigma_xx = (tau_xx - hydrostatic_pressure) / ((1.0d0 - damage) + epsilon(1.d0))
    sigma_yy = (tau_yy - hydrostatic_pressure) / ((1.0d0 - damage) + epsilon(1.d0))
    sigma_zz = (tau_zz - hydrostatic_pressure) / ((1.0d0 - damage) + epsilon(1.d0) )  !Already unscaled
    sigma_xy = tau_xy  / (1.0d0 - damage + epsilon(1.d0))
    sigma_yz = tau_yz  / (1.0d0 - damage + epsilon(1.d0))
    sigma_xz = tau_xz  / (1.0d0 - damage + epsilon(1.d0))

    ! Effective Cauchy Stress Tensor
    stress_tensor(1,1) = sigma_xx
    stress_tensor(1,2) = sigma_xy
    stress_tensor(1,3) = sigma_xz
    stress_tensor(2,1) = stress_tensor(1,2) !sigma_xy
    stress_tensor(2,2) = sigma_yy
    stress_tensor(2,3) = sigma_yz
    stress_tensor(3,1) = stress_tensor(1,3) !sigma_xz
    stress_tensor(3,2) = stress_tensor(2,3) !sigma_yz
    stress_tensor(3,3) = sigma_zz

!    LAPACK: Test symmetric matrix
!    stress_tensor(1,1) = 3.0d0
!    stress_tensor(1,2) = 2.0d0
!    stress_tensor(1,3) = 4.0d0
!    stress_tensor(2,1) = stress_tensor(1,2)
!    stress_tensor(2,2) = 0.0d0
!    stress_tensor(2,3) = 2.0d0
!    stress_tensor(3,1) = stress_tensor(1,3)
!    stress_tensor(3,2) = stress_tensor(2,3)
!    stress_tensor(3,3) = 3.0d0

    ! LAPACK: test vars
    lwork = -1
    lwmax = 1000
    info  = 0

    ! LAPACK: Query the optimal workspace.
    call DSYEV('N', 'U', 3, stress_tensor, 3, principal_stress, work, lwork, info) 
    lwork = min(lwmax,int(work(1)))

    ! LAPACK: Solve eigenproblem
    call DSYEV('N', 'U', 3, stress_tensor, 3, principal_stress, work, lwork, info)

    ! LAPACK: Check for convergence
    if (info .gt. 0) then
      print *, 'Algorithm failed to compute eigenvalues.'
    end if

    ! Maximum Principal Stress
    max_principal_stress = max(principal_stress(1),principal_stress(2),principal_stress(3))

    ! Stress invariants
    invariant_I  = sigma_xx + sigma_yy + sigma_zz 
    invariant_II_dev = 0.5d0*(tau_xx)**2.0d0 + 0.5d0*(tau_yy)**2.0d0 + 0.5d0*tau_zz**2.0d0 + (tau_xy)**2.0d0 + (tau_yz)**2.0d0 + (tau_xz)**2.0d0
    !old calc: invariant_II_dev = 0.5d0*tau_xx**2.0d0 + 0.5d0*tau_yy**2.0d0 + 0.5d0*tau_zz**2.0d0 + tau_xy**2.0d0 + tau_yz**2.0d0 + tau_xz**2.0d0


    ! Calculation of k_parameter (damage magnification parameter)
    ! As k up (down) --> damage source up (down)
    if(healing) then
      if (invariant_I >= 0) then
        k_parameter = k1_parameter * sqrt(invariant_I)
      else
        k_parameter = k1_parameter * -1.0d0 * healing_parameter * sqrt(abs(invariant_I))
      end if
    else
        k_parameter = k1_parameter * sqrt(abs(invariant_I))
    end if

    ! Stress Threshold calculation
    characSize = (FPZ_size*FZ_characSize**2.0d0)**(1.0d0/3.0d0)
    stressThreshold = stressThreshold_ref * (characSize/refLength)**(-dimSimilarity/weibullParam)

    ! Pralong and Funk (2005): Equation (25)
    psi = (alpha*max_principal_stress + beta*sqrt(abs(3.0d0*invariant_II_dev)) + (1.0d0-alpha-beta)*invariant_I) - stressThreshold
    !psi = 1.0d0 / (1.0d0 - damage)            ! FIX ME: Where did this come from in the paper?

    ! Pralong and Funk (2005): Equation (24)
    if (healing) then
      if (psi >= 0.0d0) then
        source = B_parameter * (psi**r) * ((1.0d0-damage)**(-1.0d0*k_parameter))
      else
        source = B_parameter * -1.0d0 * healing_parameter * ((-1.0d0*psi)**r) * (1.0d0-damage)**(-1.0d0*k_parameter)
      end if
    else
      source = B_parameter * (max(psi,0.0d0)**r) * ((1.0d0-damage)**(-1.0d0*k_parameter))
    end if

    ! Check for bad damage source term. 
    if (isnan(source)) then
      print *, 'Bad damage source. Setting damage source to zero.' 
      source = 0.0d0
    end if

    damage_source = source
    return
  end function damage_source
end module glide_damage

  
