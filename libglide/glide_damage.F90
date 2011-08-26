! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_damage.f90 - part of the Glimmer-CISM ice model    + 
! +  glide_damage.f90 - part of the GLIMMER ice model           +
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
! j = -1.0d0*i**(1.0d0/3.0d0)
! You should have received a copy of the GNU General Public License
! along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
!
! Glimmer-CISM is hosted on BerliOS.de:
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! new file for damage mechanics stuff

module glide_damage

  use glide_types

  private

  public :: glide_damage_init, update_damage

contains


  subroutine glide_damage_init(model)
    ! Initialize damage module

    ! To do: 
    !       1. Implement Crank-Nicholson solver (LAPACK?)
    !       2. Stress threshold calculation 
    !       3. Healing calculation
    !       4. Investigate sigma coordinates

    implicit none
    type(glide_global_type), intent(inout) :: model  !*FD Ice model parameters.    

    model%damage%sclr_damage(:,:,:) = 0.0d0

    !test advection, (up,ew,ns)
    !model%damage%sclr_damage(3,:,15) = 0.5d0
    !model%damage%sclr_damage(:,10,10)= 0.25d0
  end subroutine glide_damage_init


  subroutine update_damage(model)

    ! To do: 
    !       01. Implement Sigma coordinates.

    use glimmer_physcon, only: rhoi, grav
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
    real(dp) :: xderiv, yderiv, zderiv
    integer :: ewn, nsn, upn
    integer :: ew, ns, up
    real :: dew, dns, dt, uvel, vvel, wvel

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn-1

    allocate(damageSource(upn-1,ewn,nsn))
    damageSource = 0.0d0

    ! Step 01: Calculate dynamic function of damage
    do ew=1, ewn
      do ns = 1, nsn
        if (model%geometry%mask(ew,ns)> 0) then
          do up = 1, upn
            damageSource(up,ew,ns) = damage_source(model%stress%tau%scalar(up,ew,ns), &
                                                   model%stress%tau%xx(up,ew,ns), &
                                                   model%stress%tau%yy(up,ew,ns), &
                                                   model%stress%tau%xy(up,ew,ns), &
                                                   model%stress%tau%yz(up,ew,ns), &
                                                   model%stress%tau%xz(up,ew,ns), &
                                                   model%damage%sclr_damage(up,ew,ns), &
                                                   up)
          end do
        end if
      end do
    end do

    ! Step 02: Use FTCS scheme to solve damage advection.
    ! To Do:
    !        1. Convert derivatives to sigma coordinates
    !        2. Implement more stable scheme (Crank-Nicholson)...with LAPACK?
    dew = model%numerics%dew
    dns = model%numerics%dns
    dt  = model%numerics%dt


    do ew=1, ewn
      do ns = 1, nsn
        if (model%geometry%mask(ew,ns)<= 0) then
          model%damage%sclr_damage(:,ew,ns) = 0.0d0
        else
          do up = 1, upn

            uvel = model%velocity%uvel(up,ew,ns)
            vvel = model%velocity%vvel(up,ew,ns)
            wvel = model%velocity%wvel(up,ew,ns)

            ! Calculate x derivative
            if (ew==1) then
              xderiv = (model%damage%sclr_damage(up,2,ns)-model%damage%sclr_damage(up,1,ns))/dew
            elseif (ew == ewn) then
              xderiv = (model%damage%sclr_damage(up,ewn,ns)-model%damage%sclr_damage(up,ewn-1,ns))/dew
            else
              xderiv = (model%damage%sclr_damage(up,ew+1,ns)-2.0d0*model%damage%sclr_damage(up,ew,ns)+model%damage%sclr_damage(up,ew-1,ns))/(dew**2.0d0)
            end if 

            ! Calculate y derivative
            if (ns==1) then
              yderiv = (model%damage%sclr_damage(up,ew,2)-model%damage%sclr_damage(up,ew,1))/dns
            elseif (ew == ewn) then
              yderiv = (model%damage%sclr_damage(up,ew,nsn)-model%damage%sclr_damage(up,ewn,ns-1))/dew
            else
              yderiv = (model%damage%sclr_damage(up,ew,ns+1)-2.0d0*model%damage%sclr_damage(up,ew,ns)+model%damage%sclr_damage(up,ew,ns-1))/(dew**2.0d0)
            end if 

            ! Calculate z derivative
            if (up==1) then
              zderiv = (model%damage%sclr_damage(2,ew,ns)-model%damage%sclr_damage(1,ew,ns))/dew
            elseif (ew == ewn) then
              zderiv = (model%damage%sclr_damage(upn,ew,ns)-model%damage%sclr_damage(upn-1,ewn,ns))/dew
            else
              zderiv = (model%damage%sclr_damage(up,ew,ns+1)-2.0d0*model%damage%sclr_damage(up,ew,ns)+model%damage%sclr_damage(up-1,ew,ns))/(1.0d0**2.0d0)
            end if 

            model%damage%sclr_damage(up,ew,ns) = -1.0d0*dt*uvel*xderiv - dt*vvel*yderiv-dt*wvel*zderiv + dt*damageSource(up,ew,ns) + model%damage%sclr_damage(up,ew,ns)

            if (model%damage%sclr_damage(up,ew,ns) >= 1.0d0) then
!              model%damage%sclr_damage = 1.0d0
            end if

            if (model%damage%sclr_damage(up,ew,ns) < 0.0d0) then
!              model%damage%sclr_damage = 0.0d0
            end if
          end do
        end if
      end do
    end do

    deallocate(damageSource)
  end subroutine update_damage


  ! To Do: 1. Use LAPACK     
  !        2. Calculate stress threshold
  !        3. Investigate Sigma coordinates
  function damage_source(tau, tau_xx, tau_yy, tau_xy, tau_yz, tau_xz, damage, up)

    use glimmer_physcon, only: rhoi, grav, pi
    use glimmer_global, only: dp

    implicit none

    !--------------------------------------------------
    ! Subroutine arguments
    !--------------------------------------------------
    real (kind = dp), intent(in) :: tau, tau_xx, tau_yy, tau_xy, tau_yz, tau_xz, damage
    integer, intent(in) :: up

    real (kind = dp) :: damage_source
    !--------------------------------------------------
    ! Internal variables
    !--------------------------------------------------
    real(dp), parameter :: alpha = 0.21d0 
    real(dp), parameter :: beta = 0.63d0
    real(dp), parameter :: B_parameter = 1.7d0 * 10.0d-9 
    real(dp), parameter :: r = 0.43d0
    real(dp), parameter :: stress_threshold = 3.3d5
    real(dp), parameter :: healing_parameter = 0.4d0        !lambda_h in Pralong05
    real(dp), parameter :: k1_parameter = 3.75d0 * 10.0d-3

    real(dp) :: max_principal_stress, invariant_I, invariant_II, invariant_III, dev_invariant_II
    real(dp) :: tau_zz, sigma_xx, sigma_yy, sigma_zz, hydrostatic_pressure	
    real(dp) :: psi, k_parameter, source, phi
    integer :: ewn, nsn, upn

    real(dp) :: Eig1, Eig2, Eig3, term1, traceAvg, determinant
    real(dp), dimension(3,3) :: tempMatrix 
    integer :: idx1, idx2

!    LAPACK eigsolver vars
!    real(dp) :: stress_tensor(3,3)
!    real(dp) :: work(6)
!    real(dp) :: principal_stress(3)

    ! Begin Calculations-------------------------------------------------------------

    ! Second deviatoric stress invariant
    dev_invariant_II = tau

    !-------Calculate Maximum Principal Stress----------------------------------------
    ! Max eigenvalue of Cauchy stress tensor = maximum principal stress

    ! Could use LAPACK to get max eigenvalue of cauchy stress te
    ! call DSYEV('N', 'U', 3, stress_tensor, 3, principal_stress, work, 6, 0)

    ! zz component of deviatoric stress tensor
    tau_zz = -1.0d0 * (tau_xx + tau_yy)

    !hydrostatic_pressure = mean stress
    ! hydro pressure = rho*g*(s-z) + tau_zz(z) = one-third of mean stress
    ! sigma_ij is ij'th component of cauchy stress tensor
    ! use sigma coords?
    hydrostatic_pressure = rhoi*grav*-1.0d0*up + tau_zz
    sigma_xx = tau_xx - hydrostatic_pressure
    sigma_yy = tau_yy - hydrostatic_pressure
    sigma_zz = tau_zz - hydrostatic_pressure

    ! Stress invariants
    invariant_I = sigma_xx + sigma_yy + sigma_zz 
    invariant_II = sigma_xx*sigma_yy + sigma_yy*sigma_zz + sigma_zz*sigma_xx - tau_xy**2.0d0 - tau_yz**2.0d0 - tau_xz**2.0d0

    invariant_III = sigma_xx*sigma_yy*sigma_zz
    invariant_III = invariant_III - sigma_xx*tau_yz**2.0d0
    invariant_III = invariant_III - sigma_yy*tau_xz**2.0d0
    invariant_III = invariant_III - sigma_zz*tau_xy**2.0d0
    invariant_III = invariant_III + 2.0d0*tau_xy*tau_yz*tau_xz

    ! Eigenvalues of the stress tensor
    ! Method due to: Oliver K. Smith: Eigenvalues of a symmetric 3 × 3 matrix. Commun. ACM 4(4): 168 (1961)
    traceAvg = (sigma_xx + sigma_yy + sigma_zz)/3.0d0

    tempMatrix(1,1) = sigma_xx-traceAvg
    tempMatrix(1,2) = tau_xy
    tempMatrix(1,3) = tau_xz
    tempMatrix(2,1) = tau_xy
    tempMatrix(2,2) = sigma_yy-traceAvg
    tempMatrix(2,3) = tau_yz
    tempMatrix(3,1) = tau_xz
    tempMatrix(3,2) = tau_yz
    tempMatrix(3,3) = sigma_zz-traceAvg

    determinant = 0.5d0*(tempMatrix(1,1)*tempMatrix(2,2)*tempMatrix(3,3) + tempMatrix(1,2)*tempMatrix(2,3)*tempMatrix(3,1) + tempMatrix(1,3)*tempMatrix(2,1)*tempMatrix(3,2)- tempMatrix(1,1)*tempMatrix(2,3)*tempMatrix(3,2) - tempMatrix(1,3)*tempMatrix(2,2)*tempMatrix(3,1))

    term1 = 0.0d0
    do idx1=1, 3
      do idx2 = 2, 3
        term1 = term1 + tempMatrix(idx1,idx2)**2.0d0
      end do
    end do
    term1 = term1 / 6.0d0

    phi = (1.0d0/3.0d0)*acos(determinant/(term1**(3.0d0/2.0d0)))

    if (abs(determinant) >= abs(term1**(3.0d0/2.0d0))) then
      phi = 0.0d0
      ! Raj_print *, 'determinant is', determinant, 'term1 is', term1**(1.5d0)
    end if

    if (phi < 0.0d0) then 
      phi = phi + pi/3.0d0
      ! Raj_print *, 'increasing phi'
    end if

    Eig1 = traceAvg + 2.0d0*sqrt(term1)*cos(phi)
    Eig2 = traceAvg - sqrt(term1)*(cos(phi) + sqrt(3.0d0)*sin(phi))
    Eig3 = traceAvg - sqrt(term1)*(cos(phi) - sqrt(3.0d0)*sin(phi))

    max_principal_stress = max(Eig1,Eig2,Eig3)
    !-------End calculation of maximum principal stress of Cauchy stress tensor------

    ! k_parameter calculation
    if (invariant_I >= 0) then
      k_parameter = k1_parameter * sqrt(invariant_I)
    else
      k_parameter = k1_parameter * (-1.0d0) * healing_parameter * sqrt(abs(invariant_I))
    end if

    ! FIX ME: Stress threshold calculation here....
    ! stress_threshold = 0.0d0

    ! Equation (25) in Pralong and Funk, 2005
    psi = 1.0d0 / (1.0d0 - damage)
    psi = (psi * (alpha*max_principal_stress + beta*sqrt(abs(3.0d0*dev_invariant_II)) + (1.0d0-alpha-beta)*invariant_I)) - stress_threshold

    ! Equation (24) in Pralong and Funk, 2005          
    if (psi >=  0.0d0) then
      source = B_parameter * (psi**r)*((1.0d0-damage)**(-1.0d0*k_parameter))
    else
      source = B_parameter*((-1.0d0*healing_parameter*psi)**r)*(1.0d0-damage)**(-1.0d0*k_parameter)
    end if

    ! Check for bad damage source term. 
    if (isnan(source)) then
      !print *, 'Bad damage source. Setting damage source to zero.' 
      source = 0.0d0
    end if

    damage_source = source
    return
  end function damage_source
end module glide_damage

  
