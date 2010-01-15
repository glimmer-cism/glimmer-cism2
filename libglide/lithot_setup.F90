! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  lithot_setup.f90 - part of the GLIMMER ice model         + 
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

!> routines for setting up bedrock temperature computations
!!
!! \author Magnus Hagdorn
!! \date January 2010

module lithot_setup

  use lithot_types

contains

  !> read lithot configuration
  subroutine lithot_readconfig(litho,config)
    use glimmer_config
    implicit none
    type(lithot_type) :: litho            !< structure holding bedrock temperature configuration
    type(ConfigSection), pointer :: config !< structure holding sections of configuration file

    ! local variables
    type(ConfigSection), pointer :: section

    ! read gthf section
    call GetSection(config,section,'GTHF')
    if (associated(section)) then
       litho%do_lithot = .True.
       call GetValue(section,'num_dim',litho%num_dim)
       call GetValue(section,'nlayer',litho%nlayer)
       call GetValue(section,'surft',litho%surft)
       call GetValue(section,'rock_base',litho%rock_base)
       call GetValue(section,'numt',litho%numt)
       call GetValue(section,'rho',litho%rho_r)
       call GetValue(section,'shc',litho%shc_r)
       call GetValue(section,'con',litho%con_r)
    end if
  end subroutine lithot_readconfig

  subroutine lithot_printconfig(litho)
    use glimmer_log
    implicit none
    type(lithot_type) :: litho            !< structure holding bedrock temperature configuration

    character(len=100) :: message
    
    if (litho%do_lithot) then
       call write_log('GTHF configuration')
       call write_log('------------------')
       if (litho%num_dim.eq.1) then
          call write_log('solve 1D diffusion equation')
       else
          write(message,*) 'Wrong number of dimensions.',litho%num_dim
          call write_log(message,GM_FATAL,__FILE__,__LINE__)
       end if
       write(message,*) 'number of layers                     : ',litho%nlayer
       call write_log(message)
       write(message,*) 'initial surface temperature          : ',litho%surft
       call write_log(message)
       write(message,*) 'rock base                            : ',litho%rock_base
       call write_log(message)
       write(message,*) 'density of rock layer                : ',litho%rho_r
       call write_log(message)
       write(message,*) 'specific heat capacity of rock layer : ',litho%shc_r
       call write_log(message)
       write(message,*) 'thermal conductivity of rock layer   : ',litho%con_r
       call write_log(message)
       write(message,*) 'number of time steps for spin-up     : ',litho%numt
       call write_log(message)
       call write_log('')
    end if
  end subroutine lithot_printconfig
end module lithot_setup
