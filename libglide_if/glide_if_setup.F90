! This file configures the glide interface, which at the moment means that it
! allows the glide library used to be specified in the configuration.
#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide_if_setup
contains
  subroutine glide_if_config(config)
    use glimmer_config
    implicit none

    interface glide_if_set_lib
      subroutine glide_if_set_lib(libname)
        character(len=*), intent(in) :: libname
      end subroutine
    end interface

    type(ConfigSection), pointer :: config
    type(ConfigSection), pointer :: section

    character(len=511) libname

    call GetSection(config,section,'glide_if')
    if (associated(section)) then
      call GetValue(section,'libname',libname)
      ! this subroutine is implemented in glide_if.c
      call glide_if_set_lib(libname)
    end if
  end subroutine
end module

