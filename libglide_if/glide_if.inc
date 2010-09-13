!/* -*- mode: f90 -*- */
!/* vim: set filetype=fortran : */

!This is the glide interface to which all implementations of the ice dynamics
!must conform.  For the moment it contains the minimum set of subroutines
!necessary for the simple drivers.  In the future, this file should contain a
!more complete API, with attendant documentation.  In other words, this API is
!just a test: it is not stable.

interface glide_config
  subroutine glide_config(model,config)
    use glide_types, only: glide_global_type
    use glimmer_config, only: ConfigSection
    type(glide_global_type)      :: model
    type(ConfigSection), pointer :: config
  end subroutine
end interface

interface glide_initialise
  subroutine glide_initialise(model)
    use glide_types, only: glide_global_type
    type(glide_global_type) :: model
  end subroutine
end interface

interface glide_tstep_p1
  subroutine glide_tstep_p1(model,time)
    use glide_types, only: glide_global_type
    use glimmer_global, only: rk
    type(glide_global_type) :: model
    real(rk), intent(in)    :: time
  end subroutine
end interface

interface glide_tstep_p2
  subroutine glide_tstep_p2(model)
    use glide_types, only: glide_global_type
    type(glide_global_type) :: model
  end subroutine
end interface

interface glide_tstep_p3
  subroutine glide_tstep_p3(model)
    use glide_types, only: glide_global_type
    type(glide_global_type) :: model
  end subroutine
end interface

interface glide_finalise
  subroutine glide_finalise(model,crash)
    use glide_types, only: glide_global_type
    type(glide_global_type) :: model
    logical,optional        :: crash
  end subroutine
end interface

interface glide_nc_fillall
  subroutine glide_nc_fillall(model)
    use glide_types, only: glide_global_type
    type(glide_global_type) :: model
  end subroutine
end interface

interface spinup_lithot
  subroutine spinup_lithot(model)
    use glide_types, only: glide_global_type
    type(glide_global_type) :: model
  end subroutine
end interface

