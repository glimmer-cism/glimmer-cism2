
subroutine glide_config(model,config)
  use glide, only: glide_config_default => glide_config
  use glide_types, only: glide_global_type
  use glimmer_config, only: ConfigSection
  implicit none
  type(glide_global_type)      :: model
  type(ConfigSection), pointer :: config

  print *, 'This library prints this line and calls the default glide_config'
  call glide_config_default(model,config)
end subroutine

subroutine glide_initialise(model)
  use glide, only: glide_initialise_default => glide_initialise
  use glide_types, only: glide_global_type
  implicit none
  type(glide_global_type) :: model

  print *, 'This library prints this line and calls the default glide_initialise'
  call glide_initialise_default(model)
end subroutine

subroutine glide_tstep_p1(model,time)
  use glide, only: glide_tstep_p1_default => glide_tstep_p1
  use glide_types, only: glide_global_type
  use glimmer_global, only: rk
  implicit none
  type(glide_global_type) :: model
  real(rk), intent(in)    :: time

  print *, 'This library prints this line and calls the default glide_tstep_p1'
  call glide_tstep_p1_default(model,time)
end subroutine

subroutine glide_tstep_p2(model)
  use glide, only: glide_tstep_p2_default => glide_tstep_p2
  use glide_types, only: glide_global_type
  implicit none
  type(glide_global_type) :: model

  print *, 'This library prints this line and calls the default glide_tstep_p2'
  call glide_tstep_p2_default(model)
end subroutine

subroutine glide_tstep_p3(model)
  use glide, only: glide_tstep_p3_default => glide_tstep_p3
  use glide_types, only: glide_global_type
  implicit none
  type(glide_global_type) :: model

  print *, 'This library prints this line and calls the default glide_tstep_p3'
  call glide_tstep_p3_default(model)
end subroutine

subroutine glide_finalise(model,crash)
  use glide_stop, only: glide_finalise_default => glide_finalise
  use glide_types, only: glide_global_type
  implicit none
  type(glide_global_type) :: model
  logical,optional        :: crash

  print *, 'This library prints this line and calls the default glide_finalise'
  call glide_finalise_default(model,crash)
end subroutine

subroutine glide_nc_fillall(model)
  use glide_nc_custom, only: glide_nc_fillall_default => glide_nc_fillall
  use glide_types, only: glide_global_type
  implicit none
  type(glide_global_type) :: model

  print *, 'This library prints this line and calls the default glide_nc_fillall'
  call glide_nc_fillall_default(model)
end subroutine

subroutine spinup_lithot(model)
  use glide_lithot, only: spinup_lithot_default => spinup_lithot
  use glide_types, only: glide_global_type
  implicit none
  type(glide_global_type) :: model

  print *, 'This library prints this line and calls the default spinup_lithot'
  call spinup_lithot_default(model)
end subroutine
