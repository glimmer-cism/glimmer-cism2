! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_types.f90 - part of the GLIMMER ice model         + 
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

module glide_types

  !*FD Holds type definitions for the derived types used by each 
  !*FD instance of the ice model. Originally, each of these types
  !*FD was a module containing variables, which were used as containers
  !*FD for global variables. However, the need to allow for multiple
  !*FD ice model instances meant that the nested derived types were instituted
  !*FD instead. However, there is probably one too many levels in this scheme. 
  !*FD It would be better if the different types here were contained in the 
  !*FD higher-level instance type (\texttt{glimmer\_params}), rather than 
  !*FD the intermediate model type (\texttt{glimmer\_global\_type}). 
  !*FD 
  !*FD Note that this \emph{is} now where the defaults are defined for these
  !*FD variables.
 
  use glimmer_sparse
  use glimmer_global
  use glimmer_ncdf
  use isostasy_types
  use lithot_types
  use profile
  use glimmer_coordinates, only : coordinates_type
  use glimmer_map_types, pi_dummy=>pi
  use glide_glenflow, only: glenflow_params
  use glimmer_deriv_time, only: timeders_type
  use glide_tempFullSoln, only: type_tempFullSoln
  use glide_thckADI, only: thckADI_type
  use velo_types
  use glimmer_slap, only: slapMatrix_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_general

    !*FD Holds fundamental parameters of the ice model geometry.

    integer :: ewn = 0  !*FD The number of grid-points in the E-W direction.
    integer :: nsn = 0  !*FD The number of grid-points in the N-S direction.
    integer :: upn = 1  !*FD The number of vertical levels in the model.

    real(sp), dimension(:),GC_DYNARRAY_ATTRIB :: x0 !original x0 grid 
    real(sp), dimension(:),GC_DYNARRAY_ATTRIB :: y0 !original y0 grid
    real(sp), dimension(:),GC_DYNARRAY_ATTRIB :: x1 !original x1 grid
    real(sp), dimension(:),GC_DYNARRAY_ATTRIB :: y1 !original y1 grid
  end type glide_general
     
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_options

    !*FD Holds user options controlling the methods used in the ice-model
    !*FD integration.

    integer :: whichtemp = 1

    !*FD Method of ice temperature calculation:
    !*FD \begin{description} 
    !*FD \item[0] Set column to surface air temperature
    !*FD \item[1] Do full temperature solution (also find vertical velocity
    !*FD and apparent vertical velocity)
    !*FD \end{description}

    integer :: whichflwa = 0

    !*FD Method for calculating flow factor $A$:
    !*FD \begin{description} 
    !*FD \item[0] \emph{Patterson and Budd} relationship 
    !*FD \item[1] \emph{Patterson and Budd} relationship, 
    !*FD with temperature set to $-10^{\circ}\mathrm{C}$ 
    !*FD \item[2] Set equal to $1\times 10^{-16}\,\mathrm{yr}^{-1}
    !*FD \,\mathrm{Pa}^{-n}$
    !*FD \end{description}

    integer :: whichbwat = 2

    !*FD Basal water depth: 
    !*FD \begin{description} 
    !*FD \item[0] Calculated from local basal water balance 
    !*FD \item[1] as {\bf 0}, including constant horizontal flow 
    !*FD \item[2] Set to zero everywhere 
    !*FD \end{description}

    integer :: whichmarn = 1

    !*FD Ice thickness: 
    !*FD \begin{description} 
    !*FD \item[0] No action 
    !*FD \item[1] Set thickness to zero if floating 
    !*FD \item[2] Set thickness to zero if relaxed bedrock is more 
    !*FD than certain water depth  
    !*FD \item[3] Lose fraction of ice when edge cell
    !*FD \end{description}

    integer :: whichbtrc = 0

    !*FD Basal slip coefficient: 
    !*FD \begin{description}
    !*FD \item[0] Set equal to zero everywhere
    !*FD \item[1] Set (non--zero) constant
    !*FD \item[2] Set to (non--zero) constant where where temperature is at pressure melting point of ice, otherwise to zero
    !*FD \item[3] \texttt{tanh} function of basal water depth 
    !*FD \end{description}

    integer :: whichevol = 0

    !*FD Thickness evolution method:
    !*FD \begin{description}
    !*FD \item[0] Pseudo-diffusion approach 
    !*FD \item[2] Diffusion approach (also calculates velocities) 
    !*FD \end{description}

    integer :: whichwvel = 0

    !*FD Vertical velocities: 
    !*FD \begin{description}
    !*FD \item[0] Usual vertical integration 
    !*FD \item[1] Vertical integration constrained so that 
    !*FD upper kinematic B.C. obeyed 
    !*FD \end{description}

    integer :: whichrelaxed = 0
    !*FD relaxed topography:
    !*FD \begin{description}
    !*FD \item[0] get relaxed topo from separate variable
    !*FD \item[1] first time slice of input topo is relaxed
    !*FD \item[2] first time slice of input topo is in isostatic equilibrium
    !*FD \end{description}

    integer :: hotstart = 0
    !*FD hotstart the model
    !*FD \begin{description}
    !*FD \item[0] normal start-up
    !*FD \item[1] hotstart model from previous run
    !*FD \end{description}

    integer :: periodic_ew = 0
    !*FD \begin{description}
    !*FD \item[0] no periodic EW boundary conditions
    !*FD \item[1] periodic EW boundary conditions
    !*FD \end{description}

    integer :: which_sigma = 0
    !*FD \begin{description}
    !*FD \item[0] calculate sigma coordinates
    !*FD \item[1] sigma coordinates are given in external file
    !*FD \item[2] sigma coordinates are given in configuration file
    !*FD \end{description}

    integer :: basal_mbal = 0
    !*FD \begin{description}
    !*FD \item[0] Basal melt rate not included in continuity equation
    !*FD \item[1] Basal melt rate included in continuity equation
    !*FD \end{description}

  end type glide_options

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geometry

    !*FD Holds fields and other information relating to the
    !*FD geometry of the ice sheet and bedrock.

    real, dimension(:,:), GC_DYNARRAY_ATTRIB :: temporary0 GC_DYNARRAY_INIT
    !*FD temporary array used for masking velocity grid
    real, dimension(:,:), GC_DYNARRAY_ATTRIB :: temporary1 GC_DYNARRAY_INIT
    !*FD temporary array used for masking temperature grid

    real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: thck GC_DYNARRAY_INIT
    !*FD The thickness of the ice, divided by \texttt{thk0}.

    real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: usrf GC_DYNARRAY_INIT
    !*FD The elevation of the upper ice surface, divided by \texttt{thk0}.

    real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: lsrf GC_DYNARRAY_INIT
    !*FD The elevation of the lower ice surface, divided by \texttt{thk0}.

    real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: topg GC_DYNARRAY_INIT
    !*FD The elevation of the topography, divided by \texttt{thk0}.

    integer, dimension(:,:),GC_DYNARRAY_ATTRIB :: thkmask GC_DYNARRAY_INIT
    !*FD see glimmer_mask.f90 for possible values

    real(dp) :: ivol, iarea !*FD ice volume and ice area

  end type glide_geometry

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geomderv

    !*FD Holds the horizontal and temporal derivatives of the thickness and
    !*FD upper surface elevation, as well as the thickness on the staggered grid.

    real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: dthckdew !*FD E-W derivative of thickness.
    real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: dusrfdew !*FD E-W derivative of upper surface elevation.
    real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: dthckdns !*FD N-S derivative of thickness.
    real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: dusrfdns !*FD N-S derivative of upper surface elevation.
    real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: stagthck !*FD Thickness averaged onto the staggered grid.

  end type glide_geomderv

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_velocity

    !*FD Holds the velocity fields in 2D and 3D. At least some of these fields
    !*FD are stored on the displaced grid.

    real(dp),dimension(:,:,:),GC_DYNARRAY_ATTRIB :: uvel  !*FD 3D $x$-velocity.
    real(dp),dimension(:,:,:),GC_DYNARRAY_ATTRIB :: vvel  !*FD 3D $y$-velocity.
    real(dp),dimension(:,:,:),GC_DYNARRAY_ATTRIB :: wvel  !*FD 3D $z$-velocity.
    real(dp),dimension(:,:,:),GC_DYNARRAY_ATTRIB :: wgrd  !*FD 3D grid vertical velocity.
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: uflx  !*FD 
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: vflx  !*FD 
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: diffu !*FD 
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: total_diffu !*FD total diffusivity
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: ubas  !*FD 
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: ubas_tavg 
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: vbas  !*FD 
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: vbas_tavg  
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: bed_softness !*FD bed softness parameter
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: btrc  !*FD  basal traction
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: tau_x !*FD basal shear stress, x-dir
    real(dp),dimension(:,:)  ,GC_DYNARRAY_ATTRIB :: tau_y !*FD basal shear stress, y-dir
  end type glide_velocity

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_climate
     !*FD Holds fields used to drive the model
     real(sp),dimension(:,:),GC_DYNARRAY_ATTRIB :: acab     !*FD Annual mass balance.
     real(sp),dimension(:,:),GC_DYNARRAY_ATTRIB :: acab_tavg     !*FD Annual mass balance (time average).
     real(sp),dimension(:,:),GC_DYNARRAY_ATTRIB :: artm     !*FD Annual mean air temperature
     real(sp),dimension(:,:),GC_DYNARRAY_ATTRIB :: lati     !*FD Latitudes of model grid points
     real(sp),dimension(:,:),GC_DYNARRAY_ATTRIB :: loni     !*FD Longitudes of model grid points
     real(sp),dimension(:,:),GC_DYNARRAY_ATTRIB :: calving  !*FD Calving flux (scaled as mass balance, thickness, etc)
     real(sp) :: eus = 0.                                  !*FD eustatic sea level
  end type glide_climate

  type glide_temper

    !*FD Holds fields relating to temperature.

    real(dp),dimension(:,:,:),GC_DYNARRAY_ATTRIB :: temp !*FD Three-dimensional temperature field.
    real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: bheatflx !*FD basal heat flux
    real(dp),dimension(:,:,:),GC_DYNARRAY_ATTRIB :: flwa !*FD Glenn's $A$.
    real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: bwat !*FD Basal water depth
    real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: stagbwat !*FD Basal water depth in velo grid
    real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: stagbtemp !*FD Basal temperature on velo grid
    real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: bmlt !*FD Basal melt-rate
    real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: bmlt_tavg !*FD Basal melt-rate
    real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: bpmp !*FD Basal pressure melting point
    real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: stagbpmp !*FD Basal pressure melting point on velo grid
    
    real(sp) :: perturb = 0.0    !*FD
    real(sp) :: grid    = 0.0    !*FD
    integer  :: tpt     = 0      !*FD Pointer to time series data
    logical  :: first1  = .true. !*FD
    logical  :: newtemps = .false. !*FD new temperatures
  end type glide_temper

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_funits
    character(fname_length) :: sigfile=''                      !*FD sigma coordinates file
    character(fname_length) :: ncfile=''                       !*FD configuration file for netCDF I/O
    type(glimmer_nc_output),pointer :: out_first=>NULL()       !*FD first element of linked list defining netCDF outputs
    type(glimmer_nc_input), pointer :: in_first=>NULL()        !*FD first element of linked list defining netCDF inputs
  end type glide_funits

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_numerics

    !*FD Parameters relating to the model numerics.
    real(sp) :: tstart = 0.0      !*FD starting time
    real(sp) :: tend   = 20000.0  !*FD end time
    real(sp) :: time   =    0.0   !*FD main loop counter in years
    real(sp) :: tinc   =   20.0   !*FD time step of main loop in years 
    real(sp) :: ntem   =    1.0   !*FD temperature time step (multiplier of main time step)
    real(sp) :: nvel   =    1.0   !*FD velocity time step (multiplier of main time step)
    real(dp) :: alpha  =    0.5d0 !*FD richard suggests 1.5 - was a parameter in original
    real(dp) :: alphas =    0.5d0 !*FD was a parameter in the original
    real(dp) :: thklim =  100.0   
    real(dp) :: mlimit = -200.0d0
    real(dp) :: calving_fraction = 0.8d0
    real(dp) :: dew    =   20.0d3
    real(dp) :: dns    =   20.0d3
    real(dp) :: dt     =    0.0
    real(dp) :: dttem  =    0.0
    real(sp) :: nshlf  =    0.0

    integer  :: timecounter = 0   !*FD count time steps
    
    ! Vertical coordinate ---------------------------------------------------
                                                               
    real(dp),dimension(:),GC_DYNARRAY_ATTRIB :: sigma !*FD Sigma values for 
                                                     !*FD vertical spacing of 
                                                     !*FD model levels
    integer :: profile_period = 100            !*FD profile frequency
  end type glide_numerics


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_pcgdwk
    real(dp),dimension(4)         :: fc      = 0.0
    real(dp),dimension(6)         :: fc2     = 0.0
  end type glide_pcgdwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_thckwk
     real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: oldthck  
     real(dp),dimension(:,:),  GC_DYNARRAY_ATTRIB :: oldthck2 
     real(dp),dimension(:,:),GC_DYNARRAY_ATTRIB :: float

  end type glide_thckwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_paramets
    real(dp),dimension(5) :: bpar = (/ 2.0d0, 10.0d0, 10.0d0, 0.0d0, 1.0d0 /)
    real(dp) :: btrac_const = 0.d0 ! m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: btrac_slope = 0.0d0 ! Pa^{-1} (gets scaled during init)
    real(dp) :: btrac_max = 0.d0  !  m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: fiddle = 3.0d0    ! -
    real(dp) :: hydtim = 1000.0d0 ! yr^{-1} converted to s^{-1} and scaled, 
                                  ! 0 if no drainage = 0.0d0 * tim0 / scyr
    real(dp) :: bwat_smooth = 0.01d0 ! basal water field smoothing strength
  end type glide_paramets

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_prof_type
     integer :: geomderv
     integer :: hvelos
     integer :: ice_mask1
     integer :: temperature
     integer :: ice_evo
     integer :: ice_mask2
     integer :: isos_water
     integer :: isos
  end type glide_prof_type

  type glide_global_type
    type(glide_general)  :: general
    type(coordinates_type) :: coordinates
    type(glide_options)  :: options
    type(glide_geometry) :: geometry
    type(glide_geomderv) :: geomderv
    type(glide_velocity) :: velocity
    type(glide_climate)  :: climate
    type(glide_temper)   :: temper
    type(lithot_type)    :: lithot
    type(glide_funits)   :: funits
    type(glide_numerics) :: numerics
    type(velo_type)      :: velowk
    type(glenflow_params) :: glenflow
    type(glide_pcgdwk)   :: pcgdwk
    type(glide_thckwk)   :: thckwk
    type(glide_paramets) :: paramets
    type(profile_type)   :: profile
    type(glide_prof_type) :: glide_prof
    type(isos_type)      :: isos
    type(timeders_type) :: timederivs
    type(type_tempFullSoln) :: tempFullSoln
    type(thckADI_type) :: thckADI
  end type glide_global_type

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_GLIDE_TYPES
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_GLIDE_TYPES
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_GLIDE_TYPES
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_GLIDE_TYPES
!MH!#endif
  
  subroutine glide_allocarr(model)
    
    !*FD Allocates the model arrays, and initialises some of them to zero.
    !*FD These are the arrays allocated, and their dimensions:
    !*FD
    !*FD In \texttt{model\%temper}:
    !*FD \begin{itemize}
    !*FD \item \texttt{temp(upn,0:ewn+1,0:nsn+1))}
    !*FD \item \texttt{bheatflx(ewn,nsn))}
    !*FD \item \texttt{flwa(upn,ewn,nsn))}
    !*FD \item \texttt{bwat(ewn,nsn))}
    !*FD \item \texttt{bmlt(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%velocity}:
    !*FD \begin{itemize}
    !*FD \item \texttt{uvel(upn,ewn-1,nsn-1))}
    !*FD \item \texttt{vvel(upn,ewn-1,nsn-1))}
    !*FD \item \texttt{wvel(upn,ewn,nsn))}
    !*FD \item \texttt{wgrd(upn,ewn,nsn))}
    !*FD \item \texttt{uflx(ewn-1,nsn-1))}
    !*FD \item \texttt{vflx(ewn-1,nsn-1))}
    !*FD \item \texttt{diffu(ewn,nsn))}
    !*FD \item \texttt{btrc(ewn,nsn))}
    !*FD \item \texttt{ubas(ewn,nsn))}
    !*FD \item \texttt{vbas(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%climate}:
    !*FD \begin{itemize}
    !*FD \item \texttt{acab(ewn,nsn))}
    !*FD \item \texttt{artm(ewn,nsn))}
    !*FD \item \texttt{lati(ewn,nsn))}
    !*FD \item \texttt{loni(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%geomderv}:
    !*FD \begin{itemize}
    !*FD \item \texttt{dthckdew(ewn,nsn))}
    !*FD \item \texttt{dusrfdew(ewn,nsn))}
    !*FD \item \texttt{dthckdns(ewn,nsn))}
    !*FD \item \texttt{dusrfdns(ewn,nsn))}
    !*FD \item \texttt{dthckdtm(ewn,nsn))}
    !*FD \item \texttt{dusrfdtm(ewn,nsn))}
    !*FD \item \texttt{stagthck(ewn-1,nsn-1))}
    !*FD \end{itemize}
  
    !*FD In \texttt{model\%geometry}:
    !*FD \begin{itemize}
    !*FD \item \texttt{thck(ewn,nsn))}
    !*FD \item \texttt{usrf(ewn,nsn))}
    !*FD \item \texttt{lsrf(ewn,nsn))}
    !*FD \item \texttt{topg(ewn,nsn))}
    !*FD \item \texttt{mask(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%thckwk}:
    !*FD \begin{itemize}
    !*FD \item \texttt{olds(ewn,nsn,thckwk\%nwhich))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%numerics}:
    !*FD \begin{itemize}
    !*FD \item \texttt{sigma(upn))}
    !*FD \end{itemize}

    use glimmer_log
    use glimmer_deriv_time, only: timeders_init
    use glimmer_horizCoord, only : horizCoord_allocate
    implicit none

    type(glide_global_type),intent(inout) :: model

    integer :: ewn,nsn,upn

    ! for simplicity, copy these values...

    ewn=model%general%ewn
    nsn=model%general%nsn
    upn=model%general%upn

    ! Allocate appropriately

    allocate(model%general%x0(ewn-1))!; model%general%x0 = 0.0
    allocate(model%general%y0(nsn-1))!; model%general%y0 = 0.0
    allocate(model%general%x1(ewn))!; model%general%x1 = 0.0
    allocate(model%general%y1(nsn))!; model%general%y1 = 0.0
    allocate(model%temper%temp(upn,0:ewn+1,0:nsn+1)); model%temper%temp = 0.0
    call horizCoord_allocate(model%coordinates%ice_grid, upn, model%temper%flwa)
    call horizCoord_allocate(model%coordinates%ice_grid, model%temper%bheatflx)
    call horizCoord_allocate(model%coordinates%ice_grid, model%temper%bwat)
    call horizCoord_allocate(model%coordinates%velo_grid, model%temper%stagbwat)
    call horizCoord_allocate(model%coordinates%velo_grid, model%temper%stagbtemp)
    call horizCoord_allocate(model%coordinates%ice_grid, model%temper%bmlt)
    call horizCoord_allocate(model%coordinates%ice_grid, model%temper%bpmp)
    call horizCoord_allocate(model%coordinates%ice_grid, model%temper%bmlt_tavg)
    call horizCoord_allocate(model%coordinates%velo_grid, model%temper%stagbpmp)

    call horizCoord_allocate(model%coordinates%velo_grid, upn, model%velocity%uvel)
    call horizCoord_allocate(model%coordinates%velo_grid, upn, model%velocity%vvel)
    call horizCoord_allocate(model%coordinates%ice_grid, upn, model%velocity%wvel)
    call horizCoord_allocate(model%coordinates%ice_grid, upn, model%velocity%wgrd)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%uflx)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%vflx)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%diffu)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%total_diffu)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%bed_softness)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%btrc)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%ubas)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%ubas_tavg)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%vbas)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%vbas_tavg)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%tau_x)
    call horizCoord_allocate(model%coordinates%velo_grid, model%velocity%tau_y)

    call horizCoord_allocate(model%coordinates%ice_grid, model%climate%acab)
    call horizCoord_allocate(model%coordinates%ice_grid, model%climate%acab_tavg)
    call horizCoord_allocate(model%coordinates%ice_grid, model%climate%artm)
    call horizCoord_allocate(model%coordinates%ice_grid, model%climate%lati)
    call horizCoord_allocate(model%coordinates%ice_grid, model%climate%loni)
    call horizCoord_allocate(model%coordinates%ice_grid, model%climate%calving)

    call horizCoord_allocate(model%coordinates%velo_grid, model%geomderv%dthckdew)
    call horizCoord_allocate(model%coordinates%velo_grid, model%geomderv%dusrfdew)
    call horizCoord_allocate(model%coordinates%velo_grid, model%geomderv%dthckdns)
    call horizCoord_allocate(model%coordinates%velo_grid, model%geomderv%dusrfdns)
    call horizCoord_allocate(model%coordinates%velo_grid, model%geomderv%stagthck)
  
    call horizCoord_allocate(model%coordinates%velo_grid, model%geometry%temporary0)
    call horizCoord_allocate(model%coordinates%ice_grid, model%geometry%temporary1)
    call horizCoord_allocate(model%coordinates%ice_grid, model%geometry%thck)
    call horizCoord_allocate(model%coordinates%ice_grid, model%geometry%usrf)
    call horizCoord_allocate(model%coordinates%ice_grid, model%geometry%lsrf)
    call horizCoord_allocate(model%coordinates%ice_grid, model%geometry%topg)
    call horizCoord_allocate(model%coordinates%ice_grid, model%geometry%thkmask)

    call timeders_init(model%timederivs,ewn,nsn)

    call horizCoord_allocate(model%coordinates%ice_grid, model%thckwk%oldthck)
    call horizCoord_allocate(model%coordinates%ice_grid, model%thckwk%oldthck2)
    call horizCoord_allocate(model%coordinates%ice_grid, model%thckwk%float)

    ! If we already have sigma, don't reallocate
    if (GC_DYNARRAY_CHECK(model%numerics%sigma)) then
       if (size(model%numerics%sigma)/=upn) then
          call write_log('Wrong number of sigma levels given',GM_FATAL)
       end if
    else
       allocate(model%numerics%sigma(upn))
    endif

    ! allocate isostasy grids
    call isos_allocate(model%isos,model%coordinates)
    ! allocate geothermal heat flux grids
    call lithot_allocate(model%lithot,model%coordinates)
    ! allocate velocity work grids
    call velo_allocate(model%velowk,model%coordinates)

  end subroutine glide_allocarr

  subroutine glide_deallocarr(model)
    !*FD deallocate model arrays

    use glimmer_deriv_time, only: timeders_final
    use glide_tempFullSoln, only: destroy_tempFullSoln
    implicit none
    type(glide_global_type),intent(inout) :: model

    deallocate(model%general%x0) 
    deallocate(model%general%y0) 
    deallocate(model%general%x1) 
    deallocate(model%general%y1) 

    deallocate(model%temper%temp)
    deallocate(model%temper%flwa)
    deallocate(model%temper%bheatflx)
    deallocate(model%temper%bwat)
    deallocate(model%temper%stagbwat)
    deallocate(model%temper%stagbtemp)
    deallocate(model%temper%bmlt)
    deallocate(model%temper%bmlt_tavg)
    deallocate(model%temper%bpmp)
    deallocate(model%temper%stagbpmp)

    deallocate(model%velocity%uvel)
    deallocate(model%velocity%vvel)
    deallocate(model%velocity%wvel)
    deallocate(model%velocity%wgrd)
    deallocate(model%velocity%uflx)
    deallocate(model%velocity%vflx)
    deallocate(model%velocity%diffu)
    deallocate(model%velocity%total_diffu)
    deallocate(model%velocity%bed_softness)
    deallocate(model%velocity%btrc)
    deallocate(model%velocity%ubas)
    deallocate(model%velocity%ubas_tavg)
    deallocate(model%velocity%vbas)
    deallocate(model%velocity%vbas_tavg)
    deallocate(model%velocity%tau_x)
    deallocate(model%velocity%tau_y)

    deallocate(model%climate%acab)
    deallocate(model%climate%acab_tavg)
    deallocate(model%climate%artm)
    deallocate(model%climate%lati)
    deallocate(model%climate%loni)

    deallocate(model%geomderv%dthckdew)
    deallocate(model%geomderv%dusrfdew)
    deallocate(model%geomderv%dthckdns)
    deallocate(model%geomderv%dusrfdns)
    deallocate(model%geomderv%stagthck)
  
    deallocate(model%geometry%temporary0)
    deallocate(model%geometry%temporary1)
    deallocate(model%geometry%thck)
    deallocate(model%geometry%usrf)
    deallocate(model%geometry%lsrf)
    deallocate(model%geometry%topg)
    deallocate(model%geometry%thkmask)

    call timeders_final(model%timederivs)
    deallocate(model%thckwk%oldthck)
    deallocate(model%thckwk%oldthck2)
    deallocate(model%thckwk%float)
    deallocate(model%numerics%sigma)
    
    ! deallocate isostasy grids
    call isos_deallocate(model%isos)
    ! deallocate geothermal heat flux grids
    call lithot_deallocate(model%lithot)
    ! deallocate temperature work arrays
    call destroy_tempFullSoln(model%tempFullSoln)
    ! deallocate velo work arrays
    call velo_deallocate(model%velowk)

  end subroutine glide_deallocarr

  ! some accessor functions
  function get_dew(model)
    !*FD return scaled x node spacing
    use glimmer_paramets, only : len0
    implicit none
    real(dp) :: get_dew
    type(glide_global_type) :: model

    get_dew = model%numerics%dew * len0
  end function get_dew

  function get_dns(model)
    !*FD return scaled y node spacing
    use glimmer_paramets, only : len0
    implicit none
    real(dp) :: get_dns
    type(glide_global_type) :: model

    get_dns = model%numerics%dns * len0
  end function get_dns

  function get_tstart(model)
    !*FD return start time
    implicit none
    real(sp) :: get_tstart
    type(glide_global_type) :: model
    
    get_tstart = model%numerics%tstart
  end function get_tstart

  function get_tend(model)
    !*FD return end time
    implicit none
    real(sp) :: get_tend
    type(glide_global_type) :: model
    
    get_tend = model%numerics%tend
  end function get_tend

  function get_tinc(model)
    !*FD return time increment
    implicit none
    real(sp) :: get_tinc
    type(glide_global_type) :: model
    
    get_tinc = model%numerics%tinc
  end function get_tinc

  function get_ewn(model)
    !*FD get number of nodes in x dir
    implicit none
    integer get_ewn
    type(glide_global_type) :: model

    get_ewn = model%general%ewn
  end function get_ewn

  function get_nsn(model)
    !*FD get number of nodes in y dir
    implicit none
    integer get_nsn
    type(glide_global_type) :: model

    get_nsn = model%general%nsn
  end function get_nsn
  
  subroutine set_time(model,time)
    !*FD Set the model time counter --- useful for
    !*FD fractional year output
    implicit none
    type(glide_global_type) :: model
    real :: time

    model%numerics%time=time
  end subroutine set_time

end module glide_types

