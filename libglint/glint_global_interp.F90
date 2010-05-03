! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_global_interp.f90 - part of the GLIMMER ice model  + 
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

module glint_global_interp

  use glint_global_grid
  use glimmer_global

!lipscomb - debug
    use glimmer_paramets, only: itest, jtest, jjtest, itest_local, jtest_local, stdout  

  implicit none

contains

  subroutine global_interp (in_grid,a,out_grid,ao,in_mask,out_mask,missing,error)

    ! This subroutine does an area weighted average from one grid,
    ! on a spherical earth, to another.  Logical masks may be assigned
    ! for each grid, and only those grid boxes which are masked true
    ! on both grids will be used.  A value of amm will be assigned
    ! to all nodes of the new grid which are initially false or have
    ! no data from the old grid. The new mask will also be changed to
    ! false where no data is available.
    !
    ! Restrictions:  longitude must be the first dimension and it
    !                be monotonically increasing (west to east).
    !
    !                latitude must be the second dimension and it
    !                must be monotonic.
    !
    !                values for longitude and latitude must be in
    !                degrees.
    !
    !                arrays that wrap around must repeat longitudes
    !                with a 360 degree increment.  it will be assumed
    !                that values in the wrapped input and mask arrays
    !                will also repeat (wrapped values in these arrays
    !                will not be used).
    !
    ! input
    !
    ! integer   idl    first dimension of input a and mask.
    ! integer   il     number of grid boxes in longitude for a and mask.
    ! real      alon   longitude (deg) limits of grid boxes for a and mask.
    ! integer   jl     number of grid boxes in latitude for a and mask.
    ! real      alat   latitude (deg) limits of grid boxes for a and mask.
    ! real      a      array of input data.
    ! logical   mask   mask for input data (.false. to mask out data).
    !
    ! output
    !
    ! integer   idlo   first dimension of output ao and masko.
    ! integer   ilo    number of grid boxes in longitude for ao and masko.
    ! real      alono  longitude (deg) limits of grid boxes for ao and masko.
    ! integer   jlo    number of grid boxes in latitude for ao and masko.
    ! real      alato  latitude (deg) limits of grid boxes for ao and masko.
    ! real      ao     array of output data.
    ! logical   masko  mask for output data (.false. to mask out data).
    ! integer   ier    error indication:
    !                  (values may be summed for multiple errors)
    !                  0  no errors
    !                  1  input longitude dimension and/or length <=0.
    !                  2  output dimension and/or length <=0.
    !                  4  input latititude dimension <=0.
    !                  8  output latitude dimension <=0.
    !                 16  wrap-around on input longitude grid doesn't
    !                     repeat (+360).
    !                 32  wrap-around on output longitude grid doesn't
    !                     repeat (+360).
    !                 64  longitude of input is not monotonic increasing.
    !                128  longitude of output is not monotonic increasing.
    !                256  latitude of input is not monotonic.
    !                512  latitude of output is not monotonic.
    !               1024  input longitude wraps but doesn't repeat identically.
    !               2048  output longitude wraps but doesn't repeat identically.
    !                 -1  output mask is changed.
    !                 -2  output mask contains all false values.

    ! --------------------------------------------------------
    ! Subroutine arguments
    ! --------------------------------------------------------

    type(global_grid)                              :: in_grid
    real(rk),dimension(:,:),         intent(in)    :: a
    type(global_grid)                              :: out_grid
    real(rk),dimension(:,:),         intent(out)   :: ao
    logical, dimension(:,:),optional,intent(inout) :: in_mask
    logical, dimension(:,:),optional,intent(inout) :: out_mask
    real(rk),               optional,intent(in)    :: missing
    integer,                optional,intent(out)   :: error

    ! --------------------------------------------------------
    ! Automatic arrays
    ! --------------------------------------------------------

    logical, dimension(size(a,1) ,size(a,2))  :: mask
    logical, dimension(size(ao,1),size(ao,2)) :: masko
    real(rk),dimension(size(a,1)+1)  :: alon
    real(rk),dimension(size(a,2)+1)  :: alat
    real(rk),dimension(size(ao,1)+1) :: alono
    real(rk),dimension(size(ao,2)+1) :: alato

    ! --------------------------------------------------------
    ! Internal variables
    ! --------------------------------------------------------

    integer  :: idl,il,jl,idlo,ilo,jlo,ier
    real(rk) :: api=3.1415926536
    real(rk) :: amm,almx,almn,sgn,al,dln,almxo,almno,dlno,amnlto
    real(rk) :: amxlto,amnlt,amxlt,amnlno,amxlno,amnln,amxln,wt,avg
    real(rk) :: slatmx,wlat,slatmn,slon,slonp,slonmx,slonmn,delon
    integer  :: i,j,iil,iilo,j1,j2,jj,i1,i2,k,ii,iii,iip

    ! --------------------------------------------------------
    ! Set up array sizes and check things match up.
    ! --------------------------------------------------------

    idl=size(a,1)
    il=size(alon)-1
    jl=size(alat)-1
    idlo=size(ao,1)
    ilo=size(alono)-1
    jlo=size(alato)-1

    alon=in_grid%lon_bound
    alat=in_grid%lat_bound
    alono=out_grid%lon_bound
    alato=out_grid%lat_bound

    ! Check array sizes --------------------------------------

    if (idl/=in_grid%nx.or. &
         jl/=in_grid%ny.or. &
         idlo/=out_grid%nx.or. &
         jlo/=out_grid%ny) then
       print*,'Array size mismatch in global_interp'
       stop
    end if

    ! Deal with optional mask input --------------------------

    if (present(in_mask)) then
       mask=in_mask
    else
       mask=.true.
    endif

    if (present(out_mask)) then
       masko=out_mask
    else
       masko=.true.
    endif

    ! Set up missing value -----------------------------------

    if (present(missing)) then
       amm=missing
    else
       amm=50.
    end if
    ier=0

    ! Check that the sizes of the arrays given are sensible --

    if (idl.lt.il.or.il.le.0) ier=1
    if (idlo.lt.ilo.or.ilo.le.0) ier=ier+2
    if (jl.le.0)  ier=ier+4
    if (jlo.le.0) ier=ier+8
    if (ier.gt.0) then
       if (present(error)) error=ier
       return
    end if

    ! Check monotonic increasing input longitudes ------------

    do i=2,il
       if (alon(i).le.alon(i-1)) then
          ier=ier+64
          exit
       endif
    end do

    ! Check monotonic increasing output longitudes -----------

    do i=2,ilo 
       if (alono(i).le.alono(i-1)) then
          ier=ier+128
          exit
       endif
    end do

    ! Check monotonicity of input latitudes ------------------

    sgn=(alat(2)-alat(1))
    do j=2,jl
       if (sgn.lt.0.0) then
          if (alat(j)-alat(j-1).ge.0) then
             ier=ier+256
             exit
          endif
       else if (sgn.gt.0.0) then
          if (alat(j)-alat(j-1).le.0.0) then
             ier=ier+256
             exit
          endif
       else
          ier=ier+256
          exit
       endif
    end do

    ! Check monotonicity of output latitudes ------------------

    sgn=(alato(2)-alato(1))
    do j=2,jlo 
       if (sgn.lt.0.0) then
          if (alato(j)-alato(j-1).ge.0.0) then
             ier=ier+512
             exit
          endif
       else if (sgn.gt.0.0) then
          if (alato(j)-alato(j-1).le.0.0) then
             ier=ier+512
             exit
          endif
       else
          ier=ier+512
          exit
       endif
    end do

    ! Find wrap around of input grid, if it exists ------------

    iil=il
    almx=alon(1)
    almn=alon(1)

    do i=2,il+1
       almx=max(almx,alon(i))
       almn=min(almn,alon(i))
       al=abs(alon(i)-alon(1))-360.0
       if (abs(al).le.1.e-4) then
          iil=i-1
          exit
       else if (al.gt.0.0) then
          ier=ier+1024
          go to 12
       endif
    end do

    dln=0.0
    if (almn.lt.0.0) then
       dln=int(-almn/360.0+.001)*360.0
    else if (almn.gt.360.0) then
       dln=-int(almn/360.0+.001)*360.0
    endif
12  continue

    ! Find wrap around of output grid, if it exists -----------

    iilo=ilo
    almxo=alono(1)
    almno=alono(1)
    do i=2,ilo+1
       almxo=max(almxo,alono(i))
       almno=min(almno,alono(i))
       al=abs(alono(i)-alono(1))-360.0
       if (abs(al).le.1.e-4) then
          iilo=i-1
          exit
       else if (al.gt.0.0) then
          ier=ier+2048
          go to 15
       endif
    end do

    dlno=0.0
    if (almno.lt.0.0) then
       dlno=int(-almno/360.0+.001)*360.0
    else if (almno.gt.360.0) then
       dlno=-int(almno/360.0+.001)*360.0
    endif
15  continue

    ! Test for errors.  return if any --------------------------

    if (ier.ne.0) then
       if (present(error)) error=ier
       return
    end if

    ! The output grid needs to begin with or after the input grid.

    if (almno+dlno.lt.almn+dln) dlno=dlno+360.0

    do j=1,jlo ! loop 200 - over output latitudes
       ! find index limits in latitude to cover the new grid.
       j1=jl+1
       j2=0
       amnlto=min(alato(j),alato(j+1))
       amxlto=max(alato(j),alato(j+1))

       ! search for index limits in j.

       do jj=1,jl
          amnlt=min(alat(jj),alat(jj+1))
          amxlt=max(alat(jj),alat(jj+1))
          ! find jj limits
          if (amxlt.gt.amnlto.and.amnlt.lt.amxlto) then
             j1=min(jj,j1)
             j2=max(jj,j2)
          endif
       end do

       ! if input grid doesn't at least partially cover the
       ! output grid box, no values will be assigned.  mask out
       ! all values for the latitude.

       if (j2.lt.j1) then
          do i=1,iilo 
             ao(i,j)=amm
             if (masko(i,j)) ier=-1
             masko(i,j)=.false.
          end do
          cycle
       endif

       do i=1,iilo ! loop 100

          ! no need to compute if it is masked out.
          if (.not.masko(i,j)) cycle

          ! find index limits in longitude to cover the new grid.
          i1=3*il+1
          i2=0
          amnlno=min(alono(i),alono(i+1))+dlno
          amxlno=max(alono(i),alono(i+1))+dlno

          ! search for index limits in i.
          ! because of wrap around it is necessary to
          ! look through the data twice.
          ! the output grid longitudes have been adjusted
          ! (using dlno) such that the first longitude in
          ! the output grid is greater than the first
          ! longitude on the input grid.

          do k=0,1
             do ii=1,iil 
                amnln=min(alon(ii),alon(ii+1))+dln+k*360.0
                amxln=max(alon(ii),alon(ii+1))+dln+k*360.0
                ! find ii limits
                if (amxln.gt.amnlno.and.amnln.lt.amxlno) then
                   i1=min(ii+k*il,i1)
                   i2=max(ii+k*il,i2)
                endif
             end do
          end do

          ! if input grid doesn't partially cover the output
          ! grid box, no values will be assigned.  mask out
          ! the grid box.

          if (i2.lt.i1) then
             ao(i,j)=amm
             if (masko(i,j)) ier=-1
             masko(i,j)=.false.
             cycle
          endif

          wt=0.0
          avg=0.0

          do jj=j1,j2
             slatmx=max(alat(jj),alat(jj+1))
             slatmn=min(alat(jj),alat(jj+1))
             wlat=max(sin(min(amxlto,slatmx)*api/180.)-sin(max(amnlto,slatmn)*api/180.),0.d0)
             if (wlat.ne.0.0) then
                do iii=i1,i2
                   slon=dln
                   slonp=dln
                   if (iii.gt.iil) then
                      slon=slon+360.
                      slonp=slonp+360.
                   endif
                   ii=mod(iii-1,iil)+1
                   iip=ii+1
                   if (mask(ii,jj)) then
                      slon=slon+alon(ii)
                      slonp=slonp+alon(iip)
                      slonmx=max(slon,slonp)
                      slonmn=min(slon,slonp)
                      delon=max(min(amxlno,slonmx)-max(amnlno,slonmn),0.d0)
                      wt=wt+wlat*delon
                      avg=avg+a(ii,jj)*wlat*delon
                   endif
                end do
             endif
          end do

          if (wt.gt.0.0) then
             ao(i,j)=avg/wt
          else
             ao(i,j)=amm
             if (masko(i,j)) ier=-1
             masko(i,j)=.false.
          endif
100       continue
       end do
200    continue
    end do

    ! Finish filling the output array from wrap-around.

    if (iilo.lt.ilo) then
       do j=1,jlo 
          do i=iilo+1,ilo
             ao(i,j)=ao(i-iilo,j)
             masko(i,j)=masko(i-iilo,j)
          end do
       end do
    endif

    ! Check if output masko is all false.

    if (all(.not.masko)) ier=-2

    ! Copy outputs if necessary

    if (present(error)) error=ier
    if (present(in_mask)) in_mask=mask
    if (present(out_mask)) out_mask=masko

  end subroutine global_interp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine interp_error(err,line)

    use glimmer_log

    integer :: err,line
    character(50) :: message
    integer :: e

    write(message,'(A,I6)')'Interpolation errors at line ',line
    call write_log(message)

    e=err

    call err_check(e,2048,'output longitude wraps but doesn''t repeat identically')
    call err_check(e,1024,'input longitude wraps but doesn''t repeat identically')
    call err_check(e,512, 'latitude of output is not monotonic')
    call err_check(e,256, 'latitude of input is not monotonic')
    call err_check(e,128, 'longitude of output is not monotonic increasing')
    call err_check(e,64,  'longitude of input is not monotonic increasing')
    call err_check(e,32,  'wrap-around on output longitude grid doesn''t repeat (+360)')
    call err_check(e,16,  'wrap-around on input longitude grid doesn''t repeat (+360)')
    call err_check(e,8,   'output latitude dimension <=0')
    call err_check(e,4,   'input latitude dimension <=0')
    call err_check(e,2,   'output dimension and/or length <=0')
    call err_check(e,1,   'input longitude dimension and/or length <=0')
    stop

  end subroutine interp_error

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine err_check(e,en,out)

    use glimmer_log

    integer :: e,en
    character(*) :: out

    if (e>en) then
       call write_log(out)
       e=e-en
    end if

  end subroutine err_check

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!lipscomb mod - three new subroutines
! The next three subroutines are used to upscale and downscale fields between a global 
! grid (with multiple elevation classes per gridcell) and the local Glimmer grid.
! For coupling to the Community Climate System Model (CCSM), the global grid is the 
! land surface grid, where the surface mass balance of ice sheets is computed for 
! multiple elevation classes.

 subroutine gcm_glint_downscaling (instance,            &
                                   tsfc_g,     qice_g,  &
                                   topo_g,     gmask)
 
    use glimmer_paramets, only: thk0

    use glint_type, only: glint_instance
    use glint_interp, only: interp_to_local


    ! Downscale fields from the global grid (with multiple elevation classes)
    !  to the local grid.

    type(glint_instance), intent(inout) :: instance
    real(dp),dimension(:,:,:),intent(in) :: tsfc_g       ! Surface temperature (C)
    real(dp),dimension(:,:,:),intent(in) :: qice_g       ! Surface mass balance (m)
    real(dp),dimension(:,:,:),intent(in) :: topo_g       ! Surface elevation (m)
    integer ,dimension(:,:),  intent(in),optional :: gmask ! = 1 where global data are valid
                                                           ! = 0 elsewhere

    real(dp), parameter :: maskval = 0.0_dp    ! value written to masked out gridcells

    integer ::       &
       nec,          &      ! number of elevation classes
       i, j, n,      &      ! indices 
       nxl, nyl             ! local grid dimensions
 
    real(dp), dimension(:,:,:), allocatable ::   &
       tsfc_l,    &! interpolation of global sfc temperature to local grid
       qice_l,    &! interpolation of global mass balance to local grid
       topo_l      ! interpolation of global topography in each elev class to local grid

    real(dp) :: fact, usrf

!lipscomb - to do - Lapse rate should be set to be consistent with CLM.
!                   Read from namelist?
    real(dp), parameter :: lapse = 0.0065_dp   ! atm lapse rate, deg/m

    nec = size(qice_g,3)
    nxl = instance%lgrid%size%pt(1)
    nyl = instance%lgrid%size%pt(2)

    allocate(tsfc_l(nxl,nyl,nec))
    allocate(topo_l(nxl,nyl,nec))
    allocate(qice_l(nxl,nyl,nec))

!   Downscale global fields for each elevation class to local grid
!   Set local fields to zero where interpolation from the global grid is invalid (instance%downs%lmask = 0).

!lipscomb - to do - Is topo_g downscaled correctly?

    if (present(gmask)) then   ! set local field = maskval where the global field is masked out

       do n = 1, nec
          call interp_to_local(instance%lgrid, tsfc_g(:,:,n), instance%downs, localdp=tsfc_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
          call interp_to_local(instance%lgrid, topo_g(:,:,n), instance%downs, localdp=topo_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
          call interp_to_local(instance%lgrid, qice_g(:,:,n), instance%downs, localdp=qice_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
       enddo

    else    ! global field values are assumed to be valid everywhere

       do n = 1, nec
          call interp_to_local(instance%lgrid, tsfc_g(:,:,n), instance%downs, localdp=tsfc_l(:,:,n))
          call interp_to_local(instance%lgrid, topo_g(:,:,n), instance%downs, localdp=topo_l(:,:,n))
          call interp_to_local(instance%lgrid, qice_g(:,:,n), instance%downs, localdp=qice_l(:,:,n))
       enddo

    endif

!lipscomb - debug - Write topo in each elevation class of global cell
!lipscomb - debug
    write (stdout,*) ' ' 
    write (stdout,*) 'Interpolate fields to local grid'
    write(stdout,*) 'Global cell =', itest, jjtest
    do n = 1, nec
       write(stdout,*) n, topo_g(itest,jjtest, n)
    enddo

    do j = 1, nyl
    do i = 1, nxl
        if ( (instance%downs%xloc(i,j,1) == itest .and. instance%downs%yloc(i,j,1) == jjtest) .or.  &
             (instance%downs%xloc(i,j,2) == itest .and. instance%downs%yloc(i,j,2) == jjtest) .or.  &
             (instance%downs%xloc(i,j,3) == itest .and. instance%downs%yloc(i,j,3) == jjtest) .or.  &
             (instance%downs%xloc(i,j,4) == itest .and. instance%downs%yloc(i,j,4) == jjtest) ) then
            write(stdout,*) i, j, thk0 * instance%model%geometry%usrf(i,j)
        endif
    enddo
    enddo
    
    i = itest_local
    j = jtest_local
    write (stdout,*) ' ' 
    write (stdout,*) 'Interpolated to local cells: i, j =', i, j
    do n = 1, nec
       write (stdout,*) ' '
       write (stdout,*) 'n =', n
       write (stdout,*) 'tsfc_l =', tsfc_l(i,j,n)
       write (stdout,*) 'topo_l =', topo_l(i,j,n)
       write (stdout,*) 'qice_l =', qice_l(i,j,n)
    enddo
!lipscomb - end debug

!   Interpolate tsfc and qice to local topography using values in the neighboring 
!    elevation classes.
!   If the local topography is outside the bounds of the global elevations classes,
!    extrapolate the temperature using the prescribed lapse rate.

    do j = 1, nyl
    do i = 1, nxl

       usrf = instance%model%geometry%usrf(i,j) * thk0   ! actual sfc elevation (m)

       if (usrf <= topo_l(i,j,1)) then
          instance%acab(i,j) = qice_l(i,j,1)
          instance%artm(i,j) = tsfc_l(i,j,1) + lapse*(topo_l(i,j,1)-usrf)
       elseif (usrf > topo_l(i,j,nec)) then
          instance%acab(i,j) = qice_l(i,j,nec)
          instance%artm(i,j) = tsfc_l(i,j,nec) - lapse*(usrf-topo_l(i,j,nec))
       else
          do n = 2, nec
             if (usrf > topo_l(i,j,n-1) .and. usrf <= topo_l(i,j,n)) then
                fact = (topo_l(i,j,n) - usrf) / (topo_l(i,j,n) - topo_l(i,j,n-1)) 
                instance%acab(i,j) = fact*qice_l(i,j,n-1) + (1._dp-fact)*qice_l(i,j,n)
                instance%artm(i,j) = fact*tsfc_l(i,j,n-1) + (1._dp-fact)*tsfc_l(i,j,n)
                exit
             endif
          enddo
       endif   ! usrf

!lipscomb - debug
          if (i==itest_local .and. j==jtest_local) then
             n = 4  
             write (stdout,*) ' '
             write (stdout,*) 'Interpolated values, i, j, n =', i, j, n
             write (stdout,*) 'usrf =', usrf
             write (stdout,*) 'acab =', instance%acab(i,j)
             write (stdout,*) 'artm =', instance%artm(i,j)
             write (stdout,*) 'topo(n-1) =', topo_l(i,j,n-1)
             write (stdout,*) 'topo(n) =', topo_l(i,j,n)
             write (stdout,*) 'qice(n-1) =', qice_l(i,j,n-1)
             write (stdout,*) 'qice(n) =', qice_l(i,j,n)
             write (stdout,*) 'tsfc(n-1) =', tsfc_l(i,j,n-1)
             write (stdout,*) 'tsfc(n) =', tsfc_l(i,j,n)
             write (stdout,*) 'fact = ', (topo_l(i,j,n) - usrf) / (topo_l(i,j,n) - topo_l(i,j,n-1)) 
          endif

    enddo  ! i
    enddo  ! j

  end subroutine gcm_glint_downscaling

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_gcm_upscaled_fields(instance,    nec,      &
                                     nxl,         nyl,      &
                                     nxg,         nyg,      &
                                     gfrac,       gthck,    &
                                     gtopo,       ghflx,    &
                                     groff)

    ! Upscale fields from the local grid to the global grid (with multiple elevation classes).
    ! This subroutine is modeled on get_i_upscaled_fields in glint_type.F90.

    use glint_type, only: glint_instance
    use glimmer_paramets, only: thk0
    use glimmer_log

    ! Arguments ----------------------------------------------------------------------------
 
    type(glint_instance),     intent(in)  :: instance      ! the model instance
    integer,                  intent(in)  :: nec           ! number of elevation classes
    integer,                  intent(in)  :: nxl,nyl       ! local grid dimensions 
    integer,                  intent(in)  :: nxg,nyg       ! global grid dimensions 

    real(dp),dimension(nxg,nyg,nec),intent(out) :: gfrac   ! ice-covered fraction [0,1]
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gthck   ! ice thickness (m)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gtopo   ! surface elevation (m)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: ghflx   ! heat flux (m)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: groff   ! runoff/calving flux (kg/m^2/s)
 
    ! Internal variables ----------------------------------------------------------------------
 
    real(dp),dimension(nxl,nyl) :: local_field
    real(dp),dimension(nxl,nyl) :: local_topo

    integer :: i, j            ! indices
 
    integer :: il, jl, ig, jg

!lipscomb - to do - Pass topomax as an argument
    real(dp), dimension(nec) :: topomax   ! upper elevation limit of each class

    ! Given the value of nec, specify the upper and lower elevation boundaries of each class.
    ! Note: These must be consistent with the values in the GCM.  Better to pass as an argument.
    if (nec == 1) then
       topomax = (/ 0._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 3) then
       topomax = (/ 0._dp,  1000._dp,  2000._dp, 10000._dp, 10000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 5) then
       topomax = (/ 0._dp,   500._dp,  1000._dp,  1500._dp,  2000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 10) then
       topomax = (/ 0._dp,   200._dp,   400._dp,   700._dp,  1000._dp,  1300._dp,  &
                            1600._dp,  2000._dp,  2500._dp,  3000._dp, 10000._dp /)
    else
       write(stdout,*) 'nec =', nec
       call write_log('ERROR: Current supported values of nec (no. of elevation classes) are 1, 3, 5, or 10', &
                       GM_FATAL,__FILE__,__LINE__)
    endif

    local_topo(:,:) = thk0 * instance%model%geometry%usrf(:,:)
    
!lipscomb - debug
       ig = itest
       jg = jjtest
       il = itest_local
       jl = jtest_local
       write(stdout,*) 'In get_glc_upscaled_fields'
       write(stdout,*) 'il, jl =', il, jl
       write(stdout,*) 'ig, jg =', ig, jg
       write(stdout,*) 'nxl, nyl =', nxl,nyl
       write(stdout,*) 'nxg, nyg =', nxg,nyg
       write(stdout,*) 'topo =', local_topo(il,jl) 
       call flush(stdout)

!lipscomb - to do - Is this mask needed?
    ! temporary field

    do j = 1, nyl
    do i = 1, nxl
       if (local_topo(i,j) > 0._dp) then
          local_field(i,j) = 1._dp
       else
          local_field(i,j) = 0._dp
       endif
    enddo
    enddo

!lipscomb - debug
       il = itest_local 
       jl = jtest_local 
       write(stdout,*) 'local ifrac =', local_field(il, jl)
       write(stdout,*) 'local topo =', local_topo(il,jl)
       write(stdout,*) 'local out_mask =', instance%out_mask(il,jl)

    ! ice fraction

    call mean_to_global_mec(instance%ups,                       &
                            nxl,                nyl,            &
                            nxg,                nyg,            &
                            nec,                topomax,        &
                            local_field,        gfrac,          &
                            local_topo,         instance%out_mask)

    ! ice thickness

    local_field(:,:) = thk0 * instance%model%geometry%thck(:,:)

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                nyl,        &
                            nxg,                nyg,        &
                            nec,                topomax,    &
                            local_field,        gthck,      &
                            local_topo,         instance%out_mask)


    ! surface elevation

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_topo,          gtopo,     &
                            local_topo,          instance%out_mask)

    ! heat flux

!lipscomb - to do - Copy runoff into local_field array
    local_field(:,:) = 0._dp

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         ghflx,     &
                            local_topo,          instance%out_mask)
 
!lipscomb - to do - Copy runoff into local_field array
    local_field(:,:) = 0._dp

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         groff,     &
                            local_topo,          instance%out_mask)

!lipscomb - debug
!       write(stdout,*) ' '
!       write(stdout,*) 'global ifrac:'
!       do n = 1, nec
!          write(stdout,*) n, gfrac(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global gtopo:'
!       do n = 1, nec
!          write(stdout,*) n, gtopo(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global gthck:'
!       do n = 1, nec
!          write(stdout,*) n, gthck(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global ghflx:'
!       do n = 1, nec
!          write(stdout,*) n, ghflx(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global groff:'
!       do n = 1, nec
!          write(stdout,*) n, groff(ig, jg, n)
!       enddo

  end subroutine get_gcm_upscaled_fields

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_to_global_mec(ups,                &
                                nxl,      nyl,      &
                                nxg,      nyg,      &
                                nec,      topomax,  &
                                local,    global,   &
                                ltopo,    mask)
 
    ! Upscale from the local domain to a global domain with multiple elevation classes
    ! by areal averaging.
    !
    ! This subroutine is adapted from subroutine mean_to_global in GLIMMER.
    ! The difference is that local topography is upscaled to multiple elevation classes
    !  in each global grid cell.
    !
    ! Note: This method is not the inverse of the interp_to_local routine.
    ! Also note that each local grid cell is assumed to have the same area.
    ! It would be better to have a more sophisticated routine.
 
    use glint_interp, only: upscale
    use glimmer_log

    ! Arguments
 
    type(upscale),            intent(in)    :: ups     ! upscaling indexing data
    integer,                  intent(in)    :: nxl,nyl ! local grid dimensions 
    integer,                  intent(in)    :: nxg,nyg ! global grid dimensions 
    integer,                  intent(in)    :: nec     ! number of elevation classes 
    real(dp),dimension(0:nec),intent(in)    :: topomax ! max elevation in each class 
    real(dp),dimension(nxl,nyl),  intent(in)      :: local   ! data on local grid
    real(dp),dimension(nxg,nyg,nec),intent(out)   :: global  ! data on global grid
    real(dp),dimension(nxl,nyl),  intent(in)      :: ltopo   ! surface elevation on local grid (m)
    integer, dimension(nxl,nyl),intent(in),optional :: mask ! mask for upscaling

    ! Internal variables
 
    integer ::  &
       i, j, n,    &! indices
       ig, jg       ! indices

    integer, dimension(nxl,nyl) ::  &
        tempmask,    &! temporary mask
        gboxec        ! elevation class associated with local topography

    integer, dimension(nxg,nyg,nec) ::  &
        gnumloc       ! no. of local cells within each global cell in each elevation class


    integer :: il, jl
    real(dp) :: lsum, gsum
 
    if (present(mask)) then
       tempmask(:,:) = mask(:,:)
    else
       tempmask(:,:) = 1
    endif
 
    ! Compute global elevation class for each local grid cell
    ! Also compute number of local cells within each global cell in each elevation class

    gboxec(:,:) = 0
    gnumloc(:,:,:) = 0

    do n = 1, nec
       do j = 1, nyl
       do i = 1, nxl
          if (ltopo(i,j) >= topomax(n-1) .and. ltopo(i,j) < topomax(n)) then
             gboxec(i,j) = n
             if (tempmask(i,j)==1) then
                ig = ups%gboxx(i,j)
                jg = ups%gboxy(i,j)
                gnumloc(ig,jg,n) = gnumloc(ig,jg,n) + 1
             endif
          endif
       enddo
       enddo
    enddo

    global(:,:,:) = 0._dp

    do j = 1, nyl
    do i = 1, nxl
       ig = ups%gboxx(i,j)
       jg = ups%gboxy(i,j)
       n = gboxec(i,j)
!lipscomb - bug check
       if (n==0) then
          write(stdout,*) 'Upscaling error: local topography out of bounds'
          write(stdout,*) 'i, j, topo:', i, j, ltopo(i,j)
          write(stdout,*) 'topomax(0) =', topomax(0)
          call write_log('Upscaling error: local topography out of bounds', &
               GM_FATAL,__FILE__,__LINE__)
       endif

!lipscomb - debug
       if (i==itest_local .and. j==jtest_local) then
          write(stdout,*) ' '
          write(stdout,*) 'il, jl =', i, j
          write(stdout,*) 'ig, jg, n =', ig, jg, n
          write(stdout,*) 'Old global val =', global(ig,jg,n)
          write(stdout,*) 'local, mask =', local(i,j), tempmask(i,j)
       endif

       global(ig,jg,n) = global(ig,jg,n) + local(i,j)*tempmask(i,j)

!lipscomb - debug
       if (i==itest_local .and. j==jtest_local) then
          write(stdout,*) 'New global val =', global(ig,jg,n)
       endif

    enddo
    enddo
 
    do n = 1, nec
       do j = 1, nyg
       do i = 1, nxg
          if (gnumloc(i,j,n) /= 0) then
             global(i,j,n) = global(i,j,n) / gnumloc(i,j,n)
          else
             global(i,j,n) = 0._dp
          endif
       enddo
       enddo
    enddo

    ! conservation check

    lsum = 0._dp
    do j = 1, nyl
    do i = 1, nxl
       lsum = lsum + local(i,j)*tempmask(i,j)
    enddo    
    enddo

    gsum = 0._dp
    do n = 1, nec
    do j = 1, nyg
    do i = 1, nxg
       gsum = gsum + global(i,j,n)*gnumloc(i,j,n)
    enddo
    enddo
    enddo

!lipscomb - to do - Use a less arbitrary error threshold
    if (abs(gsum-lsum) > 1.0_dp) then 
       write(stdout,*) 'local and global sums disagree'
       write (stdout,*) 'lsum, gsum =', lsum, gsum 
       call write_log('Upscaling error: local and glocal sums disagree', &
            GM_FATAL,__FILE__,__LINE__)
    endif

  end subroutine mean_to_global_mec
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glint_global_interp
