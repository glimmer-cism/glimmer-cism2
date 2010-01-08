! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_ncio.f90 - part of the GLIMMER ice model         + 
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

#define NCO outfile%nc
#define NCI infile%nc

!> module for common netCDF I/O
!!
!! \author Magnus Hagdorn
!! \date 2004
module glimmer_ncio

  use glimmer_ncdf
  
contains
  !*****************************************************************************
  ! netCDF output
  !*****************************************************************************  
  !> open all netCDF files for output
  subroutine openall_out(out_list,coordinates,outfiles)
    use glimmer_ncdf
    use glimmer_coordinates
    implicit none
    type(glimmer_nc_output), pointer :: out_list !< first element of the output file linked list
    type(coordinates_type) :: coordinates       !< derived type holding glide coordinates
    type(glimmer_nc_output),pointer,optional :: outfiles
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>out_list
    end if

    do while(associated(oc))
       if (oc%append) then
          call glimmer_nc_openappend(oc)
       else
          call glimmer_nc_createfile(oc,coordinates)
       end if
       oc=>oc%next
    end do
  end subroutine openall_out

  !> close all netCDF files for output
  subroutine closeall_out(out_list,outfiles)
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_output), pointer :: out_list !< first element of the output file linked list
    type(glimmer_nc_output),pointer,optional :: outfiles

    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>out_list
    end if

    do while(associated(oc))
       oc=>delete(oc)
    end do
    if (.not.present(outfiles)) out_list=>NULL()
  end subroutine closeall_out

  !> open netCDF file for appending
  subroutine glimmer_nc_openappend(outfile)
    use glimmer_log
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_filenames
    use glimmer_global, only : sp
    implicit none
    type(glimmer_nc_output), pointer :: outfile !< structure containg output netCDF descriptor

    ! local variables
    integer :: status,timedimid,ntime,timeid
    real(sp),dimension(1) :: last_time
    character(len=msg_length) :: message

    ! open existing netCDF file
    status = nf90_open(process_path(NCO%filename),NF90_WRITE,NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
    call write_log_div
    write(message,*) 'Reopening file ',trim(process_path(NCO%filename)),' for output; '
    call write_log(trim(message))
    ! Find out when last time-slice was
    status = nf90_inq_dimid(NCO%id,'time',timedimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inquire_dimension(NCO%id,timedimid,len=ntime)
    call nc_errorhandle(__FILE__,__LINE__,status)
    ! Set timecounter
    outfile%timecounter=ntime+1
    write(message,*) '  Starting output at ',outfile%next_write,' and write every ',outfile%freq,' years'
    call write_log(trim(message))
    
    ! Get time varid
    status = nf90_inq_varid(NCO%id,'time',NCO%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Put dataset into define mode
    status = nf90_redef(NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)

  end subroutine glimmer_nc_openappend

  !> create a new netCDF file
  subroutine glimmer_nc_createfile(outfile,coordinates)
    use glimmer_log
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_filenames
    use glimmer_coordinates
    implicit none
    type(glimmer_nc_output), pointer :: outfile !< structure containg output netCDF descriptor
    type(coordinates_type) :: coordinates      !< derived type holding glide coordinates

    ! local variables
    integer status
    integer mapid
    character(len=msg_length) message

    ! create new netCDF file
    status = nf90_create(process_path(NCO%filename),NF90_CLOBBER,NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
    call write_log_div
    write(message,*) 'Opening file ',trim(process_path(NCO%filename)),' for output; '
    call write_log(trim(message))
    write(message,*) '  Starting output at ',outfile%next_write,' and write every ',outfile%freq,' years'
    call write_log(trim(message))
    if (outfile%end_write .lt. glimmer_nc_max_time) then
       write(message,*) '  Stop writing at ',outfile%end_write
       call write_log(trim(message))
    end if
    NCO%define_mode=.TRUE.

    ! writing meta data
    status = nf90_put_att(NCO%id, NF90_GLOBAL, 'Conventions', "CF-1.4")
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'title',trim(outfile%metadata%title))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'institution',trim(outfile%metadata%institution))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'source',trim(outfile%metadata%source))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'history',trim(outfile%metadata%history))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'references',trim(outfile%metadata%references))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'comment',trim(outfile%metadata%comment))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NF90_GLOBAL,'configuration',trim(outfile%metadata%config))
    call nc_errorhandle(__FILE__,__LINE__,status)
  
    ! defining time dimension and variable
    status = nf90_def_dim(NCO%id,'time',NF90_UNLIMITED,NCO%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    !     time -- Model time
    call write_log('Creating variable time')
    status = nf90_def_var(NCO%id,'time',NF90_FLOAT,(/NCO%timedim/),NCO%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(NCO%id, NCO%timevar, 'long_name', 'Model time')
    status = nf90_put_att(NCO%id, NCO%timevar, 'standard_name', 'time')
    status = nf90_put_att(NCO%id, NCO%timevar, 'units', 'year since 1-1-1 0:0:0')
    status = nf90_put_att(NCO%id, NCO%timevar, 'calendar', 'none')

    ! adding projection info
    if (glimmap_allocated(coordinates%projection)) then
       status = nf90_def_var(NCO%id,glimmer_nc_mapvarname,NF90_CHAR,mapid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       call glimmap_CFPutProj(NCO%id,mapid,coordinates%projection)
    end if

    ! setting the size of the level dimension
    NCO%nlevel = coordinates%sigma_grid%upn
  end subroutine glimmer_nc_createfile

  !> check if we should write to file
  subroutine glimmer_nc_checkwrite(outfile,forcewrite,time)
    use glimmer_log
    use glimmer_filenames
    use glimmer_global, only : sp
    implicit none
    type(glimmer_nc_output), pointer :: outfile !< data structure holding output file
    logical                 :: forcewrite       !< set to True when file should be written
    real(sp)                :: time             !< the next time at which output should be written

    character(len=msg_length) :: message
    integer status

    ! check if we are still in define mode and if so leave it
    if (NCO%define_mode) then
       status = nf90_enddef(NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
       NCO%define_mode = .FALSE.
    end if
    if (time.gt.NCO%processsed_time) then
       if (NCO%just_processed) then
          ! finished writing during last time step, need to increase counter...
          
          outfile%timecounter = outfile%timecounter + 1
          status = nf90_sync(NCO%id)
          call nc_errorhandle(__FILE__,__LINE__,status)
          NCO%just_processed = .FALSE.
       end if
    end if
    if (time.ge.outfile%next_write .or. (forcewrite.and.time.gt.outfile%next_write-outfile%freq)) then
       if (time.le.outfile%end_write .and. .not.NCO%just_processed) then
          call write_log_div
          write(message,*) 'Writing to file ', trim(process_path(NCO%filename)), ' at time ', time
          call write_log(trim(message))
          ! increase next_write
          outfile%next_write=outfile%next_write+outfile%freq
          NCO%processsed_time = time
          ! write time
          status = nf90_put_var(NCO%id,NCO%timevar,time,(/outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          NCO%just_processed = .TRUE.         
       end if
    end if
  end subroutine glimmer_nc_checkwrite

  !*****************************************************************************
  ! netCDF input
  !*****************************************************************************  
  !> open all netCDF files for input
  subroutine openall_in(in_list,coordinates)
    use glimmer_ncdf
    use glimmer_coordinates
    implicit none
    type(glimmer_nc_input), pointer :: in_list   !< first element of the input file linked list
    type(coordinates_type) :: coordinates       !< derived type holding glide coordinates
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ic=>in_list
    do while(associated(ic))
       call glimmer_nc_openfile(ic,coordinates)
       ic=>ic%next
    end do
  end subroutine openall_in

  !> close all netCDF files for input
  subroutine closeall_in(in_list)
    use glimmer_ncdf
    implicit none
    type(glimmer_nc_input), pointer :: in_list   !< first element of the input file linked list
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ic=>in_list
    do while(associated(ic))
       ic=>delete(ic)
    end do
    in_list=>NULL()
  end subroutine closeall_in

  !> open an existing netCDF file
  subroutine glimmer_nc_openfile(infile,coordinates)
    use glimmer_map_cfproj
    use glimmer_map_types
    use glimmer_log
    use glimmer_paramets, only: len0
    use glimmer_filenames
    use glimmer_coordinates
    implicit none
    type(glimmer_nc_input), pointer :: infile !< structure containg input netCDF descriptor
    type(coordinates_type) :: coordinates    !< derived type holding glide coordinates


    ! local variables
    integer dimsize, dimid, varid
    real, dimension(2) :: delta
    integer status    
    character(len=msg_length) message
    
    real,parameter :: small = 1.e-6

    ! open netCDF file
    status = nf90_open(process_path(NCI%filename),NF90_NOWRITE,NCI%id)
    if (status.ne.NF90_NOERR) then
       call write_log('Error opening file '//trim(process_path(NCI%filename))//': '//nf90_strerror(status),&
            type=GM_FATAL,file=__FILE__,line=__LINE__)
    end if
    call write_log_div
    call write_log('opening file '//trim(process_path(NCI%filename))//' for input')

    ! getting projection, if none defined already
    if (.not.glimmap_allocated(coordinates%projection)) coordinates%projection = glimmap_CFGetProj(NCI%id)

    ! getting time dimension
    status = nf90_inq_dimid(NCI%id, 'time', NCI%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    ! get id of time variable
    status = nf90_inq_varid(NCI%id,'time',NCI%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    
    ! getting length of time dimension and allocating memory for array containing times
    status = nf90_inquire_dimension(NCI%id,NCI%timedim,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    allocate(infile%times(dimsize))
    infile%nt=dimsize
    status = nf90_get_var(NCI%id,NCI%timevar,infile%times)

    ! setting the size of the level dimension
    NCI%nlevel = coordinates%sigma_grid%upn

    ! checking if dimensions and grid spacing are the same as in the configuration file
    ! x1
    status = nf90_inq_dimid(NCI%id,'x1',dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inquire_dimension(NCI%id,dimid,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize.ne.coordinates%ice_grid%size(1)) then
       write(message,*) 'Dimension x1 of file '//trim(process_path(NCI%filename))//' does not match with config dimension: ',&
            dimsize, coordinates%ice_grid%size(1)
       call write_log(message,type=GM_FATAL)
    end if
    status = nf90_inq_varid(NCI%id,'x1',varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_var(NCI%id,varid,delta)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (abs(delta(2)-delta(1) - coordinates%ice_grid%delta(1)*len0).gt.small) then
       write(message,*) 'deltax1 of file '//trim(process_path(NCI%filename))//' does not match with config deltax: ',&
            delta(2)-delta(1),coordinates%ice_grid%delta(1)*len0
       call write_log(message,type=GM_FATAL)
    end if

    ! y1
    status = nf90_inq_dimid(NCI%id,'y1',dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inquire_dimension(NCI%id,dimid,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize.ne.coordinates%ice_grid%size(2)) then
       write(message,*) 'Dimension y1 of file '//trim(process_path(NCI%filename))//' does not match with config dimension: ',&
            dimsize, coordinates%ice_grid%size(2)
       call write_log(message,type=GM_FATAL)
    end if
    status = nf90_inq_varid(NCI%id,'y1',varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_var(NCI%id,varid,delta)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (abs(delta(2)-delta(1) - coordinates%ice_grid%delta(2)*len0).gt.small) then
       write(message,*) 'deltay1 of file '//trim(process_path(NCI%filename))//' does not match with config deltay: ',&
            delta(2)-delta(1),coordinates%ice_grid%delta(2)*len0
       call write_log(message,type=GM_FATAL)
    end if
      
  ! Check that the number of vertical layers is the same, though it's asking for trouble
  ! to check whether the spacing is the same (don't want to put that burden on setup,
  ! plus f.p. compare has been known to cause problems here)
  status = nf90_inq_dimid(NCI%id,'level',dimid)
  ! If we couldn't find the 'level' dimension fail with a warning.
  ! We don't want to throw an error, as input files are only required to have it if they
  ! include 3D data fields.
  if (status == NF90_NOERR) then
        status = nf90_inquire_dimension(NCI%id, dimid, len=dimsize)
        call nc_errorhandle(__FILE__, __LINE__, status)
        if (dimsize.ne.coordinates%sigma_grid%upn .and. dimsize .ne. 1) then
            write(message,*) 'Dimension level of file '//trim(process_path(NCI%filename))//&
                ' does not match with config dimension: ', &
                dimsize, coordinates%sigma_grid%upn
            call write_log(message,type=GM_FATAL)
        end if
  else
        call write_log("Input file contained no level dimension.  This is not necessarily a problem.", type=GM_WARNING)
  end if
  
  end subroutine glimmer_nc_openfile

  !> check if we should read from file
  subroutine glimmer_nc_checkread(infile,time)
    use glimmer_log
    use glimmer_filenames
    use glimmer_global, only : sp
    implicit none
    type(glimmer_nc_input), pointer :: infile !< structure containg output netCDF descriptor
    real(sp) :: time                          !< time

    character(len=msg_length) :: message
    real(sp) :: sub_time

    if (infile%current_time.le.infile%nt) then
       if (.not.NCI%just_processed) then
          call write_log_div
          write(message,*) 'Reading time slice ',infile%current_time,'(',infile%times(infile%current_time),') from file ', &
               trim(process_path(NCI%filename)), ' at time ', time
          call write_log(message)
          NCI%just_processed = .TRUE.
          NCI%processsed_time = time
       end if
    end if
    if (time.gt.NCI%processsed_time) then
       if (NCI%just_processed) then
          ! finished reading during last time step, need to increase counter...
          infile%current_time = infile%current_time + 1
          NCI%just_processed = .FALSE.
       end if
    end if
  end subroutine glimmer_nc_checkread

end module glimmer_ncio

