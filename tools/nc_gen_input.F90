program nc_gen_input

  use netcdf

  implicit none

  real,parameter :: pi = 3.1415926
  real,parameter :: g = 9.81
  !real,parameter :: rhoi_rhow = 1030.d0 / 920.d0

  character(len=512) :: gen_usage = "nc_gen_input <type_code> &
       &nc_filename [params] &
       &<type_code> = [is|ls]"

  character(len=512) :: conf_shelf_params = "<fname> <nx> <ny> <hx> <hy> &
       &<slope_start_pos> <grounded_ice_thk> &
       &<ice_front_pos> <ice_front_thick> &
       &<ocean_topg> <land_topg> <k_x> & 
       &<chan_amp> <chan_init_length> & 
       &<kinbc_width> <nsswitch>"

  character(len=512) :: ghost_shelf_params = "<fname> <nx> <ny> <kinbcw> <n_level> <hx> <hy> &
       &<upstream_thk> <upstream_vel> <front_thk> &
       &<k_x> <chan_amp> <ramp_width> <landw> <ifpos>"


  character (len=2) :: type_code

  !variables to be used in all cases and in netcdf writing
  character(len=1024) :: argstr
  character (len=2056) :: fname
  integer :: nx,ny,n_level

  !data arrays
  real,dimension(:),allocatable :: xs,ys,xstag,ystag,level
  real,dimension(:,:),allocatable :: thck,topog
  real,dimension(:,:),allocatable :: kinbcmask
  real,dimension(:,:,:),allocatable :: uvelhom,vvelhom

  !get arguments from command line

  if (command_argument_count() < 1) then

     write(*,*) "Usage: ", trim(gen_usage)
     stop 1

  end if

  call get_command_argument(1,argstr)

  if (trim(argstr) == '-h' .or. trim(argstr) == '--help') then
     write(*,*) "Usage: ", trim(gen_usage)
     stop
  end if

  type_code = argstr

  if (type_code == 'is') then

     call make_ice_stream()

  else if (type_code == 'ls') then
  
     call make_linear_shelf()

  else

     write(*,*) "Unrecognized type_code: ", type_code

  end if

contains

  subroutine make_linear_shelf()

    ! local variables
    integer :: i,j,k
    integer :: ifpos,kinbcw
    integer,parameter :: zero_buf = 1
    real :: hx,hy
    real :: rhoi, rhoo
    real :: upstream_thk, if_thk, upstream_vel,inflow_a
    real :: chan_depth,otopg,ltopg,acab_per_year
    real,dimension(:),allocatable :: rand_row
    logical :: noslip
    integer,parameter :: s_marg = 2, w_marg = 1

    real :: tmp
    real :: rand_amp

    character (len=512) :: linear_shelf_params = "<fname> <nx> <ny> <n_level> <hx> <hy> <upstream_thk> &
       & <upstream_vel> <inflow_a> <ifpos> <ifthk> <ocean_depth> <kinbcw> [<noslip>]"

    if (command_argument_count() < 14) then
       write(*,*)"Incorrect number of parameters. Linear shelf requires: &
            &  ",trim(linear_shelf_params)
       stop 1
    end if

    call get_command_argument(2,argstr)
    read(argstr,'(a512)') fname
    write(*,*) 'fname ',trim(fname)

    call get_command_argument(3, argstr)
    read(argstr,'(i5)') nx
    write(*,*) 'nx:',nx

    call get_command_argument(4,argstr)
    read(argstr,'(i5)') ny
    write(*,*) 'ny:',ny

    call get_command_argument(5,argstr)
    read(argstr,'(i5)') n_level
    write(*,*) 'n_level:',n_level

    call get_command_argument(6,argstr)
    read(argstr,'(f18.12)') hx
    write(*,*) 'hx',hx

    call get_command_argument(7,argstr)
    read(argstr,'(f18.12)') hy
    write(*,*) 'hy',hy

    call get_command_argument(8,argstr)
    read(argstr,'(f18.12)') upstream_thk
    write(*,*) 'upstream_thk',upstream_thk

    call get_command_argument(9,argstr)
    read(argstr,'(f18.12)') upstream_vel
    write(*,*) 'upstream_vel',upstream_vel

    call get_command_argument(10,argstr)
    read(argstr,'(f18.12)') inflow_a
    write(*,*) 'inflow_a', inflow_a

    call get_command_argument(11,argstr)
    read(argstr,'(i5)') ifpos
    write(*,*) 'ifpos',ifpos

    call get_command_argument(12,argstr)
    read(argstr,'(f18.12)') if_thk
    write(*,*) 'if_thk', if_thk

    call get_command_argument(13,argstr)
    read(argstr,'(f18.12)') otopg
    write(*,*) 'otopg', otopg

    call get_command_argument(14,argstr)
    read(argstr,'(i5)') kinbcw
    write(*,*) 'kinbc_width',kinbcw

    if (command_argument_count() == 15) then
       call get_command_argument(15,argstr)
       read(argstr, '(l1)') noslip
    else
	noslip = .false.
    end if
    write(*,*) 'noslip', noslip

    allocate(xs(nx),ys(ny),level(n_level),xstag(nx-1),ystag(ny-1))
    allocate(topog(nx,ny),thck(nx,ny),kinbcmask(nx-1,ny-1))
    allocate(uvelhom(nx-1,ny-1,n_level), vvelhom(nx-1,ny-1,n_level))
    allocate(rand_row(nx-2*zero_buf))

    !now populate the dimension variables
    xs = (/ ( (real(i-w_marg)-0.5d0)*hx,i=1,nx ) /)
    ys = (/ ( (real(j-s_marg)-0.5d0)*hy,j=1,ny ) /)
    level = (/ ( real(i)/real((n_level-1)), i=0,(n_level-1) ) /)
    xstag = (/ ( (real(i-w_marg)*hx),i=1,nx-1 ) /)
    ystag = (/ ( (real(j-s_marg)*hy),j=1,ny-1 ) /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! define the topography, thickness, kinbcmask       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !define topg
    topog = -abs(otopg)
    
    !define thickness
    thck = 0.0
!    thck(1+zero_buf:(nx-zero_buf),(ny-kinbcw):ny) = upstream_thk 
    thck(1+zero_buf:(nx-zero_buf),1:kinbcw) = upstream_thk 
    
!    do j=ifpos,(ny-kinbcw)
    do j=(kinbcw+1),ifpos
       thck(1+zero_buf:(nx-zero_buf),j) = if_thk +  &
!              real(j-ifpos)/(ny-kinbcw+1-ifpos) * (upstream_thk-if_thk)
            real(ifpos+0.5d0-j)/(ifpos+0.5d0-kinbcw-0.5d0) * (upstream_thk-if_thk)
    end do

    ! define kinbcmask
    kinbcmask = 0

!    kinbcmask(:,(ny-6):(ny-1)) = 1 !north edge

    kinbcmask(:,1:kinbcw) = 1 !south edge
!    if (noslip) then		    
       kinbcmask(1:(1+zero_buf),:) = 1
       kinbcmask((nx-1-zero_buf):(nx-1),:) = 1
!    end if

    uvelhom = 0.d0
    vvelhom = 0.d0

    if (noslip) then
       !vvelhom(2:(nx-2),(ny-kinbcw):ny,:) = upstream_vel
       !vvelhom(:,(ny-kinbcw):(ny-1),:) = upstream_vel
	do i=2,(nx-2)
	    vvelhom(i,1:kinbcw,:) = upstream_vel * \
               real((nx-2-i)*(i-2))/real((nx/2.d0-2.d0)*(nx/2.d0-2.d0))
	end do
    else
       !vvelhom(:,ny-kinbcw:ny-1,:) = upstream_vel
       !vvelhom(:,1:kinbcw,:) = upstream_vel
       vvelhom(:,:,:) = upstream_vel
    end if

    !do j=ny-kinbcw-10,ny
    do j=2,min(10,ny)
       do i=1+zero_buf,nx-zero_buf
!          thck(i,j) = thck(i,j) + real(j-(ny-kinbcw-10))/real(ny-(ny-kinbcw-10)) & 
          thck(i,j) = thck(i,j) + real(10-j)/real(10-2) & 
                                  *inflow_a * (xs(i)-xs((nx-1)/2+1))*(xs(i)-xs((nx-1)/2+1))
       end do
    end do
    thck(:,1) = thck(:,2) 

    call write_nc_file(.true.)

    deallocate(level,xs,ys,xstag,ystag)
    deallocate(thck,topog,kinbcmask,uvelhom,vvelhom)

  end subroutine make_linear_shelf

  subroutine make_ice_stream()

    ! local variables
    integer :: i,j,k
    real :: hx,hy
    real :: rhoi, rhoo
    real :: usurf
    integer,parameter :: s_marg = 2, w_marg = 1

    real :: topog_amp
    real :: odepth
    real :: bed_depth
    real :: int_depth

    real :: surf_1 
    real :: surf_2 
    real :: surf_3
    real :: thk_final 

    real :: inflow_vel

    character (len=512) :: linear_shelf_params = "<fname> <nx> <ny> <n_level> <hx> <hy> <topog_amp> &
	                                        &<odepth> <bed_depth> <int_depth> <surf1> <surf2> <surf3> &
	                                        &<thk_final> <inflow_vel>"


    if (command_argument_count() < 16) then
       write(*,*)"Incorrect number of parameters. Linear shelf requires: &
            &  ",trim(linear_shelf_params)
       stop 1
    end if

    call get_command_argument(2,argstr)
    read(argstr,'(a512)') fname
    write(*,*) 'fname ',trim(fname)

    call get_command_argument(3, argstr)
    read(argstr,'(i5)') nx
    write(*,*) 'nx:',nx

    call get_command_argument(4,argstr)
    read(argstr,'(i5)') ny
    write(*,*) 'ny:',ny

    call get_command_argument(5,argstr)
    read(argstr,'(i5)') n_level
    write(*,*) 'n_level:',n_level

    call get_command_argument(6,argstr)
    read(argstr,'(f18.12)') hx
    write(*,*) 'hx',hx

    call get_command_argument(7,argstr)
    read(argstr,'(f18.12)') hy
    write(*,*) 'hy',hy

    call get_command_argument(8,argstr)
    read(argstr,'(f18.12)') topog_amp
    write(*,*) 'topog_amp', topog_amp
    
    call get_command_argument(9,argstr)
    read(argstr,'(f18.12)') odepth
    write(*,*) 'odepth', odepth

    call get_command_argument(10,argstr)
    read(argstr,'(f18.12)') bed_depth
    write(*,*) 'bed_depth', bed_depth

    call get_command_argument(11,argstr)
    read(argstr,'(f18.12)') int_depth
    write(*,*) 'int_depth', int_depth

    call get_command_argument(12,argstr)
    read(argstr,'(f18.12)') surf_1
    write(*,*) 'surf_1', surf_1

    call get_command_argument(13,argstr)
    read(argstr,'(f18.12)') surf_2
    write(*,*) 'surf_2', surf_2

    call get_command_argument(14,argstr)
    read(argstr,'(f18.12)') surf_3
    write(*,*) 'surf_3', surf_3

    call get_command_argument(15,argstr)
    read(argstr,'(f18.12)') thk_final
    write(*,*) 'thk_final', thk_final

    call get_command_argument(16,argstr)
    read(argstr,'(f18.12)') inflow_vel
    write(*,*) 'inflow_vel', inflow_vel

    allocate(xs(nx),ys(ny),level(n_level),xstag(nx-1),ystag(ny-1))
    allocate(topog(nx,ny),thck(nx,ny),kinbcmask(nx-1,ny-1))
    allocate(uvelhom(nx-1,ny-1,n_level), vvelhom(nx-1,ny-1,n_level))

    !now populate the dimension variables
    xs = (/ ( (real(i-w_marg)-0.5d0)*hx,i=1,nx ) /)
    ys = (/ ( (real(j-s_marg)-0.5d0)*hy,j=1,ny ) /)
    level = (/ ( real(i)/real((n_level-1)), i=0,(n_level-1) ) /)
    xstag = (/ ( (real(i-w_marg)*hx),i=1,nx-1 ) /)
    ystag = (/ ( (real(j-s_marg)*hy),j=1,ny-1 ) /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! define the topography, thickness, kinbcmask       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !define topg
    do i=1,nx
	do j=1,ny
	   topog(i,j) = (((2.0*(i-1)/(nx-1))-1.0)**4.d0) *topog_amp
	   if (j >= 3*ny/4) then
              topog(i,j) = topog(i,j) - odepth
	   elseif (j >= 2*ny/4) then
              topog(i,j) = topog(i,j) - odepth - (3.d0*ny/4.d0-j)/(ny/4.d0)*(bed_depth-odepth)
           else
              topog(i,j) = topog(i,j) - bed_depth + (2.d0*ny/4.d0-j)/(ny/2.d0)*(bed_depth-int_depth)
           end if

    	end do
    end do

    !define thickness
    thck = 0.d0

    do i=1,nx
	do j=1,int(ny/2.0)
           thck(i,j) = surf_1 - topog(i,j) + (j-1)/(ny/2.0-1)*(surf_2-surf_1)
        end do
        do j=int(ny/2.0)+1,int(3.0*ny/4.0)
           thck(i,j) = (surf_2-(j-ny/2.0)/(ny/4.0)*(surf_2-surf_3)) - topog(i,j)
        end do
        do j=int(3.0*ny/4.0)+1,ny-4
	   usurf = surf_3-(j-3.0*ny/4.0)/(ny-4.0-3.0*ny/4.0)*(surf_3 - 0.1*thk_final)
           thck(i,j) = (surf_3-topog(i,3*ny/4)) - (j-3.0*ny/4.0)/(ny-4.0-3.0*ny/4.0)*(surf_3-topog(i,3*ny/4)-thk_final)
	   if (topog(i,j)+thck(i,j) > usurf) then
	      thck(i,j) = max(0.0,usurf-topog(i,j))
	   end if
        end do
    end do

    thck = max(0.0,thck)

    ! define kinbcmask
    kinbcmask = 0
    kinbcmask(:,1:2) = 1

!    kinbcmask(:,(ny-10):ny-1) = 1
!    kinbcmask(1:10,:) = 1
!    kinbcmask(nx-10:nx,:) = 1

    vvelhom = inflow_vel
    uvelhom = 0.0

    call write_nc_file(.true.)

    deallocate(level,xs,ys,xstag,ystag)
    deallocate(thck,topog,kinbcmask,uvelhom,vvelhom)

  end subroutine make_ice_stream

  subroutine write_nc_file(write_vel)

    logical,intent(in) :: write_vel

    !local variables
    integer :: nc_id
    integer :: time_dimid,x_dimid,y_dimid,xstag_dimid,ystag_dimid,level_dimid
    integer :: x_varid,y_varid,time_varid,level_varid
    integer :: thck_varid,topog_varid,kinbcmask_varid,uvel_varid,vvel_varid
    integer :: xstag_varid,ystag_varid

    call check( nf90_create(fname, NF90_CLOBBER, nc_id) )

    ! defining dimensions
    call check( nf90_def_dim(nc_id,'time',1,time_dimid) )
    call check( nf90_def_dim(nc_id,'x1',nx,x_dimid) )
    call check( nf90_def_dim(nc_id,'y1',ny,y_dimid) )
    call check( nf90_def_dim(nc_id,'x0',nx-1,xstag_dimid) )
    call check( nf90_def_dim(nc_id,'y0',ny-1,ystag_dimid) )


    ! define variables
    call check( nf90_def_var(nc_id,'time',NF90_DOUBLE,(/time_dimid/),time_varid) )
    call check( nf90_put_att(nc_id, time_varid, 'long_name', 'time') )
    call check( nf90_put_att(nc_id, time_varid, 'units', 'seconds') )

    call check( nf90_def_var(nc_id,'x1',NF90_DOUBLE,(/x_dimid/),x_varid) )
    call check( nf90_put_att(nc_id, x_varid, 'long_name', 'Cartisian x-coordinate') )
    call check( nf90_put_att(nc_id, x_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'y1',NF90_DOUBLE,(/y_dimid/),y_varid) )
    call check( nf90_put_att(nc_id, y_varid, 'long_name', 'Cartisian y-coordinate') )
    call check( nf90_put_att(nc_id, y_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'x0',NF90_DOUBLE,(/xstag_dimid/),xstag_varid) )
    call check( nf90_put_att(nc_id, xstag_varid, 'long_name', 'Cartisian y-coordinate velocity grid') )
    call check( nf90_put_att(nc_id, xstag_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'y0',NF90_DOUBLE,(/ystag_dimid/),ystag_varid) )
    call check( nf90_put_att(nc_id, ystag_varid, 'long_name', 'Cartisian y-coordinate velocity grid') )
    call check( nf90_put_att(nc_id, ystag_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'thk',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),thck_varid) )
    call check( nf90_put_att(nc_id, thck_varid, 'long_name', 'ice thickness') )
    call check( nf90_put_att(nc_id, thck_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'topg',NF90_DOUBLE,(/x_dimid,y_dimid/),topog_varid) )
    call check( nf90_put_att(nc_id, topog_varid, 'long_name', 'topography') )
    call check( nf90_put_att(nc_id, topog_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'kinbcmask',NF90_DOUBLE,(/xstag_dimid,ystag_dimid,time_dimid/),kinbcmask_varid) )
    call check( nf90_put_att(nc_id, kinbcmask_varid, 'long_name', 'kinematic boundary condition mask') )

    if (write_vel) then
       call check( nf90_def_dim(nc_id,'level',n_level,level_dimid) )
       call check( nf90_def_var(nc_id,'level',NF90_DOUBLE,(/level_dimid/),level_varid) )
       call check( nf90_put_att(nc_id, level_varid, 'long_name', 'level') )
       call check( nf90_put_att(nc_id, level_varid, 'units', '1') )

       call check( nf90_def_var(nc_id,'uvel',NF90_DOUBLE,(/xstag_dimid,ystag_dimid,level_dimid,time_dimid/),uvel_varid) )
       call check( nf90_put_att(nc_id, uvel_varid, 'long_name', 'ice velocity in x direction') )
       call check( nf90_put_att(nc_id, uvel_varid, 'units', 'meter/year') )

       call check( nf90_def_var(nc_id,'vvel',NF90_DOUBLE,(/xstag_dimid,ystag_dimid,level_dimid,time_dimid/),vvel_varid) )
       call check( nf90_put_att(nc_id, vvel_varid, 'long_name', 'ice velocity in y direction') )
       call check( nf90_put_att(nc_id, vvel_varid, 'units', 'meter/year') )

    end if

    call check( nf90_enddef(nc_id) )

    ! populate variables

    call check( nf90_put_var(nc_id,x_varid,xs) )
    call check( nf90_put_var(nc_id,y_varid,ys) )
    call check( nf90_put_var(nc_id,xstag_varid,xstag) )
    call check( nf90_put_var(nc_id,ystag_varid,ystag) )
    call check( nf90_put_var(nc_id,time_varid, (/ 1 /) ) )


    ! write the arrays out to the netcdf file and close it
    call check( nf90_put_var(nc_id, thck_varid, thck) )
    call check( nf90_put_var(nc_id, topog_varid, topog) )
    call check( nf90_put_var(nc_id, kinbcmask_varid, kinbcmask) )
    
    if (write_vel) then
       call check( nf90_put_var(nc_id,level_varid, level) )
       call check( nf90_put_var(nc_id, uvel_varid, uvelhom) )
       call check( nf90_put_var(nc_id, vvel_varid, vvelhom) )
    end if

    call check( nf90_close(nc_id) )

  end subroutine write_nc_file

  subroutine check(status_code)

    integer,intent(in) :: status_code

    if (status_code /= 0) then
       write(*,*) 'fatal netcdf error:',nf90_strerror(status_code)
       stop 1
    end if

  end subroutine check

end program nc_gen_input
