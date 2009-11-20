! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_horizcoord.f90 - part of the GLIMMER ice model   + 
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

!> module for handling regular coordinate systems
!!
!! \author Magnus Hagdorn
!! \date June 2006
module glimmer_horizcoord


  use glimmer_global, only: dp, sp

  !> type describing coordinate systems
  type horizCoord_type
     real(kind=dp), dimension(2) :: origin   !< origin of coordinate space
     real(kind=dp), dimension(2) :: delta    !< stepsize in x and y direction
     real(kind=dp), dimension(2) :: delta_r  !< reciprocal stepsize in x and y direction
     integer, dimension(2) :: size    !< extent in x and y direction
  end type horizCoord_type
  
  !> interface of creating new coord system
  interface horizCoord_new
     module procedure horizCoord_new_real, horizCoord_new_pt
  end interface

  !> interface for allocating data for new coord system
  interface horizCoord_allocate
     module procedure horizCoord_allocate_d, horizCoord_allocate_s, horizCoord_allocate_i, horizCoord_allocate_l, &
          horizCoord_allocate_d2, horizCoord_allocate_s2, horizCoord_allocate_i2
  end interface

#ifdef DEBUG_COORDS
  character(len=msg_length), private :: message
#endif
  !NO_RESTART message
  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_COORDINATES
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_GLIMMER_COORDINATES
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_COORDINATES
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_GLIMMER_COORDINATES
!MH!#endif

  !> print horizCoord info to unit
  subroutine horizCoord_print(coord, unit)
    implicit none
    type(horizCoord_type), intent(in) :: coord  !< coordinate system
    integer,intent(in) :: unit                   !< unit to be printed to
    write(unit,*) 'Origin  ',coord%origin
    write(unit,*) 'Delta   ',coord%delta
    write(unit,*) '1/Delta ',coord%delta_r
    write(unit,*) 'Size    ',coord%size
  end subroutine horizCoord_print

  !> create new coordinate system from individual variables
  function horizCoord_new_real(ox, oy, dx, dy, sx, sy)
    implicit none
    real(kind=dp), intent(in) :: ox, oy !< coordinates of origin
    real(kind=dp), intent(in) :: dx, dy !< offsets
    integer, intent(in) :: sx, sy       !< x and y dimension
    type(horizCoord_type) :: horizCoord_new_real
    
    ! origin
    horizCoord_new_real%origin(1) = ox
    horizCoord_new_real%origin(2) = oy
    ! deltas
    horizCoord_new_real%delta(1) = dx
    horizCoord_new_real%delta(2) = dy
    horizCoord_new_real%delta_r(1) = 1.d0/dx
    horizCoord_new_real%delta_r(2) = 1.d0/dy
    ! size
    horizCoord_new_real%size(1) = sx
    horizCoord_new_real%size(2) = sy
  end function horizCoord_new_real

  !> create new coordinate system from points
  function horizCoord_new_pt(o, d, s)
    implicit none
    real(kind=dp), dimension(2), intent(in) :: o  !< coordinates of origin
    real(kind=dp), dimension(2), intent(in) :: d  !< offsets
    integer, dimension(2), intent(in) :: s !< x and y dimension
    type(horizCoord_type) :: horizCoord_new_pt

    ! origin
    horizCoord_new_pt%origin = o
    ! deltas
    horizCoord_new_pt%delta = d
    horizCoord_new_pt%delta_r(:) = 1.d0/d(:)
    ! size
    horizCoord_new_pt%size = s
  end function horizCoord_new_pt

  !> get coordinates of node
  function horizCoord_get_coord(coord,node)
    use glimmer_log
    implicit none
    type(horizCoord_type), intent(in) :: coord  !< coordinate system
    integer, dimension(2), intent(in) :: node       !< node

    real(kind=dp), dimension(2) :: horizCoord_get_coord
  
#ifdef DEBUG_COORDS
    if (.not.horizCoord_node_inside(coord,node)) then
       write(message,*) 'node (',node,') not inside coord system'
       call horizCoord_print(coord,glimmer_get_logunit())
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    end if
#endif

    horizCoord_get_coord(:) = coord%origin(:) + (node(:) - 1)*coord%delta(:)
  end function horizCoord_get_coord

  !> get index of nearest node given coords of a point
  function horizCoord_get_node(coord,point)
    use glimmer_log
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    real(kind=dp), dimension(2), intent(in) :: point      !< point
    
    integer, dimension(2) :: horizCoord_get_node
    
    horizCoord_get_node(:) = 1+floor(0.5+(point(:)-coord%origin(:))*coord%delta_r(:))
    if (horizCoord_get_node(1).eq.coord%size(1)+1) horizCoord_get_node(1) = coord%size(1)
    if (horizCoord_get_node(2).eq.coord%size(2)+1) horizCoord_get_node(2) = coord%size(2)

#ifdef DEBUG_COORDS
    if (.not.horizCoord_node_inside(coord,horizCoord_get_node)) then
       write(message,*) 'point (',point,') not inside coord system'
       call horizCoord_print(coord,glimmer_get_logunit())
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    end if
#endif
  end function horizCoord_get_node

  !> get index of lower-left node of cell into which point falls
  function horizCoord_get_llnode(coord,point)
    use glimmer_log
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    real(kind=dp), dimension(2), intent(in) :: point      !< point
    
    integer, dimension(2) :: horizCoord_get_llnode

    horizCoord_get_llnode(:) = 1+floor((point(:)-coord%origin(:))*coord%delta_r(:))
  end function horizCoord_get_llnode

  !> return true iff node is inside coord system
  function horizCoord_node_inside(coord,node)
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    integer, dimension(2), intent(in) :: node      !< node

    logical horizCoord_node_inside
    
    horizCoord_node_inside = (all(node.ge.1) .and. all(node.le.coord%size))
  end function horizCoord_node_inside

  !> return true iff point is inside coord system
  function horizCoord_point_inside(coord,point)
    use glimmer_log
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    real(kind=dp), dimension(2), intent(in) :: point      !< point
    logical horizCoord_point_inside
    integer i

    horizCoord_point_inside = .true.
    do i=1,2
       horizCoord_point_inside = (point(i).ge.coord%origin(i)) .and. &
            (point(i).le.coord%origin(i)+coord%size(i)*coord%delta(i))
       if (.not.horizCoord_point_inside) then
          exit
       end if
    end do
  end function horizCoord_point_inside
    
  !> linearise node, given coord
  function horizCoord_linearise2d(coord,node)
    use glimmer_log
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    integer, dimension(2), intent(in) :: node      !< node
    integer horizCoord_linearise2d

    horizCoord_linearise2d = -1

#ifdef DEBUG_COORDS
    if (.not.horizCoord_node_inside(coord,node)) then
       write(message,*) 'node (',node,') not inside coord system'
       call write_log(message,GM_ERROR,__FILE__,__LINE__)
       return
    end if
#endif
    
    horizCoord_linearise2d = node(1) + (node(2)-1)*coord%size(1)
  end function horizCoord_linearise2d

  !> expand linearisation
  function horizCoord_delinearise2d(coord, ind)
    use glimmer_log
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    integer, intent(in) :: ind                  !< index
    integer, dimension(2) :: horizCoord_delinearise2d

#ifdef DEBUG_COORDS
    if (ind.lt.1 .or. ind.gt.coord%size(1)*coord%size(2)) then
       write(message,*) 'index ',ind,' outside coord system'
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    end if
#endif

    horizCoord_delinearise2d(1) = mod(ind-1,coord%size(1)) + 1
    horizCoord_delinearise2d(2) = (ind-1)/coord%size(1) + 1
  end function horizCoord_delinearise2d

  !> allocate memory to pointer field
  subroutine horizCoord_allocate_d(coord, field)
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    real(kind=dp), dimension(:,:), pointer :: field !< unallocated field

    allocate(field(coord%size(1),coord%size(2)))
    field = 0.d0
  end subroutine horizCoord_allocate_d
  
  !> allocate memory to pointer field
  subroutine horizCoord_allocate_s(coord, field)
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    real(kind=sp), dimension(:,:), pointer :: field !< unallocated field

    allocate(field(coord%size(1),coord%size(2)))
    field = 0.e0
  end subroutine horizCoord_allocate_s

  !> allocate memory to pointer field
  subroutine horizCoord_allocate_i(coord, field)
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    integer, dimension(:,:), pointer :: field !< unallocated field

    allocate(field(coord%size(1),coord%size(2)))
    field = 0
  end subroutine horizCoord_allocate_i

  !> allocate memory to pointer field
  subroutine horizCoord_allocate_l(coord, field)
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    logical, dimension(:,:), pointer :: field !< unallocated field

    allocate(field(coord%size(1),coord%size(2)))
    field = .FALSE.
  end subroutine horizCoord_allocate_l

  !> allocate memory to pointer field
  subroutine horizCoord_allocate_d2(coord, nup, field)
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    integer, intent(in) :: nup !< the number of vertical points
    real(kind=dp), dimension(:,:,:), pointer :: field !< unallocated field

    allocate(field(nup,coord%size(1),coord%size(2)))
    field = 0.d0
  end subroutine horizCoord_allocate_d2

  !> allocate memory to pointer field
  subroutine horizCoord_allocate_s2(coord, nup, field)
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    integer, intent(in) :: nup !< the number of vertical points
    real(kind=sp), dimension(:,:,:), pointer :: field !< unallocated field

    allocate(field(nup,coord%size(1),coord%size(2)))
    field = 0.d0
  end subroutine horizCoord_allocate_s2

  !> allocate memory to pointer field
  subroutine horizCoord_allocate_i2(coord, nup, field)
    implicit none
    type(horizCoord_type), intent(in) :: coord !< coordinate system
    integer, intent(in) :: nup !< the number of vertical points
    integer, dimension(:,:,:), pointer :: field !< unallocated field

    allocate(field(nup,coord%size(1),coord%size(2)))
    field = 0.d0
  end subroutine horizCoord_allocate_i2

end module glimmer_horizcoord
