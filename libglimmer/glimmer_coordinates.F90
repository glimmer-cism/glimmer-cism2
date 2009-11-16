! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_coordinates.f90 - part of the GLIMMER ice model  + 
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
module glimmer_coordinates


  use glimmer_global, only: dp, sp

  !> type describing coordinate systems
  type coordsystem_type
     real(kind=dp), dimension(2) :: origin   !< origin of coordinate space
     real(kind=dp), dimension(2) :: delta    !< stepsize in x and y direction
     real(kind=dp), dimension(2) :: delta_r  !< reciprocal stepsize in x and y direction
     integer, dimension(2) :: size    !< extent in x and y direction
  end type coordsystem_type
  
  !> interface of creating new coord system
  interface coordsystem_new
     module procedure coordsystem_new_real, coordsystem_new_pt
  end interface

  !> interface for allocating data for new coord system
  interface coordsystem_allocate
     module procedure coordsystem_allocate_d, coordsystem_allocate_s, coordsystem_allocate_i, coordsystem_allocate_l, &
          coordsystem_allocate_d2, coordsystem_allocate_s2, coordsystem_allocate_i2
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

  !> print coordsystem info to unit
  subroutine coordsystem_print(coord, unit)
    implicit none
    type(coordsystem_type), intent(in) :: coord  !< coordinate system
    integer,intent(in) :: unit                   !< unit to be printed to
    write(unit,*) 'Origin  ',coord%origin
    write(unit,*) 'Delta   ',coord%delta
    write(unit,*) '1/Delta ',coord%delta_r
    write(unit,*) 'Size    ',coord%size
  end subroutine coordsystem_print

  !> create new coordinate system from individual variables
  function coordsystem_new_real(ox, oy, dx, dy, sx, sy)
    implicit none
    real(kind=dp), intent(in) :: ox, oy !< coordinates of origin
    real(kind=dp), intent(in) :: dx, dy !< offsets
    integer, intent(in) :: sx, sy       !< x and y dimension
    type(coordsystem_type) :: coordsystem_new_real
    
    ! origin
    coordsystem_new_real%origin(1) = ox
    coordsystem_new_real%origin(2) = oy
    ! deltas
    coordsystem_new_real%delta(1) = dx
    coordsystem_new_real%delta(2) = dy
    coordsystem_new_real%delta_r(1) = 1.d0/dx
    coordsystem_new_real%delta_r(2) = 1.d0/dy
    ! size
    coordsystem_new_real%size(1) = sx
    coordsystem_new_real%size(2) = sy
  end function coordsystem_new_real

  !> create new coordinate system from points
  function coordsystem_new_pt(o, d, s)
    implicit none
    real(kind=dp), dimension(2), intent(in) :: o  !< coordinates of origin
    real(kind=dp), dimension(2), intent(in) :: d  !< offsets
    integer, dimension(2), intent(in) :: s !< x and y dimension
    type(coordsystem_type) :: coordsystem_new_pt

    ! origin
    coordsystem_new_pt%origin = o
    ! deltas
    coordsystem_new_pt%delta = d
    coordsystem_new_pt%delta_r(:) = 1.d0/d(:)
    ! size
    coordsystem_new_pt%size = s
  end function coordsystem_new_pt

  !> get coordinates of node
  function coordsystem_get_coord(coord,node)
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord  !< coordinate system
    integer, dimension(2), intent(in) :: node       !< node

    real(kind=dp), dimension(2) :: coordsystem_get_coord
  
#ifdef DEBUG_COORDS
    if (.not.coordsystem_node_inside(coord,node)) then
       write(message,*) 'node (',node,') not inside coord system'
       call coordsystem_print(coord,glimmer_get_logunit())
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    end if
#endif

    coordsystem_get_coord(:) = coord%origin(:) + (node(:) - 1)*coord%delta(:)
  end function coordsystem_get_coord

  !> get index of nearest node given coords of a point
  function coordsystem_get_node(coord,point)
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    real(kind=dp), dimension(2), intent(in) :: point      !< point
    
    integer, dimension(2) :: coordsystem_get_node
    
    coordsystem_get_node(:) = 1+floor(0.5+(point(:)-coord%origin(:))*coord%delta_r(:))
    if (coordsystem_get_node(1).eq.coord%size(1)+1) coordsystem_get_node(1) = coord%size(1)
    if (coordsystem_get_node(2).eq.coord%size(2)+1) coordsystem_get_node(2) = coord%size(2)

#ifdef DEBUG_COORDS
    if (.not.coordsystem_node_inside(coord,coordsystem_get_node)) then
       write(message,*) 'point (',point,') not inside coord system'
       call coordsystem_print(coord,glimmer_get_logunit())
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    end if
#endif
  end function coordsystem_get_node

  !> get index of lower-left node of cell into which point falls
  function coordsystem_get_llnode(coord,point)
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    real(kind=dp), dimension(2), intent(in) :: point      !< point
    
    integer, dimension(2) :: coordsystem_get_llnode

    coordsystem_get_llnode(:) = 1+floor((point(:)-coord%origin(:))*coord%delta_r(:))
  end function coordsystem_get_llnode

  !> return true iff node is inside coord system
  function coordsystem_node_inside(coord,node)
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    integer, dimension(2), intent(in) :: node      !< node

    logical coordsystem_node_inside
    
    coordsystem_node_inside = (all(node.ge.1) .and. all(node.le.coord%size))
  end function coordsystem_node_inside

  !> return true iff point is inside coord system
  function coordsystem_point_inside(coord,point)
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    real(kind=dp), dimension(2), intent(in) :: point      !< point
    logical coordsystem_point_inside
    integer i

    coordsystem_point_inside = .true.
    do i=1,2
       coordsystem_point_inside = (point(i).ge.coord%origin(i)) .and. &
            (point(i).le.coord%origin(i)+coord%size(i)*coord%delta(i))
       if (.not.coordsystem_point_inside) then
          exit
       end if
    end do
  end function coordsystem_point_inside
    
  !> linearise node, given coord
  function coordsystem_linearise2d(coord,node)
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    integer, dimension(2), intent(in) :: node      !< node
    integer coordsystem_linearise2d

    coordsystem_linearise2d = -1

#ifdef DEBUG_COORDS
    if (.not.coordsystem_node_inside(coord,node)) then
       write(message,*) 'node (',node,') not inside coord system'
       call write_log(message,GM_ERROR,__FILE__,__LINE__)
       return
    end if
#endif
    
    coordsystem_linearise2d = node(1) + (node(2)-1)*coord%size(1)
  end function coordsystem_linearise2d

  !> expand linearisation
  function coordsystem_delinearise2d(coord, ind)
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    integer, intent(in) :: ind                  !< index
    integer, dimension(2) :: coordsystem_delinearise2d

#ifdef DEBUG_COORDS
    if (ind.lt.1 .or. ind.gt.coord%size(1)*coord%size(2)) then
       write(message,*) 'index ',ind,' outside coord system'
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    end if
#endif

    coordsystem_delinearise2d(1) = mod(ind-1,coord%size(1)) + 1
    coordsystem_delinearise2d(2) = (ind-1)/coord%size(1) + 1
  end function coordsystem_delinearise2d

  !> allocate memory to pointer field
  subroutine coordsystem_allocate_d(coord, field)
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    real(kind=dp), dimension(:,:), pointer :: field !< unallocated field

    allocate(field(coord%size(1),coord%size(2)))
    field = 0.d0
  end subroutine coordsystem_allocate_d
  
  !> allocate memory to pointer field
  subroutine coordsystem_allocate_s(coord, field)
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    real(kind=sp), dimension(:,:), pointer :: field !< unallocated field

    allocate(field(coord%size(1),coord%size(2)))
    field = 0.e0
  end subroutine coordsystem_allocate_s

  !> allocate memory to pointer field
  subroutine coordsystem_allocate_i(coord, field)
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    integer, dimension(:,:), pointer :: field !< unallocated field

    allocate(field(coord%size(1),coord%size(2)))
    field = 0
  end subroutine coordsystem_allocate_i

  !> allocate memory to pointer field
  subroutine coordsystem_allocate_l(coord, field)
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    logical, dimension(:,:), pointer :: field !< unallocated field

    allocate(field(coord%size(1),coord%size(2)))
    field = .FALSE.
  end subroutine coordsystem_allocate_l

  !> allocate memory to pointer field
  subroutine coordsystem_allocate_d2(coord, nup, field)
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    integer, intent(in) :: nup !< the number of vertical points
    real(kind=dp), dimension(:,:,:), pointer :: field !< unallocated field

    allocate(field(nup,coord%size(1),coord%size(2)))
    field = 0.d0
  end subroutine coordsystem_allocate_d2

  !> allocate memory to pointer field
  subroutine coordsystem_allocate_s2(coord, nup, field)
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    integer, intent(in) :: nup !< the number of vertical points
    real(kind=sp), dimension(:,:,:), pointer :: field !< unallocated field

    allocate(field(nup,coord%size(1),coord%size(2)))
    field = 0.d0
  end subroutine coordsystem_allocate_s2

  !> allocate memory to pointer field
  subroutine coordsystem_allocate_i2(coord, nup, field)
    implicit none
    type(coordsystem_type), intent(in) :: coord !< coordinate system
    integer, intent(in) :: nup !< the number of vertical points
    integer, dimension(:,:,:), pointer :: field !< unallocated field

    allocate(field(nup,coord%size(1),coord%size(2)))
    field = 0.d0
  end subroutine coordsystem_allocate_i2

end module glimmer_coordinates
