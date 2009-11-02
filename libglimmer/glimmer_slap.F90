! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_thck.f90 - part of the GLIMMER ice model         + 
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

!> A generic interface to the SLAP library, based on bits of
!! glide_thck
module glimmer_slap

  use glimmer_global, only: dp

  implicit none
  
  !> Sparse matrix type for slap library
  type slapMatrix_type
     private
     integer :: maxSize         !< Largest allowable number of elements
     integer :: numElements = 0 !< Actual number of elements
     integer :: rank            !< Rank of matrix (used for error checking)
     integer, dimension(:),allocatable :: row !< Row index of elements
     integer, dimension(:),allocatable :: col !< Column index of elements
     real(dp),dimension(:),allocatable :: val !< Values of elements
  end type slapMatrix_type

  private
  public :: slapMatrix_type, slapMatrix_init, slapMatrix_insertElement, slapMatrix_resize, slapSolve

contains

  !> Initialise sparse matrix to a given maximum size
  subroutine slapMatrix_init(matrix,rank,maxSize)

    type(slapMatrix_type), intent(inout) :: matrix
    integer,               intent(in)    :: rank
    integer,               intent(in)    :: maxSize

    matrix%maxSize     = maxSize
    matrix%rank        = rank
    matrix%numElements = 0

    ! Reallocate storage

    if (allocated(matrix%row)) deallocate(matrix%row)
    if (allocated(matrix%col)) deallocate(matrix%col)
    if (allocated(matrix%val)) deallocate(matrix%val)

    allocate(matrix%row(matrix%maxSize))
    allocate(matrix%col(matrix%maxSize))
    allocate(matrix%val(matrix%maxSize))

  end subroutine slapMatrix_init

  !-------------------------------------------------------------------
  !> Add element to sparse matrix
  subroutine slapMatrix_insertElement(matrix,value,col,row)

    use glimmer_global, only : dp

    implicit none

    type(slapMatrix_type), intent(inout) :: matrix
    integer,               intent(in)    :: row
    integer,               intent(in)    :: col
    real(dp),              intent(in)    :: value

    if (matrix%numElements==matrix%maxSize) then
       call slapMatrix_resize(matrix,int(matrix%maxSize*1.5))
    end if

    if (value /= 0.0d0) then
      matrix%numElements = matrix%numElements + 1
      matrix%val(matrix%numElements) = value
      matrix%col(matrix%numElements) = col
      matrix%row(matrix%numElements) = row
    end if

  end subroutine slapMatrix_insertElement

  !-------------------------------------------------------------------
  !> Resizes an existing matrix, retaining the already-defined elements
  subroutine slapMatrix_resize(matrix,newMaxSize)

    use glimmer_global, only : dp

    implicit none

    type(slapMatrix_type), intent(inout) :: matrix
    integer,               intent(in)    :: newMaxSize

    integer,dimension(matrix%numElements) :: tmprow
    integer,dimension(matrix%numElements) :: tmpcol
    integer,dimension(matrix%numElements) :: tmpval
    integer :: tmpN

    ! Copy exisiting contents to temporary storage
    if (matrix%numElements/=0) then
       tmprow = matrix%row
       tmpcol = matrix%col
       tmpval = matrix%val
       tmpN   = matrix%numElements
    end if

    ! Reallocate
    call slapMatrix_init(matrix,matrix%rank,newMaxSize)

    ! Copy back
    if (tmpN/=0) then
       matrix%row(1:tmpN) = tmprow
       matrix%col(1:tmpN) = tmpcol
       matrix%val(1:tmpN) = tmpval
       matrix%numElements = tmpN
    end if

  end subroutine slapMatrix_resize

!---------------------------------------------------------------------------------

  subroutine slapSolve(matrix,rhs,answ,iter,err)

    use glimmer_global, only: dp 
    use glimmer_log
    use glimmer_filenames

    implicit none

    type(slapMatrix_type),           intent(in)    :: matrix !< matrix
    real(dp),dimension(matrix%rank), intent(in)    :: rhs    !< right-hand-side of eqn
    real(dp),dimension(matrix%rank), intent(inout) :: answ   !< initial quess/final solution vector
    integer,                         intent(out)   :: iter   !< Number of iterations
    real(dp),                        intent(out)   :: err    !< Error estimate of result

    ! For call to dslucs
    real(dp), parameter :: tol   = 1.0d-12
    integer,  parameter :: isym  = 0
    integer,  parameter :: itol  = 2
    integer,  parameter :: itmax = 101
    real(dp), dimension(20 * matrix%rank) :: rwork ! Real work array
    integer,  dimension(20 * matrix%rank) :: iwork ! Integer work array
    integer :: ierr, mxnelt

    ! Variables for error handling
    character(200) :: message
    character(100) :: errfname
    integer :: lunit

    mxnelt = 20 * matrix%rank

    ! solve the problem using the SLAP package routines     

    call dslucs(matrix%rank,             &  ! n  ... order of matrix a (in)
                rhs,                     &  ! b  ... right hand side vector (in)
                answ,                    &  ! x  ... initial quess/final solution vector (in/out)
                matrix%numElements,      &  ! nelt ... number of non-zeroes in A (in)
                matrix%row,              &  ! ia  ... sparse matrix format of A (in)
                matrix%col,              &  ! ja  ... sparse matrix format of A (in)
                matrix%val,              &  ! a   ... matrix (in)
                isym,                    &  ! isym ... storage method (0 is complete) (in)
                itol,                    &  ! itol ... convergence criteria (2 recommended) (in)
                tol,                     &  ! tol  ... criteria for convergence (in)
                itmax,                   &  ! itmax ... maximum number of iterations (in)
                iter,                    &  ! iter  ... returned number of iterations (out)
                err,                     &  ! err   ... error estimate of solution (out)
                ierr,                    &  ! ierr  ... returned error message (0 is ok) (out)
                0,                       &  ! iunit ... unit for error writes during iteration (0 no write) (in)
                rwork,                   &  ! rwork ... workspace for SLAP routines (in)
                mxnelt,                  &  ! lenw
                iwork,                   &  ! iwork ... workspace for SLAP routines (in)
                mxnelt)                     ! leniw

    ! Handle errors gracefully
    if (ierr /= 0) then

       ! Acquire a file unit, and open the file
       lunit = get_free_unit()
       errfname = trim(process_path('slap_dump.txt'))
       open(lunit,file=errfname,status='unknown')

       ! Output data to file
       call dcpplt(matrix%rank,        &
                   matrix%numElements, &
                   matrix%row,         &
                   matrix%col,         &
                   matrix%val,         &
                   isym,               &
                   lunit)
       write(lunit,*) '***SLAP data ends. Matrix values follows'
       write(lunit,*) matrix%val

       ! Close unit and finish off
       close(lunit)
       write(message,*)'SLAP solution error: Data dumped to ',trim(errfname)
       call write_log(trim(message),GM_FATAL,__FILE__,__LINE__)
    end if

  end subroutine slapSolve


end module glimmer_slap
