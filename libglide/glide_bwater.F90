module glide_bwater
   use glimmer_global, only: rk, sp
   use glide_types

contains
  subroutine calcbwat(model,which,bmlt,bwat,thck,topg,btem,floater)

    use glimmer_global, only : dp 
    use glimmer_paramets, only : thk0
    use glide_thck
    use glide_temp_utils, only: calcpmptb
    use glide_grids, only: stagvarb

    implicit none

    type(glide_global_type) :: model
    integer, intent(in) :: which

    real(dp), dimension(:,:), intent(inout) :: bwat
    real(dp), dimension(:,:), intent(in) :: bmlt, thck, topg, btem
    logical, dimension(:,:), intent(in) :: floater

    real(dp), dimension(2), parameter :: &
         blim = (/ 0.00001 / thk0, 0.001 / thk0 /)

    real(dp) :: dwphidew, dwphidns, dwphi, pmpt, bave

    integer :: t_wat,ns,ew

    select case (which)
    case(0)

       do t_wat = 1, model%tempwk%nwat
          do ns = 1,model%general%nsn
             do ew = 1,model%general%ewn

                if (model%numerics%thklim < thck(ew,ns) .and. .not. floater(ew,ns)) then
                   bwat(ew,ns) = (model%tempwk%c(1) * bmlt(ew,ns) + model%tempwk%c(2) * bwat(ew,ns)) / &
                        model%tempwk%c(3)
                   if (blim(1) > bwat(ew,ns)) then
                      bwat(ew,ns) = 0.0d0
                   end if
                else
                   bwat(ew,ns) = 0.0d0
                end if

             end do
          end do
       end do

       model%tempwk%smth = 0.
       do ns = 2,model%general%nsn-1
          do ew = 2,model%general%ewn-1
             call smooth_bwat(ew-1,ew,ew+1,ns-1,ns,ns+1)
          end do
       end do
       ! apply periodic BC
       if (model%options%periodic_ew.eq.1) then
          do ns = 2,model%general%nsn-1
             call smooth_bwat(model%general%ewn-1,1,2,ns-1,ns,ns+1)
             call smooth_bwat(model%general%ewn-1,model%general%ewn,2,ns-1,ns,ns+1)
          end do
       end if

       bwat(1:model%general%ewn,1:model%general%nsn) = model%tempwk%smth(1:model%general%ewn,1:model%general%nsn)

    case(1)
       ! apply periodic BC
       if (model%options%periodic_ew.eq.1) then
          write(*,*) 'Warning, periodic BC are not implement for this case yet'
       end if
       ! ** add any melt_water

       bwat = max(0.0d0,bwat + model%numerics%dttem * bmlt)

       model%tempwk%wphi = 0.
       model%tempwk%bwatu = 0.
       model%tempwk%bwatv = 0.
       model%tempwk%fluxew = 0.
       model%tempwk%fluxns = 0.
       model%tempwk%bint = 0.


       ! ** split time evolution into steps to avoid CFL problems

       do t_wat = 1,model%tempwk%nwat

          ! ** find potential surface using paterson p112, eq 4
          ! ** if no ice then set to sea level or land surface potential
          ! ** if frozen then set high 

          do ns = 1,model%general%nsn
             do ew = 1,model%general%ewn
                if (model%numerics%thklim < thck(ew,ns) .and. .not. floater(ew,ns)) then
                   call calcpmptb(pmpt,thck(ew,ns))
                   if (btem(ew,ns) == pmpt) then
                      model%tempwk%wphi(ew,ns) = model%tempwk%c(1) * (topg(ew,ns) + bwat(ew,ns)) + model%tempwk%c(2) * thck(ew,ns)
                   else
                      model%tempwk%wphi(ew,ns) = model%tempwk%c(1) * (topg(ew,ns) + thck(ew,ns))
                   end if
                else 
                   model%tempwk%wphi(ew,ns) = max(model%tempwk%c(1) * topg(ew,ns),0.0d0)
                end if
             end do
          end do

          ! ** determine x,y components of water velocity assuming
          ! ** contstant velocity magnitude and using potential
          ! ** to determine direction

          do ns = 2,model%general%nsn-1
             do ew = 2,model%general%ewn-1
                if (thck(ew,ns) > model%numerics%thklim) then

                   dwphidew = (model%tempwk%wphi(ew+1,ns) - model%tempwk%wphi(ew-1,ns)) / model%tempwk%c(3)       
                   dwphidns = (model%tempwk%wphi(ew,ns+1) - model%tempwk%wphi(ew,ns-1)) / model%tempwk%c(4)  

                   dwphi = - model%tempwk%watvel / sqrt(dwphidew**2 + dwphidns**2)

                   model%tempwk%bwatu(ew,ns) = dwphi * dwphidew  
                   model%tempwk%bwatv(ew,ns) = dwphi * dwphidns  

                else
                   model%tempwk%bwatu(ew,ns) = 0.0d0
                   model%tempwk%bwatv(ew,ns) = 0.0d0
                end if
             end do
          end do

          ! ** use two-step law wendroff to solve dW/dt = -dF/dx - dF/dy

          ! ** 1. find fluxes F=uW

          model%tempwk%fluxew = bwat * model%tempwk%bwatu
          model%tempwk%fluxns = bwat * model%tempwk%bwatv

          ! ** 2. do 1st LW step on staggered grid for dt/2

          do ns = 1,model%general%nsn-1
             do ew = 1,model%general%ewn-1

                bave = 0.25 * sum(bwat(ew:ew+1,ns:ns+1))

                if (bave > 0.0d0) then

                   model%tempwk%bint(ew,ns) = bave - &
                        model%tempwk%c(5) * (sum(model%tempwk%fluxew(ew+1,ns:ns+1)) - sum(model%tempwk%fluxew(ew,ns:ns+1))) - &
                        model%tempwk%c(6) * (sum(model%tempwk%fluxns(ew:ew+1,ns+1)) - sum(model%tempwk%fluxns(ew:ew+1,ns)))

                else
                   model%tempwk%bint(ew,ns) = 0.0d0
                end if
             end do
          end do

          ! ** 3. find fluxes F=uW on staggered grid griven new Ws

          model%tempwk%fluxew(1:model%general%ewn-1,1:model%general%nsn-1) = model%tempwk%bint * 0.25 * &
               (model%tempwk%bwatu(1:model%general%ewn-1,1:model%general%nsn-1) + &
               model%tempwk%bwatu(2:model%general%ewn,1:model%general%nsn-1) + &
               model%tempwk%bwatu(1:model%general%ewn-1,2:model%general%nsn) + &
               model%tempwk%bwatu(2:model%general%ewn,2:model%general%nsn))
          model%tempwk%fluxns(1:model%general%ewn-1,1:model%general%nsn-1) = model%tempwk%bint * 0.25 * &
               (model%tempwk%bwatv(1:model%general%ewn-1,1:model%general%nsn-1) + &
               model%tempwk%bwatv(2:model%general%ewn,1:model%general%nsn-1) + &
               model%tempwk%bwatv(1:model%general%ewn-1,2:model%general%nsn) + &
               model%tempwk%bwatv(2:model%general%ewn,2:model%general%nsn))

          ! ** 4. finally do 2nd LW step to get back on to main grid

          do ns = 2,model%general%nsn-1
             do ew = 2,model%general%ewn-1
                if (bwat(ew,ns) > 0.0d0) then

                   bwat(ew,ns) = bwat(ew,ns) - &
                        model%tempwk%c(7) * (sum(model%tempwk%fluxew(ew,ns-1:ns)) - sum(model%tempwk%fluxew(ew-1,ns-1:ns))) - &
                        model%tempwk%c(8) * (sum(model%tempwk%fluxns(ew-1:ew,ns)) - sum(model%tempwk%fluxns(ew-1:ew,ns-1)))

                else
                   bwat(ew,ns) = 0.0d0
                end if
             end do
          end do
       end do

       where (blim(1) > bwat) 
          bwat = 0.0d0
       end where

    case default
      bwat = 0.0d0
    end select

    ! How to call the flow router.
    ! call advectflow(bwat,phi,bmlt,model%geometry%mask)

    ! now also calculate basal water in velocity coord system
    call stagvarb(model%temper%bwat, &
         model%temper%stagbwat ,&
         model%general%  ewn, &
         model%general%  nsn)

  contains
    subroutine smooth_bwat(ewm,ew,ewp,nsm,ns,nsp)
      ! smoothing basal water distrib
      implicit none
      integer, intent(in) :: ewm,ew,ewp,nsm,ns,nsp
      if (blim(2) < bwat(ew,ns)) then
         model%tempwk%smth(ew,ns) = bwat(ew,ns) + model%paramets%bwat_smooth * &
              (bwat(ewm,ns) + bwat(ewp,ns) + bwat(ew,nsm) + bwat(ew,nsp) - 4.0d0 * bwat(ew,ns))
      else 
         model%tempwk%smth(ew,ns) = bwat(ew,ns)
      end if   
    end subroutine smooth_bwat

  end subroutine calcbwat
  
  subroutine find_dt_wat(dttem,estimate,dt_wat,nwat)
    
    implicit none
    
    real(dp), intent(out) :: dt_wat
    integer, intent(out) :: nwat
    real(dp), intent(in) :: dttem, estimate
    
    nwat = int(dttem/estimate) + 1
    dt_wat = dttem / nwat

  end subroutine find_dt_wat

  subroutine flow_router(surface,input,output,mask,dx,dy)

    !*FD Routes water from input field to its destination, 
    !*FD according to a surface elevation field. The method used 
    !*FD is by Quinn et. al. (1991)

    real(sp),dimension(:,:),intent(in)  :: surface !*FD Surface elevation
    real(rk),dimension(:,:),intent(in)  :: input   !*FD Input water field
    real(rk),dimension(:,:),intent(out) :: output  !*FD Output water field
    integer, dimension(:,:),intent(in)  :: mask    !*FD Masked points
    real(rk),               intent(in)  :: dx      !*FD $x$ grid-length
    real(rk),               intent(in)  :: dy      !*FD $y$ grid-length

    ! Internal variables --------------------------------------

    integer :: nx,ny,k,nn,cx,cy,px,py,x,y
    integer, dimension(:,:),allocatable :: sorted
    real(rk),dimension(:,:),allocatable :: flats,surfcopy
    real(rk),dimension(-1:1,-1:1) :: slopes
    real(rk),dimension(-1:1,-1:1) :: dists
    logical :: flag

    ! Set up grid dimensions ----------------------------------

    nx=size(surface,1) ; ny=size(surface,2)
    nn=nx*ny

    dists(-1,:)=(/4d0,2d0*dx/dy,4d0/)
    dists(0,:)=(/2d0*dy/dx,0d0,2d0*dy/dx/)
    dists(1,:)=dists(-1,:)

    ! Allocate internal arrays and copy data ------------------

    allocate(sorted(nn,2),flats(nx,ny),surfcopy(nx,ny))
    surfcopy=surface

    ! Fill holes in data, and sort heights --------------------

    call fillholes(surfcopy,flats,mask)
    call heights_sort(surfcopy,sorted)

    ! Initialise output with input, which will then be --------
    ! redistributed -------------------------------------------

    output=input

    ! Begin loop over points, highest first -------------------

    do k=nn,1,-1
    
      ! Get location of current point -------------------------

      x=sorted(k,1)
      y=sorted(k,2)

      ! Reset flags and slope arrays --------------------------

      flag=.true.
      slopes=0.0

      ! Loop over adjacent points, and calculate slopes -------

      do cx=-1,1,1
        do cy=-1,1,1
          ! If this is the centre point, ignore
          if (cx==0.and.cy==0) continue
          ! Otherwise do slope calculation 
          px=x+cx ; py=y+cy
          if (px>0.and.px<=nx.and.py>0.and.py<=ny) then
              if (surfcopy(px,py)<surfcopy(x,y)) then
                slopes(cx,cy)=(surfcopy(x,y)-surfcopy(px,py))/dists(cx,cy)
              endif
          endif
        enddo
      enddo

      ! If there are places for the water to drain to, --------
      ! distribute it accordingly -----------------------------

      if (sum(slopes)/=0.0) then

        slopes=slopes/sum(slopes)
        do cx=-1,1
          do cy=-1,1
            px=x+cx ;py=y+cy
            if (slopes(cx,cy)/=0.0) then
              output(px,py)=output(px,py)+output(x,y)*slopes(cx,cy)
            endif
          enddo
        enddo

        ! Having distributed the water, zero the source -------

        output(x,y)=0.0

      endif

      ! End of main loop ----------------------------------------

    enddo

    ! Tidy up -------------------------------------------------

    deallocate(sorted,flats)

  end subroutine flow_router
  
!==============================================================
! Internal subroutines
!==============================================================

  subroutine fillholes(phi,flats,mask)

    implicit none

    real(rk),dimension(:,:),intent(inout) :: phi
    real(rk),dimension(:,:),intent(inout) :: flats
    integer, dimension(:,:),intent(in)    :: mask

    ! Internal variables --------------------------------------

    real(rk),allocatable,dimension(:,:) :: old_phi
    integer, allocatable,dimension(:,:) :: pool

    real(rk) :: pvs(9), max_val
    real(rk), parameter :: null = 1e+20
    integer :: flag,nx,ny,i,j

    ! ---------------------------------------------------------

    nx=size(phi,1) ; ny=size(phi,2)

    allocate(pool(nx,ny),old_phi(nx,ny))

    flag = 1

    ! ---------------------------------------------------------

    do while (flag .eq. 1)

       flag = 0

       old_phi = phi

       do i=2,nx-1
          do j=2,ny-1

             flats(i,j) = 0

             if (mask(i,j) .eq. 1) then

                if (any(old_phi(i-1:i+1,j-1:j+1) < old_phi(i,j))) then
                   pool(i,j) = 0
                else
                   pool(i,j) = 1
                end if

                if (pool(i,j) .eq. 1) then

                   flag = 1

                   pvs = (/ old_phi(i-1:i+1,j-1), old_phi(i-1:i+1,j+1), old_phi(i-1:i+1,j) /)

                   where (pvs == old_phi(i,j))
                      pvs = null
                   end where

                   max_val = minval(pvs)

                   if (max_val .ne. null) then
                      phi(i,j) = max_val
                   else
                      flag = 0
                      flats(i,j) = 1
                   end if

                end if

             end if
          end do
       end do

    end do

    deallocate(pool,old_phi)

  end subroutine fillholes

!==============================================================

  subroutine heights_sort(surface,sorted)

    real(rk),dimension(:,:) :: surface
    integer,dimension(:,:) :: sorted

    integer :: nx,ny,nn,i,j,k
    real(rk),dimension(:),allocatable :: vect
    integer,dimension(:),allocatable :: ind

    nx=size(surface,1) ; ny=size(surface,2)
    nn=size(sorted,1)

    allocate(vect(nn),ind(nn)) 

    if (nn/=nx*ny.or.size(sorted,2).ne.2) then
      print*,'Wrong dimensions'
      stop
    endif

    k=1

    do i=1,nx
      do j=1,ny
        vect(k)=surface(i,j)
        k=k+1
      enddo
    enddo

    call indexx(vect,ind)

    do k=1,nn
      sorted(k,1)=floor(real(ind(k)-1)/real(ny))+1
      sorted(k,2)=mod(ind(k)-1,ny)+1
    enddo

    do k=1,nn
      vect(k)=surface(sorted(k,1),sorted(k,2))
    enddo
    
  end subroutine heights_sort

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  ! The following two subroutines perform an index-sort of an array. 
  ! They are a GPL-licenced replacement for the Numerical Recipes routine indexx. 
  ! They are not derived from any NR code, but are based on a quicksort routine by
  ! Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
  ! in C, and issued under the GNU General Public License. The conversion to 
  ! Fortran 90, and modification to do an index sort was done by Ian Rutt.
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine indexx(array,index)

    use glimmer_log

    !*FD Performs an index sort of \texttt{array} and returns the result in
    !*FD \texttt{index}. The order of elements in \texttt{array} is unchanged.
    !*FD
    !*FD This is a GPL-licenced replacement for the Numerical Recipes routine indexx. 
    !*FD It is not derived from any NR code, but are based on a quicksort routine by
    !*FD Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    !*FD in C, and issued under the GNU General Public License. The conversion to 
    !*FD Fortran 90, and modification to do an index sort was done by Ian Rutt.

    real(rk),dimension(:) :: array !*FD Array to be indexed.
    integer, dimension(:) :: index !*FD Index of elements of \texttt{array}.
    integer :: i

    if (size(array).ne.size(index)) then
      call write_log('ERROR: INDEXX size mismatch.',GM_FATAL,__FILE__,__LINE__)
    endif

    do i=1,size(index)
       index(i)=i
    enddo

    call q_sort_index(array,index,1,size(array))

  end subroutine indexx

!==============================================================

  recursive subroutine q_sort_index(numbers,index,left,right)

    !*FD This is the recursive subroutine actually used by \texttt{indexx}. 
    !*FD
    !*FD This is a GPL-licenced replacement for the Numerical Recipes routine indexx. 
    !*FD It is not derived from any NR code, but are based on a quicksort routine by
    !*FD Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    !*FD in C, and issued under the GNU General Public License. The conversion to 
    !*FD Fortran 90, and modification to do an index sort was done by Ian Rutt.

    implicit none

    real(rk),dimension(:) :: numbers !*FD Numbers being sorted
    integer, dimension(:) :: index   !*FD Returned index
    integer :: left, right           !*FD Limit of sort region

    integer :: ll,rr
    integer :: pv_int,l_hold, r_hold,pivpos
    real(rk) :: pivot

    ll=left
    rr=right

    l_hold = ll
    r_hold = rr
    pivot = numbers(index(ll))
    pivpos=index(ll)

    do
       if (.not.(ll < rr)) exit

       do 
          if  (.not.((numbers(index(rr)) >= pivot) .and. (ll < rr))) exit
          rr=rr-1
       enddo

       if (ll.ne.rr) then
          index(ll) = index(rr)
          ll=ll+1
       endif

       do
          if (.not.((numbers(index(ll)) <= pivot) .and. (ll < rr))) exit
          ll=ll+1
       enddo

       if (ll.ne.rr) then
          index(rr) = index(ll)
          rr=rr-1
       endif
    enddo

    index(ll) = pivpos
    pv_int = ll
    ll = l_hold
    rr = r_hold
    if (ll < pv_int)  call q_sort_index(numbers, index,ll, pv_int-1)
    if (rr > pv_int)  call q_sort_index(numbers, index,pv_int+1, rr)

  end subroutine q_sort_index

end module glide_bwater
