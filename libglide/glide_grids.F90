!Helper module containing routines to move between staggered and
!unstaggered grids
#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide_grids
    use glimmer_global, only : dp
    implicit none

contains
 
  !> A special staggering algorithm that is meant to conserve mass when operating on thickness fields
  !! This incorporates Anne Le Brocq's nunatak fix and the calving front fix.
  !!
  !! \param ipvr The input thickness field (ewn x nsn)
  !! \param opvr The output (staggered) thickness field (ewn - 1 x nsn - 1)
  !! \param ewn
  !! \param nsn
  !! \param usrf   Surface elevation field, non-staggered.
  !! \param thklim Minimum thickness to enable ice dynamics
  !! \param mask   Geometry mask field (used for determining the location of shelf fronts)
  !<
  subroutine stagthickness(ipvr,opvr,ewn,nsn,usrf,thklim,mask)
    use glimmer_paramets, only: thk0
    use glide_mask
    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:) :: ipvr
    
    real(dp), intent(in), dimension(:,:) :: usrf
    real(dp) :: thklim
    integer, intent(in), dimension(:,:) :: mask
    
    integer :: ewn,nsn,ew,ns,n
    real(dp) :: tot

        do ns = 1,nsn-1
            do ew = 1,ewn-1

                !If any of our staggering points are shelf front, ignore zeros when staggering
                if (any(is_calving(mask(ew:ew+1, ns:ns+1)))) then

                !Use the "only nonzero thickness" staggering criterion for ALL marginal ice. For
                ! reasons that are not entirely clear, this corrects an error whereby the land ice 
                ! margin is defined incorrectly as existing one grid cell too far inland from where 
                ! it should be.  
                !if (any(has_ice(mask(ew:ew+1,ns:ns+1)))) then
                    n = 0
                    tot = 0
    
                    if (abs(ipvr(ew,ns)) > 1e-10) then
                        tot = tot + ipvr(ew,ns)
                        n   = n   + 1
                    end if
                    if (abs(ipvr(ew+1,ns)) > 1e-10) then
                        tot = tot + ipvr(ew+1,ns)
                        n   = n   + 1
                    end if
                    if (abs(ipvr(ew,ns+1)) > 1e-10) then
                        tot = tot + ipvr(ew,ns+1)
                        n   = n   + 1
                    end if
                    if (abs(ipvr(ew+1,ns+1)) > 1e-10) then
                        tot = tot + ipvr(ew+1,ns+1)
                        n   = n   + 1
                    end if
                    if (n > 0) then
                        opvr(ew,ns) = tot/n
                    else
                        opvr(ew,ns) = 0
                    end if
                !The following cases relate to Anne LeBroque's fix for nunataks
                !ew,ns cell is ice free:
                else if (ipvr(ew,ns) <= thklim/thk0 .and. &
                   ((usrf(ew,ns) >= usrf(ew+1,ns) .and. ipvr(ew+1,ns) >= thklim/thk0) &
                    .or. (usrf(ew,ns) >= usrf(ew,ns+1) .and. ipvr(ew,ns+1) >= thklim/thk0))) then
                        opvr(ew,ns) = 0.0

                !ew+1,ns cell is ice free:
                else if (ipvr(ew+1,ns) <= thklim/thk0 .and. &
                    ((usrf(ew+1,ns) >= usrf(ew,ns) .and. ipvr(ew,ns) >= thklim/thk0) &
                    .or. (usrf(ew+1,ns) >= usrf(ew+1,ns+1) .and. ipvr(ew+1,ns+1) >= thklim/thk0))) then
                        opvr(ew,ns) = 0.0
    
                !ew,ns+1 cell is ice free:
                else if (ipvr(ew,ns+1) <= thklim/thk0 .and. &
                    ((usrf(ew,ns+1) >= usrf(ew,ns) .and. ipvr(ew,ns) >= thklim/thk0) &
                    .or. (usrf(ew,ns+1) >= usrf(ew+1,ns+1) .and. ipvr(ew+1,ns+1) >= thklim/thk0))) then
                        opvr(ew,ns) = 0.0
    
                !ew+1,ns+1 cell is ice free:
                else if (ipvr(ew+1,ns+1) <= thklim/thk0 .and. &
                    ((usrf(ew+1,ns+1) >= usrf(ew+1,ns) .and. ipvr(ew+1,ns) >=thklim/thk0) &
                    .or. (usrf(ew+1,ns+1) >= usrf(ew,ns+1) .and. ipvr(ew,ns+1) >=thklim/thk0))) then
                        opvr(ew,ns) = 0.0
                
                !Standard Staggering   !! Not needed if only-nonzero-thickness staggering scheme is used
                else
                        opvr(ew,ns) = (ipvr(ew+1,ns) + ipvr(ew,ns+1) + &
                               ipvr(ew+1,ns+1) + ipvr(ew,ns)) / 4.0d0
                end if
  
        end do
    end do

  end subroutine stagthickness

 
  !> Moves a variable from the ice grid to the velocity grid by averaging onto the centroids
  !! \param ipvr Input variable (on the ice grid)
  !! \param opvr Output variable (on the velocity grid)
  !! \param ewn
  !! \param nsn
  !<
  subroutine stagvarb(ipvr,opvr,ewn,nsn)
    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:)  :: ipvr
    
    integer, intent(in) :: ewn,nsn

        opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                                 ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0
  end subroutine stagvarb
end module glide_grids
