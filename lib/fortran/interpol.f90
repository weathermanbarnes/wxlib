! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- Interpolations
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!@ Interpolation between regular grids
!@ 
!@ Provides for example routines to interpolate data onto a new vertical coordinate
module interpol
  use kind
  use config
  use consts !, only: nan
  !
  implicit none
contains
  !
  !@ Linear interpolation to a given index within the grid.
  !@
  !@ Interpolates all given grid indexes. The given indexes must be in acending order.
  subroutine by_idx_1d(res,nxi,nxo,dat,idxs)
    real(kind=nr), intent(in)  :: dat(nxi), idxs(nxo)
    real(kind=nr), intent(out) :: res(nxo)
    integer(kind=ni), intent(in) :: nxi,nxo
    !
    integer(kind=ni) :: ii,io
    !f2py depend(nxo) res
    ! -----------------------------------------------------------------
    !
    do io = 1_ni,nxo
       if ( isnan(idxs(io)) ) then
          res(io) = nan
       else
          ii = floor(idxs(io))
          if ( ii >= nxi ) then
             res(io) = nan
          else
             res(io) = dat(ii) + (dat(ii+1_ni)-dat(ii))*(idxs(io) - ii)
          end if
       end if
    end do
    !
  end subroutine
  !
  !@ Linear interpolation to given levels of a vertical coordinate
  !@
  !@ Interpolates all given grid indexes. The given levels must be in acending order
  !@ of the original vertical index direction
  subroutine vert_by_coord(res,nx,ny,nz,nn,dat,datv,levels)
    real(kind=nr), intent(out) :: res(nn,ny,nx)
    real(kind=nr), intent(in) :: dat(nz,ny,nx), datv(nz,ny,nx), levels(nn)
    integer(kind=ni), intent(in) :: nx,ny,nz,nn
    !f2py intent(hide,out) res
    !####f2py depend(nx,ny,nz) res
    !####f2py depend(nx,ny) dx, dy
    !
    real(kind=nr) :: idxs(nn)
    integer(kind=ni) :: i,j,ki,ko
    logical :: inverted = .false.
    ! -----------------------------------------------------------------
    !
    ! Check if the new vertical coordinate is inverted with respect to the current one
    if ( datv(2_ni,1_ni,1_ni) < datv(1_ni,1_ni,1_ni) ) then
       inverted = .true.
    end if
    !
    if ( .not. inverted ) then
       do i = 1_ni,nx
          do j = 1_ni,ny
             ki = 1_ni
             do ko = 1_ni,nn
                ! Outside input grid domain
                if ( levels(ko) < datv(1_ni,j,i) .or. levels(ko) >= datv(nz,j,i) ) then
                   idxs(ko) = nan
                ! Loop through the input grid until interpolation point is found
                else
                   do ki = ki,nz-1_ni
                      if ( datv(ki,j,i) <= levels(ko) .and. datv(ki+1_ni,j,i) > levels(ko) ) exit
                   end do
                   idxs(ko) = real(ki, nr) + (levels(ko) - datv(ki,j,i))/(datv(ki+1_ni,j,i) - datv(ki,j,i))
                end if
             end do
             !
             ! Do the actual interpolation and save the results
             call by_idx_1d(res(:,j,i), nz,nn, dat(:,j,i), idxs(:))
          end do
       end do
    ! Structurally identical, but alot of comparison operators are inverted
    else
       do i = 1_ni,nx
          do j = 1_ni,ny
             ki = 1_ni
             do ko = 1_ni,nn
                ! Outside input grid domain
                if ( levels(ko) > datv(1_ni,j,i) .or. levels(ko) <= datv(nz,j,i) ) then
                   idxs(ko) = nan
                ! Loop through the input grid until interpolation point is found
                else
                   do ki = ki,nz-1_ni
                      if ( datv(ki,j,i) >= levels(ko) .and. datv(ki+1_ni,j,i) < levels(ko) ) exit
                   end do
                   idxs(ko) = real(ki, nr) + (levels(ko) - datv(ki,j,i))/(datv(ki+1_ni,j,i) - datv(ki,j,i))
                end if
             end do
             !
             ! Do the actual interpolation and save the results
             call by_idx_1d(res(:,j,i), nz,nn, dat(:,j,i), idxs(:))
          end do
       end do
    end if
    !
  end subroutine
  !
end module interpol
