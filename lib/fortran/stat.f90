! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- statistical functions
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module stat
  use kind
  use config
  !
  implicit none
contains
  ! 
  !@ Calculation of min, max, average, standard deviation and dat(tidx)*tidx for linear trend estimation
  !@
  !@ The statistics are calculated for each grid point along the third
  !@ dimension, typically time. All statistics are calculated at once, with only
  !@ one pass through the input data array.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nt,ny,nx) and dtype float64
  !@     Input data for the statistics calculation.
  !@ tidx_start : int
  !@     Absolute index for the first index in the third direction. To be used
  !@     if trends are estimated over several subsequent calls to this function
  !@     covering a longer input data set.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nt : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Minimum value per grid point.
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Maximum value per grid point.
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Mean value per grid point.
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Standard deviation per grid point.
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     tidx*dat(tidx,::), a quantity useful for linear trend estimation
  subroutine basic(minv,maxv,mean,stddev,tprod, nx,ny,nt, dat,tidx_start)
    real(kind=nr), intent(in) :: dat(nt,ny,nx), tidx_start
    real(kind=nr), intent(out) :: mean(ny,nx), stddev(ny,nx), &
                 &                minv(ny,nx), maxv(ny,nx), tprod(ny,nx)
    integer(kind=ni) :: i,j,n, nx,ny,nt
    !f2py depend(nx,ny) mean, stddev, minv, maxv
    ! -----------------------------------------------------------------
    !
    do i=1_ni,nx
       do j=1_ni,ny
          minv  (j,i) = dat(1_ni,j,i)
          maxv  (j,i) = dat(1_ni,j,i)
          mean  (j,i) = 0._nr
          stddev(j,i) = 0._nr
          tprod (j,i) = 0._nr
          do n=1_ni,nt
             minv  (j,i) = min(dat(n,j,i),minv(j,i))
             maxv  (j,i) = max(dat(n,j,i),maxv(j,i))
             mean  (j,i) = mean(j,i)   + dat(n,j,i)
             stddev(j,i) = stddev(j,i) + dat(n,j,i)**2._nr
             tprod (j,i) = tprod(j,i) + dat(n,j,i)*(n+tidx_start)
          end do
          mean  (j,i) = mean(j,i)/nt
          stddev(j,i) = sqrt((stddev(j,i)-(2._nr*nt-1._nr)*mean(j,i)**2._nr)/(nt-1_ni))
       end do
    end do
  end subroutine
  !
  !@ Calculation of confidence intervals for the trends
  !@
  !@ Uses the squared residual as an indicator for the amount of variance
  !@ explained by the given linear trend. As :meth:`basic`, the statistics are 
  !@ calculated for each grid point along the third dimension.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nt,ny,nx) and dtype float64
  !@     Input data for the statistics calculation.
  !@ trend : np.ndarray with shape (ny,nx) and dtype float64
  !@     Estimated linear trend coefficient.
  !@ icept : np.ndarray with shape (ny,nx) and dtype float64
  !@     Offset of the trend line for tidx=0.
  !@ tidx_start : int
  !@     Absolute index for the first index in the third direction. To be used
  !@     if trends are estimated over several subsequent calls to this function
  !@     covering a longer input data set.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nt : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     The squared residual, describing the variance not explained by the linear trend.
  subroutine basic_confidence(ressq, nx,ny,nt, dat,trend,icept,tidx_start)
    real(kind=nr), intent(in) :: dat(nt,ny,nx), trend(ny,nx), icept(ny,nx), tidx_start
    real(kind=nr), intent(out) :: ressq(ny,nx)
    integer(kind=ni) :: i,j,n, nx,ny,nt
    !f2py depend(nx,ny) ressq, trend, icept
    ! -----------------------------------------------------------------
    !
    do i=1_ni,nx
       do j=1_ni,ny
          ressq(j,i) = 0._nr
          do n=1_ni,nt
             ressq(j,i) =  ressq(j,i) & 
                        + (dat(n,j,i) - icept(j,i) - trend(j,i)*(n+tidx_start))**2.0_nr
          end do
       end do
    end do
  end subroutine
  !
  !@ Calculation of min, max, average and standard deviation, respecting a weight function
  !@
  !@ The statistics are calculated for each grid point along the third
  !@ dimension, typically time. All statistics are calculated at once, with only
  !@ one pass through the input data array.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input data for the statistics calculation.
  !@ weight : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Input data weights, do not need to be normalised.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nt : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Minimum value per grid point.
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Maximum value per grid point.
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Mean value per grid point.
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Standard deviation per grid point.
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Sum of the weight function along the third dimension.
  subroutine basic_weighted(minv,maxv,mean,stddev,wsum, nx,ny,nt, dat,weight)
    real(kind=nr), intent(in)  :: dat(nt,ny,nx), weight(nt,ny,nx)
    real(kind=nr), intent(out) :: mean(ny,nx), stddev(ny,nx), wsum(ny,nx), &
                 &                minv(ny,nx), maxv(ny,nx)
    integer(kind=ni) :: i,j,n, nx,ny,nt
    !f2py depend(nx,ny) mean, stddev, wsum, minv, maxv
    !f2py depend(nx,ny,nt) weight
    ! -----------------------------------------------------------------
    !
    do i=1_ni,nx
       do j=1_ni,ny
          minv  (j,i) = dat(1_ni,j,i)
          maxv  (j,i) = dat(1_ni,j,i)
          mean  (j,i) = 0._nr
          stddev(j,i) = 0._nr
          wsum  (j,i) = 0._nr
          do n=1_ni,nt
             minv  (j,i) = min(dat(n,j,i),minv(j,i))
             maxv  (j,i) = max(dat(n,j,i),maxv(j,i))
             mean  (j,i) = mean(j,i)   + weight(n,j,i)*dat(n,j,i)
             stddev(j,i) = stddev(j,i) + weight(n,j,i)*dat(n,j,i)**2._nr
             wsum  (j,i) = wsum(j,i)   + weight(n,j,i)
          end do
          mean  (j,i) = mean(j,i)/wsum(j,i)
          stddev(j,i) = sqrt((stddev(j,i)*nt/wsum(j,i)-(2._nr*nt-1._nr)*mean(j,i)**2._nr)/(nt-1_ni))
       end do
    end do
  end subroutine
  !
  !@ Some basic binned statistics
  !@ 
  !@ Calculates the histogram using predefined bins, the most frequent value
  !@ among the bins and the median. All statistics are calculated at once, with 
  !@ only one pass through the input data array.
  !@
  !@ Parameters
  !@ ----------
  !@
  !@ dat : np.ndarray with shape (nt,ny,nx) and dtype float64
  !@     Input data for the statistics calculation.
  !@ bins : np.ndarray with shape (nbins,)
  !@     A sorted sequence of bin boundaries. For the case that ``bins(n+1) < bins(n)``,
  !@     all values outside the interval ``bins(n+1) ... bins(n)`` are put into
  !@     that bin. In case a value fits into several of the given intervals, values are 
  !@     put only in the first fitting bin.
  !@
  !@ Other parameters
  !@ ----------------
  !@
  !@ nx : int
  !@     Grid size in x-direction.
  !@ ny : int
  !@     Grid size in y-direction.
  !@ nt : int
  !@     Grid size in z- or t-direction.
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Most frequent value in the data set, estimated using the boundaries of
  !@     the bin holding most values.
  !@ np.ndarray with shape (nbins, ny,nx) and dtype float64
  !@     Histogram, following the given bins. The last values hist[nbins,::]
  !@     hold the number of values that did not fit into any of the given bins.
  !@ np.ndarray with shape (ny,nx) and dtype float64
  !@     Median value, approximated by the borders of the bin holding the median value.
  subroutine binned(mfv,hist,med, nx,ny,nt,nbin, dat,bins)
    real(kind=nr),    intent(in)  :: dat(nt,ny,nx), bins(nbin)
    real(kind=nr),    intent(out) :: mfv(ny,nx), med(ny,nx)
    integer(kind=ni), intent(out) :: hist(nbin,ny,nx)
    integer(kind=ni) :: i,j,n,b, nx,ny,nt, nbin, cnt,bmax,histmax
    logical :: lbinned
    !f2py depend(nx,ny) mfv, med
    !f2py depend(nx,ny,nbin) hist
    ! -----------------------------------------------------------------
    !
    do i=1_ni,nx
       do j=1_ni,ny
          ! Init the histogram
          do b=1_ni,nbin
             hist(b,j,i) = 0_ni
          end do
          !
          do n=1_ni,nt
             lbinned = .false.
             ! Try to fit into the defined bin intervals
             do b=1_ni,nbin-1
                ! "lower" boundary < "upper" boundary: normal interval
                if (bins(b+1_ni) >= bins(b)) then
                   if (dat(n,j,i) >= bins(b) .and. dat(n,j,i) < bins(b+1_ni)) then
                      hist(b,j,i) = hist(b,j,i) + 1_ni
                      lbinned = .true.
                      exit
                   end if
                ! "lower" boundary" > "upper" boundary: everything except the interval
                else
                   if (dat(n,j,i) >= bins(b) .or. dat(n,j,i) < bins(b+1_ni)) then
                      hist(b,j,i) = hist(b,j,i) + 1_ni
                      lbinned = .true.
                      exit
                   end if
                end if
             end do
             ! If none of the intervals in bins matched.
             if (.not. lbinned) then
                hist(nbin,j,i) = hist(nbin,j,i) + 1_ni
             end if
          end do
          !
          ! Finding the most frequent value and the median bin
          ! [OBS: only of those elements that fit into the bin intervals!]
          cnt     =  0_ni
          bmax    = -1_ni
          histmax = -1_ni
          lbinned = .false.
          do b=1_ni,nbin-1_ni
             cnt = cnt + hist(b,j,i)
             if (hist(b,j,i) > histmax) then
                histmax = hist(b,j,i)
                bmax    = b
             end if
             if (cnt >= n/2_ni .and. .not. lbinned) then
                med(j,i) = 0.5_nr*(bins(b+1_ni)+bins(b))
                lbinned = .true.
             end if
          end do
          mfv(j,i) = 0.5_nr*(bins(bmax+1_ni)+bins(bmax))
          !
       end do
    end do
  end subroutine
end module
