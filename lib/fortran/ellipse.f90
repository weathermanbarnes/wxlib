! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 	 DynLib -- highly efficient solvers for elliptical PDEs
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!@ Highly efficient solvers for elliptical PDEs, using sparse matrix inversion algorithms
!@
!@ The discretised PDEs are expressed as an operator applying the homogenous
!@ part of the equation to a given input field.
module ellipse
  use kind
  use config
  !
  implicit none
contains
  !
  !@ Generic routine for inverting a 3d operator op
  !@
  !@ The routine solves the equation::
  !@
  !@             op(res) = b 
  !@
  !@ for ``res`` using the best available method. Currently, this is the Full MultiGrid 
  !@ (FMG) method. This default is not likely to change anytime soon.
  !@ 
  !@ Parameters
  !@ ----------
  !@ 
  !@ op : function
  !@      A function or subroutine that implements the operator to be inverted, calculating
  !@      ``op(res)`` from ``res``.
  !@ b : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Inhomogeneous part of the discretised PDE ``op(res) = b``.
  !@ 
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Approximated solution of the PDE.
  !@ float64
  !@     Squared residual of the approximation
  subroutine invert_3d(x, ressq, nx,ny,nz, op, b)
    real(kind=nr), intent(in)  :: b(nz,ny,nx)
    real(kind=nr), intent(out) :: x(nz,ny,nx), ressq
    !
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) x 
    !
    interface
      subroutine op(opa, a)
        use kind
        real(kind=nr), intent(inout) :: a(:,:,:)
        real(kind=nr), intent(out) :: opa(:,:,:)
      end subroutine
    end interface
    ! -----------------------------------------------------------------
    !
    call fmg(x, ressq, nx,ny,nz, op, b, 3_ni, 2_ni, 1_ni)
    !
  end subroutine
  !
  !@ Full multigrid method for inverting a 3d operator op
  !@
  !@ The routine solves the equation::
  !@
  !@             op(res) = b 
  !@
  !@ for ``res``. BiCGstab is used as the smoother function, linear interpolation is 
  !@ used for restriction and prolongation.
  !@ 
  !@ Parameters
  !@ ----------
  !@ 
  !@ op : function
  !@      A function or subroutine that implements the operator to be inverted, calculating
  !@      ``op(res)`` from ``res``.
  !@ b : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Inhomogeneous part of the discretised PDE ``op(res) = b``.
  !@ niter_pre : int
  !@     Number of iterations of the smoother routine before restricting.
  !@ niter_post : int
  !@     Number of iterations of the smoother routine after prolongation.
  !@ nrec : int 
  !@     Number of recursions.
  !@ 
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Approximated solution of the PDE.
  !@ float64
  !@     Squared residual of the approximation
  !@ 
  !@ See Also
  !@ --------
  !@ :meth:`bicgstab`, :meth:`cg`
  subroutine fmg(x, ressq, nx,ny,nz, op, b, niter_pre, niter_post, nrec)
    real(kind=nr), intent(in)  :: b(nz,ny,nx)
    real(kind=nr), intent(out) :: x(nz,ny,nx), ressq
    integer(kind=ni), intent(in) :: niter_pre, niter_post, nrec
    !
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) x 
    !
    interface
      subroutine op(opa, a)
        use kind
        real(kind=nr), intent(inout) :: a(:,:,:)
        real(kind=nr), intent(out) :: opa(:,:,:)
      end subroutine
    end interface
    !
    integer(kind=ni), parameter :: nstep_max = 20_ni
    !
    real(kind=nr), allocatable :: xlow(:,:,:), xhigh(:,:,:), bi(:,:,:)
    integer(kind=ni) :: nxs(0_ni:nstep_max), nys(0_ni:nstep_max), nzs(0_ni:nstep_max), &
            &           n, n1, m, nstep, nxlow,nylow,nzlow
    ! -----------------------------------------------------------------
    !
    ! init
    !   [make sure this definition of the coarser grid sizes is consistent to the
    !    one in mg()!]
    nstep = 0_ni
    nxs(0_ni) = nx
    nys(0_ni) = ny
    nzs(0_ni) = nz
    nxlow = nx
    nylow = ny
    nzlow = nz
    do while ((nxlow == 2_ni .or. nylow == 2_ni) .and. nzlow == 2_ni )
       nstep = nstep + 1_ni
       if ( nstep > nstep_max ) then
          write(*,*) 'too many grids (nstep > ', nstep_max, ')'
          stop
       end if
       nxlow = min(nint(nxlow/2.0_nr), 2_ni)
       nylow = min(nint(nylow/2.0_nr), 2_ni)
       nzlow = min(nint(nzlow/2.0_nr), 2_ni)
       nxs(nstep) = nxlow
       nys(nstep) = nylow
       nzs(nstep) = nzlow
    end do
    !
    ! solution at coarsest grid
    allocate(xlow(nzs(nstep),nys(nstep),nxs(nstep)))
    xlow = 0.0_nr
    !
    do n = nstep,1_ni,-1_ni
       n1 = n-1_ni
       ! interpolate solution to higher resolution
       allocate(xhigh(nzs(n1),nys(n1),nxs(n1)), bi(nzs(n1),nys(n1),nxs(n1)) )
       call interpolate(xhigh, nxs(n1),nys(n1),nzs(n1), xlow, nxs(n),nys(n),nzs(n))
       ! interpolate forcing b to current resolution
       call interpolate(bi, nxs(n1),nys(n1),nzs(n1), b, nx,ny,nz)
       !
       ! refine resolution at higher resolution and save back to xlow
       deallocate(xlow)
       allocate(xlow(nzs(n1),nys(n1),nxs(n1)))
       do m = 1_ni,nrec
          call mg(xlow, ressq, nxs(n1),nys(n1),nzs(n1), & 
                  & op, bi, xhigh, niter_pre, niter_post, nrec)
          if ( m < nrec ) xhigh(:,:,:) = xlow(:,:,:)
       end do
       ! 
       ! prepare next iteration
       deallocate(xhigh, bi)
    end do
    !
    x(:,:,:) = xlow(:,:,:)
    deallocate(xlow)
    !
    return
  end subroutine
  !
  !@ One multigrid step
  !@
  !@ The routine is used by FMG, and is not intended to be called directly.
  recursive subroutine mg(x, ressq, nx,ny,nz, op, b, x0, niter_pre, niter_post, nrec)
    real(kind=nr), intent(in)  :: b(nz,ny,nx), x0(nz,ny,nx)
    real(kind=nr), intent(out) :: x(nz,ny,nx), ressq
    integer(kind=ni), intent(in) :: niter_pre, niter_post, nrec
    !
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) x, x0
    !f2py depend(nn) ny, nz
    !
    interface
      subroutine op(opa, a)
        use kind
        real(kind=nr), intent(inout) :: a(:,:,:)
        real(kind=nr), intent(out) :: opa(:,:,:)
      end subroutine
    end interface
    !
    real(kind=nr) :: res(nz,ny,nx), xadd(nz,ny,nx)
    real(kind=nr), allocatable :: reslow(:,:,:), zero(:,:,:), xlow(:,:,:)
    integer(kind=ni) :: nxlow, nylow, nzlow, n
    ! -----------------------------------------------------------------
    !
    ! init
    !   [make sure this definition of the coarser grid sizes is consistent to the
    !    one in fmg()!]
    nxlow = min(nint(nx/2.0_nr), 2_ni)
    nylow = min(nint(ny/2.0_nr), 2_ni)
    nzlow = min(nint(nz/2.0_nr), 2_ni)
    allocate(reslow(nzlow,nylow,nxlow), zero(nzlow,nylow,nxlow), xlow(nzlow,nylow,nxlow))
    !
    ! pre-smoothing
    call bicgstab(x, ressq, nx,ny,nz, op, b, x0, niter_pre)
    call op(res, x)
    res(:,:,:) = b - res(:,:,:)
    !
    ! restict
    call interpolate(reslow, nxlow,nylow,nzlow, res, nx,ny,nz)
    zero = 0.0
    !
    ! recursively solve on lower resolutions
    if ( (nxlow == 2_ni .or. nylow == 2_ni) .and. nzlow == 2_ni ) then
       do n = 1_ni,nrec
          call mg(xlow, ressq, nxlow,nylow,nzlow, &
                  &       op, reslow, zero, niter_pre, niter_post, nrec)
       end do
    else
       xlow = zero
    endif
    !
    ! prolongate
    call interpolate(xadd, nx,ny,nz, xlow, nxlow,nylow,nzlow)
    x(:,:,:) = x(:,:,:) + xadd(:,:,:)
    !
    ! post-smoothing
    call bicgstab(x, ressq, nx,ny,nz, op, b, x0, niter_post)
    !
    ! clean-up
    deallocate(reslow, zero, xlow)
    !
  end subroutine
  !
  !@ Interpolation helper function for interpolating between the different grids
  !@ in the fmg / mg methods.
  !@
  !@ The routine is used by FMG/MG, and is not intended to be called directly.
  subroutine interpolate(xo, nxo,nyo,nzo, xi, nxi,nyi,nzi)
    real(kind=nr), intent(in)  :: xi(nzi,nyi,nxi)
    real(kind=nr), intent(out) :: xo(nzo,nyo,nxo)
    integer(kind=ni) :: nxo,nyo,nzo, nxi,nyi,nzi
    !
    real(kind=nr) :: plane(nyi,nxi), row(nxi), x,y,z, xscale,yscale,zscale
    integer(kind=ni) :: i,j,k
    ! -----------------------------------------------------------------
    !
    xscale = real(nxi, nr)/nxo
    yscale = real(nyi, nr)/nyo
    zscale = real(nzi, nr)/nzo
    !
    z = 1.0_nr
    do k = 1_ni,nzo
       call interpolate_3d_2d(plane, nxi,nyi,nzi, xi, z, z+zscale)
       !
       y = 1.0_nr
       do j = 1_ni,nyo
          call interpolate_2d_1d(row, nxi,nyi, plane, y, y+yscale)
          !
          x = 1.0_nr
          do i = 1_ni,nxo
             call interpolate_1d_0d(xo(k,j,i), nxi, row, x, x+xscale)
             !
             x = x + xscale
          end do
          y = y + yscale
       end do
       z = z + zscale
    end do
    !
    return
  end subroutine
  !
  ! The three 1d interpolation functions below are actually one and the same
  ! function, but Fortran makes it impossible to avoid the duplication
  !
  !@ Interpolation of one plane
  !@
  !@ The routine is used by FMG/MG, and is not intended to be called directly.
  subroutine interpolate_3d_2d(res, nx,ny,nz, dat, z0, z1)
    real(kind=nr), intent(in)  :: dat(nz,ny,nx), z0, z1
    real(kind=nr), intent(out) :: res(ny,nx)
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny) res
    !
    real(kind=nr) :: overlap
    integer(kind=ni) :: k, fst, lst
    ! -----------------------------------------------------------------
    !
    res(:,:) = 0.0_nr
    !
    fst = int(z0)
    lst = min(nz, ceiling(z1))
    do k = fst,lst
       overlap = min(z1, real(k+1_ni,nr)) - max(z0, real(k,nr))
       res(:,:) = res(:,:) + overlap*dat(k,:,:)
    end do
    !
    res(:,:) = res(:,:)/(z1-z0)
    !
    return
  end subroutine
  !
  !@ Interpolation of one row
  !@
  !@ The routine is used by FMG/MG, and is not intended to be called directly.
  subroutine interpolate_2d_1d(res, nx,ny, dat, y0, y1)
    real(kind=nr), intent(in)  :: dat(ny,nx), y0, y1
    real(kind=nr), intent(out) :: res(nx)
    integer(kind=ni) :: nx,ny
    !f2py depend(nx) res
    !
    real(kind=nr) :: overlap
    integer(kind=ni) :: j, fst, lst
    ! -----------------------------------------------------------------
    !
    res(:) = 0.0_nr
    !
    fst = int(y0)
    lst = min(ny, ceiling(y1))
    do j = fst,lst
       overlap = min(y1, real(j+1_ni,nr)) - max(y0, real(j,nr))
       res(:) = res(:) + overlap*dat(j,:)
    end do
    !
    res(:) = res(:)/(y1-y0)
    !
    return
  end subroutine
  !
  !@ Interpolation of one value
  !@
  !@ The routine is used by FMG/MG, and is not intended to be called directly.
  subroutine interpolate_1d_0d(res, nx, dat, x0, x1)
    real(kind=nr), intent(in)  :: dat(nx), x0, x1
    real(kind=nr), intent(out) :: res
    integer(kind=ni) :: nx
    !
    real(kind=nr) :: overlap
    integer(kind=ni) :: i, fst, lst
    ! -----------------------------------------------------------------
    !
    res = 0.0_nr
    !
    fst = int(x0)
    lst = min(nx, ceiling(x1))
    do i = fst,lst
       overlap = min(x1, real(i+1_ni,nr)) - max(x0, real(i,nr))
       res = res + overlap*dat(i)
    end do
    !
    res = res/(x1-x0)
    !
    return
  end subroutine
  !
  !@ The BI-Conqugate Gradient STABilized (BICGSTAB) method for inverting a 3d operator op
  !@
  !@ The routine solves the equation::
  !@
  !@             op(res) = b 
  !@
  !@ for ``res``. 
  !@
  !@ Reference: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging Variant 
  !@ of Bi-CG for the Solution of Nonsymmetric Linear Systems". SIAM Journal on Scientific and 
  !@ Statistical Computing 13: 631–644.
  !@ 
  !@ Parameters
  !@ ----------
  !@ 
  !@ op : function
  !@      A function or subroutine that implements the operator to be inverted, calculating
  !@      ``op(res)`` from ``res``.
  !@ b : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Inhomogeneous part of the discretised PDE ``op(res) = b``.
  !@ x0 : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     First guess for the iteration, might be just zeros.
  !@ niter : int
  !@     Number of iterations.
  !@ 
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Approximated solution of the PDE.
  !@ float64
  !@     Squared residual of the approximation
  !@ 
  !@ See Also
  !@ --------
  !@ :meth:`fmg`, :meth:`cg`
  subroutine bicgstab(x, ressq, nx,ny,nz, op, b, x0, niter)
    real(kind=nr), intent(in)  :: b(nz,ny,nx), x0(nz,ny,nx)
    real(kind=nr), intent(out) :: x(nz,ny,nx), ressq
    integer(kind=ni), intent(in) :: niter
    !
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) x, x0
    !
    interface
      subroutine op(opa, a)
        use kind
        real(kind=nr), intent(inout) :: a(:,:,:)
        real(kind=nr), intent(out) :: opa(:,:,:)
      end subroutine
    end interface
    !
    real(kind=nr) :: res(nz,ny,nx), res0(nz,ny,nx), &
            &        ret(nz,ny,nx), opret(nz,ny,nx), &
            &        pre(nz,ny,nx), oppre(nz,ny,nx), &
            &        alpha, beta, omega, rho, rhoo
    integer(kind=ni) :: n
    ! -----------------------------------------------------------------
    !
    x   (:,:,:) = x0(:,:,:)
    call op(x, res)
    res (:,:,:) = b(:,:,:) - res(:,:,:)
    res0(:,:,:) = res(:,:,:)
    pre (:,:,:) = res(:,:,:)
    !
    rho = sum(res0(:,:,:)*res(:,:,:))
    do n = 1_ni,niter
       ! Find optimal alpha
       call op(pre, oppre)
       alpha = rho/sum(res0(:,:,:)*oppre(:,:,:))
       !
       ! Find optimal omega
       ret(:,:,:) = res(:,:,:) - alpha*oppre(:,:,:)
       call op(ret, opret)
       omega = sum(opret(:,:,:)*ret(:,:,:))/sum(opret(:,:,:)*opret(:,:,:))
       !
       ! Update solution and residual
       x  (:,:,:) = x  (:,:,:) + alpha*pre(:,:,:) + omega*ret(:,:,:)
       res(:,:,:) = ret(:,:,:) - omega*opret(:,:,:)
       !
       ! Find new search direction pre and prepare next iteration
       rhoo = rho
       rho  = sum(res(:,:,:)*res0(:,:,:))
       beta = alpha/omega * rho/rhoo
       pre(:,:,:) = res(:,:,:) + beta*(pre(:,:,:) - omega*oppre(:,:,:))
    end do
    !
    ressq = sum(res(:,:,:)*res(:,:,:))
    !
  end subroutine
  !
  !@ The Conqugate Gradient (CG) method for inverting a 3d operator op
  !@
  !@ The routine solves the equation::
  !@
  !@             op(res) = b 
  !@
  !@ for ``res``. 
  !@ 
  !@ Parameters
  !@ ----------
  !@ 
  !@ op : function
  !@      A function or subroutine that implements the operator to be inverted, calculating
  !@      ``op(res)`` from ``res``.
  !@ b : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Inhomogeneous part of the discretised PDE ``op(res) = b``.
  !@ x0 : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     First guess for the iteration, might be just zeros.
  !@ niter : int
  !@     Number of iterations.
  !@ 
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Approximated solution of the PDE.
  !@ float64
  !@     Squared residual of the approximation
  !@ 
  !@ See Also
  !@ --------
  !@ :meth:`fmg`, :meth:`cg`
  subroutine cg(x, ressq, nx,ny,nz, op, b, x0, niter)
    real(kind=nr), intent(in)  :: b(nz,ny,nx), x0(nz,ny,nx)
    real(kind=nr), intent(out) :: x(nz,ny,nx), ressq
    integer(kind=ni), intent(in) :: niter
    !
    integer(kind=ni) :: nx,ny,nz
    !f2py depend(nx,ny,nz) x, x0
    !
    interface
      subroutine op(opa, a)
        use kind
        real(kind=nr), intent(inout) :: a(:,:,:)
        real(kind=nr), intent(out) :: opa(:,:,:)
      end subroutine
    end interface
    !
    real(kind=nr) :: res(nz,ny,nx), pre(nz,ny,nx), oppre(nz,ny,nx), &
            &        alpha, beta, ressqo
    integer(kind=ni) :: n
    ! -----------------------------------------------------------------
    !
    x(:,:,:) = x0(:,:,:)
    call op(x, res)
    res(:,:,:) = b(:,:,:) - res(:,:,:)
    pre(:,:,:) = res(:,:,:)
    !
    ressqo = sum(res(:,:,:)*res(:,:,:))
    do n = 1_ni,niter
       ! Find optimal alpha
       call op(pre, oppre)
       alpha = ressqo/sum(pre(:,:,:)*oppre(:,:,:))
       !
       ! Update solution and residual
       x  (:,:,:) = x  (:,:,:) + alpha*pre(:,:,:)
       res(:,:,:) = res(:,:,:) - alpha*pre(:,:,:)
       ressq      = sum(res(:,:,:)*res(:,:,:))
       !
       ! Find new op-conjugent search direction pre
       beta = ressq/ressqo
       pre(:,:,:) = res(:,:,:) + beta*pre(:,:,:)
       !
       ! Prepare next iteration
       ressqo = ressq
    end do
    !
  end subroutine
  !
end module
