! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 		DynLib -- F90/F2py Interface to SpherePack
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Module maintained by Clemens Spensberger (csp001@uib.no)
module sphere
  use kind
  use config
  !
  implicit none
  !
  real(kind=nr), allocatable :: wshaes(:), wshses(:), wvhaes(:), wvhses(:)
  integer(kind=ni) :: last_nx=0_ni, last_ny=0_ni, &
          &           lshaes=0_ni, lshses=0_ni, lvhaes=0_ni, lvhses=0_ni
contains
  !
  !@ Spherical harmonics analysis of a scalar field
  !@
  !@ This routine is a wrapper around the SpherePack routines ``shaes`` and ``shaesi``.
  !@ The output of ``shaes`` is passed unchanged as the results of this routine.
  !@
  !@ Parameters
  !@ ----------
  !@ 
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Wind velocity component in the East direction
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Wind velocity component in the North direction
  !@
  !@ Returns
  !@ -------
  !@ 2-tuple of np.ndarray with shape (nz,ny,ny) and dtype float64
  !@     Spherical harmonics coefficients as returned by ``shaes``.
  subroutine sh_analysis(a,b, nx,ny,nz, dat)
    real(kind=nr), intent(in) :: dat(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz
    ! TODO: I tried to integrate the return into two complex arrays, but that
    !       leads to mysterious SEGFAULTS type "double free or corruption"
    real(kind=nr), intent(out) :: a(nz,ny,ny), b(nz,ny,ny)
    !
    !real(kind=nr) :: br(nz,ny,ny), bi(nz,ny,ny), cr(nz,ny,ny), ci(nz,ny,ny)
    integer(kind=ni) :: ierror, k
    !
    integer(kind=8_ni) :: lwork, ldwork    ! lwork might be larger than 32bit!
    real(kind=nr), allocatable :: work(:), dwork(:)
    !
    integer(kind=ni), parameter :: buffer = 10000_ni
    !
    !f2py depend(nz,ny) a, b
    ! -----------------------------------------------------------------
    !
    if ( nx /= last_nx .or. ny /= last_ny ) then
       call reset()
       last_nx = nx
       last_ny = ny
    end if
    !
    lwork = max( 4_ni*(ny+1_ni)**2_ni , (nz+1_ni)*ny*nx )
    ldwork = (ny + 1_ni)*2_ni ! This factor two is not in the SPHEREPACK 
    ! documentation but necessary to avoid shaesi overwriting random parts in memory.
    allocate(work(lwork), dwork(ldwork))
    !
    ! No precalculated wshaes available
    if ( lshaes == 0_ni ) then
       lshaes = ((ny+1_ni)**3_ni)/2_ni + nx + 15_ni
       allocate(wshaes(lshaes))
       call shaesi(ny,nx,wshaes,lshaes,work,lwork,dwork,ldwork,ierror)
       if ( ierror /= 0_ni ) then
          write(*,*) 'Error in shaesi, code', ierror
          stop 1
       end if
    end if
    !
    a(:,:,:) = 0.0_nr
    b(:,:,:) = 0.0_nr
    do k=1_ni,nz
       call shaes(ny,nx,0_ni,1_ni, dat(k,:,:), ny,nx, &
               & a(k,:,:),b(k,:,:), ny,ny, wshaes,lshaes,work,lwork, ierror)
    end do
    if ( ierror /= 0_ni ) then
       write(*,*) 'Error in shaes, code', ierror
       stop 1
    end if
    !
    deallocate(work, dwork)
    !
  end subroutine
  !
  !@ Spherical harmonics synthesis of a scalar field
  !@
  !@ This routine is a wrapper around the SpherePack routines ``shses`` and ``shsesi``.
  !@ The output of ``shses`` is passed unchanged as the results of this routine.
  !@
  !@ Parameters
  !@ ----------
  !@ 
  !@ nx : int
  !@     Number of grid points in zonal direction of the output arrays
  !@ a,b : np.ndarrays with shape (nz,ny,ny) and dtype float64
  !@     Spherical harmonics coefficients as returned by ``sh_analysis``
  !@
  !@ Returns
  !@ -------
  !@ np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Reconstructed field
  subroutine sh_synthesis(res, nx,ny,nz, a,b)
    real(kind=nr), intent(in) :: a(nz,ny,ny), b(nz,ny,ny)
    integer(kind=ni), intent(in) :: nx,ny,nz
    ! TODO: I tried to integrate the return into two complex arrays, but that
    !       leads to mysterious SEGFAULTS type "double free or corruption"
    real(kind=nr), intent(out) :: res(nz,ny,nx)
    !
    !real(kind=nr) :: br(nz,ny,ny), bi(nz,ny,ny), cr(nz,ny,ny), ci(nz,ny,ny)
    integer(kind=ni) :: ierror, k
    !
    integer(kind=8_ni) :: lwork, ldwork    ! lwork might be larger than 32bit!
    real(kind=nr), allocatable :: work(:), dwork(:)
    !
    !f2py depend(nz,ny) bi, cr, ci, u, v
    ! -----------------------------------------------------------------
    !
    if ( nx /= last_nx .or. ny /= last_ny ) then
       call reset()
       last_nx = nx
       last_ny = ny
    end if
    !
    lwork = max( 4_ni*(ny+1_ni)**2_ni , (nz+1_ni)*ny*nx )
    ldwork = ny + 1_ni
    allocate(work(lwork), dwork(ldwork))
    !
    ! No precalculated wshses available
    if ( lshses == 0_ni ) then
       lshses = ((ny+1_ni)**3_ni)/4_ni + nx + 15_ni
       allocate(wshses(lshses))
       call shsesi(ny,nx,wshses,lshses,work,lwork,dwork,ldwork,ierror)
       if ( ierror /= 0_ni ) then
          write(*,*) 'Error in shsesi, code', ierror
          stop 1
       end if
    end if
    !
    res(:,:,:) = 0.0_nr
    do k = 1_ni,nz
       call shses(ny,nx,0_ni,1_ni, res(k,:,:), ny,nx, &
               & a(k,:,:),b(k,:,:), ny,ny, wshses,lshses,work,lwork, ierror)
    end do
    if ( ierror /= 0_ni ) then
       write(*,*) 'Error in shses, code', ierror
       stop 1
    end if
    !
    deallocate(work, dwork)
    !
  end subroutine
  !
  !@ Spherical harmonics analysis of a vector field (u,v)
  !@
  !@ This routine is a wrapper around the SpherePack routines ``vhaes`` and ``vhaesi``.
  !@ The output of ``vhaes`` is passed unchanged as the results of this routine.
  !@
  !@ Parameters
  !@ ----------
  !@ 
  !@ u : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Vector component in the East direction
  !@ v : np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Vector component in the North direction
  !@
  !@ Returns
  !@ -------
  !@ 4-tuple of np.ndarray with shape (nz,ny,ny) and dtype float64
  !@     Spherical harmonics coefficients as returned by ``vhaes``.
  subroutine sh_analysis_vector(br,bi,cr,ci, nx,ny,nz, u,v)
    real(kind=nr), intent(in) :: u(nz,ny,nx), v(nz,ny,nx)
    integer(kind=ni), intent(in) :: nx,ny,nz
    ! TODO: I tried to integrate the return into two complex arrays, but that
    !       leads to mysterious SEGFAULTS type "double free or corruption"
    real(kind=nr), intent(out) :: br(nz,ny,ny), bi(nz,ny,ny), &
            &                     cr(nz,ny,ny), ci(nz,ny,ny)
    !
    !real(kind=nr) :: br(nz,ny,ny), bi(nz,ny,ny), cr(nz,ny,ny), ci(nz,ny,ny)
    integer(kind=ni) :: ierror, k
    !
    integer(kind=8_ni) :: lwork, ldwork    ! lwork might be larger than 32bit!
    real(kind=nr), allocatable :: work(:), dwork(:)
    !
    !f2py depend(nz,ny,nx) v
    !f2py depend(nz,ny) br, bi, cr, ci
    ! -----------------------------------------------------------------
    !
    if ( nx /= last_nx .or. ny /= last_ny ) then
       call reset()
       last_nx = nx
       last_ny = ny
    end if
    !
    lwork = max( 4_ni*(ny+1_ni)**2_ni , (2_ni*nz+1_ni)*ny*nx )
    ldwork = (ny + 1_ni)*2_ni
    allocate(work(lwork), dwork(ldwork))
    !
    ! No precalculated wvhaes available
    if ( lvhaes == 0_ni ) then
       lvhaes = ((ny+1_ni)**3_ni) + nx + 15_ni
       allocate(wvhaes(lvhaes))
       call vhaesi(ny,nx,wvhaes,lvhaes,work,lwork,dwork,ldwork,ierror)
       if ( ierror /= 0_ni ) then
          write(*,*) 'Error in vhaesi, code', ierror
          stop 1
       end if
    end if
    !
    br(:,:,:) = 0.0_nr
    bi(:,:,:) = 0.0_nr
    cr(:,:,:) = 0.0_nr
    ci(:,:,:) = 0.0_nr
    do k = 1_ni,nz
       call vhaes(ny,nx,0_ni,1_ni, -v(k,:,:),u(k,:,:), ny,nx, &
               & br(k,:,:),bi(k,:,:),cr(k,:,:),ci(k,:,:), ny,ny, wvhaes,lvhaes,work,lwork, ierror)
    end do
    if ( ierror /= 0_ni ) then
       write(*,*) 'Error in vhaes, code', ierror
       stop 1
    end if
    !
    deallocate(work, dwork)
    !
  end subroutine
  !
  !@ Spherical harmonics synthesis of a vector field (u,v)
  !@
  !@ This routine is a wrapper around the SpherePack routines ``vhses`` and ``vhsesi``.
  !@ The output of ``vhses`` is passed unchanged as the results of this routine.
  !@
  !@ Parameters
  !@ ----------
  !@ 
  !@ nx : int
  !@     Number of grid points in zonal direction of the output arrays
  !@ br, bi, cr, ci : np.ndarrays with shape (nz,ny,ny) and dtype float64
  !@     Spherical harmonics coefficients as returned by ``sh_analysis_vector``
  !@
  !@ Returns
  !@ -------
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Reconstructed vector components in East and North directions
  subroutine sh_synthesis_vector(u,v, nx,ny,nz, br,bi,cr,ci)
    real(kind=nr), intent(in) :: br(nz,ny,ny), bi(nz,ny,ny), &
            &                     cr(nz,ny,ny), ci(nz,ny,ny)
    integer(kind=ni), intent(in) :: nx,ny,nz
    ! TODO: I tried to integrate the return into two complex arrays, but that
    !       leads to mysterious SEGFAULTS type "double free or corruption"
    real(kind=nr), intent(out) :: u(nz,ny,nx), v(nz,ny,nx)
    !
    !real(kind=nr) :: br(nz,ny,ny), bi(nz,ny,ny), cr(nz,ny,ny), ci(nz,ny,ny)
    integer(kind=ni) :: ierror, k
    !
    integer(kind=8_ni) :: lwork, ldwork    ! lwork might be larger than 32bit!
    real(kind=nr), allocatable :: work(:), dwork(:)
    !
    !f2py depend(nz,ny) bi, cr, ci, u, v
    ! -----------------------------------------------------------------
    !
    if ( nx /= last_nx .or. ny /= last_ny ) then
       call reset()
       last_nx = nx
       last_ny = ny
    end if
    !
    lwork = max( 4_ni*(ny+1_ni)**2_ni , (2_ni*nz+1_ni)*ny*nx )
    ldwork = (ny + 1_ni)*2_ni
    allocate(work(lwork), dwork(ldwork))
    !
    ! No precalculated wvhses available
    if ( lvhses == 0_ni ) then
       lvhses = ((ny+1_ni)**3_ni) + nx + 15_ni
       allocate(wvhses(lvhses))
       call vhsesi(ny,nx,wvhses,lvhses,work,lwork,dwork,ldwork,ierror)
       if ( ierror /= 0_ni ) then
          write(*,*) 'Error in vhsesi, code', ierror
          stop 1
       end if
    end if
    !
    u(:,:,:) = 0.0_nr
    v(:,:,:) = 0.0_nr
    do k = 1_ni,nz
       call vhses(ny,nx,0_ni,1_ni, v(k,:,:),u(k,:,:), ny,nx, &
               & br(k,:,:),bi(k,:,:),cr(k,:,:),ci(k,:,:), ny,ny, wvhses,lvhses,work,lwork, ierror)
    end do
    if ( ierror /= 0_ni ) then
       write(*,*) 'Error in vhses, code', ierror
       stop 1
    end if
    !
    v(:,:,:) = -v(:,:,:)
    !
    deallocate(work, dwork)
    !
  end subroutine
  !
  !@ Triangular truncation in spectral space
  !@
  !@ Routine works in-place!
  !@
  !@ Parameters
  !@ ----------
  !@ 
  !@ T : int
  !@     Truncation wave number
  !@ a, b : np.ndarrays with shape (nz,ny,ny) and dtype float64
  !@     Spherical harmonics coefficients as returned by ``sh_analysis``
  !@
  !@ Returns
  !@ -------
  !@ 2-tuple of np.ndarray with shape (nz,ny,nx) and dtype float64
  !@     Truncated spherical harmonics coefficients
  subroutine trunc_triangular(ny,nz, a, b, T)
    real(kind=nr), intent(inout) :: a(nz,ny,ny), b(nz,ny,ny)
    integer(kind=ni), intent(in) :: ny,nz, T
    ! TODO: I tried to integrate the return into two complex arrays, but that
    !       leads to mysterious SEGFAULTS type "double free or corruption"
    integer(kind=ni) :: n,m,k
    !f2py depend(nz,ny) b
    ! -----------------------------------------------------------------
    !
    do n = T+2_ni,ny
       do m = 1_ni,n
          do k = 1_ni,nz
             a(k,m,n) = 0.0_nr
             b(k,m,n) = 0.0_nr
          end do
       end do
    end do
    !
  end subroutine
  !
  !@ Reset the internal precalculated work arrays
  !@
  !@ This routine is called automatically as soon as the one of the
  !@ analysis/synthesis routines recieves a grid of a different size. 
  !@ It is not neccessary to call this routine from outside this 
  !@ module.
  subroutine reset()
    ! -----------------------------------------------------------------
    !
    last_nx = 0_ni
    last_ny = 0_ni
    !
    lshaes=0_ni
    lshses=0_ni
    lvhaes=0_ni
    lvhses=0_ni
    !
    if ( allocated(wshaes) ) deallocate(wshaes)
    if ( allocated(wshses) ) deallocate(wshses)
    if ( allocated(wvhaes) ) deallocate(wvhaes)
    if ( allocated(wvhses) ) deallocate(wvhses)
  end subroutine
end module
