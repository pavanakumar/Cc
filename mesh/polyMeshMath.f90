!> The polyhedral mesh math functions module
module PolyMeshMath
  use iso_c_binding
  use PolyMeshConstant
  implicit none

 contains

  !> Cross product of vector of length 3
  pure function cross_prod(x, y) result(ret)
    implicit none
    real(kind=8),dimension(dim_),intent(in) :: x, y
    real(kind=8),dimension(dim_)            :: ret
    integer, dimension(dim_), parameter     :: c1 = (/2,3,1/), &
                                               c2 = (/3,1,2/)
    ret = x(c1) * y(c2) - y(c1) * x(c2)
  end function cross_prod

  !> Ruled face position vector
  pure function r_ruled(r, xi, eta) result(ret)
    implicit none
    real(kind=8), intent(in) :: r(dim_, quad_)
    real(kind=8), intent(in) :: xi, eta
    real(kind=8)             :: ret(dim_)
    ret = ( r(:, 1) * ( 1.0d0 - xi ) + r(:, 2) * xi ) * ( 1.0d0 - eta ) + &
          ( r(:, 4) * ( 1.0d0 - xi ) + r(:, 3) * xi ) * eta
  end function r_ruled

  !> Ruled face area element
  pure function dS_ruled(r, xi, eta) result(ret)
    implicit none
    real(kind=8), intent(in)      :: r(dim_, quad_)
    real(kind=8), intent(in)      :: eta, xi
    real(kind=8), dimension(dim_) :: ret, temp1, temp2
    temp1   = (r(:,2) - r(:,1)) * ( 1.0d0 - eta ) + (r(:,3) - r(:,4)) * eta
    temp2   = (r(:,4) - r(:,1)) * ( 1.0d0 - xi  ) + (r(:,3) - r(:,2)) * xi
    ret     = cross_prod( temp1, temp2 )
  end function dS_ruled

  !> The mesh metrics inner loop
  subroutine poly_metrics_loop( &
    iface, nnode, nface, ninternalface, ncell, &
    x, facelr, facenode, dn, fs, fc, fvol, fvolc )
    implicit none
    !> Function arguements
    integer, intent(in)         :: iface, nnode, nface, ninternalface, ncell
    real(kind=8), intent(in)    :: x(nnode, dim_)
    integer, intent(in)         :: facelr(nface, lr_), facenode(nface, lr_)
    real(kind=8), intent(out)   :: dn(nface, dim_), fs(nface), fc(nface, dim_)
    real(kind=8), intent(out)   :: fvol, fvolc(dim_)
    !> Local variables
    integer                     :: il, ir, iqdr
    real(kind=8)                :: rpos(dim_, quad_), & ! Position vectors
                                   rqdr( dim_),       & ! Quadrature position vector
                                   dSqdr(dim_),       & ! Quadrature dS
                                   magdSqdr             ! |dS|
    !> Init face values to zero
    dn(iface, :)  = 0.0d0
    fc(iface, :)  = 0.0d0
    fvol          = 0.0d0
    fvolc         = 0.0d0
    !> Face metric calculation
    rpos          = x(facenode( face1_:face4_, iface ), :)
    il            = facelr(iface, :)
    ir            = facelr(iface, :)
    !> Obtain the quadrature point values
    do iqdr = 1, 4
      rqdr        = r_ruled( rpos, qxi_(iqdr), qeta_(iqdr) )
      dSqdr       = dS_ruled( rpos, qxi_(iqdr), qeta_(iqdr) )
      magdSqdr    = sqrt( sum( dSqdr * dSqdr ) )
      dn(iface, :) = dn(iface, :) + dSqdr
      fc(iface, :) = fc(iface, :) + rqdr * magdSqdr
      fvol        = fvol  + sum( rqdr * dSqdr )
      fvolc       = fvolc + sum( rqdr * rqdr ) * dSqdr
    end do
    !> Multiply by quadrature weights and scale factor
    dn(iface, :)   = 0.250d0  * dn(iface, :)
    fs(iface)     = sqrt( sum( dn(iface, :) * dn(iface, :) ) )
    dn(iface, :)   = dn(iface, :) / fs(iface)
    fc(iface, :)   = 0.250d0 * fc(:,iface) / fs(iface)

 end subroutine poly_metrics_loop

 
  !> Calculate the mesh metrics
  !> cv : Cell volume   (scalar)
  !> cc : Cell centroid (vector*3)
  !> dn : Face normal   (vector*3)
  !> fs : Face area     (scalar)
  !> fc : Face centroid (vector*3)
  subroutine poly_metrics( &
    nnode, nface, ninternalface, ncell, &
    x, facelr, facenode, cv, cc, dn, fs, fc )
    implicit none
    integer, intent(in)         :: nnode, nface, ninternalface, ncell
    real(kind=8), intent(in)    :: x(nnode, dim_)
    integer, intent(in)         :: facelr(nface, lr_), facenode(nface, quad_)
    real(kind=8), intent(out)   :: cv(ncell), cc(ncell, dim_)
    real(kind=8), intent(out)   :: dn(nface, dim_), fs(nface), fc(nface, dim_)
    !> Local variables
    integer                     :: iface, icell, il, ir
    real(kind=8)                :: fvol, fvolc(dim_)
    !> Zero out everything
    cv = 0.0d0; cc = 0.0d0
    !> Internal face loop for metric calculation
    do iface = 1, ninternalface
      call poly_metrics_loop( iface, nnode, nface, ninternalface, ncell, &
                              x, facelr, facenode, dn, fs, fc, fvol, fvolc )
      il = facelr(iface, lcell_)
      ir = facelr(iface, rcell_)
      !> Add CC/CV contribution of internal face
      cv(il)     = cv(il)    + fvol
      cv(ir)     = cv(ir)    - fvol
      cc(il, :)  = cc(il, :) + fvolc
      cc(ir, :)  = cc(ir, :) - fvolc
    end do
    !> Boundary face loop for metric calculation
    do iface = ninternalface + 1, nface
      call poly_metrics_loop( iface, nnode, nface, ninternalface, ncell, &
                              x, facelr, facenode, dn, fs, fc, fvol, fvolc )
     il = facelr(iface, lcell_)
     !> Add CC/CV contribution of boundary face
      cv(il)       = cv(il)   + fvol
      cc(il, :)    = cc(il, :) + fvolc
    end do
    !> Do cell summation
    do icell = 1, ncell
      cv(icell)    = oneby12_ * cv(icell)
      cc(icell, :) = 0.1250d0 * cc(icell, :) / cv(icell)
    end do
  end subroutine poly_metrics


end module PolyMeshMath

