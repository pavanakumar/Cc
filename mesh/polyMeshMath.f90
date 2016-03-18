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
  pure function r_ruled(r1, r2, r3, r4, xi, eta) result(ret)
    implicit none
    real(kind=8),dimension(dim_),intent(in) :: r1, r2, r3, r4
    real(kind=8), intent(in)                :: eta, xi
    real(kind=8),dimension(dim_)            :: ret
    ret = ( r1 * ( 1.0d0 - xi ) + r2 * xi ) * ( 1.0d0 - eta ) + &
          ( r4 * ( 1.0d0 - xi ) + r3 * xi ) * eta
  end function r_ruled

  !> Ruled face area element
  pure function dS_ruled(r1, r2, r3, r4, xi, eta) result(ret)
    implicit none
    real(kind=8),dimension(dim_),intent(in) :: r1, r2, r3, r4
    real(kind=8), intent(in)                :: eta, xi
    real(kind=8),dimension(dim_)            :: ret, temp1, temp2
    temp1   = (r2 - r1) * ( 1.0d0 - eta ) + (r3 - r4) * eta
    temp2   = (r4 - r1) * ( 1.0d0 - xi  ) + (r3 - r2) * xi
    ret     = cross_prod( temp1, temp2 )
  end function dS_ruled

  !> Ruled face contribution to centroid
  pure function cc_ruled(r1, r2, r3, r4, xi, eta) result(ret)
    implicit none
    real(kind=8),dimension(dim_),intent(in) :: r1, r2, r3, r4
    real(kind=8), intent(in)                :: eta, xi
    real(kind=8),dimension(dim_)            :: temp, ret
    real(kind=8)                            :: temp_dot
    temp      = r_ruled(r1, r2, r3, r4, xi, eta)
    temp_dot  = sum( temp * temp )
    ret       = temp_dot * dS_ruled(r1, r2, r3, r4, xi, eta)
  end function cc_ruled

  !> The mesh metrics inner loop
  subroutine poly_metrics_loop( &
    iface, nnode, nface, ninternalface, ncell, &
    x, facelr, facenode, dn, fs, fc, fvol, fvolc )
    implicit none
    !> Function arguements
    integer, intent(in)         :: iface, nnode, nface, ninternalface, ncell
    real(kind=8), intent(in)    :: x(dim_, nnode)
    integer, intent(in)         :: facelr(lr_, nface), facenode(quadp1_, nface)
    real(kind=8), intent(out)   :: dn(dim_, nface), fs(nface), fc(dim_, nface)
    real(kind=8), intent(out)   :: fvol, fvolc(dim_)
    !> Local variables
    integer                     :: il, ir
    real(kind=8)                :: r1(dim_), r2(dim_), r3(dim_), r4(dim_)
    !> Face metric calculation
    r1 = x(:, facenode( 2, iface ))
    r2 = x(:, facenode( 3, iface ))
    r3 = x(:, facenode( 4, iface ))
    r4 = x(:, facenode( 5, iface ))
    il = facelr( lcell_, iface)
    ir = facelr( rcell_, iface)
    !> Face centroid
    fc(:, iface) = 0.250d0 * ( r1 + r2 + r3 + r4 )
    !> Face normal and area
    dn(:, iface) = 0.50d0 * cross_prod( (r3 - r1), (r4 - r2) )
    fs(iface)    = sqrt( sum( dn(:, iface) * dn(:, iface) ) )
    dn(:, iface) = dn(:, iface) / fs(iface)
    !> Face volume contribution
    fvol         = oneby12_ * sum( ( r1 + r2 )  * cross_prod(  (r2 + r3), (r3 + r4) ) )
    !> Face cell centroid contribution (Gauss quadrture)
    fvolc        = 0.250d0 * ( cc_ruled(r1, r2, r3, r4, xi1_, eta1_) + &
                               cc_ruled(r1, r2, r3, r4, xi1_, eta2_) + &
                               cc_ruled(r1, r2, r3, r4, xi2_, eta1_) + &
                               cc_ruled(r1, r2, r3, r4, xi2_, eta2_) ) 
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
    real(kind=8), intent(in)    :: x(dim_, nnode)
    integer, intent(in)         :: facelr(lr_, nface), facenode(quadp1_, nface)
    real(kind=8), intent(out)   :: cv(ncell), cc(dim_, ncell)
    real(kind=8), intent(out)   :: dn(dim_, nface), fs(nface), fc(dim_, nface)
    !> Local variables
    integer                     :: iface, icell, il, ir
    real(kind=8)                :: fvol, fvolc(dim_)
    !> Zero out everything
    cv = 0.0d0; cc = 0.0d0
    !> Internal face loop for metric calculation
    do iface = 1, ninternalface
      call poly_metrics_loop( iface, nnode, nface, ninternalface, ncell, &
                              x, facelr, facenode, dn, fs, fc, fvol, fvolc )
      il = facelr( lcell_, iface)
      ir = facelr( rcell_, iface)
      !> Add CC/CV contribution of internal face
      cv(il)       = cv(il)   + fvol
      cv(ir)       = cv(ir)   - fvol
      cc(:,il)     = cc(:,il) + fvolc
      cc(:,ir)     = cc(:,ir) - fvolc
    end do
    !> Boundary face loop for metric calculation
    do iface = ninternalface + 1, nface
      call poly_metrics_loop( iface, nnode, nface, ninternalface, ncell, &
                              x, facelr, facenode, dn, fs, fc, fvol, fvolc )
     il = facelr( lcell_, iface)
     !> Add CC/CV contribution of boundary face
      cv(il)       = cv(il)   + fvol
      cc(:,il)     = cc(:,il) + fvolc
    end do
    !> Do cell summation
    do icell = 1, ncell
      cc(:,icell)  = 0.50d0 * cc(:,icell) / cv(icell)
    end do
  end subroutine poly_metrics


end module PolyMeshMath

