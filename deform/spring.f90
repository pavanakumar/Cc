!! \file spring.f90
!! Spring analogy based mesh deformation
module Spring
  implicit none

  contains
  
  subroutine kmatrix( nnode, x, nedge, edgenode, ik, jk, kmat )
    implicit none
    integer, intent(in)       :: nnode
    real(kind=8), intent(in)  :: x(3, nnode)
    integer, intent(in)       :: nedge
    integer, intent(in)       :: edgenode( nedge )
    integer, intent(out)      :: ik( 3 * nnode + 1 )
    integer, intent(out)      :: jk( 9 * ( 2 * nedge + nnode ) )
    real(kind=8), intent(out) :: kmat(3, 3, 2 * nedge + nnode )

    ! Local variables
    integer :: i, ii, j, jj, icount
    integer :: ikb( nnode + 1 ), jkb( 2 * nedge + nnode )
    real(kind=8) :: invr, rij(3), kij(3,3)

    kmat = 0.0d0
    call edge_to_csr( nnode, nedge, edgenode, ikb, jkb )
    do i = 1, nnode
      do jj = ikb(i) + 1, ikb(i+1) - 1 ! first element is always the diagonal entry
        ! Construct the kij matrix using edges
        j     = jkb(jj)
        rij   = x(:, j) - x(:, i)
        invr  = 1.0d0 / sqrt( sum( rij * rij ) )
        rij   = rij * invr
        do ii = 1, 3
          kij(:, ii) = invr * rij * rij(ii)
        end do
        ! muliply and accumuate the left/right node
        kmat(:, :, jj)     = - kij
        kmat(:, :, ikb(i)) = kmat(:, :, ikb(i)) + kij ! Diagonal entry
      end do
    end do
    !write(*,*) 'ia =', ikb
    !write(*,*) 'ja = '
    !do i = 1, nnode
    !  write(*,*) jkb( ikb(i) : ikb(i+1) - 1 )
    !end do
    call bsrcsrip( nnode, 3, 2 * nedge + nnode, kmat, jkb, ikb, jk, ik )
    !write(*,*) 'ia =', ik
    !write(*,*) 'ja = '
    !do i = 1, nnode * 3
    !  write(*,*) jk( ik(i) : ik(i+1) - 1 )
    !end do
 
  end subroutine kmatrix

  subroutine kx( nnode, x, nedge, edgenode, v )
    implicit none
    integer, intent(in)         :: nnode
    real(kind=8), intent(in)    :: x(3, nnode)
    integer, intent(in)         :: nedge
    integer, intent(in)         :: edgenode( 2, nedge )
    real(kind=8), intent(inout) :: v( 3, nnode )
    ! Local variables
  end subroutine kx

  subroutine pilu( nnode, nedge, ip, jp, p, v )
    implicit none
    integer, intent(in)        :: nnode, nedge
    integer, intent(in)        :: ip( 3 * nnode + 1 )
    integer, intent(in)        :: jp( 9 * ( 2 * nedge + nnode ) )
    real(kind=8), intent(in)   :: p( 3, 3, 2 * nedge + nnode )
    real(kind=8), intent(out)  :: v( 3, nnode )
    ! Sparse matrix multiply dx = p * dx

  end subroutine pilu

end module Spring


