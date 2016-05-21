!! \file spring.f90
!! Spring analogy based mesh deformation
module Spring
  implicit none

  contains
  
  subroutine kmatrix( nnode, x, nedge, edgenode, ik, jk, k )
    implicit none
    integer, intent(in)       :: nnode
    real(kind=8), intent(in)  :: x(3, nnode)
    integer, intent(in)       :: nedge
    integer, intent(in)       :: edgenode( nedge )
    integer, intent(out)      :: ik( 3 * nnode + 1 )
    integer, intent(out)      :: jk( 9 * ( 2 * nedge + nnode ) )
    real(kind=8), intent(out) :: k(3, 3, 2 * nedge + nnode )
    ! Local variables
    integer :: i, ii, j, jj, icount
    integer :: ikb( nnode + 1 ), jkb( 2 * nedge + nnode )
    real(kind=8) :: invr, rij(3), kij(3,3)

    k = 0.0d0
    call edge_to_csr( nnode, nedge, edgenode, ikb, jkb )
    do i = 1, nnode
      do jj = ikb(i) + 1, ikb(i+1) - 1 ! first element is always the diagonal entry
        ! Construct the kij matrix using edges
        j     = jkb(jj)
        rij   = x(:,j) - x(:,i)
        invr  = 1.0d0 / sqrt( sum( rij * rij ) )
        rij   = rij * invr
        do ii = 1, 3
          kij(:, ii) = invr * rij * rij(ii)
        end do
        ! muliply and accumuate the left/right node
        k(:,:,jj)     = - kij
        k(:,:,ikb(i)) = k(:,:,ikb(i)) + kij ! Diagonal entry
        do icount = 1, 3
          write(*,*) kij(:,icount)
        end do
      end do
    end do
    !call bsrcsrnew( nnode, 9, 2 * nedge + nnode, k, jkb, ikb, jk, ik )
 
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

  !> 
  subroutine edge_to_csr( nnode, nedge, edgenode, ia, ja )
    implicit none
    integer, intent(in)         :: nnode, nedge
    integer, intent(in)         :: edgenode(2, nedge)
    integer, intent(out)        :: ia( nnode + 1 ), ja( 2 * nedge + nnode )

    !> Local subroutine variables
    integer :: ierr, i
    ia = 0
    ja = -1
    !> First pass form the sizes
    do i = 1, nedge
      ia(edgenode(1, i) + 1) = ia(edgenode(1, i) + 1 ) + 1
      ia(edgenode(2, i) + 1) = ia(edgenode(2, i) + 1 ) + 1
    end do
    ia    = ia + 1
    ia(1) = 0
    !> Form the xadj offsets (Can replace with intrisic?)
    do i = 1, nnode
      ia(i + 1) = ia(i) + ia(i + 1)
    end do
    !> Check size match
    if( ia(nnode + 1) .ne. (2 * nedge + nnode) ) then
      write(*,*) 'Error: Sizes of ia do not match 2 * nedge', &
        ia(nnode + 1), 2 * nedge + nnode
    end if
    !> Add all diagonal entries as the first entry
    do i = 1, nnode
      ia(i)     = ia(i) + 1
      ja(ia(i)) = i
    end do
    !> Second pass form the adjncy
    do i = 1, nedge
      !> Push the left face data
      ia(edgenode(1, i))     = ia(edgenode(1, i)) + 1
      ja(ia(edgenode(1, i))) = edgenode(2, i)
      !> Push the right face data
      ia(edgenode(2, i))     = ia(edgenode(2, i)) + 1
      ja(ia(edgenode(2, i))) = edgenode(1, i)
    end do
    ia = cshift(ia, nnode)
    !!! Check to see if the adjncy array is formed correctly
    if( count( ja .lt. 0 ) .ne. 0 ) then
      write(*,*) 'Error: Adjncy array from edge2nodes not constructed correctly'
      stop
    end if
    !!! Ensure the numbering is FORTRAN style (1-indexing) :(
    ia    = ia + 1
    ia(1) = 1

  end subroutine edge_to_csr


end module Spring


