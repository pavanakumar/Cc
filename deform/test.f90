program test
  use Spring
  implicit none
  integer, parameter :: nnode=5, &
                        nedge=8
  integer :: edgenode(2, nedge) = &
  reshape( &
  (/ 1, 2, &
     1, 3, &
     1, 4, &
     1, 5, &
     2, 3, &
     3, 4, &
     4, 5, &
     5, 2 /), (/ 2, nedge /) )
  integer :: iedge, ia(nnode + 1), ja( 2 * nedge + nnode ), i
  integer :: ik( 3 * nnode + 1), jk( 9 * ( 2 * nedge + nnode ) )
  real(kind=8) :: x(3, nnode) = &
  reshape( &
  (/  0.0d0,  0.0d0, 0.0d0, &
      1.0d0,  1.0d0, 0.0d0, &
      1.0d0, -1.0d0, 0.0d0, &
     -1.0d0, -1.0d0, 0.0d0, &
     -1.0d0,  1.0d0, 0.0d0 /), &
  (/ 3, nnode /) ), &
  k( 3, 3, 2 * nedge + nnode )

!  do iedge = 1, nedge
!    write(*,*) edgenode(:, iedge)
!  end do

  !call edge_to_csr( nnode, nedge, edgenode, ia, ja )
  call kmatrix( nnode, x, nedge, edgenode, ik, jk, k )

!  write(*,*) "ia = ", ik
!  write(*,*) "ja = ", jk

!  do iedge = 1, nnode
!    write(*,*) ik((iedge - 1) * 3 + 1), ik((iedge - 1) * 3 + 2), ik((iedge - 1) * 3 + 3)
!  end do

!  do iedge = 1, nnode
!    do i = ik( (iedge - 1) * 9 + 1 ), ik( iedge*9 )
!      write(*,*) jk( ik( (iedge - 1) * 3 + 1 ) : ik( iedge*3) -1 )
!    end do
!  end do


!  do iedge = 1, 2 * nedge + nnode
!    do i = 1, 3
!      write(*,*) k(:,i,iedge)
!    end do
!    write(*,*)
!  end do
 
end program test
