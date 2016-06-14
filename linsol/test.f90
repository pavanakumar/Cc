module wrap
  implicit none
  interface

  subroutine bsrcsr (n, m, na, a, ja, ia, jao, iao)
    implicit none
    integer (kind=4) :: n, m, na, ia(*), ja(*), jao(*), iao(*)
    real(kind=8) :: a(*)

  end subroutine bsrcsr


  subroutine bsrcsrnew (n, m, nzb, a, ja, ia, jao, iao)
    implicit none
    integer(kind=4) :: n, m, nzb, &
                       ia(n + 1), &
                       ja(nzb), &
                       jao(m * m * nzb), &
                       iao(m * n + 1)
    real(kind=8) :: a(*)
  end subroutine bsrcsrnew

  end interface
end module wrap


subroutine rawprint( n, a )
  implicit none
  integer :: n 
  real(kind=8) :: a(n)
  write(*,*) a
end subroutine rawprint


program test
  implicit none
  integer, parameter :: n = 3, m = 2, nnz = 5
  real(kind=8), dimension( m, m, nnz) :: & 
   a = reshape( (/ &
   111,   112,   121,  122, &
   211,   212,   221,  222, &
   311,   312,   321,  322, &
   411,   412,   421,  422, &
   511,   512,   521,  522 /), &
   (/ m, m, nnz /) )
 
!   a = reshape( (/ &
!   1,   3,   9,   11,   17, &
!   5,   7,   13,  15,   22, &
!   2,   4,   10,  12,   18, &
!   6,   8,   14,  16,   23 /), &
!   (/ m, m, 5 /) )
 
!  real(kind=8) :: ao( m * m * 5)
  integer :: ja(nnz) = reshape( (/ 1, 3, 2, 3, 1 /), (/ nnz /) ), &
             ia(n + 1) = reshape( (/ 1, 3, 5, 6 /), (/ n + 1 /) )
  integer :: jao( nnz * m * m ), iao( n * m + 1 ), nzb, i, j

  nzb = ia(n + 1) - 1

  call rawprint( m * m * nnz, a )
  call bsrcsr( n, m, nzb, a, ja, ia, jao, iao)
!  call bsrcsrnew( n, m, nzb, a, ja, ia, jao, iao )
  call rawprint( m * m * nnz, a )

end program test
