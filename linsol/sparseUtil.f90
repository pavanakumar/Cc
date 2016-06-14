! Copy raw array data ignoring the dimension
! of the array. Use with caution
subroutine rawcp( n, ain, aout )
  implicit none
  integer (kind=4), intent(in) :: n
  real(kind=8), intent(in)     :: ain(*)
  real(kind=8), intent(out)    :: aout(*)

  aout(1:n) = ain(1:n)

end subroutine rawcp

! An not so efficient in-place Block Sparse Row storage
! to Compressed Row storage matrix
subroutine bsrcsrcopy (n, m, nzb, a, ja, ia, jao, iao, ao)
  implicit none
  integer(kind=4) :: n, m, nzb, &
                     ia(n + 1), &
                     ja(nzb), &
                     jao(m * m * nzb), &
                     iao(m * n + 1)
  real(kind=8) :: a(m * m, nzb)
! Local
  integer(kind=4) :: i, istart, iend, &
                     ij, ii, irow, j, &
                     jstart, k, krow
  real(kind=8), allocatable :: ao(:)
  allocate( ao(m * m * nzb) )

  irow = 1
  krow = 1
  iao(irow) = 1      

  do ii = 1, n
    istart = ia(ii)
    iend   = ia(ii + 1) - 1
! Go row-by-row and copy contents to temp ao and then
! replace the block rows of a by ao after finishing
    do i = 1, m ! Each row in block row
      do k = istart, iend  ! Block row index
         jstart = m * (ja(k) - 1)
         do j = 1, m ! Each column in block
           ij = (j - 1) * m + i
           ao(krow) = a(ij, k)
           jao(krow) = jstart + j
           krow = krow + 1
         end do
       end do
       irow = irow + 1
       iao(irow) = krow
    end do
  end do
! copy contents
  call rawcp( m * m * nzb, ao, a )
  deallocate(ao)

end subroutine bsrcsrcopy


! An not so efficient in-place Block Sparse Row storage
! to Compressed Row storage matrix
subroutine bsrcsr (n, m, nzb, a, ja, ia, jao, iao)
  implicit none
  integer(kind=4) :: n, m, nzb, &
                     ia(n + 1), &
                     ja(nzb), &
                     jao(m * m * nzb), &
                     iao(m * n + 1)
  real(kind=8)    :: a(m * m, nzb)
! Local
  integer(kind=4) :: i, istart, iend, &
                     ij, ii, irow, j, &
                     jstart, k, krow
  real(kind=8), allocatable :: ao(:)
  allocate( ao(m * m * nzb) )

  irow = 1
  krow = 1
  iao(irow) = 1      

  do ii = 1, n
    istart = ia(ii)
    iend   = ia(ii + 1) - 1
! Go row-by-row and copy contents to temp ao and then
! replace the block rows of a by ao after finishing
    do i = 1, m ! Each row in block row
      do k = istart, iend  ! Block row index
         jstart = m * (ja(k) - 1)
         do j = 1, m ! Each column in block
           ij = (j - 1) * m + i
           ao(krow) = a(ij, k)
           jao(krow) = jstart + j
           krow = krow + 1
         end do
       end do
       irow = irow + 1
       iao(irow) = krow
    end do
  end do
! copy contents
  call rawcp( m * m * nzb, ao, a )
  deallocate(ao)

end subroutine bsrcsr


! An efficient in-place Block Sparse Row storage
! to Compressed Row storage matrix
subroutine bsrcsrip (n, m, nzb, a, ja, ia, jao, iao)
  implicit none
  integer(kind=4) :: n, m, nzb, &
                     ia(n + 1), &
                     ja(nzb), &
                     jao(m * m * nzb), &
                     iao(m * n + 1)
  real(kind=8) :: a(m * m, nzb)
! Local
  integer(kind=4) :: i, istart, iend, &
                     ij, ii, irow, j, &
                     jstart, k, krow, &
                     maxcol, kao
  real(kind=8), allocatable :: ao(:)

  krow      = 1
  irow      = 1
  iao(irow) = 1      
! Find the maximum number of non-zero column entries in a row
  maxcol    = m * m * maxval( ia( 2 : n + 1 ) - ia( 1 : n ) )
  allocate( ao(maxcol) )  
! Loop over all block rows
  do ii = 1, n
    istart = ia(ii)
    iend   = ia(ii + 1) - 1
! Go row-by-row and copy contents to ao and then
! replace the block rows of a by ao after finishing
    kao = 1
    do i = 1, m ! Each row in block row
      do k = istart, iend  ! Block row index
         jstart = m * (ja(k) - 1)
         do j = 1, m
           ij = (j - 1) * m + i
           ao(kao)   = a(ij, k)
           jao(krow) = jstart + j
           krow = krow + 1
           kao  = kao  + 1
         end do
       end do
       irow = irow + 1
       iao(irow) = krow
    end do
! Copy back to array a
    kao = 1
    do i = istart, iend
      do k = 1, m * m
        a(k, i) = ao(kao)
        kao = kao + 1
      end do
    end do
  end do

end subroutine bsrcsrip

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


!*****************************************************************************80
!
!! ATX_CR computes A'*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) W(N), the value of A'*X.
!
subroutine atx_cr ( n, nz_num, ia, ja, a, x, w )
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(ja(k1:k2)) = w(ja(k1:k2)) + a(k1:k2) * x(i)
  end do

end subroutine atx_cr 

!*****************************************************************************80
!
!! ATX_ST computes A'*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = 8 ) W(N), the value of A'*X.
!
subroutine atx_st ( n, nz_num, ia, ja, a, x, w )
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nz_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(j) = w(j) + a(k) * x(i)
  end do

end subroutine atx_st

!*****************************************************************************80
!
!! AX_CR computes A*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) W(N), the value of A*X.
!
subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
  end do

end subroutine ax_cr

!*****************************************************************************80
!
!! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!    The array UA can be used to locate the diagonal elements of the matrix.
!
!    It is assumed that every row of the matrix includes a diagonal element,
!    and that the elements of each row have been ascending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!    On output, the order of the entries of JA may have changed because of
!    the sorting.
!
!    Output, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
 implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) ua(n)

  ua(1:n) = -1

  do i = 1, n
    do k = ia(i), ia(i+1) - 1
      if ( ja(k) == i ) then
        ua(i) = k
        exit
      end if
    end do
  end do

end subroutine diagonal_pointer_cr

!*****************************************************************************80
!
!! ILU_CR_INPL computes the incomplete LU factorization of a matrix. 
!                           (In place version)
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!    FORTRAN90 (in-place) by Pavanakumar Mohanamuraly
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!    
!    Input, integer ( kind = 4 ) IW(N), work array to store non-zero entries
!
!    Output, real ( kind = 8 ) A(NZ_NUM), the ILU factorization of A.
!    Note: Original values in A are replaced by the lu factors

subroutine ilu_cr_inpl ( n, nz_num, ia, ja, ua, iw, a )
  implicit none

  integer ( kind = 4 ), intent(in)     :: n
  integer ( kind = 4 ), intent(in)     :: nz_num
  integer ( kind = 4 ), intent(in)     :: ia(n+1)
  integer ( kind = 4 ), intent(in)     :: ja(nz_num)
  integer ( kind = 4 ), intent(inout)  :: ua(n)
  integer ( kind = 4 ), intent(out)    :: iw(n)
  real ( kind = 8 ), intent(inout)     :: a(nz_num)

  ! Local variables
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) jw
  integer ( kind = 4 ) k
  real ( kind = 8 ) tl

  do i = 1, n
!
!  IW points to the nonzero entries in row I.
!
    iw(1:n) = -1

    do k = ia(i), ia(i + 1) - 1
      iw( ja(k) ) = k
    end do

    do j = ia(i), ia(i + 1) - 1
      jrow = ja(j)
      if ( i .le. jrow ) then
        exit
      end if
      tl = a(j) * a( ua(jrow) )
      a(j) = tl
      do jj = ua(jrow) + 1, ia(jrow + 1) - 1
        jw = iw( ja(jj) )
        if ( jw .ne. -1 ) then
          a(jw) = a(jw) - tl * a(jj)
        end if
      end do
    end do

    ua(i) = j

    if ( jrow .ne. i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a)' ) '  JROW ~= I'
      write ( *, '(a,i8)' ) '  JROW = ', jrow
      write ( *, '(a,i8)' ) '  I    = ', i
      stop
    end if

    if ( a(j) .eq. 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
      write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
      stop
    end if

    a(j) = 1.0D+00 / a(j)

  end do

  a( ua(1:n) ) = 1.0D+00 / a( ua(1:n) )

end subroutine ilu_cr_inpl


!*****************************************************************************80
!
!! ILU_CR computes the incomplete LU factorization of a matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
!    Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.
!
subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) iw(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) jw
  integer ( kind = 4 ) k
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) tl
  integer ( kind = 4 ) ua(n)
!
!  Copy A.
!
  l(1:nz_num) = a(1:nz_num)

  do i = 1, n
!
!  IW points to the nonzero entries in row I.
!
    iw(1:n) = -1

    do k = ia(i), ia(i+1) - 1
      iw(ja(k)) = k
    end do

    do j = ia(i), ia(i+1) - 1
      jrow = ja(j)
      if ( i <= jrow ) then
        exit
      end if
      tl = l(j) * l(ua(jrow))
      l(j) = tl
      do jj = ua(jrow) + 1, ia(jrow+1) - 1
        jw = iw(ja(jj))
        if ( jw /= -1 ) then
          l(jw) = l(jw) - tl * l(jj)
        end if
      end do
    end do

    ua(i) = j

    if ( jrow /= i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a)' ) '  JROW ~= I'
      write ( *, '(a,i8)' ) '  JROW = ', jrow
      write ( *, '(a,i8)' ) '  I    = ', i
      stop
    end if

    if ( l(j) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
      write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
      stop
    end if

    l(j) = 1.0D+00 / l(j)

  end do

  l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

end subroutine ilu_cr

!*****************************************************************************80
!
!! LUS_CR applies the incomplete LU preconditioner.
!
!  Discussion:
!
!    The linear system M * Z = R is solved for Z.  M is the incomplete
!    LU preconditioner matrix, and R is a vector supplied by the user.
!    So essentially, we're solving L * U * Z = R.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = 8 ) L(NZ_NUM), the matrix values.
!
!    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
!    of each row.
!
!    Input, real ( kind = 8 ) R(N), the right hand side.
!
!    Output, real ( kind = 8 ) Z(N), the solution of the system M * Z = R.
!
subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) ua(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) z(n)
!
!  Copy R in.
!
  w(1:n) = r(1:n)
!
!  Solve L * w = w where L is unit lower triangular.
!
  do i = 2, n
    do j = ia(i), ua(i) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
  end do
!
!  Solve U * w = w, where U is upper triangular.
!
  do i = n, 1, -1
    do j = ua(i) + 1, ia(i+1) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
    w(i) = w(i) / l(ua(i))
  end do
!
!  Copy Z out.
!
  z(1:n) = w(1:n)

end subroutine lus_cr

