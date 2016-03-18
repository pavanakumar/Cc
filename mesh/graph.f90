!>
module Graph  
  implicit none
  !>
  type :: crsGraph
    integer                                 :: nvtx = 0, nedge = 0
    integer, dimension(:), allocatable      :: xadj, adjncy
    real(kind=8), dimension(:), allocatable :: adjwgt, vvol, vsurf
  end type crsGraph

  contains

  !>
  subroutine allocate_graph( nnode, nedge, gr )
    implicit none
    integer, intent(in)         :: nnode, nedge
    type(crsGraph), intent(out) :: gr
    !> Some simple memory checks
    if( allocated(gr%xadj)   .eqv. .true. ) deallocate(gr%xadj)
    if( allocated(gr%adjncy) .eqv. .true. ) deallocate(gr%adjncy)
    if( allocated(gr%adjwgt) .eqv. .true. ) deallocate(gr%adjwgt)
    if( allocated(gr%vsurf)  .eqv. .true. ) deallocate(gr%vsurf)
    !> Allocations
    gr%nvtx = nnode
    gr%nedge = nedge
    allocate( gr%xadj(gr%nvtx + 1),    &
              gr%adjncy(2 * gr%nedge), &
              gr%adjwgt(2 * gr%nedge), &
              gr%vsurf(gr%nvtx),       &
              gr%vvol(gr%nvtx)         )
  end subroutine allocate_graph

  !>
  subroutine allocate_graph_partial( nnode, nedge, gr )
    implicit none
    integer, intent(in)         :: nnode, nedge
    type(crsGraph), intent(out) :: gr
    !> Some simple memory checks
    if( allocated(gr%xadj)   .eqv. .true. ) deallocate(gr%xadj)
    if( allocated(gr%adjncy) .eqv. .true. ) deallocate(gr%adjncy)
    if( allocated(gr%adjwgt) .eqv. .true. ) deallocate(gr%adjwgt)
    if( allocated(gr%vsurf)  .eqv. .true. ) deallocate(gr%vsurf)
    !> Allocations
    gr%nvtx  = nnode
    gr%nedge = nedge
    allocate( gr%xadj(gr%nvtx + 1), &
              gr%adjncy(2 * gr%nedge) )
  end subroutine allocate_graph_partial

  !>
  subroutine free_graph(gr)
    implicit none
    type(crsGraph), intent(inout) :: gr
    if( allocated(gr%xadj)   .eqv. .true. ) deallocate(gr%xadj)
    if( allocated(gr%adjncy) .eqv. .true. ) deallocate(gr%adjncy)
    if( allocated(gr%adjwgt) .eqv. .true. ) deallocate(gr%adjwgt)
    if( allocated(gr%vsurf)  .eqv. .true. ) deallocate(gr%vsurf)
  end subroutine free_graph

  !> 
  subroutine edge_to_graph( nnode, nedge, edgenode, gr )
    implicit none
    integer, intent(in)         :: nnode, nedge
    integer, intent(in)         :: edgenode(2, nedge)
    type(crsGraph), intent(inout) :: gr
    !> Local subroutine variables
    integer :: ierr, i
    gr%xadj   = 0
    gr%adjncy = -1
    !> First pass form the sizes
    do i = 1, nedge
      gr%xadj(edgenode(1, i)) = gr%xadj(edgenode(1, i)) + 1
      gr%xadj(edgenode(2, i)) = gr%xadj(edgenode(2, i)) + 1
    end do
    gr%xadj = cshift(gr%xadj, gr%nvtx)
    gr%xadj(1) = 0
    !> Form the xadj offsets (Can replace with intrisic?)
    do i = 1, gr%nvtx
      gr%xadj(i + 1) = gr%xadj(i) + gr%xadj(i + 1)
    end do
    !> Check size match
    if( gr%xadj(gr%nvtx + 1) .ne. (2 * nedge) ) then
      write(*,*) 'Error: Sizes of xadj do not match 2 * nedge', &
        gr%xadj(gr%nvtx + 1), 2 * nedge
    end if
    !> Second pass form the adjncy
    do i = 1, nedge
      !> Push the left face data
      gr%xadj(edgenode(1, i))            = gr%xadj(edgenode(1, i)) + 1
      gr%adjncy(gr%xadj(edgenode(1, i))) = edgenode(2, i) - 1
      !> Push the right face data
      gr%xadj(edgenode(2, i))            = gr%xadj(edgenode(2, i)) + 1
      gr%adjncy(gr%xadj(edgenode(2, i))) = edgenode(1, i) - 1
    end do
    gr%xadj    = cshift(gr%xadj, gr%nvtx)
    gr%xadj(1) = 0
    !> Check to see if the adjncy array is formed correctly
    if(count(gr%adjncy .lt. 0) .ne. 0) then
      write(*,*) 'Error: Adjncy array from edgenode not constructed correctly'
    end if

  end subroutine edge_to_graph

  !> 
  subroutine edge_data_to_graph( nnode, nedge, nbc, edgenode, &
                                 nodwgt, edgwgt, gr )
    implicit none
    integer, intent(in)           :: nnode, nedge, nbc
    integer, intent(in)           :: edgenode(2, nedge + nbc)
    real(kind=8), intent(in)      :: nodwgt(nnode), edgwgt(nedge + nbc)
    type(crsGraph), intent(inout) :: gr
    !> Local variables
    integer :: i
    !> Init values to zero
    gr%vsurf  = 0.0d0
    gr%vvol   = nodwgt
    gr%adjwgt = 0.0d0
    !> Boundary node weights vsurf
    do i = nedge + 1, nedge + nbc
      gr%vsurf(edgenode(1, i)) = gr%vsurf(edgenode(1, i)) + edgwgt(i)
    end do
    !> Add edge weights
    do i = 1, nedge
      !> Push the left face data
      gr%xadj(edgenode(1, i))            = gr%xadj(edgenode(1, i)) + 1
      gr%adjwgt(gr%xadj(edgenode(1, i))) = edgwgt(i)
      !> Push the right face data
      gr%xadj(edgenode(2, i))            = gr%xadj(edgenode(2, i)) + 1
      gr%adjwgt(gr%xadj(edgenode(2, i))) = edgwgt(i)
    end do
    gr%xadj = cshift(gr%xadj, gr%nvtx)
    gr%xadj(1) = 0
 
  end subroutine edge_data_to_graph

  !>
  subroutine graph_mg( gr, nmgpart, mgpart )
    use Wrap
    implicit none
    type(crsGraph), intent(inout) :: gr
    integer, intent(out)          :: nmgpart, mgpart( gr%nvtx )
    !> Local variables
    integer :: minsize = 1, maxsize = 10, &
               options(4) = (/ 4, 6, 128, 3 /), nmoves
    !> Call the external MG interface
    call MGridGen_f90 ( &
        gr%nvtx, gr%xadj, gr%vvol, gr%vsurf, gr%adjncy, &
        gr%adjwgt, minsize, maxsize, options, nmoves, &
        nmgpart, mgpart )
    call fix_mg_degeneracy( gr, nmgpart, mgpart )
  end subroutine graph_mg


  !> This routine checks if the MG parts are
  !> consistent and every node in the graph
  !> has more than one edge neighbour. Otherwise
  !> the node is merged with one of its edge neighbour
  !> and the part array is re-generated
  subroutine fix_mg_degeneracy( gr, nmgpart, mgpart )
    implicit none
    type(crsGraph), intent(inout) :: gr
    integer, intent(inout)        :: nmgpart, mgpart( gr%nvtx )
 
  end subroutine fix_mg_degeneracy

end module Graph

