!> Mesh module
module Constants
  real(kind=8), parameter :: oneby3_  = 0.33333333333333333333333333333333333333333d0
  real(kind=8), parameter :: oneby6_  = 0.16666666666666666666666666666666666666666d0
  real(kind=8), parameter :: oneby12_ = 0.08333333333333333333333333333333333333333d0
  !> Quadrature points (for cell-centroid)
  real(kind=8), parameter :: xi1_     = 0.78867513459481288225457439025097872782380d0
  real(kind=8), parameter :: xi2_     = 0.21132486540518711774542560974902127217620d0
  real(kind=8), parameter :: eta1_    = 0.78867513459481288225457439025097872782380d0
  real(kind=8), parameter :: eta2_    = 0.21132486540518711774542560974902127217620d0
end module Constants

!>
module Mesh
  use iso_c_binding
  use Constants
  implicit none
  !> CONSTANTS
  !> Comprises of patchstart, patchsize, patchtype, patchneigh
  integer, parameter :: pstart_ = 1, psize_ = 2, ptype_ = 3, pproc_ = 4
  integer, parameter :: dim_ = 3, lr_ = 2, lcell_ = 1, rcell_ = 2
  integer, parameter :: tri_ = 3, quad_ = 4, quadp1_ = 5
  integer, parameter :: enable_parallel_ = 0
  integer, parameter :: processor_bc_ = 0, wall_bc_ = 1, symmetry_bc_ = 2, &
                        inlet_bc_ = 3, outlet_bc_ = 4, riemann_bc_ = 5, &
                        empty_bc_ = 6
  integer, parameter :: max_cell_node_ = 8, max_cell_face_ = 6, &
                        max_cell_edge_ = 12
  !>
  type :: polyMesh
    !> Sizes
    integer :: nnode = 0, nface = 0, ninternalface = 0, &
               nedge = 0, ninternaledge = 0, &
               ncell = 0, npatch = 0, ilevel = 0, nlevel = 0
    !> Connecivity/Topology
    integer, dimension(:,:), allocatable :: facelr, cellface
    integer, dimension(:,:), allocatable :: facenode
    integer, dimension(:,:), allocatable :: edgenode
    integer, dimension(:,:), pointer     :: patchdata
    integer, dimension(:),   allocatable :: celltype
    integer, dimension(:,:), allocatable :: cellnode, celledge
    !> Metrics
    !> xyz coords, cell-center, unit-normal, face-center
    !> cell-volume, face-area
    real(kind=8), dimension(:,:), pointer     :: x
    real(kind=8), dimension(:,:), allocatable :: cc, dn, fc
    real(kind=8), dimension(:), allocatable   :: cv, fs
    real(kind=8), dimension(3)                :: xmin, xmax
    !> MPI parallel run (global numbering of local entities)
    logical                            :: parallel = .false.
    integer, dimension(:), pointer     :: cellgid, nodegid, facegid
    !> Multigrid variables
    integer, dimension(:), allocatable :: mgpart
  end type polyMesh

  !>
  type :: crsGraph
    integer                                 :: nvtx = 0, npart = 0
    integer, dimension(:), allocatable      :: xadj, adjncy
    real(kind=8), dimension(:), allocatable :: adjwgt, vvol, vsurf
  end type crsGraph

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

  !>
  subroutine create_mg_pm( nlevel, pm, ipar )
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
 
    write(*,*) "Opening reader API and reading mesh ..."
    call reader_of( nlevel, pm, ipar )
    write(*,*) "Done"
    write(*,*) "Constructing MG ..."
    call create_mg_levels( nlevel, pm, ipar )
    write(*,*) "Done"
  end subroutine create_mg_pm

  !> Read the polyMesh data from OpenFOAM reader
  !> uses parmgen and creates the multigrid levels
  subroutine reader_of( nlevel, pm, ipar )
    use Wrap
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
    !> Local variables
    integer :: ilvl, i
    !> Make the meshes aware of the multi-grid
    do ilvl = 1, nlevel
      if( ipar .eq. enable_parallel_ ) then
        pm(ilvl)%parallel = .true.
      end if
      pm(ilvl)%nlevel = nlevel
      pm(ilvl)%ilevel = ilvl
      nullify(pm(ilvl)%x)
    end do
    !> Read in the finest mesh from OpenFOAM
    call init_mesh_api( ipar )
    !> Read all sizes for allocation
    call get_pm_sizes( pm(nlevel)%nnode, pm(nlevel)%nface, &
                       pm(nlevel)%ninternalface, &
                       pm(nlevel)%nedge, pm(nlevel)%ninternaledge, &
                       pm(nlevel)%ncell, pm(nlevel)%npatch )
    !> Allocate the polyMesh
    call allocate_pm( pm(nlevel) )
    !> Read nodes
    call get_pm_nodes( pm(nlevel)%nnode, pm(nlevel)%x )
    !> Read faces
    call get_pm_faces( pm(nlevel)%nface, pm(nlevel)%ninternalface, &
                       pm(nlevel)%facelr, pm(nlevel)%facenode )
    !> Read edges
    call get_pm_edges( pm(nlevel)%nedge, pm(nlevel)%edgenode ) 
    !> Boundary condition data
    call get_pm_patches( pm(nlevel)%npatch, pm(nlevel)%patchdata )
    !> Parallel data
    if( ipar .eq. enable_parallel_ ) then
      call get_cellgid( pm(nlevel)%ncell, pm(nlevel)%cellgid )
      call get_nodegid( pm(nlevel)%nnode, pm(nlevel)%nodegid )
      call get_facegid( pm(nlevel)%nface, pm(nlevel)%facegid )
    end if
    !> Read extra data
!    call get_pm_extra( pm(nlevel)%ncell, pm(nlevel)%celltype, &
!                       pm(nlevel)%cellnode, pm(nlevel)%cellface, &
!                       pm(nlevel)%celledge )
    call mesh_metrics( pm(nlevel) ) 
    !> Close the interface
    call close_mesh_api()
  end subroutine reader_of

 !> Create the Multi-Grid levels
  !> using mgridgen wrapper
  subroutine create_mg_levels( nlevel, pm, ipar )
    use Wrap
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
    !> Local variables
    integer :: ilvl, minsize=1, maxsize=10, &
               options(4) = (/ 4, 6, 128, 3 /), nmoves
    type(crsGraph) :: gr
    real(kind=8) :: sum_fine, sum_coarse
 
    !> Multi-grid mgridgen
    do ilvl = nlevel, 2, -1
      call pm_to_graph( pm(ilvl), gr )
      call MGridGen_f90 &
      ( &
        gr%nvtx, gr%xadj, gr%vvol, gr%vsurf, gr%adjncy, &
        gr%adjwgt, minsize, maxsize, options, nmoves, &
        gr%npart, pm(ilvl)%mgpart &
      )
      !! call fix_mg_degeneracy( pm(ilvl) )
      call build_pm_coarse( pm(ilvl), pm(ilvl - 1), gr )
      call mesh_metrics( pm(ilvl - 1) ) 
      !> Check the multigrid volumes
      sum_fine   = sum(pm(ilvl)%cv)
      sum_coarse = sum(pm(ilvl - 1)%cv)
      if( abs( (sum_fine - sum_coarse) / sum_fine ) .gt. 1.0d-10 ) then
        write(*,*) "Error: In metrics", sum_fine, sum_coarse
      end if
    end do
  end subroutine create_mg_levels

  !>
  subroutine build_pm_coarse( pmf, pmc, gr )
    implicit none
    type(polyMesh) :: pmf, pmc
    type(crsGraph) :: gr
    !> Local variables
    integer, allocatable :: tmplr(:,:)
    integer :: ninternalface, nbface, iface

    !> Assign grid level infomration
    pmc%ilevel = pmf%ilevel - 1
    pmc%nlevel = pmf%nlevel
    !> Find unique internal faces in coarse mesh
    allocate(tmplr(lr_, pmf%nface))
    ninternalface = 0
    do iface = 1, pmf%nface
      tmplr(lcell_, iface) = pmf%mgpart(pmf%facelr(lcell_, iface)) + 1
      if( iface .le. pmf%ninternalface ) then
        tmplr(rcell_, iface) = pmf%mgpart(pmf%facelr(rcell_, iface)) + 1
        if( tmplr(lcell_, iface) .ne. tmplr(rcell_, iface) ) then
          ninternalface = ninternalface + 1
        end if
      end if
    end do
    !> Copy unique faces in coarse grid numbering
    pmc%nnode         = pmf%nnode
    pmc%ninternalface = ninternalface
    pmc%ncell         = gr%npart
    nbface            = pmf%nface - pmf%ninternalface
    pmc%nface         = pmc%ninternalface + nbface
    pmc%npatch        = pmf%npatch
    !> Allocate the polyMesh
    call allocate_pm( pmc )
    !> Internal faces counters
    ninternalface = 0
    pmc%facelr    = 0
    !> Copy temp face data to coarse mesh
    do iface = 1, pmf%ninternalface
      if( tmplr(lcell_, iface) .ne. tmplr(rcell_, iface) ) then
        ninternalface = ninternalface + 1
        pmc%facelr(:, ninternalface)   = tmplr(:, iface)
        pmc%facenode(:, ninternalface) = pmf%facenode(:, iface)
      end if
    end do
    !> Boundary faces
    do iface = pmf%ninternalface + 1, pmf%nface
      ninternalface = ninternalface + 1
      pmc%facelr(:, ninternalface)   = tmplr(:, iface)
      pmc%facenode(:, ninternalface) = pmf%facenode(:, iface)
    end do
    !> Assign node pointer
    pmc%x         => pmf%x
    pmc%patchdata => pmf%patchdata
    !> Nullify parallel data
    nullify(pmc%cellgid)
    nullify(pmc%facegid)
    nullify(pmc%nodegid)
    deallocate(tmplr)
  end subroutine build_pm_coarse

  !> Simpler wrapper for metrics
  subroutine mesh_metrics( pm )
    use Wrap
    implicit none
    type( polyMesh ) :: pm
    !> Tapenade diff function
    call mesh_metrics_tapenade &
           ( pm%nnode, pm%nface, &
             pm%ninternalface, pm%ncell, &
             pm%x, pm%facelr, &
             pm%facenode, pm%cv, pm%cc, &
             pm%dn, pm%fs, pm%fc )
    !> Check metrics only if fine mesh
    if( pm%nlevel .eq. pm%ilevel ) then
      call check_metrics &
           ( pm%ncell, pm%nface, pm%cv, &
             pm%cc, pm%fc, pm%fs, pm%dn )
    end if
  end subroutine mesh_metrics

  !>
  subroutine allocate_pm( pm )
    implicit none
    type(polyMesh), intent(inout) :: pm
    !> Connectivity
    allocate( pm%facelr( lr_, pm%nface ) )
    allocate( pm%facenode( quadp1_, pm%nface ) )
    allocate( pm%edgenode( lr_, pm%nedge ) )
    allocate( pm%celltype( pm%ncell ) )
    allocate( pm%cellnode( max_cell_node_, pm%ncell ), &
              pm%cellface( max_cell_face_, pm%ncell ), &
              pm%celledge( max_cell_edge_, pm%ncell ) )
    !> Metrics
    allocate( pm%cc( dim_, pm%ncell ) )
    allocate( pm%cv( pm%ncell ) )
    allocate( pm%dn( dim_, pm%nface ) )
    allocate( pm%fc( dim_, pm%nface ) )
    allocate( pm%fs( pm%nface ) )
    !> Multigrid connectivity
    allocate( pm%mgpart( pm%ncell ) )
    !> Allocate nodes only to finest level
    !> and share it across levels
    if( pm%ilevel .eq. pm%nlevel ) then
      allocate( pm%x( dim_, pm%nnode ) )
      allocate( pm%patchdata( pproc_, pm%npatch) )
      !> Parallel data
      if( pm%parallel .eqv. .true. ) then
        allocate( pm%cellgid( pm%ncell ) )
        allocate( pm%facegid( pm%nface ) )
        allocate( pm%nodegid( pm%nnode ) )
      end if
    end if
  end subroutine allocate_pm

  !>
  subroutine mesh_metrics_tapenade( &
    nnode, nface, ninternalface, ncell, &
    x, facelr, facenode, cv, cc, dn, fs, fc )
    implicit none
    integer, intent(in)         :: nnode, nface, ninternalface, ncell
    real(kind=8), intent(in)    :: x(dim_, nnode)
    integer, intent(in)         :: facelr(lr_, nface), facenode(quadp1_, nface)
    real(kind=8), intent(out) :: cv(ncell), cc(dim_, ncell)
    real(kind=8), intent(out) :: dn(dim_, nface), fs(nface), fc(dim_, nface)
    !> Local variables
    integer :: iface, inode, i, il, ir, icell
    real(kind=8) :: r1(dim_), r2(dim_), r3(dim_), r4(dim_), fvol, fvolc(dim_)

    !> Zero out everything
    cv = 0.0d0; cc = 0.0d0
    !> Internal face loop for metric calculation
    do iface = 1, ninternalface
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
      !> Add CC/CV contribution of internal face
      cv(il)       = cv(il)   + fvol
      cc(:,il)     = cc(:,il) + fvolc
      cv(ir)       = cv(ir)   - fvol
      cc(:,ir)     = cc(:,ir) - fvolc
    end do
    !> Boundary face loop for metric calculation
    do iface = ninternalface + 1, nface
      r1 = x(:, facenode( 2, iface ))
      r2 = x(:, facenode( 3, iface ))
      r3 = x(:, facenode( 4, iface ))
      r4 = x(:, facenode( 5, iface ))
      il = facelr( lcell_, iface)
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
      !> Add CC/CV contribution of boundary face
      cv(il)       = cv(il)   + fvol
      cc(:,il)     = cc(:,il) + fvolc
    end do
    !> Do cell summation
    do icell = 1, ncell
      cc(:,icell)  = 0.50d0 * cc(:,icell) / cv(icell)
    end do
  end subroutine mesh_metrics_tapenade

  !> Construct the CSR graph from the edge2nodes data, which
  !> are input to the METIS partitioner routine
  subroutine pm_to_graph( pm, gr )
    use mpi
    implicit none
    type(polyMesh), intent(in)  :: pm
    type(crsGraph), intent(out) :: gr
    !> Local subroutine variables
    integer :: ierr, i
    !> Some simple memory checks
    if( allocated(gr%xadj)   .eqv. .true. ) deallocate(gr%xadj)
    if( allocated(gr%adjncy) .eqv. .true. ) deallocate(gr%adjncy)
    if( allocated(gr%adjwgt) .eqv. .true. ) deallocate(gr%adjwgt)
    if( allocated(gr%vsurf)  .eqv. .true. ) deallocate(gr%vsurf)
    !> Allocations
    gr%nvtx = pm%ncell
    allocate( gr%xadj(gr%nvtx + 1), &
              gr%adjncy(lr_ * pm%ninternalface), &
              gr%adjwgt(lr_ * pm%ninternalface), &
              gr%vsurf(gr%nvtx), gr%vvol(gr%nvtx) )
    gr%xadj   = 0
    gr%adjncy = -1
    gr%vsurf  = 0.0d0
    gr%vvol   = pm%cv
    gr%adjwgt = -100.0d0
    !> First pass form the sizes
    do i = 1, pm%ninternalface
      gr%xadj(pm%facelr(lcell_, i)) = gr%xadj(pm%facelr(lcell_, i)) + 1
      gr%xadj(pm%facelr(rcell_, i)) = gr%xadj(pm%facelr(rcell_, i)) + 1
    end do
    gr%xadj = cshift(gr%xadj, gr%nvtx)
    gr%xadj(1) = 0
    !write(*,*) gr%xadj
    !> All boundary face surface area (attach to cell)
    do i = pm%ninternalface + 1, pm%nface
      gr%vsurf(pm%facelr(lcell_, i)) = gr%vsurf(pm%facelr(lcell_, i)) + pm%fs(i)
    end do
    !> Form the xadj offsets (Can replace with intrisic?)
    do i = 1, gr%nvtx
      gr%xadj(i+1) = gr%xadj(i) + gr%xadj(i+1)
    end do
    !> Check size match
    if( gr%xadj(gr%nvtx + 1) .ne. (lr_ * pm%ninternalface) ) then
      write(*,*) 'Error: Sizes of facelr xadj do not match ninternalface', &
        gr%xadj(gr%nvtx + 1), 2 * pm%ninternalface
    end if
    !> Second pass form the adjncy
    do i = 1, pm%ninternalface
      !> Push the left face data
      gr%xadj(pm%facelr(lcell_, i))            = gr%xadj(pm%facelr(lcell_, i)) + 1
      gr%adjncy(gr%xadj(pm%facelr(lcell_, i))) = pm%facelr(rcell_, i) - 1
      gr%adjwgt(gr%xadj(pm%facelr(lcell_, i))) = pm%fs(i)
      !> Push the right face data
      gr%xadj(pm%facelr(rcell_, i))            = gr%xadj(pm%facelr(rcell_, i)) + 1
      gr%adjncy(gr%xadj(pm%facelr(rcell_, i))) = pm%facelr(lcell_, i) - 1
      gr%adjwgt(gr%xadj(pm%facelr(rcell_, i))) = pm%fs(i)
    end do
    gr%xadj = cshift(gr%xadj, gr%nvtx)
    gr%xadj(1) = 0
    !write(*,*) gr%xadj
    !> Check to see if the adjncy array is formed correctly
    if(count(gr%adjncy .lt. 0) .ne. 0) then
      write(*,*) 'Error: Adjncy array from edge2nodes not constructed correctly'
    end if
    !> Check to see if the adjncywgt array is formed correctly
    if(count(gr%adjwgt .lt. 0.0d0) .ne. 0) then
      write(*,*) 'Error: Adjncy weight array from edge2nodes not constructed correctly'
    end if

  end subroutine pm_to_graph

end module Mesh

