module PolyMeshMG
  use iso_c_binding
  use PolyMeshMath
  use Graph
  implicit none

  !>
  type :: polyMesh
    logical :: parallel      = .false.
    !> Sizes (default values)
    integer :: nnode         = 0, &
               nface         = 0, &
               ninternalface = 0, &
               nedge         = 0, &
               ninternaledge = 0, &
               ncell         = 0, &
               npatch        = 0, &
               ilevel        = 0, &
               nlevel        = 0, &
               nmgpart       = 0
    !> Sizes MPI node based (Node based)
    integer :: mpi_nod_nxadj   = 0, &
               mpi_nod_nadjncy = 0, &
               mpi_nod_nproc   = 0, &
               mpi_nod_nsnlist = 0
    !> Connecivity/Topology
    integer, allocatable :: facelr(:,:),    &
                            cellface(:,:),  &
                            nfacenode(:),   &
                            facenode(:,:),  &
                            edgenode(:,:),  &
                            cellnode(:,:),  &
                            celledge(:,:),  &
                            celltype(:),    &
                            mgpart(:)
    !> Pointer values shared between MG levels
    integer, pointer     :: patchdata(:,:), &
                            cellgid(:),     &
                            nodegid(:),     &
                            facegid(:)
    !> Grid coordinates shared across MG levels
    real(kind=8), pointer     :: x(:,:)
    !> Metrics
    real(kind=8), allocatable :: cc(:,:), & ! Cell centroid
                                 cv(:),   & ! Cell volume
                                 dn(:,:), & ! Face unit normal
                                 fc(:,:), & ! Face centroid
                                 fs(:)      ! Face area
    !> Node based MPI schedule
    integer, allocatable :: mpi_nod_xadj(:),   &
                            mpi_nod_adjncy(:), &
                            mpi_nod_proc(:),   &
                            mpi_nod_snlist(:)
  end type polyMesh

  contains

  !> Create the multi-grid levels using MGridGen
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
                       pm(nlevel)%nfacenode, pm(nlevel)%facenode, &
                       pm(nlevel)%facelr )
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
    call get_pm_extra( pm(nlevel)%ncell, pm(nlevel)%celltype, &
                       pm(nlevel)%cellnode, pm(nlevel)%cellface, &
                       pm(nlevel)%celledge )
    call mesh_metrics( pm(nlevel) ) 
    !> Close the interface
    call free_mesh_api()
  end subroutine reader_of

  !> Create the Multi-Grid levels
  !> using mgridgen wrapper
  subroutine create_mg_levels( nlevel, pm, ipar )
    use Wrap
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
    !> Local variables
    integer :: ilvl
    type(crsGraph) :: gr
    real(kind=8) :: sum_fine, sum_coarse
 
    !> Multi-grid mgridgen
    do ilvl = nlevel, 2, -1
      call pm_to_graph( pm(ilvl), gr )
      call graph_mg( gr, pm(ilvl)%nmgpart, pm(ilvl)%mgpart )
      call build_pm_coarse( gr, pm(ilvl), pm(ilvl - 1) )
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
  subroutine build_pm_coarse( gr, pmf, pmc )
    implicit none
    type(crsGraph)       :: gr
    type(polyMesh)       :: pmf, pmc
    !> Local variables
    integer, allocatable :: tmplr(:,:)
    integer              :: ninternalface, nbface, iface

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
    pmc%ncell         = pmf%nmgpart
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
        ninternalface                   = ninternalface + 1
        pmc%facelr(:, ninternalface)    = tmplr(:, iface)
        pmc%nfacenode(ninternalface)    = pmf%nfacenode(iface)
        pmc%facenode(:, ninternalface)  = pmf%facenode(:, iface)
      end if
    end do
    !> Boundary faces
    do iface = pmf%ninternalface + 1, pmf%nface
      ninternalface                  = ninternalface + 1
      pmc%facelr(:, ninternalface)   = tmplr(:, iface)
      pmc%nfacenode(ninternalface)   = pmf%nfacenode(iface)
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
    call poly_metrics &
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
    allocate( pm%nfacenode( pm%nface ) )
    allocate( pm%facenode( quad_, pm%nface ) )
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

  !> Construct the CSR graph from the edge2nodes data, which
  !> are input to the MgridGen multi-grid routine
  subroutine pm_to_graph( pm, gr )
    implicit none
    type(polyMesh), intent(in)  :: pm
    type(crsGraph), intent(out) :: gr
    !> Local variables
    integer :: nbc
    nbc = pm%nface - pm%ninternalface
    !> Allocate memory for the graph
    call allocate_graph( pm%ncell, pm%ninternalface, gr )
    call edge_to_graph( pm%ncell, pm%ninternalface, pm%facelr, gr )
    call edge_data_to_graph( pm%ncell, pm%ninternalface, nbc, &
                             pm%facelr, pm%cv, pm%fs, gr )
  end subroutine pm_to_graph

end module PolyMeshMG

