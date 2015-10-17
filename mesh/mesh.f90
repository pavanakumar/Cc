module Mesh
  implicit none
!!!!!!!!!!!!  CONSTANTS  !!!!!!!!!!!!!
  !!! Comprises of patchstart, patchsize, patchtype, patchneigh
  integer, parameter :: _pstart = 1, _psize = 2, _ptype = 3, _pproc = 4
  integer, parameter :: _dim = 3, _lr = 2
  integer, parameter :: _tri = 3, _quad = 4, _quadp1 = 5
  integer, parameter :: _enable_parallel = 0
  integer, parameter :: _processor_bc = 0, _wall_bc = 1, _symmetry_bc = 2, &
                        _inflow_bc = 3, _outflow_bc = 4, _riemann_bc = 5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type polyMesh
    !!! Sizes
    integer :: nnode = 0, nface = 0, ninternalface = 0, &
               ncell = 0, npatch = 0, ilevel = 0, nlevel = 0
    !!! Connecivity/Topology
    integer, dimension(:,:), allocatable :: facelr
    integer, dimension(:,:), allocatable :: facenode
    integer, dimension(:,:), allocatable :: patchdata
    !!! Metrics
    !!! xyz coords, cell-center, unit-normal, face-center
    !!! cell-volume, face-area
    real(kind=8), dimension(:,:) :: x, cc, dn, fc
    real(kind=8), dimension(:)   :: cv, fs
    !!! Parallel run (global numbering of local entities)
    logical                            :: parallel = .false.
    integer, dimension(:), allocatable :: cellgid, nodegid, facegid
    !!! Multigrid variables

  end type polyMesh

  type crsGraph
    integer, dimension(:), allocatable :: xadj, adjncy, part 
  end type crsGraph

  !!! Read the polyMesh data from OpenFOAM reader
  !!! uses parmgen and creates the multigrid levels
  subroutine reader_of( nlevel, pm, ipar )
    use Wrap
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
    !!! Local variables
    integer :: ilvl, ilvl_c
    type(crsGraph) :: gr
 
    !!! Make the meshes aware of the multi-grid 
    do ilvl = 1, nlevel
      if( ipar .eq. _enable_parallel ) then
        pm(ilvl)%parallel = .true.
      end if
      pm(ilvl)%nlevel = nlevel
      pm(ilvl)%ilevel = ilvl
    end do
    !!! Read in the finest mesh from OpenFOAM    
    call init_of_mesh( ipar )
    call get_pm_sizes( pm(nlevel)%nnode, pm(nlevel)%nface, &
                       pm%(nlevel)%ninternalface, &
                       pm%(nlevel)%ncell, pm(nlevel)%npatch )
    call allocate_pm( pm(nlevel) )
    call get_pm_nodes( pm(nlevel)%nnode, pm(nlevel)%x )
    call get_pm_faces( pm(nlevel)%nface, pm(nlevel)%ninternalface, &
                       pm(nlevel)%facelr, pm(nlevel)%facenode )
    call get_pm_patches( pm(nlevel)%npatch, pm(nlevel)%patchdata )
    !!! Calculate the metrics
    call mesh_metrics( pm(nlevel) )
    !!! Multi-grid mgridgen
    do ilvl = nlevel, 2, -1
      ilvl_c = ilvl - 1
!      call allocate_graph( pm(ilvl), gr )
!      call form_graph( pm(ilvl), gr )
!      call mgridgen( gr%xadj, gr%adjncy, gr%part, pm(ilvl)%cv, &
!                      pm(ilvl)%fs, pm(ilvl)%patchdata )
!      call fine_to_coarse( pm(ilvl), pm(ilvl_c), gr%part )
!      call deallocate_graph( gr )
!      call mesh_metrics( pm(ilvl_c) ) 
    end do
    call close_of_mesh()

  end subroutine reader_of

  !!! Simpler wrapper for metrics
  subroutine mesh_metrics( pm )
    implicit none
    type( polyMesh ) :: pm
    !!! Tapenade diff function
    call mesh_metrics_tapenade &
         ( pm(nlevel)%nnode, pm(nlevel)%nface, &
           pm(nlevel)%ninternalface, pm(nlevel)%ncell, &
           pm(nlevel)%x, pm(nlevel)%facelr, &
           pm(nlevel)%facenode, pm(nlevel)%cv, pm(nlevel)%cc, &
           pm(nlevel)%dn, pm(nlevel)%fs, pm(nlevel)%fc )
  end subroutine mesh_metrics 

  subroutine allocate_pm( pm )
    implicit none
    type(polyMesh), intent(inout) :: pm
    !!! Connectivity
    allocate( pm%facelr( _lr, pm%nface ), &
              pm%facenode( _quadp1, pm%nface ), &
              pm%patchdata( _pproc, pm%npatch) )
    !!! Metrics
    allocate( pm%x( _dim, pm%nnode ), &
              pm%cc( _dim, pm%ncell ), &
              pm%cv( pm%ncell ), &
              pm%dn( _dim, pm%nface ), &
              pm%fc( _dim, pm%nface ), &
              pm%fs( pm%nface ) )
    !!! Parallel data
    if( pm%parallel .eqv. .true. ) then
      allocate( pm%cellgid( pm%ncell ), &
                pm%facegid( pm%nface ), &
                pm%nodegid( pm%nnode ) )
    end if
  end subroutine allocate_pm

  subroutine mesh_metrics_tapenade( &
    nnode, nface, ninternalface, ncell, &
    x, facelr, facenode, &
    cv, cc, dn, fs, fc )
    implicit none
    integer, intent(in)      :: nnode, nface, ninternalface, ncell
    real(kind=8), intent(in) :: x(_dim, nnode)
    integer, intent(in)      :: facelr(_lr, nface), facenode(_quadp1, nface)
    real(kind=8), intent(in) :: cv(ncell), cc(_dim, ncell)
    real(kind=8), intent(in) :: dn(_dim, nface), fs(nface), fc(_dim, nface)

  end subroutine mesh_metrics_tapenade

end module Mesh

