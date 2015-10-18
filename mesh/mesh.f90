module Mesh
  implicit none
!!!!!!!!!!!!  CONSTANTS  !!!!!!!!!!!!!
  !!! Comprises of patchstart, patchsize, patchtype, patchneigh
  integer, parameter :: pstart_ = 1, psize_ = 2, ptype_ = 3, pproc_ = 4
  integer, parameter :: dim_ = 3, lr_ = 2
  integer, parameter :: tri_ = 3, quad_ = 4, quadp1_ = 5
  integer, parameter :: enable_parallel_ = 0
  integer, parameter :: processor_bc_ = 0, wall_bc_ = 1, symmetry_bc_ = 2, &
                        inlet_bc_ = 3, outlet_bc_ = 4, riemann_bc_ = 5
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
    real(kind=8), dimension(:,:), allocatable :: x, cc, dn, fc
    real(kind=8), dimension(:), allocatable   :: cv, fs
    !!! Parallel run (global numbering of local entities)
    logical                            :: parallel = .false.
    integer, dimension(:), allocatable :: cellgid, nodegid, facegid
    !!! Multigrid variables

  end type polyMesh

  type crsGraph
    integer, dimension(:), allocatable :: xadj, adjncy, part 
  end type crsGraph

  contains

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
      if( ipar .eq. enable_parallel_ ) then
        pm(ilvl)%parallel = .true.
      end if
      pm(ilvl)%nlevel = nlevel
      pm(ilvl)%ilevel = ilvl
    end do
    !!! Read in the finest mesh from OpenFOAM    
    call init_of_mesh( ipar )
    call get_pm_sizes( pm(nlevel)%nnode, pm(nlevel)%nface, &
                       pm(nlevel)%ninternalface, &
                       pm(nlevel)%ncell, pm(nlevel)%npatch )
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
         ( pm%nnode, pm%nface, &
           pm%ninternalface, pm%ncell, &
           pm%x, pm%facelr, &
           pm%facenode, pm%cv, pm%cc, &
           pm%dn, pm%fs, pm%fc )
  end subroutine mesh_metrics 

  subroutine allocate_pm( pm )
    implicit none
    type(polyMesh), intent(inout) :: pm
    !!! Connectivity
    allocate( pm%facelr( lr_, pm%nface ), &
              pm%facenode( quadp1_, pm%nface ), &
              pm%patchdata( pproc_, pm%npatch) )
    !!! Metrics
    allocate( pm%x( dim_, pm%nnode ), &
              pm%cc( dim_, pm%ncell ), &
              pm%cv( pm%ncell ), &
              pm%dn( dim_, pm%nface ), &
              pm%fc( dim_, pm%nface ), &
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
    real(kind=8), intent(in) :: x(dim_, nnode)
    integer, intent(in)      :: facelr(lr_, nface), facenode(quadp1_, nface)
    real(kind=8), intent(in) :: cv(ncell), cc(dim_, ncell)
    real(kind=8), intent(in) :: dn(dim_, nface), fs(nface), fc(dim_, nface)

  end subroutine mesh_metrics_tapenade

end module Mesh

