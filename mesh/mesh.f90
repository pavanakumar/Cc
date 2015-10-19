module Constants
  real(kind=8), parameter :: oneby3_ = 0.33333333333333333333333333333333333333d0
  real(kind=8), parameter :: oneby6_ = 0.16666666666666666666666666666666666666666d0
  real(kind=8), parameter :: oneby12_ = 0.0833333333333333333333333333333333333333d0
end module Constants

module Mesh
  use Constants
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

  function cross_prod(x, y)
    implicit none
    real(kind=8),dimension(dim_),intent(in) :: x, y
    real(kind=8),dimension(dim_)            :: cross_prod
    integer, dimension(3), parameter        :: c1 = (/2,3,1/), &
                                               c2 = (/3,1,2/)
    cross_prod = x(c1) * y(c2) - y(c1) * x(c2)
  end function cross_prod

  subroutine create_mg_pm( nlevel, pm, ipar )
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
 
    call reader_of( nlevel, pm, ipar )
    call create_mg_levels( nlevel, pm, ipar )

  end subroutine create_mg_pm

  !!! Read the polyMesh data from OpenFOAM reader
  !!! uses parmgen and creates the multigrid levels
  subroutine reader_of( nlevel, pm, ipar )
    use Wrap
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
    !!! Local variables
    integer :: ilvl
 
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
    call close_of_mesh()
  end subroutine reader_of

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! Create the Multi-Grid levels  !!
  !!!! using mgridgen wrapper        !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine create_mg_levels( nlevel, pm, ipar )
    use Wrap
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
    !!! Local variables
    integer :: ilvl, ilvl_c
    type(crsGraph) :: gr
 
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

  end subroutine create_mg_levels

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
    integer, intent(in)         :: nnode, nface, ninternalface, ncell
    real(kind=8), intent(in)    :: x(dim_, nnode)
    integer, intent(in)         :: facelr(lr_, nface), facenode(quadp1_, nface)
    real(kind=8), intent(inout) :: cv(ncell), cc(dim_, ncell)
    real(kind=8), intent(inout) :: dn(dim_, nface), fs(nface), fc(dim_, nface)
    !!! Local variables
    integer :: iface, inode, i, il, ir
    real(kind=8) :: r1(dim_), r2(dim_), r3(dim_), r4(dim_), fvol

    !!! Zero out everything    
    cv = 0.0d0; cc = 0.0d0
    !!! Calculate the face centroids, area, and unit normal vector
    do iface = 1, nface
      r1 = x(:, facenode( 2, iface ))
      r2 = x(:, facenode( 3, iface ))
      r3 = x(:, facenode( 4, iface ))
!!! Planar face
      if( facenode(1, iface) .eq. tri_ ) then
        !!! Face centroid
        fc(:, iface) = oneby3_ * ( r1 + r2 + r3 )
        !!! Face normal and area
        dn(:, iface) = 0.50d0 * cross_prod( (r2 - r1), (r3 - r1) )
        fs(iface)    = sqrt( sum( dn(:, iface) * dn(:, iface) ) )
        dn(:, iface) = dn(:, iface) / fs(iface)
        !!! Face volume contribution
        fvol   = oneby6_ * sum( r1 * cross_prod(r2, r3) )
        il     = facelr(1, iface)
        ir     = facelr(2, iface)
        cv(il) = cv(il) + fvol
        if( iface .le. ninternalface ) cv(ir) = cv(ir) - fvol
        !!! Face cell centroid contribution

!!! Ruled face
      else if( facenode(1, iface) .eq. quad_ ) then
        r4 = x(:, facenode( 5, iface ))
        !!! Face centroid
        fc(:, iface) = 0.250d0 * ( r1 + r2 + r3 + r4 )
        !!! Face normal and area
        dn(:, iface) = 0.50d0 * cross_prod( (r3 - r1), (r4 - r2) )
        fs(iface)    = sqrt( sum( dn(:, iface) * dn(:, iface) ) )
        dn(:, iface) = dn(:, iface) / fs(iface)
        !!! Face volume contribution
        fvol   = oneby12_ * sum( (r1 + r2) * cross_prod((r1 + r2), (r3 + r4)) )
        il     = facelr(1, iface)
        ir     = facelr(2, iface)
        cv(il) = cv(il) + fvol
        if( iface .le. ninternalface ) cv(ir) = cv(ir) - fvol
        !!! Face cell centroid contribution
        
!!! Error unknow face type
      else

      end if
    end do

  end subroutine mesh_metrics_tapenade

end module Mesh

