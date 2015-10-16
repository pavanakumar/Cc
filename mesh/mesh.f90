module Mesh
  implicit none

  !!! Comprises of patchstart, patchsize, patchtype, patchneigh
  integer, parameter :: _pstart = 1, _psize = 2, _ptype = 3, _pproc = 4
  integer, parameter :: _dim = 3, _lr = 2
  integer, parameter :: _tri = 3, _quad = 4
 
  type polyMesh
    integer :: nnode, nface, ninternalface,
               ncell, npatch, ilevel, nlevel
    integer, dimension(:,:), allocatable :: facelr
    integer, dimension(:,:), allocatable :: facenode
    integer, dimension(:,:), allocatable :: patchdata
    !!! Metrics
    !!! xyz coords, cell-center, unit-normal, face-center
    real(kind=8), dimension(:,:) :: x, cc, dn, fc
    !!! cell-volume, face-area
    real(kind=8), dimension(:) :: cv, fs
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
    call mesh_metrics_tapenade&
         ( pm(nlevel)%nnode, pm(nlevel)%nface,&
           pm(nlevel)%ninternalface,&
           pm(nlevel)%x, pm(nlevel)%facelr,&
           pm(nlevel)%facenode, pm(nlevel)%cv,&
           pm(nlevel)%dn, pm(nlevel)%fs, pm(nlevel)%fc )
  end subroutine mesh_metrics 

  subroutine allocate_pm( pm )
    implicit none
    type(polyMesh), intent(inout) :: pm
    !!! Connectivity
    allocate( pm%facelr( _lr, pm%nface ),&
              pm%facenode( _quad + 1, pm%nface ), &
              pm%patchdata( _pproc, pm%npatch) )
    !!! Metrics
    allocate( pm%x( _dim, pm%nnode ), &
              pm%cc( _dim, pm%ncell ), &
              pm%cv( pm%ncell ), &
              pm%dn( _dim, pm%nface ), &
              pm%fc( _dim, pm%nface ),&
              pm%fs( pm%nface ) )
  end subroutine allocate_pm

end module Mesh

