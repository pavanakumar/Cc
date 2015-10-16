module mesh
  use iso_c_binding
  implicit none

  type polyMesh
    integer :: nnode, nface, ncell, npatch
    real(kind=8), dimension(:,:) :: node
    integer, dimension(:,:), allocatable :: facelr
    integer, dimension(:,:), allocatable :: facenode
    !!! Comprises of patchstart, patchsize, patchtype, patchneigh
    integer, dimension(:,:), allocatable :: patchdata
    !!! Optional cell data
    integer, dimension(:,:) :: cellnode !! Ordered
    integer, dimension(:,:) :: cellface !! Ordered
  end type polyMesh

  !!! Read the polyMesh data from OpenFOAM reader
  !!! uses parmgen and creates the multigrid levels
  subroutine reader_of( nlevel, pm )
    implicit none
    integer, intent(in)           :: nlevel
    type(polyMesh), intent(inout) :: pm(nlevel)



  end subroutine reader_of

end module mesh


