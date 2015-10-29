program Cc
  use Wrap
  use Mesh
  use mpi
  implicit none
  integer, parameter :: nlevel = 2
  integer :: ipar = enable_parallel_ + 1
  type(polyMesh), dimension(nlevel) :: pm
  
  call create_mg_pm( nlevel, pm, ipar )
  
  call write_cc_tecio( 1, 0, pm(nlevel)%nnode, pm(nlevel)%ncell,&
                       pm(nlevel)%nface, pm(nlevel)%ninternalface,&
                       pm(nlevel)%x, pm(nlevel)%facelr, pm(nlevel)%facenode )

end program Cc

