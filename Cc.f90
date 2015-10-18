program Cc
  use Mesh
  use mpi
  implicit none
  integer, parameter :: nlevel = 2
  integer :: ipar = enable_parallel_ + 1
  type(polyMesh), dimension(nlevel) :: pm
  
  call create_mg_pm( nlevel, pm, ipar )


end program Cc
