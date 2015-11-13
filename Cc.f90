program Cc
  use Wrap
  use Mesh
  use mpi
  implicit none
  integer, parameter :: nlevel = 3
  integer :: ipar = enable_parallel_ + 1, ierr, rank = 0, ilvl
  type(polyMesh), dimension(nlevel) :: pm

  !! Read the mesh and create coarse levels
  call create_mg_pm( nlevel, pm, ipar )
  !! Get the MPI rank if running in parallel
  if( ipar .eq. enable_parallel_ ) then
    call MPI_Comm_rank( MPI_COMM_WORLD, rank, ierr )
  end if
!  call tecio_write( rank, nlevel, pm(1:nlevel) )
  !! Write the MG levels for each rank
  do ilvl = nlevel, 1, -1
    call write_pm_tecio( ilvl, rank, pm(ilvl)%nnode, pm(ilvl)%ncell,&
                         pm(ilvl)%nface, pm(ilvl)%ninternalface,&
                         pm(ilvl)%x, pm(ilvl)%facelr, pm(ilvl)%facenode )
  end do 

end program Cc

