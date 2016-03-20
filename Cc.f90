module Cc
  use iso_c_binding
  use PolyMeshMG
  implicit none
  type(polyMesh), allocatable :: pm(:)

  contains

  !>
  subroutine init_cc( parflag, nlevel ) bind (C)
    implicit none
    integer, intent(in) :: parflag, nlevel
    !> Allocate the polymesh data-structure and read from API
    allocate( pm(nlevel) )
    call create_mg_pm( nlevel, pm, parflag )

  end subroutine init_cc

  !> 
  subroutine finalise_cc() bind(C)
    implicit none
    deallocate( pm )
  end subroutine finalise_cc

  !>
  subroutine write_cc_mesh() bind(C)
    use Wrap
    implicit none
    integer :: ilvl, nlevel, rank
    nlevel = size(pm)
    rank = 0
    do ilvl = nlevel, 1, -1
      call write_pm_tecio( ilvl, rank, pm(ilvl)%nnode, pm(ilvl)%ncell,&
                           pm(ilvl)%nface, pm(ilvl)%ninternalface,&
                           pm(ilvl)%x, pm(ilvl)%facelr, pm(ilvl)%facenode )
    end do 
  end subroutine write_cc_mesh

end module Cc

