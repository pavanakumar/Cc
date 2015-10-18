module Wrap
  use iso_c_binding
  implicit none

  interface
    
    subroutine init_of_mesh( ipar ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: ipar
    end subroutine init_of_mesh

    subroutine get_pm_sizes( nnode, nface, ninternalface, ncell, npatch ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: nnode, nface, ncell, npatch
    end subroutine get_pm_sizes

    subroutine get_pm_nodes( nnode, x ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: nnode
      real(kind=c_double), dimension(*) :: x
    end subroutine get_pm_nodes

    subroutine get_pm_faces( nface, ninternalface, facelr, facenode ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: nface, ninternalface
      integer(kind=c_int), dimension(*) :: facelr, facenode
    end subroutine get_pm_faces

    subroutine get_pm_patches( npatch, patchdata ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: npatch
      integer(kind=c_int), dimension(*) :: patchdata
    end subroutine get_pm_faces

    subroutine close_of_mesh() bind(C)
      use iso_c_binding
      implicit none
    end subroutine close_of_mesh

  end interface  

end module Wrap

