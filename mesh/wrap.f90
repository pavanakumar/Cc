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
      integer(kind=c_int) :: nnode, nface, ninternalface, ncell, npatch
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
    end subroutine get_pm_patches

    subroutine check_metrics( ncell, nface, cv, cc, fc, fs, dn ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: ncell, nface
      real(kind=c_double), dimension(*) :: cv, cc, fc, fs, dn
    end subroutine check_metrics

    subroutine close_of_mesh() bind(C)
      use iso_c_binding
      implicit none
    end subroutine close_of_mesh

    subroutine MGridGen_f90( nvtxs, xadj, vvol, vsurf,&
                             adjncy, adjwgt, minsize, maxsize,&
                             options, nmoves, nparts, part) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: nvtxs, minsize, maxsize, nmoves, nparts
      integer(kind=c_int), dimension(*) :: xadj, adjncy, part, options
      real(kind=c_double), dimension(*) :: vvol, vsurf, adjwgt
    end subroutine MGridGen_f90

    subroutine write_pm_tecio( ilevel, irank, nnode, ncell,&
                               nface, ninternalface, xyz,&
                               facelr, facenode) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: ilevel, irank, nnode, ncell, nface, ninternalface
      real(kind=c_double), dimension(*) :: xyz
      integer(kind=c_int), dimension(*) :: facelr, facenode

    end subroutine write_pm_tecio

  end interface  

end module Wrap

