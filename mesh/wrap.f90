module Wrap
  use iso_c_binding
  implicit none

  interface
    
    subroutine init_mesh_api( ipar ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: ipar
    end subroutine init_mesh_api

    subroutine get_pm_sizes( nnode, nface, ninternalface, &
                             nedge, ninternaledge, &
                             ncell, npatch ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: nnode, nface, ninternalface, &
                             nedge, ninternaledge, &
                             ncell, npatch
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

    subroutine get_pm_edges( nedge, edgenodes ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: nedge
      integer(kind=c_int), dimension(*) :: edgenodes
    end subroutine get_pm_edges

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

    subroutine close_mesh_api() bind(C)
      use iso_c_binding
      implicit none
    end subroutine close_mesh_api

    subroutine finalize_mesh_api() bind(C)
      use iso_c_binding
      implicit none
    end subroutine finalize_mesh_api

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

    subroutine sfc_perm( n, x, xmin, xmax, perm, iperm ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: n, perm(n), iperm(n)
      real(kind=c_double) :: x(3, n), xmin(3), xmax(3)
    end subroutine sfc_perm

    subroutine get_cellgid( ncell, cellgid ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: ncell, cellgid(*)
    end subroutine get_cellgid

    subroutine get_nodegid( nnode, nodegid ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: nnode, nodegid(*)
    end subroutine get_nodegid

    subroutine get_facegid( nface, facegid ) bind(C)
      use iso_c_binding
      implicit none
      integer(kind=c_int) :: nface, facegid(*)
    end subroutine get_facegid

  end interface  

end module Wrap

