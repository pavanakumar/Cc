!> The polyhedral mesh defined constants
module PolyMeshConstant
  use iso_c_binding
  implicit none
  !> Constant 1/12 hard coded for quadrature
  real(kind=8), parameter :: oneby12_ = 0.08333333333333333333333333333333333333333d0
  !> Quadrature points (for cell-centroid calculation)
  real(kind=8), parameter :: xi1_     = 0.78867513459481288225457439025097872782380d0
  real(kind=8), parameter :: xi2_     = 0.21132486540518711774542560974902127217620d0
  real(kind=8), parameter :: eta1_    = 0.78867513459481288225457439025097872782380d0
  real(kind=8), parameter :: eta2_    = 0.21132486540518711774542560974902127217620d0
  !>
  integer, parameter :: pstart_        = 1, &
                        psize_         = 2, &
                        ptype_         = 3, &
                        pproc_         = 4
  !>
  integer, parameter :: dim_           = 3, &
                        lr_            = 2, &
                        lcell_         = 1, &
                        rcell_         = 2
  !>
  integer, parameter :: tri_           = 3, &
                        quad_          = 4, &
                        quadp1_        = 5
  !>
  integer, parameter :: processor_bc_  = 0, &
                        wall_bc_       = 1, &
                        symmetry_bc_   = 2, &
                        inlet_bc_      = 3, &
                        outlet_bc_     = 4, &
                        riemann_bc_    = 5, &
                        empty_bc_      = 6
  !>
  integer, parameter :: max_cell_node_ = 8, &
                        max_cell_face_ = 6, &
                        max_cell_edge_ = 12
  !>
  integer, parameter :: enable_parallel_ = 1
end module PolyMeshConstant

