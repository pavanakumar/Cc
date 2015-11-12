!!!! Mesh module
module Constants
  real(kind=8), parameter :: oneby3_  = 0.33333333333333333333333333333333333333333d0
  real(kind=8), parameter :: oneby6_  = 0.16666666666666666666666666666666666666666d0
  real(kind=8), parameter :: oneby12_ = 0.08333333333333333333333333333333333333333d0
  !!!! Quadrature points (for cell-centroid)
  real(kind=8), parameter :: xi1_     = 0.78867513459481288225457439025097872782380d0
  real(kind=8), parameter :: xi2_     = 0.21132486540518711774542560974902127217620d0
  real(kind=8), parameter :: eta1_    = 0.78867513459481288225457439025097872782380d0
  real(kind=8), parameter :: eta2_    = 0.21132486540518711774542560974902127217620d0
end module Constants

module Mesh
  use iso_c_binding
  use Constants
  implicit none
!!!!!!!!!!!!  CONSTANTS  !!!!!!!!!!!!!
  !!! Comprises of patchstart, patchsize, patchtype, patchneigh
  integer, parameter :: pstart_ = 1, psize_ = 2, ptype_ = 3, pproc_ = 4
  integer, parameter :: dim_ = 3, lr_ = 2, lcell_ = 1, rcell_ = 2
  integer, parameter :: tri_ = 3, quad_ = 4, quadp1_ = 5
  integer, parameter :: enable_parallel_ = 0
  integer, parameter :: processor_bc_ = 0, wall_bc_ = 1, symmetry_bc_ = 2, &
                        inlet_bc_ = 3, outlet_bc_ = 4, riemann_bc_ = 5, &
                        empty_bc_ = 6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type :: polyMesh
    !!! Sizes
    integer :: nnode = 0, nface = 0, ninternalface = 0, &
               ncell = 0, npatch = 0, ilevel = 0, nlevel = 0
    !!! Connecivity/Topology
    integer, dimension(:,:), allocatable :: facelr
    integer, dimension(:,:), allocatable :: facenode
    integer, dimension(:,:), pointer     :: patchdata
    !!! Metrics
    !!! xyz coords, cell-center, unit-normal, face-center
    !!! cell-volume, face-area
    real(kind=8), dimension(:,:), pointer     :: x
    real(kind=8), dimension(:,:), allocatable :: cc, dn, fc
    real(kind=8), dimension(:), allocatable   :: cv, fs
    !!! MPI parallel run (global numbering of local entities)
    logical                            :: parallel = .false.
    integer, dimension(:), pointer     :: cellgid, nodegid, facegid
    !!! OMP parallel run data
    integer                            :: ncolour
    integer, dimension(:), allocatable :: colourxadj
    !!! Multigrid variables
    integer, dimension(:), allocatable :: mgpart
  end type polyMesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type :: crsGraph
    integer                                 :: nvtx = 0, npart = 0
    integer, dimension(:), allocatable      :: xadj, adjncy
    real(kind=8), dimension(:), allocatable :: adjwgt, vvol, vsurf
  end type crsGraph

  contains

  function cross_prod(x, y)
    implicit none
    real(kind=8),dimension(dim_),intent(in) :: x, y
    real(kind=8),dimension(dim_)            :: cross_prod
    integer, dimension(dim_), parameter        :: c1 = (/2,3,1/), &
                                               c2 = (/3,1,2/)
    cross_prod = x(c1) * y(c2) - y(c1) * x(c2)

  end function cross_prod

  function r_1234(r1, r2, r3, r4, xi, eta)
    implicit none
    real(kind=8),dimension(dim_),intent(in) :: r1, r2, r3, r4
    real(kind=8), intent(in)                :: eta, xi
    real(kind=8),dimension(dim_)            :: r_1234
    r_1234 = ( r1 * ( 1.0d0 - xi ) + r2 * xi ) * ( 1.0d0 - eta ) + &
             ( r4 * ( 1.0d0 - xi ) + r3 * xi ) * eta

  end function r_1234

  function dS_1234(r1, r2, r3, r4, xi, eta)
    implicit none
    real(kind=8),dimension(dim_),intent(in) :: r1, r2, r3, r4
    real(kind=8), intent(in)                :: eta, xi
    real(kind=8),dimension(dim_)            :: dS_1234, temp1, temp2
    temp1   = (r2 - r1) * ( 1.0d0 - eta ) + (r3 - r4) * eta
    temp2   = (r4 - r1) * ( 1.0d0 - xi  ) + (r3 - r2) * xi
    dS_1234 = cross_prod( temp1, temp2 )

  end function dS_1234

  function func_1234(r1, r2, r3, r4, xi, eta)
    implicit none
    real(kind=8),dimension(dim_),intent(in) :: r1, r2, r3, r4
    real(kind=8), intent(in)                :: eta, xi
    real(kind=8),dimension(dim_)            :: temp, func_1234
    real(kind=8)                            :: temp_dot
    temp      = r_1234(r1, r2, r3, r4, xi, eta)
    temp_dot  = sum( temp * temp )
    func_1234 = temp_dot * dS_1234(r1, r2, r3, r4, xi, eta)

  end function func_1234

  subroutine create_mg_pm( nlevel, pm, ipar )
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
 
    call reader_of( nlevel, pm, ipar )
    call create_mg_levels( nlevel, pm, ipar )

  end subroutine create_mg_pm

  !!! Read the polyMesh data from OpenFOAM reader
  !!! uses parmgen and creates the multigrid levels
  subroutine reader_of( nlevel, pm, ipar )
    use Wrap
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
    !!! Local variables
    integer :: ilvl
    !!! Make the meshes aware of the multi-grid 
    do ilvl = 1, nlevel
      if( ipar .eq. enable_parallel_ ) then
        pm(ilvl)%parallel = .true.
      end if
      pm(ilvl)%nlevel = nlevel
      pm(ilvl)%ilevel = ilvl
      nullify(pm(ilvl)%x)
    end do
    !!! Read in the finest mesh from OpenFOAM
    call init_of_mesh( ipar )
    call get_pm_sizes( pm(nlevel)%nnode, pm(nlevel)%nface, &
                       pm(nlevel)%ninternalface, &
                       pm(nlevel)%ncell, pm(nlevel)%npatch )
    call allocate_pm( pm(nlevel) )
    call get_pm_nodes( pm(nlevel)%nnode, pm(nlevel)%x )
!!! Dummy face data
    call get_pm_faces( pm(nlevel)%nface, pm(nlevel)%ninternalface, &
                       pm(nlevel)%facelr, pm(nlevel)%facenode )
    call colour_pm_faces( pm(nlevel)%ninternalface, pm(nlevel)%ncell, &
                          pm(nlevel)%facelr, pm(nlevel)%facenode, &
                          pm(nlevel)%ncolour, pm(nlevel)%colourxadj )
    call get_pm_patches( pm(nlevel)%npatch, pm(nlevel)%patchdata )
    !!! Calculate the metrics
    call mesh_metrics( pm(nlevel) )
    call close_of_mesh()
  end subroutine reader_of

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! Create the Multi-Grid levels  !!
  !!!! using mgridgen wrapper        !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine create_mg_levels( nlevel, pm, ipar )
    use Wrap
    implicit none
    integer, intent(in)           :: nlevel, ipar
    type(polyMesh), intent(inout) :: pm(nlevel)
    !!! Local variables
    integer :: ilvl, minsize=1, maxsize=10, &
               options(4) = (/ 4, 6, 128, 3 /), nmoves
    type(crsGraph) :: gr
 
    !!! Multi-grid mgridgen
    do ilvl = nlevel, 2, -1
      call pm_to_graph( pm(ilvl), gr )
      call MGridGen_f90( gr%nvtx, gr%xadj, gr%vvol, gr%vsurf, gr%adjncy, &
                         gr%adjwgt, minsize, maxsize, options, nmoves, gr%npart, &
                         pm(ilvl)%mgpart )
!      call fix_mg_degeneracy( pm(ilvl) )
      call build_pm_coarse( pm(ilvl), pm(ilvl - 1), gr )
      call mesh_metrics( pm(ilvl - 1) ) 
    end do

  end subroutine create_mg_levels

  subroutine build_pm_coarse( pmf, pmc, gr )
    implicit none
    type(polyMesh) :: pmf, pmc
    type(crsGraph) :: gr
    !! Local variables
    integer, allocatable :: tmplr(:,:)
    integer :: ninternalface, nbface, iface

    !! Assign grid level infomration
    pmc%ilevel = pmf%ilevel - 1
    pmc%nlevel = pmf%nlevel
    !! Find unique internal faces in coarse mesh
    allocate(tmplr(lr_, pmf%nface))
    ninternalface = 0
    do iface = 1, pmf%nface
      tmplr(lcell_, iface) = pmf%mgpart(pmf%facelr(lcell_, iface)) + 1
      if( iface .le. pmf%ninternalface ) then
        tmplr(rcell_, iface) = pmf%mgpart(pmf%facelr(rcell_, iface)) + 1
        if( tmplr(lcell_, iface) .ne. tmplr(rcell_, iface) ) then
          ninternalface = ninternalface + 1
        end if
      end if
    end do
    !! Copy unique faces in coarse grid numbering
    pmc%nnode         = pmf%nnode
    pmc%ninternalface = ninternalface
    pmc%ncell         = gr%npart
    nbface            = pmf%nface - pmf%ninternalface
    pmc%nface         = pmc%ninternalface + nbface
    pmc%npatch        = pmf%npatch
    !! Allocate the polyMesh
    call allocate_pm( pmc )
    !! Internal faces counters
    ninternalface = 0
    pmc%facelr    = 0
    !! Copy temp face data to coarse mesh
    do iface = 1, pmf%ninternalface
      if( tmplr(lcell_, iface) .ne. tmplr(rcell_, iface) ) then
        ninternalface = ninternalface + 1
        pmc%facelr(:, ninternalface)   = tmplr(:, iface)
        pmc%facenode(:, ninternalface) = pmf%facenode(:, iface)
      end if
    end do
    !! Boundary faces
    do iface = pmf%ninternalface + 1, pmf%nface
      ninternalface = ninternalface + 1
      pmc%facelr(:, ninternalface)   = tmplr(:, iface)
      pmc%facenode(:, ninternalface) = pmf%facenode(:, iface)
    end do
    !! Assign node pointer
    pmc%x         => pmf%x
    pmc%patchdata => pmf%patchdata
    !! Nullify parallel data
    nullify(pmc%cellgid)
    nullify(pmc%facegid)
    nullify(pmc%nodegid)
    deallocate(tmplr)

  end subroutine build_pm_coarse

  !!! Simpler wrapper for metrics
  subroutine mesh_metrics( pm )
    use Wrap
    implicit none
    type( polyMesh ) :: pm
    !!! Tapenade diff function
    call mesh_metrics_tapenade &
         ( pm%nnode, pm%nface, &
           pm%ninternalface, pm%ncell, &
           pm%x, pm%facelr, &
           pm%facenode, pm%cv, pm%cc, &
           pm%dn, pm%fs, pm%fc )
    !!! Check metrics only if fine mesh
    if( pm%nlevel .eq. pm%ilevel ) then
      call check_metrics &
           ( pm%ncell, pm%nface, pm%cv, &
             pm%cc, pm%fc, pm%fs, pm%dn )
    end if

  end subroutine mesh_metrics 

  subroutine allocate_pm( pm )
    implicit none
    type(polyMesh), intent(inout) :: pm
    !!! Connectivity
    allocate( pm%facelr( lr_, pm%nface ) )
    allocate( pm%facenode( quadp1_, pm%nface ) )
    !!! Metrics
    allocate( pm%cc( dim_, pm%ncell ) )
    allocate( pm%cv( pm%ncell ) )
    allocate( pm%dn( dim_, pm%nface ) )
    allocate( pm%fc( dim_, pm%nface ) )
    allocate( pm%fs( pm%nface ) )
    !!! Multigrid connectivity
    allocate( pm%mgpart( pm%ncell ) )
    !!! Allocate nodes only to finest level
    !!! and share it across levels
    if( pm%ilevel .eq. pm%nlevel ) then
      allocate( pm%x( dim_, pm%nnode ) )
      allocate( pm%patchdata( pproc_, pm%npatch) )
      !!! Parallel data
      if( pm%parallel .eqv. .true. ) then
        allocate( pm%cellgid( pm%ncell ) )
        allocate( pm%facegid( pm%nface ) )
        allocate( pm%nodegid( pm%nnode ) )
      end if
    end if

  end subroutine allocate_pm

  subroutine mesh_metrics_tapenade( &
    nnode, nface, ninternalface, ncell, &
    x, facelr, facenode, &
    cv, cc, dn, fs, fc )
    implicit none
    integer, intent(in)         :: nnode, nface, ninternalface, ncell
    real(kind=8), intent(in)    :: x(dim_, nnode)
    integer, intent(in)         :: facelr(lr_, nface), facenode(quadp1_, nface)
    real(kind=8), intent(inout) :: cv(ncell), cc(dim_, ncell)
    real(kind=8), intent(inout) :: dn(dim_, nface), fs(nface), fc(dim_, nface)
    !!! Local variables
    integer :: iface, inode, i, il, ir, icell
    real(kind=8) :: r1(dim_), r2(dim_), r3(dim_), r4(dim_), fvol, fvolc(dim_)

    !!! Zero out everything    
    cv = 0.0d0; cc = 0.0d0
    !!! Calculate the face centroids, area, and unit normal vector
    do iface = 1, nface
      r1 = x(:, facenode( 2, iface ))
      r2 = x(:, facenode( 3, iface ))
      r3 = x(:, facenode( 4, iface ))
      il = facelr( lcell_, iface)
      ir = facelr( rcell_, iface)
!!! Planar face
      if( facenode(1, iface) .eq. tri_ ) then
        !!! Face centroid
        fc(:, iface) = oneby3_ * ( r1 + r2 + r3 )
        !!! Face normal and area
        dn(:, iface) = 0.50d0 * cross_prod( (r2 - r1), (r3 - r1) )
        fs(iface)    = sqrt( sum( dn(:, iface) * dn(:, iface) ) )
        dn(:, iface) = dn(:, iface) / fs(iface)
        !!! Face cell centroid contribution
        fvolc = sum(r1 * r1) + sum(r2 * r2) + sum(r3 * r3) + &
                sum(r1 * r2) + sum(r2 * r3) + sum(r1 * r3)
        fvolc = fvolc * oneby12_ * cross_prod((r2 - r1), (r3 - r1))
!!! Ruled face
      else if( facenode(1, iface) .eq. quad_ ) then
        r4 = x(:, facenode( 5, iface ))
        !!! Face centroid
        fc(:, iface) = 0.250d0 * ( r1 + r2 + r3 + r4 )
        !!! Face normal and area
        dn(:, iface) = 0.50d0 * cross_prod( (r3 - r1), (r4 - r2) )
        fs(iface)    = sqrt( sum( dn(:, iface) * dn(:, iface) ) )
        dn(:, iface) = dn(:, iface) / fs(iface)
        !!! Face cell centroid contribution (Gauss quadrture)
        fvolc = 0.250d0 * ( func_1234(r1, r2, r3, r4, xi1_, eta1_) + &
                            func_1234(r1, r2, r3, r4, xi1_, eta2_) + &
                            func_1234(r1, r2, r3, r4, xi2_, eta1_) + &
                            func_1234(r1, r2, r3, r4, xi2_, eta2_) )
!!! Error unknow face type
      else
        write(*,*) "Something terribly wrong ... Cannot have more than 4 nodes for a face"
        stop
      end if
!!! Add CC/CV contribution of face
      !!!! Face volume contribution
      fvol     = oneby3_ * sum( fc(:, iface) * dn(:, iface) ) * fs(iface)
      cv(il)   = cv(il)   + fvol
      cc(:,il) = cc(:,il) + fvolc
      if( iface .le. ninternalface ) then
        cv(ir)   = cv(ir)   - fvol
        cc(:,ir) = cc(:,ir) - fvolc
      end if
    end do
    do icell = 1, ncell
      cc(:,icell) = 0.50d0 * cc(:,icell) / cv(icell)
    end do

  end subroutine mesh_metrics_tapenade

  !!! Construct the CSR graph from the edge2nodes data, which
  !!! are input to the METIS partitioner routines 
  subroutine pm_to_graph( pm, gr )
    use mpi
    implicit none
    type(polyMesh), intent(in)  :: pm
    type(crsGraph), intent(out) :: gr
    !!! Local subroutine variables
    integer :: ierr, i
    !!! Some simple memory checks
    if( allocated(gr%xadj) .eqv. .true.   ) deallocate(gr%xadj)
    if( allocated(gr%adjncy) .eqv. .true. ) deallocate(gr%adjncy)
    if( allocated(gr%adjwgt) .eqv. .true. ) deallocate(gr%adjwgt)
    if( allocated(gr%vsurf) .eqv. .true.  ) deallocate(gr%vsurf)
    !!! Allocations
    gr%nvtx = pm%ncell
    allocate( gr%xadj(gr%nvtx + 1), &
              gr%adjncy(lr_ * pm%ninternalface), &
              gr%adjwgt(lr_ * pm%ninternalface), &
              gr%vsurf(gr%nvtx), gr%vvol(gr%nvtx) )
    gr%xadj   = 0
    gr%adjncy = -1
    gr%vsurf  = 0.0d0
    gr%vvol   = pm%cv
    !!! First pass form the sizes
    do i = 1, pm%ninternalface
      gr%xadj(pm%facelr(lcell_, i)) = gr%xadj(pm%facelr(lcell_, i)) + 1
      gr%xadj(pm%facelr(rcell_, i)) = gr%xadj(pm%facelr(rcell_, i)) + 1
    end do
    gr%xadj = cshift(gr%xadj, gr%nvtx)
    gr%xadj(1) = 0
    !write(*,*) gr%xadj
    !!! All boundary face surface area (attach to cell)
    do i = pm%ninternalface + 1, pm%nface
      gr%vsurf(pm%facelr(lcell_, i)) = gr%vsurf(pm%facelr(lcell_, i)) + pm%fs(i)
    end do
    !!! Form the xadj offsets (Can replace with intrisic?)
    do i = 1, gr%nvtx
      gr%xadj(i+1) = gr%xadj(i) + gr%xadj(i+1)
    end do
    !!! Check size match
    if( gr%xadj(gr%nvtx + 1) .ne. (lr_ * pm%ninternalface) ) then
      write(*,*) 'Error: Sizes of facelr xadj do not match ninternalface', &
        gr%xadj(gr%nvtx + 1), 2 * pm%ninternalface
    end if
    !!! Second pass form the adjncy
    do i = 1, pm%ninternalface
      !!! Push the left face data
      gr%xadj(pm%facelr(lcell_, i)) = gr%xadj(pm%facelr(lcell_, i)) + 1
      gr%adjncy(gr%xadj(pm%facelr(lcell_, i))) = pm%facelr(rcell_, i) - 1
      gr%adjwgt(gr%xadj(pm%facelr(lcell_, i))) = pm%fs(i)
      !!! Push the right face data
      gr%xadj(pm%facelr(rcell_, i)) = gr%xadj(pm%facelr(rcell_, i)) + 1
      gr%adjncy(gr%xadj(pm%facelr(rcell_, i))) = pm%facelr(lcell_, i) - 1
      gr%adjwgt(gr%xadj(pm%facelr(rcell_, i))) = pm%fs(i)
    end do
    gr%xadj = cshift(gr%xadj, gr%nvtx)
    gr%xadj(1) = 0
    !write(*,*) gr%xadj
    !!! Check to see if the adjncy array is formed correctly
    if(count(gr%adjncy .lt. 0) .ne. 0) then
      write(*,*) 'Error: Adjncy array from edge2nodes not constructed correctly'
    end if

  end subroutine pm_to_graph

  subroutine cell_face_xadj_adjncy(ninternalface, ncell, facelr, xadj, adjncy)
    implicit none
    integer, intent(in)  :: ninternalface, ncell
    integer, intent(in)  :: facelr(lr_, ninternalface)
    integer, intent(out) :: xadj(ncell + 1), &
                            adjncy(lr_ * ninternalface)
    !--------------------------------------------------------
    integer, dimension(ncell) :: counter
    integer :: iface, lcell, rcell, icell

    ! first sweep: count how many faces are connected to each cell
    counter = 0
    do iface = 1, ninternalface
      lcell = facelr(lcell_, iface)
      rcell = facelr(rcell_, iface)
      counter(lcell) = counter(lcell) + 1
      counter(rcell) = counter(rcell) + 1
    end do

    ! use the counter to initialise xadj
    xadj(1) = 1
    do icell = 1, ncell
      xadj(icell + 1) = xadj(icell) + counter(icell)
    end do
    ! reset the counter, we use it again
    counter = 0

    ! fill the adjncy array. For each cell, we fill
    ! the section that starts at the position marked by xadj
    do iface = 1, ninternalface
      lcell = facelr(lcell_, iface)
      rcell = facelr(rcell_, iface)
      adjncy(xadj(lcell) + counter(lcell)) = iface
      adjncy(xadj(rcell) + counter(rcell)) = iface
      counter(lcell) = counter(lcell) + 1
      counter(rcell) = counter(rcell) + 1
    end do
    
  end subroutine cell_face_xadj_adjncy

  subroutine face_colouring(ninternalface, ncell, facelr, ncolour, colour)
    implicit none
    integer, intent(in)    :: ninternalface, ncell
    integer, intent(inout) :: facelr(lr_, ninternalface), colour(ninternalface)
    integer, intent(out)   :: ncolour
    !--------------------------------------------------------
    integer :: iface, lcell, rcell, ileft, iright
    integer,dimension(lr_ * ninternalface) :: adjncy
    integer,dimension(ncell + 1) :: xadj
    integer :: mycolour, ncolor
    logical :: free
    !! Form the face-cell xadj and adjncy
    call cell_face_xadj_adjncy(ninternalface, ncell, facelr, xadj, adjncy)
    !! Init colours to zero
    colour = 0
    !! Loop over all faces
    do iface = 1, ninternalface
      lcell = facelr(lcell_,iface)
      rcell = facelr(rcell_,iface)
      ! go through the colors and check if we can still assign them to this face
      mycolour = 1
      do
        ! start assuming that we can use this color
        free = .true.
        ! check if the color is already taken by any other face connected to left cell
        do ileft = xadj(lcell), xadj(lcell + 1) - 1
          if(mycolour .eq. colour( adjncy(ileft) ) )then
            free = .false.
          end if
        end do
        ! check if the color is already taken by any other face connected to right cell
        do iright = xadj(rcell), xadj(rcell + 1) - 1
          if(mycolour .eq. colour( adjncy(iright) )) then
            free = .false.
          end if
        end do
        ! if the color wasn't used by any other face, we can still assign it.
        ! exit here and don't change c any further
        if(free) then
          exit
        ! if it isn't free any more, we try the next color. Go back
        ! and do another cycle.
        else
          mycolour = mycolour + 1
          cycle
        end if
      end do
      ! if the loop terminated, it found a free color and stored it in mycolour. We assign
      ! this color to i^{th} face now.
      colour(iface) = mycolour
    end do
    !! Total number of colours
    ncolour = maxval(colour)
    !! Check if everything went well
    if( minval(colour) .eq. 0 ) then
      write(*,*) "Error: in colour module setting ncolour to zero"
      ncolour = 0
    end if

  end subroutine face_colouring

  subroutine colour_pm_faces( ninternalface, ncell, facelr, &
                              facenode, ncolour, colourxadj )
    implicit none 
    integer, intent(in) :: ninternalface, ncell
    integer, intent(inout) :: facelr(lr_, ninternalface), &
                              facenode(quadp1_, ninternalface)
    integer, intent(out) :: ncolour
    integer, intent(out), allocatable :: colourxadj(:)
    !!! Local vars
    integer :: iface, icolour
    integer, allocatable :: colour(:), dummylr(:,:), &
                            dummynode(:,:)
!!! Colouring of faces
    allocate(dummylr(lr_, ninternalface), &
             dummynode(quadp1_, ninternalface), &
             colour(ninternalface))
!!! Copy all internal faces to dummy arrays to calculate colouring
    dummylr(:, 1:ninternalface) = facelr(:, 1:ninternalface)
    dummynode(:, 1:ninternalface) = facenode(:, 1:ninternalface)
    call face_colouring( ninternalface, ncell, dummylr, ncolour, colour )
    write(*,*) "Total colours = ", ncolour
!!! Construct colour xadj array
    allocate(colourxadj(ncolour + 1))
    colourxadj = 0
    do iface = 1, ninternalface
      colourxadj(colour(iface) + 1) = colourxadj(colour(iface) + 1) + 1
    end do
    colourxadj(1) = 1
    do icolour = 1, ncolour
      colourxadj(icolour + 1) = colourxadj(icolour) + colourxadj(icolour + 1)
    end do
!!! Order faces by colour
    do iface = 1, ninternalface
      icolour = colour(iface)
      facelr(:, colourxadj(icolour))   = dummylr(:, iface)
      facenode(:, colourxadj(icolour)) = dummynode(:, iface)
      colourxadj(icolour) = colourxadj(icolour) + 1
    end do
    colourxadj = cshift(colourxadj, ncolour)
    colourxadj(1) = 1   
    deallocate( dummylr, dummynode, colour )
!!!!
  end subroutine colour_pm_faces

end module Mesh

