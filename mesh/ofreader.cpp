#include "ofreader.h"

/******************************************
   Opens the mesh
   1. Serial   if ipar != _enable_parallel
   2. Parallel if ipar == _enable_parallel
*******************************************/
void init_mesh_api( int *ipar ) { 
  char *arg[] = 
  {
    const_cast<char *>("ddd\0"),
    const_cast<char *>("-parallel\0")
  };
  /// Init the OpenFOAM mesh
  if( *ipar == _enable_parallel )
    global_of_mesh = new __c_ofreader( 2, arg );
  else
    global_of_mesh = new __c_ofreader( 1, arg );
}

/*************************************
   Frees the mesh and cleares memory
   Note: This will not close the Time
   object, which will call MPI_Finalize
   so use close_mesh_api() to do
   complete clean up of API interface.
***************************************/
void free_mesh_api() {
  delete global_of_mesh;
}

/*************************************
   Closes the mesh and cleares memory
   Note: This will also close the Time
   object, which will call MPI_Finalize
   so use with caution.
***************************************/
void close_mesh_api() {
  delete global_of_mesh;
}

/*****************************************
  Read sizes from mesh and pass it to Cc
******************************************/
void get_pm_sizes( int *nnode, int *nface, int *ninternalface,
                   int *nedge, int *ninternaledge, int *ncell, 
                   int *npatch ) {
  *nnode = global_of_mesh->mesh()->nPoints();
  *nface = global_of_mesh->mesh()->nFaces();
  *ninternalface = global_of_mesh->mesh()->nInternalFaces();
  *nedge = global_of_mesh->mesh()->nEdges();
  *ninternaledge = global_of_mesh->mesh()->nInternalEdges();
  *ncell = global_of_mesh->mesh()->nCells();
  *npatch = global_of_mesh->mesh()->boundaryMesh().size();
}

/*****************************************
  Read nodes from mesh and pass it to Cc
******************************************/
void get_pm_nodes( int *nnode, double *x ) {
  for( int i = 0; i < *nnode; ++i ) {
    for( int j = 0; j < 3; ++j ) {
      x[ i * 3 + j ] = global_of_mesh->mesh()->points()[i][j];
    }
  }
}

/*****************************************
  Read faces from mesh and pass it to Cc
******************************************/
void get_pm_faces( int *nface, int *ninternalface,
                   int *nfacenode, int *facenode,
                   int *facelr ) {
  int size = 0, right;
  const int _quad = 4, _left = 0, _right = 1;
  ///// Form facenode
  for( int i = 0; i < *nface; ++i ) {
    size = global_of_mesh->mesh()->faces()[i].size();
    nfacenode[i] = size;
    /// Face nodes loop
    for( int j = 0; j < size; ++j )
      facenode[ i * _quad + j ] = global_of_mesh->mesh()->faces()[i][j] + 1; // +1 FORTRAN indexing
    ///// Fix facenode for planar faces     0   1   2   3
    ///// (last node repeats (for example) [13, 22, 37, 37]
    if( size == 3 ) // Planar face
      facenode[ i * _quad + _quad - 1 ] = facenode[ i * _quad + _quad - 2 ];
    /// Face L/R
    if( i < *ninternalface )
      right = global_of_mesh->mesh()->faceNeighbour()[i] + 1; // +1 FORTRAN indexing
    else
      right = 0; // The right face for a boundary is an invalid index
    facelr[ i * 2 + _left  ] = global_of_mesh->mesh()->faceOwner()[i] + 1; // +1 FORTRAN indexing
    facelr[ i * 2 + _right ] = right;
  }
}

/************************************************************
   Get the mesh edge data from OpenFOAM to Cc
*************************************************************/
void get_pm_edges( int *nedge, int *edgenode ) {
  const Foam::edgeList &edgL = global_of_mesh->mesh()->edges();
  for( int i = 0; i < *nedge; ++i ) {
    edgenode[ i * 2 + 0 ] = edgL[i][0] + 1;  /* +1 FORTRAN indexing */
    edgenode[ i * 2 + 1 ] = edgL[i][1] + 1; /* +1 FORTRAN indexing */
  }
}

/******************************************************
  Read bc patches from mesh and pass it to Cc
  Some constants:
  ---------------
  0 = patchstart
  1 = patchsize
  2 = patchtype
  3 = patchneigh (proc ID) - only for pocessor patches
********************************************************/
void get_pm_patches( int *npatch, int *patchdata ) {
  Foam::wordList basicTypes =
        global_of_mesh->mesh()->boundaryMesh().types();
  Foam::wordList physicalTypes =
        global_of_mesh->mesh()->boundaryMesh().physicalTypes();
  for( int i=0; i<*npatch; ++i ) {
    /// Get patch name for debugging purpose
    Foam::word patchName =
          global_of_mesh->mesh()->boundaryMesh()[i].name(); 
    //// Fix the patch values
    patchdata[ i * 4 + 0 ] = global_of_mesh->mesh()->boundaryMesh()[i].start() + 1; // +1 FORTRAN indexing
    patchdata[ i * 4 + 1 ] = global_of_mesh->mesh()->boundaryMesh()[i].size();
    patchdata[ i * 4 + 3 ] = 0; // Default value is 0 for processor neighbour
    ////  BC type  /////
    if( basicTypes[i] == "wall" )
      patchdata[ i * 4 + 2 ] = _wall_bc;
    else if( basicTypes[i] == "symmetryPlane" )
      patchdata[ i * 4 + 2 ] = _symmetry_bc;
    else if( basicTypes[i] == "empty" )
      patchdata[ i * 4 + 2 ] = _empty_bc;
    else if( basicTypes[i] == "patch" ) {  /// specially treat inlet/outlet/riemann
      if( physicalTypes[i] == "riemann" ) {
        patchdata[ i * 4 + 2 ] = _riemann_bc;
      }else if( physicalTypes[i] == "inlet" ) {
        patchdata[ i * 4 + 2 ] = _inlet_bc;
      }else if( physicalTypes[i] == "outlet" ) {
        patchdata[ i * 4 + 2 ] = _outlet_bc;
      }else if( physicalTypes[i] == "symmetry" ) {
        patchdata[ i * 4 + 2 ] = _symmetry_bc;
      }else if( basicTypes[i] == "empty" ) {
        patchdata[ i * 4 + 2 ] = _empty_bc;
      }else {
        Foam::Info << "Unknown physical patch type \"" << physicalTypes[i] << "\"\n";
        if( physicalTypes[i] == "" )
          Foam::Info << "Did you forget to specify a physicalType "
                     << "in constant\\polyMesh\\boundary for patch \""
                     << patchName << "\"?\n";
        Foam::FatalError.exit();
      }
    }
    else if( basicTypes[i] == "processor" ) {
      patchdata[ i * 4 + 2 ] = _processor_bc;
      processorPolyPatch &myProcPatch = 
        const_cast<processorPolyPatch &>
        (
          refCast<const processorPolyPatch>
            ( global_of_mesh->mesh()->boundaryMesh()[i] )
        );
      /////  Neighbour processor number
      patchdata[ i * 4 + 3 ] = myProcPatch.neighbProcNo();
    }
    else {
      Foam::Info << "Unsupported bc type " << basicTypes[i] << "\n";
      Foam::FatalError.exit();
    }
  }
}

/************************************************************
   Get the wall-distance values
*************************************************************/
void get_pm_walldist( int *ncell, double *walldist ) {

}

/************************************************************
   Check the mesh metrics if they are calculated correctly
*************************************************************/
void check_metrics
( 
  int *ncell, int *nface, double *cv,
  double *cc, double *fc, double *fs,
  double *dn ) {
//  std::cerr << "Epsilon = " << std::numeric_limits<double>::epsilon() << "\n";
//  std::cerr << "Round err  = " << std::numeric_limits<double>::round_error() << "\n";
  const double eps = 1.0e-10;//std::numeric_limits<double>::epsilon() * 10.0;
#if 0
  Foam::Info << "Face normal\n";
  for( int i=0; i<*nface; ++i ) {
    Foam::Info << i << "  "
               << global_of_mesh->mesh()->Sf()[i][0] << "  "
               << global_of_mesh->mesh()->Sf()[i][1] << "  "
               << global_of_mesh->mesh()->Sf()[i][2] << "  "
               << dn[ i * 3 + 0 ] * fs[i] << "  "
               << dn[ i * 3 + 1 ] * fs[i] << "  "
               << dn[ i * 3 + 2 ] * fs[i] << "\n";
  }
#endif
#if 1
  for( int i=0; i<*ncell; ++i ) {
    double cv_chk = std::abs( (cv[i] - global_of_mesh->mesh()->V()[i]) / cv[i] );
    if( cv_chk > eps ) {
      Foam::Info << "volume error  " << i << " = " << cv_chk << "\n";
      Foam::Info << "Error: Metrics not calculated correctly (cell volume) : "
                 << cv_chk << "  " << cv[i] << " != " << global_of_mesh->mesh()->V()[i] << "\n";
    }
  }
#endif
#if 0
    for( int j=0; j<3; ++j ) {
      cv_chk = std::abs( (cc[i * 3 + j] - global_of_mesh->mesh()->C()[i][j]) );// / cc[i * 3 + j] );
      if( cv_chk > eps ) {
        Foam::Info << "Error: Metrics not calculated correctly (cell centroid) : "
                   << cv_chk << " " << cc[i * 3 + j]  << "\n";
      }
    }
#endif
#if 0
    std::cout << "volume   " << i << " = " << cv[i] << "  " << global_of_mesh->mesh()->V()[i] << "\n";
    std::cout << "centroid " << i << " = " << cc[i*3] << "  " << cc[i*3+1] << "  " << cc[i*3+2] << "  "
              << global_of_mesh->mesh()->C()[i][0] << "  " << global_of_mesh->mesh()->C()[i][1] << "  "
              << global_of_mesh->mesh()->C()[i][2] << "\n";
  }
#endif
}

/************************************************************
   Get the cell gids from OpenFOAM parallel mesh
*************************************************************/
void get_cellgid( int *ncell, int *cellgid ) {

}

/************************************************************
   Get the node gids from OpenFOAM parallel mesh
*************************************************************/
void get_nodegid( int *nnode, int *nodegid ) {

}

/************************************************************
   Get the face gids from OpenFOAM parallel mesh
*************************************************************/
void get_facegid( int *nface, int *facegid ) {

}

/************************************************************************
   Get the extra data required for forming the dual mesh
   ncell    : Total number of cells in the mesh (for check)
   celltype : Value ==> Type of cell (refer to OpenFOAM cellShape class)
              0 ==> unknown
              1 ==> tet
              2 ==> hex
              3 ==> prism
              4 ==> pyr
              5 ==> wedge
              6 ==> tetWedge
   cellnode : Node forming the cell (size==8)
              (refer to cellShape class for convention)
   cellface : Faces forming the cell (size==6)
              (refer to cellShape class for convention)
   celledge : Edges forming the cell (size==12)
              (refer to cellShape class for convention)
***********************************************************************/
void get_pm_extra( int *ncell, int *celltype, int *cellnode,
                   int *cellface, int *celledge ) {
  /// Constants defining the array limits
  const int MAX_CELL_NODE = 8, MAX_CELL_FACE = 6, MAX_CELL_EDGE = 12;

  /// Allowed cell types
  const cellModel& tet      = *(cellModeller::lookup("tet"));
  const cellModel& hex      = *(cellModeller::lookup("hex"));
  const cellModel& prism    = *(cellModeller::lookup("prism"));
  const cellModel& pyr      = *(cellModeller::lookup("pyr"));
  const cellModel& wedge    = *(cellModeller::lookup("wedge"));
  const cellModel& tetWedge = *(cellModeller::lookup("tetWedge"));
  /// Get the cell shapes
  const cellShapeList& cellShapes =
        global_of_mesh->mesh()->cellShapes();
  ///
  forAll(cellShapes, cellI) {
    /// Easy access
    const cellModel& model = cellShapes[cellI].model();
    const labelList edges = cellShapes[cellI].meshEdges
                             ( global_of_mesh->mesh()->edges(), 
                               global_of_mesh->mesh()->cellEdges()[cellI] );
    const labelList faces  = cellShapes[cellI].meshFaces
                             ( global_of_mesh->mesh()->faces(),
                               global_of_mesh->mesh()->cells()[cellI] );
    /// 
    for( int i = 0; i < cellShapes[cellI].size(); ++i )
      cellnode[ cellI * MAX_CELL_NODE + i] = cellShapes[cellI][i] + 1; // +1 for FORTRAN indexing
    for( int i = 0; i < faces.size(); ++i )
      cellface[ cellI * MAX_CELL_FACE + i] = faces[i] + 1; // +1 for FORTRAN indexing
    for( int i = 0; i < edges.size(); ++i )
      celledge[ cellI * MAX_CELL_EDGE + i] = edges[i] + 1; // +1 for FORTRAN indexing
    /// Check for correct cell model and trf the data
    if     ( model == tet )      celltype[cellI] = 1;
    else if( model == hex )      celltype[cellI] = 2;
    else if( model == prism )    celltype[cellI] = 3;
    else if( model == pyr )      celltype[cellI] = 4;
    else if( model == wedge )    celltype[cellI] = 5;
    else if( model == tetWedge ) celltype[cellI] = 6;
    else                         celltype[cellI] = 0;
  }
}

