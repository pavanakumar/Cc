#include "ofreader.h"

/******************************************
   Opens the mesh
   1. Serial   if ipar != _enable_parallel
   2/ Parallel if ipar == _enable_parallel
*******************************************/
void init_of_mesh( int *ipar ) { 
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
   Closes the mesh and cleares memory
   Note: This will also close the Time
   object, which will call MPI_Finalize
   so use with caution.
***************************************/
void close_of_mesh() {
  delete global_of_mesh;
}
/*****************************************
  Read sizes from mesh and pass it to Cc
******************************************/
void get_pm_sizes( int *nnode, int *nface, int *ninternalface,
                   int *ncell, int *npatch ) {
  *nnode = global_of_mesh->mesh()->points().size();
  *ncell = global_of_mesh->mesh()->cells().size();
  *nface = global_of_mesh->mesh()->faces().size();
  *ninternalface = global_of_mesh->mesh()->owner().size();
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
                   int *facelr, int *facenodes ) {
  int size = 0;
  ///// Form facenodes
  for( int i=0; i<*nface; ++i ) {
    size = global_of_mesh->mesh()->faces()[i].size();
    facenodes[ i * 5 ] = size;
    for( int j=1; j <= size; ++j ) { /// Loop over number of face nodes
      facenodes[ i * 5 + j ] = global_of_mesh->mesh()->faces()[i][j-1] + 1; // +1 FORTRAN indexing
    }
  }
  /// Form facelr for all internal cells
  for( int i=0; i<*ninternalface; ++i ) {
    facelr[ i * 2 ]     = global_of_mesh->mesh()->owner()[i] + 1; // +1 FORTRAN indexing
    facelr[ i * 2 + 1 ] = global_of_mesh->mesh()->neighbour()[i] + 1; // +1 FORTRAN indexing
  } 
  /// Patch faces only have owner and no neighbour cell
  int npatch = global_of_mesh->mesh()->boundaryMesh().size();
  for( int i=0; i<npatch; ++i ) {
    int start = global_of_mesh->mesh()->boundaryMesh()[i].start();
    int size  = global_of_mesh->mesh()->boundaryMesh()[i].size();
    for( int j=0; j<size; ++j ) {
      facelr[ ( start + j ) * 2 ] = global_of_mesh->mesh()->boundaryMesh()[i].faceCells()[j] + 1; // +1 FORTRAN indexing
      facelr[ ( start + j ) * 2 + 1] = 0; /// Boundary patches have 0 as their neighbour cell ID (invalid)
    }
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
    else if( basicTypes[i] == "patch" ) {  /// specially treat inlet/outlet/riemann
      if( physicalTypes[i] == "riemann" ) {
        patchdata[ i * 4 + 2 ] = _riemann_bc;
      }else if( physicalTypes[i] == "inlet" ) {
        patchdata[ i * 4 + 2 ] = _inlet_bc;
      }else if( physicalTypes[i] == "outlet" ) {
        patchdata[ i * 4 + 2 ] = _outlet_bc;
      }else if( physicalTypes[i] == "symmetry" ) {
        patchdata[ i * 4 + 2 ] = _symmetry_bc;
      }else {
        Foam::Info << "Unknown physical patch type \"" << physicalTypes[i] << "\"\n";
        if( physicalTypes[i] == "" )
          Foam::Info << "Did you forget to specify a physicalType in constant\\polyMesh\\boundary for patch \"" << patchName << "\"?\n";
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
   Check the mesh metrics if they are calculated correctly
*************************************************************/
void check_metrics
( 
  int *ncell, int *nface, double *cv,
  double *cc, double *fc, double *fs,
  double *dn ) {
  std::cerr << "Epsilon = " << std::numeric_limits<double>::epsilon() << "\n";
  std::cerr << "Round err  = " << std::numeric_limits<double>::round_error() << "\n";
  const double eps = 1.0e-10;//std::numeric_limits<double>::epsilon() * 10.0;
  for( int i=0; i<*ncell; ++i ) {
    double cv_chk = std::abs( (cv[i] - global_of_mesh->mesh()->V()[i]) / cv[i] );
    if( cv_chk > eps ) {
      Foam::Info << "Error: Metrics not calculated correctly (cell volume) : "
                 << cv_chk << "\n";
    }
    for( int j=0; j<3; ++j ) {
      cv_chk = std::abs( (cc[i * 3 + j] - global_of_mesh->mesh()->C()[i][j]) );// / cc[i * 3 + j] );
      if( cv_chk > eps ) {
        Foam::Info << "Error: Metrics not calculated correctly (cell centroid) : "
                   << cv_chk << " " << cc[i * 3 + j]  << "\n";
      }
    }
  }
}


