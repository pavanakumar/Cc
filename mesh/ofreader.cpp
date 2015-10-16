#include "of_reader.h"

/*************************************
   Opens the mesh
   1. Serial   if ipar != __enable_par
   2/ Parallel if ipar == __enable_par
***************************************/
void init_of_mesh( int *ipar ) { 
  char *arg[] = 
  {
    const_cast<char *>("ddd\0"),
    const_cast<char *>("-parallel\0")
  };
  /// Init the OpenFOAM mesh
  if( *ipar == __enable_par )
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
void get_pm_nodes( int *nNodes, double *x ) {
  for( int i = 0; i < *nNodes; ++i ) {
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
      facenodes[ i * 5 + j ] = global_of_mesh->mesh()->faces()[i][j] + 1; // +1 FORTRAN indexing
    }
  }
  /// Form facelr for all internal cells
  for( int i=0; i<*ninternalface; ++i ) {
    facelr[ i * 2 ]     = global_of_mesh->mesh()->owner()[i] + 1; // +1 FORTRAN indexing
    facelr[ i * 2 + 1 ] = global_of_mesh->mesh()->neighbour()[i] + 1; // +1 FORTRAN indexing
  } 
  /// Patch faces only have owner and no neighbour cell
  for( int i=0; i<*npatch; ++i ) {
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
  for( int i=0; i<*npatch; ++i ) {
    patchdata[ i * 4 + 0 ] = global_of_mesh->mesh()->boundaryMesh()[i].start();
    patchdata[ i * 4 + 1 ] = global_of_mesh->mesh()->boundaryMesh()[i].size();
    //// To do ...
//    patchdata[ i * 4 + 2 ] = global_of_mesh->mesh()->boundaryMesh()[i].size();
//    if( global_of_mesh->mesh()->boundaryMesh()[i].type() == "processor" )
//      patchdata[ i * 4 + 3 ] = global_of_mesh->mesh()->boundaryMesh()[i].size();
  }
}

