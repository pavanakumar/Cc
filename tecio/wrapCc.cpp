#include "wrapCc.h"

void write_cc_tecio
(
  int *ilvl, int *rank,
  int *nnode, int *ncell, int *nface,
  int *ninternalface, double *xyz,
  int *facelr, int *facenode
){
  size_t face_node_adj_size = 0;
  std::stringstream cat;
  cat << "prefix_" << *ilvl << "_" << *rank << ".plt";
  tec360::tecplotGridFile tec_file( cat.str().c_str() );
  /// Set sizes
  tec_file.SetNumNodes( *nnode );
  tec_file.SetNumCells( *ncell );
  tec_file.SetNumFaces( *nface );
  std::cout << "nnode = "  << *nnode
            << " ncell = " << *ncell
            << " nface = " << *nface
            << " ninternalface = " 
            << *ninternalface << "\n";
  /// Set XYZ coordinates
  for( int i = 0; i < *nnode; ++i ) {
    tec_file.X()[i] = xyz[ i * 3 + 0 ];
    tec_file.Y()[i] = xyz[ i * 3 + 1 ];
    tec_file.Z()[i] = xyz[ i * 3 + 2 ];
  }
  /// Set L/R faces
  for( int i = 0; i < *nface; ++i ) {
    tec_file.FaceLeftCell()[i]   = facelr[ i * 2 ];
    tec_file.FaceRightCell()[i]  =
        ( i < *ninternalface ) ? facelr[ i * 2 + 1 ] : 0;
    tec_file.FaceNodeCounts()[i] = facenode[ i * 5 + 0 ];
    face_node_adj_size          += facenode[ i * 5 + 0 ];
  }
  std::cout << "face_node_adj_size = " << face_node_adj_size << "\n";

  /// Face-node adjncy data
  tec_file.SetNumFaceNodes( face_node_adj_size ); // Adj size
  size_t count = 0;
  for( int i = 0; i < *nface; ++i )
    for( int j = 1; j <= facenode[ i * 5 + 0 ]; ++j )
      tec_file.FaceNodes()[count++] = facenode[ i * 5 + j ];

  std::cout << "count = " << count << "\n";
  tec_file.Write();
  tec_file.Close();
}

