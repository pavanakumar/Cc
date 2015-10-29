/*! \file tecplotGridFile.h
**
**/
#ifndef TECPLOT_GRID_H

#define TECPLOT_GRID_H

#include "tecplotCommon.h"

namespace tec360 {

  class tecplotGridFile {

    public:
    tecplotGridFile();

    tecplotGridFile( const char *file_name );

    ~tecplotGridFile();

    void Open( const char *file_name );

    void Close();
  
    void Write();

    void SetNumFaces( size_t num_faces );

    void SetNumCells( size_t num_cells );

    void SetNumNodes( size_t num_nodes );

    void SetNumFaceNodes( size_t num_face_nodes );

    void SetVarList( const char *list );

    INTEGER4 *FaceLeftCell();

    INTEGER4 *FaceRightCell();

    INTEGER4 *FaceNodes();

    std::vector<INTEGER4> &FaceNodesVector();

    INTEGER4 *FaceNodeCounts();

    double *X();

    double *Y();

    double *Z();

    void DryRun();

    private:
    INTEGER4 _file_type, _debug;
    bool _is_open;
    std::string _zone_name, _var_list, _file_name, _scratch_dir;
    std::string _title;
    INTEGER4 _z_type, _n_nodes, _n_cells, _n_faces;
    INTEGER4 _i_max, _j_max, _k_max;
    double   _sol_time;
    INTEGER4 _strand_id, _parent_zone, _data_packing, _n_face_conns, _f_n_mode;
    INTEGER4 *_pasiv_var_arr, *_val_loc_arr, *_var_shr_arr;
    INTEGER4 _shr_conn;
    INTEGER4 _n_face_nodes, _n_patch_conns, _n_patch_items;
    std::vector<double> _x, _y, _z;
    INTEGER4 _is_double;
    std::vector<INTEGER4> _face_node_counts, _face_nodes;
    std::vector<INTEGER4> _left_cell, _right_cell;

    void WriteGrid();

    void WriteZoneInfo();

    void WriteNodeInfo();

    void WriteFaceInfo();

  };

  /******* Class implementation **********/
  
  tecplotGridFile
  ::tecplotGridFile()
  : _file_type( GRID ), _debug(0), _is_open(false),
    _zone_name("POLY"), _var_list("X Y Z"),
    _scratch_dir("./"), _title("FUSE"),
    _z_type(FEPOLYHEDRON),
    _n_nodes(0), _n_cells(0), _n_faces(0),
    _i_max(0), _j_max(0),
    _k_max(0), _sol_time(0.0), _strand_id(0),
    _parent_zone(0), _data_packing(BLOCK), _n_face_conns(0),
     _f_n_mode(0),
    _pasiv_var_arr(0), _val_loc_arr(0), _var_shr_arr(0),
    _shr_conn(0), _n_face_nodes(0), _n_patch_conns(0),
    _n_patch_items(0), _is_double(1)
  { /* empty */ }

  
  tecplotGridFile
  ::tecplotGridFile( const char *file_name )
  : _file_type( GRID ), _debug(0), _is_open(false),
    _zone_name("POLY"), _var_list("X Y Z"), 
    _scratch_dir("./"),  _title("FUSE"), _z_type(FEPOLYHEDRON),
    _n_nodes(0), _n_cells(0), 
    _n_faces(0), _i_max(0), _j_max(0),
    _k_max(0), _sol_time(0.0), _strand_id(0),
    _parent_zone(0), _data_packing(BLOCK), 
    _n_face_conns(0), _f_n_mode(0), 
    _pasiv_var_arr(0), _val_loc_arr(0), _var_shr_arr(0),
     _shr_conn(0), _n_face_nodes(0),
    _n_patch_conns(0), 
    _n_patch_items(0), _is_double(1)
  {
    Open( file_name );
  }
 
  void tecplotGridFile
  ::Open( const char *file_name )
  {
    assert
    ( 
      TECINI112
      ( 
        const_cast<char *>( _title.c_str() ),
        const_cast<char *>( _var_list.c_str() ),
        const_cast<char*>( file_name ),
        const_cast<char*>( _scratch_dir.c_str() ),
        &_file_type, &_debug, &_is_double
      ) == 0 
    );
    _is_open = true;
  }

  tecplotGridFile
  ::~tecplotGridFile()
  {
    Close();
  }

  
  void tecplotGridFile
  ::Close()
  {
    if( _is_open ) {
      assert
      ( 
        TECEND112() 
        == 0 
      );
      _is_open = false;
    }
  }

  
  void tecplotGridFile
  ::Write()
  { 
      assert( _is_open );
      WriteGrid();
  }

  void tecplotGridFile
  ::WriteGrid()
  { 
    WriteZoneInfo();
    WriteNodeInfo();
    WriteFaceInfo();
  }

  void tecplotGridFile
  ::WriteZoneInfo()
  {
    assert
    (
      TECZNE112
      (
        const_cast<char *>( _zone_name.c_str() ),
        &_z_type,
        &_n_nodes,
        &_n_cells,
        &_n_faces,
        &_i_max,
        &_j_max,
        &_k_max,
        &_sol_time,
        &_strand_id,
        &_parent_zone,
        &_data_packing,
        &_n_face_conns,
        &_f_n_mode,
        &_n_face_nodes,
        &_n_patch_conns,
        &_n_patch_items,
        _pasiv_var_arr,
        _val_loc_arr,
        _var_shr_arr,
        &_shr_conn
      ) == 0 
    );
  }

  
  void tecplotGridFile
  ::WriteNodeInfo()
  {
    assert
    (
      TECDAT112
      (
        &_n_nodes, &_x[0], &_is_double
      ) == 0 
    );
    assert
    (
      TECDAT112
      (
        &_n_nodes, &_y[0], &_is_double
      ) == 0 
    );
    assert
    (
      TECDAT112
      (
        &_n_nodes, &_z[0], &_is_double
      ) == 0 
    );
  }

  
  void tecplotGridFile
  ::WriteFaceInfo()
  {
    assert
    (
      TECPOLY112
      ( 
        &_face_node_counts[0], /* The face node counts array */
        &_face_nodes[0],      /* The face nodes array */
        &_left_cell[0],  /* The left elements array  */
        &_right_cell[0], /* The right elements array  */
        NULL,       /* No boundary connection counts */
        NULL,       /* No boundary connection elements */
        NULL
      ) == 0 
    );      /* No boundary connection zones */
  }

/*** Set functions ****/
  
  void tecplotGridFile
  ::SetNumFaces( size_t num_faces )
  {
    _n_faces = num_faces; 
    _face_node_counts.resize(_n_faces);
    _left_cell.resize(_n_faces);
    _right_cell.resize(_n_faces);
  }
  
  
  void tecplotGridFile
  ::SetNumCells( size_t num_cells )
  { 
    _n_cells = num_cells;
  }

  
  void tecplotGridFile
  ::SetNumNodes( size_t num_nodes )
  {
    _n_nodes = num_nodes;
    _x.resize( _n_nodes );
    _y.resize( _n_nodes );
    _z.resize( _n_nodes );
  }

  
  void tecplotGridFile
  ::SetNumFaceNodes( size_t num_face_nodes )
  {
    _n_face_nodes = num_face_nodes;
    _face_nodes.resize( _n_face_nodes );
  }

  
  void tecplotGridFile
  ::SetVarList( const char *list )
  {
    _var_list.clear();
    _var_list = std::string(list);
  }

/*** Data access function ****/

  
  INTEGER4 *tecplotGridFile
  ::FaceLeftCell()
  {
    return &_left_cell[0];  
  }

  
  INTEGER4 *tecplotGridFile
  ::FaceRightCell()
  {
    return &_right_cell[0];  
  }

  
  INTEGER4 *tecplotGridFile
  ::FaceNodes()
  {
    return &_face_nodes[0];
  }

  
  std::vector<INTEGER4> &tecplotGridFile
  ::FaceNodesVector()
  {
    return _face_nodes;
  }


  
  INTEGER4 *tecplotGridFile
  ::FaceNodeCounts()
  {
    return &_face_node_counts[0];  
  }

  
  double *tecplotGridFile
  ::X()
  {
    return &_x[0];  
  }

  
  double *tecplotGridFile
  ::Y()
  {
    return &_y[0];  
  }

  
  double *tecplotGridFile
  ::Z()
  {
    return &_z[0];  
  }


/***** Dont look at this dry run stuff its just a test code *****/

  
  void tecplotGridFile
  ::DryRun()
  {
    SetNumNodes(5);
    SetNumCells(1);
    SetNumFaces(5);
    SetNumFaceNodes(16);

    _x[0] = 0;
    _y[0] = 0;
    _z[0] = 0;

    _x[1] = 1;
    _y[1] = 1;
    _z[1] = 2;

    _x[2] = 2;
    _y[2] = 0;
    _z[2] = 0;

    _x[3] = 2;
    _y[3] = 2;
    _z[3] = 0;

    _x[4] = 0;
    _y[4] = 2;
    _z[4] = 0;

    _face_node_counts[0] = 3;
    _face_node_counts[1] = 3;
    _face_node_counts[2] = 3;
    _face_node_counts[3] = 3;
    _face_node_counts[4] = 4;
   
    /* Face Nodes for Face 1 */
    _face_nodes[0]  = 1;
    _face_nodes[1]  = 2;
    _face_nodes[2]  = 3;

    /* Face Nodes for Face 2 */
    _face_nodes[3]  = 3;
    _face_nodes[4]  = 2;
    _face_nodes[5]  = 4;

    /* Face Nodes for Face 3 */
    _face_nodes[6]  = 5;
    _face_nodes[7]  = 2;
    _face_nodes[8]  = 4;

    /* Face Nodes for Face 4 */
    _face_nodes[9]  = 1;
    _face_nodes[10] = 2;
    _face_nodes[11] = 5;

    /* Face Nodes for Face 5 */
    _face_nodes[12] = 1;
    _face_nodes[13] = 5;
    _face_nodes[14] = 4;
    _face_nodes[15] = 3;

    _left_cell[0]  = 1;
    _left_cell[1]  = 1;
    _left_cell[2]  = 0;
    _left_cell[3]  = 0;
    _left_cell[4]  = 0;

    _right_cell[0]  = 0;
    _right_cell[1]  = 0;
    _right_cell[2]  = 1;
    _right_cell[3]  = 1;
    _right_cell[4]  = 1;

    WriteGrid();
  }  

} // End of tecplot namespace

#endif

