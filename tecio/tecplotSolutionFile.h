/*! \file tecplotSolution.h
**
**/
#ifndef TECPLOT_SOLUTION_H

#define TECPLOT_SOLUTION_H

#include "tecplotCommon.h"

namespace tec360 {

  class tecplotSolutionFile {

    public:
    tecplotSolutionFile();

    tecplotSolutionFile( const char *file_name );

    ~tecplotSolutionFile();

    void Open( const char *file_name );

    void Close();
  
    void SetNumFaces( size_t num_faces );

    void SetNumCells( size_t num_cells );

    void SetNumNodes( size_t num_nodes );

    void SetVarList( const char *list );

    double *Buf();
   
    void SetNumComp( unsigned num_comp );

    void SetNumComp( unsigned num_comp, valueLocation loc );

    void WriteHeader();

    void WriteZoneInfo();

    void WriteVar();

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
    std::vector<double> _buf;
    INTEGER4 _is_double;
    unsigned _num_comp;

  };

  /******* Class implementation **********/
  
  tecplotSolutionFile
  ::tecplotSolutionFile()
  : _file_type( SOLUTION ), _debug(0), _is_open(false),
    _zone_name("POLY"), _var_list(""), _scratch_dir("./"),
    _title("FUSE"), _z_type(FEPOLYHEDRON),
    _n_nodes(0), _n_cells(0),
    _n_faces(0), _i_max(0), _j_max(0),
    _k_max(0), _sol_time(0.0), _strand_id(0),
    _parent_zone(0), _data_packing(BLOCK), _n_face_conns(0),
    _f_n_mode(0), 
    _pasiv_var_arr(0), _val_loc_arr(0), _var_shr_arr(0),
    _shr_conn(0), _n_face_nodes(0), _n_patch_conns(0), 
    _n_patch_items(0), _is_double(1), _num_comp(1)
  { /* empty */ }

  
  tecplotSolutionFile
  ::tecplotSolutionFile( const char *file_name )
  : _file_type( SOLUTION ), _debug(0), _is_open(false),
    _zone_name("POLY"), _var_list(""),
    _file_name( file_name ), _scratch_dir("./"),
    _title("FUSE"), _z_type(FEPOLYHEDRON),
    _n_nodes(0), _n_cells(0),
    _n_faces(0), _i_max(0), _j_max(0),
    _k_max(0), _sol_time(0.0), _strand_id(0),
    _parent_zone(0), _data_packing(BLOCK), _n_face_conns(0),
    _f_n_mode(0), 
    _pasiv_var_arr(0), _val_loc_arr(0), _var_shr_arr(0),
    _shr_conn(0), _n_face_nodes(0), _n_patch_conns(0),
    _n_patch_items(0), _is_double(1), 
    _num_comp(1)
  {
    SetNumComp( _num_comp );
  }

  void tecplotSolutionFile 
  ::Open( const char *file_name )
  {
    _file_name = std::string(file_name);
  }

  void tecplotSolutionFile 
  ::SetNumComp( unsigned num_comp )
  {
    _num_comp = num_comp;
    if( _val_loc_arr != 0 )
      delete[] _val_loc_arr;
    _val_loc_arr = new INTEGER4[ _num_comp ];
    for( unsigned i = 0; i < _num_comp; ++i )
      _val_loc_arr[i] = CELL_CENTER;
  }

  void tecplotSolutionFile 
  ::SetNumComp( unsigned num_comp, valueLocation loc )
  {
    _num_comp = num_comp;
    if( _val_loc_arr != 0 )
      delete[] _val_loc_arr;
    _val_loc_arr = new INTEGER4[ _num_comp ];
    for( unsigned i = 0; i < _num_comp; ++i )
      _val_loc_arr[i] = loc;
  }

  void tecplotSolutionFile
  ::WriteHeader()
  {
    assert
    ( 
#ifdef TECIO_2008
      TECINI111
#else
      TECINI112
#endif
      ( 
        const_cast<char *>( _title.c_str() ),
        const_cast<char *>( _var_list.c_str() ),
        const_cast<char*>( _file_name.c_str() ),
        const_cast<char*>( _scratch_dir.c_str() ),
        &_file_type, &_debug, &_is_double
      ) == 0
    );
    _is_open = true;
  }

  tecplotSolutionFile
  ::~tecplotSolutionFile()
  {
    Close();
  }
  
  void tecplotSolutionFile
  ::Close()
  {
    if( _is_open ) {
      assert
      ( 
#ifdef TECIO_2008
        TECEND111()
#else
        TECEND112() 
#endif
        == 0 
      );
      _is_open = false;
    }
  }

  void tecplotSolutionFile
  ::WriteZoneInfo()
  {
    assert
    (
#ifdef TECIO_2008
      TECZNE111
#else
      TECZNE112
#endif
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

  void tecplotSolutionFile
  ::WriteVar()
  {
    assert
    (
#ifdef TECIO_2008
      TECDAT111
#else
      TECDAT112
#endif
      (
        &_n_cells, &_buf[0], &_is_double
      ) == 0 
    );
  }
  
/*** Set functions ****/
  
  void tecplotSolutionFile
  ::SetNumFaces( size_t num_faces )
  {
    _n_faces = num_faces; 
  }
  
  
  void tecplotSolutionFile
  ::SetNumCells( size_t num_cells )
  { 
    _n_cells = num_cells;
    _buf.resize( _n_cells );
  }

  
  void tecplotSolutionFile
  ::SetNumNodes( size_t num_nodes )
  {
    _n_nodes = num_nodes;
  }

  
  void tecplotSolutionFile
  ::SetVarList( const char *list )
  {
    _var_list.clear();
    _var_list = std::string(list);
  }

/*** Data access function ****/

  double *tecplotSolutionFile
  ::Buf()
  {
    return &_buf[0];  
  }

} // End of tecplot namespace

#endif

