#ifndef __CC_OPENFOAM_READER__

#define __CC_OPENFOAM_READER__

#include "fvCFD.H"
#include "cellShape.H"
#include "cellModeller.H"
#include "argList.H"
#include "mpi.h"

#include <set>
#include <map>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>
#include <vector>

/**********************************************************
              Constants from FORTRAN code
**********************************************************/
const int _enable_parallel = 0;
const int _processor_bc    = 0,
          _wall_bc         = 1,
          _symmetry_bc     = 2,
          _inlet_bc        = 3,
          _outlet_bc       = 4,
          _riemann_bc      = 5,
          _empty_bc        = 6;
/**********************************************************
   Function prototypes for the FORTRAN wrapper
**********************************************************/
extern "C" {
  void init_mesh_api( int *ipar );
  void close_mesh_api ();
  void free_mesh_api  ();
  void get_pm_sizes   ( int *nnode, int *nface, int *ninternalface,
                        int *nedge, int *ninternaledge, int *ncell, 
                        int *npatch );
  void get_pm_nodes   ( int *nNodes, double *x );
  void get_pm_faces   ( int *nface, int *ninternalface,
                        int *nfacenode, int *facenode,
                        int *facelr );
  void get_pm_edges   ( int *nedge, int *edgenode );
  void get_pm_patches ( int *npatch, int *patchdata );
  void check_metrics  ( int *ncell, int *nface,
                        double *cv, double *cc,
                        double *fc, double *fs, double *dn );
  void get_cellgid    ( int *ncell, int *cellgid );
  void get_nodegid    ( int *nnode, int *nodegid );
  void get_facegid    ( int *nface, int *facegid );
  void get_pm_extra   ( int *ncell, int *celltype, int *cellnode,
                        int *cellface, int *celledge );
  void get_pm_walldist( int *ncell, double *walldist );
  void get_pm_shared_node_size
  (
    int *nprocs, int *nxadj,
    int *nadjncy, int *nsnlist
  );
  void get_pm_shared_node_shed
  (
    int *procs, int *xadj,
    int *adjncy, int *snlist
  );
};

/**************************************
   OpenFOAM reader - C wrapper class
***************************************/
struct __c_ofreader {
  ////// Constructor
  __c_ofreader( int nargs, char *args[] ) {
    /////// Commandline args parser
    _args = new Foam::argList( nargs, args );
    if ( !_args->checkRootCase() )
      Foam::FatalError.exit();
    ////// Runtime object
    Foam::Info<< "Create time\n" << Foam::endl;
    _time = new Foam::Time( Foam::Time::controlDictName, *_args );
    ////// Mesh object
    Foam::Info << "Create mesh for time = " << _time->timeName() << Foam::nl << Foam::endl;
    _mesh = new Foam::fvMesh
    (
       Foam::IOobject
       (
         Foam::fvMesh::defaultRegion,
         _time->timeName(),
         *_time,
         Foam::IOobject::MUST_READ
       )
    );
    if(Pstream::parRun()) shared_node_mpi();
  }

  void shared_node_mpi() {
    read_processor_node_gid();
    form_shared_node_list();
    form_shared_node_proc_list();
    form_shared_node_map();
  }

  void read_processor_node_gid() {
    _node_gid = new labelIOList
    (
      IOobject
      (
        "pointProcAddressing",
        _mesh->pointsInstance(),
        polyMesh::meshSubDir,
        *(_mesh),
        IOobject::READ_IF_PRESENT
      ),
      labelList(0)
    );
    forAll( *_node_gid, inode )  (*_node_gid)[inode]++; // +1 for FORTRAN indexing
  }

  void form_shared_node_list() {
    std::set<int> shared_node_set;
    Foam::wordList basicTypes = _mesh->boundaryMesh().types();
    forAll( basicTypes, itype ) {
      if( basicTypes[itype] == "processor" ) {
        auto myProcPatch = const_cast<processorPolyPatch &>
        (
          refCast<const processorPolyPatch>( _mesh->boundaryMesh()[itype] )
        );
        int begin = _mesh->boundaryMesh()[itype].start(),
            end   = begin + _mesh->boundaryMesh()[itype].size();
        for( int i = begin; i < end; ++i )
          for( int j = 0; j < _mesh->faces()[i].size(); ++j )
            shared_node_set.insert( (*_node_gid)[ _mesh->faces()[i][j] ] );
      }
    }
    _shared_node_list.resize( shared_node_set.size() );
    std::copy( shared_node_set.begin(), shared_node_set.end(), _shared_node_list.begin() );
  }

  void form_shared_node_proc_list() {
    std::vector<int> globalSize( Pstream::nProcs() ),
                     displ( Pstream::nProcs() + 1 ),
                     globalSharedList;
    int localSize;
    localSize = _shared_node_list.size();
    MPI_Allgather( &localSize, 1, MPI_INT, &globalSize[0], 1, MPI_INT, MPI_COMM_WORLD );
    displ[0] = 0;
    for( int i = 0; i < Pstream::nProcs(); ++i )
      displ[i+1] = globalSize[i] + displ[i];
    MPI_Datatype gatherType;
    MPI_Type_contiguous( localSize, MPI_INT, &gatherType );
    MPI_Type_commit( &gatherType );

    globalSharedList.resize( displ[Pstream::nProcs()] );
    MPI_Allgatherv( &_shared_node_list[0], 1, gatherType,
                    &globalSharedList[0], &globalSize[0], &displ[0],
                    MPI_INT, MPI_COMM_WORLD );
    std::map< Foam::label, std::set<int> > shared_node_proc;
    for( int i = 0; i < Pstream::nProcs(); ++i )
      for( int j = displ[i]; j < displ[i+1]; ++j )
        shared_node_proc[ globalSharedList[j] ].insert(i);
    /// Clear memory
    globalSize.clear();
    displ.clear();
    globalSharedList.clear();
    for( auto i = _shared_node_list.begin(); i != _shared_node_list.end(); ++i )
      if( shared_node_proc.find( *i ) != shared_node_proc.end() )
         _shared_node_proc_list[*i] = shared_node_proc[*i];
  }

  void form_shared_node_map() {
    forAll( (*_node_gid), inode )
      _node_gid_map[ (*_node_gid)[inode] ] = inode;
    for( auto i  = _shared_node_proc_list.begin(); i != _shared_node_proc_list.end(); ++i )
      for( auto j = i->second.begin(); j != i->second.end(); ++j )
        if( *j != Pstream::myProcNo() )
          _proc_shared_node_map[*j].insert(i->first);
    _proc_xadj.push_back(0);
    for( auto i  = _proc_shared_node_map.begin();
              i != _proc_shared_node_map.end(); ++i ) {
      _proc_list.push_back( i->first );
      _proc_xadj.push_back( i->second.size() );
      _proc_adjncy.insert( _proc_adjncy.end(), i->second.begin(), i->second.end() );
    }
    for( int i = 0; i < _proc_xadj.size() - 1; ++i )
      _proc_xadj[i + 1] += _proc_xadj[i];
    /// FORTRAN +1 indexing
    for( auto i = _proc_xadj.begin(); i != _proc_xadj.end(); ++i ) (*i)++;
    for( auto i = _proc_adjncy.begin(); i != _proc_adjncy.end(); ++i )
      (*i) = _node_gid_map[ (*i) ] + 1;
  }

  /////// Destructor
  ~__c_ofreader(){
    close_mesh();
    close_args();
    close_time();
  }
  /////// Access functions
  Foam::fvMesh *mesh() { return _mesh; }
  std::vector<int> &proc_list() { return _proc_list; }
  std::vector<int> &proc_xadj() { return _proc_xadj; }
  std::vector<int> &proc_adjncy() { return _proc_adjncy; }
  std::vector<int> &proc_snlist() { return _shared_node_list; }

  /////// Destory functions
  void close_mesh() { if( _mesh != nullptr ) delete _mesh; _mesh = nullptr; }
  void close_args() { if( _args != nullptr ) delete _args; _args = nullptr; }
  void close_time() { if( _time != nullptr ) delete _time; _time = nullptr; }

  private:
  /////// Member data
  Foam::argList     *_args;
  Foam::Time        *_time;
  Foam::fvMesh      *_mesh;
  /// Shared node list
  labelIOList        *_node_gid;
  std::map< Foam::label, Foam::label> _node_gid_map;
  std::vector<int>   _shared_node_list;
  std::map< Foam::label, std::set<int> > _shared_node_proc_list;
  std::map< int, std::set<int> > _proc_shared_node_map;
  std::vector<int>   _proc_list, _proc_xadj, _proc_adjncy; 
};

/// The global variable storing the OpenFOAM C wrapper object
__c_ofreader *global_of_mesh;
 
#endif

