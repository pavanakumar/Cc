#ifndef __CC_OPENFOAM_READER__

#define __CC_OPENFOAM_READER__

#include "fvCFD.H"
#include "argList.H"

#include <set>
#include <map>
#include <sstream>
#include <fstream>

/**********************************************************
              Constants from FORTRAN code
**********************************************************/
const int _enable_parallel = 0;
const int _processor_bc = 0, _wall_bc = 1, _symmetry_bc = 2,
          _inlet_bc = 3, _outlet_bc = 4, _riemann_bc = 5;
/**********************************************************
   Function prototypes for the FORTRAN wrapper
**********************************************************/
extern "C" {
  void init_of_mesh( int *ipar );
  void close_of_mesh();
  void get_pm_sizes( int *nnode, int *nface, int *ninternalface, int *ncell, int *npatch );
  void get_pm_nodes( int *nNodes, double *x );
  void get_pm_faces( int *nface, int *ninternalface, int *facelr, int *facenodes );
  void get_pm_patches( int *npatch, int *patchdata );
  void check_metrics( int *ncell, int *nface, double *cv, double *cc, double *fc, double *fs, double *dn );
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
  }
  /////// Destructor
  ~__c_ofreader(){
    close_mesh();
    close_args();
    close_time();
  }
  /////// Access functions
  Foam::fvMesh *mesh() { return _mesh; }
  /////// Destory functions
  void close_mesh() { if( _mesh != nullptr ) delete _mesh; _mesh = nullptr; }
  void close_args() { if( _args != nullptr ) delete _args; _args = nullptr; }
  void close_time() { if( _time != nullptr ) delete _time; _time = nullptr; }

  private:
  /////// Member data
  Foam::argList     *_args;
  Foam::Time        *_time;
  Foam::fvMesh      *_mesh;
 
};

/// The global variable storing the OpenFOAM C wrapper object
__c_ofreader *global_of_mesh;
 
#endif

