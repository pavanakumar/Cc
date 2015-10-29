/*! \file tecplotCommon.h
**
**/
#ifndef TECPLOT_COMMON_H

#define TECPLOT_COMMON_H

#include "TECIO.h"
#include "MASTER.h"

#include<cassert>

namespace tec360 {

  /// Determine the file type
  enum fileType { FULL=0, GRID=1, SOLUTION=2 };

  /// Cell types
  enum cellType { FEPOLYHEDRON=7 };

  /// Data packing
  enum dataPacking { POINT=0, BLOCK=1 };

  /// Value location
  enum valueLocation { CELL_CENTER=0, NODAL=1 };


} // End of tecplot namespace

#endif

