#ifndef __WRAP_TECIO_CC_H__

#define __WRAP_TECIO_CC_H__

#include "tecplotGridFile.h"
#include <string>
#include <sstream>
#include <iostream>

extern "C" {
  void write_cc_tecio
  (
    int *ilvl, int *rank, int *nnode,
    int *ncell, int *nface,
    int *ninternalface, double *xyz,
    int *facelr, int *facenode
  );
}


#endif
