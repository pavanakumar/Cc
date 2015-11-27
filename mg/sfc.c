#include<stdio.h>
#include<limits.h>
#include<stdlib.h>

struct record_t {
  int index;
  unsigned key;
};

int comp_record_t
(
  const void *x,
  const void *y
)
{
  const struct record_t *a = (const struct record_t *)x;
  const struct record_t *b = (const struct record_t *)y;
  if( a->key < b->key ) return -1;
  else if( a->key == b->key ) return 0;
  else return 1;
}

unsigned morton10( unsigned x, unsigned y, unsigned z ) {
  /* Interleave x */
  x = (x | (x << 16)) & 0x030000FF;
  x = (x | (x <<  8)) & 0x0300F00F;
  x = (x | (x <<  4)) & 0x030C30C3;
  x = (x | (x <<  2)) & 0x09249249;
  /* Interleave y */
  y = (y | (y << 16)) & 0x030000FF;
  y = (y | (y <<  8)) & 0x0300F00F;
  y = (y | (y <<  4)) & 0x030C30C3;
  y = (y | (y <<  2)) & 0x09249249;
  /* Interleave z */
  z = (z | (z << 16)) & 0x030000FF;
  z = (z | (z <<  8)) & 0x0300F00F;
  z = (z | (z <<  4)) & 0x030C30C3;
  z = (z | (z <<  2)) & 0x09249249;
  /* return value */
  return x | (y << 1) | (z << 2); 
}

void sfc_perm
(
  int *n, double *x, double xmin[3],
  double xmax[3], int *perm, int *iperm
)
{
  unsigned i, j, intx[3], tenbit = (1<<10)-1;
  double eps[3] = { (double)tenbit, (double)tenbit, (double)tenbit };
  /* Create the records for sorting */
  struct record_t *records = (struct record_t *)malloc( sizeof(struct record_t) * (*n) );
  /* Precalculate the multiplication factor */
  for( i=0; i<3; ++i ) eps[i] /= ( xmin[i] - xmax[i] );
  /* Find the pixel xyz and interleave to get morton code */ 
  for( i=0; i<(*n); ++i) {
    for( j=0; j<3; ++j )
      intx[j] = (unsigned)( ( x[3 * i + j] - xmin[j] ) * eps[j] );
    records[i].index = i; /* FORTRAN indexing */
    records[i].key = morton10( intx[0], intx[1], intx[2] );
  }
  /* Sort the records */
  qsort( records, *n, sizeof(struct record_t), comp_record_t );
  /* Form the perm array */
  for( i=0; i<*n; ++i )
    perm[i] = records[i].index;
  /* Form the iperm array */
  for( i=0; i<*n; ++i )
    iperm[ perm[i] ] = i;
  /* FORTRAN indexing */
  for( i=0; i<*n; ++i ) {
    perm[i]++;
    iperm[i]++;
  }
  free(records);
}

