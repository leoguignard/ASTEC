/*************************************************************************
 * chunks.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Fri Nov 16 22:45:44 CET 2012
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */

#ifndef _CHUNKS_H_
#define _CHUNKS_H_

#ifdef __cplusplus
extern "C" {
#endif





#include <stdlib.h>
#include <stdio.h>

typedef enum {
  _NO_PARALLELISM_,
  _DEFAULT_SCHEDULING_,
  _DYNAMIC_SCHEDULING_,
  _DYNAMIC_ONE_SCHEDULING_,
  _GUIDED_SCHEDULING_,
  _STATIC_SCHEDULING_
} schedulingType;

typedef  int (*_chunk_callfunction)(void *parameter, size_t first, size_t last);

/* defines 
   - one interval [first,last],
   - parameters to be passed to the procedure
     of type _chunk_callfunction()
     Typically, it is a struct containing the parameters
   - the returned value
*/
typedef struct {
  size_t first;
  size_t last;
  void *parameters;
  int ret;
} typeChunk;

typedef struct {
  typeChunk *data;
  int n_allocated_chunks;
  schedulingType scheduling;
} typeChunks;



extern void setVerboseInChunks( int v );
extern void incrementVerboseInChunks(  );
extern void decrementVerboseInChunks(  );
extern void setDebugInChunks( int v );
extern void incrementDebugInChunks(  );
extern void decrementDebugInChunks(  );

extern void setMaxChunks( int c );
extern int getMaxChunks( );
extern void setOpenMPScheduling( schedulingType s );
extern schedulingType getOpenMPScheduling( );



/************************************************************
 *
 * chunks processing
 *
 ************************************************************/

extern int processChunks( _chunk_callfunction ftn, typeChunks *chunks, char *from );



/************************************************************
 *
 * chunks construction
 *
 ************************************************************/

extern int buildChunks( typeChunks *chunks, size_t first, size_t last, char *from );

extern int allocBuildOneChunk( typeChunks *chunks, size_t first, size_t last );
extern int buildEqualChunks( typeChunk *chunk, size_t first, size_t last, int n );
extern int allocBuildEqualChunks( typeChunks *chunks, size_t first, size_t last, int n );
extern int allocBuildInequalChunks( typeChunks *chunks, size_t first, size_t last, 
			     double fraction,   /* fraction of data to be put in one bucket */
			     int chunks_bucket, /* number of chunks for one bucket */
			     size_t minimal_chunk_size, 
			     int maximal_chunk_number );



/************************************************************
 *
 * management
 *
 ************************************************************/

extern void initChunks( typeChunks *chunks );
extern void freeChunks( typeChunks *chunks );
extern int allocChunks( typeChunks *chunks, int n );
extern void printChunks( FILE *theFile, typeChunks *chunks, char *s );





#ifdef __cplusplus
}
#endif

#endif
