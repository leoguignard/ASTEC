/*************************************************************************
 * chunks.c -
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <chunks.h>
#ifdef _OPENMP
#include <omp.h>
#endif



static int _verbose_ = 1;
static int _debug_ = 0;

void setVerboseInChunks( int v )
{
  _verbose_ = v;
}

void incrementVerboseInChunks(  )
{
  _verbose_ ++;
}

void decrementVerboseInChunks(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}

void setDebugInChunks( int v )
{
  _debug_ = v;
}

void incrementDebugInChunks(  )
{
  _debug_ ++;
}

void decrementDebugInChunks(  )
{
  _debug_ --;
  if ( _debug_ < 0 ) _debug_ = 0;
}




#ifdef _OPENMP
static int _max_chunks_ = 100;
#else
static int _max_chunks_ = 1;
#endif

static int _min_elements_ = 100;

static schedulingType _scheduling_ = _DYNAMIC_SCHEDULING_;





void setMaxChunks( int c )
{
  if ( c >= 1 )
    _max_chunks_ = c;
}


int getMaxChunks( )
{
  return( _max_chunks_ );
}


void setOpenMPScheduling( schedulingType s )
{
  _scheduling_ = s;
}


schedulingType getOpenMPScheduling( )
{
  return( _scheduling_ );
}




/************************************************************
 *
 * chunks processing
 *
 ************************************************************/



int processChunks( _chunk_callfunction ftn, typeChunks *chunks, char *from )
{
  int n;

  if ( chunks->scheduling != _NO_PARALLELISM_ && chunks->n_allocated_chunks > 1 ) {

    switch( chunks->scheduling ) {

    default :
    case _DEFAULT_SCHEDULING_ :
      if ( _debug_ ) fprintf( stderr, "%s: default openmp scheduling\n", from );    
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( n=0; n<chunks->n_allocated_chunks; n++ ) {	
	if ( _debug_ >=2 ) {
	  fprintf( stderr, "%s: processing chunk #%d/%d", from, n+1, chunks->n_allocated_chunks );
#ifdef _OPENMP
	  fprintf( stderr, " attributed to thread #%d", omp_get_thread_num() );
#endif
	  fprintf( stderr, "\n" );
	}
	chunks->data[n].ret = (*ftn)( chunks->data[n].parameters, chunks->data[n].first, chunks->data[n].last );
      }
      break;

    case _DYNAMIC_SCHEDULING_ :
      if ( _debug_ ) fprintf( stderr, "%s: dynamic openmp scheduling\n", from );    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for ( n=0; n<chunks->n_allocated_chunks; n++ ) {	
	if ( _debug_ >=2 ) {
	  fprintf( stderr, "%s: processing chunk #%d/%d", from, n+1, chunks->n_allocated_chunks );
#ifdef _OPENMP
	  fprintf( stderr, " attributed to thread #%d", omp_get_thread_num() );
#endif
	  fprintf( stderr, "\n" );
	}
	chunks->data[n].ret = (*ftn)( chunks->data[n].parameters, chunks->data[n].first, chunks->data[n].last );
      }
      break;

    case _DYNAMIC_ONE_SCHEDULING_ :
      if ( _debug_ ) fprintf( stderr, "%s: (dynamic,1) openmp scheduling\n", from );    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
      for ( n=0; n<chunks->n_allocated_chunks; n++ ) {	
	if ( _debug_ >=2 ) {
	  fprintf( stderr, "%s: processing chunk #%d/%d", from, n+1, chunks->n_allocated_chunks );
#ifdef _OPENMP
	  fprintf( stderr, " attributed to thread #%d", omp_get_thread_num() );
#endif
	  fprintf( stderr, "\n" );
	}
	chunks->data[n].ret = (*ftn)( chunks->data[n].parameters, chunks->data[n].first, chunks->data[n].last );
      }
      break;

    case _GUIDED_SCHEDULING_ :
      if ( _debug_ ) fprintf( stderr, "%s: (guided) openmp scheduling\n", from );    
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
      for ( n=0; n<chunks->n_allocated_chunks; n++ ) {	
	if ( _debug_ >=2 ) {
	  fprintf( stderr, "%s: processing chunk #%d/%d", from, n+1, chunks->n_allocated_chunks );
#ifdef _OPENMP
	  fprintf( stderr, " attributed to thread #%d", omp_get_thread_num() );
#endif
	  fprintf( stderr, "\n" );
	}
	chunks->data[n].ret = (*ftn)( chunks->data[n].parameters, chunks->data[n].first, chunks->data[n].last );
      }
      break;

    case _STATIC_SCHEDULING_ :
      if ( _debug_ ) fprintf( stderr, "%s: (static) openmp scheduling\n", from );    
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for ( n=0; n<chunks->n_allocated_chunks; n++ ) {	
	if ( _debug_ >=2 ) {
	  fprintf( stderr, "%s: processing chunk #%d/%d", from, n+1, chunks->n_allocated_chunks );
#ifdef _OPENMP
	  fprintf( stderr, " attributed to thread #%d", omp_get_thread_num() );
#endif
	  fprintf( stderr, "\n" );
	}
	chunks->data[n].ret = (*ftn)( chunks->data[n].parameters, chunks->data[n].first, chunks->data[n].last );
      }
      break;
    }

  }
  else {
    
    if ( _debug_ ) 
      fprintf( stderr, "%s: sequential loop\n", from );    
    for ( n=0; n<chunks->n_allocated_chunks; n++ ) {	
      if ( _debug_ >= 2 ) {
	fprintf( stderr, "%s: processing chunk #%d/%d", from, n+1, chunks->n_allocated_chunks );
	fprintf( stderr, " in a sequential way" );
	fprintf( stderr, "\n" );
      }
      chunks->data[n].ret = (*ftn)( chunks->data[n].parameters, chunks->data[n].first, chunks->data[n].last );
    }

  }



  /* check the returned values in case of error
   */
  for ( n=0; n<chunks->n_allocated_chunks; n++ ) {
    if ( chunks->data[n].ret != 1 ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: error when computing chunk #%d\n", from, n );
      }
      return( -1 );
    }
  }

  return( 1 );
}



/************************************************************
 *
 * chunks construction (ad hoc function)
 *
 ************************************************************/

int buildChunks( typeChunks *chunks, size_t first, size_t last, char *from )
{
  char *proc = "buildChunks";
  size_t size = last-first+1;
  int n =  _max_chunks_;


  chunks->scheduling = _scheduling_;

  /* only one chunk
   */
  if ( n == 1 ) {
    if ( allocBuildOneChunk( chunks, first, last ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error when allocating one chunk\n", proc );
      return( -1 );
    }
    return( 1 );
  }


  /* parallelism is not possible
     => one chunk
  */
#ifdef _OPENMP
  if ( omp_get_max_threads() == 1 ) {
    if ( allocBuildOneChunk( chunks, first, last ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error when building one chunk\n", proc );
      return( -1 );
    }
    return( 1 );
  }
#endif

  

  /* what to do when there are too few elements ?
   */
  if ( size < n*_min_elements_ ) {
    n = sqrt( (double)size );
#ifdef _OPENMP
    if ( size < omp_get_max_threads() ) {
      n = size;
    }
    else {
      if ( n < omp_get_max_threads() ) 
	n = omp_get_max_threads();
    }
#endif    
  }

  if ( allocBuildEqualChunks( chunks, first, last, n ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when building %d chunks (call from %s)\n", proc, n, from );
    return( -1 );
  }
  
  if ( _debug_ ) {
    fprintf( stderr, "\n" );
    if ( from != (char*)NULL )
      fprintf( stderr, "%s: has computed %d chunks from '%s'\n", proc, chunks->n_allocated_chunks, from );
    else
      fprintf( stderr, "%s: has computed %d chunks from '%s'\n", proc, chunks->n_allocated_chunks, from );
  }
  
  return( 1 );
}





/************************************************************
 *
 * chunks construction (generic functions)
 *
 ************************************************************/


int allocBuildOneChunk( typeChunks *chunks, size_t first, size_t last )
{
  char *proc = "oneChunk";

  if ( allocChunks( chunks, 1 ) <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when allocating one chunk\n", proc );
    return( -1 );
  }
  chunks->data[0].first = first;
  chunks->data[0].last = last;
  return( 1 );
}



int buildEqualChunks( typeChunk *chunk, size_t first, size_t last, int n )
{
  char *proc = "buildEqualChunks";
  size_t totalSize = last+1-first;
  size_t f = first;
  size_t i, chunkSize;
  
  if ( n <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: negative number of chunks\n", proc );
    return( -1 );
  }
    
  chunkSize = totalSize / n;

  if ( chunkSize <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: negative or null chunk size (too many chunks ?)\n", proc );
    return( -1 );
  }

  for ( i=0; i<n; i++, f+=chunkSize ) {
    /* chunk of size chunkSize+1
     */
    if ( chunkSize * (n-i) < totalSize ) {
      totalSize -= chunkSize + 1;
      chunk[ i ].first = f;
      chunk[ i ].last  = f+(chunkSize+1)-1;
      f++;
    }
    /* chunk of size chunkSize
     */
    else {
      totalSize -= chunkSize;
      chunk[ i ].first = f;
      chunk[ i ].last  = f+chunkSize-1;
    }
  }

  return( 1 );
}



int allocBuildEqualChunks( typeChunks *chunks, size_t first, size_t last, int nchunks )
{
  char *proc = "allocBuildEqualChunks";
  size_t size = last - first + 1;
  int n = nchunks;
  
  if ( n <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: negative number of chunks\n", proc );
    return( -1 );
  }

  if ( size < n ) {
    n = size;
  }

  if ( allocChunks( chunks, n ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when allocating chunks\n", proc );
    return( -1 );
  }
  
  if ( buildEqualChunks( &(chunks->data[0]), first, last, n ) != 1 ) {
    freeChunks( chunks );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when calculating chunks\n", proc );
    return( -1 );
  }
  
  return( 1 );
}



int allocBuildInequalChunks( typeChunks *chunks, size_t first, size_t last, 
			     double fraction,   /* fraction of data to be put in one bucket */
			     int chunks_bucket, /* number of chunks for one bucket */
			     size_t minimal_chunk_size, 
			     int maximal_chunk_number )
{
  char *proc = "allocBuildInequalChunks";
  size_t totalSize = last+1-first;
  size_t f = first;
  size_t i, chunkSize;
  
  int maxchunks;
  int n, nchunks;
  
  if ( maximal_chunk_number <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: negative maximal number of chunks\n", proc );
    return( -1 );
  }
  
  if ( chunks_bucket <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: negative number of chunks in a bucket\n", proc );
    return( -1 );
  }
  
  if ( fraction < 0.0 || 1.0 <= fraction ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: fraction should be in ]0,1[\n", proc );
    return( -1 );
  }
  


  /* calculating the number of chunks
   */
  for ( nchunks=0, maxchunks=maximal_chunk_number, totalSize=last+1-first; totalSize > 0 && maxchunks > 0;  ) {

    /* put the remaining data in the remaining chunks
       -> equal repartition in maxchunks chunks
    */
    if ( maxchunks <= chunks_bucket ) {
      nchunks += maxchunks;
      totalSize = 0;
      if ( _debug_ ) fprintf( stderr, "(1) add %d chunks -> %d\n", maxchunks, nchunks );
      continue;
    }
    
    /* the remaining data fits in one bucket
       -> equal repartition in chunks_bucket chunks
    */
    if ( totalSize / chunks_bucket < minimal_chunk_size ) {
      nchunks += chunks_bucket;
      totalSize = 0;
      if ( _debug_ ) fprintf( stderr, "(2) add %d chunks -> %d\n", chunks_bucket, nchunks );
      continue;
    }
    
    /* the considered fraction of data in the bucket
       yield too small chunks
       -> equal repartition in min( totalSize / minimal_chunk_size, maxchunks)
    */
    if ( (fraction * totalSize) / chunks_bucket < minimal_chunk_size ) {
      n = totalSize / minimal_chunk_size;
      if ( n > maxchunks ) n = maxchunks;
      nchunks += n;
      totalSize = 0;
      if ( _debug_ ) fprintf( stderr, "(3) add %d chunks -> %d\n", n, nchunks );
      continue;
    }
    
    /* after filling the bucket, there will be too much
       remaining data
       -> equal repartition in  maxchunks
    */
    if ( (maxchunks <= 2*chunks_bucket) &&
	 (fraction * totalSize) / chunks_bucket < ((1-fraction) * totalSize) / (maxchunks-chunks_bucket) ) {
      nchunks += maxchunks;
      totalSize = 0;
      if ( _debug_ ) fprintf( stderr, "(4) add %d chunks -> %d\n", maxchunks, nchunks );
      continue;
    }

    /* generic case
     */
    chunkSize = (fraction * totalSize) / chunks_bucket;
    
    nchunks += chunks_bucket;
    maxchunks -= chunks_bucket;
    totalSize -= chunks_bucket * chunkSize;
    
    if ( _debug_ ) fprintf( stderr, "(5) add %d chunks -> %d\n", chunks_bucket, nchunks );
  }
  
  
  
  /* allocating the chunks
   */
  if ( allocChunks( chunks, nchunks ) <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when allocating chunks\n", proc );
    return( -1 );
  }




  /* filling chunks
   */
  for ( nchunks=0, maxchunks=maximal_chunk_number, totalSize=last+1-first, f=first; totalSize > 0 && maxchunks > 0;  ) {
    
    /* put the remaining data in the remaining chunks
       -> equal repartition in maxchunks chunks
    */
    if ( maxchunks <= chunks_bucket ) {
      if ( buildEqualChunks( &(chunks->data[nchunks]), f, last, maxchunks ) != 1 ) {
	freeChunks( chunks );
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when calculating chunks\n", proc );
	return( -1 );
      }
      totalSize = 0;
      continue;
    }
    
    /* the remaining data fits in one bucket
       -> equal repartition in chunks_bucket chunks
    */
    if ( totalSize / chunks_bucket < minimal_chunk_size ) {
      if ( buildEqualChunks( &(chunks->data[nchunks]), f, last, chunks_bucket ) != 1 ) {
	freeChunks( chunks );
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when calculating chunks\n", proc );
	return( -1 );
      }
      totalSize = 0;
      continue;
    }
    
    /* the considered fraction of data in the bucket
       yield too small chunks
       -> equal repartition in min( totalSize / minimal_chunk_size, maxchunks)
    */
    if ( (fraction * totalSize) / chunks_bucket < minimal_chunk_size ) {
      n = totalSize / minimal_chunk_size;
      if ( n > maxchunks ) n = maxchunks;
      if ( buildEqualChunks( &(chunks->data[nchunks]), f, last, n ) != 1 ) {
	freeChunks( chunks );
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when calculating chunks\n", proc );
	return( -1 );
      }
      totalSize = 0;
      continue;
    }
    
    /* after filling the bucket, there will be too much
       remaining data
       -> equal repartition in  maxchunks
    */
    if ( (maxchunks <= 2*chunks_bucket) &&
	 (fraction * totalSize) / chunks_bucket < ((1-fraction) * totalSize) / (maxchunks-chunks_bucket) ) {
      if ( buildEqualChunks( &(chunks->data[nchunks]), f, last, maxchunks ) != 1 ) {
	freeChunks( chunks );
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when calculating chunks\n", proc );
	return( -1 );
      }
      totalSize = 0;
      continue;
    }

    /* generic case
     */
    chunkSize = (fraction * totalSize) / chunks_bucket;
    for (i=0; i<chunks_bucket; i++, f+=chunkSize, totalSize-=chunkSize, nchunks++, maxchunks-- ) {
      chunks->data[ nchunks ].first = f;
      chunks->data[ nchunks ].last  = f+chunkSize-1;
    }
    
  }
  
  return( 1 );
}





/************************************************************
 *
 * management
 *
 ************************************************************/



void initChunk( typeChunk *chunk )
{
  chunk->first = 1;
  chunk->last = 0;
  chunk->parameters = (void*)NULL;
  chunk->ret = 0;
}



void initChunks( typeChunks *chunks )
{
  chunks->data = (typeChunk*)NULL;
  chunks->n_allocated_chunks = 0;
  chunks->scheduling = _DEFAULT_SCHEDULING_;
}



void freeChunks( typeChunks *chunks )
{
  if ( chunks->data != (typeChunk*)NULL ) free( chunks->data );
  initChunks( chunks );
}



int allocChunks( typeChunks *chunks, int n )
{
  char *proc = "allocChunks";
  int i;
  
  if ( n <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: negative or null number of chunks\n", proc );
    return( -1 );
  }
    
  chunks->data = (typeChunk*)malloc( n *sizeof(typeChunk) );
  if ( chunks->data == (typeChunk*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when allocating %d chunks\n", proc, n );
    return( -1 );
  }

  for ( i=0; i<n; i++ )
    initChunk( &(chunks->data[i] ) );

  chunks->n_allocated_chunks = n;

  return( 1 );
}



void printChunks( FILE *theFile, typeChunks *chunks, char *s )
{
  char *proc = "printChunks";
  int i;
  FILE *f = theFile;

  if ( f == NULL ) f = stdout;
  
  if ( s != (char *)NULL )
    fprintf( f, "%s: information on '%s'\n", proc, s );
  
  if ( chunks->n_allocated_chunks <= 0 || chunks->data == (typeChunk*)NULL ) {
    fprintf( f, "empty chunks\n" );
    return;
  }
    
  for ( i=0; i<chunks->n_allocated_chunks; i++ )
    fprintf( f, "#%3d [%12lu %12lu] = %12lu\n", 
	     i, chunks->data[i].first, chunks->data[i].last, chunks->data[i].last + 1 - chunks->data[i].first );
  fprintf( f, "\n" );
}



