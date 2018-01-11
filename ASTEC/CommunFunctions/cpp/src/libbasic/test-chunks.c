/*************************************************************************
 * test-chunks.c -
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <chunks.h>

int main ( int argc, char *argv[] )
{
  int t, ntests = 1;
  int bucket, maxchunks;
  double fraction;
  size_t minsize = 1000;
  size_t first = 0;
  size_t last;
  
  typeChunks chunks;

  initChunks( &chunks );
  srandom( time( 0 ) );

  
  for ( t=0; t<ntests; t++ ) {
    last = random();

    fraction = 0.25 + 0.75 * ((double)random()/(double)RAND_MAX);
    bucket = 1 + 14 * ((double)random()/(double)RAND_MAX);
    minsize = last / 100 + (last/100)*2.0*((double)random()/(double)RAND_MAX);
    maxchunks = 5 + 25 * ((double)random()/(double)RAND_MAX);

    if ( 0 ) {
      fraction  = 0.462016;
      bucket    = 1;
      minsize   = 2171099;
      maxchunks = 8;
      last = 105426069;
    }
    if ( 0 ) {
      fraction  = 0.316970;
      bucket    = 12;
      minsize   = 15422445;
      maxchunks = 22;
      last = 1237880908;
    }


    fprintf( stdout, "\n" );
    fprintf( stdout, "========================================\n" );
    fprintf( stdout, "==== fraction  = %f\n", fraction );
    fprintf( stdout, "==== bucket    = %d\n", bucket );
    fprintf( stdout, "==== minsize   = %lu\n", minsize );
    fprintf( stdout, "==== maxchunks = %d\n", maxchunks );
    fprintf( stdout, "====        [%lu - %lu]\n", first, last );
    fprintf( stdout, "========================================\n" );

    fprintf( stdout, "\n" );
    fprintf( stdout, "---- Equal repartition\n" );

    if ( allocBuildEqualChunks( &chunks, first, last, maxchunks ) != 1 ) {
      fprintf( stderr, "error when calculating equal repartition\n" );
    }
    printChunks( stdout, &chunks, (char*)NULL );
    freeChunks( &chunks );

    fprintf( stdout, "\n" );
    fprintf( stdout, "---- Inequal repartition\n" );

    if ( allocBuildInequalChunks( &chunks, first, last, fraction, bucket, minsize, maxchunks ) != 1 ) {
      fprintf( stderr, "error when calculating equal repartition\n" );
    }
    printChunks( stdout, &chunks, (char*)NULL );
    freeChunks( &chunks );

    fprintf( stdout, "\n" );


  }


  exit( 0 );
}

