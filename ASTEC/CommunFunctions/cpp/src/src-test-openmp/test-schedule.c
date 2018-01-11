#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <time.h>
#include <sys/time.h>

static double _GetTime();
static double _GetClock();


int main (int argc, char const *argv[])
{
  int t, ntest = 10;
  int i, nelts = 100000000;

  double *elts = NULL;
  double sum;

  double sumtimeuser;
  double sumtimeproc;
  double sumratio;

  double timeproc;
  double timeuser;
  double ratio;


  int chunksize = 0;

#ifdef _OPENMP

  elts = (double*)malloc( nelts*sizeof(double) );
  if ( elts == NULL ) {
    fprintf( stderr, "pb allocation\n" );
    exit( 0 );
  }
  srandom( time(0) );

  for ( i=0; i<nelts; i++ )
    elts[i] = (double)random()/(double)(RAND_MAX);

#define _BEG { \
    timeuser = (- _GetTime()); \
    timeproc = (- _GetClock()); \
    sum = 0.0; \
}

#define _END { \
    timeproc += _GetClock(); \
    timeuser += _GetTime(); \
    ratio = timeproc / timeuser; \
    \
    fprintf( stderr, "test #%2d: user = %7.3f, proc = %7.3f, ratio = %7.5f --- sum = %g, average = %g\n", \
	     t, timeuser, timeproc, ratio, sum, sum/(double)nelts );	\
    \
    sumtimeuser += timeuser; \
    sumtimeproc += timeproc; \
    sumratio += ratio; \
}

  sumtimeuser = sumtimeproc = sumratio = 0.0;
  for ( t=0; t<ntest; t ++ ) {
    _BEG
    for ( i=0; i<nelts; i++ ) {
      sum += elts[i];
    }
    _END
  }
  fprintf( stderr, "sequential\n" );
  fprintf( stderr, "average time user = %7.3f\n", sumtimeuser / (double)ntest );
  fprintf( stderr, "average time proc = %7.3f\n", sumtimeproc / (double)ntest );
  fprintf( stderr, "average ratio     = %7.5f\n", sumratio / (double)ntest );
  fprintf( stderr, "\n" );




  sumtimeuser = sumtimeproc = sumratio = 0.0;
  for ( t=0; t<ntest; t ++ ) {
    _BEG
#pragma omp parallel for reduction(+:sum)
    for ( i=0; i<nelts; i++ ) {
      sum += elts[i];
    }
    _END
  }
  fprintf( stderr, "parallel\n" );
  fprintf( stderr, "average time user = %7.3f\n", sumtimeuser / (double)ntest );
  fprintf( stderr, "average time proc = %7.3f\n", sumtimeproc / (double)ntest );
  fprintf( stderr, "average ratio     = %7.5f\n", sumratio / (double)ntest );
  fprintf( stderr, "\n" );



  sumtimeuser = sumtimeproc = sumratio = 0.0;
  for ( t=0; t<ntest; t ++ ) {
    _BEG
#pragma omp parallel for schedule(static) reduction(+:sum)
    for ( i=0; i<nelts; i++ ) {
      sum += elts[i];
    }
    _END
  }
  fprintf( stderr, "schedule(static)\n" );
  fprintf( stderr, "average time user = %7.3f\n", sumtimeuser / (double)ntest );
  fprintf( stderr, "average time proc = %7.3f\n", sumtimeproc / (double)ntest );
  fprintf( stderr, "average ratio     = %7.5f\n", sumratio / (double)ntest );
  fprintf( stderr, "\n" );



  sumtimeuser = sumtimeproc = sumratio = 0.0;
  for ( t=0; t<ntest; t ++ ) {
    _BEG
#pragma omp parallel for schedule(dynamic) reduction(+:sum)
    for ( i=0; i<nelts; i++ ) {
      sum += elts[i];
    }
    _END
  }
  fprintf( stderr, "schedule(dynamic)\n" );
  fprintf( stderr, "average time user = %7.3f\n", sumtimeuser / (double)ntest );
  fprintf( stderr, "average time proc = %7.3f\n", sumtimeproc / (double)ntest );
  fprintf( stderr, "average ratio     = %7.5f\n", sumratio / (double)ntest );
  fprintf( stderr, "\n" );



 sumtimeuser = sumtimeproc = sumratio = 0.0;
  for ( t=0; t<ntest; t ++ ) {
    _BEG
#pragma omp parallel for schedule(guided) reduction(+:sum)
    for ( i=0; i<nelts; i++ ) {
      sum += elts[i];
    }
    _END
  }
  fprintf( stderr, "schedule(guided)\n" );
  fprintf( stderr, "average time user = %7.3f\n", sumtimeuser / (double)ntest );
  fprintf( stderr, "average time proc = %7.3f\n", sumtimeproc / (double)ntest );
  fprintf( stderr, "average ratio     = %7.5f\n", sumratio / (double)ntest );
  fprintf( stderr, "\n" );


  



  for ( chunksize=nelts/10; chunksize>=1000; chunksize /= 10 ) {
    
    sumtimeuser = sumtimeproc = sumratio = 0.0;
    for ( t=0; t<ntest; t ++ ) {
      _BEG
#pragma omp parallel for  schedule(static,chunksize) reduction(+:sum)
      for ( i=0; i<nelts; i++ ) {
	sum += elts[i];
      }
      _END
    }
    
    fprintf( stderr, "schedule (static,%d)\n", chunksize );
    fprintf( stderr, "average time user = %7.3f\n", sumtimeuser / (double)ntest );
    fprintf( stderr, "average time proc = %7.3f\n", sumtimeproc / (double)ntest );
    fprintf( stderr, "average ratio     = %7.5f\n", sumratio / (double)ntest );
    fprintf( stderr, "\n" );
  }

#else
  printf("openmp non present\n");
#endif

  return( 0 );
}




static double _GetTime() 
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}

static double _GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}




