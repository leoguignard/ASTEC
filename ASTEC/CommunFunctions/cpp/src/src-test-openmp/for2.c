#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static double _GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}

static double _GetTime() 
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}


int main (int argc, char const *argv[]){

  int n;
  int nthreads;

  double clock_init;
  double clock_exit;

  double time_init;
  double time_exit;

#ifdef _OPENMP

  fprintf( stderr, "number of processors available = %d\n", omp_get_num_procs() );
  fprintf( stderr, "number of threads in the current team = %d\n", omp_get_num_threads() );
  fprintf( stderr, "maximum number of threads that could be used = %d\n", omp_get_max_threads() );

  for ( nthreads=1; nthreads<=8; nthreads++ ) {
    printf( "nthreads = %d\n", nthreads);

    clock_init = _GetClock();
    time_init = _GetTime();

    /* omp_set_num_threads(nthreads);*/

#pragma omp parallel for schedule(static,1) num_threads(nthreads)
    for(n=0;n<8;n++) {
      sleep( 10 );
      printf("  Element %d traité par le thread %d \n",n,omp_get_thread_num());
    }
    
    clock_exit = _GetClock();
    time_exit = _GetTime();

    fprintf( stderr, "elapsed time = %f\n", clock_exit - clock_init );
    fprintf( stderr, "elapsed time = %f\n", time_exit - time_init );
    printf( "\n" );

  }
    
#else
  printf("openmp non present\n");
#endif

  return 0;

}
