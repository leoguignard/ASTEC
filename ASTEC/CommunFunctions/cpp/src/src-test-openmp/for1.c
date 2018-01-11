#include <stdio.h>

#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main (int argc, char const *argv[]){

  int n;

#ifdef _OPENMP

  #pragma omp parallel for
  for(n=0;n<8;n++){
    printf("Element %d traité par le thread %d \n",n,omp_get_thread_num());
  }

#else
  printf("openmp non present\n");
#endif

  return 0;

}
