#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* cf
   http://www.unixgarden.com/index.php/gnu-linux-magazine/decouverte-de-la-programmation-parallele-avec-openmp
   https://computing.llnl.gov/tutorials/openMP/
   http://openmp.org/wp/resources/
   compilation avec 
   gcc -fopenmp helloworld-openmp.c 
*/

 int main (int argc, char *argv[]) {
   int th_id, nthreads;

#ifdef _OPENMP

   /*Définition d'une région parallèle */
   #pragma omp parallel private(th_id)
   {
     /*Obtenir le numéro du thread en question*/
     th_id = omp_get_thread_num();
     printf("Hello World from thread %d\n", th_id);
     /*Synchroniser avec barrier */
     #pragma omp barrier
     /*Test si on est avec le thread principal (master thread)*/
     if ( th_id == 0 ) {
       nthreads = omp_get_num_threads();
       printf("There are %d threads\n",nthreads);
     }
   }


#else
  printf("openmp non present\n");
#endif


   return 0;
 }
