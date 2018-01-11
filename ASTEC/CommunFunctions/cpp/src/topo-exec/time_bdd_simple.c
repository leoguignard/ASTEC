#include <vt_bdd.h>
#include <vt_bdd2.h>
#include <vt_bdd3.h>
#include <stdio.h>
#include <time.h>

#define _FALSE_ 0
#define _TRUE_  1

#define NEIGHBORHOOD_SIZE 27
#define PERMUTATION_SIZE 26
/*---------------- Static functions ----------------*/
static int  Next_Permutation( int *tab, int l );
static void Print_Neighborhood( int neighbors[NEIGHBORHOOD_SIZE] );

/*---
  on teste un 3x3x3 voisinage pour compter
  les occurrences des differentes configurations
  topologiques
---*/

int main( int argc, char *argv[] )
{
    int permutations[PERMUTATION_SIZE];
    int i, nb;

    long int time_0 = 0;
    long int time_1 = 0;
    long int time_2 = 0;
    long int time_3 = 0;
    long int time_4 = 0;
    double cps=CLOCKS_PER_SEC;
    double clocks_to_microsecs = 1000000.0 / cps;


    time_0 = clock();

    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
    nb = 0;
    do {
	nb ++;
    } while ( Next_Permutation( permutations, PERMUTATION_SIZE ) == _TRUE_ );

    time_1 = clock();

    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
    nb = 0;
    do {

      (void)_VT_IsSimple( permutations);

	nb ++;
    } while ( Next_Permutation( permutations, PERMUTATION_SIZE ) == _TRUE_ );

    time_2 = clock();

    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
    nb = 0;
    do {

      (void)IsSimple( permutations);

	nb ++;
    } while ( Next_Permutation( permutations, PERMUTATION_SIZE ) == _TRUE_ );

    time_3 = clock();

    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
    nb = 0;
    do {

      (void)NewIsSimple( permutations);

	nb ++;
    } while ( Next_Permutation( permutations, PERMUTATION_SIZE ) == _TRUE_ );

    time_4 = clock();

    printf("----- RESULTS : times for computation -----\n");
    printf("      # of simple points with bdds\n");
    printf("      %8d tested configurations \n", nb);
    printf(" CLOCKS_PER_SEC = %g \n", cps );
    
    printf(" t0 = %ld \n", time_0 );
    printf(" t1 = %ld (after generation of all neighborhoods)\n", time_1 );
    printf(" t2 = %ld (after 1st treatment of all neighborhoods)\n", time_2 );
    printf(" t3 = %ld (after 2nd treatment of all neighborhoods)\n", time_3 );
    printf(" t4 = %ld (after 3rd treatment of all neighborhoods)\n", time_4 );
    printf(" Total time for neighborhoods' management,    t1-t0 = %ld clocks\n\n", 
	   time_1 - time_0 );

    printf(" Fonction 1.\n");
    printf(" Total time for all treatements,              t2-t1 = %ld clocks\n", 
	   time_2 - time_1 );
    printf(" Total time for all function calls, (t2-t1)-(t1-t0) = %ld clocks\n", 
	   time_2 - time_1 - time_1 + time_0);
    printf(" Average time for all function calls\n");
    printf("                         ((t2-t1)-(t1-t0))/%d = %g clocks\n", 
	   nb, (double)(time_2 - time_1 - time_1 + time_0)/(double)nb );
      printf("                                                    = %g microsecs\n\n", 
	   (double)(time_2 - time_1 - time_1 + time_0)/(double)nb * clocks_to_microsecs );
  
    printf(" Fonction 2.\n");
    printf(" Total time for all treatements,              t3-t2 = %ld clocks\n", 
	   time_3 - time_2 );
    printf(" Total time for all function calls, (t3-t2)-(t1-t0) = %ld clocks\n", 
	   time_3 - time_2 - time_1 + time_0);
    printf(" Average time for all function calls\n");
    printf("                         ((t3-t2)-(t1-t0))/%d = %g clocks\n", 
	   nb, (double)(time_3 - time_2 - time_1 + time_0)/(double)nb );
    printf("                                                    = %g microsecs\n\n", 
	   (double)(time_3 - time_2 - time_1 + time_0)/(double)nb * clocks_to_microsecs );
    
    printf(" Fonction 3.\n");
    printf(" Total time for all treatements,              t4-t3 = %ld clocks\n", 
	   time_4 - time_3 );
    printf(" Total time for all function calls, (t4-t3)-(t1-t0) = %ld clocks\n", 
	   time_4 - time_3 - time_1 + time_0);
    printf(" Average time for all function calls\n");
    printf("                         ((t4-t3)-(t1-t0))/%d = %g clocks\n", 
	   nb, (double)(time_4 - time_3 - time_1 + time_0)/(double)nb );
    printf("                                                    = %g microsecs\n\n", 
	   (double)(time_4 - time_3 - time_1 + time_0)/(double)nb * clocks_to_microsecs );
    
    return( 0 ); 
}

static int Next_Permutation( int *tab, int l )
{
    int i, c;
    
    i = 0;
    c = _FALSE_;

    while ( (i < l) && (c == _FALSE_) ) {
        if ( tab[i] == 0 ) {
            tab[i] = 1;
            c = _TRUE_;
        } else {
            tab[i] = 0;
            i++;
        }
    }
    return( c );
}

static void Print_Neighborhood( int neighbors[NEIGHBORHOOD_SIZE] )
{
    printf("\n");
    printf("%d %d %d  -  ", neighbors[ 0], neighbors[ 1], neighbors[ 2]);
    printf("%d %d %d  -  ", neighbors[ 9], neighbors[10], neighbors[11]);
    printf("%d %d %d\n",    neighbors[18], neighbors[19], neighbors[20]);

    printf("%d %d %d  -  ", neighbors[ 3], neighbors[ 4], neighbors[ 5]);
    printf("%d %d %d  -  ", neighbors[12], neighbors[13], neighbors[14]);
    printf("%d %d %d\n",    neighbors[21], neighbors[22], neighbors[23]);

    printf("%d %d %d  -  ", neighbors[ 6], neighbors[ 7], neighbors[ 8]);
    printf("%d %d %d  -  ", neighbors[15], neighbors[16], neighbors[17]);
    printf("%d %d %d\n",    neighbors[24], neighbors[25], neighbors[26]);
    printf("\n");
}
