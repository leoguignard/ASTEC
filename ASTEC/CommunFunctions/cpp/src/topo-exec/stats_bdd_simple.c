#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */

#include <vt_bdd2.h>

/*-------------------- statistics --------------------*/
#define TESTS_MAX 260
int TESTS_nb;
static int TESTS_array[TESTS_MAX + 1];
/*----------------------------------------------------*/


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
    char *base = basename( argv[0] );
    char filename[80];
    int permutations[PERMUTATION_SIZE];
    int i, nb, n ,np = _TRUE_;;


    for ( i=0 ; i <= TESTS_MAX; i++ ) TESTS_array[i] = 0;
    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
    nb = 0;


    do {
      for ( n=0; (n <100000) && (np == _TRUE_); n ++ ) {
	
	TESTS_nb = 0;
	(void)IsSimple( permutations);
	if ( (TESTS_nb > TESTS_MAX) ) {
	  fprintf(stderr, "TESTS_MAX is too small (should be at least %d)\n", TESTS_nb );
	  exit( 0 );
	} else {
	  TESTS_array[ TESTS_nb ] ++;
	}

	nb ++;
	np = Next_Permutation( permutations, PERMUTATION_SIZE );
      }
      fprintf(stderr, " %08d / 67108864 tested configurations\r", nb );
    } while ( np == _TRUE_ );

    printf("----- RESULTS: stats of computation -----\n");
    printf("      # of  simple points with bdds\n");
    printf("      %8d tested configurations \n", nb);

    {
      int min, max;
      int n;
      long int vmax, sum;

      for ( i=0 ; (i <= TESTS_MAX) && (TESTS_array[i]==0); i++ );
      min = i;
      for ( i=TESTS_MAX ; (i >= 0) && (TESTS_array[i]==0); i-- );
      max = i;
      sum = n = 0;
      vmax = TESTS_array[0];
      for ( i=0 ; i <= TESTS_MAX; i++ ) {
	sum += i * TESTS_array[i];
	n +=  TESTS_array[i];
	if ( vmax < TESTS_array[i] ) vmax = TESTS_array[i];
      }
      printf("      #minimum of tests = %d\n", min );
      printf("      #maximum of tests = %d\n", max );
      printf("      maximum of #tests = %ld\n", vmax );
      printf("      #average of tests = %ld / %d = %g\n", sum, n, (double)sum / (double)n );
    }
      
    {
      FILE *f, *fopen();
      f = fopen( strcat( strcpy( filename, base ), ".TESTS.RESULTS" ), "w" );
      for ( i=0 ; i <= TESTS_MAX; i++ ) {
	fprintf(f, " %3d %8d\n", i, TESTS_array[i] );
      }
      fclose( f );
    }
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
