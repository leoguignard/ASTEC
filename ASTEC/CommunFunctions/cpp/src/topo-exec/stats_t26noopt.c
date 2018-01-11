#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */

#include <vt_t26noopt.h>

/*-------------------- statistics --------------------*/
#define TESTS_MAX 260
int TESTS_nb, EQUIV_nb;
static int TESTS_array[TESTS_MAX + 1];
static int EQUIV_array[TESTS_MAX + 1];
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
    char *base = (char*)basename( argv[0] );
    char filename[80];
    int permutations[PERMUTATION_SIZE];
    int neighbors[NEIGHBORHOOD_SIZE];
    int i, nb, n ,np = _TRUE_;


    fprintf( stderr, " Histogram files are %s \n",
	     strcat( strcpy( filename, base ), ".TESTS.RESULTS" ) );
    fprintf( stderr, "                 and %s \n",
	     strcat( strcpy( filename, base ), ".EQUIV.RESULTS" ) );

    for ( i=0 ; i <= TESTS_MAX; i++ ) TESTS_array[i] = EQUIV_array[i] = 0;
    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
    nb = 0;


    do {
      for ( n=0; (n <100000) && (np == _TRUE_); n ++ ) {
	/*--- the neigborhood is filled.
	  
	  the neighborhood is numbered as follows 
	  
	  0  1  2  -   9 10 11  -  18 19 20
	  3  4  5  -  12 13 14  -  21 22 23
	  6  7  8  -  15 16 17  -  24 25 26
	  
	  the central point (#13) does not matter.
	  ---*/
        for (i=0; i<13; i++)  neighbors[i] = permutations[i];
        neighbors[13] = 1;
        for (i=14; i<NEIGHBORHOOD_SIZE; i++) neighbors[i] = permutations[i-1];
	
	TESTS_nb = EQUIV_nb = 0;
	(void)Compute_NotOptT26( neighbors );
	if ( (TESTS_nb > TESTS_MAX) || (EQUIV_nb > TESTS_MAX ) ) {
	  fprintf(stderr, "TESTS_MAX is too small (should be at least %d or %d)\n", TESTS_nb, EQUIV_nb );
	  exit( 0 );
	} else {
	  TESTS_array[ TESTS_nb ] ++;
	  EQUIV_array[ EQUIV_nb ] ++;
	}
	
	nb ++;
	np = Next_Permutation( permutations, PERMUTATION_SIZE );
      }
      fprintf(stderr, " %08d / 67108864 tested configurations\r", nb );
    } while ( np == _TRUE_ );

    printf("----- RESULTS: stats of computation -----\n");
    printf("      # of 26-components in the 26-neighborhood\n");
    printf("      %8d tested configurations \n", nb);

    {
      int min, max;
      int n;
      long int vmax, sum;

      for ( i=0 ; (i <= TESTS_MAX) && (TESTS_array[i]==0); i++ );
      min = i;
      for ( i=TESTS_MAX ; (i >= 0) && (TESTS_array[i]==0); i-- );
      max = i;
      sum = 0;
      vmax = TESTS_array[0];
      for ( i=0 ; i <= TESTS_MAX; i++ ) {
	sum += i * TESTS_array[i];
	if ( vmax < TESTS_array[i] ) vmax = TESTS_array[i];
      }
      printf("      #minimum of tests = %d\n", min );
      printf("      #maximum of tests = %d\n", max );
      printf("      maximum of #tests = %ld\n", vmax );
      printf("      #average of tests = %ld / %d = %g\n", sum, nb, (double)sum / (double)nb );

      for ( i=0 ; (i <= TESTS_MAX) && (EQUIV_array[i]==0); i++ );
      min = i;
      for ( i=TESTS_MAX ; (i >= 0) && (EQUIV_array[i]==0); i-- );
      max = i;
      sum = n = 0;
      vmax = EQUIV_array[1];
      for ( i=1 ; i <= TESTS_MAX; i++ ) {
	sum += i * EQUIV_array[i];
	n += EQUIV_array[i];
	if ( vmax < EQUIV_array[i] ) vmax = EQUIV_array[i];
      }
      printf("      #minimum of equivalences = %d\n", min );
      printf("      #maximum of equivalences = %d\n", max );
      printf("      maximum of #equivalences = %ld\n", vmax );
      printf("      #average of equivalences = %ld / %d = %g\n", sum, nb, (double)sum / (double)nb );
      printf("      #cases where equivalences are needed = %8d\n", n );
      printf("                     percentage / %8d = %6.2f\n", nb, (double)n*100.0/(double)nb );
    }
    {
      FILE *f, *fopen();
      f = fopen( strcat( strcpy( filename, base ), ".TESTS.RESULTS" ), "w" );
      for ( i=0 ; i <= TESTS_MAX; i++ ) {
	fprintf(f, " %3d %8d\n", i, TESTS_array[i] );
      }
      fclose( f );
    }
    {
      FILE *f, *fopen();
      f = fopen( strcat( strcpy( filename, base ), ".EQUIV.RESULTS" ), "w" );
      for ( i=1 ; i <= TESTS_MAX; i++ ) {
	fprintf(f, " %3d %8d\n", i, EQUIV_array[i] );
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
