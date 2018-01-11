#include <stdio.h>
#include <stdlib.h>

#include <vt_t06t26.h>
#include <vt_t06.h>
#include <vt_t06noopt.h>


#define _FALSE_ 0
#define _TRUE_  1

#define NEIGHBORHOOD_SIZE 27
#define PERMUTATION_SIZE 26
#define ARRAY_SIZE 7

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
    int neighbors[NEIGHBORHOOD_SIZE];
    register int i, nb, n ,np = _TRUE_;
    int results_array[ARRAY_SIZE];
    register int t_06_1, t_06_2;
    int t_06_3, t_26;


    /*--- initialization ---*/
    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
    neighbors[13] = 1;
    for ( i=0 ; i < 9; i++ )  results_array[i] = 0;
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
	
	t_06_1 = Compute_NotOptT06( neighbors );
	t_06_2 = Compute_T06( neighbors );
	Compute_T06_and_T26( neighbors, &t_06_3, &t_26 );
	
	if ( (t_06_1 != t_06_2) || (t_06_1 != t_06_3) ) {
	  fprintf(stderr," Erreur : ");
	  fprintf(stderr," t_06_1 = %d -", t_06_1 );
	  fprintf(stderr," t_06_2 = %d -", t_06_2 );
	  fprintf(stderr," t_06_3 = %d \n", t_06_3 );
	  Print_Neighborhood( neighbors );
	  exit( 0 );
	}
	
	results_array[ t_06_1 ] ++;
	
	nb ++;
	np = Next_Permutation( permutations, PERMUTATION_SIZE );
      }
      fprintf(stderr, " %08d / 67108864 tested configurations\r", nb );
    } while ( np == _TRUE_ );

    printf( "----- RESULTS : # of  6-connected components -----\n");
    printf( " %8d tested configurations \n", nb);
    for ( i=0 ; i < ARRAY_SIZE; i++ )
	printf( " %8d configurations have %d  6-connected components\n", results_array[i], i );

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
