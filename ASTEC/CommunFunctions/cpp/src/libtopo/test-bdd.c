#include <stdio.h>
#include <string.h>

#include <vt_bdd.h>


/*-------------------- statistics --------------------*/
/*----------------------------------------------------*/


#define _FALSE_ 0
#define _TRUE_  1

#define NEIGHBORHOOD_SIZE 27
#define PERMUTATION_SIZE 26
/*---------------- Static functions ----------------*/
static int  Next_Permutation( int *tab, int l );

#ifdef _UNUSED_
static void Print_Neighborhood( int neighbors[NEIGHBORHOOD_SIZE] );
#endif

/*---
  on teste un 3x3x3 voisinage pour compter
  les occurrences des differentes configurations
  topologiques
---*/

int main( int argc, char *argv[] )
{
    int permutations[PERMUTATION_SIZE];
    int i, nb, n ,np = _TRUE_;

    int s=0;

    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
    nb = 0;


    do {
      for ( n=0; (n <100000) && (np == _TRUE_); n ++ ) {
	
	if ( VT_IsSimple( permutations) == 1 ) s++;

	nb ++;
	np = Next_Permutation( permutations, PERMUTATION_SIZE );
      }
      fprintf(stderr, " %08d / 67108864 tested configurations\r", nb );
    } while ( np == _TRUE_ );

    printf("----- RESULTS: stats of computation -----    \n");
    printf("      # of simple points with bdds = %d\n", s);
    printf("      %8d tested configurations \n", nb);

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

#ifdef _UNUSED_
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
#endif
