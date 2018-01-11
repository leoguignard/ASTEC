#include <vt_t06t26.h>
#include <stdio.h>

#define _FALSE_ 0
#define _TRUE_  1


#define NEIGHBORHOOD_SIZE 9
#define PERMUTATION_SIZE 8
/*---------------- Static functions ----------------*/
static int  Next_Permutation( int *tab, int l );
static void Print_Neighborhood( int neighbors[NEIGHBORHOOD_SIZE] );
static void Print_Neighborhood_And_Ts( int neighbors[NEIGHBORHOOD_SIZE], int t4, int t8 );



/*---
  on teste un 3x3x3 voisinage pour compter
  les occurrences des differentes configurations
  topologiques
---*/

void main( int argc, char *argv[] )
{
    int permutations[PERMUTATION_SIZE];
    int neighbors[NEIGHBORHOOD_SIZE];
    int i, j, nb;
    int results_array[9][7];
    int t_08, t_04;

    /*--- initialization ---*/
    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
    for ( i=0 ; i < 9; i++ )  
      for ( j=0 ; j < 7; j++ )  
	results_array[i][j] = 0;

    nb = 0;

    do {
        
        /*--- the neigborhood is filled.
	---*/
        for (i=0; i<4; i++)  neighbors[i] = permutations[i];
        neighbors[4] = 1;
        for (i=5; i<NEIGHBORHOOD_SIZE; i++) neighbors[i] = permutations[i-1];
	
	Compute_T04_and_T08( neighbors, &t_04, &t_08 );

	results_array[ t_08 ][ t_04 ] ++;

	Print_Neighborhood_And_Ts( neighbors, t_04, t_08 );

	nb ++;
    
} while ( Next_Permutation( permutations, PERMUTATION_SIZE ) == _TRUE_ );

    printf( "----- RESULTS -----\n");
    printf( " %8d tested configurations \n", nb);

    for ( j=0 ; j < 7; j++ )  
      for ( i=0 ; i < 9; i++ )
	if ( results_array[ i ][ j ] != 0 ) {
	  printf( " t08 = %2d - t04 = %2d - #configurations = %8d\n", i,j,results_array[ i ][ j ] );
	}

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

static void Print_Neighborhood_And_Ts( int neighbors[NEIGHBORHOOD_SIZE], int t4, int t8 )
{
    printf("%d %d %d\n", neighbors[ 0], neighbors[ 1], neighbors[ 2]);
    printf("%d %d %d  -  t4=%d t8=%d\n", neighbors[ 3], neighbors[4], neighbors[5],t4,t8);
    printf("%d %d %d\n",    neighbors[6], neighbors[7], neighbors[8]);
    printf("\n");
}
