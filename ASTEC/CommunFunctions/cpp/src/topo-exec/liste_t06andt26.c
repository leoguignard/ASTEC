#include <vt_t06andt26.h>
#include <stdio.h>

#define _TEST_TIME_

#if defined(_TEST_TIME_)
#include <time.h>
#endif

#define _FALSE_ 0
#define _TRUE_  1

/*---------------- Static functions ----------------*/
static int  Next_Permutation( int *tab, int l );
static void Print_Neighborhood( int neighbors[27] );




/*---
  on teste un 3x3x3 voisinage pour compter
  les occurrences des differentes configurations
  topologiques
---*/

void main( int argc, char *argv[] )
{
    int permutations[26];
    int neighbors[27];
    int i, j, nb;
    int results_array[9][7];
    int t_26, t_06;
#if defined(_TEST_TIME_)
    long int t0, t1, t2;
#endif

    /*--- initialization ---*/
    for ( i=0 ; i < 9; i++ )  
      for ( j=0 ; j < 7; j++ )  
	results_array[i][j] = 0;


#if defined(_TEST_TIME_)
    t0 = clock ();
#endif
    for ( i=0 ; i < 26; i++ ) permutations[i] = 0;
    nb = 0;
    do {
        
        /*--- the neigborhood is filled.

              the neighborhood is numbered as follows 

	      0  1  2  -   9 10 11  -  18 19 20
	      3  4  5  -  12 13 14  -  21 22 23
	      6  7  8  -  15 16 17  -  24 25 26

	      the central point (#13) does not matter.
	---*/
        for (i=0; i<13; i++)  neighbors[i] = permutations[i];
        neighbors[13] = 1;
        for (i=14; i<27; i++) neighbors[i] = permutations[i-1];
	
	Compute_T06_and_T26( neighbors, &t_06, &t_26 );


	nb ++;
#ifndef _TEST_TIME_
	results_array[ t_26 ][ t_06 ] ++;
	if ( nb % 100000 == 0 )
	  fprintf(stderr, " %8d tested configurations\r", nb );
#endif
    } while ( Next_Permutation( permutations, 26 ) == _TRUE_ );
    
#if defined(_TEST_TIME_)
    t1 = clock ();
    for ( i=0 ; i < 26; i++ ) permutations[i] = 0;
    nb = 0;
    do {
      for (i=0; i<13; i++)  neighbors[i] = permutations[i];
      neighbors[13] = 1;
      for (i=14; i<27; i++) neighbors[i] = permutations[i-1];
      nb ++;
    } while ( Next_Permutation( permutations, 26 ) == _TRUE_ );
    t2 = clock ();
#endif


    printf( "----- RESULTS -----\n");
    printf( " %8d tested configurations \n", nb);

#ifndef _TEST_TIME_
    for ( j=0 ; j < 7; j++ )  
      for ( i=0 ; i < 9; i++ )
	if ( results_array[ i ][ j ] != 0 ) {
	  printf( " t26 = %2d - t06 = %2d - #configurations = %8d\n", i,j,results_array[ i ][ j ] );
	}
#endif

#if defined(_TEST_TIME_)
    printf( " - times in microseconds \n");
    printf( " t0 = %ld  -  t1 = %ld  -  t2 = %ld\n", t0, t1, t2 );
    printf( " Total time = %ld  -  total management time = %ld\n", t1-t0, t2-t1 );
    printf( " Average total time = %g  - average routine time = %g\n",
	    (double)(t1-t0)/(double)nb, (double)(2*t1-t0-t2)/(double)nb );
    printf( " - times in seconds \n");
    printf( " Average total time = %g  - average routine time = %g\n",
	    (double)(t1-t0)/(double)nb/(double)CLOCKS_PER_SEC, 
	    (double)(2*t1-t0-t2)/(double)nb/(double)CLOCKS_PER_SEC );
#endif

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

static void Print_Neighborhood( int neighbors[27] )
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
