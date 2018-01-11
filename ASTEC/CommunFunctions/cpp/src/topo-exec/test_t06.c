#include <vt_t06.h>

#include <stdio.h>
#include <time.h>

#define _FALSE_ 0
#define _TRUE_  1

/*---------------- Static functions ----------------*/
static int  Next_Permutation( int *tab, int l );
static void Print_Neighborhood( int neighbors[27] );


void main( int argc, char *argv[] )
{
    int permutations[26];
    int neighbors[27];
    int i, nb;
    int results_array[7];
    int t_06_1, t_06_2, t_06_3, t_26;
#if defined(_TEST_TIME_)
    long int time_1 = 0;
    long int time_2 = 0;
    long int time_3 = 0;
#endif


    /*--- initialization ---*/
    for ( i=0 ; i < 26; i++ ) permutations[i] = 0;
    for ( i=0 ; i < 9; i++ )  results_array[i] = 0;
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
	

#if defined(_TEST_TIME_)
	time_1 -= clock();
#endif
	t_06_1 = VT_ComputeCbarre( neighbors );
#if defined(_TEST_TIME_)
	time_1 += clock();
#endif

#if defined(_TEST_TIME_)
	time_2 -= clock();
#endif
	t_06_2 = Compute_T06( neighbors );
#if defined(_TEST_TIME_)
	time_2 += clock();
#endif

#if defined(_TEST_TIME_)
	time_3 -= clock();
#endif
	Compute_T06_and_T26( neighbors, &t_06_3, &t_26 );
#if defined(_TEST_TIME_)
	time_3 += clock();
#endif
	
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
	if ( nb % 100000 == 0 )
	  fprintf(stderr, " %8d tested configurations\r", nb );
    } while ( Next_Permutation( permutations, 26 ) == _TRUE_ );

    printf( "----- RESULTS : # of  6-connected components -----\n");
    printf( " %8d tested configurations \n", nb);
    for ( i=0 ; i < 7; i++ )
	printf( " %8d configurations have %d  6-connected components\n", results_array[i], i );
    
#if defined(_TEST_TIME_)
    printf( "Total time for method #1 = %d ms, average time = %f \n",
	    time_1, (double)time_1/(double)nb );
    
    printf( "Total time for method #2 = %d ms, average time = %f \n",
	    time_2, (double)time_2/(double)nb );
    
    printf( "Total time for method #3 = %d ms, average time = %f \n",
	    time_3, (double)time_3/(double)nb );
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
