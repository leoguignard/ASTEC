#include <vt_t26.h>
#include <vt_t26noopt.h>
#include <stdio.h>
#include <time.h>      /* for clock() */

#define _FALSE_ 0
#define _TRUE_  1

#define NEIGHBORHOOD_SIZE 27
#define PERMUTATION_SIZE 26
/*---------------- Static functions ----------------*/
static int  Next_Permutation( int *tab, int l );

/*---
  on teste un 3x3x3 voisinage pour compter
  les occurrences des differentes configurations
  topologiques
---*/

int main( int argc, char *argv[] )
{
    int permutations[PERMUTATION_SIZE];
    int neighbors[NEIGHBORHOOD_SIZE];
    int i, nb;

    int clock_0 = 0;
    int clock_1 = 0;
    int clock_2 = 0;
    int clock_3 = 0;
    double cps=CLOCKS_PER_SEC;
    double clocks_to_microsecs = 1000000.0 / cps;


    clock_0 = clock();

    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
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
        for (i=14; i<NEIGHBORHOOD_SIZE; i++) neighbors[i] = permutations[i-1];
	nb ++;
    } while ( Next_Permutation( permutations, PERMUTATION_SIZE ) == _TRUE_ );

    clock_1 = clock();

    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
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
        for (i=14; i<NEIGHBORHOOD_SIZE; i++) neighbors[i] = permutations[i-1];

	(void)Compute_T26( neighbors );

	nb ++;
    } while ( Next_Permutation( permutations, PERMUTATION_SIZE ) == _TRUE_ );

    clock_2 = clock();

    for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = 0;
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
        for (i=14; i<NEIGHBORHOOD_SIZE; i++) neighbors[i] = permutations[i-1];

	(void)Compute_NotOptT26( neighbors );

	nb ++;
    } while ( Next_Permutation( permutations, PERMUTATION_SIZE ) == _TRUE_ );

    clock_3 = clock();

    printf("----- RESULTS : times for computation -----\n");
    printf("      # of 26-components in the 26-neighborhood\n");
    printf("      %8d tested configurations \n", nb);
    printf(" CLOCKS_PER_SEC = %g \n\n", cps );

    printf(" c0 = %12d \n", clock_0 );
    printf(" c1 = %12d (after generation of all neighborhoods)\n", clock_1 );
    printf(" c2 = %12d (after treatment of all neighborhoods)\n", clock_2 );
    printf(" c3 = %12d (after not opt. treatment of all neighborhoods)\n", clock_3 );
    printf(" Total time for neighborhoods' management,    c1-c0 = %d clock ticks\n", 
	   clock_1 - clock_0 );
    printf("                                                    = %g microsecs\n\n", 
	   (clock_1 - clock_0) * clocks_to_microsecs );
    
    printf(" OPT.\n");
    printf(" Total time for all treatements,              c2-c1 = %d clock ticks\n", 
	   clock_2 - clock_1 );
    printf(" Total time for all function calls, (c2-c1)-(c1-c0) = %d clock ticks\n", 
	   clock_2 - clock_1 - clock_1 + clock_0);
    printf(" Average time for all function calls\n");
    printf("                         ((c2-c1)-(c1-c0))/%d = %g clock ticks\n", 
	   nb, ((double)clock_2 - clock_1 - clock_1 + clock_0)/(double)nb );
    printf("                                                    = %g microsecs\n\n",
	   ((double)clock_2 - clock_1 - clock_1 + clock_0)/(double)nb * clocks_to_microsecs );

    
    printf(" NOT OPT.\n");
    printf(" Total time for all treatements,              c3-c2 = %d clock ticks\n", 
	   clock_3 - clock_2 );
    printf(" Total time for all function calls, (c3-c2)-(c1-c0) = %d clock ticks\n", 
	   clock_3 - clock_2 - clock_1 + clock_0);
    printf(" Average time for all function calls\n");
    printf("                         ((c3-c2)-(c1-c0))/%d = %g clock ticks\n", 
	   nb, (double)(clock_3 - clock_2 - clock_1 + clock_0)/(double)nb );
    printf("                                                    = %g microsecs\n\n",
	   ((double)clock_3 - clock_2 - clock_1 + clock_0)/(double)nb * clocks_to_microsecs );

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
