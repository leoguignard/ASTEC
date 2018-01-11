#include <vt_t06.h>
#include <vt_t06noopt.h>
#include <stdio.h>
#include <sys/times.h> /* for times() */
#include <time.h>      /* for clock() */

#define _FALSE_ 0
#define _TRUE_  1

#define NEIGHBORHOOD_SIZE 27
#define PERMUTATION_SIZE 26
/*---------------- Static functions ----------------*/
static int  Next_Permutation( int *tab, int l );
static void Print_Neighborhood( int neighbors[NEIGHBORHOOD_SIZE] );
static void PrintTimes( struct tms *ts );

/*---
  on teste un 3x3x3 voisinage pour compter
  les occurrences des differentes configurations
  topologiques
---*/

void main( int argc, char *argv[] )
{
    int permutations[PERMUTATION_SIZE];
    int neighbors[NEIGHBORHOOD_SIZE];
    int i, nb;

    int time_0 = 0, clock_0 = 0;
    int time_1 = 0, clock_1 = 0;
    int time_2 = 0, clock_2 = 0;
    int time_3 = 0, clock_3 = 0;
    struct tms tms_0, tms_1, tms_2, tms_3;
    double cps=CLOCKS_PER_SEC;
    double clocks_to_microsecs = 1000000.0 / cps;


    time_0 = times( &tms_0 );
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

    time_1 = times( &tms_1 );
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

	(void)Compute_T06( neighbors );

	nb ++;
    } while ( Next_Permutation( permutations, PERMUTATION_SIZE ) == _TRUE_ );

    time_2 = times( &tms_2 );
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

	(void)Compute_NotOptT06( neighbors );

	nb ++;
    } while ( Next_Permutation( permutations, PERMUTATION_SIZE ) == _TRUE_ );

    time_3 = times( &tms_3 );
    clock_3 = clock();

    printf("----- RESULTS : times for computation -----\n");
    printf("      # of  6-components in the 18-neighborhood\n");
    printf("      %8d tested configurations \n", nb);
    printf(" CLOCKS_PER_SEC = %g \n\n", cps );

    printf(" c0 = %12d, t0 = %12d \n", clock_0, time_0 );
    PrintTimes( &tms_0 );
    printf(" c1 = %12d, t1 = %12d (after generation of all neighborhoods)\n", clock_1, time_1 );
    PrintTimes( &tms_1 );
    printf(" c2 = %12d, t2 = %12d (after treatment of all neighborhoods)\n", clock_2, time_2 );
    PrintTimes( &tms_2 );
    printf(" c3 = %12d, t3 = %12d (after not opt. treatment of all neighborhoods)\n", clock_3, time_3 );
    PrintTimes( &tms_3 );
    printf(" Total time for neighborhoods' management,    t1-t0 = %d clock ticks\n", 
	   time_1 - time_0 );
    printf(" Total time for neighborhoods' management,    c1-c0 = %g microsecs\n\n", 
	   (clock_1 - clock_0) * clocks_to_microsecs );
    
    printf(" OPT.\n");
    printf(" Total time for all treatements,              t2-t1 = %d clock ticks\n", 
	   time_2 - time_1 );
    printf(" Total time for all function calls, (t2-t1)-(t1-t0) = %d clock ticks\n", 
	   time_2 - time_1 - time_1 + time_0);
    printf(" Average time for all function calls\n");
    printf("                         ((t2-t1)-(t1-t0))/%d = %g clock ticks\n", 
	   nb, ((double)time_2 - time_1 - time_1 + time_0)/(double)nb );
    printf("                         ((c2-c1)-(c1-c0))/%d = %g microsecs\n\n", nb,
	   ((double)clock_2 - clock_1 - clock_1 + clock_0)/(double)nb * clocks_to_microsecs );

    
    printf(" NOT OPT.\n");
    printf(" Total time for all treatements,              t3-t2 = %d clock ticks\n", 
	   time_3 - time_2 );
    printf(" Total time for all function calls, (t3-t2)-(t1-t0) = %d clock ticks\n", 
	   time_3 - time_2 - time_1 + time_0);
    printf(" Average time for all function calls\n");
    printf("                         ((t3-t2)-(t1-t0))/%d = %g clock ticks\n", 
	   nb, (double)(time_3 - time_2 - time_1 + time_0)/(double)nb );
    printf("                         ((c2-c1)-(c1-c0))/%d = %g microsecs\n\n", nb,
	   ((double)clock_3 - clock_2 - clock_1 + clock_0)/(double)nb * clocks_to_microsecs );
}

static void PrintTimes( struct tms *ts )
{
  printf(" [User time], [System time], [User time, children], [System time, children]\n" );
  printf("   %8d ,     %8d ,             %8d ,               %8d\n", 
	 (int)ts->tms_utime, (int)ts->tms_stime, (int)ts->tms_cutime, (int)ts->tms_cstime );
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
