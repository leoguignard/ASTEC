#include <vt_greynumbers.h>
#include <stdio.h>

#define _FALSE_ 0
#define _TRUE_  1


#define NEIGHBORHOOD_SIZE 9
#define PERMUTATION_SIZE 8
#define PERMUTATION_SIZE_MINUS_ONE 7
/*---------------- Static functions ----------------*/
static int  Next_Permutation( int *tab, int l );
static void ImprimeResults( FILE *f );


/*---
  on teste un 3x3x3 voisinage pour compter
  les occurrences des differentes configurations
  topologiques
---*/

long int nb;
int perminit[PERMUTATION_SIZE] = {0,0,0,0,0,0,0,0};
int permutations[PERMUTATION_SIZE];
long int results_array[7][7][9][9];

void main( int argc, char *argv[] )
{
  int neighbors[NEIGHBORHOOD_SIZE];
  int np = _TRUE_, n, m, i, j, k, l;
  int tm, tmm, tp ,tpp;

  char *filename = "liste_allts_2D.RESULTS";
  FILE *f, *fopen();


  /*--- initialization ---*/
  for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = perminit[i];
  neighbors[4] = 1;

    for ( i=0 ; i < 7; i++ )  
    for ( j=0 ; j < 7; j++ )  
    for ( k=0 ; k < 9; k++ )  
    for ( l=0 ; l < 9; l++ )  
	results_array[i][j][k][l] = 0;
    
    nb = 0;

    do {
      for ( m = 0; (m < 1000) && (np == _TRUE_); m++ ) {
	for ( n = 0; (n < 1000000) && (np == _TRUE_); n++ ) {
	  
	  /*--- the neigborhood is filled.
	    ---*/
	  for (i=0; i<4; i++)  neighbors[i] = permutations[i];
	  for (i=5; i<NEIGHBORHOOD_SIZE; i++) neighbors[i] = permutations[i-1];
	
	  Compute_allTs_2D( neighbors, &tp, &tpp, &tm, &tmm );
	
	  results_array[ tmm ][ tm ][ tp ][ tpp ] ++; 
	  nb ++;

	  if ( neighbors[0] >= 1 ) {
	    /* on compte 0 ..... X */
	    results_array[ tmm ][ tm ][ tp ][ tpp ] ++; 
	    nb ++;
	    /* on fait X .... 1 */
	    neighbors[PERMUTATION_SIZE] = 1;
	    Compute_allTs_2D( neighbors, &tp, &tpp, &tm, &tmm );
	    results_array[ tmm ][ tm ][ tp ][ tpp ] ++; 
	    nb ++;
	    if ( neighbors[0] >= 2 ) {
	      /* si c'etait 2 .... 1 
		 on compte aussi 1 .... 2 */
	      results_array[ tmm ][ tm ][ tp ][ tpp ] ++; 
	      nb ++;
	      /* on fait 2 .... 2 */
	      neighbors[PERMUTATION_SIZE] = 2;
	      Compute_allTs_2D( neighbors, &tp, &tpp, &tm, &tmm );
	      results_array[ tmm ][ tm ][ tp ][ tpp ] ++; 
	      nb ++;
	    }
	    neighbors[PERMUTATION_SIZE] = 0;
	  }

	  np = Next_Permutation( permutations, PERMUTATION_SIZE_MINUS_ONE );
	}
	fprintf(stderr, " %04ld / 6561 tested configurations\r", nb );
      }

      f = fopen( filename, "w" );
      ImprimeResults( f );
      fclose( f );
      
    } while ( np == _TRUE_ );

    f = fopen( filename, "w" );
    fprintf( f, "----- RESULTS -----\n");
    fprintf( f, " %13ld tested configurations \n", nb);
    ImprimeResults( f );
    fclose( f );
}

static int Next_Permutation( int *tab, int l )
{
    int i, c;
    
    i = 0;
    c = _FALSE_;

    while ( (i < l) && (c == _FALSE_) ) {
        if ( tab[i] <= 1 ) {
            tab[i] ++;
            c = _TRUE_;
        } else {
            tab[i] = 0;
            i++;
        }
    }
    return( c );
}

static void ImprimeResults( FILE *f )
{
  int i,j,k,l,n;
  long int total;

  fprintf(f, " %04ld / 6561 tested configurations\r", nb );
  fprintf(f, " first permutation was:\n" );
  fprintf(f," {");
  for ( n = 0; n<PERMUTATION_SIZE; n++ ) {
    fprintf(f,"%d", perminit[n] );
    if ( n < PERMUTATION_SIZE-1 ) fprintf(f,",");
  }
  fprintf(f,"}\n");
  fprintf(f, " next permutation is:\n" );
  fprintf(f," {");
  for ( n = 0; n<PERMUTATION_SIZE; n++ ) {
    fprintf(f,"%d", permutations[n] );
    if ( n < PERMUTATION_SIZE-1 ) fprintf(f,",");
  }
  fprintf(f,"}\n");
  
  total = 0;
  for ( i=0 ; i < 7; i++ )  
    for ( j=0 ; j < 7; j++ )  
      for ( k=0 ; k < 9; k++ )  
	for ( l=0 ; l < 9; l++ )  
	  if ( results_array[ i ][ j ][ k ][ l ] != 0 ) { 
	    total += results_array[ i ][ j ][ k ][ l ];
	    fprintf( f, " T-- = %d  ", i );
	    fprintf( f, " T- = %d  ", j );
	    fprintf( f, " T+ = %d  ", k );
	    fprintf( f, " T++ = %d  ", l );
	    fprintf( f, " #configurations = %13ld ", results_array[ i ][ j ][ k ][ l ] );
	    
	    if ( i==0 ) {
	      fprintf( f, "minimal" );
	      if ( j == 0 ) fprintf( f, " -> well" );
	      if ( l==0 ) fprintf( f, "+maximal=interior");
	      if ( (l==1) && (j==1) ) fprintf( f, "+constructible");
	      if ( l>1 ) fprintf( f, "+convergent");
	    } else if ( l==0 ) {
	      fprintf( f, "maximal" );
	      if ( k == 0 ) fprintf( f, " -> peak" );
	      if ( (i==1) && (k==1) ) fprintf( f, "+destructible");
	      if ( i>1 ) fprintf( f, "+divergent");
	    } else if ( (i==1) && (k==1) ) {
	      fprintf( f, "destructible");
	      if ( (l==1) && (j==1) ) fprintf( f, "+constructible=simple side");
	      if ( l>1 ) fprintf( f, "+convergent");
	    } else if ( (l==1) && (j==1) ) {
	      fprintf( f, "constructible");
	      if ( i>1 ) fprintf( f, "+divergent");
	    } else if ( l>1 ) {
	      fprintf( f, "convergent");
	      if ( i>1 ) fprintf( f, "+divergent=saddle point");
	    } else if ( i>1 ) {
	      fprintf(f, "divergent");
	    }
	    fprintf( f, "\n" );
	  }
  
  fprintf( f,"                                                         -------------\n");
  fprintf( f,"                                                         %013ld\n", total);
  
}