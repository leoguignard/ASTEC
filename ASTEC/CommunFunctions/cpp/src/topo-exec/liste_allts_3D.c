#include <vt_greynumbers.h>
#include <stdio.h>

#define _FALSE_ 0
#define _TRUE_  1

#define _STRLENGTH_ 256

#define NEIGHBORHOOD_SIZE 27
#define PERMUTATION_SIZE 26
#define PERMUTATION_SIZE_MINUS_ONE 25
#define DIMN 1000000
#define DIMM 300




long long int nb, nc;
int perminit[PERMUTATION_SIZE];
int permutations[PERMUTATION_SIZE];
typedef long long int tab_result[7][7][9][9];
tab_result results_array;


/*---------------- Static functions ----------------*/
static int  Next_Permutation( int *tab, int l );

static void ReadFile( char *name, 
	      long long int *nConfig,
	      int *perm_deb,
	      int *perm_fin,
	      tab_result results );

static void ImprimeResults( char *name, 
	      long long int nConfig,
	      int *perm_deb,
	      int *perm_fin,
	      tab_result results );


/*---
  on teste un 3x3x3 voisinage pour compter
  les occurrences des differentes configurations
  topologiques
---*/


void main( int argc, char *argv[] )
{
  int neighbors[NEIGHBORHOOD_SIZE];
  int np = _TRUE_, n, m, i, j, k, l;
  int tm, tmm, tp ,tpp;

  char *filename = "liste_allts_3D.RESULTS";




  /*--- initialization ---*/
  for ( i=0 ; i < PERMUTATION_SIZE; i++ ) permutations[i] = perminit[i] = 0;
  neighbors[13] = 1;

    for ( i=0 ; i < 7; i++ )  
    for ( j=0 ; j < 7; j++ )  
    for ( k=0 ; k < 9; k++ )  
    for ( l=0 ; l < 9; l++ )  
	results_array[i][j][k][l] = 0;
    
    nb = nc = 0;

    /* lecture du fichier 
     */

    ReadFile( filename, &nb, perminit, permutations, results_array );
#if defined(_LINUX_)
	fprintf(stderr, " read %013Ld / 2541865828329 tested configurations", nb );
	fprintf(stderr, "\n" );
#else
	fprintf(stderr, " read %013ld / 2541865828329 tested configurations", nb );
	fprintf(stderr, "\n" );
#endif
	





    do {
      for ( m = 0; (m < DIMM) && (np == _TRUE_); m++ ) {
	for ( n = 0; (n < DIMN) && (np == _TRUE_); n++, nc++ ) {
	  
	  /*--- the neigborhood is filled.
	    ---*/
	  for (i=0; i<13; i++)  neighbors[i] = permutations[i];
	  for (i=14; i<NEIGHBORHOOD_SIZE; i++) neighbors[i] = permutations[i-1];
	
	  Compute_allTs_3D( neighbors, &tp, &tpp, &tm, &tmm );
	
	  results_array[ tmm ][ tm ][ tp ][ tpp ] ++; 
	  nb ++;

	  if ( neighbors[0] >= 1 ) {
	    /* on compte 0 ..... X */
	    results_array[ tmm ][ tm ][ tp ][ tpp ] ++; 
	    nb ++;
	    /* on fait X .... 1 */
	    neighbors[PERMUTATION_SIZE] = 1;
	    Compute_allTs_3D( neighbors, &tp, &tpp, &tm, &tmm );
	    results_array[ tmm ][ tm ][ tp ][ tpp ] ++; 
	    nb ++;
	    if ( neighbors[0] >= 2 ) {
	      /* si c'etait 2 .... 1 
		 on compte aussi 1 .... 2 */
	      results_array[ tmm ][ tm ][ tp ][ tpp ] ++; 
	      nb ++;
	      /* on fait 2 .... 2 */
	      neighbors[PERMUTATION_SIZE] = 2;
	      Compute_allTs_3D( neighbors, &tp, &tpp, &tm, &tmm );
	      results_array[ tmm ][ tm ][ tp ][ tpp ] ++; 
	      nb ++;
	    }
	    neighbors[PERMUTATION_SIZE] = 0;
	  }

	  np = Next_Permutation( permutations, PERMUTATION_SIZE_MINUS_ONE );
	}
#if defined(_LINUX_)
	fprintf(stderr, " %013Ld / 2541865828329 tested configurations", nb );
	fprintf(stderr, " (%012Ld/847288609443 perms) --- next writing %3d/%3d", nc, m, DIMM );
	fprintf(stderr, "\r" );
#else
	fprintf(stderr, " %013ld / 2541865828329 tested configurations", nb );
	fprintf(stderr, " (%012ld/847288609443 perms) --- next writing %3d/%3d", nc, m, DIMM );
	fprintf(stderr, "\r" );
#endif
      }


      ImprimeResults( filename, nb, perminit, permutations, results_array );
      
    } while ( np == _TRUE_ );



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


































static void ReadFile( char *name, 
	      long long int *nConfig,
	      int *perm_deb,
	      int *perm_fin,
	      tab_result results )
{
  FILE *f, *fopen();
  char str[_STRLENGTH_];
  char *tmp;
  int i, j, k, l, s;
  long long int n;
  
  f = fopen( name, "r" );

  fgets( str, _STRLENGTH_, f );
  sscanf( str, " %013Ld / 2541865828329 tested configurations", nConfig );

  fgets( str, _STRLENGTH_, f );
  fgets( str, _STRLENGTH_, f );
  tmp = str; tmp += 2;
  i = 0;
  while( i < PERMUTATION_SIZE ) {
     sscanf( tmp,"%d", &(perm_deb[i]) );
     tmp += 2;
     i++;
  }

  fgets( str, _STRLENGTH_, f );
  fgets( str, _STRLENGTH_, f );
  tmp = str; tmp += 2;
  i = 0;
  while( i < PERMUTATION_SIZE ) {
     sscanf( tmp,"%d", &(perm_fin[i]) );
     tmp += 2;
     i++;
  }

  while ( fgets( str, _STRLENGTH_, f ) ) {
    s = sscanf( str, " T-- = %d   T- = %d   T+ = %d   T++ = %d   #configurations = %13ld ", &i, &j, &k, &l, &n );
    if ( s == 5 ) {
      results[ i ][ j ][ k ][ l ] = n;
    }
  }

  fclose( f );
}



















static void ImprimeResults( char *name, 
	      long long int nConfig,
	      int *perm_deb,
	      int *perm_fin,
	      tab_result results )
{
  int i,j,k,l,n;
  long long int total;
  FILE *f, *fopen();

  if ( name == (char*)NULL ) f = stdout;
  else f = fopen( name, "w" );


#if defined(_LINUX_)
  fprintf(f, " %013Ld / 2541865828329 tested configurations\n", nConfig );
#else
  fprintf(f, " %013ld / 2541865828329 tested configurations\n", nConfig );
#endif
  fprintf(f, " first permutation was:\n" );
  fprintf(f," {");
  for ( n = 0; n<PERMUTATION_SIZE; n++ ) {
    fprintf(f,"%d", perm_deb[n] );
    if ( n < PERMUTATION_SIZE-1 ) fprintf(f,",");
  }
  fprintf(f,"}\n");
  fprintf(f, " next permutation is:\n" );
  fprintf(f," {");
  for ( n = 0; n<PERMUTATION_SIZE; n++ ) {
    fprintf(f,"%d", perm_fin[n] );
    if ( n < PERMUTATION_SIZE-1 ) fprintf(f,",");
  }
  fprintf(f,"}\n");
  
  total = 0;
  for ( i=0 ; i < 7; i++ )  
    for ( j=0 ; j < 7; j++ )  
      for ( k=0 ; k < 9; k++ )  
	for ( l=0 ; l < 9; l++ )  
	  if ( results[ i ][ j ][ k ][ l ] != 0 ) { 
	    total += results[ i ][ j ][ k ][ l ];
	    fprintf( f, " T-- = %d  ", i );
	    fprintf( f, " T- = %d  ", j );
	    fprintf( f, " T+ = %d  ", k );
	    fprintf( f, " T++ = %d  ", l );
#if defined(_LINUX_)
	    fprintf( f, " #configurations = %13Ld ", results[ i ][ j ][ k ][ l ] );
#else
	    fprintf( f, " #configurations = %13ld ", results[ i ][ j ][ k ][ l ] );
#endif

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
#if defined(_LINUX_)
  fprintf( f,"                                                         %013Ld\n", total);
#else
  fprintf( f,"                                                         %013ld\n", total);
#endif

  if ( name != (char*)NULL ) fclose( f );
 
}
