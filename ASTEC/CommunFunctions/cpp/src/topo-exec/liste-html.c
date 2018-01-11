#include <stdio.h>

#define _STRLENGTH_ 256
#define PERMUTATION_SIZE 26

int perminit[PERMUTATION_SIZE];
int permutations[PERMUTATION_SIZE];

/*
  long long int results_array[7][7][9][9];
*/
typedef long long int tab_result[7][7][9][9];
tab_result results_array;

long long int nb, nc;



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

static void ImprimeHtmlTable( char *name, 
	      long long int nConfig,
	      int *perm_deb,
	      int *perm_fin,
	      tab_result results );






void main( int argc, char *argv[] )
{
  int n, f, i, j, k, l;

  char name[ _STRLENGTH_ ];
  int p1[PERMUTATION_SIZE];
  int p2[PERMUTATION_SIZE];
  long long int conf;
  tab_result res;
  

  nb = 0;
  for ( n = 0; n<PERMUTATION_SIZE; n++ ) 
    perminit[n] = permutations[n] = 0;
  for ( i=0 ; i < 7; i++ )  
  for ( j=0 ; j < 7; j++ )  
  for ( k=0 ; k < 9; k++ )  
  for ( l=0 ; l < 9; l++ )  
    results_array[ i ][ j ][ k ][ l ] = 0;


  conf = 0;
  for ( n = 0; n<PERMUTATION_SIZE; n++ ) 
    p1[n] = p2[n] = 0;
  for ( i=0 ; i < 7; i++ )  
  for ( j=0 ; j < 7; j++ )  
  for ( k=0 ; k < 9; k++ )  
  for ( l=0 ; l < 9; l++ )  
    res[ i ][ j ][ k ][ l ] = 0;

  sprintf( name, "liste_allts_3D.RESULTS" );
  ReadFile( name, &conf, p1, p2, res );

  fprintf( stderr, " %s : %013ld configurations\n", name, conf );

  ImprimeHtmlTable( "3Dnumbers.html", conf, p1, p2, res );
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
  int verb = 1;
  
  f = fopen( name, "r" );

  fgets( str, _STRLENGTH_, f );
  if ( verb ) fprintf( stderr,"#1 : %s", str );

  sscanf( str, " %013Ld / 2541865828329 tested configurations", nConfig );

  fgets( str, _STRLENGTH_, f );
  if ( verb ) fprintf( stderr,"#2 : %s", str );
  fgets( str, _STRLENGTH_, f );
  if ( verb ) fprintf( stderr,"#3 : %s", str );

  tmp = str; tmp += 2;
  i = 0;
  while( i < PERMUTATION_SIZE ) {
     sscanf( tmp,"%d", &(perm_deb[i]) );
     tmp += 2;
     i++;
  }

  fgets( str, _STRLENGTH_, f );
  if ( verb ) fprintf( stderr,"#4 : %s", str );
  fgets( str, _STRLENGTH_, f );
  if ( verb ) fprintf( stderr,"#5 : %s", str );

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
    } else {
      if ( verb ) fprintf( stderr,"does not understand : %s", str );
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








static void ImprimeHtmlTable( char *name, 
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

  
  fprintf( f, "\n" );
  fprintf( f, "<center>\n" );
  fprintf( f, "<TABLE align=center border=1 cellspacing=1 cellpadding=1>\n" );
  fprintf( f, "<TR align=center>\n" );
  if ( 0 ) {
    fprintf( f, "<TH width=150> T04--\n" );
    fprintf( f, "<TH width=150> T04-\n" );
    fprintf( f, "<TH width=150> T08+\n" );
    fprintf( f, "<TH width=150> T08++\n" );
  } else {
    fprintf( f, "<TH width=150> T06--\n" );
    fprintf( f, "<TH width=150> T06-\n" );
    fprintf( f, "<TH width=150> T26+\n" );
    fprintf( f, "<TH width=150> T26++\n" );
  }
  fprintf( f, "<TH> number of configurations\n" );
  fprintf( f, "<TH> configuration name (if any)\n" );
  fprintf( f, "</TR>\n" );

  
  total = 0;
  for ( i=0 ; i < 7; i++ )  
    for ( j=0 ; j < 7; j++ )  
      for ( k=0 ; k < 9; k++ )  
	for ( l=0 ; l < 9; l++ )  
	  if ( results[ i ][ j ][ k ][ l ] != 0 ) { 
	    total += results[ i ][ j ][ k ][ l ];

	    fprintf( f, "<TR>" );
	    fprintf( f, " <TD align=center>%d", i );
	    fprintf( f, " <TD align=center>%d", j );
	    fprintf( f, " <TD align=center>%d", k );
	    fprintf( f, " <TD align=center>%d", l );
	    fprintf( f, "\n" );

	    fprintf( f, "<TD align=right>%ld\n", results[ i ][ j ][ k ][ l ] );
	    fprintf( f, "<TD align=left>" );

    
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
	    fprintf( f, "</TR>\n" );
	  }

  
  fprintf( f, "<TR>\n" );
  fprintf( f, "<TD align=right colspan=4><b>total (3<sup>26</sup>)</b>\n" );
  fprintf( f, "<TD align=right>%ld\n",total );
  fprintf( f, "</TR>\n" );
  fprintf( f, "<CAPTION align=bottom>\n" );
  fprintf( f, "<EM>Detailled 3D grey valued topological classification, with the number of\n" );
  fprintf( f, "configurations (neighborhood) for each set of values.\n" );
  fprintf( f, "</EM>\n" );
  fprintf( f, "</caption>\n" );
  fprintf( f, "</TABLE>\n" );
  fprintf( f, "</center>\n" );

  if ( name != (char*)NULL ) fclose( f );
 
}
