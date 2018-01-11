#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <YelTr3D.h>

#define STD_STR_LENGTH 256

static int _verbose_ = 1;
static int _debug_ = 1;





/* allow to sort in ascending order
 */

int compareDouble( const void *d1, const void *d2 ) {
  if ( *((double*)d1) >  *((double*)d2) ) {
    return( 1 );
  }
  else if  ( *((double*)d1) <  *((double*)d2) ) {
    return( -1 );
  }
  return( 0 );
}





void initTriangulation( typeTriangulation *t )
{
  t->nVertex = 0;
  t->vertex = (typeVertex*)NULL;
  t->nTriangle = 0;
  t->triangle = (typeTriangle*)NULL;
}



void freeTriangulation( typeTriangulation *t )
{
  if ( t->vertex != (typeVertex*)NULL )
    free( t->vertex );
  
  if (  t->triangle != (typeTriangle*)NULL )
    free( t->triangle );

  initTriangulation( t );
}



int readTriangulation( typeTriangulation *triangulation,
		       char *filename )
{
  char *proc = "readTriangulation";
  FILE *f;
  char *str = (char*)NULL;
  int i;
  int readingheader = 1;
  int waituntilend;
  int l=0;

  int nv, nt;
  typeVertex *v = (typeVertex*)NULL;
  typeTriangle *t = (typeTriangle*)NULL;



  f = fopen( filename, "r" );
  if ( f == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open file %s\n", proc, filename );
    return( -1 );
  }

  str = (char*)malloc( STD_STR_LENGTH * sizeof( char ) );
  if ( str == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate string\n", proc );
    fclose( f );
    return( -1 );
  }

  /* skip header (should be of 10 lines)
   */
  while ( readingheader ) {
    
    l++;
    if ( fgets( str, STD_STR_LENGTH, f ) == NULL ) {
	if ( _verbose_ ) {
	  fprintf( stderr, "%s: unexpected EOF when parsing %s\n", proc, filename );
	}
	free( str );
	fclose( f );
	return( -1 );
    }
  
    /* comment */
    if ( str[0] == '#' ) continue;
    
    /* skip containers CONTENT and TRIANGULATION3D
     */
    if ( strncmp( str, "CONTENT", 7 ) == 0 
	 || strncmp( str, "TRIANGULATION3D", 15 ) == 0 ) {
      waituntilend = 1;
      while ( waituntilend ) {

	l++;
	if ( fgets( str, STD_STR_LENGTH, f ) == NULL ) {
	  if ( _verbose_ ) {
	    fprintf( stderr, "%s: unexpected EOF when parsing %s\n", proc, filename );
	  }
	  free( str );
	  fclose( f );
	  return( -1 );
	}

	if ( strncmp( str, "END", 3 ) == 0 ) 
	  waituntilend = 0;
      }
      continue;
    }
      
    readingheader = 0;

  }
  


  /* header has be skiped
   */

  
  
  /* get number of vertices
   */
  if ( sscanf( str, "%d", &nv ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to read #vertices when parsing %s\n", proc, filename );
    }
    free( str );
    fclose( f );
    return( -1 );
  }
  


  if ( _verbose_ || _debug_ )
       fprintf( stderr, "%s: '%s' has %d vertices\n", proc, filename, nv );


  
  v = (typeVertex*)malloc( nv * sizeof(typeVertex) );
  if ( v == (typeVertex*)NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate vertex list (#vertices = %d)\n", proc, nv );
    }
    free( str );
    fclose( f );
    return( -1 );
  }
  
  for ( i=0; i<nv; i ++ ) {

    l++;
    if ( fgets( str, STD_STR_LENGTH, f ) == NULL ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: unexpected EOF when parsing %s\n", proc, filename );
      }
      free( v );
      free( str );
      fclose( f );
      return( -1 );
    }

    if ( sscanf( str, "%lf %lf %lf", &(v[i].x),  &(v[i].y), &(v[i].z) ) != 3 ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: unable to read vertex[%d] when parsing %s\n", proc, i, filename );
      }
      free( v );
      free( str );
      fclose( f );
      return( -1 );
    }
    
  }


  /* get number of triangles
   */
  
  l++;
  if ( fgets( str, STD_STR_LENGTH, f ) == NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unexpected EOF when parsing %s\n", proc, filename );
    }
      free( v );
      free( str );
      fclose( f );
      return( -1 );
  }

  if ( sscanf( str, "%d", &nt ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to read #vertices when parsing %s\n", proc, filename );
    }
    free( v );
    free( str );
    fclose( f );
    return( -1 );
  }



  if ( _verbose_ || _debug_ )
    fprintf( stderr, "%s: '%s' has %d triangles\n", proc, filename, nt );


  
  t = (typeTriangle*)malloc( nt * sizeof(typeTriangle) );
  if ( t == (typeTriangle*)NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate triangle list (#triangles = %d)\n", proc, nt );
    }
    free( v );
    free( str );
    fclose( f );
    return( -1 );
  }
  
  for ( i=0; i<nt; i ++ ) {

    l++;
    if ( fgets( str, STD_STR_LENGTH, f ) == NULL ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: unexpected EOF when parsing %s\n", proc, filename );
      }
      free( t );
      free( v );
      free( str );
      fclose( f );
      return( -1 );
    }

    /* fprintf( stderr, "line %d: %s\n", l, str ); */

    if ( sscanf( str, "%d %d %d %d %d %d", &(t[i].vertex[0]), &(t[i].vertex[1]), 
		 &(t[i].vertex[2]), &(t[i].neighbor[0]), &(t[i].neighbor[1]), 
		 &(t[i].neighbor[2])  ) != 6 ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: unable to read triangle[%d] when parsing %s\n", proc, i, filename );
      } 
      free( t );
      free( v );
      free( str );
      fclose( f );
      return( -1 );
    }
    
  }

  free( str );
  fclose( f );

  triangulation->nVertex = nv;
  triangulation->vertex = v;
  triangulation->nTriangle = nt;
  triangulation->triangle = t;

  return( 1 );
}



int printTriangulation( typeTriangulation *triangulation,
		       char *filename )
{
  char *proc = "printTriangulation";
  FILE *f;
  int i;

  f = fopen( filename, "w" );
  if ( f == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open file %s\n", proc, filename );
    return( -1 );
  }

  fprintf( f, "%d \n", triangulation->nVertex );
  
  for ( i=0; i<triangulation->nVertex; i++ ) {
    fprintf( f, "%13.10f %13.10f %13.10f \n", triangulation->vertex[i].x,
	     triangulation->vertex[i].y, triangulation->vertex[i].z );
  }
  
  fprintf( f, "%d \n", triangulation->nTriangle );

  for ( i=0; i<triangulation->nTriangle; i++ ) {
    fprintf( f, "%d %d %d %d %d %d \n", triangulation->triangle[i].vertex[0], 
	     triangulation->triangle[i].vertex[1], triangulation->triangle[i].vertex[2], 
	     triangulation->triangle[i].neighbor[0], triangulation->triangle[i].neighbor[1], 
	     triangulation->triangle[i].neighbor[2] );
  }

  fclose( f );

  return( 1 );
}



static double _smallestDistanceToVertex( typeVertex *v,
					 typeTriangulation *t )
{
  char *proc = "_smallestDistanceToVertex";
  double x, y, z;
  double min, d, dx, dy, dz;
  int i;

  if ( t->nVertex <= 0 || t->vertex == (typeVertex*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: empty triangulation ?!\n", proc );
    return( -1.0 );  
  }
  
  x = v->x;
  y = v->y;
  z = v->z;

  dx = t->vertex[0].x - x;
  dy = t->vertex[0].y - y;
  dz = t->vertex[0].z - z;

  min = dx*dx + dy*dy + dz*dz;
  
  for ( i=1; i<t->nVertex; i++ ) {
    dx = t->vertex[i].x - x;
    d = dx*dx;
    if ( d > min ) continue;
    dy = t->vertex[i].y - y;
    d += dy*dy;
    if ( d > min ) continue;
    dz = t->vertex[i].z - z;
    d += dz*dz;
    if ( d > min ) continue;
    min = d;
  }

  return( sqrt(min) );
}



int computeDistancesWithClosestVertex( typeTriangulation *t1, 
				       typeTriangulation *t2,
				       double *d )
{
  int i;

  for ( i=0; i<t1->nVertex; i++ ) {
    d[i] = _smallestDistanceToVertex( &(t1->vertex[i]), t2 );
  }

  return( 1 );

}



int computeDistancesWithIdenticalVertex( typeTriangulation *t1, 
					 typeTriangulation *t2,
					 double *d )
{
  char *proc = "computeDistancesWithIdenticalVertex";
  int i;
  double dx, dy, dz;

  if ( t1->nVertex != t2->nVertex ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: triangulations do not have the same number of vertices\n", proc );
    return( -1 );
  }

  for ( i=0; i<t1->nVertex; i++ ) {
    dx = t1->vertex[i].x -  t2->vertex[i].x;
    dy = t1->vertex[i].y -  t2->vertex[i].y;
    dz = t1->vertex[i].z -  t2->vertex[i].z;
    d[i] = sqrt( dx*dx + dy*dy + dz*dz );
  }

  return( 1 );

}


void printReport( double *d, int nd, char *t1, char *t2, char *desc )
{
  int i, middle1, middle2;
  double s, m, stddev;

  printf( "\n" );
  printf( "# %s\n", desc ); 
  printf( "# from vertices of triangulation '%s'\n", t1 );
  printf( "#   to vertices of triangulation '%s'\n", t2 );
 
  qsort( (void*)d, nd, sizeof(double), &compareDouble );

  s = 0;
  for ( i=0; i<nd; i++ )
    s += d[i];
  m = s/nd;

  s = 0;
  for ( i=0; i<nd; i++ )
    s += (d[i]-m)*(d[i]-m);
  stddev = sqrt( s/nd );

  middle1 = (nd-1)/2;
  middle2 = nd/2;

  if ( 0 ) printf( "middle %d %d\n", middle1, middle2 );

  printf( "minimum distance = %f\n", d[0] );
  if ( middle1 == middle2 )
    printf( "median  distance = %f\n", d[middle1] );
  else 
    printf( "median  distance = [ %f , %f ]\n", d[middle1], d[middle2] );
  printf( "average distance = %f - standard deviation = %f\n", m, stddev );
  printf( "maximum distance = %f\n", d[nd-1] );
  
  printf( "\n" );

}
