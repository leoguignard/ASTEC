

#include <stdlib.h>
#include <stdio.h>
#include <time.h>


#include <vt_jacobien.h>

#define TINY 1e-12

static void PrintMat4x4( double *m )
{
  register int i;
  for (i=0;i<16;i++ ) {
    printf( "   %12g", m[i] );
    if ( (i+1)%4 == 0 ) printf( "\n");
  }
}

static void InverseMat4x4( double *matrice, double *inv )
{
  register int i, j, k;
  int kmax, rang = 4;
  register double c, max;
  double mat [16];
  
  for (i=0; i<16; i++ ) {
    mat[i] = matrice[i] ;
    inv[i] = 0.0;
  }
  inv[0] = inv[5] = inv[10] = inv[15] = 1.0;
  
  for ( j=0; j<4; j++ ) {
    if ( (mat[j*4+j] > (-TINY)) && (mat[j*4+j] < TINY) ) {
      /* recherche du plus grand element non nul sur la colonne j */
      kmax = j;
      max = 0.0;
      for (k=j+1; k<4; k++ ) {
	c = ( mat[k*4+j] > 0.0 ) ? mat[k*4+j] : (-mat[k*4+j]) ;
	if ( (c > TINY) && (c > max) ) { max = c; kmax = k; }
      }
      if ( kmax == j ) {
	/* la ligne est nulle */
	rang --;
      } else {
	/* sinon, on additionne */
	for ( i=0; i<4; i++ ) {
	  mat[j*4+i] += mat[kmax*4+i];
	  inv[j*4+i] += inv[kmax*4+i];
	}
      }
    }
    if ( (mat[j*4+j] < (-TINY)) || (mat[j*4+j] > TINY) ) {
      /* les autres lignes */
      for (k=0; k<4; k++) {
	if ( k != j ) {
	  c = mat[k*4 + j] / mat[j*4 + j];
	  for ( i=0; i<4; i++ ) {
	    mat[k*4 + i] -= c * mat[j*4 + i];
	    inv[k*4 + i] -= c * inv[j*4 + i];
	  }
	}
      }
      /* la ligne */
      c = mat[j*4 + j];
      for ( i=0; i<4; i++ ) {
	mat[j*4 + i] /= c;
	inv[j*4 + i] /= c;
      }
    } else {
      /* cas singulier */
      /*
      for (k=0; k<j; k++) {
	if ( (mat[k*4+k] < (-TINY)) || (mat[k*4+k] > TINY) ) {
	  c = mat[k*4+j] / mat[k*4+k];
	  for ( i=0; i<4; i++ ) {
	    mat[4*i+j] -= c * mat[4*i+k];
	    inv[4*i+j] -= c * inv[4*i+k];
	  }
	}
      }
      */
    }
    /*
    printf("\n--- MATRICE %d rang = %d ---\n",j,rang);
    PrintMat4x4( mat );
    printf("--- INVERSE %d ---\n",j);
    PrintMat4x4( inv );
    {
    double prd[16];
      Mat4x4ByMat4x4( matrice, inv, prd );
      printf("--- PRODUIT %d ---\n",j);
      PrintMat4x4( prd );
    }
    printf("\n");
    */
  }
  if ( rang < 4 ) {
    printf("--- MATRICE %d rang = %d ---\n",j,rang);
    PrintMat4x4( mat );
  }
}

static void Mat4x4ByMat4x4( double *m1, double *m2, double *prod )
{
  prod[0] = m1[0]*m2[0] + m1[1]*m2[4] + m1[2]*m2[8] + m1[3]*m2[12];
  prod[1] = m1[0]*m2[1] + m1[1]*m2[5] + m1[2]*m2[9] + m1[3]*m2[13];
  prod[2] = m1[0]*m2[2] + m1[1]*m2[6] + m1[2]*m2[10] + m1[3]*m2[14];
  prod[3] = m1[0]*m2[3] + m1[1]*m2[7] + m1[2]*m2[11] + m1[3]*m2[15];

  prod[4] = m1[4]*m2[0] + m1[5]*m2[4] + m1[6]*m2[8] + m1[7]*m2[12];
  prod[5] = m1[4]*m2[1] + m1[5]*m2[5] + m1[6]*m2[9] + m1[7]*m2[13];
  prod[6] = m1[4]*m2[2] + m1[5]*m2[6] + m1[6]*m2[10] + m1[7]*m2[14];
  prod[7] = m1[4]*m2[3] + m1[5]*m2[7] + m1[6]*m2[11] + m1[7]*m2[15];

  prod[8] = m1[8]*m2[0] + m1[9]*m2[4] + m1[10]*m2[8] + m1[11]*m2[12];
  prod[9] = m1[8]*m2[1] + m1[9]*m2[5] + m1[10]*m2[9] + m1[11]*m2[13];
  prod[10] = m1[8]*m2[2] + m1[9]*m2[6] + m1[10]*m2[10] + m1[11]*m2[14];
  prod[11] = m1[8]*m2[3] + m1[9]*m2[7] + m1[10]*m2[11] + m1[11]*m2[15];

  prod[12] = m1[12]*m2[0] + m1[13]*m2[4] + m1[14]*m2[8] + m1[15]*m2[12];
  prod[13] = m1[12]*m2[1] + m1[13]*m2[5] + m1[14]*m2[9] + m1[15]*m2[13];
  prod[14] = m1[12]*m2[2] + m1[13]*m2[6] + m1[14]*m2[10] + m1[15]*m2[14];
  prod[15] = m1[12]*m2[3] + m1[13]*m2[7] + m1[14]*m2[11] + m1[15]*m2[15];
}

/********************************************************************************
 *
 *
 *
 ********************************************************************************/

void JCB_InitRandom( int seed )
{
  if ( seed < 0 ) srandom( time(0) );
  else            srandom( seed );
}



/* return a value between 0 and 1
 */
double JCB_Random()
{
  double r=((double)random()) / ((double)RAND_MAX);
  return( r );
}

static void _initMatrix( double *mat )
{
  int i;
  for ( i=0; i<16; i++ ) mat[i] = 0;
  mat[ 0] = 1;
  mat[ 5] = 1;
  mat[10] = 1;
  mat[15] = 1;
}


void JCB_PrintMatrix( FILE *f, double *mat, char *str )
{
  int i, j;
  if ( str != NULL ) 
    fprintf( f, "%s\n", str );
  else
    fprintf( f, "\n" );
  for ( i=0; i<4; i++ ) {
    for ( j=0; j<4; j++ ) 
      fprintf( f, " %8.5g ", mat[i*4+j] );
    fprintf( f, "\n" );
  }
  fprintf( f, "\n" );
}



void JCB_RndTranslation( double *mat, double tmin, double tmax )
{
  int i;
  for ( i=0; i<3; i++ )
    mat[i*4+3] = tmin + (tmax-tmin)*JCB_Random();
}

void JCB_RndTranslationMatrix( double *mat, double tmin, double tmax )
{

  _initMatrix( mat );
  JCB_RndTranslation( mat, tmin, tmax );
}

void JCB_RndScaleMatrix( double *mat, double smin, double smax, double tmin, double tmax )
{

  _initMatrix( mat );
  JCB_RndTranslation( mat, tmin, tmax );
  mat[0] = mat[5] = mat[10] = smin + (smax-smin)*JCB_Random();
}

void JCB_RndSimilitudeMatrix( double *mat, double smin, double smax, double tmin, double tmax )
{
  int i;
  _initMatrix( mat );
  JCB_RndTranslation( mat, tmin, tmax );
  for ( i=0; i<3; i++ )
    mat[i*4+i] = smin + (smax-smin)*JCB_Random();
}

void JCB_RndAffineMatrix( double *mat, double amin, double amax, double tmin, double tmax )
{
  int i, j;
  _initMatrix( mat );
  JCB_RndTranslation( mat, tmin, tmax );
  for ( i=0; i<3; i++ )
    for ( j=0; j<3; j++ )     
      mat[i*4+j] = amin + (amax-amin)*JCB_Random();
}

/********************************************************************************
 *
 *
 *
 ********************************************************************************/

void JCB_initListPoint3D( typeListPoint3D *l )
{
  l->n = 0;
  l->pts = NULL;
}

int JCB_allocListPoint3D( typeListPoint3D *l, int n )
{
  typePoint3D *p;
 
  if ( n <= 0 ) return( -1 );
  p = (typePoint3D*)malloc( n*sizeof( typePoint3D ));
  if ( p == NULL ) return( -2 );
  l->n = n;
  l->pts = p;
  return( 1 );
}

void JCB_freeListPoint3D( typeListPoint3D *l )
{
  if ( l->pts != NULL ) free( l->pts );
  JCB_initListPoint3D( l );
}

void JCB_PrintListPoint3D( FILE *f, typeListPoint3D *l , char *str )
{
  int i;
  if ( str != NULL ) 
    fprintf( f, "%s\n", str );
  else
    fprintf( f, "\n" );

  for ( i=0; i<l->n; i++ ) {
    fprintf( f, "#%2d %8.5g  %8.5g  %8.5g\n", i, l->pts[i].x, l->pts[i].y, l->pts[i].z );
  }
  fprintf( f, "\n" );
}

int JCB_get6Neighborhood( typeListPoint3D *l )
{
  int i;
  if ( JCB_allocListPoint3D( l, 7 ) <= 0 ) {
    fprintf( stderr, "allocation failed\n" );
    return( -1 );
  }
  for ( i=0; i<7; i++ )
    l->pts[i].x = l->pts[i].y = l->pts[i].z = 0;
  l->pts[1].x =  1;
  l->pts[2].x = -1;
  l->pts[3].y =  1;
  l->pts[4].y = -1;
  l->pts[5].z =  1;
  l->pts[6].z = -1;
  return( 1 );
}

void JCB_transformListPoint3D( typeListPoint3D *theList, typeListPoint3D *resList, double *mat )
{
  int n;
  double x, y, z;


  for ( n=0; n<theList->n; n++ ) {
    x = mat[ 0] * theList->pts[n].x + mat[ 1] * theList->pts[n].y + mat[ 2] * theList->pts[n].z + mat[ 3];
    y = mat[ 4] * theList->pts[n].x + mat[ 5] * theList->pts[n].y + mat[ 6] * theList->pts[n].z + mat[ 7];
    z = mat[ 8] * theList->pts[n].x + mat[ 9] * theList->pts[n].y + mat[10] * theList->pts[n].z + mat[11];
    resList->pts[n].x = x;
    resList->pts[n].y = y;
    resList->pts[n].z = z;
  }
}

/********************************************************************************
 *
 *
 *
 ********************************************************************************/

static void JCB_ComputeCorrelationMatrix( double *mat, typeListPoint3D *r, typeListPoint3D *l )
{
  int i;
  for ( i=0; i<16; i++ ) mat[i]=0.0;
  
  for ( i=0; i<l->n; i++ ) {
    mat[ 0] += r->pts[i].x * l->pts[i].x;
    mat[ 1] += r->pts[i].x * l->pts[i].y;
    mat[ 2] += r->pts[i].x * l->pts[i].z;
    mat[ 3] += r->pts[i].x;

    mat[ 4] += r->pts[i].y * l->pts[i].x;
    mat[ 5] += r->pts[i].y * l->pts[i].y;
    mat[ 6] += r->pts[i].y * l->pts[i].z;
    mat[ 7] += r->pts[i].y;

    mat[ 8] += r->pts[i].z * l->pts[i].x;
    mat[ 9] += r->pts[i].z * l->pts[i].y;
    mat[10] += r->pts[i].z * l->pts[i].z;
    mat[11] += r->pts[i].z;

    mat[12] += l->pts[i].x;
    mat[13] += l->pts[i].y;
    mat[14] += l->pts[i].z;
    mat[15] += 1;

  }
}


/********************************************************************************
 *
 *
 *
 ********************************************************************************/


void JCB_TestTranslation( int n )
{
  double tmin = -4;
  double tmax =  8;
  int i;
  double mat[16];
  double mxx[16];
  double ixx[16];
  double myx[16];
  double prd[16];
  
  typeListPoint3D theList;
  typeListPoint3D resList;

  JCB_initListPoint3D( &theList );
  JCB_initListPoint3D( &resList );
  (void)JCB_get6Neighborhood( &theList );
  (void)JCB_allocListPoint3D( &resList, theList.n );

  for ( i=0; i<n; i++ ) {
    JCB_RndTranslationMatrix( mat, tmin, tmax );
    JCB_PrintMatrix( stdout, mat, NULL );

    JCB_transformListPoint3D( &theList, &resList, mat );
    
    JCB_ComputeCorrelationMatrix( mxx, &theList, &theList );
    JCB_ComputeCorrelationMatrix( myx, &resList, &theList );
    
    InverseMat4x4( mxx, ixx );
    Mat4x4ByMat4x4( myx, mxx, prd );
    JCB_PrintMatrix( stdout, mat, "Produit" );
    fprintf( stdout, "Jacobien = %g\n", JCB_ComputeJacobien( &resList, &theList ) );
  }

}


void JCB_TestSimilitude( int n )
{
  double smin =  0.5;
  double smax =  1.5;
  double tmin = -4;
  double tmax =  8;
  int i;
  double mat[16];
  double mxx[16];
  double ixx[16];
  double myx[16];
  double prd[16];
  
  typeListPoint3D theList;
  typeListPoint3D resList;

  JCB_initListPoint3D( &theList );
  JCB_initListPoint3D( &resList );
  (void)JCB_get6Neighborhood( &theList );
  (void)JCB_allocListPoint3D( &resList, theList.n );

  for ( i=0; i<n; i++ ) {
    JCB_RndSimilitudeMatrix( mat, smin, smax, tmin, tmax );
    JCB_PrintMatrix( stdout, mat, NULL );

    JCB_transformListPoint3D( &theList, &resList, mat );
    
    JCB_ComputeCorrelationMatrix( mxx, &theList, &theList );
    JCB_ComputeCorrelationMatrix( myx, &resList, &theList );
    
    InverseMat4x4( mxx, ixx );
    Mat4x4ByMat4x4( myx, mxx, prd );
    JCB_PrintMatrix( stdout, mat, "Produit" );
    fprintf( stdout, "Jacobien = %g\n", JCB_ComputeJacobien( &resList, &theList ) );
  }

}


void JCB_TestAffine( int n )
{
  double amin =  0.5;
  double amax =  1.5;
  double tmin = -4;
  double tmax =  8;
  int i;
  double mat[16];
  double mxx[16];
  double ixx[16];
  double myx[16];
  double prd[16];
  
  typeListPoint3D theList;
  typeListPoint3D resList;

  JCB_initListPoint3D( &theList );
  JCB_initListPoint3D( &resList );
  (void)JCB_get6Neighborhood( &theList );
  (void)JCB_allocListPoint3D( &resList, theList.n );

  for ( i=0; i<n; i++ ) {
    JCB_RndAffineMatrix( mat, amin, amax, tmin, tmax );
    JCB_PrintMatrix( stdout, mat, NULL );

    JCB_transformListPoint3D( &theList, &resList, mat );
    
    JCB_ComputeCorrelationMatrix( mxx, &theList, &theList );
    JCB_ComputeCorrelationMatrix( myx, &resList, &theList );
    
    InverseMat4x4( mxx, ixx );
    Mat4x4ByMat4x4( myx, mxx, prd );
    JCB_PrintMatrix( stdout, mat, "Produit" );
    fprintf( stdout, "Jacobien = %g\n", JCB_ComputeJacobien( &resList, &theList ) );
  }

}


/********************************************************************************
 *
 *
 *
 ********************************************************************************/

/* points resultats, points initiaux
 */
double JCB_ComputeJacobien( typeListPoint3D *resList, typeListPoint3D *theList )
{
  double mxx[16];
  double myx[16];
  double dxx, dyx;

  JCB_ComputeCorrelationMatrix( mxx, theList, theList );
  JCB_ComputeCorrelationMatrix( myx, resList, theList );

  dxx = mxx[ 0] * mxx[ 5] * mxx[10]
    + mxx[ 1] * mxx[ 6] * mxx[ 8]
    + mxx[ 2] * mxx[ 4] * mxx[ 9];
  dyx = myx[ 0] * myx[ 5] * myx[10]
    + myx[ 1] * myx[ 6] * myx[ 8]
    + myx[ 2] * myx[ 4] * myx[ 9];
  return( dyx / dxx );
}
