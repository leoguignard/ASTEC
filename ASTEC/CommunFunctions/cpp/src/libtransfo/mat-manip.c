#include <stdlib.h>

/******************************************
 Procedures de manipulation de matrices 4x4 
*******************************************/
int MatRead( char *name, double mat[4][4] ) 
{
  FILE *file;
  int status;
  
  if ( name == NULL || name[0] == '\0' ) return( 0 );

  file = fopen( name, "r" );
    
  if ( file == NULL ) {
    if ( 0 ) fprintf( stderr, "Erreur a l'ouverture de '%s'\n", name );
    return( 0 );
  }
    
  status = fscanf(file , "(\nO8\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf", 
		  &mat[0][0], &mat[0][1], &mat[0][2], &mat[0][3], 
		  &mat[1][0], &mat[1][1], &mat[1][2], &mat[1][3], 
		  &mat[2][0], &mat[2][1], &mat[2][2], &mat[2][3], 
		  &mat[3][0], &mat[3][1], &mat[3][2], &mat[3][3]);
  if ( status != 16 ){
    if ( 0 ) fprintf( stderr, "Erreur a la lecture de '%s'\n", name );
    return( -1 );
  }
  
  fclose( file );
  return( 1 );
}


void MatPrintf( FILE *f, double mat[4][4] ) 
{
  fprintf(f , "(\nO8\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n)\n", 
	  mat[0][0], mat[0][1], mat[0][2], mat[0][3], 
	  mat[1][0], mat[1][1], mat[1][2], mat[1][3], 
	  mat[2][0], mat[2][1], mat[2][2], mat[2][3], 
	  mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
  
}


int MatWrite( char *name, double mat[4][4] ) 
{
  FILE *file;
  
  if ( name == NULL || name[0] == '\0' ) return( 0 );

  file = fopen( name, "w" );
  
  if ( file == NULL ) {
    if ( 0 ) fprintf( stderr, "Erreur a l'ouverture de '%s'\n", name );
    return( 0 );
  }
  
  MatPrintf( file, mat );

  fclose( file );
  return( 1 );
}





void MatIdentity( double mat[4][4] )
{
  int i, j;

  for (i=0; i<4; i++)
  for (j=0; j<4; j++)
    mat[i][j] = 0.0;
  for (i=0; i<4; i++)
    mat[i][i] = 1.0;
}


/* Copy mat1 into mat2 */ 
void MatCopy(double mat1[4][4], double mat2[4][4] )
{
  int i, j;

  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      mat2[i][j] = mat1[i][j];
}

/* Left-hand multiply mat1 with mat2 (mat1*mat2) and put the result into mat3 */
void MatMult(double mat1[4][4], double mat2[4][4], double mat3[4][4])
{
  int i, j, k;

  for (i=0; i<4; i++)
    for (j=0; j<4; j++){
      mat3[i][j] = 0;
      for (k=0; k<4; k++)
	mat3[i][j] += mat1[i][k]*mat2[k][j];
    }
}

/* Inverse mat1 into mat2 */
#define TINY 1e-12
int MatInverse(double mat1[4][4], double mat2[4][4])
{
  register int i, j, k;
  int kmax, rang = 4;
  register double c, max;
  double mat[4][4];

  MatPrintf( stdout, mat1 );

  for (i=0; i<4; i++ ) 
    for (j=0; j<4; j++) {
      mat[i][j] = mat1[i][j] ;
      mat2[i][j] = 0.0;
    }
  mat2[0][0] = mat2[1][1] = mat2[2][2] = mat2[3][3] = 1.0;
  
  for ( j=0; j<4; j++ ) {
    if ( (mat[j][j] > (-TINY)) && (mat[j][j] < TINY) ) {
      /* recherche du plus grand element non nul sur la colonne j */
      kmax = j;
      max = 0.0;
      for (k=j+1; k<4; k++ ) {
        c = ( mat[j][k] > 0.0 ) ? mat[j][k] : (-mat[j][k]) ;
        if ( (c > TINY) && (c > max) ) { max = c; kmax = k; }
      }
      if ( kmax == j ) {
        /* la ligne est nulle */
        rang --;
      } else {
        /* sinon, on additionne */
        for ( i=0; i<4; i++ ) {
          mat[i][j] += mat[i][kmax];
          mat2[i][j] += mat2[i][kmax];
        }
      }
    }
    if ( (mat[j][j] < (-TINY)) || (mat[j][j] > TINY) ) {
      /* les autres lignes */
      for (k=0; k<4; k++) {
        if ( k != j ) {
          c = mat[j][k] / mat[j][j];
          for ( i=0; i<4; i++ ) {
            mat[i][k] -= c * mat[i][j];
            mat2[i][k] -= c * mat2[i][j];
          }
        }
      }
      /* la ligne */
      c = mat[j][j];
      for ( i=0; i<4; i++ ) {
        mat[i][j] /= c;
        mat2[i][j] /= c;
      }
    }
  }
  MatPrintf( stdout, mat2 );

  return( rang );
}
