#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*--- aide ---*/
static char *usage = "2 usages:\n-inv trsf1 [inverse trsf name (default: inv)]\n-n [nb trsfs] trsf1 trsf2 ... [result transf name (default: compotrsfs.trsf)]";
static char *detail ="\t Compose transformations";


/******************************************
 Procedures de manipulation de matrices 4x4 
*******************************************/

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
  return( rang );
}


/* Programme principal */

int main(int argc, char *argv[]){
  
  FILE *file1;
  double Trsfs[500][4][4], Trsf_tmp[4][4], Trsf_res[4][4];
  int i, nb_trsfs;
  char la_trsf[500];
  int status;
  
  /*--- y-a-t'il le bon nombre d'arguments ? 
    s'il n'y en a pas, on affiche l'aide complete ---*/
  if ( argc == 1 ) {
    printf("\n%s\n%s\n\n", usage, detail );
    exit(0);
  }
  
  /* On verifie les arguments */
  if ( (strcmp(argv[1], "-n") != 0) && (strcmp(argv[1], "-inv") != 0) ){
    printf( "\nWrong argument (-inv OR -n nb_trsfs)\n%s\n\n", usage );
    exit(2);
  }
  
  /* Cas 1 : pour inverser une matrice */
  if (strcmp(argv[1], "-inv") == 0){
    i = 0;
    file1 = fopen(argv[2+i], "r");
    fprintf(stderr, "Lecture en cours de  %s\n", argv[2+i]);
    status = fscanf(file1, "(\nO8\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf", 
	   &Trsf_tmp[0][0], &Trsf_tmp[0][1], &Trsf_tmp[0][2], &Trsf_tmp[0][3], 
	   &Trsf_tmp[1][0], &Trsf_tmp[1][1], &Trsf_tmp[1][2], &Trsf_tmp[1][3], 
	   &Trsf_tmp[2][0], &Trsf_tmp[2][1], &Trsf_tmp[2][2], &Trsf_tmp[2][3], 
	   &Trsf_tmp[3][0], &Trsf_tmp[3][1], &Trsf_tmp[3][2], &Trsf_tmp[3][3]);
    if ( status != 16 ){
      fprintf(stderr, "Erreur dans la lecture de %d\n", i+1 );
      exit(2);
    }
    fclose(file1);
    
    MatInverse(Trsf_tmp, Trsf_res);

    /* On ecrit la transformation dans un fichier */
    if ( argc == 4 )
      sprintf(la_trsf, "%s.trsf", argv[argc-1]);
    else 
      strcpy(la_trsf, "inv.trsf");
    file1=fopen(la_trsf, "w");
    if ( file1 == NULL ) {
      fprintf( stderr, "Unable to open '%s'\n", la_trsf );
    }
    fprintf(file1,"(\nO8\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n)\n",
	    Trsf_res[0][0],  Trsf_res[0][1], Trsf_res[0][2], Trsf_res[0][3], 
	    Trsf_res[1][0],  Trsf_res[1][1], Trsf_res[1][2], Trsf_res[1][3], 
	    Trsf_res[2][0],  Trsf_res[2][1], Trsf_res[2][2], Trsf_res[2][3], 
	    Trsf_res[3][0],  Trsf_res[3][1], Trsf_res[3][2], Trsf_res[3][3]);
    fclose(file1);
  }
  /* Cas 2 : pour composer des matrices */
  else if (strcmp(argv[1], "-n") == 0){
    /* Nombre de transformations */
    nb_trsfs = atoi(argv[2]);
    
    if ( argc == (nb_trsfs+4) )
      strcpy(la_trsf, argv[argc-1]);
    else 
      strcpy(la_trsf, "compotrsfs.trsf");
    
    /* On lit les matrices et on les stocke dans Trsfs */
    for (i=0; i<nb_trsfs; i++){
      file1 = fopen(argv[3+i], "r");
      if ( file1 == NULL ) {
	fprintf( stderr, "Unable to open '%s'\n", la_trsf );
      }
      fprintf(stderr, "Lecture en cours de  %s\n", argv[3+i]);
      status = fscanf(file1, "(\nO8\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf", 
		      &Trsfs[i][0][0], &Trsfs[i][0][1], &Trsfs[i][0][2], &Trsfs[i][0][3], 
		      &Trsfs[i][1][0], &Trsfs[i][1][1], &Trsfs[i][1][2], &Trsfs[i][1][3], 
		      &Trsfs[i][2][0], &Trsfs[i][2][1], &Trsfs[i][2][2], &Trsfs[i][2][3], 
		      &Trsfs[i][3][0], &Trsfs[i][3][1], &Trsfs[i][3][2], &Trsfs[i][3][3]);
      if ( status != 16 ){
	fprintf(stderr, "Erreur dans la lecture de %d\n", i+1 );
	exit(2);
      }
      fclose(file1);
    }
    
    /* Calcul de la matrice de transformation = trn*...*tr2*tr1 */
    Trsf_res[0][0] = 1.0; Trsf_res[0][1] = 0.0; Trsf_res[0][2] = 0.0;
    Trsf_res[0][3] = 0.0;
    Trsf_res[1][0] = 0.0; Trsf_res[1][1] = 1.0; Trsf_res[1][2] = 0.0;
    Trsf_res[1][3] = 0.0;
    Trsf_res[2][0] = 0.0; Trsf_res[2][1] = 0.0; Trsf_res[2][2] = 1.0;
    Trsf_res[2][3] = 0.0;
    Trsf_res[3][0] = 0.0; Trsf_res[3][1] = 0.0; Trsf_res[3][2] = 0.0;
    Trsf_res[3][3] = 1.0;
    
    for (i=0; i<nb_trsfs; i++){
      MatMult(Trsfs[i], Trsf_res, Trsf_tmp);
      MatCopy(Trsf_tmp, Trsf_res);
    }
    
    /* On ecrit la transformation dans un fichier */
    file1=fopen(la_trsf, "w");
    if ( file1 == NULL ) {
      fprintf( stderr, "Unable to open '%s'\n", la_trsf );
    }
    fprintf(file1, "(\nO8\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n)\n",
	    Trsf_res[0][0],  Trsf_res[0][1], Trsf_res[0][2], Trsf_res[0][3], 
	    Trsf_res[1][0],  Trsf_res[1][1], Trsf_res[1][2], Trsf_res[1][3], 
	    Trsf_res[2][0],  Trsf_res[2][1], Trsf_res[2][2], Trsf_res[2][3], 
	    Trsf_res[3][0],  Trsf_res[3][1], Trsf_res[3][2], Trsf_res[3][3]);
    fclose(file1);
  }
    return 0;
}

