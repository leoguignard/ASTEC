#include<epidaure.h>
#include<epidaure.ee>

/* Definition de type */
typedef struct Monpoint {
  double valx, valy, error;
} point;

/*--- Definition des fonctions statiques ---*/
#ifndef NO_PROTO
static void E_ErrorParse( char *str, int l );
#else 
static void E_ErrorParse();
#endif

static  int errorcompare( point *pt1, point *pt2)
{
  if ((*pt1).error > (*pt2).error)
    return (1);
  if ((*pt1).error < (*pt2).error)
    return (-1);
  return (0);
}

static  double dabs( double val)
{
  if (val > 0)
    return val;
  else
    return -val;
}

/*--- aide ---*/
static char *usage = "[image1-in] [image2-in] [name image out(2)] [max background in image1 (often 10)] [max background in image2 (often 10)] [ratio of rejected point at each iteration] [error on a parameter] [error on b parameter]";
static char *detail ="\t\n";


#if defined(_ANSI_)
void main( int argc, char *argv[] )
#else
void main( argc, argv )
int argc;
char *argv[];
#endif
{
    /*--- t_Image est une structure image ---*/
    t_Image theIm1, theIm2, resIm;
    int i, j, k, iteration;
    long cpt;
    double x, y, moy_x, moy_y, moy_xy, moy_yy, moy_xx;
    double max_im;
    double n, n_best, max1, max2, a, b, A, a_err, b_err, ratio, old_a, old_b;
    point *bonpoints;

    /*--- y-a-t'il le bon nombre d'arguments ? 
          s'il n'y en a pas, on affiche l'aide complete ---*/
    if ( argc == 1 ) E_ErrorParse( "", TRUE );
    if ( argc < 9 )  E_ErrorParse( " not enough arguments", FALSE );
    if ( argc > 9 )  E_ErrorParse( " too many arguments", FALSE );

    max1 = atof(argv[4]);
    max2 = atof(argv[5]);
    ratio = atof(argv[6]);
    a_err = atof(argv[7]);
    b_err = atof(argv[8]);    

    /*--- initialisation de la structure image entree ---*/
    T_Init3Image( &theIm1, 0, 0, 0, TYPEIM_UNKNOWN, argv[1] );
    T_Init3Image( &theIm2, 0, 0, 0, TYPEIM_UNKNOWN, argv[2] );

    /*--- lecture de l'image entree ---*/
    T_ReadInrim( &theIm1 );
    T_ReadInrim( &theIm2 );

    puts("Images-in read");

    
    if (theIm1.type != theIm2.type)
      puts("Images must have same format");

    max_im = 0;
    if (theIm1.type == TYPEIM_256)
      max_im = 255.0;
    else if (theIm1.type == TYPEIM_16B)
      max_im = 65535.0; 
    else if (theIm1.type == TYPEIM_S16B)
      max_im = 37767.0; 

    fprintf( stderr, "Image 1 '%25s' %d x %d x %d fond=%f\n",
	     theIm1.fi.name, theIm1.tx, theIm1.ty, theIm1.tz, max1 );
    fprintf( stderr, "Image 2 '%25s' %d x %d x %d fond=%f\n",
	     theIm2.fi.name, theIm2.tx, theIm2.ty, theIm2.tz, max2 );

    /*--- initialisation de l'image sortie :
          on lui donne les memes dimensions 
          que l'image d'entree et on lui donne
	  le type sur un octet ---*/
    T_Init3Image( &resIm, theIm1.tx, theIm1.ty, theIm1.tz, theIm1.type, argv[3] );

    /*--- allocation de l'image sortie ---*/
    T_AllocImage( &resIm );

    /* Allocation du tableau de points */
    bonpoints = malloc(theIm1.tx*theIm1.ty*theIm1.tz*sizeof(point));

    puts("Image-out allocated");

    puts("Computation of an initialization of a and b coefficients");

    moy_x = 0.0;
    moy_y = 0.0;
    moy_xy = 0.0;
    moy_xx = 0.0;
    moy_yy = 0.0;
    n=0.0;
    for ( i = 0; i < theIm1.tx; i++ )
      for ( j = 0; j < theIm1.ty; j++ )
	for ( k = 0; k < theIm1.tz; k++ )
	  {
	    /*
	    y = max_im * T_Get3Point( i, j, k, &theIm1 );
	    x = max_im * T_Get3Point( i, j, k, &theIm2 );
	    */
	    y = (double)T_GetIPoint( i, j, k, &theIm1 );
	    x = (double)T_GetIPoint( i, j, k, &theIm2 );
	    
	    /* on ne tient pas compte des voxels du fond */
	    if ((max2 < x) && (max1 < y))
	      {
		/*
		fprintf( stdout, "i=%d/%d j=%d/%d k=%d/%d ",
			 i, theIm1.tx, j, theIm1.ty, k, theIm1.tz );
		fprintf( stdout, "(%f<%f)&&(%f<%f) \n", max1,x,max2,y );
		*/
		/* On ajoute le point a la liste de points */
		bonpoints[(long)n].valx = x;
		bonpoints[(long)n].valy = y;
		moy_x += x;
		moy_y += y;
		moy_xy += x*y;
		moy_xx += x*x;
		moy_yy += y*y;
		n++;
	      }
	  }

    if ( n < 1.0 ) {
      puts("no points in images" );
      exit( -1 );
    }

    moy_x  /= (double)n;
    moy_y  /= (double)n;
    moy_xx /= (double)n;
    moy_xy /= (double)n;
    moy_yy /= (double)n;

    if ( 0 )fprintf( stderr, " %g %g %g %g %g\n", moy_x, moy_y, moy_xy, moy_xx, moy_yy );

    /* Calcul de A */
    A = 0.5 * 
      ( (moy_yy - moy_y*moy_y) - (moy_xx - moy_x*moy_x) ) /
      ( moy_x*moy_y - moy_xy );

    /* Calcul du coeff dir a */
    if (moy_xy - moy_x*moy_y > 0)
      a = -A + sqrt(A*A+1);
    else
      a = -A - sqrt(A*A+1);

    /* Calcul de l'ordonnee a l'origine */
    b = moy_y - a*moy_x;
    b = 0.0;

    printf("Initial affine bias : y = %f * x + %f\n", a, b);

    /* Iteration 0 pour le moment */
    iteration = 0;

    /* On calcule sur combien de points on fait faire l'estimation robuste */
    n_best = n * (1 - ratio) ;
    printf("We will use only %ld on %ld points\n", (long)n_best, (long)n );

    /* On initialise de maniere a forcer la premiere iteration */
    old_a = a - a_err - 1;
    old_b = b - b_err - 1;

    /* Boucle de traitements */
    /*    while ( (dabs(a-old_a) > a_err) || (dabs(b-old_b) > b_err) )  */
    while ( (dabs(a-old_a) > a_err) ) 
      {	
	iteration++;
	printf("Iteration no %d\n", iteration);
	printf("a error = %f\nb error = %f\n", dabs(a-old_a), dabs(b-old_b));
	old_a = a;
	old_b = b;

	/* Calcul des erreurs */
	for ( cpt = 0 ; cpt < (long)n; cpt++) 
	  {
	    bonpoints[cpt].error = abs ( max_im * bonpoints[cpt].valy - 
		       ( a * max_im * bonpoints[cpt].valx + b ) ) /
		     (1 + b * b );
	  }
	
	/* Trie des erreurs */
	qsort( bonpoints, (long)n, sizeof(point), errorcompare );
	
	/* Estimation sur les "meilleurs" points */
	moy_x = 0.0;
	moy_y = 0.0;
	moy_xy = 0.0;
	moy_xx = 0.0;
	moy_yy = 0.0;
	n = 0.0;
	for ( cpt = 0; cpt < (long)n_best; cpt++ )
	  {
	    x = bonpoints[cpt].valx;
	    y = bonpoints[cpt].valy;
	    moy_x = (n * moy_x + x) / (n+1);
	    moy_y = (n * moy_y + y) / (n+1);
	    moy_xy = (n * moy_xy + x*y) / (n+1);
	    moy_xx = (n * moy_xx + x*x) / (n+1);
	    moy_yy = (n * moy_yy + y*y) / (n+1);
	    n += 1;
	  }
	
	/* Calcul de A */
	A = 0.5 * 
	  ( (moy_yy - moy_y*moy_y) - (moy_xx - moy_x*moy_x) ) /
	  ( moy_x*moy_y - moy_xy );
	
	/* Calcul du coeff dir a */
	if (moy_xy - moy_x*moy_y > 0)
	  a = -A + sqrt(A*A+1);
	else
	  a = -A - sqrt(A*A+1);
	
	/* Calcul de l'ordonnee a l'origine */
	b = moy_y - a*moy_x;
	b = 0.0;
	
	printf("Affine bias at iteration %d : y = %f * x + %f\n", iteration, a, b);
	
      }
    
    puts("Computation of image-out");
    
    /* Calcul de l'image de sortie */
    for ( i = 0; i < theIm1.tx; i++ )
      for ( j = 0; j < theIm1.ty; j++ )
	for ( k = 0; k < theIm1.tz; k++ )
	  {
	    /*
	    x = T_Get3Point( i, j, k, &theIm2 );
	    */
	    x = (double)T_GetIPoint( i, j, k, &theIm2 );
	    y = a*x + b;
	    if ( y < 0 ) y = 0;
	    /* on ne corrige pas les voxels du fond */
	    T_Set3Point( i, j, k, &resIm, y/max_im, (int)MODE_SET );
	  }
    
    puts("Image-out computed");

    /*--- la multiplication est realisee :
          on ecrit l'image sortie ---*/
    T_WriteInrim( &resIm );
/*     T_WriteInrim( &resIm_real ); */
 
    puts("Image-out written");

    /*--- on desalloue les images ---*/
    T_FreeImage( &theIm1 );
    T_FreeImage( &theIm2 );
    T_FreeImage( &resIm );
/*     T_FreeImage( &resIm_real ); */

    puts("Memory free");
}

/* Fonction d'erreur.

   Si une erreur est trouvee lors de l'examen des
   arguments, on fait appel a cette fonction (locale).

   Le message contenu dans "str" est un message explicatif
   de l'erreur.

   Si "flag" est a TRUE, un message d'aide plus complet
   est affiche.
*/
#if defined(_ANSI_)
static void E_ErrorParse( char *str, int flag )
#else
static void E_ErrorParse( str, flag )
char *str;
int flag;
#endif
{

        (void)fprintf(stderr,"Usage %s : %s\n",program, usage);
        if ( flag == TRUE ) (void)fprintf(stderr,"%s",detail);
        T_Err( str, CNULL );
}
