/*************************************************************************
 * compute-speed.c -
 *
 * $Id: compute-speed.c,v 1.2 2005/06/24 13:01:05 greg Exp $
 *
 * Copyright (c) INRIA 2005
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * 
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>

typedef struct local_par {
  vt_names names;
  char imx[STRINGLENGTH];
  char imy[STRINGLENGTH];
  char imz[STRINGLENGTH];
  char matlab[STRINGLENGTH];
  double time_step;
  double pixel_size;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-x %s] [-y %s] [-z %s] [-ts %lf] [-ps %lf]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t options-de-type : -o 1    : unsigned char\n\
\t                   -o 2    : unsigned short int\n\
\t                   -o 2 -s : short int\n\
\t                   -o 4 -s : int\n\
\t                   -r      : float\n\
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n\
 $Revision: 1.2 $ $Date: 2005/06/24 13:01:05 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, *imx, *imy, *imz;
  float ***bufx, ***bufy, ***bufz;
  unsigned char ***buf;
  double *m, *ec;
  double vx, vy, vz, vn, va;
  int dimz, x, y, z, n;
  double sum;

  FILE *f, *fopen();
  int i, fd;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);

  imx = _VT_Inrimage( par.imx );
  if ( imx == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input X image\n", 0);
  imy = _VT_Inrimage( par.imy );
  if ( imy == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input Y image\n", 0);
  imz = _VT_Inrimage( par.imz );
  if ( imz == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input Z image\n", 0);
  
  if ( image->type != UCHAR ||
       imx->type != FLOAT ||
       imy->type != FLOAT ||
       imz->type != FLOAT ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_FreeImage( imx );
    VT_Free( (void**)&imx );
    VT_FreeImage( imy );
    VT_Free( (void**)&imy );
    VT_FreeImage( imy );
    VT_Free( (void**)&imy );
    VT_ErrorParse("image types not handled\n", 0);
  }


  buf = (unsigned char***)(image->array);
  bufx = (float***)(imx->array);
  bufy = (float***)(imy->array);
  bufz = (float***)(imz->array);

  
  m = (double*)malloc( image->dim.z * sizeof(double) );
  ec = (double*)malloc( image->dim.z * sizeof(double) );

  dimz = image->dim.z;

  for ( z=0; z<image->dim.z; z++ ) {

    sum = 0;
    n = 0;

    if ( _VT_VERBOSE_ ) {
      fprintf( stderr, "%3d", z );
    }

    for ( y=0; y<image->dim.y; y++ )
    for ( x=0; x<image->dim.x; x++ ) {

      if ( buf[z][y][x] == 0 ) continue;
      n++;

      vx = bufx[z][y][x];      
      vy = bufy[z][y][x];
      vz = bufz[z][y][x];
      va = sqrt( vx*vx + vy*vy  + vz*vz  );
      vn = sqrt( vx*vx + vy*vy  );

      vz /= vn;

      /* fprintf( stdout, "[%f/%f %f] ", vn,va,vz ); */
      
      /* vitesse dans le plan 
	 = produit scalaire du gradient avec z 
	 = vz
      */

      sum += vz ;
      
    }
 
    m[z] = sum / n;

    if ( _VT_VERBOSE_ ) {
      fprintf( stderr, " %5d %f", n, m[z] );
    }

    sum = 0;

    for ( y=0; y<image->dim.y; y++ )
    for ( x=0; x<image->dim.x; x++ ) {

      if ( buf[z][y][x] == 0 ) continue;

      vx = bufx[z][y][x];      
      vy = bufy[z][y][x];
      vz = bufz[z][y][x];
      vn = sqrt( vx*vx + vy*vy  );
 
      vz /= vn;

      sum += ( m[z] - vz )*( m[z] - vz );
      
    }

    ec[z] = sqrt( sum / n );

    if ( _VT_VERBOSE_ ) {
      fprintf( stderr, " %f\n", ec[z] );
    }

  }

    
  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_Free( (void**)&image );

  VT_FreeImage( imx );
  VT_Free( (void**)&imx );
  VT_FreeImage( imy );
  VT_Free( (void**)&imy );
  VT_FreeImage( imz );
  VT_Free( (void**)&imz );


  if ( par.names.out[0] != '>' && par.names.out[0] != '\0' )
    sprintf( par.names.ext, "%s.raw", par.names.out );
  else
    sprintf( par.names.ext, "foo.raw" );
    
  fd = creat( par.names.ext, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  if ( write( fd, m, dimz*sizeof( double ) ) == -1 ) 
    fprintf( stderr, "error when writing 'measure' values in %s\n", par.names.ext );
  if ( write( fd, ec, dimz*sizeof( double ) ) == -1 ) 
    fprintf( stderr, "error when writing 'standard deviation' values in %s\n", par.names.ext );
  close( fd );

  i = strlen( par.names.out )-1;
  for ( ; i >= 0 && par.names.out[i] != '/' ; i-- ) 
    ;

  if ( par.names.out[0] != '>' && par.names.out[0] != '\0' )
    sprintf( par.names.ext, "%s.m", par.names.out );
  else
    sprintf( par.names.ext, "foo.m"  );

  f = fopen( par.names.ext, "w" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "%%\n" );
  fprintf( f, "%%" );
  for ( x=0; x<argc; x++ ) fprintf( f, " %s", argv[x] );
  fprintf( f, "\n" );
  fprintf( f, "%%\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );

  fprintf( f, "echo off\n" );
  if ( par.names.out[0] != '>' && par.names.out[0] != '\0' )
    fprintf( f, "fid = fopen('%s.raw', 'r' );\n", &(par.names.out[i+1]) );
  else
    fprintf( f, "fid = fopen('foo.raw', 'r' );\n" );

  fprintf( f, "\n" );
  fprintf( f, "xgrad_%s=0:%lf:%lf;\n", par.matlab, par.time_step, par.time_step*(dimz-1) );
  fprintf( f, "vgrad_%s = fread( fid, [%d], 'float%lu' );\n", 
	    par.matlab, dimz, 8*sizeof( double ) );
  fprintf( f, "ecgrad_%s = fread( fid, [%d], 'float%lu' );\n", 
	    par.matlab, dimz, 8*sizeof( double ) );
  fprintf( f, "fclose( fid );\n" );
  fprintf( f, "\n" );
  fprintf( f, "mgrad_%s  = mean(vgrad_%s);\n", par.matlab, par.matlab );
  fprintf( f, "rmgrad_%s = trimmean(vgrad_%s, 25);\n", par.matlab, par.matlab );
  fprintf( f, "\n" );
  fprintf( f, "disp(sprintf('Vitesse moyenne          = %%f', mgrad_%s ));\n", par.matlab );
  fprintf( f, "disp(sprintf('Vitesse moyenne (robuste)= %%f', rmgrad_%s ));\n", par.matlab );
  fprintf( f, "\n" );
  
  fprintf( f, "figure;\n" );
  fprintf( f, "hold on;\n" );
  fprintf( f, "errorbar( xgrad_%s, vgrad_%s, ecgrad_%s );\n", 
	   par.matlab, par.matlab, par.matlab );
  fprintf( f, "plot( xgrad_%s, ones(size(xgrad_%s))*mgrad_%s, 'g' );\n", 
	   par.matlab, par.matlab, par.matlab );
  fprintf( f, "plot( xgrad_%s, ones(size(xgrad_%s))*rmgrad_%s, 'r' );\n", 
	   par.matlab, par.matlab, par.matlab );
  fprintf( f, "legend('mesures','ecart-type','moyenne','moyenne robuste');\n" );
  fprintf( f, "hold off;\n" );
  fprintf( f, "\n" );
  fclose( f );

  free( m );
  free( ec );
  

  return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status;
  char text[STRINGLENGTH];
  
  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
  i = 1; nb = 0;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( argv[i][1] == '\0' ) {
	if ( nb == 0 ) {
	  /*--- standart input ---*/
	  strcpy( par->names.in, "<" );
	  nb += 1;
	}
      }
      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
	VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
	_VT_VERBOSE_ = 1;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
      }

      
      else if ( strcmp ( argv[i], "-x" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -x...\n", 0 );
	strncpy( par->imx, argv[i], STRINGLENGTH );  
      }
      else if ( strcmp ( argv[i], "-y" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -y...\n", 0 );
	strncpy( par->imy, argv[i], STRINGLENGTH );  
      }
      else if ( strcmp ( argv[i], "-z" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -y...\n", 0 );
	strncpy( par->imz, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-matlab" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -matlab...\n", 0 );
	strncpy( par->matlab, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-ts" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -ts...\n", 0 );
	status = sscanf( argv[i], "%lf", &(par->time_step) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -ts...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-ps" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -ps...\n", 0 );
	status = sscanf( argv[i], "%lf", &(par->pixel_size) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -ps...\n", 0 );
      }


      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }

      /*--- option inconnue ---*/
      else {
	sprintf(text,"unknown option %s\n",argv[i]);
	VT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( nb == 0 ) { 
	strncpy( par->names.in, argv[i], STRINGLENGTH );  
	nb += 1;
      }
      else if ( nb == 1 ) {
	strncpy( par->names.out, argv[i], STRINGLENGTH );  
	nb += 1;
      }
      else 
	VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }
  
  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
    strcpy( par->names.out, ">" );  /* standart output */
  
}






static void VT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}








static void VT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );
  par->imx[0] = '\0';
  par->imy[0] = '\0';
  par->imz[0] = '\0';
  par->matlab[0] = '\0';
  par->time_step = 1;
  par->pixel_size = 1;
}
