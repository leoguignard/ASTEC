/*************************************************************************
 * connexe.c -
 *
 * $Id: connexe.c,v 1.4 2003/06/20 09:05:09 greg Exp $
 *
 * Copyright (c) INRIA 1999
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * ?
 *
 * ADDITIONS, CHANGES
 *
 * - Tue Apr 11 19:07:13 MET DST 2000, Gregoire Malandain
 *   propagation de la taille du voxel
 */

#include <vt_common.h>

#include <vt_connexe.h>

#define _VT_CONNECTED  1
#define _VT_HYSTERESIS 2
#define _VT_SEEDPT     3
#define _VT_SEEDSIM    4

typedef struct local_par {
  vt_names names;
  vt_connexe cpar;
  char matlab_name[STRINGLENGTH];
  char matlab_file[STRINGLENGTH];
  char ascii_file[STRINGLENGTH];
  double time_step;
  double pixel_size;
} local_par;

/*------- Definition des fonctions statiques ----------*/
#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
#else 
static void VT_Parse();
static void VT_ErrorParse();
static void VT_InitParam();
#endif

static char *usage = "[image-in]\n\
\t [-matlab_file %s] [-matlab_name %s] [-ts %f] [-ps %f]\n\
\t [-ascii_file %s]\n";

static char *detail = "\
\t if 'image-in' is equal to '-', we consider stdin\n\
\t if 'image-out' is not specified, we consider stdout\n";

static char program[STRINGLENGTH];

typedef struct {
  int xmin, xmax;
  int ymin, ymax;
  int s;
  int v;
} bb2D;

#ifndef NO_PROTO
int main( int argc, char *argv[] )
#else
int main( argc, argv )
int argc;
char *argv[];
#endif
{
  local_par par;
  vt_image *image;
  vt_image imtmp;
  int retour;
  double *size = NULL;
  int dimy, dimz;
  
  FILE *f, *fopen();
  int i, x, z, fd;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  /* parametres */
  par.cpar.dim = VT_2D;
  par.cpar.type_connexite = N04;
  par.cpar.type_output = VT_GREY;
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  VT_InverseImage( image );
  
  VT_InitFromImage( &imtmp, image, "foo.inr", USHORT );
  if ( VT_AllocImage( &imtmp ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate auxiliary image\n", 0);
  }
  
  retour = VT_ConnectedComponents( image, &imtmp, 1, &(par.cpar) );
  VT_FreeImage( image );
  VT_Free( (void**)&image );
  
  dimy = imtmp.dim.y;
  dimz = imtmp.dim.z;
  size = (double*)malloc((imtmp.dim.z)*sizeof(double) );
  
  {
    int y, j, s, max;
    unsigned short int ***buf = (unsigned short int ***)imtmp.array;
    bb2D *boxes;
    
    s = imtmp.dim.x*imtmp.dim.y;
    
    for ( z=0; z<imtmp.dim.z; z++ ) {
      
      size[z] = 0;
      
      max = 0;
      for ( y=0; y<imtmp.dim.y; y++ ) 
	for ( x=0; x<imtmp.dim.x; x++ ) {
	  if ( buf[z][y][x] > max ) max = buf[z][y][x];
	}
      
      boxes = (bb2D*)malloc( (max+1)*sizeof(bb2D) );
      for ( i=0;i<=max;i++ ) {
	boxes[i].xmin = imtmp.dim.x - 1;
	boxes[i].xmax = 0;
	boxes[i].ymin = imtmp.dim.y - 1;
	boxes[i].ymax = 0;
	boxes[i].s = 0;
	boxes[i].v = 1;
      }
      
      for ( y=0; y<imtmp.dim.y; y++ ) 
	for ( x=0; x<imtmp.dim.x; x++ ) {
	  i = buf[z][y][x];
	  if ( i == 0 ) continue;
	  if ( boxes[i].xmin > x ) boxes[i].xmin = x;
	  if ( boxes[i].xmax < x ) boxes[i].xmax = x;
	  if ( boxes[i].ymin > y ) boxes[i].ymin = y;
	  if ( boxes[i].ymax < y ) boxes[i].ymax = y;
	  boxes[i].s ++;
	}
      
      j = 1;
      for ( i=1;i<=max;i++ ) {
	if ( boxes[j].s < boxes[i].s ) j = i;
      }
      
      if ( max > 1 ) {
	fprintf( stdout, "slice #%3d, %2d components, max=#%d\n", z, max, j );
      }
      
      size[z] = boxes[j].s;
      
      free( boxes );
    }
  }

	
  
  VT_FreeImage( &imtmp );

  
  
  /* fichier ascii
   */
  if ( par.ascii_file[0] != '\0' ) {
    f = fopen( par.ascii_file, "w" );
    for ( z=0; z<dimz; z++ ) 
      fprintf( f, "%d\n", (int)size[z] );
    fclose(f);
  }


  /* par coherence vis-a-vis d'avant, on divise par dimy
   */
  for ( z=0; z<dimz; z++ ) 
    size[z] /= (double)dimy;
	

  /* fichiers matlab
   */
  if ( par.matlab_file[0] != '>' && par.matlab_file[0] != '\0' ) {

    /* fichier raw
     */

    sprintf( par.names.ext, "%s.raw", par.matlab_file );

    fd = creat( par.names.ext, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  
    if ( write( fd, size, dimz*sizeof( double ) ) == -1 ) 
      fprintf( stderr, "error when writing in %s\n", par.names.ext );

    for ( z=1; z<dimz; z++ ) 
      size[z-1] -= size[z];
    size[dimz-1] = size[dimz-2];
    
    if ( write( fd, size, dimz*sizeof( double ) ) == -1 ) 
      fprintf( stderr, "error when writing in %s\n", par.names.ext );
    
    close( fd );

    /* fichier m
     */
    
    i = strlen( par.matlab_file )-1;
    for ( ; i >= 0 && par.matlab_file[i] != '/' ; i-- ) 
      ;

    sprintf( par.names.ext, "%s.m", par.matlab_file );
    
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
    fprintf( f, "fid = fopen('%s.raw', 'r' );\n", &(par.matlab_file[i+1]) );


    fprintf( f, "xaire_%s=0:%lf:%lf;\n", par.matlab_name, par.time_step, par.time_step*(dimz-1) );
    fprintf( f, "aire_%s = fread( fid, [%d], 'float%lu' );\n", 
	     par.matlab_name, dimz, 8*sizeof( double ) );
    fprintf( f, "daire_%s = fread( fid, [%d], 'float%lu' );\n", 
	     par.matlab_name, dimz, 8*sizeof( double ) );
    fprintf( f, "fclose( fid );\n" );
    
    fprintf( f, "\n" );
    fprintf( f, "mdaire_%s = mean( daire_%s );\n", par.matlab_name, par.matlab_name );
    fprintf( f, "rmdaire_%s = trimmean( daire_%s, 25 );\n", par.matlab_name, par.matlab_name );
    fprintf( f, "cfitaire_%s = polyfit( xaire_%s', aire_%s, 1 );\n",
	     par.matlab_name, par.matlab_name, par.matlab_name );
    fprintf( f, "rcfitaire_%s = robustfit( xaire_%s, aire_%s );\n",
	     par.matlab_name, par.matlab_name, par.matlab_name );
    fprintf( f, "\n" );
    
    fprintf( f, "disp(sprintf('Demi-pente   moyenne          = %%f', abs(cfitaire_%s(1)/2) ));\n", 
	     par.matlab_name );
    fprintf( f, "disp(sprintf('Demi-pente   moyenne (robuste)= %%f', abs(rcfitaire_%s(2)/2) ));\n", 
	     par.matlab_name );
    fprintf( f, "disp(sprintf('Demi-vitesse moyenne          = %%f', mdaire_%s/2 ));\n", par.matlab_name );
    fprintf( f, "disp(sprintf('Demi-vitesse moyenne (robuste)= %%f', rmdaire_%s/2 ));\n", par.matlab_name );
    
    fprintf( f, "\n" );
   
    
    fprintf( f, "figure;\n" );
    fprintf( f, "hold on;\n" );
    fprintf( f, "plot( xaire_%s, aire_%s );\n", par.matlab_name, par.matlab_name );
    fprintf( f, "plot( xaire_%s, xaire_%s*cfitaire_%s(1) + cfitaire_%s(2), 'g' );\n", 
	     par.matlab_name, par.matlab_name, par.matlab_name, par.matlab_name );
    fprintf( f, "plot( xaire_%s, xaire_%s*rcfitaire_%s(2) + rcfitaire_%s(1), 'r' );\n", 
	     par.matlab_name, par.matlab_name, par.matlab_name, par.matlab_name );
    fprintf( f, "legend('mesures','moindres carres','moindres carres robustes');\n" );
    fprintf( f, "hold off;\n" );
    
    fprintf( f, "\n" );
    fprintf( f, "figure;\n" );
    fprintf( f, "hold on;\n" );
    fprintf( f, "t_%s=0:%lf:%lf;\n", par.matlab_name, par.time_step, par.time_step*(dimz-1) );
    fprintf( f, "plot( xaire_%s, daire_%s );\n", par.matlab_name, par.matlab_name );
    fprintf( f, "plot( xaire_%s, ones(size(xaire_%s))*mdaire_%s, 'g' );\n", 
	     par.matlab_name, par.matlab_name, par.matlab_name );
    fprintf( f, "plot( xaire_%s, ones(size(xaire_%s))*rmdaire_%s, 'r' );\n", 
	     par.matlab_name, par.matlab_name, par.matlab_name );
    fprintf( f, "legend('mesures','moyenne','moyenne robuste');\n" );
    fprintf( f, "hold off;\n" );
    
    
    
    fclose( f );
  }

	free( size );

	return( 1 );
}

#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], local_par *par )
#else
static void VT_Parse( argc, argv, par )
int argc;
char *argv[];
local_par *par;
#endif
{
    int i, nb, status;
    int o=0, s=0, r=0;
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
		par->cpar.verbose = 1;
	    }
	    else if ( strcmp ( argv[i], "-D" ) == 0 ) {
		_VT_DEBUG_ = 1;
	    }


      else if ( strcmp ( argv[i], "-matlab_name" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -matlab_name...\n", 0 );
	strncpy( par->matlab_name, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-matlab_file" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -matlab_file...\n", 0 );
	strncpy( par->matlab_file, argv[i], STRINGLENGTH );  
      }

       else if ( strcmp ( argv[i], "-ascii_file" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -ascii_file...\n", 0 );
	strncpy( par->ascii_file, argv[i], STRINGLENGTH );  
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


	    /*--- lecture du type de l'image de sortie ---*/
	    else if ( strcmp ( argv[i], "-r" ) == 0 ) {
		r = 1;
	    }
	    else if ( strcmp ( argv[i], "-s" ) == 0 ) {
		s = 1;
	    }
	    else if ( strcmp ( argv[i], "-o" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -o...\n", 0 );
		status = sscanf( argv[i],"%d",&o );
		if ( status <= 0 ) VT_ErrorParse( "parsing -o...\n", 0 );
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

#ifndef NO_PROTO
static void VT_ErrorParse( char *str, int flag )
#else
static void VT_ErrorParse( str, flag )
char *str;
int flag;
#endif
{
	(void)fprintf(stderr,"Usage : %s %s\n",program, usage);
        if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
        (void)fprintf(stderr,"Erreur : %s",str);
        exit(0);
}

#ifndef NO_PROTO
static void VT_InitParam( local_par *par )
#else
static void VT_InitParam( par )
local_par *par;
#endif
{
	VT_Names( &(par->names) );
	VT_Connexe( &(par->cpar) );
  par->matlab_file[0] = '\0';
  par->matlab_name[0] = '\0';
  par->ascii_file[0] = '\0';
  par->time_step = 1;
  par->pixel_size = 1;
}
