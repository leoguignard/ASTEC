/*************************************************************************
 * parcellisation.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2008
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

#include <parcelling.h>

typedef struct local_par {
  vt_names names;
  char seedfile[STRINGLENGTH];
  char seedout[STRINGLENGTH];
  char seedimg[STRINGLENGTH];
  int nparcels;
  int type;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-dist %s] [-iterations|-i %d] [-force] [-seeds %s] [-sout %s] [-wi %s]\n\
\t [-parcels|-p %d] [-random-seed|-rs %ld]\n\
\t [-inv] [-swap] [-v] [-nv] [-D] [-help] [options-de-type]";

static char *detail = "\
\t -parcels|-p %d: number of parcels\n\
\t -dist %s: output distance image\n\
\t -iterations|-i %d: maximal number of iterations\n\
\t -force: force iterative computation of parcel center\n\
\t -seeds: list of seed points\n\
\t         each line is of the form 'X Y Z'\n\
\t -sout %s: write list of seeds computed in %s\n\
\t -random-seed %ld: specify random seed\n\
\t -wi %s: seeds written in binary image file %s\n\
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
 $Revision: 1.5 $ $Date: 2000/08/16 16:31:56 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;
  vt_image imres;
  vt_image imdist;
  int theDim[3];

  FILE *fopen(),*fseeds;
  int n;
  int nallocated = 10000;
  int **theSeeds = NULL;
  int *seeds = NULL;
  int ret;
  int inSeeds = 0;

  /*--- initialisation des parametres ---*/
  _VT_VERBOSE_ = 1;
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );

  theDim[0] = image->dim.x;
  theDim[1] = image->dim.y;
  theDim[2] = image->dim.z;



  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  if ( par.nparcels < 256 && par.seedfile[0] == '\0')
    VT_InitFromImage( &imres, image, par.names.out, UCHAR );
  else 
    VT_InitFromImage( &imres, image, par.names.out, USHORT );

  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  
  
  VT_Image( &imdist );
  if ( par.names.ext[0] != '\0' ) {
    VT_InitFromImage( &imdist, image, par.names.ext, USHORT );
    if ( VT_AllocImage( &imdist ) != 1 ) {
      VT_FreeImage( &imres );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to allocate distance image\n", 0);
    }
  }

  if ( par.seedfile[0] != '\0' ) {
    inSeeds = 1;
    fseeds = fopen( par.seedfile, "r" );
    if ( fseeds == NULL ) {
      if ( par.names.ext[0] != '\0' ) 
	VT_FreeImage( &imdist );
      VT_FreeImage( &imres );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("error when opening seeds file\n", 0);
    }
    
    seeds = (int*)malloc( 3*nallocated * sizeof( int ) );
    theSeeds = (int**)malloc( nallocated * sizeof( int* ) );
    if ( seeds == NULL || theSeeds == NULL ) {
      if ( seeds != NULL ) free( seeds );
      if ( theSeeds != NULL ) free( theSeeds );
      if ( par.names.ext[0] != '\0' ) 
	VT_FreeImage( &imdist );
      VT_FreeImage( &imres );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("error when allocating seeds arrays\n", 0);
    }

    for ( n=0;n<nallocated;n++ ) {
      theSeeds[n] = &(seeds[3*n]);
    }
    
    n = 0;
    while ( (ret = fscanf( fseeds, "%d %d %d\n", &theSeeds[n][0], &theSeeds[n][1], &theSeeds[n][2]) ) != EOF ) {
      if ( ret != 3 ) {
	free( seeds );
	free( theSeeds );
	if ( par.names.ext[0] != '\0' ) 
	  VT_FreeImage( &imdist );
	VT_FreeImage( &imres );
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("error in reading seeds file\n", 0);
      }
      n++;
    }
    fclose( fseeds );
    
    fprintf( stderr, "%s: read %d seeds in file '%s'\n", program, n, par.seedfile );
    
    if ( n == 0 ) {
      free( seeds );
      free( theSeeds );
      if ( par.names.ext[0] != '\0' ) 
	VT_FreeImage( &imdist );
      VT_FreeImage( &imres );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("empty seeds file ?\n", 0);
    }
    par.nparcels = n;
  }
  else {
    if( par.seedimg[0] != '\0' || par.seedout[0] != '\0') {
      seeds = (int*)malloc( 3*par.nparcels * sizeof( int ) );
      theSeeds = (int**)malloc( par.nparcels * sizeof( int* ) );
      if ( seeds == NULL || theSeeds == NULL ) {
        if ( seeds != NULL ) free( seeds );
        if ( theSeeds != NULL ) free( theSeeds );
        if ( par.names.ext[0] != '\0' )
          VT_FreeImage( &imdist );
        VT_FreeImage( &imres );
        VT_FreeImage( image );
        VT_Free( (void**)&image );
        VT_ErrorParse("error when allocating seeds arrays\n", 0);
      }

      for ( n=0;n<par.nparcels;n++ ) {
        theSeeds[n] = &(seeds[3*n]);
      }
    }
  }

  if ( par.names.ext[0] != '\0' ) {
    if ( parcelling( image->buf, image->type, 
		     theSeeds, par.nparcels,
		     imres.buf, imres.type, 
		     imdist.buf, imdist.type, 
		     theDim, inSeeds ) != 1 ) {
      VT_FreeImage( &imdist );
      VT_FreeImage( &imres );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      if ( seeds != NULL ) free( seeds );
      if ( theSeeds != NULL ) free( theSeeds );
      VT_ErrorParse("error when processing\n", 0);
    }
  }
  else {
    if ( parcelling( image->buf, image->type, 
                     theSeeds, par.nparcels,
		     imres.buf, imres.type, 
		     NULL, TYPE_UNKNOWN, 
		     theDim, inSeeds ) != 1 ) {
      VT_FreeImage( &imres );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      if ( seeds != NULL ) free( seeds );
      if ( theSeeds != NULL ) free( theSeeds );
      VT_ErrorParse("error when processing\n", 0);
    }	   
  }






  /*--- ecriture de l'image resultat ---*/


  /* faudrait ecrire les points 
   */

  if ( (seeds != NULL) || ( theSeeds != NULL )) {
    if ( par.seedout[0] != '\0') {
      fseeds = fopen( par.seedout, "w" );

      if ( fseeds == NULL ) {
        if ( par.names.ext[0] != '\0' )
          VT_FreeImage( &imdist );
        VT_FreeImage( &imres );
        VT_FreeImage( image );
        VT_Free( (void**)&image );
        free( seeds );
        free( theSeeds );
        VT_ErrorParse("error when writing seeds output file\n", 0);
      }

      for ( n=0; n<par.nparcels; n++ ) {
        fprintf( fseeds, "%d %d %d\n", theSeeds[n][0], theSeeds[n][1], theSeeds[n][2] );
      }
    }
    else {
      for ( n=0; n<par.nparcels; n++ ) {
        fprintf( stdout, "#%03d: %d %d %d\n", n, theSeeds[n][0], theSeeds[n][1], theSeeds[n][2] );
      }
    }
    free( seeds );
    free( theSeeds );
  }

  if ( par.names.ext[0] != '\0' ) {
    if ( VT_WriteInrimage( &imdist ) == -1 ) {
      VT_FreeImage( &imdist );
      VT_FreeImage( &imres );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to write output distance image\n", 0);
    }
    VT_FreeImage( &imdist );
  }

  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( &imres );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( &imres );
  VT_FreeImage( image );
  VT_Free( (void**)&image );
  return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status;
  int o=0, s=0, r=0;
  int iterations;
  char text[STRINGLENGTH];
  long int randomseed;

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
	_VT_VERBOSE_ ++ ;
	parcelling_setverbose();
      }
      else if ( strcmp ( argv[i], "-nv" ) == 0 ) {
	_VT_VERBOSE_ = 0;
	parcelling_setnoverbose();
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
      }

      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }





      /*---  ---*/

      else if ( strcmp ( argv[i], "-dist" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -dist...\n", 0 );
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
      }
      
      
      
      else if ( strcmp ( argv[i], "-seeds" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -seeds...\n", 0 );
        strncpy( par->seedfile, argv[i], STRINGLENGTH );
      }

      else if ( strcmp ( argv[i], "-sout" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -sout...\n", 0 );
        strncpy( par->seedout, argv[i], STRINGLENGTH );
      }

      else if ( strcmp ( argv[i], "-wi" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -wi...\n", 0 );
        strncpy( par->seedimg, argv[i], STRINGLENGTH );
      }

      
      
      else if ( strcmp ( argv[i], "-parcels" ) == 0 || strcmp ( argv[i], "-p" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -parcels...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->nparcels) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -parcels...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-random-seed" ) == 0 || strcmp ( argv[i], "-rs" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -random-seed...\n", 0 );
	status = sscanf( argv[i],"%ld",&randomseed );
	if ( status <= 0 ) VT_ErrorParse( "parsing -random-seed...\n", 0 );
	parcelling_setRandomSeed( randomseed );
      }

      else if ( strcmp ( argv[i], "-iterations" ) == 0 || strcmp ( argv[i], "-i" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -iterations...\n", 0 );
	status = sscanf( argv[i],"%d",&iterations );
	if ( status <= 0 ) VT_ErrorParse( "parsing -iterations...\n", 0 );
	parcelling_setNumberOfIterations( iterations );
      }

      else if ( strcmp ( argv[i], "-force" ) == 0 ) {
	parcelling_ForceExactCenterCalculation();
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
  
  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 1) && (r == 0) )  par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = SINT;
  if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
  /* if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program); */
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
  par->nparcels = 100;
  par->type = TYPE_UNKNOWN;
  par->seedfile[0] = '\0';
  par->seedout[0] = '\0';
  par->seedimg[0] = '\0';
}
