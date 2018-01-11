/*************************************************************************
 * membrane.c -
 *
 * $Id: membrane.c,v 2.0 2013/10/22 14:22:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 * 
 * CREATION DATE: 
 * 2013/06/17
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <mt_membrane3D.h>
#include <connexe.h>
#include <time.h>

#include <chunks.h>

typedef enum {
  MULTI_SCALE,
  SINGLE_SCALE
} enumComputation;

typedef struct local_par {
  vt_names names;
  int type;

  int writeImages;
  enumComputation typeComputation;
  enumMode mode;

  enumStructureColor structureColor;

  double scale1;
  double scale2;
  int nbscales;

  double zfact;  
  
  int dimension;
  int hsp;

} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-init %lf] [-last %lf] [-zfact %lf] [-nb %d] [-wi] [-2D]\n\
\t [-single] [-wi] [-black|-noir|-white|-blanc] [-modelbased|-acme]\n\
\t [-hsp %d]\n\
\t [-parallel|-no-parallel] [-max-chunks %d]\n\
\t [-parallel-scheduling|-ps default|static|dynamic-one|dynamic|guided]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t [-modelbased|-acme] : selectionne le mode de calcul de la reponse\n\
\t -init %lf : echelle (sigma) fine\n\
\t -last %lf : echelle (sigma) grossiere\n\
\t -zfact %lf: rapport d'échelle entre l'axe z et les axes x et y\n\
\t -nb   %d  : nombre d'echelles (interpolation logarithmique)\n\
\t -hsp  %d  : half size plane (pour le cas modelbased : integration sur \n\
\t des petits plans plutot que seulement sur 2 points\n\
\t -wi : ecrit toutes les images intermediaires\n\
\t       theta, phi : angles definissant le vecteur 3D\n\
\t       rep : reponse 3D multi-echelle\n\
\t       scale : echelle de la reponse maximale\n\
\t -single : une seule echelle de calcul (spécifiee dans -init)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;
  vt_3Dimres imsRes;
  int flag_3D = 1;

  clock_t start, stop;
  double elapsed;

  start = clock();


  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    MT_ErrorParse("unable to read input image\n", 0);
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );
  

  if ( par.dimension == 2 || image->dim.z == 1 )
    flag_3D = 0;
  
  if ( VT_Alloc3DImres( &imsRes, image, par.names.out, flag_3D ) != 1 ) {
    MT_ErrorParse("unable to allocate response images\n", 0 );
  }


  /*--- calculs ---*/

  switch ( par.typeComputation ) {
  
  case SINGLE_SCALE :

      par.scale2 = par.scale1;
      par.nbscales = 1;
    
  case MULTI_SCALE :
        
    if ( flag_3D ) {
      // 3D case
      if ( MT_Compute3DMultiScale( image, &imsRes, par.scale1,
			   par.scale2, par.nbscales, par.zfact,
			   par.structureColor, par.mode, par.hsp ) != 1 ) {
	MT_ErrorParse("unable to compute response\n", 0 );
      }
      if ( par.writeImages ) VT_Write3DImres( &imsRes, flag_3D );
      
      MT_Compute3DExtrema( &imsRes, par.zfact, &(imsRes.imTheta) );
      sprintf( imsRes.imTheta.name, "%s.ext.inr", par.names.out );
      VT_WriteInrimage( &(imsRes.imTheta) );

    }
    else {
      // 2D case
      if ( MT_Compute2DMultiScale( image, &imsRes, par.scale1,
				   par.scale2, par.nbscales, 
				   par.structureColor ) != 1 ) {
	MT_ErrorParse("unable to compute response\n", 0 );
      }
      if ( par.writeImages ) VT_Write3DImres( &imsRes, flag_3D );

      MT_Compute2DExtrema( &imsRes, &(imsRes.imTheta) );
      sprintf( imsRes.imTheta.name, "%s.ext.inr", par.names.out );
      VT_WriteInrimage( &(imsRes.imTheta) );
    }
    
    break;

  default :
    {
    	MT_ErrorParse("incorrect typeComputation selected\n", 0 );
    }
  }
  


  
  /*--- liberations memoires ---*/
  VT_Free3DImres( &imsRes );
  VT_FreeImage( image );
  VT_Free( (void**)&image );

  stop = clock();
  elapsed = (double)(stop-start)/CLOCKS_PER_SEC;
  fprintf(stdout, "Elapsed time : \t%.1fn\n", elapsed);

  return( 1 );
}








static void MT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status, tmp;
  char text[STRINGLENGTH];
  
  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) MT_ErrorParse("\n", 0 );
  
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
	MT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
	_VT_VERBOSE_ = 1;
	VT_IncrementVerboseInVtTube3D(  );
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
	VT_IncrementDebugInVtTube3D(  );
      }

      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }



      
      else if ( strcmp ( argv[i], "-single" ) == 0 ) {
	par->typeComputation = SINGLE_SCALE;
      }

      else if ( strcmp ( argv[i], "-acme" ) == 0 ) {
        par->mode = ACME;
      }

      else if ( strcmp ( argv[i], "-modelbased" ) == 0 ) {
	par->mode = MODELBASED;
      }

      
      else if ( strcmp ( argv[i], "-wi" ) == 0 ) {
	par->writeImages = 1;
      }



      else if ( strcmp ( argv[i], "-black" ) == 0 
		|| strcmp ( argv[i], "-noir" ) == 0 ) {
	par->structureColor = _BLACK_;
      }
      else if ( strcmp ( argv[i], "-white" ) == 0 
		|| strcmp ( argv[i], "-blanc" ) == 0 ) {
	par->structureColor = _WHITE_;
      }
  


      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
        par->dimension = 2;
      }



      /* Parametres de calcul */
      else if ( strcmp ( argv[i], "-init" ) == 0 ) {
	i += 1;
	if ( i >= argc)    MT_ErrorParse( "parsing -init...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->scale1) );
	if ( status <= 0 ) MT_ErrorParse( "parsing -init...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-last" ) == 0 ) {
	i += 1;
	if ( i >= argc)    MT_ErrorParse( "parsing -last...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->scale2) );
	if ( status <= 0 ) MT_ErrorParse( "parsing -last...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-nb" ) == 0 ) {
	i += 1;
	if ( i >= argc)    MT_ErrorParse( "parsing -nb...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->nbscales) );
	if ( status <= 0 ) MT_ErrorParse( "parsing -nb...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-zfact" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -last...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->zfact) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -last...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-hsp" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -hsp...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->hsp) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -hsp...\n", 0 );
      }

      /* parallelism
       */
      else if ( strcmp ( argv[i], "-parallel" ) == 0 ) {
        setMaxChunks( 100 );
      }

      else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
        setMaxChunks( 1 );
      }

      else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-max-chunks", 0 );
        status = sscanf( argv[i], "%d", &tmp );
        if ( status <= 0 ) MT_ErrorParse( "-max-chunks", 0 );
        if ( tmp >= 1 ) setMaxChunks( tmp );
      }

      else if ( strcmp ( argv[i], "-parallel-scheduling" ) == 0 ||
                ( strcmp ( argv[i], "-ps" ) == 0 && argv[i][3] == '\0') ) {
        i ++;
        if ( i >= argc)    MT_ErrorParse( "-parallel-scheduling", 0 );
        if ( strcmp ( argv[i], "default" ) == 0 ) {
          setOpenMPScheduling( _DEFAULT_SCHEDULING_ );
        }
        else if ( strcmp ( argv[i], "static" ) == 0 ) {
          setOpenMPScheduling( _STATIC_SCHEDULING_ );
        }
        else if ( strcmp ( argv[i], "dynamic-one" ) == 0 ) {
          setOpenMPScheduling( _DYNAMIC_ONE_SCHEDULING_ );
        }
        else if ( strcmp ( argv[i], "dynamic" ) == 0 ) {
          setOpenMPScheduling( _DYNAMIC_SCHEDULING_ );
        }
        else if ( strcmp ( argv[i], "guided" ) == 0 ) {
          setOpenMPScheduling( _GUIDED_SCHEDULING_ );
        }
        else {
          MT_ErrorParse( "-parallel-scheduling", 0 );
        }
      }



      /*--- option inconnue ---*/
      else {
	sprintf(text,"unknown option %s\n",argv[i]);
	MT_ErrorParse(text, 0);
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
	MT_ErrorParse("too much file names when parsing\n", 0 );
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







static void MT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}








static void MT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );
  par->type = TYPE_UNKNOWN;

  par->writeImages = 0;
  par->typeComputation = 0;
  par->mode = 0;

  par->structureColor = _WHITE_;

  par->scale1 = 1.0;
  par->scale2 = 1.0;
  par->nbscales = 1;
  par->zfact = 1.0;
  
  par->dimension = 3;
  par->hsp = 0;

}
