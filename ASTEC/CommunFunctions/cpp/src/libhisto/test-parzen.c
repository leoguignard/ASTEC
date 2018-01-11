/*************************************************************************
 * minimum.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 1999
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
#include <vt_sliceHisto.h>



typedef enum {
  _NONE_,
  _GAUSSIANS_,
  _JOINTHISTO_,
  _ENTROPY_,
  _SSD_,
  _SAD_,
  _CORRELATION_
} enumComputation;





typedef struct local_par {
  vt_names names;
  double sigma;
  int zref;
  int type;
  enumComputation typeComputation;
} local_par;







/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-matlab %s]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

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
 $Revision: 1.2 $ $Date: 2001/09/28 17:09:55 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;
  int **theHisto;
  double **thePDF;
  int i, iref;

  int intensity_max = 0;

  typeProbabilite *newPDF;
  int newn;
  double binsize = 1.0;
  

  int fd = 0;
  FILE *fp = NULL;
  char *longname;


  /*--- initialisation des parametres ---*/
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
  

  theHisto = _GetSlicesHisto( image, &intensity_max );
  for (i=0; i<image->dim.z; i++ )
    theHisto[i][0] = 0;
  thePDF = _GetPDFFromHisto( theHisto, image->dim.z, intensity_max+1 );






  if ( par.zref <= 0 || par.zref >= image->dim.z )
    iref = _GetReferenceSlice( theHisto, image->dim.z, intensity_max+1 );
  else
    iref = par.zref;
  printf( "coupe de reference = %d\n", iref );


  if ( par.names.ext[0] != '\0' ) {
    int i, k;
    
    longname = malloc( strlen(par.names.ext)+ 10 );
    sprintf( longname, "%s.raw", par.names.ext );
    fd = creat( longname, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
    sprintf( longname, "%s.m", par.names.ext );
    fp = fopen( longname, "w" );
    free( longname );

    fprintf( fp, "\n" );
    k = strlen( par.names.ext );
    for ( i = k-1; i >= 0 && par.names.ext[i] != '/' ; i-- )
      ;
    fprintf( fp, "fid = fopen('%s.raw', 'r' );\n", &(par.names.ext[i+1]) );
    fprintf( fp, "\n" );
  }

  /* on a besoin de ca (ecriture de l'histo complet
   */
  _Print2DSlicesPDFForMatlab( fd, fp, thePDF, image->dim.z, intensity_max+1, 0 );
  /* pour le plot c'est iref+1, qui correspond
     a notre iref a nous
  */
  _PrintOneSlicePDFForMatlab( fd, fp, NULL, thePDF, iref+1, intensity_max+1, 0 );


  printf("\n sigma = %f \n", par.sigma );
  newPDF = _buildOneNewPDF( binsize, intensity_max+1, &newn );
  _newPDFFromPDF( newPDF, newn, thePDF[iref], intensity_max+1, par.sigma );
  _PrintOnePDFForMatlab( fd, fp, newPDF, newn, iref );



  if ( par.names.ext[0] != '\0' ) {
    fprintf( fp, "\n" );
    fprintf( fp, "fclose(fid);\n" );
    fclose( fp );
    close( fd );
  }




  free( thePDF);
  free( theHisto );
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
      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }





      else if ( strcmp ( argv[i], "-mode" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -mode...\n", 0 );
	if ( strcmp ( argv[i], "gauss" ) == 0 ) {
	  par->typeComputation = _GAUSSIANS_;
	} 
	else if ( strcmp ( argv[i], "joint" ) == 0 ) {
	  par->typeComputation = _JOINTHISTO_;
	} 
	else if ( strcmp ( argv[i], "entropy" ) == 0 ||
		  strcmp ( argv[i], "ent" ) == 0 ) {
	  par->typeComputation = _ENTROPY_;
	} 
	else if ( strcmp ( argv[i], "ssd" ) == 0 ) {
	  par->typeComputation = _SSD_;
	} 
	else if ( strcmp ( argv[i], "sad" ) == 0 ) {
	  par->typeComputation = _SAD_;
	} 
	else if ( strcmp ( argv[i], "corr" ) == 0 || strcmp ( argv[i], "cor" ) == 0 ) {
	  par->typeComputation = _CORRELATION_;
	} 
	else {
	  VT_ErrorParse( "unknown mode...\n", 0 );
	}
	  
      }

      else if ( strcmp ( argv[i], "-matlab" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -matlab...\n", 0 );
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
      }



      else if ( strcmp ( argv[i], "-zref" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -zref...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->zref) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -zref...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-sigma" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -sigma...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->sigma) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -sigma...\n", 0 );
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
  par->type = TYPE_UNKNOWN;
  par->sigma = 5.0;
  par->zref = -1;
  par->typeComputation = _NONE_;

}
