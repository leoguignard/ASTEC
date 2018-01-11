/*************************************************************************
 * inr2interfile.c
 *
 * $Id: inr2interfile.c,v 1.2 2000/05/03 17:41:32 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue May  2 11:50:15 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>



#define descNumber       66
#define descValueLength 256

static char descValue[descNumber][descValueLength];

static char descInBuffer[4097];
static char descOutBuffer[4097];

static char *descInterfile[] = { "INTERFILE := ",
			  "imaging modality := ",
			  "originating system := ",
			  "version of keys := ",
			  "date of keys := ",
			  "conversion program := ",
			  "program author := ",
			  "program version := ",
			  "program date := ",
			  "GENERAL DATA := ",
			  "original institution := ",
			  "contact person := ",
			  "data description := ",
			  "data starting block := ",
			  "name of data file := ",
			  "patient name := ",
			  "patient ID := ",
			  "patient dob := ",
			  "patient sex := ",
			  "study ID := ",
			  "exam type := ",
			  "data compression := ",
			  "data encode := ",
			  "GENERAL IMAGE DATA := ",
			  "type of data := ",
			  "total number of images := ",
			  "study date := ",
			  "study time := ",
			  "imagedata byte order := ",
			  "number of energy windows := ",
			  "energy window [1] := ",
			  "flood corrected := ",
			  "decay corrected := ",
			  "SPECT STUDY (general) := ",
			  "number of images/energy window := ",
			  "matrix size [1] := ",
			  "matrix size [2] := ",
			  "number format := ",
			  "number of bytes per pixel := ",
			  "number of projections := ",
			  "extent of rotation := ",
			  "process status := ",
			  "time per projection (sec) := ",
			  "study duration (sec) := ",
			  "scaling factor (mm/pixel) [1] := ",
			  "scaling factor (mm/pixel) [2] := ",
			  "maximum pixel count := ",
			  "patient orientation := ",
			  "patient rotation := ",
			  "SPECT STUDY (reconstructed data) := ",
			  "method of reconstruction := ",
			  "number of slices := ",
			  "number of reference frame := ",
			  "slice orientation := ",
			  "slice thickness (pixels) := ",
			  "center-center slice separation (pixels) := ",
			  "filter name := ",
			  "filter parameters := ",
			  "z-axis filter := ",
			  "attenuation correction coefficient/cm := ",
			  "method of attenuation correction := ",
			  "scatter corrected := ",
			  "method of scatter correction := ",
			  "oblique reconstruction := ",
			  "oblique orientation := ",
			  "END OF INTERFILE := " };

typedef struct local_par { 
  vt_names names;
  int type;
} local_par;






/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );








static char *usage = "[image-in] [image-out]\n\
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
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n";

static char program[STRINGLENGTH];





int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;
  int i, j, size;
  int ifd, max = 0;
  char *theBuf;
  char tmpBuf[STRINGLENGTH];

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);


  for ( i=0; i<4097;i++ )
    descOutBuffer[i] = '\0';

  for ( i=0; i<descNumber; i++ ) 
  for (j=0; j<descValueLength; j++ )
    descValue[i][j] = '\0';
    
    
  ifd = VT_ROpen( par.names.ext );
  if ( ifd == -1 ) {
    VT_ErrorParse("Unable to open model file for reading",0 );
  }
	
  size = VT_Read( ifd, descInBuffer, 4096 );
  if ( size != 4096 ) {
    printf( " erreur lors de la lecture de '%s'\n", par.names.ext );
    exit( 0 );
  }
  VT_Close( ifd );



  theBuf = descInBuffer;
  for ( i=0; i<descNumber; i++ ) {
    if ( strncmp( theBuf, descInterfile[i], strlen( descInterfile[i] ) ) == 0 ) {
      theBuf += strlen( descInterfile[i] );
      j = 0;
      do {
	descValue[i][j] = *theBuf;
	theBuf ++;
	j++;
      } while( *theBuf != '\n' );
      theBuf ++;
    }
  }
  


  switch ( image->type ) {
  default :
    break;
  case USHORT :
    {
      unsigned short int *buf = image->buf;
      max = buf[0];
      for ( i=1; i<image->dim.x*image->dim.y*image->dim.z; i++ )
	if ( max < buf[i] ) max = buf[i];
    }
  }



  theBuf = descOutBuffer;
  for ( i=0; i<descNumber; i++ ) {
    if ( max > 0 &&
	 strncmp( descInterfile[i], "maximum pixel count := ", 
		  strlen( "maximum pixel count := " ) ) == 0 ) {
      sprintf( tmpBuf, "%s%d\r\n", descInterfile[i], max );
    }
    else if ( strncmp( descInterfile[i], "study ID := ",
		       strlen( "study ID := " ) ) == 0 ) {
      sprintf( tmpBuf, "%sHMPAORECAL\r\n", descInterfile[i] );
    }
    else if ( strncmp( descInterfile[i], "name of data file := ",
		       strlen( "name of data file := " ) ) == 0 ) {
      for ( j = strlen( par.names.out )-1; j >= 0 && par.names.out[j] != '/' ; j-- )
	;
      sprintf( tmpBuf, "%s%s\r\n", descInterfile[i], &(par.names.out[j+1]) );
    }
    else {
      for (j=0; j<STRINGLENGTH; j++ ) tmpBuf[j] = '\0';
      sprintf( tmpBuf, "%s%s\n", descInterfile[i], descValue[i] );
    }
    sprintf( theBuf, "%s", tmpBuf );
    theBuf += strlen( tmpBuf );
  }
  descOutBuffer[4095] = 032;

  ifd = VT_WOpen( par.names.out );
  VT_Write( ifd, descOutBuffer, 4096 );
  size = VT_SizeImage( image );
  (void)VT_Write( ifd, (char*)(image->buf), size );
  VT_Close( ifd );
    
  return( 1 );
}















static void VT_Parse( int argc, char *argv[], local_par *par )
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


      else if ( strcmp ( argv[i], "-mod" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -mod...\n", 0 );
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
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
}
