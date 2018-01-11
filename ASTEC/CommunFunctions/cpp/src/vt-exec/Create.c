/*************************************************************************
 * Create.c -
 *
 * $Id: Create.c,v 1.3 2001/04/13 18:12:12 greg Exp $
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

typedef struct local_par {
  vt_names names;
  vt_4vpt dim;
  vt_fpt voxel;
  vt_fpt offset;
    int type;
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

static char *usage = "[image-out] [-x %d] [-y %d] [-z %d] [-v %d]\n\
[-vs %f %f %f] [-tx %f] [-ty %f] [[-tz %f] options-de-type]\n\
\t [-w] [-D] [-help]";

static char *detail = "\
\t si 'image-out' est absent, on prendra stdout\n\
\t -x %d      : taille de l'image selon X (1 par defaut)\n\
\t -y %d      : taille de l'image selon Y (1 par defaut)\n\
\t -z %d      : taille de l'image selon Z (1 par defaut)\n\
\t -v %d      : taille de l'image selon V (1 par defaut)\n\
\t -vs %f %f %f : voxel size\n\
\t options-de-type : -o 1    : unsigned char\n\
\t                   -o 1 -s :   signed char\n\
\t                   -o 2    : unsigned short int\n\
\t                   -o 2 -s :   signed short int\n\
\t                   -o 4 -s : int\n\
\t                   -r      : float\n\
\t si aucune de ces options n'est presente, on prend le type 'unsigned char'\n\
\n\
 $Revision: 1.3 $ $Date: 2001/04/13 18:12:12 $ $Author: greg $\n";

static char program[STRINGLENGTH];

int main( int argc, char *argv[] )
{
	local_par par;
	vt_image image;
	char header[257];
	int fd ,nb;

	/*--- initialisation des parametres ---*/
	VT_InitParam( &par );
	/*--- lecture des parametres ---*/
	VT_Parse( argc, argv, &par );

	/*--- verification des parametres ---*/
	if ( (par.dim.x <= 0) || (par.dim.y <= 0) || (par.dim.z <= 0) || (par.dim.v <= 0) )
	    VT_ErrorParse("unable to deal with negative dimensions\n", 0);
	
	/*--- initialisation de l'image auxiliaire ---*/
	VT_Image( &image );
	VT_InitVImage( &image, par.names.out, par.dim.v, par.dim.x, par.dim.y, par.dim.z, par.type );
	image.siz.x = par.voxel.x;
	image.siz.y = par.voxel.y;
	image.siz.z = par.voxel.z;
	VT_SetImageOffset( &image, par.offset.x, par.offset.y, par.offset.z );
	/*--- creation de l'header ---*/
	VT_FillInrimHeader( header, &image );
	/*--- opening the file ---*/
	fd = VT_WOpen( image.name );
	if ( fd == -1 ) {
	    VT_ErrorParse("Unable to open file for writing", 0);
	}
	/*--- writing the header ---*/
	nb = VT_Write( fd, header, 256 );
	if ( nb == -1 ) {
	    VT_Close( fd );
	    VT_ErrorParse("error when writing the header", 0);
	}
	if ( nb < 256 ) {
	    VT_Close( fd );
	    VT_ErrorParse("not enough space left for writing the header", 0);
	}

	VT_Close( fd );
	
	exit( 0 );
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
		VT_ErrorParse("parsing - \n", 1);
	    }
	    /*--- options generales ---*/
	    else if ( strcmp ( argv[i], "-help" ) == 0 ) {
		VT_ErrorParse("\n", 1);
	    }
	    else if ( strcmp ( argv[i], "-w" ) == 0 ) {
		_VT_VERBOSE_ = 1;
	    }
	    else if ( strcmp ( argv[i], "-D" ) == 0 ) {
		_VT_DEBUG_ = 1;
	    }
	    /*--- dimension de l'image ---*/
	    else if ( strcmp ( argv[i], "-x" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -x...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->dim.x) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -x...\n", 0 );
	    }
	    else if ( strcmp ( argv[i], "-y" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -y...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->dim.y) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -y...\n", 0 );
	    }
	    else if ( strcmp ( argv[i], "-z" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -z...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->dim.z) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -z...\n", 0 );
	    }
	    else if ( strcmp ( argv[i], "-tx" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -tx...\n", 0 );
		status = sscanf( argv[i],"%f",&(par->offset.x) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -tx...\n", 0 );
	    }
	    else if ( strcmp ( argv[i], "-ty" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -ty...\n", 0 );
		status = sscanf( argv[i],"%f",&(par->offset.y) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -ty...\n", 0 );
	    }
	    else if ( strcmp ( argv[i], "-tz" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -tz...\n", 0 );
		status = sscanf( argv[i],"%f",&(par->offset.z) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -tz...\n", 0 );
	    }
	    else if ( strcmp ( argv[i], "-v" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -v...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->dim.v) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -v...\n", 0 );
	    }

	    else if ( strcmp ( argv[i], "-vs" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -vs..\n", 0 );
		status = sscanf( argv[i],"%f",&(par->voxel.x) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -vs...\n", 0 );
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -vs..\n", 0 );
		status = sscanf( argv[i],"%f",&(par->voxel.y) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -vs...\n", 0 );
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -vs..\n", 0 );
		status = sscanf( argv[i],"%f",&(par->voxel.z) );
		if ( status <= 0 ) VT_ErrorParse( "parsing -vs...\n", 0 );
	    }

	    /*--- lecture du type de l'image ---*/
	    else if ( strcmp ( argv[i], "-f" ) == 0 ) {
	      /*--- rien : fixed ---*/
	    }
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
	    else {
		sprintf(text,"unknown option %s\n",argv[i]);
		VT_ErrorParse(text, 0);
	    }
	}
	else if ( argv[i][0] != 0 ) {
	    if ( nb == 0 ) { 
		strncpy( par->names.out, argv[i], STRINGLENGTH );  
		nb += 1;
	    }
	    else 
		VT_ErrorParse("too much file names when parsing\n", 0 );
	}
	i += 1;
    }
    if (nb == 0) {
	strcpy( par->names.out, ">" );  /* standart output */
    }
    
    /*--- type de l'image resultat ---*/
    if ( (o == 1) && (s == 1) && (r == 0) )  par->type = SCHAR;
    if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
    if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
    if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
    if ( (o == 4) && (s == 1) && (r == 0) )  par->type = SINT;
    if ( ((o == 0) || (o == 4)) && (r == 1) )  par->type = FLOAT;
    if ( (o == 8) && (r == 1) )  par->type = DOUBLE;
    if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program);
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
	par->dim.x = par->dim.y = par->dim.z = par->dim.v = 1;
	par->voxel.x = par->voxel.y = par->voxel.z = 1.0;
	par->offset.x = par->offset.y = par->offset.z = 0.0;
	par->type = UCHAR;
}
