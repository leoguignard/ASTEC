/*************************************************************************
 * meristemeFormationAxes.c -
 *
 * $Id: meristemeFormationAxes.c,v 1.0 2014/04/24 18:40:00 gael Exp $
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/04/24
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <mt_membrane3D.h>
#include <vt_meristemeFormationAxes.h>

#include <time.h>





typedef struct local_par {
  vt_names names;
  int r;
  double alpha;
  int tensor;
  int writeImages;
  double scale;
  int NangleIter;
} local_par;



/*------- Definition des fonctions statiques ----------*/
static void MT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void MT_InitParam( local_par *par );



static char *usage = "[image-in] [image-out]\n\
\t [-r %d] [-tensor | -t] [-alpha | -a %f [-nangles|-n %d]]\n\
\t [-wi] [-v] [-D] [-help]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -r %d : rayon de vote (toute l'image si non specifie)\n\
\t -alpha %d : semi-angle du cone votant (default=0)\n\
\t -nangles %d : nombre d'iterations pour calculer des directions de champs\n\
\t -tensor | -t : vote de type tenseur\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  char name[STRINGLENGTH];
  local_par par;
  vt_image *imageIn;
  vt_image *imageTht;
  vt_image *imagePhi;
  vt_image imres;
  vt_image *theIm;
  vt_image* imagesIn[3];

  int dim[3];
  int pt[3];
  bufferType t;

  int flag_3D = 1;
  int marker = 0;
  
  char prefix[STRINGLENGTH];
  int i, j, k;

  float ***theTht=NULL, ***thePhi=NULL;
  float ***theIn32=NULL;
  unsigned char ***theInU8=NULL;
  double coef, tht, phi;
  
  clock_t start, stop;
  double elapsed;

  start = clock();


  /*--- initialisation des parametres ---*/
  MT_InitParam( &par );
  //fprintf(stdout, "alpha = %lf \ttensor = %d \tr = %d\n", par.alpha, par.tensor, par.r);

  /*--- lecture des parametres ---*/
  MT_Parse( argc, argv, &par );
  
  /*--- lecture des images d'entree ---*/
  sprintf( prefix, "%s", par.names.in );
  i=0;
  while (prefix[i]!='\0')
  {
    i++;
  }
  if (i>4 && prefix[i-4]=='.' &&
      ((prefix[i-3]=='h' && prefix[i-2]=='d' && prefix[i-1]=='r') ||
       (prefix[i-3]=='i' && prefix[i-2]=='n' && prefix[i-1]=='r')))
  {
    prefix[i-4]='\0';
    sprintf( name, "%s", par.names.in );
    marker = 1;
  }
  else
  {
    sprintf( name, "%s.inr", prefix );
  }
  imageIn = _VT_Inrimage( name );
  if ( imageIn == (vt_image*)NULL )
  {
	if ( marker == 1 )
	  MT_ErrorParse("unable to read input image A\n", 0);	  
    sprintf( name, "%s.hdr", prefix );	
	imageIn = _VT_Inrimage( name );
	if ( imageIn == (vt_image*)NULL )
	  MT_ErrorParse("unable to read input image A\n", 0);	  
  }
  if(_VT_VERBOSE_)
	fprintf(stdout, "Image A : %s\n", name);
  switch ( imageIn->type ) {
  case UCHAR:
    theInU8 = (unsigned char ***)imageIn->array;
    break;
  case FLOAT:
    theIn32 = (float ***)imageIn->array;
    break;
  default:
    MT_ErrorParse("image Phi type not handled yet\n", 0);		   
  }
	  
  if ( imageIn->dim.z == 1 )
    flag_3D = 0;

  sprintf( name, "%s.theta.inr", prefix );
  imageTht = _VT_Inrimage( name );
  if ( imageTht == (vt_image*)NULL )
  {
     fprintf(stderr, "%s unreadable\n", name);
     sprintf( name, "%s.theta.hdr", prefix );
     imageTht = _VT_Inrimage( name );
     if ( imageTht == (vt_image*)NULL )
     {
		fprintf(stderr, "%s unreadable\n", name);
	    j=0;
        i=0;
        while (prefix[i]!='\0')
        {
          if (prefix[i]=='.') j = i;
          i++;
        }
        if (j>0)
        {
          prefix[j]='\0';
		  sprintf( name, "%s.theta.inr", prefix );
		  imageTht = _VT_Inrimage( name );
		  if ( imageTht == (vt_image*)NULL )
		  {
			fprintf(stderr, "%s unreadable\n", name);
			sprintf( name, "%s.theta.hdr", prefix );
			imageTht = _VT_Inrimage( name );
			if ( imageTht == (vt_image*)NULL )
			{
			  VT_FreeImage( imageIn );
			  fprintf(stderr, "%s unreadable\n", name);
			  MT_ErrorParse("unable to read input image B\n", 0);		   
			}
		  }
	    }
	}
  }	
  if(_VT_VERBOSE_)
	fprintf(stdout, "Image d'entree theta : %s\n", name);
  if (imageTht->type != FLOAT) {
	VT_FreeImage( imageIn );
	VT_FreeImage( imageTht );
	MT_ErrorParse("image Theta type not handled yet\n", 0);		   
    
  }
  theTht = (float ***)imageTht->array;

  if (flag_3D == 1)
  {
	sprintf( name, "%s.phi.inr", prefix );
	imagePhi = _VT_Inrimage( name );
	if ( imagePhi == (vt_image*)NULL )
	{
	  fprintf(stderr, "%s unreadable\n", name);
	  sprintf( name, "%s.phi.hdr", prefix );
	  imagePhi = _VT_Inrimage( name );
	  if ( imagePhi == (vt_image*)NULL )
	  {
		VT_FreeImage( imageIn );
		VT_FreeImage( imageTht );
		fprintf(stderr, "%s unreadable\n", name);
		MT_ErrorParse("unable to read input image C\n", 0);		   
	  }
	}
	if(_VT_VERBOSE_)
	  fprintf(stdout, "Image d'entree phi : %s\n", name);
	if (imagePhi->type != FLOAT) {
	  VT_FreeImage( imageIn );
	  VT_FreeImage( imageTht );
	  VT_FreeImage( imagePhi );
	  MT_ErrorParse("image Phi type not handled yet\n", 0);		   
	}
	thePhi = (float ***)imagePhi->array;
  }

  // Calculs
  if(par.tensor == 0)
  {
	/*--- AJOUT DE SCALAIRES ---*/

	/*--- Output initialization ---*/
	sprintf( name, "%s", par.names.out );
	t =  FLOAT;
	
	VT_Image( &imres );

	VT_InitVImage( &imres, name, 1, 
		     imageIn->dim.x, imageIn->dim.y, imageIn->dim.z, 
		     t );
		     
	imres.siz.x = imageIn->siz.x;
	imres.siz.y = imageIn->siz.y;
	imres.siz.z = imageIn->siz.z;

	if ( VT_AllocImage( &imres ) != 1 ) {
	  VT_FreeImage( imageIn );
      VT_FreeImage( imageTht );
	  if ( flag_3D != 0 )
		VT_FreeImage( imagePhi );
	  MT_ErrorParse("unable to allocate output image\n", 0);
    }

	theIm = &imres;
	dim[0]=theIm->dim.x;
	dim[1]=theIm->dim.y;
	dim[2]=theIm->dim.z;

	/*--- calculs ---*/
	if(_VT_VERBOSE_)
      fprintf(stdout, "Calculs...\n");

	if ( flag_3D == 1 ) {
	  // 3D case
	  
	  if(1 || par.alpha==0.0)	//TEMPORAIRE
	  //if( par.alpha==0.0)	//TEMPORAIRE
      for (k=0 ; k < imageIn->dim.z ; k++)
      for (j=0 ; j < imageIn->dim.y ; j++)
      for (i=0 ; i < imageIn->dim.x ; i++) 
	  {
		if ( _VT_VERBOSE_ && i==0 && j==0 && (k*100/imageIn->dim.z % 5 == 0) && ((k-1)*100/imageIn->dim.z % 5 != 0))
		  fprintf(stdout, "Slice %d\ton %d\t(%d %%)\n", k+1, (int)imageIn->dim.z, (int)(k*100/imageIn->dim.z));
		  
		switch (imageIn->type) {
		case UCHAR:
	      coef=(double) (theInU8[k][j][i]); 
	      break;
		case FLOAT:
	      coef=(double) (theIn32[k][j][i]); 
	      break;
		default:
		  VT_FreeImage( &imres );
		  VT_FreeImage( imageIn );
		  VT_FreeImage( imageTht );
		  VT_FreeImage( imagePhi );
	      MT_ErrorParse( "input image type not handled yet\n", 0 );
		}
		coef = (coef > 0.0) ? 1.0 : 0.0;
		
		if (coef <= 0.0) 
		  continue;

	    tht=theTht[k][j][i];
	    phi=thePhi[k][j][i];
	  
		pt[0]=i;
		pt[1]=j;
		pt[2]=k;
		
		if (par.alpha == 0.0) {
		  if ( addLineToBuf(theIm->buf, dim, t, pt, tht, phi, par.r) == 0 )
	      {
			VT_FreeImage( &imres );
			VT_FreeImage( imageIn );
			VT_FreeImage( imageTht );
			VT_FreeImage( imagePhi );
			MT_ErrorParse( "Error while adding a line\n", 0 );
		  }
	    }
		else {
		  if ( addConeToBuf(theIm->array, dim, t, pt, tht, phi, (double) par.r, par.alpha) == 0 ) 
	      {
			VT_FreeImage( &imres );
			VT_FreeImage( imageIn );
			VT_FreeImage( imageTht );
			VT_FreeImage( imagePhi );
			MT_ErrorParse( "Error while adding a cone\n", 0 );
		  }
		}

	  }
	  else // Conic vote with voting fields
	  {
		imagesIn[0]=imageIn;
		imagesIn[1]=imageTht;
		if (flag_3D == 1)
		  imagesIn[2]=imagePhi;
			
		if( VT_conicVote(imagesIn, &imres, (double) par.r, par.alpha, par.NangleIter ) == 0 )
		{
		  VT_FreeImage( &imres );
		  VT_FreeImage( imageIn );
		  VT_FreeImage( imageTht );
		  VT_FreeImage( imagePhi );
		  MT_ErrorParse( "Error while cone voting\n", 0 );
		}
	  }
    
	}
	else {
      // 2D case : TODO
	  VT_FreeImage( imageIn );
	  VT_FreeImage( imageTht  );
	  VT_FreeImage( &imres  );
      MT_ErrorParse( "2D case not handled yet\n", 0 );
	}
	
	
	/*--- ecriture de l'image resultat ---*/
	if ( VT_WriteInrimage( &imres ) == -1 ) {
      VT_FreeImage( &imres );
      MT_ErrorParse("unable to write output image\n", 0);
	}
	VT_FreeImage( &imres );
  }
  else
  {
	/*--- TENSEURS---*/
	vt_3Dtensor theTensor;
	mt_2Dtensor theTensor2D;

	/* Ecriture des adresses dans l'input imagesIn de la fonction de TV */
	imagesIn[0]=imageIn;
	imagesIn[1]=imageTht;
	if (flag_3D == 1)
	  imagesIn[2]=imagePhi;
  
	/*--- tensor voting ---*/
	//TV 3D
	fprintf(stdout, "Entree dans la fonction de tensor voting\n");
	if (flag_3D == 1)
	{
      //if (MT_Compute3DTensorVotingTest(&theTensor, imagesIn, 3, (double) par.r,
      //	par.alpha, LINE, par.names.out, par.writeImages) != 1)
      if ( MT_Cumul3DTensorLine(&theTensor, imagesIn, par.r, par.names.out) != 1)
      {
		MT_ErrorParse("problem while computing the tensor voting\n", 0 );
        VT_FreeImage( imageTht );
        VT_FreeImage( imagePhi );
		VT_FreeImage( imageIn );
		return(-1);
      }
	}
	else  //TODO
	{
      VT_FreeImage( imageTht );
      VT_FreeImage( imageIn );
      MT_ErrorParse("2D Tensor mode not implemented yet\n", 0);
      /*
      if (MT_Compute2DTensorVoting(&theTensor2D, imagesIn, par.scale, 0,
        4, TVCLASSIC, LINE,
        par.names.out, par.writeImages) != 1)
      {
		MT_ErrorParse("problem while computing the tensor voting\n", 0 );
        VT_FreeImage( imageTht );
		VT_FreeImage( imageIn );
		return(-1);
      }
      */
	}
	fprintf(stdout, "Retour dans la fonction principale\n");

	fprintf(stdout, "Ecriture des images...\n");
	if (flag_3D == 1)
      VT_Write3Dtensor( &theTensor);
	else
      VT_Write2DTensor( &theTensor2D );



	if (flag_3D == 1)
      VT_Free3Dtensor(&theTensor);
	else
      VT_Free2DTensor(&theTensor2D);

  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( imageIn );
  VT_FreeImage( imageTht  );
  if (flag_3D ==1)
	VT_FreeImage( imagePhi  );
  
  

  stop = clock();
  elapsed = (double)(stop-start)/CLOCKS_PER_SEC;

  if(_VT_VERBOSE_)
    fprintf(stdout, "Elapsed time : \t%.1fs\n", elapsed);

  return( 1 );
}








static void MT_Parse( int argc,
                      char *argv[],
                      local_par *par )
{
  int i, nb, status;
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
      else if ( strcmp ( argv[i], "-wi" ) == 0 ) {
        par->writeImages = 1;
      }


      /*--- traitement eventuel de l'image d'entree ---*/



      /*--- Parametres d'input/output ---*/


      /* Parametres de calcul */

      else if ( strcmp ( argv[i], "-r" ) == 0 ||
		(strcmp ( argv[i], "-r" ) == 0 ) ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -r...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->r) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -r...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-t" ) == 0 ||  
				strcmp ( argv[i], "-tensor" ) == 0 ) {

        par->tensor = 1;
      }

      else if ( strcmp ( argv[i], "-alpha" ) == 0 ||
		(strcmp ( argv[i], "-a" ) == 0 ) ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -alpha...\n", 0 );
		status = sscanf( argv[i],"%lf",&(par->alpha) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -alpha...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-nangles" ) == 0 ||
		(strcmp ( argv[i], "-n" ) == 0 ) ) {
		i += 1;
		if ( i >= argc)    MT_ErrorParse( "parsing -nangles...\n", 0 );
		status = sscanf( argv[i],"%d",&(par->NangleIter) );
		if ( status <= 0 ) MT_ErrorParse( "parsing -nangles...\n", 0 );
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
  {
    strcpy( par->names.out, ">" );  /* standart output */
  }

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
  par->r = -1;
  par->tensor = 0;
  par->writeImages = 0;
  par->alpha = 0.0;
  par->NangleIter = 3;
}


