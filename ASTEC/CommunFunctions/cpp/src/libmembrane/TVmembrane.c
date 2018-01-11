/*************************************************************************
 * TVmembrane.c -
 *
 * $Id: TVmembrane.c,v 2.0 2013/10/22 14:22:00 gael Exp $
 *
 *
 * AUTHOR:
 * Gaël Michelin (gael.michelin@inria.fr)
 * 
 * CREATION DATE: 
 * 2013/06/24
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <mt_membrane3D.h>
#include <parcelling.h>

#include <time.h>
#include <sys/time.h>

#include <chunks.h>


typedef enum {
  _UNKNOWN_ /* unknown mode */,
  _RANDOM_SAMPLING_,
  _REGULAR_SAMPLING_
} enumSampleMode;


typedef struct local_par {
  vt_names names;
  int type;
  shapeExtraction shape;
  enumTVmode TVmode;
  double scale;
  double zfact;
  int niter;   // pour sparse voting
  int nangles; // nanglesiter
  int nsticks; // nombre de sticks pour le calcul des plate fields
  int writeImages;
  hessianMode initHessian;
  int print_time;
  double sample;
  int power;
  enumSampleMode sampleMode;
  int itermax;
//  double theCoefficient;

//  double theTensorLargeCoefficient;
//  double theTensorSmallCoefficient;
//  double theTensorMultiCoefficient;
  int dimension;

} local_par;





/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void MT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static double _GetTime();
static double _GetClock();
static char *_BaseName( char *p );
static int _count_points( void *theBuf, bufferType theBufType, int *theDim );





static char *usage = "[image-in] [image-out] [-wi] [-2D]\n\
\t [-hessian [plane|line]] [-shape m|v|b|mv|mb|vb|mvb]\n\
\t [-scale %lf] [-zfact %lf] \n\
\t [-niter %d] [-nangles %d] [-nsticks %d] [-cftv|-tvclassic]\n\
\t [-sample %lf [-power]] [-random|-regular [-i %d]] \n\
\t [-parallel|-no-parallel] [-max-chunks %d]\n\
\t [-parallel-scheduling|-ps default|static|dynamic-one|dynamic|guided]\n\
\t [-inv] [-swap] [-v] [-D] [-time] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -time : ecrit le temps d'execution du code dans le fichier <time.txt>\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t options-de-type : -o 1    : unsigned char\n\
\t                   -o 2    : unsigned short int\n\
\t                   -o 2 -s : short int\n\
\t                   -o 4 -s : int\n\
\t                   -r      : float\n\
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\t -2D : slice par slice\n\
\t -membrane : detecteur de structures planes (extr. maxima + hessian init)\n\
\t -vessel : detecteur de structures lineiques (extr. maxima + hessian init)\n\
\t -scale %lf : echelle du tensor voting\n\
\t -zfact %lf : facteur de resolution selon l'axe des z\n\
\t -niter %d : nombre d'iterations de SPARSE voting\n\
\t -nangles %d : nombre d'angles pour la discretisation des champs de \
tenseurs\n\
\t -nsticks %d : nombre d'angles pour calculer le champ de tenseurs plate\n\
\t -hessian : recuperation des informations directionnelles donnees par\
la matrice hessienne calculee lors de la 1ere binarisation\
(donner le prefixe ou le binaire en argument, necessite imtheta (+ imphi en 3D)\
avec le meme prefixe que le binaire, formats supportes .inr et .hdr)\n\
\t -hessian plane|line : plane si les directions sont des normales\n\
\t aux structures (plans, defaut) ou line pour des tangentes aux structures\n\
\t (lignes)\n\
\t -shape {m,v,b} : extrait les maxima de {membranes, vaisseaux, balles}\n\
\t -sample %lf : echantillonne les votants de l'image initiale selon le \
coefficient \n\
\t -power : l'argument de -sample est alors converti en puissance pour \n\
l'echantillonnage : %age=10^-argument\n\
\t -random : sampling aleatoire\n\
\t -regular [-i %d]: sampling regulier (type Voronoi) \n\
(-i %d : nb max d'iterations)\n\
\t -cftv : champs de tenseurs de type Closed Form Tensor Voting\n\
\t -tvclassic : champs de tenseurs de type Medioni\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  int flag_3D = 1;
  local_par par;
  vt_image *imageBin;
  vt_image *imageTht=NULL;
  vt_image *imagePhi=NULL;
  vt_image* imagesIn[3];
  vt_3Dtensor theTensor;
  mt_2Dtensor theTensor2D;
  vt_image theExtrema;
  
  int t = (int) UCHAR;
  char name[256];
  char prefix[256];
  clock_t start, stop;
  double elapsed;

  start = clock();
  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;
  double time_sampling;
  double clock_sampling;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  

  /*--- lecture de l'image d'entree ---*/
  if (par.initHessian == NONE)
  {
    //imageBin=malloc(sizeof(vt_image*));
    sprintf( name, "%s", par.names.in );
    imageBin = _VT_Inrimage( name );
    if ( imageBin == (vt_image*)NULL )
    {
      //format .hdr
      int i=0;
      while (name[i]!='\0')
      {
        i++;
      }
      name[i++] = '.';
      name[i++] = 'h';
      name[i++] = 'd';
      name[i++] = 'r';
      name[i]  = '\0';

      imageBin = _VT_Inrimage( name );
      if ( imageBin == (vt_image*)NULL )
      {
        //format .inr
        name[i-3] = 'i';
        name[i-2] = 'n';
        fprintf(stdout, "i=%d\n", i);
        fprintf(stdout, "name=%s\n", name);

        imageBin = _VT_Inrimage( name );
        if ( imageBin == (vt_image*)NULL )
        {
          fprintf(stderr, "%s unreadable\n", par.names.in );
          MT_ErrorParse("unable to read input imageA\n", 0);
        }
      }
    }
    fprintf(stdout, "Image d'entree : %s\n", name);
    if ( par.dimension == 2 || imageBin->dim.z == 1 )
      flag_3D = 0;
  }
  else
  {
    sprintf( prefix, "%s", par.names.in );
    int i=0;
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
    }
    else
    {
      sprintf( name, "%s.hdr", prefix );
    }

    fprintf(stdout, "Lecture des fichiers %s...\n", prefix);
    imageBin = _VT_Inrimage( name );
    if ( imageBin == (vt_image*)NULL )
    {
      int i=0;
      while (name[i]!='\0')
      {
        i++;
      }
      name[i-3] = 'i';
      name[i-2] = 'n';

      imageBin = _VT_Inrimage( name );
      if ( imageBin == (vt_image*)NULL )
      {
        fprintf(stderr, "%s unreadable\n", name);
        MT_ErrorParse("unable to read input imageB\n", 0);
      }
    }
    fprintf(stdout, "Image d'entree binaire : %s\n", name);

    if ( par.dimension == 2 || imageBin->dim.z == 1 )
      flag_3D = 0;


    sprintf( name, "%s.theta.hdr", prefix );
    imageTht = _VT_Inrimage( name );
    if ( imageTht == (vt_image*)NULL )
    {
      int i=0;
      while (name[i]!='\0')
      {
        i++;
      }
      name[i-3] = 'i';
      name[i-2] = 'n';
      imageTht = _VT_Inrimage( name );
      if ( imageTht == (vt_image*)NULL )
      {
          int j=0;
          i=0;
          while (prefix[i]!='\0')
          {
            if (prefix[i]=='.') j = i;
            i++;
          }
          if (j>0)
          {
            prefix[j]='\0';
            sprintf( name, "%s.theta.hdr", prefix );

            imageTht = _VT_Inrimage( name );
            if ( imageTht == (vt_image*)NULL )
            {
              i=0;
              while (name[i]!='\0')
              {
                i++;
              }
              name[i-3] = 'i';
              name[i-2] = 'n';
              imageTht = _VT_Inrimage( name );
              if ( imageTht == (vt_image*)NULL )
              {
                  VT_FreeImage( imageBin );
                  fprintf(stderr, "%s unreadable\n", name);
                  MT_ErrorParse("unable to read input imageC\n", 0);
              }
            }
          }

      }
    }
    fprintf(stdout, "Image d'entree theta : %s\n", name);

    if (flag_3D == 1)
    {

      sprintf( name, "%s.phi.hdr", prefix );
      imagePhi = _VT_Inrimage( name );
      if ( imagePhi == (vt_image*)NULL )
      {
        int i=0;
        while (name[i]!='\0')
        {
          i++;
        }
        name[i-3] = 'i';
        name[i-2] = 'n';
        imagePhi = _VT_Inrimage( name );
        if ( imagePhi == (vt_image*)NULL )
        {
          VT_FreeImage( imageTht );
          VT_FreeImage( imageBin );
          fprintf(stderr, "%s unreadable\n", name);
          MT_ErrorParse("unable to read input imageD\n", 0);
        }
      }
      fprintf(stdout, "Image d'entree phi : %s\n", name);
    }
  }
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( imageBin );
  if ( par.names.swap == 1 ) VT_SwapImage( imageBin );



  /* Echantillonnage de l'image */
  if (par.power == 1)
  {
    par.sample = exp((par.sample>=0 ? - par.sample : par.sample)*log(10));
  }

  if ( par.sample != 1.0 )
  {
    time_sampling = _GetTime();
    clock_sampling = _GetClock();
    if ( _VT_VERBOSE_ )
      fprintf(stdout, "Sampling step : sampling coefficient = %f\n", par.sample);
    if (par.sample > 0 && par.sample < 1)
    {
      if (flag_3D == 1)
      {
        if (par.sampleMode == _RANDOM_SAMPLING_)
        {
          if (MT_SampleBin(imageBin, par.sample) != 1 )
          {
            if (par.initHessian != NONE)
            {
              VT_FreeImage( imageTht );
              VT_FreeImage( imagePhi );
            }
            VT_FreeImage( imageBin );
            MT_ErrorParse("unexpected error while sampling binary image\n", 0);
          }
        }
        if (par.sampleMode == _REGULAR_SAMPLING_)
        {
          int nparcels;
          int **theSeeds = NULL;
          int *seeds = NULL;
          int theDim[3];
          int i,j,k,n;

          theDim[0] = imageBin->dim.x;
          theDim[1] = imageBin->dim.y;
          theDim[2] = imageBin->dim.z;

          nparcels =(int) (((double)_count_points(imageBin->buf, imageBin->type, theDim)) * par.sample);

          //fprintf(stdout, "nparcels=%d\n", nparcels);

          seeds = (int*)malloc( 3*nparcels * sizeof( int ) );
          theSeeds = (int**)malloc( nparcels * sizeof( int* ) );
          if ( seeds == NULL || theSeeds == NULL ) {
            if ( seeds != NULL ) free( seeds );
            if ( theSeeds != NULL ) free( theSeeds );
            VT_FreeImage( imageBin );
            if (par.initHessian != NONE)
            {
              VT_FreeImage( imageTht );
              if ( flag_3D == 1)
				VT_FreeImage( imagePhi );
            }
            MT_ErrorParse("error when allocating seeds arrays\n", 0);
          }

          for ( n=0;n<nparcels;n++ ) {
            theSeeds[n] = &(seeds[3*n]);
          }

          if ( parcelling( imageBin->buf, imageBin->type,
                           theSeeds, nparcels,
                           NULL, TYPE_UNKNOWN,
                           NULL, TYPE_UNKNOWN,
                           theDim, 0 ) != 1 ) {
            VT_FreeImage( imageBin );
            if (par.initHessian != NONE)
            {
              VT_FreeImage( imageTht );
              if (flag_3D==1)
				VT_FreeImage( imagePhi );
            }
            if ( seeds != NULL ) free( seeds );
            if ( theSeeds != NULL ) free( theSeeds );
            MT_ErrorParse("error when processing\n", 0);
          }
          switch(imageBin->type)
          {
          case UCHAR :
            {
              u8 ***array = (unsigned char ***)imageBin->array;
              for (i=0;i<theDim[0];i++)
              for (j=0;j<theDim[1];j++)
              for (k=0;k<theDim[2];k++)
              {
                array[k][j][i]= (unsigned char) 0 ;
              }
              for (n=0;n<nparcels;n++)
              {
                array[theSeeds[n][2]][theSeeds[n][1]][theSeeds[n][0]]=(unsigned char) 255;
              }
              break;
            }
          case FLOAT :
            {
              r32 ***array = (float ***)imageBin->array;
              for (n=0;n<nparcels;n++)
              {
                array[theSeeds[n][2]][theSeeds[n][1]][theSeeds[n][0]]*= -1;
              }
              for (i=0;i<theDim[0];i++)
              for (j=0;j<theDim[1];j++)
              for (k=0;k<theDim[2];k++)
              {
                array[k][j][i]= (array[k][j][i]>0) ? 0.0 : -array[k][j][i];
              }
              break;
            }
          default:
            VT_FreeImage( imageBin );
            if (par.initHessian != NONE)
            {
              VT_FreeImage( imageTht );
              if (flag_3D==1)
				VT_FreeImage( imagePhi );
            }
            if ( seeds != NULL ) free( seeds );
            if ( theSeeds != NULL ) free( theSeeds );
            MT_ErrorParse("error: image extension not implemented yet\n", 0);
          }
          if ( seeds != NULL ) free( seeds );
          if ( theSeeds != NULL ) free( theSeeds );
        }
      }
      else
      {
        if (MT_SampleBin2D(imageBin, par.sample) != 1 )
        {
          if (par.initHessian != NONE)
          {
            VT_FreeImage( imageTht );
          }
          VT_FreeImage( imageBin );
          MT_ErrorParse("unexpected error while sampling binary image\n", 0);
        }
      }
      if (1 || par.writeImages == 1)
      {
        char tmp[256];
        sprintf(tmp, "%s.sample.inr", par.names.out);
        fprintf(stdout, "Ecriture de %s ...\n", tmp);
        VT_WriteInrimageWithName(imageBin, tmp);
      }
    }
    else
    {
      if (par.initHessian != NONE)
      {
        VT_FreeImage( imageTht );
        if (flag_3D == 1)
          VT_FreeImage( imagePhi );
      }
      VT_FreeImage( imageBin );
      MT_ErrorParse("incorrect sampling parameter : 0 < expected value <= 1\n", 0);
    }

    time_exit = _GetTime();
    clock_exit = _GetClock();
    if ( par.print_time ) {
      char Ftime[256];
      sprintf(Ftime, "%s.timesampling.txt", par.names.out);
      fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName(program), time_exit - time_sampling );
      fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_sampling );
      fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_sampling)/(time_exit - time_sampling) );
      FILE *fichier = fopen (Ftime, "a" );
      if (fichier == NULL)
      {
        perror (Ftime);
      }
      else
      {
        int i;

        fprintf(fichier, "\n");
        for(i=0;i<argc;i++)
          fprintf( fichier, "%s ", argv[i]);
        fprintf(fichier, "\n");
        fprintf(fichier, "%s: elapsed (real) time = %f\n", _BaseName(program), time_exit - time_sampling );
        fprintf(fichier, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_sampling );
        fprintf(fichier, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_sampling)/(time_exit - time_sampling) );

        fclose(fichier);
      }
    }


  }


  /* Ecriture des adresses dans l'input imagesIn de la fonction de TV */
  imagesIn[0]=imageBin;
  imagesIn[1]=imageTht;
  if (flag_3D == 1)
    imagesIn[2]=imagePhi;
  
  /*--- tensor voting ---*/
  //TV 3D
  fprintf(stdout, "Entree dans la fonction de tensor voting\n");
  if (flag_3D == 1)
  {
    if (MT_Compute3DTensorVoting(&theTensor, imagesIn, par.scale, par.zfact, par.niter,
      par.nangles,par.nsticks, par.TVmode, par.initHessian,
      par.names.out, par.writeImages) != 1)
    {
      MT_ErrorParse("problem while computing the tensor voting\n", 0 );
      if (par.initHessian != NONE)
      {
        VT_FreeImage( imageTht );
        VT_FreeImage( imagePhi );
      }
      VT_FreeImage( imageBin );
      return(-1);
    }
  }
  else
  {
    if (MT_Compute2DTensorVoting(&theTensor2D, imagesIn, par.scale, par.niter,
        par.nangles, par.TVmode, par.initHessian,
        par.names.out, par.writeImages) != 1)
    {
      MT_ErrorParse("problem while computing the tensor voting\n", 0 );
      if (par.initHessian != NONE)
      {
        VT_FreeImage( imageTht );
      }
      VT_FreeImage( imageBin );
      return(-1);
    }
  }
  fprintf(stdout, "Retour dans la fonction principale\n");

  fprintf(stdout, "Ecriture des images...\n");
  if (flag_3D == 1)
    VT_Write3Dtensor( &theTensor );
  else
    VT_Write2DTensor( &theTensor2D );
  
  //Surfaces extraction
  if(par.shape == MEMBRANE || par.shape == MV || par.shape == MB || 
	 par.shape == MVB ) {
	  t = (int) FLOAT;
	  if (par.shape == MEMBRANE)
		sprintf( name, "%s.ext.inr", par.names.out );
	  else
		sprintf( name, "%s.plane.inr", par.names.out );
	  VT_InitImage( &theExtrema, name, imageBin->dim.x, imageBin->dim.y,
				  imageBin->dim.z, t );
	  if ( VT_AllocImage( &(theExtrema) ) != 1 ) {
		MT_ErrorParse("problem while allocating theExtrema\n", 0 );
		if (par.initHessian != NONE)
		{
		  VT_FreeImage( imageTht );
		  if (flag_3D ==1)
			VT_FreeImage( imagePhi );
		}
		VT_FreeImage( imageBin );
		if (flag_3D == 1)
		  VT_Free3Dtensor(&theTensor);
		else
		  VT_Free2DTensor(&theTensor2D);
		return( -1 );
	  }

	  fprintf(stdout, "Extraction des plans centraux issus du tensor voting...\n");
	  if (flag_3D == 1)
		MT_ComputeTensorSurfaceExtrema( &theTensor, &theExtrema, par.zfact );
	  else
		MT_Compute2DTensorLineicExtrema( &theTensor2D, &theExtrema);

	  fprintf(stdout, "Ecriture du resultat dans %s...\n", name);
	  VT_WriteInrimage( &theExtrema);
	  VT_FreeImage(&theExtrema);
  }

  //Lines extraction
  if(par.shape == VESSEL || par.shape == MV || par.shape == VB || 
	 par.shape == MVB ) {
	  t = (int) FLOAT;
	  sprintf( name, "%s.line.inr", par.names.out );
	  VT_InitImage( &theExtrema, name, imageBin->dim.x, imageBin->dim.y,
				  imageBin->dim.z, t );
	  if ( VT_AllocImage( &(theExtrema) ) != 1 ) {
		MT_ErrorParse("problem while allocating theExtrema\n", 0 );
		if (par.initHessian != NONE)
		{
		  VT_FreeImage( imageTht );
		  if (flag_3D ==1)
			VT_FreeImage( imagePhi );
		}
		VT_FreeImage( imageBin );
		if (flag_3D == 1)
		  VT_Free3Dtensor(&theTensor);
		else
		  VT_Free2DTensor(&theTensor2D);
		return( -1 );
	  }

	  fprintf(stdout, "Extraction des plans centraux issus du tensor voting...\n");
	  if (flag_3D == 1)
		MT_ComputeTensorLineExtrema( &theTensor, &theExtrema, par.zfact );
	  else
		MT_Compute2DTensorLineicExtrema( &theTensor2D, &theExtrema);

	  fprintf(stdout, "Ecriture du resultat dans %s...\n", name);
	  VT_WriteInrimage( &theExtrema);
	  VT_FreeImage(&theExtrema);
  }


  //Balls extraction
  if(par.shape == BALL || par.shape == MB || par.shape == VB || 
	 par.shape == MVB ) {
	  t = (int) FLOAT;
	  sprintf( name, "%s.ball.inr", par.names.out );
	  VT_InitImage( &theExtrema, name, imageBin->dim.x, imageBin->dim.y,
				  imageBin->dim.z, t );
	  if ( VT_AllocImage( &(theExtrema) ) != 1 ) {
		MT_ErrorParse("problem while allocating theExtrema\n", 0 );
		if (par.initHessian != NONE)
		{
		  VT_FreeImage( imageTht );
		  if (flag_3D ==1)
			VT_FreeImage( imagePhi );
		}
		VT_FreeImage( imageBin );
		if (flag_3D == 1)
		  VT_Free3Dtensor(&theTensor);
		else
		  VT_Free2DTensor(&theTensor2D);
		return( -1 );
	  }

	  fprintf(stdout, "Extraction des plans centraux issus du tensor voting...\n");
	  if (flag_3D == 1)
		MT_ComputeTensorBallExtrema( &theTensor, &theExtrema, par.zfact );
	  else
		MT_Compute2DTensorBallExtrema( &theTensor2D, &theExtrema);

	  fprintf(stdout, "Ecriture du resultat dans %s...\n", name);
	  VT_WriteInrimage( &theExtrema);
	  VT_FreeImage(&theExtrema);
  }




  /*--- liberations memoires ---*/
  VT_FreeImage( imageBin );
  if (par.initHessian != NONE)
  {
    VT_FreeImage( imageTht );
    if (flag_3D ==1)
      VT_FreeImage( imagePhi );
  }
  if (flag_3D == 1)
    VT_Free3Dtensor(&theTensor);
  else
    VT_Free2DTensor(&theTensor2D);

  fprintf(stdout, "Fin des operations\n");

  stop = clock();
  elapsed = (double)(stop-start)/CLOCKS_PER_SEC;
  fprintf(stdout, "Elapsed time : \t%.1fs\n", elapsed);

  time_exit = _GetTime();
  clock_exit = _GetClock();
  if ( par.print_time ) {
    char Ftime[256];
    sprintf(Ftime, "%s.time.txt", par.names.out);
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName(program), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
    FILE *fichier = fopen (Ftime, "a" );
    if (fichier == NULL)
    {
      perror (Ftime);
    }
    else
    {
      int i;

      fprintf(fichier, "\n");
      for(i=0;i<argc;i++)
        fprintf( fichier, "%s ", argv[i]);
      fprintf(fichier, "\n");
      fprintf(fichier, "%s: elapsed (real) time = %f\n", _BaseName(program), time_exit - time_init );
      fprintf(fichier, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
      fprintf(fichier, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );



      fclose(fichier);
    }
  }


  return( 1 );
}








static void VT_Parse( int argc, 
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
      }
      else if ( strcmp ( argv[i], "-time" ) == 0 ) {
        par->print_time = 1;
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


      /* TV params */
      else if ( strcmp ( argv[i], "-scale" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -scale...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->scale) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -scale...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-zfact" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -zfact...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->zfact) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -zfact...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-niter" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -niter...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->niter) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -niter...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-nangles" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -nangles...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->nangles) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -nangles...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-nsticks" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -nsticks...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->nsticks) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -nsticks...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-sample" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -sample...\n", 0 );
        status = sscanf( argv[i],"%lf",&(par->sample) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -sample...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-power" ) == 0 ) {
        par->power = 1;
      }

      else if ( strcmp ( argv[i], "-random" ) == 0 ) {
        par->sampleMode = _RANDOM_SAMPLING_;
      }

      else if ( strcmp ( argv[i], "-regular" ) == 0 ) {
        par->sampleMode = _REGULAR_SAMPLING_;
      }

      else if ( strcmp ( argv[i], "-iterations" ) == 0 || strcmp ( argv[i], "-i" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -iterations...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->itermax) );
        if ( status <= 0 ) MT_ErrorParse( "parsing -iterations...\n", 0 );
        parcelling_setNumberOfIterations( par->itermax );
      }


      else if ( strcmp ( argv[i], "-cftv" ) == 0 ) {
        par->TVmode = CFTV;
      }

      else if ( strcmp ( argv[i], "-tvclassic" ) == 0 ) {
        par->TVmode = TVCLASSIC;
      }


      else if ( strcmp ( argv[i], "-hessian" ) == 0 ) {
        par->initHessian = NONE;
        i += 1;
        if ( i < argc)    {
          if ( strcmp(argv[i],"plane") == 0 )
			par->initHessian = PLANE;
          if ( strcmp(argv[i],"line") == 0 )
			par->initHessian = LINE;
	    }
	    if (par->initHessian == NONE) {
		  par->initHessian = PLANE;
		  i -= 1;
	    }

      }

      else if ( strcmp ( argv[i], "-shape" ) == 0 ) {
        i += 1;
        if ( i >= argc)    MT_ErrorParse( "parsing -shape...\n", 0 );
        if ( strcmp(argv[i],"m") == 0 )
		  par->shape = MEMBRANE;
        else if ( strcmp(argv[i],"v") == 0 )
		  par->shape = VESSEL;
        else if ( strcmp(argv[i],"b") == 0 )
		  par->shape = BALL;
        else if ( strcmp(argv[i],"mv") == 0 )
		  par->shape = MV;
        else if ( strcmp(argv[i],"mb") == 0 )
		  par->shape = MB;
        else if ( strcmp(argv[i],"vb") == 0 )
		  par->shape = VB;
        else if ( strcmp(argv[i],"mvb") == 0 )
		  par->shape = MVB;
        else
          MT_ErrorParse( "parsing -shape...\n", 0 );

      }

      else if ( strcmp ( argv[i], "-wi" ) == 0 ) {
        par->writeImages = 1;
      }

      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
        par->dimension = 2;
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








static void VT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );
  par->type = TYPE_UNKNOWN;

  par->writeImages = 0;
  par->initHessian = NONE;

  par->TVmode = 0;
  par->nangles = 1;
  par->niter = 1;
  par->scale = 1;
  par->zfact = 1;
  par->nsticks = 36;
  par->print_time = 0;
  par->dimension = 3;
  par->sample = 1.0;
  par->power=0;
  par->sampleMode = _RANDOM_SAMPLING_;
  par->itermax = 0;
  par->shape = MEMBRANE;
}


static double _GetTime()
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}

static double _GetClock()
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}

static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}

static int _count_points( void *theBuf, bufferType theBufType, int *theDim )
{
  char *proc = "_count_points";
  int i, v = theDim[0]*theDim[1]*theDim[2];
  int n = 0;

  switch ( theBufType ) {

  default :
    if ( _VT_VERBOSE_ ) {
      fprintf( stderr, "%s: such label type not handled yet\n", proc );
    }
    return( -1 );

  case UCHAR :
    {
      u8 *buf = (u8*)theBuf;
      for ( i=0; i<v; i++ )
        if ( buf[i] ) n++;
    }
    break;

  case USHORT :
    {
      u16 *buf = (u16*)theBuf;
      for ( i=0; i<v; i++ )
        if ( buf[i] ) n++;
    }
    break;

  case FLOAT :
    {
      r32 *buf = (r32*)theBuf;
      for ( i=0; i<v; i++ )
        if ( buf[i] != 0 ) n++;
    }
    break;

  }

  return( n );
}
