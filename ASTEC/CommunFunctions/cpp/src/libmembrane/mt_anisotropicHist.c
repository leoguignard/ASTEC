/*
 * mt_anisotropicHist.c -
 *
 * $Id: mt_anisotropicHist.c,v 1.0 2013/10/22 17:08:34 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2013/10/22
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <mt_anisotropicHist.h>
#include <transfo.h>

static int _verbose_ = 1;

int MT_ComputeRepAnisotropic(vt_Rep_Angles imageIn,
                                    vt_image *imRepAnisotropes)
{
    char *proc = "MT_ComputeRepAnisotropic";
    int i,j,k;
    double n[3];

    float rep;
    float ***theTht, ***thePhi;
    float ***theRep;

    vt_image *repX;
    vt_image *repY;
    vt_image *repZ;
    float ***theRepX, ***theRepY, ***theRepZ;

    char name[STRINGLENGTH], prefix[STRINGLENGTH];


    if(1)   // Init imResAnisotropes[]
    {
      sprintf(prefix, "%s", imageIn.imrep->name);
      i=0;
      while(prefix[i] != '\0')
        i++;
      while(prefix[i] != '.')
        i--;
      prefix[i]='\0';
      sprintf(name,"%s.X.inr",prefix);
      VT_InitFromImage( imRepAnisotropes, (imageIn.imrep), name,  (int)FLOAT );
      sprintf(name,"%s.Y.inr",prefix);
      VT_InitFromImage( imRepAnisotropes+1, (imageIn.imrep), name,  (int)FLOAT );
      sprintf(name,"%s.Z.inr",prefix);
      VT_InitFromImage( imRepAnisotropes+2, (imageIn.imrep), name,  (int)FLOAT );

      if ( VT_AllocImage( imRepAnisotropes ) != 1 ) {
          fprintf(stderr, "%s: Unable to allocate 1st anisotropic response image\n", proc);
        return( 0 );
      }

      repX=imRepAnisotropes;

      if ( VT_AllocImage( imRepAnisotropes+1 ) != 1 ) {
          fprintf(stderr, "%s: Unable to allocate 2nd anisotropic response image\n", proc);
        VT_FreeImage( repX );
        return( 0 );
      }

      repY=imRepAnisotropes+1;

      if ( VT_AllocImage( imRepAnisotropes+2 ) != 1 ) {
          fprintf(stderr, "%s: Unable to allocate 3rd anisotropic response image\n", proc);
        VT_FreeImage( repX );
        VT_FreeImage( repY );
        return( 0 );
      }

      repZ=imRepAnisotropes+2;

      theRepX=(float***)repX->array;
      theRepY=(float***)repY->array;
      theRepZ=(float***)repZ->array;
    }

    theTht=(float***)imageIn.imtheta->array;
    thePhi=(float***)imageIn.imphi->array;
    theRep=(float***)imageIn.imrep->array;

    for(i=0;i<imageIn.imrep->dim.z;i++)
    for(j=0;j<imageIn.imrep->dim.y;j++)
    for(k=0;k<imageIn.imrep->dim.x;k++)
    {

      SphericalAnglesToUnitVector( (double)(theTht[i][j][k]), (double)(thePhi[i][j][k]), n);
      rep=theRep[i][j][k];

      theRepX[i][j][k] = rep*n[0]*n[0];
      theRepY[i][j][k] = rep*n[1]*n[1];
      theRepZ[i][j][k] = rep*n[2]*n[2];
    }


    return (1);
}

int MT_ComputeHisto3D( vt_Rep_Angles imageIn,
                             vt_anisotropicHisto *histo, vt_histmode mode,
                             vt_histparam inithistparam, int normalize, 
                             double lmin, int initlmin, 
                             double lmax, int initlmax, 
                             int wi,
                             vt_image *imResAnisotropes)
{
  char *proc = "MT_ComputeHisto3D";

  int i,j,k,ind,memX,memY,memZ;
  double n[3];
  float rep;
  float ***theTht, ***thePhi;
  float ***theRep;

  double *bufX;
  double *bufY;
  double *bufZ;

  double *bufnx;
  double *bufny;
  double *bufnz;

  vt_image *repX;
  vt_image *repY;
  vt_image *repZ;
  float ***theRepX, ***theRepY, ***theRepZ;
  char name[STRINGLENGTH], prefix[STRINGLENGTH];

  double minX=0, minY=0, minZ=0, maxX=0, maxY=0, maxZ=0;
  vt_histDescriptor *descX, *descY, *descZ;
  descX=&(histo->descriptorX);
  descY=&(histo->descriptorY);
  descZ=&(histo->descriptorZ);

  int dimBuf=imageIn.imrep->dim.x*imageIn.imrep->dim.y*imageIn.imrep->dim.z;
  bufX=malloc(dimBuf*sizeof(double));
  bufY=malloc(dimBuf*sizeof(double));
  bufZ=malloc(dimBuf*sizeof(double));

  if((mode == _HIST_N_) || (mode == _HIST_N_2_))
  {
    bufnx=malloc(dimBuf*sizeof(double));
    bufny=malloc(dimBuf*sizeof(double));
    bufnz=malloc(dimBuf*sizeof(double));
  }

  if(wi == 1 )
  {
    sprintf(prefix, "%s", imageIn.imrep->name);
    i=0;
    while(prefix[i] != '\0')
      i++;
    while(prefix[i] != '.')
      i--;
    prefix[i]='\0';
    sprintf(name,"%s.X.inr",prefix);
    VT_InitFromImage( imResAnisotropes, (imageIn.imrep), name,  (int)FLOAT );
    sprintf(name,"%s.Y.inr",prefix);
    VT_InitFromImage( imResAnisotropes+1, (imageIn.imrep), name,  (int)FLOAT );
    sprintf(name,"%s.Z.inr",prefix);
    VT_InitFromImage( imResAnisotropes+2, (imageIn.imrep), name,  (int)FLOAT );

    if ( VT_AllocImage( imResAnisotropes ) != 1 ) {
      free(bufX); bufX=NULL;
      free(bufY); bufY=NULL;
      free(bufZ); bufZ=NULL;
      if((mode == _HIST_N_) || (mode == _HIST_N_2_))
      {
        free(bufnx); bufnx=NULL;
        free(bufny); bufny=NULL;
        free(bufnz); bufnz=NULL;
      }
      return( -1 );
    }
    repX=imResAnisotropes;
    if ( VT_AllocImage( imResAnisotropes+1 ) != 1 ) {
      VT_FreeImage( repX );
      free(bufX); bufX=NULL;
      free(bufY); bufY=NULL;
      free(bufZ); bufZ=NULL;
      if((mode == _HIST_N_) || (mode == _HIST_N_2_))
      {
        free(bufnx); bufnx=NULL;
        free(bufny); bufny=NULL;
        free(bufnz); bufnz=NULL;
      }
      return( -1 );
    }
    repY=imResAnisotropes+1;

    if ( VT_AllocImage( imResAnisotropes+2 ) != 1 ) {
      VT_FreeImage( repX );
      VT_FreeImage( repY );
      free(bufX); bufX=NULL;
      free(bufY); bufY=NULL;
      free(bufZ); bufZ=NULL;
      if((mode == _HIST_N_) || (mode == _HIST_N_2_))
      {
        free(bufnx); bufnx=NULL;
        free(bufny); bufny=NULL;
        free(bufnz); bufnz=NULL;
      }
      return( -1 );
    }
    repZ=imResAnisotropes+2;

    theRepX=(float***)repX->array;
    theRepY=(float***)repY->array;
    theRepZ=(float***)repZ->array;
  }

  theTht=(float***)imageIn.imtheta->array;
  thePhi=(float***)imageIn.imphi->array;
  theRep=(float***)imageIn.imrep->array;

  memX=memY=memZ=0;
  for(i=0;i<imageIn.imrep->dim.z;i++)
  for(j=0;j<imageIn.imrep->dim.y;j++)
  for(k=0;k<imageIn.imrep->dim.x;k++)
  {
    ind=k+j*imageIn.imrep->dim.x+i*imageIn.imrep->dim.x*imageIn.imrep->dim.y;
    SphericalAnglesToUnitVector( (double)(theTht[i][j][k]), (double)(thePhi[i][j][k]), n);
    rep=theRep[i][j][k];
    switch (mode) {
    case _HIST_PROJ_ :
      bufX[ind]=((double)rep)*fabs((double)n[0]);//(n[0]>=0 ? n[0] : -n[0]);
      bufY[ind]=((double)rep)*fabs((double)n[1]);//(n[1]>=0 ? n[1] : -n[1]);
      bufZ[ind]=((double)rep)*fabs((double)n[2]);//(n[2]>=0 ? n[2] : -n[2]);
      break;
    case _HIST_PROJ_2_ :
      bufX[ind] = bufY[ind] = bufZ[ind] = 0;
      if ((fabs(n[0]) >= fabs(n[1])) && (fabs(n[0]) >= fabs(n[2])))
        bufX[ind] = rep;
      if ((fabs(n[1]) >= fabs(n[0])) && (fabs(n[1]) >= fabs(n[2])))
        bufY[ind] = rep;
      if ((fabs(n[2]) >= fabs(n[0])) && (fabs(n[2]) >= fabs(n[1])))
        bufZ[ind] = rep;
      break;
    case _HIST_N_ :
      bufX[ind] = bufY[ind] = bufZ[ind] = rep;
      bufnx[ind]= fabs(n[0]);
      bufny[ind]= fabs(n[1]);
      bufnz[ind]= fabs(n[2]);
      break;
    case _HIST_N_2_ :
      bufX[ind] = bufY[ind] = bufZ[ind] = rep;
      bufnx[ind]= n[0]*n[0];
      bufny[ind]= n[1]*n[1];
      bufnz[ind]= n[2]*n[2];
      break;
    }
    if(wi == 1)
    {
      theRepX[i][j][k] = (float) (bufX[ind]);
      theRepY[i][j][k] = (float) (bufY[ind]);
      theRepZ[i][j][k] = (float) (bufZ[ind]);
    }

    // MAJ des minima et maxima de réponses X, Y, Z
    if(bufX[ind]>0)
    {
      if (memX==0)
      {
        memX=1;
        minX=bufX[ind];
        maxX=bufX[ind];
      }
      else
      {
        if(bufX[ind]<minX)   minX=bufX[ind];
        if(bufX[ind]>maxX)   maxX=bufX[ind];
      }
    }
    if(bufY[ind]>0)
    {
      if (memY==0)
      {
        memY=1;
        minY=bufY[ind];
        maxY=bufY[ind];
      }
      else
      {
        if(bufY[ind]<minY)   minY=bufY[ind];
        if(bufY[ind]>maxY)   maxY=bufY[ind];
      }
    }
    if(bufZ[ind]>0)
    {
      if (memZ==0)
      {
        memZ=1;
        minZ=bufZ[ind];
        maxZ=bufZ[ind];
      }
      else
      {
        if(bufZ[ind]<minZ)   minZ=bufZ[ind];
        if(bufZ[ind]>maxZ)   maxZ=bufZ[ind];
      }
    }
  }

  MT_InitDescriptor(descX, minX, maxX, dimBuf, inithistparam);
  MT_InitDescriptor(descY, minY, maxY, dimBuf, inithistparam);
  MT_InitDescriptor(descZ, minZ, maxZ, dimBuf, inithistparam);


  if(VT_AllocAnisotropicHisto( histo) != 1)
  {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate anisotropic histogram\n", proc );
    free(bufX); bufX=NULL;
    free(bufY); bufY=NULL;
    free(bufZ); bufZ=NULL;
    if((mode == _HIST_N_) || (mode == _HIST_N_2_))
    {
      free(bufnx); bufnx=NULL;
      free(bufny); bufny=NULL;
      free(bufnz); bufnz=NULL;
    }
    if(wi==1)
    {
      VT_FreeImage( repX );
      VT_FreeImage( repY );
      VT_FreeImage( repZ );
    }
    return( -1 );
  }

  histo->descriptorX.total=MT_ComputeHistoN(bufX, bufnx, dimBuf, mode, histo->histX, *descX);
  histo->descriptorY.total=MT_ComputeHistoN(bufY, bufny, dimBuf, mode, histo->histY, *descY);
  histo->descriptorZ.total=MT_ComputeHistoN(bufZ, bufnz, dimBuf, mode, histo->histZ, *descZ);

  free(bufX); bufX=NULL;
  free(bufY); bufY=NULL;
  free(bufZ); bufZ=NULL;
  if((mode == _HIST_N_) || (mode == _HIST_N_2_))
  {
    free(bufnx); bufnx=NULL;
    free(bufny); bufny=NULL;
    free(bufnz); bufnz=NULL;
  }


  if(normalize == 1)
  {
	int i_min, i_max, i_tot, inc;
	double *H;
	double x, delta, Hmax;
//	int axe;
	
	/*--- X axis ---*/
	H=(double*) histo->histX;
	if(initlmin==0)
	  lmin=descX->lower;
	if(initlmax==0)
	  lmax=descX->upper;
	i_tot=descX->ncell;
	i_min=0;
	delta=descX->precision;
	x=descX->lower+0.5*delta;
	while (i_min<i_tot && x<lmin)
	{
	  x+=delta;
	  i_min++;
	}
	if (i_min>=i_tot) {
	  fprintf(stderr, "%s: Error while normalizing X histogram", proc);
	  return(0);
	}
	i_max=i_min;
	while (i_max<i_tot && x<lmax)
	{
	  x+=delta;
	  i_max++;
	}
	Hmax=H[i_min];
	for (inc=i_min+1 ; inc<i_max ; inc++)
	{
	  if (Hmax<H[inc])
	    Hmax=H[inc];
	}
	if (Hmax>0)
	  for (inc=0 ; inc<i_tot ; inc++) 
		H[inc]/=Hmax;
	
	/*--- Y axis ---*/
	H=(double*) histo->histY;
	if(initlmin==0)
	  lmin=descY->lower;
	if(initlmax==0)
	  lmax=descY->upper;
	i_tot=descY->ncell;
	i_min=0;
	delta=descY->precision;
	x=descY->lower+0.5*delta;
	while (i_min<i_tot && x<lmin)
	{
	  x+=delta;
	  i_min++;
	}
	if (i_min>=i_tot) {
	  fprintf(stderr, "%s: Error while normalizing Y histogram", proc);
	  return(0);
	}
	i_max=i_min;
	while (i_max<i_tot && x<lmax)
	{
	  x+=delta;
	  i_max++;
	}
	Hmax=H[i_min];
	for (inc=i_min+1 ; inc<i_max ; inc++)
	{
	  if (Hmax<H[inc])
	    Hmax=H[inc];
	}
	if (Hmax>0)
	  for (inc=0 ; inc<i_tot ; inc++) 
		H[inc]/=Hmax;
	
	/*--- Z axis ---*/
	H=(double*) histo->histZ;
	if(initlmin==0)
	  lmin=descZ->lower;
	if(initlmax==0)
	  lmax=descZ->upper;
	i_tot=descZ->ncell;
	i_min=0;
	delta=descZ->precision;
	x=descZ->lower+0.5*delta;
	while (i_min<i_tot && x<lmin)
	{
	  x+=delta;
	  i_min++;
	}
	if (i_min>=i_tot) {
	  fprintf(stderr, "%s: Error while normalizing Z histogram", proc);
	  return(0);
	}
	i_max=i_min;
	while (i_max<i_tot && x<lmax)
	{
	  x+=delta;
	  i_max++;
	}
	Hmax=H[i_min];
	for (inc=i_min+1 ; inc<i_max ; inc++)
	{
	  if (Hmax<H[inc])
	    Hmax=H[inc];
	}
	if (Hmax>0)
	  for (inc=0 ; inc<i_tot ; inc++) 
		H[inc]/=Hmax;
  }

  return(1);
}

double MT_ComputeHistoN(double *buf, double *bufn, int dim, vt_histmode mode,
                            double *hist, vt_histDescriptor desc)
{
  double upper = desc.upper, lower = desc.lower;
  int ncell=desc.ncell;
  int index;
  int i;
  double sum=0;

  for(i=0;i<dim;i++)
  {
    index=(int)((buf[i]-lower)/(upper-lower)*ncell);
    if(index<0)
    {
      continue;
    }
    if(index>=ncell)
    {
      if(buf[i]==upper)
      {
        switch (mode) {
        default :
        case _HIST_PROJ_ :
        case _HIST_PROJ_2_ :
          sum++;
          hist[ncell-1]++;
          break;
        case _HIST_N_ :
        case _HIST_N_2_ :
          sum+=bufn[i];
          hist[ncell-1]+=bufn[i];
          break;
        }
      }
      continue;
    }
    switch (mode) {
    default :
    case _HIST_PROJ_ :
    case _HIST_PROJ_2_ :
      sum++;
      hist[index]++;
      break;
    case _HIST_N_ :
    case _HIST_N_2_ :
      sum+=bufn[i];
      hist[index]+=bufn[i];
      break;
    }
  }

  return(sum);
}

void MT_InitHistParam ( vt_histparam *histparam)
{
  histparam->ncell = 0;
  histparam->initncell = 0;
  histparam->min = 0;
  histparam->initmin = 0;
  histparam->max = 0;
  histparam->initmax = 0;
  histparam->precision = 0;
  histparam->initprecision = 0;
}


int  MT_InitDescriptor ( vt_histDescriptor *desc, double l, double u, int nelt,
    vt_histparam param)
{
  char *proc = "MT_InitDescriptor";

  double precision;
  double upper, lower;
  int ncell;

  desc->total=0;

  if(param.initncell != 0 && param.ncell <= 0) {
    fprintf(stderr,"%s: unable to initialise descriptor, wrong ncell initialisation\n",proc);
    return(-1);
  }
  if(param.initprecision != 0 && param.initprecision <= 0 ) {
    fprintf(stderr,"%s: unable to initialise descriptor, wrong precision initialisation\n",proc);
    return(-1);
  }
  if(param.initmin != 0 && param.initmax <= 0 && param.min>=param.max) {
    fprintf(stderr,"%s: unable to initialise descriptor, initialisation with min>=max\n",proc);
    return(-1);
  }

  if (param.initmin == 0) {
    if(param.initmax == 0 ) {
      if (param.initprecision == 0) {
        if( param.initncell == 0) {
          upper=u;
          lower=l;
          precision=(upper-lower)/(double)(nelt-1);
          ncell=(int)((upper-lower)/precision)+1;
        }
        else{ // initcell != 0
          ncell=param.ncell;
          upper=u;
          lower=l;
          precision=(upper-lower)/(double)(ncell);
        }
      }
      else // initprecision != 0
      {
        if(param.initncell == 0) {
          upper=u;
          lower=l;
          precision=param.precision;
          ncell=(int)((upper-lower)/precision)+1;
        }
        else { // initncell != 0
    		  fprintf(stderr,"%s: initmin=%f, initmax=%f, initncell=%d, initprecision=%f\n", proc, param.min, param.max, param.initncell, param.precision);
          fprintf(stderr,"%s: unable to initialise descriptor, bad parameters initialisation A\n",proc);
          return(-1);
        }
      }
    }
    else { // initmax != 0
      upper = param.max;
      if (param.initprecision == 0) {
        if( param.initncell == 0) {
          lower=l;
          precision=(upper-lower)/(double)(nelt-1);
          ncell=(int)((upper-lower)/precision)+1;
        }
        else{ // initcell != 0
          lower=l;
          ncell=param.ncell;
          precision=(upper-lower)/(double)(ncell);
        }
      }
      else // initprecision != 0
      {
        if(param.initncell == 0) {
          lower=l;
          precision=param.precision;
          ncell=(int)((upper-lower)/precision)+1;
        }
        else { // initncell != 0
          precision=param.precision;
          ncell=param.ncell;
          lower=upper-(double)ncell*precision;
        }
      }
    }
  }
  else { // initmin != 0
    lower = param.min;
    if(param.initmax == 0 ) {
      if (param.initprecision == 0) {
        if( param.initncell == 0) {
          upper=u;
          precision=(upper-lower)/(double)(nelt-1);
          ncell=(int)((upper-lower)/precision)+1;
        }
        else{ // initcell != 0
          upper=u;
          ncell=param.ncell;
          precision=(upper-lower)/(double)(ncell);
        }
      }
      else // initprecision != 0
      {
        if(param.initncell == 0) {
          upper=u;
          precision=param.precision;
          ncell=(int)((upper-lower)/precision)+1;
        }
        else { // initncell != 0
          precision=param.precision;
          ncell=param.ncell;
          upper=lower+(double)ncell*precision;
        }
      }
    }
    else { // initmax != 0
      upper=param.max;
      if (param.initprecision == 0) {
        if( param.initncell == 0) {
          precision=(upper-lower)/(double)(nelt-1);
          ncell=(int)((upper-lower)/precision)+1;
        }
        else{ // initcell != 0
          ncell=param.ncell;
          precision=(upper-lower)/(double)(ncell);
        }
      }
      else // initprecision != 0
      {
        if(param.initncell == 0) {
          precision=param.precision;
          ncell=(int)((upper-lower)/precision)+1;
        }
        else { // initncell != 0
      			fprintf(stderr,"%s: initmin=%f, initmax=%f, initncell=%d, initprecision=%f\n", proc, param.min, param.max, param.initncell, param.precision);
            fprintf(stderr,"%s: unable to initialise descriptor, bad parameters initialisation B\n",proc);
            return(-1);
        }

      }
    }

  }

  desc->upper=upper;
  desc->lower=lower;
  desc->precision= precision;
  desc->ncell=ncell;
  desc->total=0;

  return(1);
}

int  MT_InitDescriptorWithNCell ( vt_histDescriptor *desc, double l, double u, int ncell)
{
  desc->upper=u;
  desc->lower=l;
  desc->ncell=ncell;
  desc->precision=(u-l)/(double)(ncell);
  desc->total=0;
  return(1);
}




int MT_WriteHisto3D(vt_anisotropicHisto histo, int cptLevenberg, vt_LM_type *LM_types, int nbLeven, double *parLevenOutX, double *parLevenOutY, double *parLevenOutZ, double **thresholds,
                    double *ratiosIn, int nratio, double *thX, double *thY, double *thZ, double *MSF, int nMSF, double *MT, int nMT, vt_binmode *binarise, int nbin, double **ratiosOut, char *filename, char *prefix)
{
  int i,j,n,axe;
  double ratio;
  double *parLevenOut;
  char word[STRINGLENGTH];
  /* Generation du fichier */

  FILE *fichier;
  if (filename != NULL && filename[0] == '>')
    fichier = stdout;
  else
    fichier = fopen (filename, "w" );

  if (fichier == NULL)
  {
    perror (filename);
  }
  else
  {
    fprintf(fichier, "%-10s\t\t%s%s", "Histogram", "File ", prefix);

    /* histX */
    fprintf(fichier, "\n\nX projection\n\n");
    /*    -- description */
    fprintf(fichier, "\t%-14s\n", "X description :");
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "lower", histo.descriptorX.lower);
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "upper", histo.descriptorX.upper);
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "delta", histo.descriptorX.precision);
    fprintf(fichier, "\t\t%-7s : %d\n", "ncell", histo.descriptorX.ncell);
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "total", histo.descriptorX.total);
    /*    -- histogramme */
    fprintf(fichier, "\t%-14s\n\t", "X histogram : ");
    for (i=0;i<histo.descriptorX.ncell;i++)
      fprintf(fichier,"\t%10f", histo.histX[i]);

    /* histY */
    fprintf(fichier, "\n\nY projection\n\n");
    /*    -- description */
    fprintf(fichier, "\t%-14s\n", "Y description :");
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "lower", histo.descriptorY.lower);
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "upper", histo.descriptorY.upper);
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "delta", histo.descriptorY.precision);
    fprintf(fichier, "\t\t%-7s : %d\n", "ncell", histo.descriptorY.ncell);
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "total", histo.descriptorY.total);
    /*    -- histogramme */
    fprintf(fichier, "\t%-14s\n\t", "Y histogram : ");
    for (i=0;i<histo.descriptorY.ncell;i++)
      fprintf(fichier,"\t%10f", histo.histY[i]);

    /* histZ */
    fprintf(fichier, "\n\nZ projection\n\n");
    /*    -- description */
    fprintf(fichier, "\t%-14s\n", "Z description :");
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "lower", histo.descriptorZ.lower);
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "upper", histo.descriptorZ.upper);
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "delta", histo.descriptorZ.precision);
    fprintf(fichier, "\t\t%-7s : %d\n", "ncell", histo.descriptorZ.ncell);
    fprintf(fichier, "\t\t%-7s : %7.2f\n", "total", histo.descriptorZ.total);
    /*    -- histogramme */
    fprintf(fichier, "\t%-14s\n\t", "Z histogram : ");
    for (i=0;i<histo.descriptorZ.ncell;i++)
      fprintf(fichier,"\t%10f", histo.histZ[i]);
    fprintf(fichier, "\n");

    /* Fitting levenberg */
    if(cptLevenberg == 1)
    {
      for (axe = 0 ; axe < 3 ; axe ++) {
        j=0;
        switch (axe) {
        case 0 :
          parLevenOut = parLevenOutX;
          fprintf(fichier, "\nAjustement Axe X :\n");
          break;
        case 1 :
          parLevenOut = parLevenOutY;
          fprintf(fichier, "\nAjustement Axe Y :\n");
          break;
        case 2 :
          parLevenOut = parLevenOutZ;
          fprintf(fichier, "\nAjustement Axe Z :\n");
          break;
        default :
          fprintf(stderr, "Valeur d'axe inattendue (axe=%d)\n",axe);
          return(-1);
        }
        for ( i=0 ; i<nbLeven ; i++)
        {
          switch ( LM_types[i] ) {
          case _GAUSS_ :
            fprintf( fichier, "Gauss            : amplitude = %f moyenne = %f sigma = %f\n",
                parLevenOut[j+0], parLevenOut[j+1], parLevenOut[j+2] );
            j+=3;
            break;
          case _RAYLEIGH_ :
            fprintf( fichier, "Rayleigh         : amplitude = %f moyenne = %f sigma = %f\n",
                parLevenOut[j+0], parLevenOut[j+1], parLevenOut[j+2]);
            j+=3;
            break;
          case _RAYLEIGH_POS_ :
            fprintf( fichier, "RayleighPos      : amplitude = %f moyenne = %f sigma = %f\n",
                parLevenOut[j+0], parLevenOut[j+1], parLevenOut[j+2]);
            j+=3;
            break;
          case _RAYLEIGH_CENTERED_ :
            fprintf( fichier, "RayleighCentered : amplitude = %f sigma = %f\n",
                parLevenOut[j+0], parLevenOut[j+1]);
            j+=2;
            break;
          case _Xn_ :
            fprintf( fichier, "Xn          : coeff = %f exposant = %f\n",
                parLevenOut[j+0], parLevenOut[j+1]);
            j+=2;
            break;
          default :
            fprintf(stderr,"Erreur, valeur inattendue dans LM_types[i]\n");
            return(-1);
          }
        }
      }
      fprintf(fichier, "\n");
    }

    /* anisotropic thresholds */
    if (nratio>0 && (nMSF==0 && nMT==0))
    {
      fprintf(fichier, "\n\nAnisotropic thresholds\n");
      for (i=0;i<nratio;i++)
      {
        ratio=ratiosIn[i];
        fprintf(fichier, "\n\t%-14s%f\n", "Global rate : ", ratio);
        fprintf(fichier, "\t%-14s%f\n", "X threshold : ", thresholds[0][i]);
        fprintf(fichier, "\t%-14s%f\n", "Y threshold : ", thresholds[1][i]);
        fprintf(fichier, "\t%-14s%f\n", "Z threshold : ", thresholds[2][i]);
        for(n=0;n<nbin;n++)
        {
          switch (binarise[n]) {
          default :
          case _NO_BIN_ :
            sprintf(word,"%s rate : ","NoBin ");
            break;
          case _BIN_PROJ_ :
            sprintf(word,"%s rate : ","Proj  ");
            break;
          case _BIN_DOT_ :
            sprintf(word,"%s rate : ","Dot   ");
            break;
          case _BIN_AND_ :
            sprintf(word,"%s rate : ","And   ");
            break;
          case _BIN_OR_ :
            sprintf(word,"%s rate : ","Or    ");
            break;
          }
          fprintf(fichier, "\t%-14s%f\n", word, ratiosOut[i][n]);
        }
      }
    }

    if (nratio==0 && (nMSF>0 || nMT>0))
    {
      if (nMSF>0 && nMT>0)
      {
        fprintf(stderr, "nMSF > 0 and nMT > 0 at the same time...\n");
        return(-1);
      }
      if (nMSF>0)
        fprintf(fichier, "\n\nMSF Anisotropic thresholds\n");
      if (nMT>0)
        fprintf(fichier, "\n\nMT Anisotropic thresholds\n");
      for (i=0;i<nMSF+nMT;i++)
      {
        if (nMSF>0)
          fprintf(fichier, "\n\t%-14s%f\n", "Global rate : ", MSF[i]);
        if (nMT>0)
          fprintf(fichier, "\n\t%-14s%f\n", "Global rate : ", MT[i]);
        fprintf(fichier, "\t%-14s%f\n", "X threshold : ", thX[i]);
        fprintf(fichier, "\t%-14s%f\n", "Y threshold : ", thY[i]);
        fprintf(fichier, "\t%-14s%f\n", "Z threshold : ", thZ[i]);
      }
    }

  }
  return(1);
}


int MT_ThresholdsFromHisto( double **thres,
                             vt_anisotropicHisto histo,
                             double *rates,
                             double nrates)
{
  int i,n;
  int sum;
  double rate;
  double low, delta;
  int total;

  for(i=0;i<3;i++)
    thres[i]=malloc(nrates*sizeof(double));

  for (n=0;n<nrates;n++)
  {
    rate=rates[n];
    if (rate<0 || rate >1)
    {
      for (i=0;i<3;i++)
      {
        free(thres[i]);
        thres[i]=NULL;
      }
      return(-1);
    }
    // histX
    low=histo.descriptorX.lower;
    delta=histo.descriptorX.precision;
    total=histo.descriptorX.total;
    i=0;
    sum=histo.histX[0];
    while(sum<total*(1-rate))
    {
      sum+=histo.histX[++i];
    }
    thres[0][n]=low+delta*i;

    // histY
    low=histo.descriptorY.lower;
    delta=histo.descriptorY.precision;
    total=(double)histo.descriptorY.total;
    i=0;
    sum=histo.histY[0];
    while(sum<total*(1-rate))
    {
      sum+=histo.histY[++i];
    }
    thres[1][n]=low+delta*i;

    // histZ
    low=histo.descriptorZ.lower;
    delta=histo.descriptorZ.precision;
    total=(double)histo.descriptorZ.total;
    i=0;
    sum=histo.histZ[0];
    while(sum<total*(1-rate))
    {
      sum+=histo.histZ[++i];
    }
    thres[2][n]=low+delta*i;
  }
  return(1);
}

int MT_Binarise(vt_Rep_Angles imIn, vt_image *imResAnisotropes,
                 vt_image *imBin, char *name,
                 double thX, double thY, double thZ, vt_binmode bin,
                 double *r)
{
  char *proc = "MT_Binarise";

  int i;
  double thres;
  double n[3];
  unsigned char binX,binY,binZ;

  float *theTht, *thePhi;
  float *theResX, *theResY, *theResZ, *theRes;
  unsigned char *theBin;
  int total=0;

  VT_InitFromImage( imBin, imResAnisotropes, name,  (int)UCHAR );
  if ( VT_AllocImage( imBin ) != 1 ) {
    fprintf(stderr, "%s: unable to allocate binary image\n", proc);
    return( -1 );
  }

  theResX = (float *)imResAnisotropes[0].buf;
  theResY = (float *)imResAnisotropes[1].buf;
  theResZ = (float *)imResAnisotropes[2].buf;
  theRes  = (float *)imIn.imrep->buf;
  theTht  = (float *)imIn.imtheta->buf;
  thePhi  = (float *)imIn.imphi->buf;
  theBin  = (unsigned char *) imBin->buf;

  for (i=0;i<imIn.imrep->dim.x*imIn.imrep->dim.y*imIn.imrep->dim.z; i++)
  {
    SphericalAnglesToUnitVector( theTht[i], thePhi[i], n);
    switch (bin) {
    case _BIN_PROJ_:
      binX=binY=binZ=(unsigned char) 1;
      if((fabs(n[0])>=fabs(n[1])) && (fabs(n[0])>=fabs(n[2])))
        binX = (theResX[i] > thX ? (unsigned char) 1 : (unsigned char) 0);
      if((fabs(n[1])>=fabs(n[0])) && (fabs(n[1])>=fabs(n[2])))
        binY = (theResY[i] > thY ? (unsigned char) 1 : (unsigned char) 0);
      if((fabs(n[2])>=fabs(n[0])) && (fabs(n[2])>=fabs(n[1])))
        binZ = (theResZ[i] > thZ ? (unsigned char) 1 : (unsigned char) 0);
      theBin[i] = ((unsigned char)255)*binX*binY*binZ;
      break;
    case _BIN_DOT2_:
      thres = n[0]*n[0]*thX+n[1]*n[1]*thY+n[2]*n[2]*thZ;
      theBin[i]=(theRes[i] >thres ? (unsigned char) 255 : (unsigned char) 0);
      break;
    case _BIN_DOT_ :
      thres = fabs(n[0])*thX+fabs(n[1])*thY+fabs(n[2])*thZ;
      theBin[i]=(theRes[i] >thres ? (unsigned char) 255 : (unsigned char) 0);
      break;
    case _BIN_OR_ :
      theBin[i]=(theResZ[i] > thZ ? (unsigned char) 255 :
      (theResY[i] > thY ? (unsigned char) 255 :
      (theResX[i] > thX ? (unsigned char) 255 : (unsigned char) 0)));
      break;
    case _BIN_AND_ :
      theBin[i]=(theResZ[i] > thZ ? (theResY[i] > thY ?
          (theResX[i] > thX ? (unsigned char) 255 : (unsigned char) 0) :
          (unsigned char) 0) : (unsigned char) 0);
      break;
    case _NO_BIN_ :
    default:
      VT_FreeImage(imBin);
      fprintf(stderr, "%s: bad binmode argument\n", proc);
      return(-1);
    }
    if (theBin[i] != '\0')
      total++;
  }
  *r = ((double)total)/((double)i);
  return(1);
}



int  VT_AllocRep_AnglesWithImageRep( vt_Rep_Angles *par)
{
  char *name;
  char prefix[256];
  int i=0;
  sprintf(prefix, "%s", par->imrep->name);

  vt_image *imrep = par->imrep;
  vt_image *imtheta=par->imtheta;

  name = imtheta->name;
  while(prefix[i] != '\0')
    i++;
  if(i>4)
    i -=4 ;
  prefix[i]='\0';
  fprintf(stdout, "STRINGLENGTH : %d\n", (int)STRINGLENGTH);
  sprintf( name, "imtheta");

  fprintf(stdout, "Image d'entree : %s\n", imrep->name);
  fprintf(stdout, "dimx : %lu\tdimy :%lu\tdimz : %lu\n", imrep->dim.x, imrep->dim.y, imrep->dim.z);

  fprintf(stdout, "STRINGLENGTH : %d\n", (int)STRINGLENGTH);
  sprintf( imtheta->name, "imtheta");

  fprintf(stdout, "STRINGLENGTH : %d\n", (int)STRINGLENGTH);
  imtheta->type = (int) DOUBLE;
  imtheta->cpu = MY_CPU;
  imtheta->buf = (void*)NULL;
  imtheta->array = (void***)NULL;
  imtheta->nuser = 0;
  imtheta->user = NULL;




  if ( VT_AllocImage( (par->imtheta) ) != 1 ) {
    VT_FreeImage( (par->imrep) );
    return( -1 );
  }
  sprintf( name, "imphi");

  fprintf(stdout, "Image d'entree : %s\n", name);

  VT_InitFromImage( (par->imphi), (par->imrep), name,  (int)DOUBLE );
  if ( VT_AllocImage( (par->imphi) ) != 1 ) {
    VT_FreeImage( (par->imrep) );
    VT_FreeImage( (par->imtheta) );
    return( -1 );
  }
  return(1);
}


void VT_FreeRep_Angles ( vt_Rep_Angles *par )
{
  VT_FreeImage( (par->imrep) );
  VT_FreeImage( (par->imtheta) );
  VT_FreeImage( (par->imphi) );
}

int  VT_AllocAnisotropicHisto ( vt_anisotropicHisto *par)
{
  int dim,i;
  double *h;
  dim=par->descriptorX.ncell;
  par->histX=malloc(dim*sizeof(double));
  h=par->histX;
  for (i=0;i<dim;i++) h[i]=0;
  dim=par->descriptorY.ncell;
  par->histY=malloc(dim*sizeof(double));
  h=par->histY;
  for (i=0;i<dim;i++) h[i]=0;
  dim=par->descriptorZ.ncell;
  par->histZ=malloc(dim*sizeof(double));
  h=par->histZ;
  for (i=0;i<dim;i++) h[i]=0;
  return(1);
}


void VT_FreeAnisotropicHisto ( vt_anisotropicHisto *par )
{
  free(par->histX); par->histX = NULL;
  free(par->histY); par->histY = NULL;
  free(par->histZ); par->histZ = NULL;
}

