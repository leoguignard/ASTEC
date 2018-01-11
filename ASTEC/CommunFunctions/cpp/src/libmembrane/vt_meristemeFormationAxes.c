/*************************************************************************
 * vt_meristemeFormationAxes.c -
 *
 * $Id: vt_meristemeFormationAxes.c,v 1.0 2014/06/16 16:07:34 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 * 
 * CREATION DATE: 
 * 2014/06/16 
 *
 * ADDITIONS, CHANGES
 *
 *
 */
 

#include <vt_common.h>
#include <transfo.h>
#include <mt_membrane3D.h>

#include <vt_meristemeFormationAxes.h>

#define FLTZERO 1e-8

static int _verbose_ = 1;

int VT_conicVote(vt_image **imagesIn, vt_image *imres, double r, double alpha, int NangleIter )
{
  char *proc = "VT_conicVote";
  double **angles;
  int Nangles;
  
  int dim[3];

  int dimFields[3];
  unsigned char ****fields;

  int i, j;
  int x,y,z;
  int yy,zz;
  int n;
  
  vt_image imbin=*imagesIn[0];
  vt_image imtht=*imagesIn[1];;
  vt_image imphi=*imagesIn[2];;
  unsigned char ***bin;
  float ***tht, ***phi, ***array;

  dim[0]=imbin.dim.x;
  dim[1]=imbin.dim.y;
  dim[2]=imbin.dim.z;

  dimFields[0]=2*(int)r+1;
  dimFields[1]=2*(int)r+1;
  dimFields[2]=2*(int)r+1;

  switch (imbin.type) {
  case UCHAR:
    bin=(unsigned char ***)imbin.array;
    break;
  default:
    fprintf(stderr, "%s: input image type not handled yet\n", proc);
    return( 0 );
  }

  switch (imtht.type) {
  case FLOAT:
    tht=(float ***)imtht.array;
    break;
  default:
    fprintf(stderr, "%s: input image type not handled yet\n", proc);
    return( 0 );
  }

  switch (imphi.type) {
  case FLOAT:
    phi=(float ***)imphi.array;
    break;
  default:
    fprintf(stderr, "%s: input image type not handled yet\n", proc);
    return( 0 );
  }

  switch (imres->type) {
  case FLOAT:
    array=(float ***)imres->array;
    break;
  default:
    fprintf(stderr, "%s: output image type not handled yet\n", proc);
    return( 0 );
  }


  /* Calcul des angles répartis de façon homogène sur la boule unite */
  if (_verbose_)
    fprintf(stdout, "%s: computating angles...\n", proc);

  Nangles = MT_Compute3DAngles(&angles, NangleIter);
  if (Nangles <= 0) {        
    if (_verbose_)
      fprintf(stderr, "%s: error in computating angles\n", proc);
    return( 0 );
  }
  
  /* Allocation des champs de vote coniques */
  if(_verbose_) fprintf(stdout, "%s: allocating conic fields...\n", proc);

  fields = malloc(Nangles*sizeof(unsigned char***));
  if (fields==NULL)
  {
	for (n=0;n<Nangles;n++){
      free(angles[n]);
      angles[n]=NULL;
    }
    fprintf(stderr, "%s: error in allocating fields\n", proc);
    return(0);
  }
  for(i=0; i<Nangles; i++)
  {
    fields[i] = malloc(dimFields[2]*sizeof(unsigned char**));
	if (fields[i]==NULL)
	{
	  for (n=0;n<Nangles;n++){
        free(angles[n]);
        angles[n]=NULL;
      }
      for (j=0; j<i; j++){
		free(fields[j]);
		fields[j]=NULL;
	  }
	  free(fields);
	  fields=NULL;
      fprintf(stderr, "%s: error in allocating fields\n", proc);
      return(0);
    }
    for(z=0; z<dimFields[2]; z++)
    {
      fields[i][z] = malloc(dimFields[1]*sizeof(unsigned char*));
	  if (fields[i][z]==NULL)
	  {
		for (n=0;n<Nangles;n++){
          free(angles[n]);
          angles[n]=NULL;
		}
		for (j=0; j<=i; j++){
		  for (zz=0; zz<(j==i) ? z : dimFields[2] ; zz++){
		    free(fields[j][zz]); fields[j][zz]=NULL;
		  }
		  free(fields[j]); fields[j]=NULL;
	    }
		free(fields);
		fields=NULL;
		fprintf(stderr, "%s: error in allocating fields\n", proc);
		return(0);
	  }
      for(y=0; y<dimFields[1]; y++)
      {
        fields[i][z][y] = malloc(dimFields[0]*sizeof(unsigned char));
		if(fields[i][z][y]==NULL)
		{
		  for (n=0;n<Nangles;n++){
			free(angles[n]);
			angles[n]=NULL;
		  }
		  for (j=0; j<=i; j++){
			for (zz=0; zz<(j==i) ? z : dimFields[2] ; zz++){
			  for (yy=0; yy<(j==i && zz==z) ? y : dimFields[1] ; yy++) {
				free(fields[j][zz][yy]); fields[j][zz][yy]=NULL;
			  }
			  free(fields[j][zz]); fields[j][zz]=NULL;
			}
			free(fields[j]); fields[j]=NULL;
	      }
		  free(fields);
		  fields=NULL;
		  fprintf(stderr, "%s: error in allocating fields\n", proc);
		  return(0);
		}
      }
    }
  }
  
  /* Calcul des champs */
  if(_verbose_) fprintf(stdout, "%s: computing conic voting fields...\n", proc);
  for (i=0; i<Nangles; i++)
	VT_ComputeConicField3D(fields[i], angles[i], dimFields, r, alpha);
  
  /* Accumulation des votes */
  if(_verbose_) fprintf(stdout, "%s: computing votes...\n", proc);
  for (z=0; z<dim[2]; z++) 
  { 
	//fprintf(stdout, "z=%d on %d\n", z, dim[2]); 
	if ( _verbose_ && (z*100/dim[2] % 10 == 0) && ((z-1)*100/dim[2] % 10 != 0))
	  fprintf(stdout, "Slice %d\ton %d\t(%d %%)\n", z+1, dim[2], (z*100/dim[2]));

	for (y=0; y<dim[1]; y++)
	for (x=0; x<dim[0]; x++)
	{
	  if(bin[z][y][x]==(char)0) continue;
	  n = VT_nearestAngle(tht[z][y][x], phi[z][y][x], angles, Nangles);
	  VT_AddConicField(array, dim, fields[n], dimFields, (float)bin[z][y][x], x,y,z );

    }
  
  }
  

  /* Libération mémoire */
  if(_verbose_) fprintf(stdout, "%s: free memory...\n", proc);
  for (n=0;n<Nangles;n++){
	free(angles[n]);
	angles[n]=NULL;
  }
  for (j=0; j<Nangles; j++){
	for (zz=0; zz<dimFields[2] ; zz++){
	  for (yy=0; yy<dimFields[1] ; yy++) {
		free(fields[j][zz][yy]); //fields[j][zz][yy]=NULL;
	  }
	  free(fields[j][zz]); //fields[j][zz]=NULL;
	}
	free(fields[j]); //fields[j]=NULL;
  }
  free(fields);
  fields=NULL;

  
  return(1);
}

void VT_AddConicField(float ***array, int *dim, unsigned char ***field, int *dimField, float coef, int x, int y, int z )
{
  
  float val;

  int xf0, yf0, zf0;
  //int xf1, yf1, zf1;
  int x0, y0, z0, x1, y1, z1;
  int i, j, k;
  int ifield, jfield, kfield;


  int hsf[3];

  int xmax = dim[0];
  int ymax = dim[1];
  int zmax = dim[2];

  hsf[0]=dimField[0]/2;
  hsf[1]=dimField[1]/2;
  hsf[2]=dimField[2]/2;


  if (hsf[0]<=x)
  {
    xf0 = 0;
    x0 = x-hsf[0];
  }
  else
  {
    x0 = 0;
    xf0 = hsf[0]-x;
  }
  if (hsf[1]<=y)
  {
    yf0 = 0;
    y0 = y-hsf[1];
  }
  else
  {
    y0 = 0;
    yf0 = hsf[1]-y;
  }
  if (hsf[2]<=z)
  {
    zf0 = 0;
    z0 = z-hsf[2];
  }
  else
  {
    z0 = 0;
    zf0 = hsf[2]-z;
  }

  //xf0 = (hsf[0]<=x) ? 0 : hsf[0]-x;
  //yf0 = (hsf[1]<=y) ? 0 : hsf[1]-y;
  //zf0 = (hsf[2]<=z) ? 0 : hsf[2]-z;

  //xf1 = (x+hsf[0]<xmax) ? dimField[0] : (xmax-1-x)+hsf[0]+1;
  //yf1 = (y+hsf[1]<ymax) ? dimField[1] : (ymax-1-y)+hsf[1]+1;
  //zf1 = (z+hsf[2]<zmax) ? dimField[2] : (zmax-1-z)+hsf[2]+1;

  x1 = (x+hsf[0]<xmax) ? x+hsf[0]+1 : xmax;
  y1 = (y+hsf[1]<ymax) ? y+hsf[1]+1 : ymax;
  z1 = (z+hsf[2]<zmax) ? z+hsf[2]+1 : zmax;

  kfield = zf0;
  for (k=z0; k<z1;k++)
  {
    jfield = yf0;
    for (j=y0; j<y1; j++)
    {
      ifield = xf0;
      for (i=x0; i<x1;i++)
      {
		val=(float)field[kfield][jfield][ifield];
		if (val>0)
		{
		  val*=coef;
		  array[k][j][i]+=val;
		}
        ifield++;
      }
      jfield++;
    }
    kfield++;
  }
  return;
}

int VT_nearestAngle(double theta, double phi, double **angles, int Nangles)
{
  int n, N=0;
  double dot, DOT=0;

  double v[3], v0[3];
  SphericalAnglesToUnitVector( theta, phi, v0);

  for (n=0;n<Nangles;n++)
  {
    SphericalAnglesToUnitVector( angles[n][0], angles[n][1], v);
    dot = fabs(v[0]*v0[0]+v[1]*v0[1]+v[2]*v0[2]);
    if (dot<DOT)
      continue;
    DOT = dot;
    N = n;
    if (1-DOT < FLTZERO)
      break;
  }

  return (N);
}

void VT_ComputeConicField3D(unsigned char ***field, double *angles, int *dimField, double r, double alpha)
{
  int O[3];
  int x,y,z;
  double norme, norme2;
  double dot;
  double cosalpha;
  double n[3];
  double r2=r*r;
    
  O[0]=(int)(0.5*(dimField[0]-1));
  O[1]=(int)(0.5*(dimField[1]-1));
  O[2]=(int)(0.5*(dimField[2]-1));
  
  cosalpha=cos(alpha); 
  if (cosalpha<0) cosalpha*=-1;
  SphericalAnglesToUnitVector( angles[0], angles[1], n);


  for (z=0 ; z<dimField[2] ; z++)
  for (y=0 ; y<dimField[1] ; y++)
  for (x=0 ; x<dimField[0] ; x++)
  {
	norme2=(x-O[0])*(x-O[0])+(y-O[1])*(y-O[1])+(z-O[2])*(z-O[2]);
	if (norme2>r2) continue;
	norme=sqrt(norme2);
	dot=(x-O[0])*n[0]+(y-O[1])*n[1]+(z-O[2])*n[2];
	if (dot<0) dot=-dot;
	if (dot<cosalpha*norme) continue;
	field[z][y][x]=(char)1;
  }

  return;
}


int addConeToBuf(void ***array, int *dim, bufferType t, int *pt, double tht, double phi, double r, double alpha)
{
  char *proc = "addConeToBuf";
  

  //int dimFields[3];
  

  //double **angles;
  //int Nangles;

  //double theCoeff[3];
  //theCoeff[0] = theCoeff[1] = scale;
  //theCoeff[2] = zfact*scale;
  //dimFields[0] = dimFields[1] = dimFields[2] = (int)(2*r+1);

  int x,y,z;
  double norme, norme2, dot;
  double cosalpha;
  double n[3];
  double r2=r*r;
  
  float*** res32;
  //unsigned char ***res8;
  

  switch (t)
  {
  //case UCHAR :
  //  res8 = (unsigned char***)array;
  //  break;
  case FLOAT :
    res32 = (float***)array;
    break;
  default:
    fprintf(stderr, "%s: array type not handled yet...\n", proc);
    return(0);
  }
  
  cosalpha=cos(alpha); 
  if (cosalpha<0) cosalpha*=-1;
  SphericalAnglesToUnitVector( tht, phi, n);
  
  int xmin, xmax, ymin, ymax, zmin, zmax;
  xmin=(pt[0]<r)? 0: pt[0]-r;
  ymin=(pt[1]<r)? 0: pt[1]-r;
  zmin=(pt[2]<r)? 0: pt[2]-r;
  xmax=(pt[0]>=dim[0]-r)? dim[0]-1: pt[0]+r;
  ymax=(pt[1]>=dim[1]-r)? dim[1]-1: pt[1]+r;
  zmax=(pt[2]>=dim[2]-r)? dim[2]-1: pt[2]+r;
  
  //fprintf(stdout, "pt     = \t[%d %d %d]\t", pt[0], pt[1], pt[2]);
  //fprintf(stdout, "dim    = \t[%d %d %d]\n", dim[0], dim[1], dim[2]);
  //fprintf(stdout, "bornes = \t{%d %d}\t{%d %d}\t{%d %d}\n", xmin, xmax, ymin, ymax, zmin, zmax);
  
  for (z=zmin ; z<zmax ; z++)
  for (y=ymin ; y<ymax ; y++)
  for (x=xmin ; x<xmax ; x++)
  {
	//fprintf(stdout, "{x y z} = {%d %d %d}\t", x, y, z);
	norme2=(x-pt[0])*(x-pt[0])+(y-pt[1])*(y-pt[1])+(z-pt[2])*(z-pt[2]);
	if (norme2>r2) continue;
	norme=sqrt(norme2);
	dot=(x-pt[0])*n[0]+(y-pt[1])*n[1]+(z-pt[2])*n[2];
	if (dot<0) dot=-dot;
	if (dot<cosalpha*norme) continue;
	res32[z][y][x]+=1;
  }
  
  return ( 1 );
}

int addLineToBuf(void *buf, 	// buffer array
				 int *dim, 		// buffer dimensions
				 bufferType t, 	// buffer type
				 int *pt, 		// origin point
				 double tht, 	// angle 1
				 double phi, 	// angle 2
				 int r			// rayon de "vote"
				)
{
  char *proc = "addLineToBuf";

  int Pt1[3],Pt2[3];
  int i, j;
  double n[3];
  double lambdas[6];
  double tmp;
  double maxdim;
  float *theBuf;
  
  //size_t size;
  int ix, ixy;
  int dx, dy, dz;  
  int ax, ay, az;
  int sx, sy, sz;
  int x, y, z;
  int xd, yd, zd;

  int r2=r*r;
  double d2;

  
  // n : director vector of the line
  SphericalAnglesToUnitVector(tht, phi, n);
  
  // Computing the limit points Pt1, Pt2 of the image
  maxdim=(double)(dim[0]<dim[1]) ? ((dim[1]<dim[2]) ? dim[2] : dim[1]) : ((dim[0]<dim[2]) ? dim[2] : dim[0]);
  lambdas[0]=(n[0]==0) ? 2*maxdim : (0-(double)pt[0])/n[0];
  lambdas[1]=(n[1]==0) ? 2*maxdim : (0-(double)pt[1])/n[1];
  lambdas[2]=(n[2]==0) ? 2*maxdim : (0-(double)pt[2])/n[2];
  lambdas[3]=(n[0]==0) ? 2*maxdim : ((double)dim[0]-(double)pt[0]-1)/n[0];
  lambdas[4]=(n[1]==0) ? 2*maxdim : ((double)dim[1]-(double)pt[1]-1)/n[1];
  lambdas[5]=(n[2]==0) ? 2*maxdim : ((double)dim[2]-(double)pt[2]-1)/n[2];
  
    
  /* Sort the lambdas
   * 
  procédure tri_insertion
   * 
   */
  for (i=1; i<6; i++) {
	tmp=lambdas[i];
	j=i;
	while (j>0 && lambdas[j-1]>tmp) {
	  lambdas[j]=lambdas[j-1];
	  j--;
	}
	lambdas[j]=tmp;
  }
  
  // Pt1 
  i=0;
  while((pt[0]+lambdas[i]*n[0]+0.5 < 0 || pt[0]+lambdas[i]*n[0]+0.5 >= dim[0] ||
		 pt[1]+lambdas[i]*n[1]+0.5 < 0 || pt[1]+lambdas[i]*n[1]+0.5 >= dim[1] ||
		 pt[2]+lambdas[i]*n[2]+0.5 < 0 || pt[2]+lambdas[i]*n[2]+0.5 >= dim[2] ) && i<6 )
	i++;
  if (i>=6)
  {
	fprintf(stderr,"%s: error while computing the limit points: 1st point not found\n", proc);
	return(0);
  }
  //i1=i;
  Pt1[0]= pt[0] + (int) (lambdas[i]*n[0]+0.5); 
  Pt1[1]= pt[1] + (int) (lambdas[i]*n[1]+0.5); 
  Pt1[2]= pt[2] + (int) (lambdas[i]*n[2]+0.5); 
  
  // Pt2
  i=5;
  while((pt[0]+lambdas[i]*n[0]+0.5 < 0 || pt[0]+lambdas[i]*n[0]+0.5 >= dim[0] ||
		 pt[1]+lambdas[i]*n[1]+0.5 < 0 || pt[1]+lambdas[i]*n[1]+0.5 >= dim[1] ||
		 pt[2]+lambdas[i]*n[2]+0.5 < 0 || pt[2]+lambdas[i]*n[2]+0.5 >= dim[2] ) && i>=0 )
	i--;
  if (i<0)
  {
	fprintf(stderr,"%s: error while computing the limit points: 2nd point not found\n", proc);
	return(0);
  }
  //i2=i;
  Pt2[0]= pt[0] + (int) (lambdas[i]*n[0]+0.5); 
  Pt2[1]= pt[1] + (int) (lambdas[i]*n[1]+0.5); 
  Pt2[2]= pt[2] + (int) (lambdas[i]*n[2]+0.5); 

  
  
  if(r>0)				// Computing Pt1, Pt2 so that voting distance < r  
  {
	// Pt1 :
	
	d2=(Pt1[0]-pt[0])*(Pt1[0]-pt[0])+(Pt1[1]-pt[1])*(Pt1[1]-pt[1])+(Pt1[2]-pt[2])*(Pt1[2]-pt[2]);
	
	if (r2<d2)
	{
	  if(n[0]*(Pt1[0]-pt[0])+n[1]*(Pt1[1]-pt[1])+n[2]*(Pt1[2]-pt[2])>0)
	  {
		Pt1[0]=(int) (pt[0] + r*n[0] + 0.5);
		Pt1[1]=(int) (pt[1] + r*n[1] + 0.5);
		Pt1[2]=(int) (pt[2] + r*n[2] + 0.5);
	  }
	  else
	  {
		Pt1[0]=(int) (pt[0] - r*n[0] + 0.5);
		Pt1[1]=(int) (pt[1] - r*n[1] + 0.5);
		Pt1[2]=(int) (pt[2] - r*n[2] + 0.5);
	  }
	}
	
	// Pt2 :

	d2=(Pt2[0]-pt[0])*(Pt2[0]-pt[0])+(Pt2[1]-pt[1])*(Pt2[1]-pt[1])+(Pt2[2]-pt[2])*(Pt2[2]-pt[2]);
	
	if (r2<d2)
	{
	  if(n[0]*(Pt2[0]-pt[0])+n[1]*(Pt2[1]-pt[1])+n[2]*(Pt2[2]-pt[2])>0)
	  {
		Pt2[0]=(int) (pt[0] + r*n[0] + 0.5);
		Pt2[1]=(int) (pt[1] + r*n[1] + 0.5);
		Pt2[2]=(int) (pt[2] + r*n[2] + 0.5);
	  }
	  else
	  {
		Pt2[0]=(int) (pt[0] - r*n[0] + 0.5);
		Pt2[1]=(int) (pt[1] - r*n[1] + 0.5);
		Pt2[2]=(int) (pt[2] - r*n[2] + 0.5);
	  }
	}
  }


  //double longueur=sqrt((Pt1[0]-Pt2[0])*(Pt1[0]-Pt2[0])+(Pt1[1]-Pt2[1])*(Pt1[1]-Pt2[1])+(Pt1[2]-Pt2[2])*(Pt1[2]-Pt2[2]));
  //if (longueur<(double)2*r-2 )
	//fprintf(stderr, "longueur = %lf", longueur);
  //Bresenham algorithm
  switch (t) {
  default:
	if ( _VT_VERBOSE_ )
	  fprintf(stderr,"%s: such image type not handled yet\n", proc);
	return(0);
  case FLOAT:
	theBuf = (float *)buf;
	float v;
	
	//COPY
	
	ix = dim[0];
	ixy = dim[0] * dim[1];
	//size = (size_t)dim[0] * (size_t)dim[1] * (size_t)dim[2];

	dx = Pt2[0] - Pt1[0];
	dy = Pt2[1] - Pt1[1];
	dz = Pt2[2] - Pt1[2];

	ax = ABS(dx) << 1;
	ay = ABS(dy) << 1;
	az = ABS(dz) << 1;

	sx = ZSGN(dx);
	sy = ZSGN(dy);
	sz = ZSGN(dz);

	x = Pt1[0];
	y = Pt1[1];
	z = Pt1[2];

	v=1.0;
	
    /* x dominant
     */
    if (ax >= MAX(ay, az)) {
      yd = ay - (ax >> 1);
      zd = az - (ax >> 1);
      for (;;) {
		theBuf[z*ixy + y*ix + x] += v;
		if (x == Pt2[0]) {
		  break;
		}
	    if (yd >= 0) {
	      y += sy;
	      yd -= ax;
	    }
	    if (zd >= 0) {
	      z += sz;
	      zd -= ax;
	    }
		x += sx;
		yd += ay;
		zd += az;
      }
    }
    /* y dominant
     */
    else if (ay >= MAX(ax, az)) {
      xd = ax - (ay >> 1);
      zd = az - (ay >> 1);
      for (;;) {
		theBuf[z*ixy + y*ix + x] += v;
		if (y == Pt2[1]) {
	      break;
		}
		if (xd >= 0) {
	      x += sx;
	      xd -= ay;
		}
		if (zd >= 0) {
	      z += sz;
	      zd -= ay;
		}
		y += sy;
		xd += ax;
		zd += az;
      }
    }
    /* z dominant 
     */
    else if (az >= MAX(ax, ay)) {
      xd = ax - (az >> 1);
      yd = ay - (az >> 1);
      for (;;) {
		theBuf[z*ixy + y*ix + x] += v;
		if (z == Pt2[2]) {
	      break;
		}
		if (xd >= 0) {
	      x += sx;
	      xd -= az;
		}
	    if (yd >= 0) {
		  y += sy;
	      yd -= az;
		}
		z += sz;
		xd += ax;
		yd += ay;
      }
    }

	//END COPY
	
	break;
  }
  
  
  return(1);
}


int extractLabelFromLine(void *buf, int *dim, bufferType t, int *pt, double tht, double phi, int r, int *label, int *p1, int *p2)
{
  char *proc = "extractLabelFromLine";

  int Pt1[3],Pt2[3];
  int i, j;
  double n[3];
  double lambdas[6];
  double tmp;
  double maxdim;
  unsigned char *theBufU8=NULL;
  unsigned short int *theBufU16=NULL;
  
  //size_t size;
  int ix, ixy;
  int dx, dy, dz;  
  int ax, ay, az;
  int sx, sy, sz;
  int x, y, z;
  int xd, yd, zd;

  int r2=r*r;
  double d2;

  int l;
  
  int p1marker=0;

  
  // n : director vector of the line
  SphericalAnglesToUnitVector(tht, phi, n);
  
  // Computing the limit points Pt1, Pt2 of the image
  maxdim=(double)(dim[0]<dim[1]) ? ((dim[1]<dim[2]) ? dim[2] : dim[1]) : ((dim[0]<dim[2]) ? dim[2] : dim[0]);
  lambdas[0]=(n[0]==0) ? 2*maxdim : (0-(double)pt[0])/n[0];
  lambdas[1]=(n[1]==0) ? 2*maxdim : (0-(double)pt[1])/n[1];
  lambdas[2]=(n[2]==0) ? 2*maxdim : (0-(double)pt[2])/n[2];
  lambdas[3]=(n[0]==0) ? 2*maxdim : ((double)dim[0]-(double)pt[0]-1)/n[0];
  lambdas[4]=(n[1]==0) ? 2*maxdim : ((double)dim[1]-(double)pt[1]-1)/n[1];
  lambdas[5]=(n[2]==0) ? 2*maxdim : ((double)dim[2]-(double)pt[2]-1)/n[2];
  
    
  /* Sort the lambdas
   * 
  procédure tri_insertion
   * 
   */
  for (i=1; i<6; i++) {
	tmp=lambdas[i];
	j=i;
	while (j>0 && lambdas[j-1]>tmp) {
	  lambdas[j]=lambdas[j-1];
	  j--;
	}
	lambdas[j]=tmp;
  }
  
  // Pt1 
  i=0;
  while((pt[0]+lambdas[i]*n[0]+0.5 < 0 || pt[0]+lambdas[i]*n[0]+0.5 >= dim[0] ||
		 pt[1]+lambdas[i]*n[1]+0.5 < 0 || pt[1]+lambdas[i]*n[1]+0.5 >= dim[1] ||
		 pt[2]+lambdas[i]*n[2]+0.5 < 0 || pt[2]+lambdas[i]*n[2]+0.5 >= dim[2] ) && i<6 )
	i++;
  if (i>=6)
  {
	fprintf(stderr,"%s: error while computing the limit points: 1st point not found\n", proc);
	return(0);
  }
  //i1=i;
  Pt1[0]= pt[0] + (int) (lambdas[i]*n[0]+0.5); 
  Pt1[1]= pt[1] + (int) (lambdas[i]*n[1]+0.5); 
  Pt1[2]= pt[2] + (int) (lambdas[i]*n[2]+0.5); 
  
  // Pt2
  i=5;
  while((pt[0]+lambdas[i]*n[0]+0.5 < 0 || pt[0]+lambdas[i]*n[0]+0.5 >= dim[0] ||
		 pt[1]+lambdas[i]*n[1]+0.5 < 0 || pt[1]+lambdas[i]*n[1]+0.5 >= dim[1] ||
		 pt[2]+lambdas[i]*n[2]+0.5 < 0 || pt[2]+lambdas[i]*n[2]+0.5 >= dim[2] ) && i>=0 )
	i--;
  if (i<0)
  {
	fprintf(stderr,"%s: error while computing the limit points: 2nd point not found\n", proc);
	return(0);
  }
  //i2=i;
  Pt2[0]= pt[0] + (int) (lambdas[i]*n[0]+0.5); 
  Pt2[1]= pt[1] + (int) (lambdas[i]*n[1]+0.5); 
  Pt2[2]= pt[2] + (int) (lambdas[i]*n[2]+0.5); 

  
  
  if(r>0)				// Computing Pt1, Pt2 so that voting distance < r  
  {
	// Pt1 :
	
	d2=(Pt1[0]-pt[0])*(Pt1[0]-pt[0])+(Pt1[1]-pt[1])*(Pt1[1]-pt[1])+(Pt1[2]-pt[2])*(Pt1[2]-pt[2]);
	
	if (r2<d2)
	{
	  if(n[0]*(Pt1[0]-pt[0])+n[1]*(Pt1[1]-pt[1])+n[2]*(Pt1[2]-pt[2])>0)
	  {
		Pt1[0]=(int) (pt[0] + r*n[0] + 0.5);
		Pt1[1]=(int) (pt[1] + r*n[1] + 0.5);
		Pt1[2]=(int) (pt[2] + r*n[2] + 0.5);
	  }
	  else
	  {
		Pt1[0]=(int) (pt[0] - r*n[0] + 0.5);
		Pt1[1]=(int) (pt[1] - r*n[1] + 0.5);
		Pt1[2]=(int) (pt[2] - r*n[2] + 0.5);
	  }
	}
	
	// Pt2 :

	d2=(Pt2[0]-pt[0])*(Pt2[0]-pt[0])+(Pt2[1]-pt[1])*(Pt2[1]-pt[1])+(Pt2[2]-pt[2])*(Pt2[2]-pt[2]);
	
	if (r2<d2)
	{
	  if(n[0]*(Pt2[0]-pt[0])+n[1]*(Pt2[1]-pt[1])+n[2]*(Pt2[2]-pt[2])>0)
	  {
		Pt2[0]=(int) (pt[0] + r*n[0] + 0.5);
		Pt2[1]=(int) (pt[1] + r*n[1] + 0.5);
		Pt2[2]=(int) (pt[2] + r*n[2] + 0.5);
	  }
	  else
	  {
		Pt2[0]=(int) (pt[0] - r*n[0] + 0.5);
		Pt2[1]=(int) (pt[1] - r*n[1] + 0.5);
		Pt2[2]=(int) (pt[2] - r*n[2] + 0.5);
	  }
	}
  }


  //double longueur=sqrt((Pt1[0]-Pt2[0])*(Pt1[0]-Pt2[0])+(Pt1[1]-Pt2[1])*(Pt1[1]-Pt2[1])+(Pt1[2]-Pt2[2])*(Pt1[2]-Pt2[2]));
  //if (longueur<(double)2*r-2 )
	//fprintf(stderr, "longueur = %lf", longueur);
  //Bresenham algorithm
  switch (t) {
  default:
	if ( _VT_VERBOSE_ )
	  fprintf(stderr,"%s: such image type not handled yet\n", proc);
	return(0);
  case UCHAR:
  case SCHAR:
	  theBufU8 = (unsigned char *)buf;
//	float v;
	  break;
	case USHORT:
	case SSHORT:
	  theBufU16 = (unsigned short int *)buf;
	  break;
	}
	//COPY
	
	ix = dim[0];
	ixy = dim[0] * dim[1];
	//size = (size_t)dim[0] * (size_t)dim[1] * (size_t)dim[2];

	dx = Pt2[0] - Pt1[0];
	dy = Pt2[1] - Pt1[1];
	dz = Pt2[2] - Pt1[2];

	ax = ABS(dx) << 1;
	ay = ABS(dy) << 1;
	az = ABS(dz) << 1;

	sx = ZSGN(dx);
	sy = ZSGN(dy);
	sz = ZSGN(dz);

	x = Pt1[0];
	y = Pt1[1];
	z = Pt1[2];
	
	

//	v=1.0;
	l=0;
    /* x dominant
     */
    if (ax >= MAX(ay, az)) {
      yd = ay - (ax >> 1);
      zd = az - (ax >> 1);
      for (;;) {
    switch(t) {
    case UCHAR:
    case SCHAR:
		  if ((int)theBufU8[z*ixy + y*ix + x] != 0)
		  {
		    if(p1marker==0) {
		      p1marker=1;
		      p1[0]=x; p1[1]=y; p1[2]=z;
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }
		    else {
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }
		    
		    if(l==0) l=(int)theBufU8[z*ixy + y*ix + x];
		    else {
		      if ((int)theBufU8[z*ixy + y*ix + x]!=l)
		        l=-1;
		    }
		  }
		  break;
	  case USHORT:
    case SSHORT:
		  if ((int)theBufU16[z*ixy + y*ix + x] != 0)
		  {
		    if(p1marker==0) {
		      p1marker=1;
		      p1[0]=x; p1[1]=y; p1[2]=z;
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }
		    else {
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }

		    if(l==0) l=(int)theBufU16[z*ixy + y*ix + x];
		    else {
		      if ((int)theBufU16[z*ixy + y*ix + x]!=l)
		        l=-1;
		    }
		  }
		  break;
	  default:
    	if ( _VT_VERBOSE_ )
	    fprintf(stderr,"%s: such image type not handled yet\n", proc);
	    return(0);
		}
    if(l<0) break;
    
		if (x == Pt2[0]) {
		  break;
		}
	    if (yd >= 0) {
	      y += sy;
	      yd -= ax;
	    }
	    if (zd >= 0) {
	      z += sz;
	      zd -= ax;
	    }
		x += sx;
		yd += ay;
		zd += az;
      }
    }
    /* y dominant
     */
    else if (ay >= MAX(ax, az)) {
      xd = ax - (ay >> 1);
      zd = az - (ay >> 1);
      for (;;) {
    switch(t) {
    case UCHAR:
    case SCHAR:
		  if ((int)theBufU8[z*ixy + y*ix + x] != 0)
		  {
		    if(p1marker==0) {
		      p1marker=1;
		      p1[0]=x; p1[1]=y; p1[2]=z;
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }
		    else {
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }

		    if(l==0) l=(int)theBufU8[z*ixy + y*ix + x];
		    else {
		      if ((int)theBufU8[z*ixy + y*ix + x]!=l)
		        l=-1;
		    }
		  }
		  break;
	  case USHORT:
    case SSHORT:
		  if ((int)theBufU16[z*ixy + y*ix + x] != 0)
		  {
		    if(p1marker==0) {
		      p1marker=1;
		      p1[0]=x; p1[1]=y; p1[2]=z;
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }
		    else {
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }

		    if(l==0) l=(int)theBufU16[z*ixy + y*ix + x];
		    else {
		      if ((int)theBufU16[z*ixy + y*ix + x]!=l)
		        l=-1;
		    }
		  }
		  break;
	  default:
    	if ( _VT_VERBOSE_ )
	    fprintf(stderr,"%s: such image type not handled yet\n", proc);
	    return(0);
		}
    if(l<0) break;

		if (y == Pt2[1]) {
	      break;
		}
		if (xd >= 0) {
	      x += sx;
	      xd -= ay;
		}
		if (zd >= 0) {
	      z += sz;
	      zd -= ay;
		}
		y += sy;
		xd += ax;
		zd += az;
      }
    }
    /* z dominant 
     */
    else if (az >= MAX(ax, ay)) {
      xd = ax - (az >> 1);
      yd = ay - (az >> 1);
      for (;;) {
    switch(t) {
    case UCHAR:
    case SCHAR:
		  if ((int)theBufU8[z*ixy + y*ix + x] != 0)
		  {
		    if(p1marker==0) {
		      p1marker=1;
		      p1[0]=x; p1[1]=y; p1[2]=z;
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }
		    else {
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }

		    if(l==0) l=(int)theBufU8[z*ixy + y*ix + x];
		    else {
		      if ((int)theBufU8[z*ixy + y*ix + x]!=l)
		        l=-1;
		    }
		  }
		  break;
	  case USHORT:
    case SSHORT:
		  if ((int)theBufU16[z*ixy + y*ix + x] != 0)
		  {
		    if(p1marker==0) {
		      p1marker=1;
		      p1[0]=x; p1[1]=y; p1[2]=z;
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }
		    else {
		      p2[0]=x; p2[1]=y; p2[2]=z;
		    }

		    if(l==0) l=(int)theBufU16[z*ixy + y*ix + x];
		    else {
		      if ((int)theBufU16[z*ixy + y*ix + x]!=l)
		        l=-1;
		    }
		  }
		  break;
	  default:
    	if ( _VT_VERBOSE_ )
	    fprintf(stderr,"%s: such image type not handled yet\n", proc);
	    return(0);
		}
    if(l<0) break;

		if (z == Pt2[2]) {
	      break;
		}
		if (xd >= 0) {
	      x += sx;
	      xd -= az;
		}
	    if (yd >= 0) {
		  y += sy;
	      yd -= az;
		}
		z += sz;
		xd += ax;
		yd += ay;
      }
    }

	//END COPY
   if(p1marker==0) {
     p1[0]=pt[0]; p1[1]=pt[1]; p1[2]=pt[2];
     p2[0]=pt[0]; p2[1]=pt[1]; p2[2]=pt[2];
   }
	
  *label=l;
  return(1);
}

int extractLabelFromCone(void ***array, int *dim, bufferType t, int *pt, double tht, double phi, double r, double alpha, int *label, int *p1, int *p2)
{
  char *proc = "extractLabelFromCone";
  
  
  //int dimFields[3];
  

  //double **angles;
  //int Nangles;

  //double theCoeff[3];
  //theCoeff[0] = theCoeff[1] = scale;
  //theCoeff[2] = zfact*scale;
  //dimFields[0] = dimFields[1] = dimFields[2] = (int)(2*r+1);

  int x,y,z,i;
  double norme, norme2, dot;
  double cosalpha;
  double n[3];
  double r2=r*r;
  
  unsigned char*** resU8=NULL;
  unsigned short int*** resU16=NULL;
  
  int l;
  double d1=0, d2=0;
  //unsigned char ***res8;
  

  switch (t)
  {
  //case UCHAR :
  //  res8 = (unsigned char***)array;
  //  break;
  case UCHAR :
  case SCHAR :
    resU8 = (unsigned char***)array;
    break;
  case USHORT :
  case SSHORT :
    resU16 = (unsigned short int***)array;
    break;
  default:
    fprintf(stderr, "%s: array type not handled yet...\n", proc);
    return(0);
  }
  
  cosalpha=cos(alpha); 
  if (cosalpha<0) cosalpha*=-1;
  SphericalAnglesToUnitVector( tht, phi, n);
  
  int xmin, xmax, ymin, ymax, zmin, zmax;
  xmin=(pt[0]<r)? 0: pt[0]-r;
  ymin=(pt[1]<r)? 0: pt[1]-r;
  zmin=(pt[2]<r)? 0: pt[2]-r;
  xmax=(pt[0]>=dim[0]-r)? dim[0]-1: pt[0]+r;
  ymax=(pt[1]>=dim[1]-r)? dim[1]-1: pt[1]+r;
  zmax=(pt[2]>=dim[2]-r)? dim[2]-1: pt[2]+r;
  
  //fprintf(stdout, "pt     = \t[%d %d %d]\t", pt[0], pt[1], pt[2]);
  //fprintf(stdout, "dim    = \t[%d %d %d]\n", dim[0], dim[1], dim[2]);
  //fprintf(stdout, "bornes = \t{%d %d}\t{%d %d}\t{%d %d}\n", xmin, xmax, ymin, ymax, zmin, zmax);
  
  l=0;
  for (i=0; i<3; i++) {p1[i]=pt[i]; p2[i]=pt[i];};

  
  for (z=zmin ; z<zmax ; z++)
  {
    for (y=ymin ; y<ymax ; y++)
    {
      for (x=xmin ; x<xmax ; x++)
      {
      	//fprintf(stdout, "{x y z} = {%d %d %d}\t", x, y, z);
      	norme2=(x-pt[0])*(x-pt[0])+(y-pt[1])*(y-pt[1])+(z-pt[2])*(z-pt[2]);
      	if (norme2>r2) continue;
      	norme=sqrt(norme2);
      	dot=(x-pt[0])*n[0]+(y-pt[1])*n[1]+(z-pt[2])*n[2];
      	if (dot<0) dot=-dot;
      	if (dot<cosalpha*norme) continue;
        switch (t) {
        case UCHAR:
        case SCHAR:
          if ((int)resU8[z][y][x]!=0) {
            if ((x-pt[0])*n[0]+(y-pt[1])*n[1]+(z-pt[2])*n[2]>0) {
              if(norme>d1) {
                p1[0]=x;
                p1[1]=y;
                p1[2]=z;
                d1=norme;
              }
            }
            else {
              if(norme>d2) {
                p2[0]=x;
                p2[1]=y;
                p2[2]=z;
                d2=norme;
              }
            }
          }
          if ((int)resU8[z][y][x]==0) continue;
          if ((int)resU8[z][y][x]==0) continue;
          l=(l==0)? (int)resU8[z][y][x]:-1;
          break;
        case USHORT:
        case SSHORT:
          if ((int)resU16[z][y][x]!=0) {
            if ((x-pt[0])*n[0]+(y-pt[1])*n[1]+(z-pt[2])*n[2]>0) {
              if(norme>d1) {
                p1[0]=x;
                p1[1]=y;
                p1[2]=z;
                d1=norme;              }
            }
            else {
              if(norme>d2) {
                p2[0]=x;
                p2[1]=y;
                p2[2]=z;
                d2=norme;              
              }
            }
          }
          if ((int)resU16[z][y][x]==0) continue;
          if ((int)resU16[z][y][x]==l) continue;
          l=(l==0)? (int)resU16[z][y][x]:-1;
          break;
        default:
          fprintf(stderr, "%s: array type not handled yet...\n", proc);
          return(0);
        }  
      	if (l<0) break;
      }
      if (l<0) break;
    }
    if (l<0) break;
  }
  if(d1==0)
  {
    p1[0]=p2[0];
    p1[1]=p2[1];
    p1[2]=p2[2];
  }
  if(d2==0)
  {
    p2[0]=p1[0];
    p2[1]=p1[1];
    p2[2]=p1[2];
  }

  *label=l;
  return(1);
}

int setLineToLabel(void *buf, int *dim, bufferType t, int *Pt1, int *Pt2, int label)
{
  char *proc = "setLineToLabel";

//  int Pt1[3],Pt2[3];
  //int i, j;
  //double n[3];
  //double lambdas[6];
  //double tmp;
  //double maxdim;
  unsigned char *theBufU8=NULL;
  unsigned short int *theBufU16=NULL;
  
  //size_t size;
  int ix, ixy;
  int dx, dy, dz;  
  int ax, ay, az;
  int sx, sy, sz;
  int x, y, z;
  int xd, yd, zd;

  //int r2=r*r;
  //double d2;
  

  //double longueur=sqrt((Pt1[0]-Pt2[0])*(Pt1[0]-Pt2[0])+(Pt1[1]-Pt2[1])*(Pt1[1]-Pt2[1])+(Pt1[2]-Pt2[2])*(Pt1[2]-Pt2[2]));
  //if (longueur<(double)2*r-2 )
	//fprintf(stderr, "longueur = %lf", longueur);
  //Bresenham algorithm
  switch (t) {
  default:
	if ( _VT_VERBOSE_ )
	  fprintf(stderr,"%s: such image type not handled yet\n", proc);
	return(0);
  case UCHAR:
  case SCHAR:
  	theBufU8 = (unsigned char *)buf;
	  break;
	case USHORT:
  case SSHORT:
  	theBufU16 = (unsigned short int *)buf;
	  break;
	}
	//float v;
	
	//COPY
	
	ix = dim[0];
	ixy = dim[0] * dim[1];
	//size = (size_t)dim[0] * (size_t)dim[1] * (size_t)dim[2];

	dx = Pt2[0] - Pt1[0];
	dy = Pt2[1] - Pt1[1];
	dz = Pt2[2] - Pt1[2];

	ax = ABS(dx) << 1;
	ay = ABS(dy) << 1;
	az = ABS(dz) << 1;

	sx = ZSGN(dx);
	sy = ZSGN(dy);
	sz = ZSGN(dz);

	x = Pt1[0];
	y = Pt1[1];
	z = Pt1[2];

	//v=1.0;
	
    /* x dominant
     */
    if (ax >= MAX(ay, az)) {
      yd = ay - (ax >> 1);
      zd = az - (ax >> 1);
      for (;;) {
    switch (t) {
    case UCHAR:
    case SCHAR:
      if ((int)theBufU8[z*ixy + y*ix + x] == 0)
    		theBufU8[z*ixy + y*ix + x] = (unsigned char) label;
      if ((int)theBufU8[z*ixy + y*ix + x] != label)
        theBufU8[z*ixy + y*ix + x] = (unsigned char) 255;
      break;
    case USHORT:
    case SSHORT:
      if ((int)theBufU16[z*ixy + y*ix + x] == 0)
    		theBufU16[z*ixy + y*ix + x] = (unsigned short int) label;
      if ((int)theBufU16[z*ixy + y*ix + x] != label)
        theBufU16[z*ixy + y*ix + x] = (unsigned short int) 32767;
      break;
    default:
    	if ( _VT_VERBOSE_ )
    	  fprintf(stderr,"%s: such image type not handled yet\n", proc);
    	return(0);
    }
		if (x == Pt2[0]) {
		  break;
		}
	    if (yd >= 0) {
	      y += sy;
	      yd -= ax;
	    }
	    if (zd >= 0) {
	      z += sz;
	      zd -= ax;
	    }
		x += sx;
		yd += ay;
		zd += az;
      }
    }
    /* y dominant
     */
    else if (ay >= MAX(ax, az)) {
      xd = ax - (ay >> 1);
      zd = az - (ay >> 1);
      for (;;) {
    switch (t) {
    case UCHAR:
    case SCHAR:
      if ((int)theBufU8[z*ixy + y*ix + x] == 0)
    		theBufU8[z*ixy + y*ix + x] = (unsigned char) label;
      if ((int)theBufU8[z*ixy + y*ix + x] != label)
        theBufU8[z*ixy + y*ix + x] = (unsigned char) 255;
      break;
    case USHORT:
    case SSHORT:
      if ((int)theBufU16[z*ixy + y*ix + x] == 0)
    		theBufU16[z*ixy + y*ix + x] = (unsigned short int) label;
      if ((int)theBufU16[z*ixy + y*ix + x] != label)
        theBufU16[z*ixy + y*ix + x] = (unsigned short int) 32767;
      break;
    default:
    	if ( _VT_VERBOSE_ )
    	  fprintf(stderr,"%s: such image type not handled yet\n", proc);
    	return(0);
    }
		if (y == Pt2[1]) {
	      break;
		}
		if (xd >= 0) {
	      x += sx;
	      xd -= ay;
		}
		if (zd >= 0) {
	      z += sz;
	      zd -= ay;
		}
		y += sy;
		xd += ax;
		zd += az;
      }
    }
    /* z dominant 
     */
    else if (az >= MAX(ax, ay)) {
      xd = ax - (az >> 1);
      yd = ay - (az >> 1);
      for (;;) {
    switch (t) {
    case UCHAR:
    case SCHAR:
      if ((int)theBufU8[z*ixy + y*ix + x] == 0)
    		theBufU8[z*ixy + y*ix + x] = (unsigned char) label;
      if ((int)theBufU8[z*ixy + y*ix + x] != label)
        theBufU8[z*ixy + y*ix + x] = (unsigned char) 255;
      break;
    case USHORT:
    case SSHORT:
      if ((int)theBufU16[z*ixy + y*ix + x] == 0)
    		theBufU16[z*ixy + y*ix + x] = (unsigned short int) label;
      if ((int)theBufU16[z*ixy + y*ix + x] != label)
        theBufU16[z*ixy + y*ix + x] = (unsigned short int) 32767;
      break;
    default:
    	if ( _VT_VERBOSE_ )
    	  fprintf(stderr,"%s: such image type not handled yet\n", proc);
    	return(0);
    }
		if (z == Pt2[2]) {
	      break;
		}
		if (xd >= 0) {
	      x += sx;
	      xd -= az;
		}
	    if (yd >= 0) {
		  y += sy;
	      yd -= az;
		}
		z += sz;
		xd += ax;
		yd += ay;
      }
    }

	//END COPY
	
	//break;
  //}
  
  
  return(1);
}

int setConeToLabel(void ***array, int *dim, bufferType t, int *p1, int *p2, double tht, double phi, double alpha, int l)
{
  char *proc = "setConeToLabel";
  

  //int dimFields[3];
  

  //double **angles;
  //int Nangles;

  //double theCoeff[3];
  //theCoeff[0] = theCoeff[1] = scale;
  //theCoeff[2] = zfact*scale;
  //dimFields[0] = dimFields[1] = dimFields[2] = (int)(2*r+1);

  int x,y,z;
  double norme, norme2, dot;
  double cosalpha;
  double n[3];
  double r=sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
  double r2=r*r;
  double v[3];
  v[0]=p2[0]-p1[0];
  v[1]=p2[1]-p1[1];
  v[2]=p2[2]-p1[2];
    
  unsigned char*** resU8=NULL;
  unsigned short int*** resU16=NULL;
  

  //unsigned char ***res8;
  

  switch (t)
  {
  //case UCHAR :
  //  res8 = (unsigned char***)array;
  //  break;
  case UCHAR :
  case SCHAR :
    resU8 = (unsigned char***)array;
    break;
  case USHORT :
  case SSHORT :
    resU16 = (unsigned short int***)array;
    break;
  default:
    fprintf(stderr, "%s: array type not handled yet...\n", proc);
    return(0);
  }
  
  cosalpha=cos(alpha); 
  if (cosalpha<0) cosalpha*=-1;
  SphericalAnglesToUnitVector( tht, phi, n);
  
  int xmin, xmax, ymin, ymax, zmin, zmax;
  xmin=(p1[0]<r)? 0: p1[0]-r;
  ymin=(p1[1]<r)? 0: p1[1]-r;
  zmin=(p1[2]<r)? 0: p1[2]-r;
  xmax=(p1[0]>=dim[0]-r)? dim[0]-1: p1[0]+r;
  ymax=(p1[1]>=dim[1]-r)? dim[1]-1: p1[1]+r;
  zmax=(p1[2]>=dim[2]-r)? dim[2]-1: p1[2]+r;
  
  //fprintf(stdout, "pt     = \t[%d %d %d]\t", pt[0], pt[1], pt[2]);
  //fprintf(stdout, "dim    = \t[%d %d %d]\n", dim[0], dim[1], dim[2]);
  //fprintf(stdout, "bornes = \t{%d %d}\t{%d %d}\t{%d %d}\n", xmin, xmax, ymin, ymax, zmin, zmax);
  
  for (z=zmin ; z<zmax ; z++)
  {
    for (y=ymin ; y<ymax ; y++)
    {
      for (x=xmin ; x<xmax ; x++)
      {
        if(v[0]*(x-p1[0])+v[1]*(y-p1[1])+v[2]*(z-p1[2])<0) continue;
      	//fprintf(stdout, "{x y z} = {%d %d %d}\t", x, y, z);
      	norme2=(x-p1[0])*(x-p1[0])+(y-p1[1])*(y-p1[1])+(z-p1[2])*(z-p1[2]);
      	if (norme2>r2) continue;
      	norme=sqrt(norme2);
      	dot=(x-p1[0])*n[0]+(y-p1[1])*n[1]+(z-p1[2])*n[2];
      	if (dot<0) dot=-dot;
      	if (dot<cosalpha*norme) continue;
        switch (t) {
        case UCHAR:
        case SCHAR:
          if ((int)resU8[z][y][x]==0) resU8[z][y][x]=(unsigned char) l;
          if ((int)resU8[z][y][x]!=l) resU8[z][y][x]=(unsigned char) 255;
          break;
        case USHORT:
        case SSHORT:
          if ((int)resU16[z][y][x]==0) resU16[z][y][x]=(unsigned short int) l;
          if ((int)resU16[z][y][x]!=l) resU16[z][y][x]=(unsigned short int) 32767;
          break;
        default:
          fprintf(stderr, "%s: array type not handled yet...\n", proc);
          return(0);
        }  
      }
    }
  }
  
  return(1);
}
