/*************************************************************************
 * mt_membrane3D.c -
 *
 * $Id: mt_membrane3D.c,v 2.0 2013/10/22 14:22:34 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 * 
 * CREATION DATE: 
 * 2013/06/20 
 *
 * ADDITIONS, CHANGES
 *
 *
 */
 


#include <mt_membrane3D.h>
#include <transfo.h>
#include <eigens.h>
#include <linearFiltering.h>

#define FLTZERO 1e-8

/* find maximum of a and b */
#define MAX(a,b) (((a)>(b))?(a):(b))

/* absolute value of a */
#define ABS(a) (((a)<0) ? -(a) : (a))

/* take sign of a, either -1, 0, or 1 */
#define ZSGN(a) (((a)<0) ? -1 : (a)>0 ? 1 : 0)

static int _borderLength_[3] = { 20, 20, 20 };
static int _verbose_ = 1;
//static int _debug_ = 1;



int MT_Compute3DMultiScale( vt_image *theIm,
			    vt_3Dimres *imsRes,
			    double scale1,
			    double scale2,
			    int nbscales,
			    double zfact,
			    enumStructureColor color,
			    enumMode mode,
			    int hsp)
{
  char *proc = "MT_Compute3DMultiScale";
  int s, nbs = nbscales;
  double scale;

  vt_3Dimages ims3D;
  vt_2Dimauxs ims2D;

  double tau = 1;


  float theCoeff[3];

  int slice;

  typeResponseInSlice aux2D;

  int x, y, z;
  float ***maxRep = (float***)imsRes->imRep.array;
  float ***maxTht = (float***)imsRes->imTheta.array;
  float ***maxPhi = (float***)imsRes->imPhi.array;
  float ***maxScl = (float***)imsRes->imScale.array;

  float ***theRep;
  float ***theTht;
  float ***thePhi;

  double newRep;



  if ( VT_Alloc3Dimages( &ims3D, theIm, proc ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate 3D auxiliary images\n", proc );
      return( -1 );
    }
  }


  theRep = (float***)ims3D.imzz.array;
  theTht = (float***)ims3D.tmp0.array;
  thePhi = (float***)ims3D.tmp1.array;


  if ( VT_Alloc2Dimauxs( &ims2D, theIm, proc ) != 1 ) {
    VT_Free3DImages( &ims3D );
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate 2D auxiliary images\n", proc );
      return( -1 );
    }
  }






  if ( scale2 <= scale1 ) nbs = 1;

  for ( s = 0; s < nbs; s ++ ) {

    if ( s == 0 ) {
      scale = scale1;
    } else {
      scale = exp( ((double)(nbs-1-s)*log(scale1) + (double)(s)*log(scale2)) /
		   (double)(nbs-1) );
    }

    if ( _verbose_ ) {
      fprintf( stdout, "... processing scale #%d, sigma = %f\n", s, scale );
    }



    theCoeff[0] = theCoeff[1] = scale;
    theCoeff[2] = zfact * scale;

    /* Calcul des images 3D suivantes
       tmp0 : . . 0
       tmp1 : . . 1
       imx  : 1 0 0
       imy  : 0 1 0
       imz  : 0 0 1
       imzz : 0 0 2
    */

    if ( _verbose_ ) {
      fprintf( stdout, "\tPartial filtering of the whole 3D image...\n");
    }


    if ( MT_Filter3Dimages( theIm, &ims3D, theCoeff ) != 1 ) {
      VT_Free2Dimauxs( &ims2D );
      VT_Free3DImages( &ims3D );
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to 3D filter\n", proc );
    	return( -1 );
    }



    aux2D.theX  = (float***)(ims3D.imx.array);
    aux2D.theY  = (float***)(ims3D.imy.array);
    aux2D.theZ  = (float***)(ims3D.imz.array);



    for ( slice = 0; slice < theIm->dim.z; slice ++ ) {

        if ( _verbose_ && (theIm->dim.z<=30 || (slice%(theIm->dim.z/30) ==0))) {
          fprintf( stdout, "\tTreating slice #%d/%d\t(%d%%)\n",
		   slice, (int)theIm->dim.z, (int)(100*(float)slice/(float)theIm->dim.z));
        }


      /* Calcul des images 2D suivantes
	 imxx : 2 0 0
	 imyy : 0 2 0
	 imxy : 1 1 0
	 imxz : 1 0 1
	 imyz : 0 1 1
      */

      if ( MT_Filter2Dimauxs( &ims3D, &ims2D, theCoeff, slice ) != 1 ) {
	VT_Free2Dimauxs( &ims2D );
	VT_Free3DImages( &ims3D );
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to 2D filter\n", proc );
	return( -1 );
	
      }

      aux2D.theXX = ((float***)(ims2D.imxx.array))[0];
      aux2D.theYY = ((float***)(ims2D.imyy.array))[0];
      aux2D.theXY = ((float***)(ims2D.imxy.array))[0];

      aux2D.theXZ = ((float***)(ims2D.imxz.array))[0];
      aux2D.theYZ = ((float***)(ims2D.imyz.array))[0];

      aux2D.theZZ = ((float***)(ims3D.imzz.array))[slice];
      
      aux2D.theRep   = theRep[slice];
      aux2D.theTheta = theTht[slice];
      aux2D.thePhi   = thePhi[slice];

      MT_Compute3DresponseInSlice( &aux2D,  
				   theIm->dim.x, theIm->dim.y,                   
				   theIm->dim.z, slice, tau, theCoeff, color,
				   mode, hsp);
      

    } /* for ( slice = 0; slice < theIm->dim.z; slice ++ ) */

    if ( _verbose_ ) {
      fprintf( stdout, "\tMultiscale conversion of scale %f...\n", scale );
    }


    /* on recopie le resultat dans les images resultats
     */
    
    if ( s == 0 ) {
      for ( z = 0; z < theIm->dim.z; z ++ )
      for ( y = 0; y < theIm->dim.y; y ++ )
      for ( x = 0; x < theIm->dim.x; x ++ ) {
	if ( theRep[z][y][x] > 0.0 ) {

	  maxRep[z][y][x] = theRep[z][y][x];

	  maxTht[z][y][x] = theTht[z][y][x];
	  maxPhi[z][y][x] = thePhi[z][y][x];
	  maxScl[z][y][x] = scale * theIm->siz.x;
	} else {
	  maxRep[z][y][x] = maxTht[z][y][x] = maxPhi[z][y][x] = 0.0;
	  maxScl[z][y][x] = 0.0;
	}
      }
    } else {
      for ( z = 0; z < theIm->dim.z; z ++ )
      for ( y = 0; y < theIm->dim.y; y ++ )
      for ( x = 0; x < theIm->dim.x; x ++ ) {
        newRep = theRep[z][y][x];

        if ( newRep > maxRep[z][y][x] ) {
	  maxRep[z][y][x] = newRep;
	  maxTht[z][y][x] = theTht[z][y][x];
	  maxPhi[z][y][x] = thePhi[z][y][x];
	  maxScl[z][y][x] = scale * theIm->siz.x;
	}
      }
    }

  } /* for ( s = 0; s < nbs; s ++ ) */


  VT_Free3DImages( &ims3D );
  VT_Free2Dimauxs( &ims2D );
  return( 1 );
}
			   












/* calcul des extrema 3D pour des membranes
 */

void MT_Compute3DExtrema( vt_3Dimres *imsRes,
                          double zfact,
			  vt_image *imExt )
{
  int x, y, z;
  
  float ***theExt = (float***)(imExt->array);
  float ***theRep = (float***)(imsRes->imRep.array);
  float ***theTht = (float***)(imsRes->imTheta.array);
  float ***thePhi = (float***)(imsRes->imPhi.array);

  double v1[3], v2[3], v3[3];
  
  double rep, rn[2];

  double theta, phi;

  double vx=0.0, vy=0.0, vz=0.0;
  double norm;
  

  double rx, ry, rz;
  double dx, dy, dz;
  int j, ix, iy, iz;
  double dxdy, dxdz, dydz, dxdydz;
  double v4, v5, v6;




  for ( z=0; z<imExt->dim.z; z++ ) 
  for ( y=0; y<imExt->dim.y; y++ ) 
  for ( x=0; x<imExt->dim.x; x++ ) {

    theta = (double)theTht[z][y][x];
    phi   = (double)thePhi[z][y][x];
    
    theExt[z][y][x] = 0.0;
    if ( theRep[z][y][x] <= 0.0 ) continue;

    rep = theRep[z][y][x];


    /* v3 est le vecteur normal a la membrane
       v1 et v2 sont deux du plan tangent
    */
    SphericalAnglesToUnitsVectors( theta, phi, v3, v2, v1 );
    
    for (j=0;j<2;j++) {
      switch(j) {
      default:
      case 0 :
  vx = zfact*v3[0];    vy = zfact*v3[1];    vz = v3[2]; break;
      case 1 : 
  vx = -zfact*v3[0];  vy = -zfact*v3[1];  vz = -v3[2];  break;
      }
      norm=sqrt(vx*vx+vy*vy+vz*vz);
      vx/=norm; vy/=norm; vz/=norm;
      
      rx = x + vx;
      ix = (int)rx;
      if ( rx <= 0.0 || ix >= imExt->dim.x-1 ) {
	j = 3;
	continue;
      }
      ry = y + vy;
      iy = (int)ry;
      if ( ry <= 0.0 || iy >= imExt->dim.y-1 ) {
	j = 3;
	continue;
      }
      rz = z + vz;
      iz = (int)rz;
      if ( rz <= 0.0 || iz >= imExt->dim.z-1 ) {
	j = 3;
	continue;
      }      
      

      dx = rx - ix;
      dy = ry - iy;
      dz = rz - iz;

      dxdy = dx*dy;
      dxdz = dx*dz;
      dydz = dy*dz;
      dxdydz = dxdy*dz;

      v6 = dxdz-dxdydz;
      v5 = dxdy-dxdydz;
      v4 = dx-dxdy-v6;

      rn[j] = 0;
      rn[j] += dxdydz        * theRep[iz+1][iy+1][ix+1];
      rn[j] += (dydz-dxdydz) * theRep[iz+1][iy+1][ix  ];
      rn[j] += v6            * theRep[iz+1][iy  ][ix+1];
      rn[j] += (dz-dydz-v6)  * theRep[iz+1][iy  ][ix  ];
      rn[j] += v5            * theRep[iz  ][iy+1][ix+1];
      rn[j] += (dy-dydz-v5)  * theRep[iz  ][iy+1][ix  ];
      rn[j] += v4            * theRep[iz  ][iy  ][ix+1];
      rn[j] += (1-dy-dz+dydz-v4) * theRep[iz  ][iy][ix];

      if ( rn[j] > rep ) {
        j = 3;
        continue;
      }

    } /* fin de la boucle suivant la normale à la membrane */

    if ( rn[0] == rep && rn[1] == rep ) {
      continue;
    }

    if ( j < 3 ) theExt[z][y][x] = rep;
  } //x,y,z
}




















void MT_Compute3DresponseInSlice( typeResponseInSlice *aux,
				  int dimx,
				  int dimy,
				  int dimz,
				  int slice,
				  double tau,
				  float *theCoeff,
				  enumStructureColor color,
				  enumMode mode,
				  int hsp)
{
  float **theXX = aux->theXX;
  float **theYY = aux->theYY;
  float **theZZ = aux->theZZ;
  float **theXY = aux->theXY;
  float **theXZ = aux->theXZ;
  float **theYZ = aux->theYZ;
  
  float ***theX = aux->theX;
  float ***theY = aux->theY;
  float ***theZ = aux->theZ;

  float **theRep   = aux->theRep;
  float **theTheta = aux->theTheta;
  float **thePhi   = aux->thePhi;

  int x, y;

  double HESSIEN[9];
  double VALPROP[3];
  double VECPROP[9];

  double v[3];

  //double rapportXZ = theCoeff[0]/theCoeff[2];
  //double rapportXY = theCoeff[0]/theCoeff[1];
  //double norm;

//  case MODELBASED :
  //int n;
  int a, a1, a2;
  double val;
  double qx, qy, qz;
  double vx, vy, vz;
  double wx, wy, wz;
  double wwx, wwy, wwz;
  double gx, gy, gz;
  double rx, ry, rz;
  double dx, dy, dz;
  int ix, iy, iz;
  double dxdy, dxdz, dydz, dxdydz;
  double v4, v5, v6;

  int nbPosPts;
  int nbNegPts;
  int nbValPts;

  double sumValPts;
  double sumPosPts;
  double sumNegPts;

  double radii[3];
  radii[0]=tau*(double)theCoeff[0];
  radii[1]=tau*(double)theCoeff[1];
  radii[2]=tau*(double)theCoeff[2];
  //double sigma = theCoeff[0];
  //double radius = tau*sigma;

//  case ACME :
  double A, B, S;
  double T0, T1, T2, T3;
  double gamma = 1, alpha = 0.5, beta = 0.5;
  double c = 0.01;
  
  
  double theta, phi;
  
  double R;


  for ( y = 0; y < dimy; y ++ )
  for ( x = 0; x < dimx; x ++ ) {

    /* calcul des valeurs et vecteurs propres
     */

    HESSIEN[0] =              (double)(theXX[y][x]*theCoeff[0]*theCoeff[0]);
    HESSIEN[1] = HESSIEN[3] = (double)(theXY[y][x]*theCoeff[0]*theCoeff[1]);
    HESSIEN[2] = HESSIEN[6] = (double)(theXZ[y][x]*theCoeff[0]*theCoeff[2]);
    HESSIEN[4] =              (double)(theYY[y][x]*theCoeff[1]*theCoeff[1]);
    HESSIEN[5] = HESSIEN[7] = (double)(theYZ[y][x]*theCoeff[1]*theCoeff[2]);
    HESSIEN[8] =              (double)(theZZ[y][x]*theCoeff[2]*theCoeff[2]);

    
    theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0;



    if ( _ComputeEigensOfSymetricSquareMatrix( HESSIEN, VALPROP, VECPROP, 3 )
         != 1 ) {
      theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0.0;
      continue;
    }
    if ( _SortEigensInAbsIncreasingOrder( VALPROP, VECPROP, 3 ) != 1 ) {
      theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0.0;
      continue;
    }



    switch( color ) {
    default:
    case _WHITE_ :

      /* le point appartient a un tube blanc sur fond
	 noir si les deux plus grandes sont negatives
	 et beaucoup plus grandes que la plus petite
	 
	 On a fabs(valprop[0]) <= fabs(valprop[1]) <= fabs(valprop[2])
	 
	 On veut :
	 1.        valprop[2] < 0
	 2.        fabs(valprop[1]) * 2   <=   fabs(valprop[2]) si MODELBASED
	           fabs(valprop[2]) <=  FLTZERO si ACME
      */
      if ( VALPROP[2] >= 0 ) {
	theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0.0;
	continue;
      }
      break;

    case _BLACK_ :

      /* le point appartient a un tube noir sur fond
	 blanc si les deux plus grandes sont positives
	 et beaucoup plus grandes que la plus petite
	 
	 On a fabs(VALPROP[0]) <= fabs(VALPROP[1]) <= fabs(VALPROP[2])
	 
	 On veut :
	 1.        VALPROP[2] > 0
	 2.        fabs(VALPROP[1]) * 2   <=   fabs(VALPROP[2]) si MODELBASED
	           fabs(VALPROP[2]) <=  FLTZERO si ACME
      */
      if ( VALPROP[2] <= 0 ) {
	theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0.0;
	continue;
      }
      break;

    } 

    /* OK, ici on regarde le point 
       le vecteur normal a la membrane est celui
       associe a la plus grande valeur propre
       ie vecprop[2,5,8]
       les deux autres sont vecprop[1,4,7] et vecprop[0,3,6]
    */
    
    switch (mode) {
    default:
    case MODELBASED :

      //if ( fabs(VALPROP[2]) <= 2 * fabs(VALPROP[1]) ) {
      //  theRep[y][x] = theTheta[y][x] = thePhi[y][x] = 0.0;
      //  continue;
      //}


      nbValPts = nbPosPts = nbNegPts = 0;
      sumValPts = sumPosPts = sumNegPts = 0;


      for (a=-1 ; a<=1 ; a+=2)       // de part et d'autre de la membrane...
      {

        vx=(double)a*VECPROP[2];
        vy=(double)a*VECPROP[5];
        vz=(double)a*VECPROP[8];


        //qx=radius*vx+x;
        //qy=radius*vy+y;
        //qz=radius*vz+slice;
        qx=radii[0]*vx+x;
        qy=radii[1]*vy+y;
        qz=radii[2]*vz+slice;

        for (a1=-hsp;a1<=hsp;a1++)      // on se deplace selon v1...
        {

          wx=(double)a1*VECPROP[0];
          wy=(double)a1*VECPROP[3];
          wz=(double)a1*VECPROP[6];

          for (a2=-hsp;a2<=hsp;a2++)      // puis selon v2...
          {
            //wx+=(double)a2*vecprop[1];
            //wy+=(double)a2*vecprop[4];
            //wz+=(double)a2*vecprop[7];
            wwx=wx+(double)a2*VECPROP[1];
            wwy=wy+(double)a2*VECPROP[4];
            wwz=wz+(double)a2*VECPROP[7];

            rx=theCoeff[0]*wwx+qx;
            ry=theCoeff[1]*wwy+qy;
            rz=theCoeff[2]*wwz+qz;

            if (rx<=0.0 || ry<=0.0 || rz<=0.0) {
              continue;
            }
                                          // et on somme les reponses
            ix=(int)rx;
            iy=(int)ry;
            iz=(int)rz;

            if (ix>=dimx-1 || iy>=dimy-1 || iz>=dimz-1) {
              continue;
            }
        
            dx = rx-ix;
            dy = ry-iy;
            dz = rz-iz;

            dxdy = dx*dy;
            dxdz = dx*dz;
            dydz = dy*dz;
            dxdydz = dxdy*dz;

            v6 = dxdz-dxdydz;
            v5 = dxdy-dxdydz;
            v4 = dx-dxdy-v6;
        
            gx = 0;
            gx += dxdydz        * theX[iz+1][iy+1][ix+1];
            gx += (dydz-dxdydz) * theX[iz+1][iy+1][ix  ];
            gx += v6            * theX[iz+1][iy  ][ix+1];
            gx += (dz-dydz-v6)  * theX[iz+1][iy  ][ix  ];
            gx += v5            * theX[iz  ][iy+1][ix+1];
            gx += (dy-dydz-v5)  * theX[iz  ][iy+1][ix  ];
            gx += v4            * theX[iz  ][iy  ][ix+1];
            gx += (1-dy-dz+dydz-v4) * theX[iz  ][iy][ix];
            gx *= theCoeff[0];

            gy = 0;
            gy += dxdydz        * theY[iz+1][iy+1][ix+1];
            gy += (dydz-dxdydz) * theY[iz+1][iy+1][ix  ];
            gy += v6            * theY[iz+1][iy  ][ix+1];
            gy += (dz-dydz-v6)  * theY[iz+1][iy  ][ix  ];
            gy += v5            * theY[iz  ][iy+1][ix+1];
            gy += (dy-dydz-v5)  * theY[iz  ][iy+1][ix  ];
            gy += v4            * theY[iz  ][iy  ][ix+1];
            gy += (1-dy-dz+dydz-v4) * theY[iz  ][iy][ix];
            gy *= theCoeff[1];

            gz = 0;
            gz += dxdydz        * theZ[iz+1][iy+1][ix+1];
            gz += (dydz-dxdydz) * theZ[iz+1][iy+1][ix  ];
            gz += v6            * theZ[iz+1][iy  ][ix+1];
            gz += (dz-dydz-v6)  * theZ[iz+1][iy  ][ix  ];
            gz += v5            * theZ[iz  ][iy+1][ix+1];
            gz += (dy-dydz-v5)  * theZ[iz  ][iy+1][ix  ];
            gz += v4            * theZ[iz  ][iy  ][ix+1];
            gz += (1-dy-dz+dydz-v4) * theZ[iz  ][iy][ix];
            gz *= theCoeff[2];

            val = - (gx*vx+gy*vy+gz*vz);

            nbValPts++;
            sumValPts += val;
        
            if(val>0) {
              nbPosPts++;
              sumPosPts+=val;
            }
            if(val<0) {
              nbNegPts++;
              sumNegPts-=val;
            }
          }
        }
      }
      if (color==_BLACK_) {
        int nb;
        double sum;
        sumValPts = (-sumValPts);
        sum = sumPosPts;
        sumPosPts = sumNegPts;
        sumNegPts = sum;

        nb = nbPosPts;
        nbPosPts = nbNegPts;
        nbNegPts = nb;
      }


      if (sumValPts<=0.0) {
        continue;
      }
      
      R = sumValPts/(double)nbValPts;
      break;
      
    case ACME :
      /*  Calcul de la reponse de type ACME */


      /* Calcul de A, B et S */
      A = fabs(VALPROP[1]/VALPROP[2]);
      B = fabs(sqrt(fabs(VALPROP[0]*VALPROP[1]))/VALPROP[2]);
      S = sqrt(VALPROP[0]*VALPROP[0]+VALPROP[1]*VALPROP[1]
               +VALPROP[2]*VALPROP[2]);

      T0 = 1 - exp(-S*S/(2*gamma*gamma));       // Foreground vs. Background
      T1 = exp(-A*A/(2*alpha*alpha));           // Plane vs. Rod
      T2 = exp(-B*B/(2*beta*beta));             // Plane vs. Ball
      T3 = exp(-2*c*c/(VALPROP[2]*VALPROP[2])); // Smooth plane

      /* Reponse */
      R = T0*T1*T2*T3;

      break;
    }
    
    //vecteur normal au plan: v = vecprop associe a la grande VALPROP de hessien
    v[0] = VECPROP[2];
    v[1] = VECPROP[5];
    v[2] = VECPROP[8];

    theRep[y][x] = R;
    UnitVectorToSphericalAngles( v, &theta, &phi );
    theTheta[y][x] = theta;
    thePhi[y][x]   = phi;

  } // fin du  for(y) for(x)...

}





int MT_Filter3Dimages( vt_image *im, vt_3Dimages *par, float *theCoeffs )
{
  char *proc = "MT_Filter3Dimages";
  filterType theFilter = GABOR_YOUNG_2002;
  derivativeOrder theDerivatives[3] = {NODERIVATIVE,NODERIVATIVE,NODERIVATIVE};
  int theDim[3];
  int i;
  typeFilteringCoefficients filter[3];
  for (i=0; i<3; i++)
  {
    initFilteringCoefficients(&filter[i]);
    filter[i].type=theFilter;
    filter[i].coefficient=theCoeffs[i];
  }

  theDim[0] = im->dim.x;
  theDim[1] = im->dim.y;
  theDim[2] = im->dim.z;

  
  theDerivatives[0] = NODERIVATIVE;
  theDerivatives[1] = NODERIVATIVE;
  theDerivatives[2] = DERIVATIVE_0;

    

  if ( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (.,.,0)\n" );

  if ( _verbose_ >= 2 ) {
    VT_PrintImage( stderr, im, "input image" );
    VT_PrintImage( stderr, &(par->tmp0), "par->tmp0" );
  }

  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering( im->buf, im->type,
				par->tmp0.buf, par->tmp0.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (.,.,0)\n", proc );

    return( -1 );
  }



  theDerivatives[0] = DERIVATIVE_0;
  theDerivatives[1] = DERIVATIVE_1;
  theDerivatives[2] = NODERIVATIVE;

  if( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (0,1,.)\n" );

  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering( par->tmp0.buf, par->tmp0.type,
				par->imy.buf, par->imy.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,1,0)\n", proc );

    return( -1 );
  }



  theDerivatives[0] = DERIVATIVE_1;
  theDerivatives[1] = DERIVATIVE_0;
  theDerivatives[2] = NODERIVATIVE;

  if( _verbose_ >= 2 ) fprintf( stderr, " 3D filtering (1,0,.)\n" );

  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering( par->tmp0.buf, par->tmp0.type,
				par->imx.buf, par->imx.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (1,0,0)\n", proc );

    return( -1 );
  }



  theDerivatives[0] = NODERIVATIVE;
  theDerivatives[1] = NODERIVATIVE;
  theDerivatives[2] = DERIVATIVE_1;

  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering( im->buf, im->type,
				par->tmp1.buf, par->tmp1.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (.,.,1)\n", proc );

    return( -1 );
  }



  theDerivatives[0] = DERIVATIVE_0;
  theDerivatives[1] = DERIVATIVE_0;
  theDerivatives[2] = NODERIVATIVE;

  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering( par->tmp1.buf, par->tmp1.type,
				par->imz.buf, par->imz.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,0,1)\n", proc );

    return( -1 );
  }



  theDerivatives[0] = DERIVATIVE_0;
  theDerivatives[1] = DERIVATIVE_0;
  theDerivatives[2] = DERIVATIVE_2;

  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering( im->buf, im->type,
				par->imzz.buf, par->imzz.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error in computating derivatives (0,0,2)\n", proc );
  
    return( -1 );
  }
  
  return( 1 );
}


int MT_Filter2Dimauxs( vt_3Dimages *ims, vt_2Dimauxs *par,
		       float *theCoeffs, int slice )
{
  char *proc = "MT_Filter2Dimauxs";
  filterType theFilter = GABOR_YOUNG_2002;
  derivativeOrder theDerivatives[3] = {NODERIVATIVE,NODERIVATIVE,NODERIVATIVE};
  int theDim[3];
  
  float ***theTmp0 = (float ***)(ims->tmp0.array);
  float ***theTmp1 = (float ***)(ims->tmp1.array);

  if ( ims->tmp0.type != FLOAT || ims->tmp1.type != FLOAT ) {
    return( -1 );
  }

  int i;
  typeFilteringCoefficients filter[3];
  for (i=0; i<3; i++)
  {
    initFilteringCoefficients(&filter[i]);
    filter[i].type=theFilter;
    filter[i].coefficient=theCoeffs[i];
  }

  theDim[0] = ims->imzz.dim.x;
  theDim[1] = ims->imzz.dim.y;
  theDim[2] = 1;


  theDerivatives[0] = DERIVATIVE_2;
  theDerivatives[1] = DERIVATIVE_0;

  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering((void*)(&(theTmp0[slice][0][0])), ims->tmp0.type,
				par->imxx.buf, par->imxx.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr,
               "%s: error in computating derivatives (2,0,0) in slice %d\n",
	       proc, slice );
    return( -1 );
  }
  
  
  theDerivatives[0] = DERIVATIVE_0;
  theDerivatives[1] = DERIVATIVE_2;
  
  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering((void*)(&(theTmp0[slice][0][0])), ims->tmp0.type,
				par->imyy.buf, par->imyy.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr,
               "%s: error in computating derivatives (0,2,0) in slice %d\n",
	       proc, slice );
    return( -1 );
  }
  
  
  theDerivatives[0] = DERIVATIVE_1;
  theDerivatives[1] = DERIVATIVE_1;
  
  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering((void*)(&(theTmp0[slice][0][0])), ims->tmp0.type,
				par->imxy.buf, par->imxy.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr,
               "%s: error in computating derivatives (1,1,0) in slice %d\n",
	       proc, slice );
    return( -1 );
  }
  
  
  theDerivatives[0] = DERIVATIVE_1;
  theDerivatives[1] = DERIVATIVE_0;
  
  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering((void*)(&(theTmp1[slice][0][0])), ims->tmp1.type,
				par->imxz.buf, par->imxz.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr,
               "%s: error in computating derivatives (1,0,1) in slice %d\n",
	       proc, slice );
    return( -1 );
  }
  
  
  theDerivatives[0] = DERIVATIVE_0;
  theDerivatives[1] = DERIVATIVE_1;
  
  for (i=0; i<3; i++)
  {
    filter[i].derivative=theDerivatives[i];
  }

  if ( separableLinearFiltering((void*)(&(theTmp1[slice][0][0])), ims->tmp1.type,
				par->imyz.buf, par->imyz.type,
				theDim, _borderLength_,
				filter ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr,
               "%s: error in computating derivatives (0,1,1) in slice %d\n",
	       proc, slice );
    return( -1 );
  }
  
  
  
  return( 1 );
}











///////////////// TENSOR VOTING ///////////////////////////////////////////////


int MT_SampleBin(vt_image *imageBin, double sample)
{
  char *proc = "MT_SampleBin";

  int nind = 0;
  int *ind;
  int i, j;
  int n0;
  int tmp;

  unsigned char* bufU8=NULL;
  float* bufFLT=NULL;

  if (sample < 0.0 || sample > 1.0)
  {
    fprintf(stderr, "%s : bad sample argument, expected to be between 0 and 1 (sample=%f)", proc, sample);
    return(-1);
  }

  switch (imageBin->type)
  {
  case UCHAR :
    bufU8=(unsigned char*)imageBin->buf;
    for (i=0; i<(imageBin->dim.x*imageBin->dim.y*imageBin->dim.z); i++)
      if(bufU8[i] != '\0')
        nind ++;
    break;
  case FLOAT :
    bufFLT=(float*)imageBin->buf;
    for (i=0; i<(imageBin->dim.x*imageBin->dim.y*imageBin->dim.z); i++)
      if(bufFLT[i] > 0.0)
        nind ++;
    break;
  default:
    fprintf(stderr, "%s : image extension not implemented yet\n", proc);
    return(-1);
  }

  ind=malloc(nind*sizeof(int));
  j=0;

  switch (imageBin->type)
  {
  case UCHAR :
    for (i=0; i<(imageBin->dim.x*imageBin->dim.y*imageBin->dim.z); i++)
      if(bufU8[i] != '\0')
        ind[j++] = i;
    break;
  case FLOAT :
    for (i=0; i<(imageBin->dim.x*imageBin->dim.y*imageBin->dim.z); i++)
      if(bufFLT[i] > 0.0)
        ind[j++] = i;
    break;
  default:
    free(ind);
    ind=(int*)NULL;
    fprintf(stderr, "%s : image extension not implemented yet\n", proc);
    return(-1);
  }

  n0 = (int) ( ((double) nind)*(1.0-sample) );

  fprintf(stdout, "nind=%d\tn0=%d\n", nind, n0);

  srand(time(NULL));

  for (i=0; i<n0; i++)
  {
    j=rand_(nind-i);
    tmp=ind[j];
    ind[j]=ind[nind-i-1];
    ind[nind-i-1]=tmp;
    switch (imageBin->type)
    {
    case UCHAR :
      bufU8[tmp]='\0';
      break;
    case FLOAT :
      bufFLT[tmp]=0.0;
      break;
    default:
      free(ind);
      ind=(int*)NULL;
      fprintf(stderr, "%s : image extension not implemented yet\n", proc);
      return(-1);
    }
  }

  if(0 && _verbose_)
  {
    fprintf(stdout, "indices supprimes : [   ");
    for (i=0; i<n0; i++)
      fprintf(stdout, "%d   ", ind[nind-n0+i]);
    fprintf(stdout, "]\n");
  }

  free(ind);
  ind=(int*)NULL;

  return (1);
}

int MT_Compute3DTensorVoting( vt_3Dtensor *theTensor,
			    vt_image **imBin,
			    double scale,
			    double zfact,
			    int  Niter,
                            int  NanglesIter,
                            int Nsticks,
                            enumTVmode mode,
                            hessianMode initHessian,
                            char *parName,
                            int writeImages)
{

  char *proc = "MT_Compute3DTensorVoting";
  int hsf;

  int dimFields[3];
  

  double **angles;
  int Nangles;

  double theCoeff[3];
  theCoeff[0] = theCoeff[1] = scale;
  theCoeff[2] = zfact*scale;
  hsf = ceil(sqrt(-pow(scale,2)*log(0.01)));  
  dimFields[0] = dimFields[1] = 2*hsf+1;
  dimFields[2] = (int)(2*zfact*hsf+1);

  int x,y,z;
  int dimx = imBin[0]->dim.x;
  int dimy = imBin[0]->dim.y;
  int dimz = imBin[0]->dim.z;

  float ***theXX;
  float ***theYY;
  float ***theZZ;
  float ***theXY;
  float ***theYZ;
  float ***theXZ;

  float ***theVP1;
  float ***theVP2;
  float ***theVP3;

  float ***theTheta1;
  float ***thePhi1;
  float ***theTheta2;
  float ***thePhi2;
  float ***theTheta3;
  float ***thePhi3;

  unsigned char ***theZEROS;

  float ***imTht=NULL;
  float ***imPhi=NULL;

  float*** bin=NULL;
  unsigned char ***bin0;
  

  switch (imBin[0]->type)
  {
  case UCHAR :
    bin0 = (unsigned char***)imBin[0]->array;
    bin=malloc(dimz*sizeof(float**));
    for (z=0;z<dimz;z++)
    {
      bin[z]=malloc(dimy*sizeof(float*));
      for (y=0;y<dimy;y++)
      {
        bin[z][y]=malloc(dimx*sizeof(float));
        for (x=0;x<dimx;x++)
          bin[z][y][x] = bin0[z][y][x] > 0 ? (float) 1 : (float) 0;
      }
    }
    break;
  case FLOAT :
    bin = (float***)imBin[0]->array;
    break;
  default:
    break;
  }


  if (initHessian!=NONE)
  {
    imTht = (float***)imBin[1]->array;
    imPhi = (float***)imBin[2]->array;
  }

  int n,m,i;
  //ImageType t = DOUBLE; 

  vt_3Dtensor *sfields;
  vt_3Dtensor bfield;
  vt_3Dtensor *pfields;
  
  double lambda1, lambda2, lambda3;
  
  enumVote typeVote;
  char name[256];

  //float hessien[9];
  double theta, phi;
  double v1[3],v2[3],v3[3];

  float val;

  /* Calcul des angles répartis de façon homogène sur la boule unite */
  if (_verbose_)
    fprintf(stdout, "%s: computating angles...\n", proc);

  Nangles = MT_Compute3DAngles(&angles, NanglesIter);
  if (Nangles <= 0) {        
    if (_verbose_)
      fprintf(stderr, "%s: error in computating angles\n", proc);
    return( -1 );
  }
  

  
  /* Calcul des champs de tenseurs Stick/Plate/Ball */
  if (initHessian != LINE && _verbose_)
    fprintf(stdout, "%s: computating stick fields...\n", proc);
  
  if (initHessian != LINE && MT_Compute3DStickFields( &sfields, dimFields, 
          angles, theCoeff, Nangles, mode) != 1)
  {
    for (n=0;n<Nangles;n++){
      free(angles[n]);
      angles[n]=NULL;
    }
    free(angles);
    angles = NULL;
    switch (imBin[0]->type)
    {
    case UCHAR :
      for(z=0;z<dimz;z++)
      {
        for(y=0;y<dimy;y++)
        {
          free(bin[z][y]);
          bin[z][y]=NULL;
        }
        free(bin[z]);
        bin[z]=NULL;
      }
      free(bin);
      bin=NULL;
      break;
    default:
      break;
    }
    return (-1);
  }

  if (initHessian != PLANE && _verbose_)
    fprintf(stdout, "%s: computating plate fields...\n", proc);

  if (initHessian != PLANE && MT_Compute3DPlateFields( &pfields, dimFields, theCoeff,
      angles, Nangles, Nsticks, mode) != 1)
  {
    for (n=0;n<Nangles;n++){
      free(angles[n]);
      angles[n]=NULL;
      if (initHessian != LINE)
		VT_Free3Dtensor(&sfields[n]);
    }
    if (initHessian != LINE)
	  free(sfields);
    sfields=NULL;
    free(angles);
    angles = (double**)NULL;    
    switch (imBin[0]->type)
    {
    case UCHAR :
      for(z=0;z<dimz;z++)
      {
        for(y=0;y<dimy;y++)
        {
          free(bin[z][y]);
          bin[z][y]=NULL;
        }
        free(bin[z]);
        bin[z]=NULL;
      }
      free(bin);
      bin=NULL;
      break;
    default:
      break;
    }
    return (-1);
  }  

  
  if (initHessian == NONE && _verbose_)
    fprintf(stdout, "%s: computating ball fields...\n", proc);

  if (initHessian == NONE && MT_Compute3DBallFieldFromStick( &bfield, sfields, dimFields,
          Nangles ) != 1)
  {
    for (n=0;n<Nangles;n++){
      free(angles[n]);
      angles[n]=NULL;
      VT_Free3Dtensor(sfields+n);
      VT_Free3Dtensor(pfields+n);
    }
    switch (imBin[0]->type)
    {
    case UCHAR :
      for(z=0;z<dimz;z++)
      {
        for(y=0;y<dimy;y++)
        {
          free(bin[z][y]);
          bin[z][y]=NULL;
        }
        free(bin[z]);
        bin[z]=NULL;
      }
      free(bin);
      bin=NULL;
      break;
    default:
      break;
    }
    free(sfields);
    sfields=NULL;
    free(pfields);
    pfields=NULL;
    free(angles);
    angles = NULL;
    return (-1);
  }  
  
  /*  Tensor image allocation */
  if (_verbose_)
    fprintf(stdout, "%s: allocating tensor image...\n", proc);
  
  sprintf( name, "%s", parName );
  if (VT_Alloc3DtensorFromImage( theTensor, (imBin[0]), name ) !=1)
  {
    if (_verbose_)
      fprintf(stderr, "%s: allocation of tensor image failed\n", proc);

    for (n=0;n<Nangles;n++){
      free(angles[n]);
      angles[n]=NULL;
      if (initHessian != LINE)
		VT_Free3Dtensor(sfields+n);
      if (initHessian != PLANE)
        VT_Free3Dtensor(pfields+n);
    }
    if (initHessian != LINE)
	  free(sfields);
    sfields=NULL;
    if (initHessian == NONE)
    {
      VT_Free3Dtensor(&bfield);
    }
    if (initHessian != PLANE )
    {
	  free(pfields);
      pfields=NULL;
	}
    free(angles);
    angles = NULL;
    switch (imBin[0]->type)
    {
    case UCHAR :
      for(z=0;z<dimz;z++)
      {
        for(y=0;y<dimy;y++)
        {
          free(bin[z][y]);
          bin[z][y]=NULL;
        }
        free(bin[z]);
        bin[z]=NULL;
      }
      free(bin);
      bin=NULL;
      break;
    default:
      break;
    }

    return (-1);

  }

  theXX = (float ***)theTensor->imxx.array;
  theYY = (float ***)theTensor->imyy.array;
  theZZ = (float ***)theTensor->imzz.array;
  theXY = (float ***)theTensor->imxy.array;
  theYZ = (float ***)theTensor->imyz.array;
  theXZ = (float ***)theTensor->imxz.array;

  theVP1 = (float ***)theTensor->imvp1.array;
  theVP2 = (float ***)theTensor->imvp2.array;
  theVP3 = (float ***)theTensor->imvp3.array;

  theTheta1 = (float ***)theTensor->imtheta1.array;
  thePhi1 = (float ***)theTensor->imphi1.array;
  theTheta2 = (float ***)theTensor->imtheta2.array;
  thePhi2 = (float ***)theTensor->imphi2.array;
  theTheta3 = (float ***)theTensor->imtheta3.array;
  thePhi3 = (float ***)theTensor->imphi3.array;

  theZEROS = (unsigned char ***)theTensor->iszero.array;


  /* Initialisation de l'image tenseur a partir d'imBin */
  if (_verbose_)
    fprintf(stdout,"%s: initialisation du tenseur a partir d'imBin...\n", proc);

  n=0;
  m=0;

  for (z=0;z<dimz;z++)
  {
    if (_verbose_ && (dimz<=30 || z%(dimz/30)==0))
      fprintf(stdout,"\tFrame #%d/%d (%d%%)... \n",
          z+1, dimz, (int)(100*(float)z/(float)dimz));
    for (y=0;y<dimy;y++)
    for (x=0;x<dimx;x++)
    {
      val = bin[z][y][x];
      if ( val < FLTZERO)
      {
        /*  le voxel(x,y,z) ne vote pas ni ne reçoit de votes  */
        theXX[z][y][x] = 0;
        theYY[z][y][x] = 0;
        theZZ[z][y][x] = 0;
        theXY[z][y][x] = 0;
        theXZ[z][y][x] = 0;
        theYZ[z][y][x] = 0;
        theZEROS[z][y][x] = (unsigned char)1;
        n++;
        continue;
      }
      if (initHessian==NONE)
      {
      /* Tenseur Ball */
      theXX[z][y][x] = val;
      theYY[z][y][x] = val;
      theZZ[z][y][x] = val;
      theXY[z][y][x] = 0;
      theXZ[z][y][x] = 0;
      theYZ[z][y][x] = 0;

      theVP1[z][y][x] = val;
      theVP2[z][y][x] = val;
      theVP3[z][y][x] = val;
      }
      if (initHessian==PLANE)
      {
        /* Tenseur Stick */
        theTheta3[z][y][x] = imTht[z][y][x];
        thePhi3[z][y][x]   = imPhi[z][y][x];

        theVP1[z][y][x] = 0;
        theVP2[z][y][x] = 0;
        theVP3[z][y][x] = val;
      }
      if (initHessian==LINE)
      {
        /* Tenseur Plate */
        SphericalAnglesToUnitsVectors( imTht[z][y][x], imPhi[z][y][x], v1, v2, v3 );
		UnitVectorToSphericalAngles(v2,&theta,&phi);
        theTheta2[z][y][x] = theta;
        thePhi2[z][y][x]   = phi;
        UnitVectorToSphericalAngles(v3,&theta,&phi);
        theTheta3[z][y][x] = theta;
        thePhi3[z][y][x]   = phi;

        theVP1[z][y][x] = 0;
        theVP2[z][y][x] = val;
        theVP3[z][y][x] = val;
      }
      theZEROS[z][y][x] = (unsigned char)0;
      m++;
    }
  }

  if (initHessian != NONE)
    MT_InitTensorFromAngles(theTensor);


  if (_verbose_)
  {
    fprintf(stdout, "%s : nombre de zeros = %d \tnombre de non-zeros = %d \n",
        proc, n,m);
  }

  /* Niter etapes de votes SPARSE + 1 etape de votes DENSE */

  typeVote = SPARSE;

  if (initHessian != NONE)
    Niter=0;

  for (i=0;i<=Niter;i++)
  {

    if (writeImages)
      {
      if (_verbose_)
      {
        sprintf( name, "%s.eparse.%d", parName,i );
        fprintf(stdout,
        "Ecriture du champ de tenseurs eparse #%d dans %s.inr...\n", i, name);
      }
      VT_Write3DtensorWithName( theTensor, name );
    }
    if (_verbose_)
      fprintf(stdout, "%s: etape de vote : iteration #%d/%d...\n",
        proc, i+1, Niter+1);

    /* Vote DENSE si derniere iteration */
    if (i == Niter)
    {
        typeVote = DENSE;
    }


    if (_verbose_) {
      if ( typeVote == SPARSE) fprintf(stdout, "\tSPARSE voting...\n");
      else fprintf(stdout, "\tDENSE voting...\n");
	}
	
    MT_ReinitTensorImage(theTensor);

    /* Cumul des votes (seuls imxx,yy,zz,xy,xz,yz et iszero sont affectes) */
    for (z=0;z<dimz;z++)
    {
      if (_verbose_)
        if (dimz<=30 || z%(dimz/30)==0)
          fprintf(stdout, "\t\tTreating frame #%d/%d (%d%%)...\n",
              z+1, dimz, (int)(100*(float)z/(float)dimz));

      for (y=0;y<dimy;y++)
      for (x=0;x<dimx;x++)
      {
        if ( theZEROS[z][y][x] != 0)
        {
          /*  le voxel(x,y,z) ne vote pas  */
          continue;
        }

        lambda1 = theVP1[z][y][x];
        lambda2 = theVP2[z][y][x];
        lambda3 = theVP3[z][y][x];

        if (lambda3-lambda2 > FLTZERO)  // STICKS
        {
          theta = theTheta3[z][y][x];
          phi = thePhi3[z][y][x];
          n = nearestAngle(theta, phi, angles, Nangles);
          MT_AddField(theTensor, sfields[n], lambda3-lambda2,x,y,z, typeVote);
        }
        //if (initHessian == 0)
        //{
          if (lambda2-lambda1 > FLTZERO)  // PLATES
          {
            theta = theTheta1[z][y][x];
            phi = thePhi1[z][y][x];
            n = nearestAngle(theta, phi, angles, Nangles);
            MT_AddField(theTensor, pfields[n], lambda2-lambda1,x,y,z, typeVote);
          }
          if (lambda1 > FLTZERO)          // BALLS
            MT_AddField(theTensor, bfield, lambda1,x,y,z, typeVote);
        //}
      }
    }

    /* Conversion de la hessienne en vp, angles, etc... */
    if (_verbose_)
      fprintf(stdout, "\tConversion des tenseurs en valeurs/vecteurs propres, \
          etc...\n");

    for (z=0;z<dimz;z++)
    for (y=0;y<dimy;y++)
    for (x=0;x<dimx;x++)
    {
      if ( theZEROS[z][y][x] == 1)
      {
        /*  le voxel(x,y,z) contient un tenseur nul */
        continue;
      }

      if (setTensor(theTensor, x,y,z, typeVote) !=1)
      {
        if (_verbose_)
          fprintf(stderr, "%s: probleme lors du calcul des parametres des \
              tenseurs", proc);
        for (n=0;n<Nangles;n++){
          free(angles[n]);
          angles[n]=NULL;
          if (initHessian != LINE)
			VT_Free3Dtensor(sfields+n);
          if (initHessian != PLANE)
            VT_Free3Dtensor(pfields+n);
        }
		if (initHessian != LINE)
          free(sfields);
        sfields=NULL;
        if (initHessian != PLANE)
        {
          free(pfields);
          pfields=NULL;
        }
        if (initHessian == NONE )
		{
          VT_Free3Dtensor(&bfield);
		}

        free(angles);
        angles = NULL;
        VT_Free3Dtensor(theTensor);
        switch (imBin[0]->type)
        {
        case UCHAR :
          for(z=0;z<dimz;z++)
          {
            for(y=0;y<dimy;y++)
            {
              free(bin[z][y]);
              bin[z][y]=NULL;
            }
            free(bin[z]);
            bin[z]=NULL;
          }
          free(bin);
          bin=NULL;
          break;
        default:
          break;
        }
        return(-1);
      }

    }

  } // Niter


  if (_verbose_)
    fprintf(stdout, "%s: Liberation de la memoire allouee...\n", proc);

  for (n=0;n<Nangles;n++){
    free(angles[n]);
    angles[n]=NULL;
    if (initHessian != LINE)
	  VT_Free3Dtensor(sfields+n);
    if (initHessian != PLANE)
      VT_Free3Dtensor(pfields+n);
  }
  if (initHessian != LINE)
    free(sfields);
  sfields=NULL;
  if (initHessian != PLANE)
  {
    free(pfields);
    pfields=NULL;
  }
  if (initHessian == NONE )
  {
    VT_Free3Dtensor(&bfield);
  }

  free(angles);
  angles = NULL;
  switch (imBin[0]->type)
  {
  case UCHAR :
    for(z=0;z<dimz;z++)
    {
      for(y=0;y<dimy;y++)
      {
        free(bin[z][y]);
        bin[z][y]=NULL;
      }
      free(bin[z]);
      bin[z]=NULL;
    }
    free(bin);
    bin=NULL;
    break;
  default:
    break;
  }

  
  return(1);
}

////////////////////////////////////////////////////////////////////////////////
/// \brief MT_Cumul3DTensorLine
/// \param theTensor
/// \param imBin
/// \param r
/// \param name
/// \return
///
int MT_Cumul3DTensorLine(vt_3Dtensor *theTensor, vt_image **imBin, double r, char *name)
{
    char *proc = "MT_Cumul3DTensorLine";

    int x,y,z;
    int dimx = imBin[0]->dim.x;
    int dimy = imBin[0]->dim.y;
    int dimz = imBin[0]->dim.z;
    int dim[3];
    dim[0]=dimx; dim[1]=dimy; dim[2]=dimz;

    float stick[9];
    double n[3];
    int pt[3];

    double tht, phi;
    float val;

    float ***theXX;
    float ***theYY;
    float ***theZZ;
    float ***theXY;
    float ***theYZ;
    float ***theXZ;

    float ***theVP1;
    float ***theVP2;
    float ***theVP3;

    /*
    float ***theTheta1;
    float ***thePhi1;
    float ***theTheta2;
    float ***thePhi2;
    float ***theTheta3;
    float ***thePhi3;
    */

    unsigned char ***theZEROS;

    float ***imTht;
    float ***imPhi;

    float*** bin=NULL;
    unsigned char ***bin0;


    switch (imBin[0]->type)
    {
    case UCHAR :
      bin0 = (unsigned char***)imBin[0]->array;
      bin=malloc(dimz*sizeof(float**));
      for (z=0;z<dimz;z++)
      {
        bin[z]=malloc(dimy*sizeof(float*));
        for (y=0;y<dimy;y++)
        {
          bin[z][y]=malloc(dimx*sizeof(float));
          for (x=0;x<dimx;x++)
            bin[z][y][x] = bin0[z][y][x] > 0 ? (float) 1 : (float) 0;
        }
      }
      break;
    case FLOAT :
      bin = (float***)imBin[0]->array;
      break;
    default:
      break;
    }
    imTht = (float***)imBin[1]->array;
    imPhi = (float***)imBin[2]->array;

    /*  Tensor image allocation */
    if (_verbose_)
      fprintf(stdout, "%s: allocating tensor image...\n", proc);

    if (VT_Alloc3DtensorFromImage( theTensor, (imBin[0]), name ) !=1)
    {
      if (_verbose_)
        fprintf(stderr, "%s: allocation of tensor image failed\n", proc);
      switch (imBin[0]->type)
      {
      case UCHAR :
        for(z=0;z<dimz;z++)
        {
          for(y=0;y<dimy;y++)
          {
            free(bin[z][y]);
            bin[z][y]=NULL;
          }
          free(bin[z]);
          bin[z]=NULL;
        }
        free(bin);
        bin=NULL;
        break;
      default:
        break;
      }
      return (-1);
    }

    theXX = (float ***)theTensor->imxx.array;
    theYY = (float ***)theTensor->imyy.array;
    theZZ = (float ***)theTensor->imzz.array;
    theXY = (float ***)theTensor->imxy.array;
    theYZ = (float ***)theTensor->imyz.array;
    theXZ = (float ***)theTensor->imxz.array;

    theVP1 = (float ***)theTensor->imvp1.array;
    theVP2 = (float ***)theTensor->imvp2.array;
    theVP3 = (float ***)theTensor->imvp3.array;

    /*
    theTheta1 = (float ***)theTensor->imtheta1.array;
    thePhi1 = (float ***)theTensor->imphi1.array;
    theTheta2 = (float ***)theTensor->imtheta2.array;
    thePhi2 = (float ***)theTensor->imphi2.array;
    theTheta3 = (float ***)theTensor->imtheta3.array;
    thePhi3 = (float ***)theTensor->imphi3.array;
    */

    theZEROS = (unsigned char ***)theTensor->iszero.array;

    /* INIT */
    for (z=0;z<dimz;z++)
    for (y=0;y<dimy;y++)
    for (x=0;x<dimx;x++)
    {
        /*  le voxel(x,y,z) ne vote pas ni ne reçoit de votes  */
        theXX[z][y][x] = 0;
        theYY[z][y][x] = 0;
        theZZ[z][y][x] = 0;
        theXY[z][y][x] = 0;
        theXZ[z][y][x] = 0;
        theYZ[z][y][x] = 0;
        theVP1[z][y][x] = 0;
        theVP2[z][y][x] = 0;
        theVP3[z][y][x] = 0;
        theZEROS[z][y][x] = (unsigned char)1;
    }

    /* CUMUL LIGNES */
    for (z=0;z<dimz;z++)
    for (y=0;y<dimy;y++)
    for (x=0;x<dimx;x++)
    {
        val = bin[z][y][x];
        if(val<=0) continue;
        pt[0]=x; pt[1]=y; pt[2]=z;
        tht=imTht[z][y][x];
        phi=imPhi[z][y][x];
        SphericalAnglesToUnitVector( tht, phi, n);
        stick[0]=n[0]*n[0];
        stick[1]=n[0]*n[1];
        stick[2]=n[0]*n[2];
        stick[3]=n[1]*n[0];
        stick[4]=n[1]*n[1];
        stick[5]=n[1]*n[2];
        stick[6]=n[2]*n[0];
        stick[7]=n[2]*n[1];
        stick[8]=n[2]*n[2];

        if (    MT_AddLineToTensorImage3D(theTensor->imxx.buf, dim, theTensor->imxx.type, pt, tht, phi, r, val*stick[0])!=1 ||
                MT_AddLineToTensorImage3D(theTensor->imxy.buf, dim, theTensor->imxy.type, pt, tht, phi, r, val*stick[1])!=1 ||
                MT_AddLineToTensorImage3D(theTensor->imxz.buf, dim, theTensor->imxz.type, pt, tht, phi, r, val*stick[2])!=1 ||
                MT_AddLineToTensorImage3D(theTensor->imyy.buf, dim, theTensor->imyy.type, pt, tht, phi, r, val*stick[4])!=1 ||
                MT_AddLineToTensorImage3D(theTensor->imyz.buf, dim, theTensor->imyz.type, pt, tht, phi, r, val*stick[5])!=1 ||
                MT_AddLineToTensorImage3D(theTensor->imzz.buf, dim, theTensor->imzz.type, pt, tht, phi, r, val*stick[8])!=1 ||
                MT_ResetLineOfTensorImage3D(theTensor->iszero.buf, dim, theTensor->iszero.type, pt, tht, phi, r)!=1 )
        {
            if (_verbose_){
              fprintf(stderr, "%s: probleme lors de l'addition des tenseurs\n", proc);
              fprintf(stderr, "stick[%d][%d][%d]=\n%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f\n",
                    z, y, x,
                    stick[0], stick[1], stick[2], stick[3], stick[4], stick[5], stick[6], stick[7], stick[8]);
            }

            VT_Free3Dtensor(theTensor);
            switch (imBin[0]->type)
            {
            case UCHAR :
              for(z=0;z<dimz;z++)
              {
                for(y=0;y<dimy;y++)
                {
                  free(bin[z][y]);
                  bin[z][y]=NULL;
                }
                free(bin[z]);
                bin[z]=NULL;
              }
              free(bin);
              bin=NULL;
              break;
            default:
              break;
            }
            return(-1);

        }
    }

    /* Calcul VP, angles */
    /* Conversion de la hessienne en vp, angles, etc... */
    if (_verbose_)
      fprintf(stdout, "\tConversion des tenseurs en valeurs/vecteurs propres, \
          etc...\n");

    for (z=0;z<dimz;z++)
    for (y=0;y<dimy;y++)
    for (x=0;x<dimx;x++)
    {
      if ( theZEROS[z][y][x] == (unsigned char)1)
      {
        /*  le voxel(x,y,z) contient un tenseur nul */
        continue;
      }
      if (setTensor(theTensor, x,y,z, DENSE) !=1)
      {
        if (_verbose_){
          fprintf(stderr, "%s: probleme lors du calcul des parametres des tenseurs\n", proc);
          fprintf(stderr, "setTensor #%d\n",x+y*dimx+z*dimy*dimx);
          fprintf(stderr, "Hessian[%d][%d][%d]=\n%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f\n",
                z, y, x,
                theXX[z][y][x],theXY[z][y][x],theXZ[z][y][x],
                theXY[z][y][x],theYY[z][y][x],theYZ[z][y][x],
                theXZ[z][y][x],theYZ[z][y][x],theZZ[z][y][x]);
        }
        VT_Free3Dtensor(theTensor);
        switch (imBin[0]->type)
        {
        case UCHAR :
          for(z=0;z<dimz;z++)
          {
            for(y=0;y<dimy;y++)
            {
              free(bin[z][y]);
              bin[z][y]=NULL;
            }
            free(bin[z]);
            bin[z]=NULL;
          }
          free(bin);
          bin=NULL;
          break;
        default:
          break;
        }
        return(-1);
      }
    }

    /* Write tensor image */
    //sprintf(name, "%s", parName);
    //VT_Write3DtensorWithName( theTensor, name );

    /* Free memory */
    //VT_Free3Dtensor(theTensor);
    switch (imBin[0]->type)
    {
    case UCHAR :
      for(z=0;z<dimz;z++)
      {
        for(y=0;y<dimy;y++)
        {
          free(bin[z][y]);
          bin[z][y]=NULL;
        }
        free(bin[z]);
        bin[z]=NULL;
      }
      free(bin);
      bin=NULL;
      break;
    default:
      break;
    }

    return(1);
}


int MT_AddLineToTensorImage3D(void *buf, 	// buffer array
                 int *dim, 		// buffer dimensions
                 bufferType t, 	// buffer type
                 int *pt, 		// origin point
                 double tht, 	// angle 1
                 double phi, 	// angle 2
                 int r,         // rayon de "vote"
                 double val     // value to be added
                )
{
  char *proc = "MT_AddLineToTensorImage3D";

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

    v=(float) val;

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

int MT_ResetLineOfTensorImage3D(void *buf, int *dim, bufferType t, int *pt, double tht, double phi, int r)
{
    char *proc = "MT_ResetLineOfTensorImage3D";

    int Pt1[3],Pt2[3];
    int i, j;
    double n[3];
    double lambdas[6];
    double tmp;
    double maxdim;
    unsigned char *theBuf;

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
    case UCHAR:
    case SCHAR:
      theBuf = (unsigned char *)buf;

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


      /* x dominant
       */
      if (ax >= MAX(ay, az)) {
        yd = ay - (ax >> 1);
        zd = az - (ax >> 1);
        for (;;) {
          theBuf[z*ixy + y*ix + x] = (unsigned char) 0;
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
          theBuf[z*ixy + y*ix + x] = (unsigned char) 0;
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
          theBuf[z*ixy + y*ix + x] = (unsigned char) 0;
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

int MT_Compute3DTensorVotingTest(vt_3Dtensor *theTensor, vt_image **imBin, int NanglesIter, double r,
		double alpha, hessianMode initHessian, char *parName, int writeImages)
{

  char *proc = "MT_Compute3DTensorVotingTest";
  int hsf;

  int dimFields[3];
  

  double **angles;
  int Nangles;

  hsf = r;  
  dimFields[0] = dimFields[1] = dimFields[2] = 2*hsf+1;
//  dimFields[2] = (int)(2*zfact*hsf+1);

  int x,y,z;
  int dimx = imBin[0]->dim.x;
  int dimy = imBin[0]->dim.y;
  int dimz = imBin[0]->dim.z;

  float ***theXX;
  float ***theYY;
  float ***theZZ;
  float ***theXY;
  float ***theYZ;
  float ***theXZ;

  float ***theVP1;
  float ***theVP2;
  float ***theVP3;

  float ***theTheta1;
  float ***thePhi1;
  //float ***theTheta2;
  //float ***thePhi2;
  //float ***theTheta3;
  //float ***thePhi3;

  unsigned char ***theZEROS;

  float ***imTht=NULL;
  float ***imPhi=NULL;

  float*** bin=NULL;
  unsigned char ***bin0;
  

  switch (imBin[0]->type)
  {
  case UCHAR :
    bin0 = (unsigned char***)imBin[0]->array;
    bin=malloc(dimz*sizeof(float**));
    for (z=0;z<dimz;z++)
    {
      bin[z]=malloc(dimy*sizeof(float*));
      for (y=0;y<dimy;y++)
      {
        bin[z][y]=malloc(dimx*sizeof(float));
        for (x=0;x<dimx;x++)
          bin[z][y][x] = bin0[z][y][x] > 0 ? (float) 1 : (float) 0;
      }
    }
    break;
  case FLOAT :
    bin = (float***)imBin[0]->array;
    break;
  default:
    break;
  }


  if (initHessian!=NONE)
  {
    imTht = (float***)imBin[1]->array;
    imPhi = (float***)imBin[2]->array;
  }

  int n,m;
  //ImageType t = DOUBLE; 

  vt_3Dtensor *tfields;
  
  double lambda1, lambda2;//, lambda3;
  
  char name[256];

  //float hessien[9];
  double theta, phi;
  //double v1[3],v2[3],v3[3];

  float val;

  /* Calcul des angles répartis de façon homogène sur la boule unite */
  if (_verbose_)
    fprintf(stdout, "%s: computating angles...\n", proc);

  Nangles = MT_Compute3DAngles(&angles, NanglesIter);
  if (Nangles <= 0) {        
    if (_verbose_)
      fprintf(stderr, "%s: error in computating angles\n", proc);
    return( -1 );
  }
  

  
  /* Calcul des champs de tenseurs Test */
  if (initHessian != NONE && _verbose_)
    fprintf(stdout, "%s: computating test fields...\n", proc);
  
  if (initHessian != NONE && MT_Compute3DTestFields( &tfields, dimFields, 
          angles, Nangles, r, alpha) != 1)
  {
    for (n=0;n<Nangles;n++){
      free(angles[n]);
      angles[n]=NULL;
    }
    free(angles);
    angles = NULL;
    switch (imBin[0]->type)
    {
    case UCHAR :
      for(z=0;z<dimz;z++)
      {
        for(y=0;y<dimy;y++)
        {
          free(bin[z][y]);
          bin[z][y]=NULL;
        }
        free(bin[z]);
        bin[z]=NULL;
      }
      free(bin);
      bin=NULL;
      break;
    default:
      break;
    }
    return (-1);
  }
  
  /*  Tensor image allocation */
  if (_verbose_)
    fprintf(stdout, "%s: allocating tensor image...\n", proc);
  
  sprintf( name, "%s", parName );
  if (VT_Alloc3DtensorFromImage( theTensor, (imBin[0]), name ) !=1)
  {
    if (_verbose_)
      fprintf(stderr, "%s: allocation of tensor image failed\n", proc);

    for (n=0;n<Nangles;n++){
      free(angles[n]);
      angles[n]=NULL;
      if (initHessian != NONE)
		VT_Free3Dtensor(tfields+n);
    }
    if (initHessian != NONE)
	  free(tfields);
    tfields=NULL;
    free(angles);
    angles = NULL;
    switch (imBin[0]->type)
    {
    case UCHAR :
      for(z=0;z<dimz;z++)
      {
        for(y=0;y<dimy;y++)
        {
          free(bin[z][y]);
          bin[z][y]=NULL;
        }
        free(bin[z]);
        bin[z]=NULL;
      }
      free(bin);
      bin=NULL;
      break;
    default:
      break;
    }
    return (-1);
  }

  theXX = (float ***)theTensor->imxx.array;
  theYY = (float ***)theTensor->imyy.array;
  theZZ = (float ***)theTensor->imzz.array;
  theXY = (float ***)theTensor->imxy.array;
  theYZ = (float ***)theTensor->imyz.array;
  theXZ = (float ***)theTensor->imxz.array;

  theVP1 = (float ***)theTensor->imvp1.array;
  theVP2 = (float ***)theTensor->imvp2.array;
  theVP3 = (float ***)theTensor->imvp3.array;

  theTheta1 = (float ***)theTensor->imtheta1.array;
  thePhi1 = (float ***)theTensor->imphi1.array;
  //theTheta2 = (float ***)theTensor->imtheta2.array;
  //thePhi2 = (float ***)theTensor->imphi2.array;
  //theTheta3 = (float ***)theTensor->imtheta3.array;
  //thePhi3 = (float ***)theTensor->imphi3.array;

  theZEROS = (unsigned char ***)theTensor->iszero.array;

  /* Initialisation de l'image tenseur a partir d'imBin */
  if (_verbose_)
    fprintf(stdout,"%s: initialisation du tenseur a partir d'imBin...\n", proc);

  n=0;
  m=0;

  for (z=0;z<dimz;z++)
  {
    if (_verbose_ && (dimz<=30 || z%(dimz/30)==0))
      fprintf(stdout,"\tFrame #%d/%d (%d%%)... \n",
          z+1, dimz, (int)(100*(float)z/(float)dimz));
    for (y=0;y<dimy;y++)
    for (x=0;x<dimx;x++)
    {
      val = bin[z][y][x];
      if ( val < FLTZERO)
      {
        /*  le voxel(x,y,z) ne vote pas ni ne reçoit de votes  */
        theXX[z][y][x] = 0;
        theYY[z][y][x] = 0;
        theZZ[z][y][x] = 0;
        theXY[z][y][x] = 0;
        theXZ[z][y][x] = 0;
        theYZ[z][y][x] = 0;
        theZEROS[z][y][x] = (unsigned char)1;
        n++;
        continue;
      }
      if (initHessian==NONE)
      {
      /* Tenseur Ball */
      theXX[z][y][x] = val;
      theYY[z][y][x] = val;
      theZZ[z][y][x] = val;
      theXY[z][y][x] = 0;
      theXZ[z][y][x] = 0;
      theYZ[z][y][x] = 0;

      theVP1[z][y][x] = val;
      theVP2[z][y][x] = val;
      theVP3[z][y][x] = val;
      }
      if (initHessian!=NONE)
      {
        /* Tenseur Test */
        theTheta1[z][y][x] = imTht[z][y][x];
        thePhi1[z][y][x]   = imPhi[z][y][x];

        theVP1[z][y][x] = 0;
        theVP2[z][y][x] = val;
        theVP3[z][y][x] = val;
      }
      theZEROS[z][y][x] = (unsigned char)0;
      m++;
    }
  }

  if (initHessian != NONE)
    MT_InitTensorTest(theTensor);

  if (writeImages)
  {
    if (_verbose_)
    {
      sprintf( name, "%s.init", parName );
      fprintf(stdout,
      "Ecriture du champ de tenseurs initial dans %s.inr...\n", name);
    }
    VT_Write3DtensorWithName( theTensor, name );
  }


  if (_verbose_)
  {
    fprintf(stdout, "%s : nombre de zeros = %d \tnombre de non-zeros = %d \n",
        proc, n,m);
  }

  /* Votes DENSE */

  if (_verbose_)
    fprintf(stdout, "\tDENSE voting...\n");

  MT_ReinitTensorImage(theTensor);

  /* Cumul des votes (seuls imxx,yy,zz,xy,xz,yz et iszero sont affectes) */
  for (z=0;z<dimz;z++)
  {
    if (_verbose_)
      if (dimz<=30 || z%(dimz/30)==0)
        fprintf(stdout, "\t\tTreating frame #%d/%d (%d%%)...\n",
            z+1, dimz, (int)(100*(float)z/(float)dimz));
    for (y=0;y<dimy;y++)
    for (x=0;x<dimx;x++)
    {
      if ( theZEROS[z][y][x] != 0)
      {
        /*  le voxel(x,y,z) ne vote pas  */
        continue;
      }

      lambda1 = theVP1[z][y][x];
      lambda2 = theVP2[z][y][x];
      //lambda3 = theVP3[z][y][x];

      if (lambda2-lambda1 > FLTZERO)  
      {
        theta = theTheta1[z][y][x];
        phi = thePhi1[z][y][x];
        n = nearestAngle(theta, phi, angles, Nangles);
        MT_AddField(theTensor, tfields[n], lambda2, x, y, z, DENSE);
      }
    }
  }

  /* Conversion de la hessienne en vp, angles, etc... */
  if (_verbose_)
    fprintf(stdout, "\tConversion des tenseurs en valeurs/vecteurs propres, \
        etc...\n");

  for (z=0;z<dimz;z++)
  for (y=0;y<dimy;y++)
  for (x=0;x<dimx;x++)
  {
    if ( theZEROS[z][y][x] == 1)
    {
      /*  le voxel(x,y,z) contient un tenseur nul */
      continue;
    }
    if (setTensor(theTensor, x,y,z, DENSE) !=1)
    {
      if (_verbose_){
        fprintf(stderr, "%s: probleme lors du calcul des parametres des tenseurs\n", proc);
	    fprintf(stderr, "setTensor #%d\n",x+y*dimx+z*dimy*dimx); 
	    fprintf(stderr, "Hessian[%d][%d][%d]=\n%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f\n", 
			  z, y, x,
			  theXX[z][y][x],theXY[z][y][x],theXZ[z][y][x],
			  theXY[z][y][x],theYY[z][y][x],theYZ[z][y][x],
			  theXZ[z][y][x],theYZ[z][y][x],theZZ[z][y][x]);
	  }
      for (n=0;n<Nangles;n++){
        free(angles[n]);
        angles[n]=NULL;
        if (initHessian != NONE)
		  VT_Free3Dtensor(tfields+n);
      }
	  if (initHessian != NONE) {
        free(tfields);
		tfields=NULL;
	  }

      free(angles);
      angles = NULL;
      VT_Free3Dtensor(theTensor);
      switch (imBin[0]->type)
      {
      case UCHAR :
        for(z=0;z<dimz;z++)
        {
          for(y=0;y<dimy;y++)
          {
            free(bin[z][y]);
            bin[z][y]=NULL;
          }
          free(bin[z]);
          bin[z]=NULL;
        }
        free(bin);
        bin=NULL;
        break;
      default:
        break;
      }
      return(-1);
    }

  }


  if (_verbose_)
    fprintf(stdout, "%s: Liberation de la memoire allouee...\n", proc);

  for (n=0;n<Nangles;n++){
    free(angles[n]);
    angles[n]=NULL;
    if (initHessian != NONE)
	  VT_Free3Dtensor(tfields+n);
  }
  if (initHessian != NONE){
    free(tfields);
	tfields=NULL;
  }
  free(angles);
  angles = NULL;
  switch (imBin[0]->type)
  {
  case UCHAR :
    for(z=0;z<dimz;z++)
    {
      for(y=0;y<dimy;y++)
      {
        free(bin[z][y]);
        bin[z][y]=NULL;
      }
      free(bin[z]);
      bin[z]=NULL;
    }
    free(bin);
    bin=NULL;
    break;
  default:
    break;
  }

  
  return(1);
}

////////////////////////////////////////////////////////////////////////////////

void MT_ComputeTensorSurfaceExtrema( vt_3Dtensor *imTensor,
                          vt_image *imExt, double zfact )
{
  int x, y, z;

  float ***theExt = (float***)(imExt->array);

  //float ***theVP1 = (float***)(imTensor->imvp1.array);
  float ***theVP2 = (float***)(imTensor->imvp2.array);
  float ***theVP3 = (float***)(imTensor->imvp3.array);
  float ***theTht1 = (float***)(imTensor->imtheta1.array);
  float ***theTht2 = (float***)(imTensor->imtheta2.array);
  float ***theTht3 = (float***)(imTensor->imtheta3.array);
  float ***thePhi1 = (float***)(imTensor->imphi1.array);
  float ***thePhi2 = (float***)(imTensor->imphi2.array);
  float ***thePhi3 = (float***)(imTensor->imphi3.array);
  unsigned char ***theZeros = (unsigned char***)(imTensor->iszero.array);



  double v1[3], v2[3], v3[3];

  double rep, rn[2];

  double theta1, phi1, theta2, phi2, theta3, phi3;

  double vx=0.0, vy=0.0, vz=0.0;
  double norm;

  double rx, ry, rz;
  double dx, dy, dz;
  int j, ix, iy, iz;
  double dxdy, dxdz, dydz, dxdydz;
  double v4, v5, v6;


  for ( z=0; z<imExt->dim.z; z++ )
  for ( y=0; y<imExt->dim.y; y++ )
  for ( x=0; x<imExt->dim.x; x++ ) {

    theta3 = (double)theTht3[z][y][x];
    phi3   = (double)thePhi3[z][y][x];
    theta2 = (double)theTht2[z][y][x];
    phi2   = (double)thePhi2[z][y][x];
    theta1 = (double)theTht1[z][y][x];
    phi1   = (double)thePhi1[z][y][x];

    theExt[z][y][x] = 0.0;
    if ( theZeros[z][y][x] != 0 ) continue;

    rep = (double)(theVP3[z][y][x]-theVP2[z][y][x]);



    /* v3 est le vecteur normal a la membrane
       v1 et v2 sont deux vecteurs du plan tangent
    */
    SphericalAnglesToUnitVector( theta3, phi3, v3);
    SphericalAnglesToUnitVector( theta2, phi2, v2);
    SphericalAnglesToUnitVector( theta1, phi1, v1);

    for (j=0;j<2;j++) {
      switch(j) {
      default:
      case 0 :
  vx = zfact*v3[0];    vy = zfact*v3[1];    vz = v3[2]; break;
      case 1 :
  vx = -zfact*v3[0];  vy = -zfact*v3[1];  vz = -v3[2];  break;
      }
      norm=sqrt(vx*vx+vy*vy+vz*vz);
      vx/=norm; vy/=norm; vz/=norm;

      rx = x + vx;
      ix = (int)rx;
      if ( rx <= 0.0 || ix >= imExt->dim.x-1 ) {
        j = 3;
        continue;
      }
      ry = y + vy;
      iy = (int)ry;
      if ( ry <= 0.0 || iy >= imExt->dim.y-1 ) {
        j = 3;
        continue;
      }
      rz = z + vz;
      iz = (int)rz;
      if ( rz <= 0.0 || iz >= imExt->dim.z-1 ) {
        j = 3;
        continue;
      }


      dx = rx - ix;
      dy = ry - iy;
      dz = rz - iz;

      dxdy = dx*dy;
      dxdz = dx*dz;
      dydz = dy*dz;
      dxdydz = dxdy*dz;

      v6 = dxdz-dxdydz;
      v5 = dxdy-dxdydz;
      v4 = dx-dxdy-v6;

      rn[j] = 0;
      rn[j]+= dxdydz       *(theVP3[iz+1][iy+1][ix+1]-theVP2[iz+1][iy+1][ix+1]);
      rn[j]+= (dydz-dxdydz)*(theVP3[iz+1][iy+1][ix  ]-theVP2[iz+1][iy+1][ix  ]);
      rn[j]+= v6           *(theVP3[iz+1][iy  ][ix+1]-theVP2[iz+1][iy  ][ix+1]);
      rn[j]+= (dz-dydz-v6) *(theVP3[iz+1][iy  ][ix  ]-theVP2[iz+1][iy  ][ix  ]);
      rn[j]+= v5           *(theVP3[iz  ][iy+1][ix+1]-theVP2[iz  ][iy+1][ix+1]);
      rn[j]+= (dy-dydz-v5) *(theVP3[iz  ][iy+1][ix  ]-theVP2[iz  ][iy+1][ix  ]);
      rn[j]+= v4           *(theVP3[iz  ][iy  ][ix+1]-theVP2[iz  ][iy  ][ix+1]);
      rn[j]+= (1-dy-dz+dydz-v4)*(theVP3[iz  ][iy][ix]-theVP2[iz  ][iy  ][ix  ]);

      if ( rn[j] > rep ) {
        j = 3;
        continue;
      }

    } /* fin de la boucle suivant la normale à la membrane */

    if ( rn[0] == rep && rn[1] == rep ) {
      continue;
    }

    if ( j < 3 ) theExt[z][y][x] = rep;
  } //x,y,z
}


void MT_ComputeTensorLineExtrema( vt_3Dtensor *imTensor,
                          vt_image *imExt, double zfact )
{
  int x, y, z;
  
  float ***theExt = (float***)(imExt->array);
  float ***theVP1 = (float***)(imTensor->imvp1.array);
  float ***theVP2 = (float***)(imTensor->imvp2.array);
  float ***theTht = (float***)(imTensor->imtheta1.array);
  float ***thePhi = (float***)(imTensor->imphi1.array);

  double v1[3], v2[3], v3[3];
  
  double rep, r;

  double theta, phi;

  double vx=0.0, vy=0.0, vz=0.0;
  
  double c = sqrt(2.0)/2.0;

  double rx, ry, rz;
  double dx, dy, dz;
  int i, ix, iy, iz;
  double dxdy, dxdz, dydz, dxdydz;
  double v4, v5, v6;





  for ( z=0; z<imExt->dim.z; z++ ) 
  for ( y=0; y<imExt->dim.y; y++ ) 
  for ( x=0; x<imExt->dim.x; x++ ) {

    theta = (double)theTht[z][y][x];
    phi   = (double)thePhi[z][y][x];
    
    theExt[z][y][x] = 0.0;

    rep = theVP2[z][y][x]-theVP1[z][y][x];
	if ( rep <= 0.0 ) continue;


    /* v1 est le vecteur directeur du vaisseau
       v2 et v3 sont deux vecteurs orthogonaux
    */
    SphericalAnglesToUnitsVectors( theta, phi, v1, v2, v3 );
    
    for ( i=0; i<8; i++ ) {

      switch( i ) {
      default :
      case 0 :
	vx = v2[0];    vy = v2[1];    vz = v2[2];    break;
      case 1 :
	vx = -v2[0];   vy = -v2[1];   vz = -v2[2];   break;
      case 2 :
	vx = v3[0];    vy = v3[1];    vz = v3[2];    break;
      case 3 :
	vx = -v3[0];   vy = -v3[1];   vz = -v3[2];   break;
      case 4 :
	vx =  c * v2[0] + c * v3[0];
	vy =  c * v2[1] + c * v3[1];
	vz =  c * v2[2] + c * v3[2];
	break;
      case 5 :
	vx = -c * v2[0] + c * v3[0];
	vy = -c * v2[1] + c * v3[1];
	vz = -c * v2[2] + c * v3[2];
	break;
       case 6 :
	vx =  c * v2[0] - c * v3[0];
	vy =  c * v2[1] - c * v3[1];
	vz =  c * v2[2] - c * v3[2];
	break;
      case 7 :
	vx = -c * v2[0] - c * v3[0];
	vy = -c * v2[1] - c * v3[1];
	vz = -c * v2[2] - c * v3[2];
	break;
      }

      rx = x + vx;
      ix = (int)rx;
      if ( rx <= 0.0 || ix >= imExt->dim.x-1 ) {
	i = 10;
	continue;
      }
      ry = y + vy;
      iy = (int)ry;
      if ( ry <= 0.0 || iy >= imExt->dim.y-1 ) {
	i = 10;
	continue;
      }
      rz = z + vz;
      iz = (int)rz;
      if ( rz <= 0.0 || iz >= imExt->dim.z-1 ) {
	i = 10;
	continue;
      }


      /* si on veut que tous les points autour
	 aient une reponse positive
	 C'est tres (trop?) restrictif
      */
      if ( 0 ) {
	if ( theVP2[iz+1][iy+1][ix+1]-theVP1[iz+1][iy+1][ix+1] <= 0.0 ||
	     theVP2[iz+1][iy+1][ix  ]-theVP1[iz+1][iy+1][ix  ] <= 0.0 ||
	     theVP2[iz+1][iy  ][ix+1]-theVP1[iz+1][iy  ][ix+1] <= 0.0 ||
	     theVP2[iz+1][iy  ][ix  ]-theVP1[iz+1][iy  ][ix  ] <= 0.0 ||
	     theVP2[iz  ][iy+1][ix+1]-theVP1[iz  ][iy+1][ix+1] <= 0.0 ||
	     theVP2[iz  ][iy+1][ix  ]-theVP1[iz  ][iy+1][ix  ] <= 0.0 ||
	     theVP2[iz  ][iy  ][ix+1]-theVP1[iz  ][iy  ][ix+1] <= 0.0 ||
	     theVP2[iz  ][iy  ][ix  ]-theVP1[iz  ][iy  ][ix  ] <= 0.0 ) {
	  i = 10;
	  continue;
	}
      }
      
      

      dx = rx - ix;
      dy = ry - iy;
      dz = rz - iz;

      dxdy = dx*dy;
      dxdz = dx*dz;
      dydz = dy*dz;
      dxdydz = dxdy*dz;

      v6 = dxdz-dxdydz;
      v5 = dxdy-dxdydz;
      v4 = dx-dxdy-v6;

      r = 0;
      r += dxdydz        * (theVP2[iz+1][iy+1][ix+1]-theVP1[iz+1][iy+1][ix+1]);
      r += (dydz-dxdydz) * (theVP2[iz+1][iy+1][ix  ]-theVP1[iz+1][iy+1][ix  ]);
      r += v6            * (theVP2[iz+1][iy  ][ix+1]-theVP1[iz+1][iy  ][ix+1]);
      r += (dz-dydz-v6)  * (theVP2[iz+1][iy  ][ix  ]-theVP1[iz+1][iy  ][ix  ]);
      r += v5            * (theVP2[iz  ][iy+1][ix+1]-theVP1[iz  ][iy+1][ix+1]);
      r += (dy-dydz-v5)  * (theVP2[iz  ][iy+1][ix  ]-theVP1[iz  ][iy+1][ix  ]);
      r += v4            * (theVP2[iz  ][iy  ][ix+1]-theVP1[iz  ][iy  ][ix+1]);
      r += (1-dy-dz+dydz-v4) * (theVP2[iz ][iy ][ix ]-theVP1[iz ][iy  ][ix  ]);
      

      
      if ( r >= rep ) {
	i = 10;
	continue;
      }
      
    } /* fin de la boucle sur les 4 directions
       */
    if ( i < 10 ) theExt[z][y][x] = rep;
  }
}


void MT_ComputeTensorBallExtrema( vt_3Dtensor *imTensor,
                          vt_image *imExt, double zfact )
{
	//TODO
  int x, y, z;
  
  float ***theExt = (float***)(imExt->array);
  float ***theVP1 = (float***)(imTensor->imvp1.array);

  float rep;

  int i;

  int dx[]={-1, -1, -1, -1, -1, -1, -1, -1, -1, 
			 0,  0,  0,  0,  0,  0,  0,  0,
			 1,  1,  1,  1,  1,  1,  1,  1,  1};

  int dy[]={-1, -1, -1,  0,  0,  0,  1,  1,  1, 
			-1, -1, -1,  0,  0,  1,  1,  1,
			-1, -1, -1,  0,  0,  0,  1,  1,  1};

  int dz[]={-1,  0,  1, -1,  0,  1, -1,  0,  1, 
			-1,  0,  1, -1,  1, -1,  0,  1,
			-1,  0,  1, -1,  0,  1, -1,  0,  1};

  for ( z=1; z<imExt->dim.z-1; z++ ) 
  for ( y=1; y<imExt->dim.y-1; y++ ) 
  for ( x=1; x<imExt->dim.x-1; x++ ) {
	rep=theVP1[z][y][x];
	
	if (rep<=0.0) continue;
	
	
	for (i=0; i<26; i++) {
	  if ( rep > theVP1[z+dz[i]][y+dy[i]][x+dx[i]] )
		continue;
	  i=27;
	}
	if ( i == 26 )
	  theExt[z][y][x] = rep;
  }
}







void MT_ReinitTensorImage(vt_3Dtensor *theTensor)
{
  /* Reinitialise seulement les champs XX, YY, ZZ, XY, XZ et YZ */
  int x,y,z;
  int dimx = theTensor->imxx.dim.x;
  int dimy = theTensor->imxx.dim.y;
  int dimz = theTensor->imxx.dim.z;

  float ***theXX = (float ***)theTensor->imxx.array;
  float ***theYY = (float ***)theTensor->imyy.array;
  float ***theZZ = (float ***)theTensor->imzz.array;
  float ***theXY = (float ***)theTensor->imxy.array;
  float ***theYZ = (float ***)theTensor->imyz.array;
  float ***theXZ = (float ***)theTensor->imxz.array;

  for (z=0; z<dimz; z++)
  for (y=0; y<dimy; y++)
  for (x=0; x<dimx; x++)
  {
    theXX[z][y][x] = 0;
    theYY[z][y][x] = 0;
    theZZ[z][y][x] = 0;
    theXY[z][y][x] = 0;
    theXZ[z][y][x] = 0;
    theYZ[z][y][x] = 0;
  }
  return;
}

void MT_InitTensorFromAngles(vt_3Dtensor *theTensor)
{
  int x,y,z;
  int dimx = theTensor->imxx.dim.x;
  int dimy = theTensor->imxx.dim.y;
  int dimz = theTensor->imxx.dim.z;

  float ***theXX = (float ***)theTensor->imxx.array;
  float ***theYY = (float ***)theTensor->imyy.array;
  float ***theZZ = (float ***)theTensor->imzz.array;
  float ***theXY = (float ***)theTensor->imxy.array;
  float ***theYZ = (float ***)theTensor->imyz.array;
  float ***theXZ = (float ***)theTensor->imxz.array;
  float ***theTheta3 = (float ***)theTensor->imtheta3.array;
  float ***thePhi3 = (float ***)theTensor->imphi3.array;
  float ***theVP3 = (float ***)theTensor->imvp3.array;
  float ***theTheta2 = (float ***)theTensor->imtheta2.array;
  float ***thePhi2 = (float ***)theTensor->imphi2.array;
  float ***theVP2 = (float ***)theTensor->imvp2.array;
  float ***theTheta1 = (float ***)theTensor->imtheta1.array;
  float ***thePhi1 = (float ***)theTensor->imphi1.array;
  float ***theVP1 = (float ***)theTensor->imvp1.array;

  unsigned char ***zeros = (unsigned char ***)theTensor->iszero.array;

  double v[3], theta, phi;

  for (z=0; z<dimz; z++)
  for (y=0; y<dimy; y++)
  for (x=0; x<dimx; x++)
  {
    if (zeros[z][y][x] == 1)
      continue;

    theta = (double)theTheta3[z][y][x];
    phi = (double)thePhi3[z][y][x];
    SphericalAnglesToUnitVector(theta, phi, v);
    theXX[z][y][x] = (float)v[0]*v[0]*theVP3[z][y][x];
    theYY[z][y][x] = (float)v[1]*v[1]*theVP3[z][y][x];
    theZZ[z][y][x] = (float)v[2]*v[2]*theVP3[z][y][x];
    theXY[z][y][x] = (float)v[0]*v[1]*theVP3[z][y][x];
    theXZ[z][y][x] = (float)v[0]*v[2]*theVP3[z][y][x];
    theYZ[z][y][x] = (float)v[1]*v[2]*theVP3[z][y][x];

	if (theVP2[z][y][x] <= FLTZERO)
	  continue;

    theta = (double)theTheta2[z][y][x];
    phi = (double)thePhi2[z][y][x];
    SphericalAnglesToUnitVector(theta, phi, v);
    theXX[z][y][x] += (float)v[0]*v[0]*theVP2[z][y][x];
    theYY[z][y][x] += (float)v[1]*v[1]*theVP2[z][y][x];
    theZZ[z][y][x] += (float)v[2]*v[2]*theVP2[z][y][x];
    theXY[z][y][x] += (float)v[0]*v[1]*theVP2[z][y][x];
    theXZ[z][y][x] += (float)v[0]*v[2]*theVP2[z][y][x];
    theYZ[z][y][x] += (float)v[1]*v[2]*theVP2[z][y][x];

	if (theVP1[z][y][x] <= FLTZERO)
	  continue;

    theta = (double)theTheta1[z][y][x];
    phi = (double)thePhi1[z][y][x];
    SphericalAnglesToUnitVector(theta, phi, v);
    theXX[z][y][x] += (float)v[0]*v[0]*theVP1[z][y][x];
    theYY[z][y][x] += (float)v[1]*v[1]*theVP1[z][y][x];
    theZZ[z][y][x] += (float)v[2]*v[2]*theVP1[z][y][x];
    theXY[z][y][x] += (float)v[0]*v[1]*theVP1[z][y][x];
    theXZ[z][y][x] += (float)v[0]*v[2]*theVP1[z][y][x];
    theYZ[z][y][x] += (float)v[1]*v[2]*theVP1[z][y][x];
  }

  return;
}


void MT_InitTensorTest(vt_3Dtensor *theTensor)	
{
  int x,y,z;
  int dimx = theTensor->imxx.dim.x;
  int dimy = theTensor->imxx.dim.y;
  int dimz = theTensor->imxx.dim.z;

  float ***theXX = (float ***)theTensor->imxx.array;
  float ***theYY = (float ***)theTensor->imyy.array;
  float ***theZZ = (float ***)theTensor->imzz.array;
  float ***theXY = (float ***)theTensor->imxy.array;
  float ***theYZ = (float ***)theTensor->imyz.array;
  float ***theXZ = (float ***)theTensor->imxz.array;
  float ***theVP3 = (float ***)theTensor->imvp3.array;
  float ***theTheta1 = (float ***)theTensor->imtheta1.array;
  float ***thePhi1 = (float ***)theTensor->imphi1.array;
  
  unsigned char ***zeros = (unsigned char ***)theTensor->iszero.array;

  double v[3], theta, phi;

  for (z=0; z<dimz; z++)
  for (y=0; y<dimy; y++)
  for (x=0; x<dimx; x++)
  {
    if (zeros[z][y][x] == 1)
      continue;

    theta = (double)theTheta1[z][y][x];
    phi = (double)thePhi1[z][y][x];
    SphericalAnglesToUnitVector(theta, phi, v);
    
    theXX[z][y][x] = (float)(1 - v[0]*v[0])*theVP3[z][y][x];
    theYY[z][y][x] = (float)(1 - v[1]*v[1])*theVP3[z][y][x];
    theZZ[z][y][x] = (float)(1 - v[2]*v[2])*theVP3[z][y][x];
    theXY[z][y][x] = (float)- v[0]*v[1]*theVP3[z][y][x];
    theXZ[z][y][x] = (float)- v[0]*v[2]*theVP3[z][y][x];
    theYZ[z][y][x] = (float)- v[1]*v[2]*theVP3[z][y][x];

  }

  return;
}

int nearestAngle(double theta, double phi, double **angles, int Nangles)
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



int MT_AddField(vt_3Dtensor *tensorImg, vt_3Dtensor field,
    double coef, int x, int y, int z, enumVote typeVote)
{
  if (coef < FLTZERO)
    return(1);


  float ***fieldXX = (float***)(field.imxx.array);
  float ***fieldYY = (float***)(field.imyy.array);
  float ***fieldZZ = (float***)(field.imzz.array);
  float ***fieldXY = (float***)(field.imxy.array);
  float ***fieldXZ = (float***)(field.imxz.array);
  float ***fieldYZ = (float***)(field.imyz.array);

  unsigned char ***fieldZeros = (unsigned char***)(field.iszero.array);


  float ***theXX = (float***)(tensorImg->imxx.array);
  float ***theYY = (float***)(tensorImg->imyy.array);
  float ***theZZ = (float***)(tensorImg->imzz.array);
  float ***theXY = (float***)(tensorImg->imxy.array);
  float ***theXZ = (float***)(tensorImg->imxz.array);
  float ***theYZ = (float***)(tensorImg->imyz.array);

  unsigned char ***theZeros = (unsigned char***)(tensorImg->iszero.array);


  int xf0, yf0, zf0;
  int x0, y0, z0, x1, y1, z1;
  int i, j, k;
  int ifield, jfield, kfield;

  int sf[3]={(int)field.imxx.dim.x,(int)field.imxx.dim.y,(int)field.imxx.dim.z};

  int hsf[3];

  int xmax = tensorImg->imxx.dim.x;
  int ymax = tensorImg->imxx.dim.y;
  int zmax = tensorImg->imxx.dim.z;

  hsf[0]=sf[0]/2;
  hsf[1]=sf[1]/2;
  hsf[2]=sf[2]/2;


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
        if (fieldZeros[kfield][jfield][ifield]==0)
          if (typeVote == DENSE || theZeros[k][j][i]==0)
          {
            theZeros[k][j][i]=0;
            theXX[k][j][i]+=(float)coef*fieldXX[kfield][jfield][ifield];
            theYY[k][j][i]+=(float)coef*fieldYY[kfield][jfield][ifield];
            theZZ[k][j][i]+=(float)coef*fieldZZ[kfield][jfield][ifield];
            theXY[k][j][i]+=(float)coef*fieldXY[kfield][jfield][ifield];
            theXZ[k][j][i]+=(float)coef*fieldXZ[kfield][jfield][ifield];
            theYZ[k][j][i]+=(float)coef*fieldYZ[kfield][jfield][ifield];
          }
        ifield++;
      }
      jfield++;
    }
    kfield++;
  }

  return(1);
}


int setTensor(vt_3Dtensor *tensorImg, int x, int y, int z, enumVote mode)
{
  char *proc = "setTensor";
  float ***theXX = (float***)(tensorImg->imxx.array);
  float ***theYY = (float***)(tensorImg->imyy.array);
  float ***theZZ = (float***)(tensorImg->imzz.array);
  float ***theXY = (float***)(tensorImg->imxy.array);
  float ***theXZ = (float***)(tensorImg->imxz.array);
  float ***theYZ = (float***)(tensorImg->imyz.array);

  float ***theVP1 = (float***)(tensorImg->imvp1.array);
  float ***theVP2 = (float***)(tensorImg->imvp2.array);
  float ***theVP3 = (float***)(tensorImg->imvp3.array);


  float ***theTht1 = (float***)(tensorImg->imtheta1.array);
  float ***theTht2 = (float***)(tensorImg->imtheta2.array);
  float ***theTht3 = (float***)(tensorImg->imtheta3.array);

  float ***thePhi1 = (float***)(tensorImg->imphi1.array);
  float ***thePhi2 = (float***)(tensorImg->imphi2.array);
  float ***thePhi3 = (float***)(tensorImg->imphi3.array);

  unsigned char ***theZeros = (unsigned char ***) (tensorImg->iszero.array);

  double hessien[9];
  double valprop[3], vecprop[9];
  double v[3], theta, phi;

  float oldVPsum, newVPsum;

  hessien[0] = (double)theXX[z][y][x];
  hessien[1] = hessien[3] = (double)theXY[z][y][x];
  hessien[2] = hessien[6] = (double)theXZ[z][y][x];
  hessien[4] = (double)theYY[z][y][x];
  hessien[5] = hessien[7] = (double)theYZ[z][y][x];
  hessien[8] = (double)theZZ[z][y][x];

  if ( _ComputeEigensOfSymetricSquareMatrix( hessien, valprop, vecprop, 3 )
       != 1 ) {
    if (_verbose_)
      fprintf(stderr, "%s: error in computing eigens\n", proc);
    return( -1 );
  }
  if ( _SortEigensInAbsIncreasingOrder( valprop, vecprop, 3 ) != 1 ) {

    if (_verbose_)
      fprintf(stderr, "%s: error in sorting eigens\n", proc);
    return( -1 );
  }

  /* Les valeurs et vecteurs propres sont tries par fabs croissante */
  if (mode == SPARSE) {
    oldVPsum=theVP1[z][y][x]+theVP2[z][y][x]+theVP3[z][y][x];
    newVPsum=(float)(valprop[0]+valprop[1]+valprop[2]);
  }

  theVP1[z][y][x]=(float)valprop[0];
  v[0]=vecprop[0]; v[1]=vecprop[3]; v[2]=vecprop[6];
  UnitVectorToSphericalAngles( v, &theta, &phi );
  theTht1[z][y][x]=(float)theta;
  thePhi1[z][y][x]=(float)phi;

  theVP2[z][y][x]=(float)valprop[1];
  v[0]=vecprop[1]; v[1]=vecprop[4]; v[2]=vecprop[7];
  UnitVectorToSphericalAngles( v, &theta, &phi );
  theTht2[z][y][x]=(float)theta;
  thePhi2[z][y][x]=(float)phi;

  theVP3[z][y][x]=(float)valprop[2];
  v[0]=vecprop[2]; v[1]=vecprop[5]; v[2]=vecprop[8];
  UnitVectorToSphericalAngles( v, &theta, &phi );
  theTht3[z][y][x]=(float)theta;
  thePhi3[z][y][x]=(float)phi;

  if (valprop[2]>FLTZERO)
  {
    theZeros[z][y][x] = (unsigned char)0;
    if (mode == SPARSE)
    {
        /* On "normalise" : sum(vp avant vote)=sum(vp apres vote) */
        theVP1[z][y][x]=theVP1[z][y][x]*oldVPsum/newVPsum;
        theVP2[z][y][x]=theVP2[z][y][x]*oldVPsum/newVPsum;
        theVP3[z][y][x]=theVP3[z][y][x]*oldVPsum/newVPsum;

        theXX[z][y][x]=theXX[z][y][x]*oldVPsum/newVPsum;
        theYY[z][y][x]=theYY[z][y][x]*oldVPsum/newVPsum;
        theZZ[z][y][x]=theZZ[z][y][x]*oldVPsum/newVPsum;
        theXY[z][y][x]=theXY[z][y][x]*oldVPsum/newVPsum;
        theXZ[z][y][x]=theXZ[z][y][x]*oldVPsum/newVPsum;
        theYZ[z][y][x]=theYZ[z][y][x]*oldVPsum/newVPsum;
    }

  }
  else
  {
      theZeros[z][y][x] = (unsigned char)1;
  }

  return(1);
}


/////////////////////////////////////

int MT_Compute3DTestFields( vt_3Dtensor **tfields, int *dimFields, 
          double **angles, int Nangles, double r, double alpha)
{

  int n,m;
  char *proc = "MT_Compute3DTestFields";
  double angle[2];
  char name[256];  
  int t = (int) FLOAT; 
  vt_3Dtensor *field;

  (*tfields) =  malloc(Nangles*sizeof(vt_3Dtensor));

  if((*tfields)==NULL)
  {
    if (_verbose_)
      fprintf(stderr, "%s: allocation de (*tfield) echouee", proc);
    return(-1);
  }

  for (n=0; n<Nangles; n++) {
    sprintf( name, "%s.%d", "testfield", n );
    field = (*tfields+n);
    if ( VT_Alloc3Dtensor(field, name, dimFields[0], dimFields[1], dimFields[2],
        t ) != 1 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: unable to allocate 3D tensor images\n", proc );
      for (m=0;m<n;m++){
        field = (*tfields+m);
        VT_Free3Dtensor( field );
      }
      free((*tfields));
      tfields[0]=NULL;
      return( -1 );
      
    }
  }
  
  for (n=0; n<Nangles; n++) {
    angle[0]=angles[n][0];
    angle[1]=angles[n][1];
    field = (*tfields+n);
    if (MT_Compute3DTestField(field, dimFields, angle, r, alpha) != 1)
    {
      if (_verbose_)
        fprintf(stderr, "%s: error in computating test field #%d\n", proc, n);
      for (m=0;m<Nangles;m++)
      {
        field = (*tfields+m);
        VT_Free3Dtensor( field );
      }
      free((*tfields));
      (*tfields)=NULL;
      return( -1 );
    }
  }

  return( 1 );
}

/////////////////////////////////

int MT_Compute3DStickFields( vt_3Dtensor **sfields, int *dimFields, 
          double **angles, double *theCoeffs, int Nangles, enumTVmode mode)
{

  int n,m;
  char *proc = "MT_Compute3DStickFields";
  double angle[2];
  char name[256];  
  int t = (int) FLOAT; 
  vt_3Dtensor *field;

  (*sfields) =  malloc(Nangles*sizeof(vt_3Dtensor));

  if((*sfields)==NULL)
  {
    if (_verbose_)
      fprintf(stderr, "%s: allocation de (*sfield) echouee", proc);
    return(-1);
  }

  for (n=0; n<Nangles; n++) {
    sprintf( name, "%s.%d", "stickfield", n );
    field = (*sfields+n);
    if ( VT_Alloc3Dtensor(field, name, dimFields[0], dimFields[1], dimFields[2],
        t ) != 1 ) {
      if ( _verbose_ ) 
        fprintf( stderr, "%s: unable to allocate 3D tensor images\n", proc );
      for (m=0;m<n;m++){
        field = (*sfields+m);
        VT_Free3Dtensor( field );
      }
      free((*sfields));
      sfields[0]=NULL;
      return( -1 );
      
    }
  }
  
  for (n=0; n<Nangles; n++) {
    angle[0]=angles[n][0];
    angle[1]=angles[n][1];
    field = (*sfields+n);
    if (MT_Compute3DStickField(field, dimFields, angle, theCoeffs, mode) != 1)
    {
      if (_verbose_)
        fprintf(stderr, "%s: error in computating stick field #%d\n", proc, n);
      for (m=0;m<Nangles;m++)
      {
        field = (*sfields+m);
        VT_Free3Dtensor( field );
      }
      free((*sfields));
      (*sfields)=NULL;
      return( -1 );
    }
  }

  return( 1 );
}



///////////////////

int MT_Compute3DTestField( vt_3Dtensor *tfield, int *dimFields, 
          double *angle, double ray, double Alpha)
{

  int x,y,z;
  double pi = 3.14159265358979323846;

  int pt[3] = {dimFields[0]/2, dimFields[1]/2, dimFields[2]/2};
  double r2=ray*ray;
  double norme2;
  double n;
  double n1[3], m1[3];
  double r[3];
  double alpha, beta;
  double cosalpha;
  double dot;
  
  float ***theXX = (float***)(tfield->imxx.array);
  float ***theYY = (float***)(tfield->imyy.array);
  float ***theZZ = (float***)(tfield->imzz.array);
  float ***theXY = (float***)(tfield->imxy.array);
  float ***theXZ = (float***)(tfield->imxz.array);
  float ***theYZ = (float***)(tfield->imyz.array);

  float ***theVP1 = (float***)(tfield->imvp1.array);
  float ***theVP2 = (float***)(tfield->imvp2.array);
  float ***theVP3 = (float***)(tfield->imvp3.array);


  float ***theTht1 = (float***)(tfield->imtheta1.array);
  float ***theTht2 = (float***)(tfield->imtheta2.array);
  float ***theTht3 = (float***)(tfield->imtheta3.array);

  float ***thePhi1 = (float***)(tfield->imphi1.array);
  float ***thePhi2 = (float***)(tfield->imphi2.array);
  float ***thePhi3 = (float***)(tfield->imphi3.array);
  
  unsigned char ***theZeros = (unsigned char***)(tfield->iszero.array);
  
  SphericalAnglesToUnitVector( angle[0], angle[1], n1);
  
  cosalpha = cos(Alpha);
  if (cosalpha < 0) cosalpha=-cosalpha;
  
  for (z=0;z<dimFields[2];z++)
  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++) 
  {
	r[0]=(double)x-pt[0];
	r[1]=(double)y-pt[1];
	r[2]=(double)z-pt[2];
	
	norme2=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
	if (norme2<=r2) 
	{
	  //norme=sqrt(norme2);
	  dot=r[0]*n1[0]+r[1]*n1[1]+r[2]*n1[2];
	  //if (dot<0) dot=-dot;
	  if (dot*dot>=cosalpha*cosalpha*norme2) 
	  {
		//fprintf(stdout, "normale=[%f\t%f\t%f], \tvect=[%f\t%f\t%f], \tangle=%f\n",
		//		n1[0],n1[1],n1[2],r[0],r[1],r[2], acos(dot/sqrt(norme2)));
		if(norme2>0){
		  //m1[0]=n1[0]-2*dot/norme2*r[0]; // m1 : symetrique de n1 par rapport 
		  //m1[1]=n1[1]-2*dot/norme2*r[1]; // au plan passant par le centre de
		  //m1[2]=n1[2]-2*dot/norme2*r[2]; // la sphere tangente et par le
										 // milieu du segment reliant les 2 
										 // points consideres
		  m1[0]=r[0];
		  m1[1]=r[1];
		  m1[2]=r[2];
		}
		else{
		  m1[0]=n1[0];
		  m1[1]=n1[1];
		  m1[2]=n1[2];
		}
		n=sqrt(m1[0]*m1[0]+m1[1]*m1[1]+m1[2]*m1[2]);
		m1[0]/=n;
		m1[1]/=n;
		m1[2]/=n;
		
		UnitVectorToSphericalAngles( m1, &alpha, &beta );
		theXX[z][y][x] = 1 - m1[0]*m1[0];
		theYY[z][y][x] = 1 - m1[1]*m1[1];
		theZZ[z][y][x] = 1 - m1[2]*m1[2];
		theXY[z][y][x] = - m1[0]*m1[1];
		theXZ[z][y][x] = - m1[0]*m1[2];
		theYZ[z][y][x] = - m1[1]*m1[2];

		theVP1[z][y][x]= 0;
		theVP2[z][y][x]= 1;
		theVP3[z][y][x]= 1;            

		theTht1[z][y][x]= alpha;            
		theTht2[z][y][x]= alpha;
		theTht3[z][y][x]= alpha+pi/(double)2;
		thePhi1[z][y][x]= beta;    
		thePhi2[z][y][x]= beta-pi/(double)2;
		thePhi3[z][y][x]= pi/(double)2;

		theZeros[z][y][x]= 0;            
		
		continue;
	  }
    }
    
    /* angle superieur a alpha ou rayon superieur a ray -> on ne vote pas */
    theXX[z][y][x] = 0;
    theYY[z][y][x] = 0;
    theZZ[z][y][x] = 0;
    theXY[z][y][x] = 0;
    theXZ[z][y][x] = 0;
    theYZ[z][y][x] = 0;
    theVP1[z][y][x]= 0;
    theVP2[z][y][x]= 0;
    theVP3[z][y][x]= 0;            
    theTht1[z][y][x]= angle[0]+pi/(double)2;
    theTht2[z][y][x]= angle[0];
    theTht3[z][y][x]= angle[0];            
    thePhi1[z][y][x]= pi/(double)2;
    thePhi2[z][y][x]= angle[1]-pi/(double)2;
    thePhi3[z][y][x]= angle[1];                  
    theZeros[z][y][x]= 1;                  
    
  } /* for(x,y,z) */
  
  return( 1 );
}

///////////////////




int MT_Compute3DStickField( vt_3Dtensor *sfield, int *dimFields, 
          double *angle, double *theCoeffs, enumTVmode mode)
{

  int x,y,z;

  int mid[3] = {dimFields[0]/2, dimFields[1]/2, dimFields[2]/2};
  double r[3];
  double n[3];
  double a[3];
  double rr[6]; // rr = r*transpose(r) matrice symetrique 3x3 = 6 elts
  double v[3];
  double sinteta;
  double fact=0;
  double alpha, beta;
  double s, kappa, l;
  double SUM = 0;
  
  double scale = theCoeffs[0];
  double zfact = theCoeffs[2]/theCoeffs[0];
  double pi = 3.14159265358979323846;
  double c = -16*log(0.1)*(scale-1)/pow(pi,2);
  
  float ***theXX = (float***)(sfield->imxx.array);


  float ***theYY = (float***)(sfield->imyy.array);
  float ***theZZ = (float***)(sfield->imzz.array);
  float ***theXY = (float***)(sfield->imxy.array);
  float ***theXZ = (float***)(sfield->imxz.array);
  float ***theYZ = (float***)(sfield->imyz.array);

  float ***theVP1 = (float***)(sfield->imvp1.array);
  float ***theVP2 = (float***)(sfield->imvp2.array);
  float ***theVP3 = (float***)(sfield->imvp3.array);


  float ***theTht1 = (float***)(sfield->imtheta1.array);
  float ***theTht2 = (float***)(sfield->imtheta2.array);
  float ***theTht3 = (float***)(sfield->imtheta3.array);

  float ***thePhi1 = (float***)(sfield->imphi1.array);
  float ***thePhi2 = (float***)(sfield->imphi2.array);
  float ***thePhi3 = (float***)(sfield->imphi3.array);
  
  unsigned char ***theZeros = (unsigned char***)(sfield->iszero.array);
  
  SphericalAnglesToUnitVector( angle[0], angle[1], n);
  
  for (z=0;z<dimFields[2];z++)
  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++) {

    r[0]=(double) x-mid[0];
    r[1]=(double) y-mid[1];
    r[2]=((double) z-mid[2]) / zfact;
    l = sqrt(pow(r[0],2) + pow(r[1],2) + pow(r[2],2));

    if (l == 0)
    {
      /*  Ici le point qui reçoit le vote est le point votant */
      theXX[z][y][x] = n[0]*n[0];
      theYY[z][y][x] = n[1]*n[1];
      theZZ[z][y][x] = n[2]*n[2];
      theXY[z][y][x] = n[0]*n[1];
      theXZ[z][y][x] = n[0]*n[2];
      theYZ[z][y][x] = n[1]*n[2];
      theVP1[z][y][x]= 0;
      theVP2[z][y][x]= 0;
      theVP3[z][y][x]= 1;            
      theTht1[z][y][x]= angle[0]+pi/(double)2;
      theTht2[z][y][x]= angle[0];
      theTht3[z][y][x]= angle[0];            
      thePhi1[z][y][x]= pi/(double)2;
      thePhi2[z][y][x]= angle[1]-pi/(double)2;
      thePhi3[z][y][x]= angle[1];    
      theZeros[z][y][x]= 0;            
      SUM += 1;
      continue;
    }

    r[0]=r[0]/l;
    r[1]=r[1]/l;
    r[2]=r[2]/l;

    
    if (mode == TVCLASSIC && fabs(r[0]*n[0]+r[1]*n[1]+r[2]*n[2])>sqrt(2)/2)
    {
      /* angle inferieur a pi/4 entre la normale en le point votant et 
         le vecteur reliant les deux points -> on ne vote pas */
      theXX[z][y][x] = 0;
      theYY[z][y][x] = 0;
      theZZ[z][y][x] = 0;
      theXY[z][y][x] = 0;
      theXZ[z][y][x] = 0;
      theYZ[z][y][x] = 0;
      theVP1[z][y][x]= 0;
      theVP2[z][y][x]= 0;
      theVP3[z][y][x]= 0;            
      theTht1[z][y][x]= angle[0]+pi/(double)2;
      theTht2[z][y][x]= angle[0];
      theTht3[z][y][x]= angle[0];            
      thePhi1[z][y][x]= pi/(double)2;
      thePhi2[z][y][x]= angle[1]-pi/(double)2;
      thePhi3[z][y][x]= angle[1];                  
      theZeros[z][y][x]= 1;                  
      continue;
    }      

    /* Calcul du vecteur v normal a l'arc de cercle reliant le point
       votant et le point recevant le vote */
    rr[0]=r[0]*r[0]; //rxx
    rr[1]=r[1]*r[1]; //ryy
    rr[2]=r[2]*r[2]; //rzz
    rr[3]=r[0]*r[1]; //rxy
    rr[4]=r[0]*r[2]; //rxz
    rr[5]=r[1]*r[2]; //ryz
    
    a[0] = rr[0]*n[0]+rr[3]*n[1]+rr[4]*n[2];
    a[1] = rr[3]*n[0]+rr[1]*n[1]+rr[5]*n[2];
    a[2] = rr[4]*n[0]+rr[5]*n[1]+rr[2]*n[2];

    v[0] = n[0]-2*a[0];
    v[1] = n[1]-2*a[1]; 
    v[2] = n[2]-2*a[2];       

    /* Calcul du facteur de décroissance en le pixel [x,y,z] */
    sinteta = r[0]*n[0]+r[1]*n[1]+r[2]*n[2];

    switch (mode) {
    default:
      fprintf(stderr, "MT_Compute3DStickField: unexpected value of input enumTVmode");
      return(0);
    case TVCLASSIC :
      s = (sinteta == 0) ? l : l*asin(sinteta)/sinteta;
      kappa = 2*sinteta/l;
      fact = exp(-(pow(s,2)+c*pow(kappa,2))/pow(scale,2));
      break;
    case CFTV :
      fact = exp(-pow(l/scale,2))*(1-pow(sinteta,2)); 
      break;
    }
    
    if (fact < 0.01) fact = 0;
    
    UnitVectorToSphericalAngles( v, &alpha, &beta );

    theXX[z][y][x] = v[0]*v[0]*fact;
    theYY[z][y][x] = v[1]*v[1]*fact;
    theZZ[z][y][x] = v[2]*v[2]*fact;
    theXY[z][y][x] = v[0]*v[1]*fact;
    theXZ[z][y][x] = v[0]*v[2]*fact;
    theYZ[z][y][x] = v[1]*v[2]*fact;

    theVP1[z][y][x]= 0;
    theVP2[z][y][x]= 0;
    theVP3[z][y][x]= fact;            

    theTht1[z][y][x]= alpha+pi/(double)2;
    theTht2[z][y][x]= alpha;
    theTht3[z][y][x]= alpha;            
    thePhi1[z][y][x]= pi/(double)2;
    thePhi2[z][y][x]= beta-pi/(double)2;
    thePhi3[z][y][x]= beta;    

    theZeros[z][y][x]= (fact>0) ? 0 : 1;            
    
    SUM += fact;
  } /* for(x,y,z) */
  for (z=0;z<dimFields[2];z++)
  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++)
  {
    theVP3[z][y][x] = theVP3[z][y][x]/(float)SUM;
    theXX[z][y][x] = theXX[z][y][x]/(float)SUM;
    theYY[z][y][x] = theYY[z][y][x]/(float)SUM;
    theZZ[z][y][x] = theZZ[z][y][x]/(float)SUM;
    theXY[z][y][x] = theXY[z][y][x]/(float)SUM;
    theXZ[z][y][x] = theXZ[z][y][x]/(float)SUM;
    theYZ[z][y][x] = theYZ[z][y][x]/(float)SUM;
  }
  
  return( 1 );
}



int MT_Compute3DBallFieldFromStick( vt_3Dtensor *bfield, vt_3Dtensor *sfields, 
          int *dimFields, int Nangles)
{
  char *proc = "MT_Compute3DBallFieldFromStick";
  int n, x, y, z;
  double hessien[9];
  double valprop[3];
  double vecprop[9];
  double v[3];
  double SUM = 0;
  double theta, phi;

  float ***sXX;
  float ***sYY;
  float ***sZZ;
  float ***sXY;
  float ***sXZ;
  float ***sYZ;

  float ***theXX;
  float ***theYY;
  float ***theZZ;
  float ***theXY;
  float ***theXZ;
  float ***theYZ;

  float ***theVP1;
  float ***theVP2;
  float ***theVP3;


  float ***theTht1;
  float ***theTht2;
  float ***theTht3;

  float ***thePhi1;
  float ***thePhi2;
  float ***thePhi3;

  unsigned char ***theZeros;

  char name[256];
  int t = (int) FLOAT;

  sprintf( name, "%s", "ballfield" );
  if ( VT_Alloc3Dtensor( bfield, name, dimFields[0], dimFields[1], dimFields[2],
      t ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate 3D tensor images\n", proc );
    return( -1 );

  }

  theXX = (float***)(bfield->imxx.array);
  theYY = (float***)(bfield->imyy.array);
  theZZ = (float***)(bfield->imzz.array);
  theXY = (float***)(bfield->imxy.array);
  theXZ = (float***)(bfield->imxz.array);
  theYZ = (float***)(bfield->imyz.array);

  theVP1 = (float***)(bfield->imvp1.array);
  theVP2 = (float***)(bfield->imvp2.array);
  theVP3 = (float***)(bfield->imvp3.array);


  theTht1 = (float***)(bfield->imtheta1.array);
  theTht2 = (float***)(bfield->imtheta2.array);
  theTht3 = (float***)(bfield->imtheta3.array);

  thePhi1 = (float***)(bfield->imphi1.array);
  thePhi2 = (float***)(bfield->imphi2.array);
  thePhi3 = (float***)(bfield->imphi3.array);
  
  theZeros = (unsigned char***)(bfield->iszero.array);
  
  for (z=0;z<dimFields[2];z++)
  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++) {
    theXX[z][y][x] = 0;
    theYY[z][y][x] = 0;
    theZZ[z][y][x] = 0;
    theXY[z][y][x] = 0;
    theXZ[z][y][x] = 0;
    theYZ[z][y][x] = 0;
    theZeros[z][y][x]=1;
  }  

  for (n=0;n<Nangles;n++)
  {
    vt_3Dtensor* sfield = sfields+n;
    sXX = (float ***)sfield->imxx.array;
    sYY = (float ***)sfield->imyy.array;
    sZZ = (float ***)sfield->imzz.array;
    sXY = (float ***)sfield->imxy.array;
    sYZ = (float ***)sfield->imyz.array;
    sXZ = (float ***)sfield->imxz.array;

    for (z=0;z<dimFields[2];z++)
    for (y=0;y<dimFields[1];y++)
    for (x=0;x<dimFields[0];x++) {
      
      theXX[z][y][x] = theXX[z][y][x] + sXX[z][y][x];
      theYY[z][y][x] = theYY[z][y][x] + sYY[z][y][x];
      theZZ[z][y][x] = theZZ[z][y][x] + sZZ[z][y][x];
      theXY[z][y][x] = theXY[z][y][x] + sXY[z][y][x];
      theXZ[z][y][x] = theXZ[z][y][x] + sXZ[z][y][x];
      theYZ[z][y][x] = theYZ[z][y][x] + sYZ[z][y][x];
    }  
  }
    
  for (z=0;z<dimFields[2];z++)
  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++) {
    hessien[0] = theXX[z][y][x];
    hessien[1] = hessien[3] = theXY[z][y][x];
    hessien[2] = hessien[6] = theXZ[z][y][x];
    hessien[4] = theYY[z][y][x];
    hessien[5] = hessien[7] = theYZ[z][y][x];
    hessien[8] = theZZ[z][y][x];

    if ( _ComputeEigensOfSymetricSquareMatrix( hessien, valprop, vecprop, 3 )
        != 1 ) {
      if (_verbose_)
        fprintf(stderr, "%s: error in computing eigens\n", proc);
      VT_Free3Dtensor(bfield);
      return( -1 );
    }
    if ( _SortEigensInAbsIncreasingOrder( valprop, vecprop, 3 ) != 1 ) {
    
      if (_verbose_)
        fprintf(stderr, "%s: error in sorting eigens\n", proc);
      VT_Free3Dtensor(bfield);
      return( -1 );
    }

    /* Les valeurs et vecteurs propres sont tries par fabs croissante */
    theVP1[z][y][x]=valprop[0];
    v[0]=vecprop[0]; v[1]=vecprop[3]; v[2]=vecprop[6];
    UnitVectorToSphericalAngles( v, &theta, &phi );
    theTht1[z][y][x]=theta;
    thePhi1[z][y][x]=phi;

    theVP2[z][y][x]=valprop[1];
    v[0]=vecprop[1]; v[1]=vecprop[4]; v[2]=vecprop[7];
    UnitVectorToSphericalAngles( v, &theta, &phi );
    theTht2[z][y][x]=theta;
    thePhi2[z][y][x]=phi;

    theVP3[z][y][x]=valprop[2];
    v[0]=vecprop[2]; v[1]=vecprop[5]; v[2]=vecprop[8];
    UnitVectorToSphericalAngles( v, &theta, &phi );
    theTht3[z][y][x]=theta;
    thePhi3[z][y][x]=phi;

    if (valprop[2]>0)
    {
      theZeros[z][y][x] = 0;
      SUM += valprop[2];
    }
  }  

  for (z=0;z<dimFields[2];z++)
  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++)
  {
    theVP1[z][y][x] = theVP1[z][y][x]/(float)SUM;
    theVP2[z][y][x] = theVP2[z][y][x]/(float)SUM;
    theVP3[z][y][x] = theVP3[z][y][x]/(float)SUM;
    theXX[z][y][x] = theXX[z][y][x]/(float)SUM;
    theYY[z][y][x] = theYY[z][y][x]/(float)SUM;
    theZZ[z][y][x] = theZZ[z][y][x]/(float)SUM;
    theXY[z][y][x] = theXY[z][y][x]/(float)SUM;
    theXZ[z][y][x] = theXZ[z][y][x]/(float)SUM;
    theYZ[z][y][x] = theYZ[z][y][x]/(float)SUM;

  }
  return( 1 );
}




int MT_Compute3DPlateFields( vt_3Dtensor **pfields, int *dimFields,
    double *theCoeffs, double **angles, int Nangles, int Nsticks,
    enumTVmode mode)
{

  /* On alloue *pfields,
   * On crée et alloue un array de stickfields de longueur Nsticks,
   * Pour chaque angle,
   * - on alloue le
   * - calcul de Nstick angles discrétisant le demi-cercle du plan normal à
   *   l'angle
   * - calcul des stick fields correspondants à chacun de ces Nstick angles
   * - on somme ces Nstick stickfields et on obtient le platefield correspondant
   * Fin Pour,
   * On libère les stickfields,
   * On retourne 1.
   */

  int n,m;
  char *proc = "MT_Compute3DPlateFields";
  double angle[2];
  char name[256];
  int t = (int) FLOAT;
  vt_3Dtensor *field;

  double **anglesPlane;
  vt_3Dtensor *sfields;

  /* Allocation des champs de tenseurs plate */
  (*pfields) =  malloc(Nangles*sizeof(vt_3Dtensor));

  if((*pfields)==NULL)
  {
    if (_verbose_)
      fprintf(stderr, "%s: allocation de (*pfield) echouee", proc);
    return(-1);
  }


  for (n=0; n<Nangles; n++) {
    sprintf( name, "%s.%d", "platefield", n );
    field = (*pfields+n);
    if ( VT_Alloc3Dtensor(field, name, dimFields[0], dimFields[1], dimFields[2],
        t ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate 3D tensor images\n", proc );
      for (m=0;m<n;m++){
        field = (*pfields+m);
        VT_Free3Dtensor( field );
      }
      free((*pfields));
      (*pfields)=NULL;
      return( -1 );

    }
  }

  /* Allocation des angles des champs de tenseurs stick normaux a
     l'angle du tenseur plate */
  anglesPlane = malloc(Nsticks*sizeof(double*));
  if (anglesPlane == NULL)
  {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate table of angles of plane \
            normal to vector\n", proc );
      for (n=0;n<Nangles;n++){
        field = (*pfields+n);
        VT_Free3Dtensor( field );
      }
      free((*pfields));
      (*pfields)=NULL;
      return( -1 );
  }
  for (n=0;n<Nsticks;n++)
  {
     anglesPlane[n] = malloc(2*sizeof(double));
     if (anglesPlane[n]==NULL)
     {
       if ( _verbose_ )
         fprintf( stderr, "%s: unable to allocate angle #%d of table of angles \
             of plane normal to vector\n", proc, n);
       for (n=0;n<Nangles;n++){
         field = (*pfields+n);
         VT_Free3Dtensor( field );
       }
       free((*pfields));
       (*pfields)=NULL;
       for (m=0;m<n;m++){
         free(anglesPlane[m]);
         anglesPlane[m]=NULL;
       }
       free(anglesPlane);
       anglesPlane=NULL;
       return( -1 );
     }
  }


  /* Pour chaque angle de la demi-sphere : */
  for (n=0; n<Nangles; n++) {
    angle[0]=angles[n][0];
    angle[1]=angles[n][1];
    field = (*pfields+n);


    // Calcul des Nangles champs de tenseurs plate
    if (MT_ComputeNormalAngles(&anglesPlane, Nsticks, angle[0], angle[1]) !=1)
    {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to compute normal angles of angle #%d = \
              (%f ; %f)\n", proc, n, angle[0], angle[1]);
        for (m=0;m<Nangles;m++){
          field = (*pfields+m);
          VT_Free3Dtensor( field );
        }
        free((*pfields));
        (*pfields)=NULL;
        for (m=0;m<Nsticks;m++){
          free(anglesPlane[m]);
          anglesPlane[m]=NULL;
        }
        free(anglesPlane);
        anglesPlane=NULL;

        return( -1 );
    }


    // Calcul des stick tensor fields de directions normales a la direction
    // du plate tensor
    if(MT_Compute3DStickFields( &sfields, dimFields, anglesPlane, theCoeffs,
        Nsticks, mode) !=1)
    {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to compute stick fields of normal \
              angles of angle #%d = (%f ; %f)\n", proc, n, angle[0], angle[1]);
        for (m=0;m<Nangles;m++){
          field = (*pfields+m);
          VT_Free3Dtensor( field );
        }
        free((*pfields));
        (*pfields)=NULL;
        for (m=0;m<Nsticks;m++){
          free(anglesPlane[m]);
          anglesPlane[m]=NULL;
        }
        free(anglesPlane);
        anglesPlane=NULL;

        return( -1 );
    }


    // Calcul du plate tensor field a partir des stick tensor fields
    if(MT_Compute3DPlateFieldFromSticks(field, sfields, dimFields, Nsticks, n)
        !=1)
    {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to compute plate field of angle #%d = \
              (%f ; %f)\n", proc, n, angle[0], angle[1]);
        for (m=0;m<Nangles;m++){
          field = (*pfields+m);
          VT_Free3Dtensor( field );
        }
        free((*pfields));
        (*pfields)=NULL;
        for (m=0;m<Nsticks;m++){
          free(anglesPlane[m]);
          anglesPlane[m]=NULL;
          VT_Free3Dtensor( sfields+m);
        }
        free(sfields);
        sfields=NULL;
        free(anglesPlane);
        anglesPlane=NULL;
        return( -1 );
    }

    // Liberation des stick fields (alloues dans MT_Compute3DStickFields)
    for (m=0;m<Nsticks;m++){
      VT_Free3Dtensor( sfields+m);
    }
    free(sfields);
    sfields=(vt_3Dtensor *)NULL;

  }

  // Liberation de la memoire allouee aux anglesPlanes
  for (m=0;m<Nsticks;m++){
    free(anglesPlane[m]);
    anglesPlane[m]=NULL;
  }
  free(anglesPlane);
  anglesPlane=NULL;

  return( 1 );
}


int MT_Compute3DPlateFieldFromSticks( vt_3Dtensor *pfield,
    vt_3Dtensor *sfields, int *dimFields, int Nsticks, int N)
{
  

  // Faire la somme des sticks, tout simplement (cf Compute3DBallFieldFromStick)

  char *proc = "MT_Compute3DPlateFieldFromSticks";
  int n, x, y, z;
  double hessien[9];
  double valprop[3];
  double vecprop[9];
  double v[3];

  double theta, phi;
  double SUM=0;

  float ***sXX;
  float ***sYY;
  float ***sZZ;
  float ***sXY;
  float ***sXZ;
  float ***sYZ;

  float ***theXX;
  float ***theYY;
  float ***theZZ;
  float ***theXY;
  float ***theXZ;
  float ***theYZ;

  float ***theVP1;
  float ***theVP2;
  float ***theVP3;


  float ***theTht1;
  float ***theTht2;
  float ***theTht3;

  float ***thePhi1;
  float ***thePhi2;
  float ***thePhi3;

  unsigned char ***theZeros;



  theXX = (float***)(pfield->imxx.array);
  theYY = (float***)(pfield->imyy.array);
  theZZ = (float***)(pfield->imzz.array);
  theXY = (float***)(pfield->imxy.array);
  theXZ = (float***)(pfield->imxz.array);
  theYZ = (float***)(pfield->imyz.array);

  theVP1 = (float***)(pfield->imvp1.array);
  theVP2 = (float***)(pfield->imvp2.array);
  theVP3 = (float***)(pfield->imvp3.array);


  theTht1 = (float***)(pfield->imtheta1.array);
  theTht2 = (float***)(pfield->imtheta2.array);
  theTht3 = (float***)(pfield->imtheta3.array);

  thePhi1 = (float***)(pfield->imphi1.array);
  thePhi2 = (float***)(pfield->imphi2.array);
  thePhi3 = (float***)(pfield->imphi3.array);

  theZeros = (unsigned char***)(pfield->iszero.array);

  for (z=0;z<dimFields[2];z++)
  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++) {
    theXX[z][y][x] = 0;
    theYY[z][y][x] = 0;
    theZZ[z][y][x] = 0;
    theXY[z][y][x] = 0;
    theXZ[z][y][x] = 0;
    theYZ[z][y][x] = 0;
    theZeros[z][y][x]=1;
  }

  for (n=0;n<Nsticks;n++)
  {
    vt_3Dtensor* sfield = sfields+n;
    sXX = (float ***)sfield->imxx.array;
    sYY = (float ***)sfield->imyy.array;
    sZZ = (float ***)sfield->imzz.array;
    sXY = (float ***)sfield->imxy.array;
    sYZ = (float ***)sfield->imyz.array;
    sXZ = (float ***)sfield->imxz.array;

    for (z=0;z<dimFields[2];z++)
    for (y=0;y<dimFields[1];y++)
    for (x=0;x<dimFields[0];x++) {

      theXX[z][y][x] = theXX[z][y][x] + sXX[z][y][x];
      theYY[z][y][x] = theYY[z][y][x] + sYY[z][y][x];
      theZZ[z][y][x] = theZZ[z][y][x] + sZZ[z][y][x];
      theXY[z][y][x] = theXY[z][y][x] + sXY[z][y][x];
      theXZ[z][y][x] = theXZ[z][y][x] + sXZ[z][y][x];
      theYZ[z][y][x] = theYZ[z][y][x] + sYZ[z][y][x];
    }
  }

  for (z=0;z<dimFields[2];z++)
  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++) {
    hessien[0] = theXX[z][y][x];
    hessien[1] = hessien[3] = theXY[z][y][x];
    hessien[2] = hessien[6] = theXZ[z][y][x];
    hessien[4] = theYY[z][y][x];
    hessien[5] = hessien[7] = theYZ[z][y][x];
    hessien[8] = theZZ[z][y][x];

    if ( _ComputeEigensOfSymetricSquareMatrix( hessien, valprop, vecprop, 3 )
        != 1 ) {
      if (_verbose_)
        fprintf(stderr, "%s: error in computing eigens\n", proc);
      VT_Free3Dtensor(pfield);
      return( -1 );
    }
    if ( _SortEigensInAbsIncreasingOrder( valprop, vecprop, 3 ) != 1 ) {

      if (_verbose_)
        fprintf(stderr, "%s: error in sorting eigens\n", proc);
      VT_Free3Dtensor(pfield);
      return( -1 );
    }

    /* Les valeurs et vecteurs propres sont tries par fabs croissante */
    theVP1[z][y][x]=valprop[0];
    v[0]=vecprop[0]; v[1]=vecprop[3]; v[2]=vecprop[6];
    UnitVectorToSphericalAngles( v, &theta, &phi );
    theTht1[z][y][x]=theta;
    thePhi1[z][y][x]=phi;

    theVP2[z][y][x]=valprop[1];
    v[0]=vecprop[1]; v[1]=vecprop[4]; v[2]=vecprop[7];
    UnitVectorToSphericalAngles( v, &theta, &phi );
    theTht2[z][y][x]=theta;
    thePhi2[z][y][x]=phi;

    theVP3[z][y][x]=valprop[2];
    v[0]=vecprop[2]; v[1]=vecprop[5]; v[2]=vecprop[8];
    UnitVectorToSphericalAngles( v, &theta, &phi );
    theTht3[z][y][x]=theta;
    thePhi3[z][y][x]=phi;

    if (valprop[2]>0)
    {
      theZeros[z][y][x] = 0;
      SUM += valprop[2];
    }
  }

  for (z=0;z<dimFields[2];z++)
  for (y=0;y<dimFields[1];y++)
  for (x=0;x<dimFields[0];x++)
  {
    theVP1[z][y][x] = theVP1[z][y][x]/(float)SUM;
    theVP2[z][y][x] = theVP2[z][y][x]/(float)SUM;
    theVP3[z][y][x] = theVP3[z][y][x]/(float)SUM;
    theXX[z][y][x] = theXX[z][y][x]/(float)SUM;
    theYY[z][y][x] = theYY[z][y][x]/(float)SUM;
    theZZ[z][y][x] = theZZ[z][y][x]/(float)SUM;
    theXY[z][y][x] = theXY[z][y][x]/(float)SUM;
    theXZ[z][y][x] = theXZ[z][y][x]/(float)SUM;
    theYZ[z][y][x] = theYZ[z][y][x]/(float)SUM;

  }


  return( 1 );
}



int MT_ComputeNormalAngles(double ***anglesPlane, int Nsticks,
    double theta, double phi)
{
  double v0[3], v1[3], v2[3], v[3];
  double angle, th, ph;
  double pi = 3.14159265358979323846;
  int i;

  SphericalAnglesToUnitsVectors( theta, phi, v0, v1, v2 );
  angle = pi/Nsticks;
  for (i = 0; i<Nsticks; i++)
  {
      v[0]=cos(i*angle)*v1[0]+sin(i*angle)*v2[0];
      v[1]=cos(i*angle)*v1[1]+sin(i*angle)*v2[1];
      v[2]=cos(i*angle)*v1[2]+sin(i*angle)*v2[2];
      UnitVectorToSphericalAngles( v, &th, &ph );
      (*anglesPlane)[i][0]=th;
      (*anglesPlane)[i][1]=ph;
  }


  return(1);
}

int MT_Compute3DAngles(double ***angles, int Niter)
{ 
  char *proc = "MT_Compute3DAngles";
  int i,j;
  int Nangles, indA;


  double theta, phi;
  double v[3];
  int nP, nT;
  triangle *T;
  point *P;  
  
  Nangles = MT_RecTriangulation(&P, &T, &nT, &nP, Niter);
  if( Nangles <0)
  {
    if (_verbose_)
      fprintf(stderr, "%s: error in computing recursive triangulation\n", proc);
    free(T);
    free(P);
    T=NULL;
    P=NULL;
    return(-1);
  }


  (*angles) = malloc(sizeof(double*)*Nangles);
  if((*angles)==NULL)
  {
    if (_verbose_)
      fprintf(stderr, "%s: error in allocating angles\n", proc);
    free(T);
    free(P);
    T=NULL;
    P=NULL;
    return(-1);
  }

  for (i=0;i<Nangles;i++)
  {
    (*angles)[i]= malloc(sizeof(double)*2);
    if((*angles)[i]==NULL) {
      if (_verbose_)
        fprintf(stderr, "%s: error in allocating angle #%d\n", proc, i);
      free(T);
      free(P);
      T=NULL;
      P=NULL;
      for(j=0;j<i;j++){
        free((*angles)[j]);
        (*angles)[j]=NULL;
      }
      free((*angles));
      (*angles)=NULL;
      return(-1);
    }
  }
  indA = 0;
  for (i=0; i<nP; i++)
  {
    v[0]=P[i].x;
    v[1]=P[i].y;
    v[2]=P[i].z;


    if (v[2]<-FLTZERO) continue;
    if (v[2]<FLTZERO && v[1]<-FLTZERO) continue;
    if (v[2]<FLTZERO && v[1]<FLTZERO && v[0]<-FLTZERO) continue;
    if (indA == Nangles)
    {
      if (_verbose_)
        fprintf(stderr, "%s: number of angles computed higher than expected \
            (more than %d)\n", proc, Nangles);
      free(T);
      free(P);
      T=NULL;
      P=NULL;
      for(j=0;j<Nangles;j++){
        free((*angles)[j]);
        (*angles)[j]=NULL;
      }
      free((*angles));
      (*angles)=NULL;
      return(-1);
    }
    UnitVectorToSphericalAngles( v, &theta, &phi );
    (*angles)[indA][0] = theta;
    (*angles)[indA][1] = phi;
    indA++;
  }  

  if (indA != Nangles)
  {
    if (_verbose_)
      fprintf(stderr, "%s: unexpected number of angles computed \
          (%d instead of %d)\n", proc, indA, Nangles);
    free(T);
    free(P);
    T=NULL;
    P=NULL;
    for(i=0;i<Nangles;i++){
      free((*angles)[i]);
      (*angles)[i]=NULL;
    }
    free((*angles));
    (*angles)=NULL;
    return(-1);
  }

  free(T);
  free(P);
  T=NULL;
  P=NULL;
  
  /* Eventuel qsorting pour que les angles soient bien tries...
     Non fonctionnel pr l'instant */
  //qsort(*angles, Nangles, sizeof(double*), cmpAngles);

  return( Nangles );

}


int addNewPointToAngles(double ***angles, point P, int ind)
{
  char *proc = "addNewPointToAngles";
  int i;
  double theta, phi;
  double dot;

  double v[3], v0[3];
  v0[0]=P.x; v0[1]=P.y; v0[2]=P.z;
  for (i=0; i<ind; i++)
  {
    SphericalAnglesToUnitVector( (*angles)[i][0], (*angles)[i][1], v);;
    if (fabs(fabs(v[0])-fabs(v0[0])) > FLTZERO)
      continue;
    if (fabs(fabs(v[1])-fabs(v0[1])) > FLTZERO)
      continue;
    if (fabs(fabs(v[2])-fabs(v0[2])) > FLTZERO)
      continue;
    dot = v[0]*v0[0]+v[1]*v0[1]+v[2]*v0[2];
    if (fabs(fabs(dot)-1)>FLTZERO)
      continue;
    return(0);
  }

  UnitVectorToSphericalAngles( v, &theta, &phi );

  double **temp = NULL;
  temp = realloc((*angles), sizeof(double*)*ind+1);
  if (temp==NULL)
  {
    if (_verbose_)
     fprintf(stderr,"%s: Reallocation impossible\n", proc);
   return(-1);
  }
  (*angles) = temp;

  (*angles)[ind]= malloc(sizeof(double)*2);
  if((*angles)[ind]==NULL) {
    if (_verbose_)
      fprintf(stderr, "%s: error in allocating angle[%d]\n", proc, ind);
    return(-1);
  }

  (*angles)[ind][0]=theta;
  (*angles)[ind][0]=phi;

  return (1);
}


void midPoint(point *mid, point a, point b)
{
  double norm;
  double v[3], va[3], vb[3];
  va[0]=a.x; va[1]=a.y; va[2]=a.z;
  vb[0]=b.x; vb[1]=b.y; vb[2]=b.z;

  v[0]=(va[0]+vb[0]);
  v[1]=(va[1]+vb[1]);
  v[2]=(va[2]+vb[2]);
  norm = pow(pow(v[0],2)+pow(v[1],2)+pow(v[2],2), 0.5);
  v[0]=v[0]/norm;
  v[1]=v[1]/norm;
  v[2]=v[2]/norm;
  mid->x = v[0];
  mid->y = v[1];
  mid->z = v[2];

  return;
}

void newTriangle(triangle *t, point a, point b, point c)
{
  t->a = a;
  t->b = b;
  t->c = c;
  return;
}

void setTriangle(triangle **T, triangle t, int n)
{
  (*T)[n] = t;
  return;
}

int addPoint(point **P, point p, int n, int n0)
{
  int i;
  for (i=n0;i<n;i++)
  {
    if (isequalPoint((*P)[i],p)==1)
      return(0);
  }
  (*P)[n]=p;

  return(1);
}

int isequalPoint(point a, point b)
{

  double diff;
  double xa = a.x, xb = b.x;
  diff = xa-xb;
  if (diff > FLTZERO || -diff > FLTZERO)
    return(0);
  double ya = a.y, yb = b.y;
  diff = ya-yb;
  if (diff > FLTZERO || -diff > FLTZERO)
    return(0);
  double za = a.z, zb = b.z;
  diff = za-zb;
  if (diff > FLTZERO || -diff > FLTZERO)
    return(0);

  return(1);
}



int MT_RecTriangulation(point **P, triangle **T, int *nT, int *nP, int niter)
{
  /* Points et triangles de la sphere unite */
  /* La fonction retourne le nombre de points calcules appartenant a la demi
     sphere du demi espace (z>0|z==0&y>0|z==0&y==0&x>0) */
  char *proc = "MT_RecTriangulation";
  double v[3];
  
  if (niter < 1){
    if (_verbose_)
      fprintf(stderr, "%s: bad input argument (niter=%d)",proc, niter);   
    return(-1);
  }    
  if (niter == 1)
  {

    *nT = 8;
    *nP = 6;


    (*T) = malloc((*nT)*sizeof(triangle));
    if ((*T) == NULL){
      if (_verbose_)
        fprintf(stderr, "%s: allocation failed",proc);   
      return(-1);
    }
    (*P) = malloc((*nP)*sizeof(point));
    if ((*P) == NULL)
    {
      if (_verbose_)
        fprintf(stderr, "%s: allocation failed",proc);   
      free((*T));
      (*T) = NULL;
      return(-1);
    }

    /* Iteration 1 : 6 points */
    v[0] = 1; v[1] = 0; v[2] = 0;

    (*P)[0].x = v[0];
    (*P)[0].y = v[1];
    (*P)[0].z = v[2];

    v[0] = 0; v[1] = 1; v[2] = 0;

    (*P)[1].x = v[0];
    (*P)[1].y = v[1];
    (*P)[1].z = v[2];


    v[0] = 0; v[1] = 0; v[2] = 1;

    (*P)[2].x = v[0];
    (*P)[2].y = v[1];
    (*P)[2].z = v[2];


    v[0] = -1; v[1] = 0; v[2] = 0;

    (*P)[3].x = v[0];
    (*P)[3].y = v[1];
    (*P)[3].z = v[2];

    v[0] = 0; v[1] = -1; v[2] = 0;

    (*P)[4].x = v[0];
    (*P)[4].y = v[1];
    (*P)[4].z = v[2];

    v[0] = 0; v[1] = 0; v[2] = -1;
    (*P)[5].x = v[0];
    (*P)[5].y = v[1];
    (*P)[5].z = v[2];

    /* Iteration 1 : 8 triangles */
    (*T)[0].a = (*P)[0];
    (*T)[0].b = (*P)[1];
    (*T)[0].c = (*P)[2];
  
    (*T)[1].a = (*P)[1];
    (*T)[1].b = (*P)[3];
    (*T)[1].c = (*P)[2];
  
    (*T)[2].a = (*P)[3];
    (*T)[2].b = (*P)[4];
    (*T)[2].c = (*P)[2];
  
    (*T)[3].a = (*P)[4];
    (*T)[3].b = (*P)[0];
    (*T)[3].c = (*P)[2];
  
    (*T)[4].a = (*P)[1];
    (*T)[4].b = (*P)[0];
    (*T)[4].c = (*P)[5];
  
    (*T)[5].a = (*P)[3];
    (*T)[5].b = (*P)[1];
    (*T)[5].c = (*P)[5];
  
    (*T)[6].a = (*P)[4];
    (*T)[6].b = (*P)[3];
    (*T)[6].c = (*P)[5];
  
    (*T)[7].a = (*P)[0];
    (*T)[7].b = (*P)[4];
    (*T)[7].c = (*P)[5];

    
    return(3);  // seuls (1,0,0), (0,1,0) et (0,0,1) sont dans la demi-sphere
  } /* niter == 1 */
 
  triangle *oldT;
  int oT;
  int oP;
  int nbT = *nT;
  int nbP = *nP;
  int i;
  int indT;
  int indP;
  int res = MT_RecTriangulation(P, &oldT, &oT, &oP, niter-1);
  
  if (res < 0)
  {
    if (_verbose_)
      fprintf(stderr, "%s: error in computing recursive triangulation #%d\n",
          proc, niter-1);
    free((*P));
    (*P)=NULL;
    free(oldT);
    oldT=NULL;
    free((*T));
    T=NULL;
    return(-1);     
  }
  


  nbT = (int)pow(2,2*niter+1);        //nombre de triangles a la nieme iteration
  nbP = (int)6+4*(pow(4, niter-1)-1); //nombre de points a la nieme iteration

  (*T) = malloc(sizeof(triangle)*nbT);//il faudra liberer oldT en fin de fct
  if((*T)==NULL)
  {
    if (_verbose_)
      fprintf(stderr,"%s: Allocation impossible\n", proc);
    free((*P));
    free(oldT);
    (*P)=NULL;
    oldT=NULL; 
  }  
  point *temp = NULL;
  temp = realloc((*P), sizeof(point)*nbP);  
  if (temp==NULL)
  {
    if (_verbose_)
     fprintf(stderr,"%s: Reallocation impossible\n", proc);
   free((*P));
   free(oldT);
   free((*T));
   (*P)=NULL;
   oldT=NULL;
   (*T)=NULL;
   return(-1);
  }
  
  (*P) = temp;
    
  indT = 0; // indice auquel on placera chaque nouveau triangle dans T
  indP = oP;// indice auquel on placera chaque nouveau point dans P
  
  
  
  for (i=0; i<oT; i++) {
    triangle t = oldT[i];
    point a = t.a, b=t.b, c=t.c;

    point mab, mac, mbc;
    triangle tn1, tn2, tn3, tn4;
    
    midPoint(&mab, a, b);
    midPoint(&mac, a, c);
    midPoint(&mbc, b, c);
    
    newTriangle(&tn1, a, mac, mab);
    newTriangle(&tn2, b, mbc, mab);
    newTriangle(&tn3, c, mac, mbc);
    newTriangle(&tn4, mbc, mac, mab);

        
    setTriangle(T, tn1, indT++);
    setTriangle(T, tn2, indT++);
    setTriangle(T, tn3, indT++);
    setTriangle(T, tn4, indT++);
    
    if(addPoint(P, mab, indP, oP) == 1)
    {
      if (mab.z>FLTZERO || (mab.z>=-FLTZERO && mab.y>FLTZERO) ||
          (mab.z>=-FLTZERO && mab.y>=-FLTZERO && mab.x>FLTZERO))
        res++;
      indP++;
    }

    if(addPoint(P, mac, indP, oP) == 1)
    {
      if (mac.z>FLTZERO || (mac.z>=-FLTZERO && mac.y>FLTZERO) ||
          (mac.z>=-FLTZERO && mac.y>=-FLTZERO && mac.x>FLTZERO))
        res++;
      indP++;
    }

    if(addPoint(P, mbc, indP, oP) == 1)
    {
      if (mbc.z>FLTZERO || (mbc.z>=-FLTZERO && mbc.y>FLTZERO) ||
          (mbc.z>=-FLTZERO && mbc.y>=-FLTZERO && mbc.x>FLTZERO))
        res++;
      indP++;
    }

  }

  
  if (indP != nbP)
  {
    fprintf(stderr, "%s: nb de sommets inattendu :\nnbP=%d, indP=%d, op=%d\n",
        proc,nbP,indP,oP);
    if (indT != nbT)
      fprintf(stderr,"%s: nb de facettes inattendu :\nnbT=%d, indT=%d\n",
          proc,nbT,indT);
    free(oldT);
    oldT=NULL;
    return(-1);
  }
  if (indT != nbT)
  {
    fprintf(stderr,"%s: nb de facettes inattendu :\nnbT=%d, indT=%d\n",
        proc,nbT,indT);
    free(oldT);
    oldT=NULL;
    return(-1);
  }

  *nT = nbT;
  *nP = nbP;

  free(oldT);
  oldT=NULL;
  
  return(res);
}

int cmpAngles(double *angle1, double *angle2)// pour le qsort -> non fonctionnel
{
  if (angle1[0]-angle2[0]<-FLTZERO)
    return(-1);
  if (angle1[0]-angle2[0]<FLTZERO && angle1[1]-angle2[1]<-FLTZERO)
    return(-1);
  if (angle1[0]-angle2[0]<FLTZERO && angle1[1]-angle2[1]<FLTZERO)
    return(0);
  return(1);
}



////////////////// ALLOCATIONS ////////////////////////////////////////////////


int VT_Alloc3DtensorFromImage ( vt_3Dtensor *par, vt_image *im,
    char *genericname )
{
  char name[256];

  sprintf( name, "%s.imxx.inr", genericname );
  VT_InitFromImage( &(par->imxx), im, name,  (int)FLOAT );
  sprintf( name, "%s.imxy.inr", genericname );
  VT_InitFromImage( &(par->imxy), im, name,  (int)FLOAT );
  sprintf( name, "%s.imyy.inr", genericname );
  VT_InitFromImage( &(par->imyy), im, name,  (int)FLOAT );
  sprintf( name, "%s.imzz.inr", genericname );
  VT_InitFromImage( &(par->imzz), im, name,  (int)FLOAT );
  sprintf( name, "%s.imxz.inr", genericname );
  VT_InitFromImage( &(par->imxz), im, name,  (int)FLOAT );
  sprintf( name, "%s.imyz.inr", genericname );
  VT_InitFromImage( &(par->imyz), im, name,  (int)FLOAT );



  sprintf( name, "%s.imvp1.inr", genericname );
  VT_InitFromImage( &(par->imvp1),  im, name,  (int)FLOAT );
  sprintf( name, "%s.imvp2.inr", genericname );
  VT_InitFromImage( &(par->imvp2),  im, name,  (int)FLOAT );
  sprintf( name, "%s.imvp3.inr", genericname );
  VT_InitFromImage( &(par->imvp3),  im, name,  (int)FLOAT );

  sprintf( name, "%s.imtheta1.inr", genericname );
  VT_InitFromImage( &(par->imtheta1),  im, name,  (int)FLOAT );
  sprintf( name, "%s.imtheta2.inr", genericname );
  VT_InitFromImage( &(par->imtheta2),  im, name,  (int)FLOAT );
  sprintf( name, "%s.imtheta3.inr", genericname );
  VT_InitFromImage( &(par->imtheta3),  im, name,  (int)FLOAT );
  sprintf( name, "%s.imphi1.inr", genericname );
  VT_InitFromImage( &(par->imphi1),  im, name,  (int)FLOAT );
  sprintf( name, "%s.imphi2.inr", genericname );
  VT_InitFromImage( &(par->imphi2),  im, name,  (int)FLOAT );
  sprintf( name, "%s.imphi3.inr", genericname );
  VT_InitFromImage( &(par->imphi3),  im, name,  (int)FLOAT );
  sprintf( name, "%s.iszero.inr", genericname );
  VT_InitFromImage( &(par->iszero),  im, name,  (int)UCHAR );

  if ( VT_AllocImage( &(par->imxx) ) != 1 ) {
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imxy) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imyy) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imzz) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imxz) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imyz) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imvp1) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imvp2) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imvp3) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imtheta1) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imtheta2) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imtheta3) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imphi1) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    VT_FreeImage( &(par->imtheta3) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imphi2) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    VT_FreeImage( &(par->imtheta3) );
    VT_FreeImage( &(par->imphi1) );        
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imphi3) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    VT_FreeImage( &(par->imtheta3) );
    VT_FreeImage( &(par->imphi1) );
    VT_FreeImage( &(par->imphi2) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->iszero) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    VT_FreeImage( &(par->imtheta3) );
    VT_FreeImage( &(par->imphi1) );
    VT_FreeImage( &(par->imphi2) );
    VT_FreeImage( &(par->imphi3) );
    return( -1 );
  }


  return( 1 );
} 



int VT_Alloc3Dtensor ( vt_3Dtensor *par /* image to be initialized */,
 		  char *genericname /* image generic name */,
		  int dimx /* X dimension */,
		  int dimy /* Y dimension */,
		  int dimz /* Z dimension */,
		  int type /* image type  */ )
{
  char name[256];

  sprintf( name, "%s.imxx.inr", genericname );
  VT_InitVImage( &(par->imxx), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imxy.inr", genericname );
  VT_InitVImage( &(par->imxy), name, (int)1, dimx, dimy, dimz, type );  
  sprintf( name, "%s.imyy.inr", genericname );
  VT_InitVImage( &(par->imyy), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imzz.inr", genericname );
  VT_InitVImage( &(par->imzz), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imxz.inr", genericname );
  VT_InitVImage( &(par->imxz), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imyz.inr", genericname );
  VT_InitVImage( &(par->imyz), name, (int)1, dimx, dimy, dimz, type );



  sprintf( name, "%s.imvp1.inr", genericname );
  VT_InitVImage( &(par->imvp1), name, (int)1, dimx, dimy, dimz, type );  
  sprintf( name, "%s.imvp2.inr", genericname );
  VT_InitVImage( &(par->imvp2), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imvp3.inr", genericname );
  VT_InitVImage( &(par->imvp3), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imtheta1.inr", genericname );
  VT_InitVImage( &(par->imtheta1), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imtheta2.inr", genericname );
  VT_InitVImage( &(par->imtheta2), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imtheta3.inr", genericname );
  VT_InitVImage( &(par->imtheta3), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imphi1.inr", genericname );
  VT_InitVImage( &(par->imphi1), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imphi2.inr", genericname );
  VT_InitVImage( &(par->imphi2), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.imphi3.inr", genericname );
  VT_InitVImage( &(par->imphi3), name, (int)1, dimx, dimy, dimz, type );
  sprintf( name, "%s.iszero.inr", genericname );
  VT_InitVImage( &(par->iszero), name, (int)1, dimx, dimy, dimz, (int)UCHAR );

  if ( VT_AllocImage( &(par->imxx) ) != 1 ) {
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imxy) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imyy) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imzz) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imxz) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imyz) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imvp1) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    return( -1 );
  }
  if ( VT_AllocImage( &(par->imvp2) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imvp3) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imtheta1) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imtheta2) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imtheta3) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imphi1) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    VT_FreeImage( &(par->imtheta3) );
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imphi2) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    VT_FreeImage( &(par->imtheta3) );
    VT_FreeImage( &(par->imphi1) );        
    return( -1 );
  }

  if ( VT_AllocImage( &(par->imphi3) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    VT_FreeImage( &(par->imtheta3) );
    VT_FreeImage( &(par->imphi1) );
    VT_FreeImage( &(par->imphi2) );
    return( -1 );
  }
  
  if ( VT_AllocImage( &(par->iszero) ) != 1 ) {
    VT_FreeImage( &(par->imxx) );
    VT_FreeImage( &(par->imxy) );
    VT_FreeImage( &(par->imyy) );
    VT_FreeImage( &(par->imzz) );
    VT_FreeImage( &(par->imxz) );
    VT_FreeImage( &(par->imyz) );
    VT_FreeImage( &(par->imvp1) );
    VT_FreeImage( &(par->imvp2) );
    VT_FreeImage( &(par->imvp3) );
    VT_FreeImage( &(par->imtheta1) );
    VT_FreeImage( &(par->imtheta2) );
    VT_FreeImage( &(par->imtheta3) );
    VT_FreeImage( &(par->imphi1) );
    VT_FreeImage( &(par->imphi2) );
    VT_FreeImage( &(par->imphi3) );
    return( -1 );
  }

  return( 1 );
} 



void VT_Free3Dtensor( vt_3Dtensor *par )
{
  VT_FreeImage( &(par->imxx) );
  VT_FreeImage( &(par->imxy) );
  VT_FreeImage( &(par->imyy) );
  VT_FreeImage( &(par->imzz) );
  VT_FreeImage( &(par->imxz) );
  VT_FreeImage( &(par->imyz) );

  VT_FreeImage( &(par->imvp1) );
  VT_FreeImage( &(par->imvp2) );
  VT_FreeImage( &(par->imvp3) );
  VT_FreeImage( &(par->imtheta1) );
  VT_FreeImage( &(par->imtheta2) );
  VT_FreeImage( &(par->imtheta3) );
  VT_FreeImage( &(par->imphi1) );
  VT_FreeImage( &(par->imphi2) );
  VT_FreeImage( &(par->imphi3) );
  VT_FreeImage( &(par->iszero) );
}



void VT_Write3Dtensor( vt_3Dtensor *par )
{
  VT_WriteInrimage( &(par->imxx) );
  VT_WriteInrimage( &(par->imxy) );
  VT_WriteInrimage( &(par->imyy) );
  VT_WriteInrimage( &(par->imzz) );
  VT_WriteInrimage( &(par->imxz) );
  VT_WriteInrimage( &(par->imyz) );

  VT_WriteInrimage( &(par->imvp1) );
  VT_WriteInrimage( &(par->imvp2) );
  VT_WriteInrimage( &(par->imvp3) );
  VT_WriteInrimage( &(par->imtheta1) );
  VT_WriteInrimage( &(par->imtheta2) );
  VT_WriteInrimage( &(par->imtheta3) );
  VT_WriteInrimage( &(par->imphi1) );
  VT_WriteInrimage( &(par->imphi2) );
  VT_WriteInrimage( &(par->imphi3) );
  VT_WriteInrimage( &(par->iszero) );
}

void VT_Write3DtensorWithName( vt_3Dtensor *par, char *genericname)
{
  char name[256];
  sprintf( name, "%s.imxx.inr", genericname );
  VT_WriteInrimageWithName( &(par->imxx), name );
  sprintf( name, "%s.imxy.inr", genericname );
  VT_WriteInrimageWithName( &(par->imxy), name );
  sprintf( name, "%s.imyy.inr", genericname );
  VT_WriteInrimageWithName( &(par->imyy), name );
  sprintf( name, "%s.imzz.inr", genericname );
  VT_WriteInrimageWithName( &(par->imzz), name );
  sprintf( name, "%s.imxz.inr", genericname );
  VT_WriteInrimageWithName( &(par->imxz), name );
  sprintf( name, "%s.imyz.inr", genericname );
  VT_WriteInrimageWithName( &(par->imyz), name );

  sprintf( name, "%s.imvp1.inr", genericname );
  VT_WriteInrimageWithName( &(par->imvp1), name );
  sprintf( name, "%s.imvp2.inr", genericname );
  VT_WriteInrimageWithName( &(par->imvp2), name );
  sprintf( name, "%s.imvp3.inr", genericname );
  VT_WriteInrimageWithName( &(par->imvp3), name );
  sprintf( name, "%s.imtheta1.inr", genericname );
  VT_WriteInrimageWithName( &(par->imtheta1), name );
  sprintf( name, "%s.imtheta2.inr", genericname );
  VT_WriteInrimageWithName( &(par->imtheta2), name );
  sprintf( name, "%s.imtheta3.inr", genericname );
  VT_WriteInrimageWithName( &(par->imtheta3), name );
  sprintf( name, "%s.imphi1.inr", genericname );
  VT_WriteInrimageWithName( &(par->imphi1), name );
  sprintf( name, "%s.imphi2.inr", genericname );
  VT_WriteInrimageWithName( &(par->imphi2), name );
  sprintf( name, "%s.imphi3.inr", genericname );
  VT_WriteInrimageWithName( &(par->imphi3), name );
  sprintf( name, "%s.iszero.inr", genericname );
  VT_WriteInrimageWithName( &(par->iszero), name );
}

