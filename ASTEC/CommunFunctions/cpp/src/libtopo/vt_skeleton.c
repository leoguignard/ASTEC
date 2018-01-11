/*************************************************************************
 * vT_skeleton.c
 *
 * $Id: vt_skeleton.c,v 1.1 2000/07/26 07:25:02 greg Exp $
 *
 * Copyright (c) INRIA 1999
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue Jul 25 21:39:22 MEST 2000
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <vt_skeleton.h>



#define _NORME2D_(n,z,y,x) n  = ((int)vx[z][y][x]) * ((int)vx[z][y][x]); \
			   n += ((int)vy[z][y][x]) * ((int)vy[z][y][x]);
#define _NORME3D_(n,z,y,x) n  = ((int)vx[z][y][x]) * ((int)vx[z][y][x]); \
			   n += ((int)vy[z][y][x]) * ((int)vy[z][y][x]); \
			   n += ((int)vz[z][y][x]) * ((int)vz[z][y][x])



#define _SCALAIRE2D_(c,b,a) scalaire  = ((double)vx[z][y][x]) * ((double)vx[c][b][a]); \
			    scalaire += ((double)vy[z][y][x]) * ((double)vy[c][b][a]); \
			    scalaire /= sqrt( (double)(norme * n) ); \
                            angle = acos( scalaire ) * 57.2957795
#define _SCALAIRE3D_(c,b,a) scalaire  = ((double)vx[z][y][x]) * ((double)vx[c][b][a]); \
			    scalaire += ((double)vy[z][y][x]) * ((double)vy[c][b][a]); \
			    scalaire += ((double)vz[z][y][x]) * ((double)vz[c][b][a]); \
			    scalaire /= sqrt( (double)(norme * n) ); \
                            angle = acos( scalaire ) * 57.2957795
/* 57.2957795 = 180.0 / 3.1415926536 */




#define _DIST2D_(c,b,a) dx = (double)vx[c][b][a] - (double)vx[z][y][x]; \
                        dy = (double)vy[c][b][a] - (double)vy[z][y][x]; 

#define _DIST3D_(c,b,a) dx = (double)vx[c][b][a] - (double)vx[z][y][x]; \
                        dy = (double)vy[c][b][a] - (double)vy[z][y][x]; \
                        dz = (double)vz[c][b][a] - (double)vz[z][y][x]






#define _CHOIX_ switch ( typeCoefficient ) { \
		case  _DISTANCE_ : \
		  r = dist; \
		  break; \
		case _PRODUIT_ : \
		  r = angle * dist; \
		  break; \
		case _PRODUIT_LOG_ : \
		  if ( dist > 0.0 ) r = angle * log( (double)dist ); \
		  else r = 0.0; \
		  break; \
		case _LOG_PRODUIT_ : \
		  if ( angle > 0.0 ) r = log( (double)angle ) * dist ; \
		  else r = 0.0; \
		  break; \
		case _PRODUIT_2_ : \
		  if ( (dist > 0.0) ) \
		    r = angle * exp( (double)( .2630344061 * log( (double)dist ) ) ); \
		  else r = 0.0; \
		  break; \
		case _PRODUIT_3_ : \
		  if ( (dist > 0.0) ) \
		    r = angle * exp( (double)( .3690702462 * log( (double)dist ) ) ); \
		  else r = 0.0; \
		  break; \
		case _PRODUIT_4_ : \
		  if ( (dist > 0.0) ) \
		    r = angle * sqrt( (double)dist ); \
		  else r = 0.0; \
		  break; \
		case _PRODUIT_5_ : \
		  if ( (dist > 0.0) ) \
		    r = angle * exp( (double)( 1.2223924213 * log( (double)dist ) ) ); \
		  else r = 0.0; \
		  break; \
		case _PRODUIT_6_ : \
		  if ( (dist > 0.0) ) \
		    r = angle * exp( (double)( 1.4150374992 * log( (double)dist ) ) ); \
		  else r = 0.0; \
		  break; \
		case _PRODUIT_7_ : \
		  if ( (dist > 0.0) ) \
		    r = angle * exp( (double)( 1.5 * log( (double)dist ) ) ); \
		  else r = 0.0; \
		  break; \
		case _ANGLE_ : \
		default : \
		  r = angle; \
		  break; \
		}

                  




typedef struct {
  int x;
  int y;
  int z;
} typeNeighbor;





int VT_Compute2DSkeletonCoefficient( vt_image *imx,
				      vt_image *imy,
				      vt_image *imc,
				      vt_image *imn,
				      enumCoefficient typeCoefficient )
{
  char *proc = "VT_Compute2DSkeletonCoefficient";
  int x, y, z;
  int dimx, dimy, dimz;

  unsigned char ***theNeighbor = (unsigned char ***)imn->array;
  float ***theResult = (float ***)imc->array;

  double r, scalaire, angle, dist;
  double dx, dy;
  int norme, n;

  int j, nbNeighbors;
  typeNeighbor theNeighbors[4];

  
  nbNeighbors = 4;
  theNeighbors[0].x = -1;   theNeighbors[0].y =  0;
  theNeighbors[1].x =  1;   theNeighbors[1].y =  0;
  theNeighbors[2].x =  0;   theNeighbors[2].y = -1;
  theNeighbors[3].x =  0;   theNeighbors[3].y =  1;






  if ( imx->type != imy->type ) {
    fprintf( stderr, "%s: can not deal with different vector image type\n", proc );
    return( -1 );
  }


  dimx = imx->dim.x;
  dimy = imx->dim.y;
  dimz = imx->dim.z;



  for ( z = 0; z < dimz; z ++ )
  for ( y = 0; y < dimy; y ++ )
  for ( x = 0; x < dimx; x ++ ) {
    theNeighbor[z][y][x] = _HHH_;
    theResult[z][y][x] = 0.0;
  }


  
  switch( imx->type ) {
  case SCHAR :
    {
      s8 *** vx = (s8 ***)imx->array;
      s8 *** vy = (s8 ***)imy->array;

      for ( z = 0; z < dimz; z ++ )
      for ( y = 0; y < dimy; y ++ )
      for ( x = 0; x < dimx; x ++ ) {
	if ( (vx[z][y][x] == 0) && (vy[z][y][x] == 0) ) continue;
	_NORME2D_(norme,z,y,x);

	/* check if it is a possible maximal ball
	   norme = |(vx[M],  vy[M])|
	   n     = |(vx[M+u],vy[M+u])|
	   
	   on teste si    sqrt(norme)+1 <= sqrt(n)
	   si oui, ce n'est un centre d'une boule maximale
	*/
	
	for ( j=0; j<nbNeighbors; j++ ) {
	  if ( x+theNeighbors[j].x < 0 ||  x+theNeighbors[j].x >= dimx ||
	       y+theNeighbors[j].y < 0 ||  y+theNeighbors[j].y >= dimy ) continue;
	  if ( vx[z][y+theNeighbors[j].y][x+theNeighbors[j].x] == 0 &&
	       vy[z][y+theNeighbors[j].y][x+theNeighbors[j].x] == 0 ) continue;
	  _NORME2D_(n, z, y+theNeighbors[j].y, x+theNeighbors[j].x);
	  if ( sqrt( (double)norme ) + 1.0 <= sqrt( (double)n ) ) continue;
	}

	
	
	/*  si c'est un centre
	    on calcule la valeur
	 */
	for ( j=0; j<nbNeighbors; j++ ) {
	  if ( x+theNeighbors[j].x < 0 ||  x+theNeighbors[j].x >= dimx ||
	       y+theNeighbors[j].y < 0 ||  y+theNeighbors[j].y >= dimy ) continue;
	  if ( vx[z][y+theNeighbors[j].y][x+theNeighbors[j].x] == 0 &&
	       vy[z][y+theNeighbors[j].y][x+theNeighbors[j].x] == 0 ) continue;
	  _NORME2D_(n, z, y+theNeighbors[j].y, x+theNeighbors[j].x);
	  if ( (theNeighbors[j].x == -1 && (norme + 2.0 * (double)vx[z][y][x] + 1.0) <= n) ||
	       (theNeighbors[j].x ==  1 && (norme - 2.0 * (double)vx[z][y][x] + 1.0) <= n) ||
	       (theNeighbors[j].y == -1 && (norme + 2.0 * (double)vy[z][y][x] + 1.0) <= n) ||
	       (theNeighbors[j].y ==  1 && (norme - 2.0 * (double)vy[z][y][x] + 1.0) <= n) )
	    continue;
	  _SCALAIRE2D_(z, y+theNeighbors[j].y, x+theNeighbors[j].x);
	  if ( (angle < 91.0) && (norme == 1) && (n == 1) ) angle = 0.0;
	  _DIST2D_(z, y+theNeighbors[j].y, x+theNeighbors[j].x);
	  dx += theNeighbors[j].x;
	  dy += theNeighbors[j].y;
	  dist = sqrt( dx*dx + dy*dy  );
	  _CHOIX_
	  if ( r > theResult[z][y][x] ) { 
	    theResult[z][y][x] = (r32)( r );
	    if ( theNeighbors[j].x == -1 ) {
	      theNeighbor[z][y][x] = _HHB_;
	    } else if ( theNeighbors[j].x ==  1 ) {
	      theNeighbor[z][y][x] = _HHF_;
	    } else if ( theNeighbors[j].y == -1 ) {
	      theNeighbor[z][y][x] = _HBH_;
	    } else if ( theNeighbors[j].y ==  1 ) {
	      theNeighbor[z][y][x] = _HFH_;	 
	    }
	  }
	  
	}

      }
    }
    break;

  case SSHORT :
    {
      s16 *** vx = (s16 ***)imx->array;
      s16 *** vy = (s16 ***)imy->array;

      for ( z = 0; z < dimz; z ++ )
      for ( y = 0; y < dimy; y ++ )
      for ( x = 0; x < dimx; x ++ ) {
	if ( (vx[z][y][x] == 0) && (vy[z][y][x] == 0) ) continue;
	_NORME2D_(norme,z,y,x);

	/* check if it is a possible maximal ball
	   norme = |(vx[M],  vy[M])|
	   n     = |(vx[M+u],vy[M+u])|
	   
	   on teste si    sqrt(norme)+1 <= sqrt(n)
	   si oui, ce n'est un centre d'une boule maximale
	*/
	
	for ( j=0; j<nbNeighbors; j++ ) {
	  if ( x+theNeighbors[j].x < 0 ||  x+theNeighbors[j].x >= dimx ||
	       y+theNeighbors[j].y < 0 ||  y+theNeighbors[j].y >= dimy ) continue;
	  if ( vx[z][y+theNeighbors[j].y][x+theNeighbors[j].x] == 0 &&
	       vy[z][y+theNeighbors[j].y][x+theNeighbors[j].x] == 0 ) continue;
	  _NORME2D_(n, z, y+theNeighbors[j].y, x+theNeighbors[j].x);
	  if ( sqrt( (double)norme ) + 1.0 <= sqrt( (double)n ) ) continue;
	}

	
	
	/*  si c'est un centre
	    on calcule la valeur
	 */
	for ( j=0; j<nbNeighbors; j++ ) {
	  if ( x+theNeighbors[j].x < 0 ||  x+theNeighbors[j].x >= dimx ||
	       y+theNeighbors[j].y < 0 ||  y+theNeighbors[j].y >= dimy ) continue;
	  if ( vx[z][y+theNeighbors[j].y][x+theNeighbors[j].x] == 0 &&
	       vy[z][y+theNeighbors[j].y][x+theNeighbors[j].x] == 0 ) continue;
	  _NORME2D_(n, z, y+theNeighbors[j].y, x+theNeighbors[j].x);
	  if ( (theNeighbors[j].x == -1 && (norme + 2.0 * (double)vx[z][y][x] + 1.0) <= n) ||
	       (theNeighbors[j].x ==  1 && (norme - 2.0 * (double)vx[z][y][x] + 1.0) <= n) ||
	       (theNeighbors[j].y == -1 && (norme + 2.0 * (double)vy[z][y][x] + 1.0) <= n) ||
	       (theNeighbors[j].y ==  1 && (norme - 2.0 * (double)vy[z][y][x] + 1.0) <= n) )
	    continue;
	  _SCALAIRE2D_(z, y+theNeighbors[j].y, x+theNeighbors[j].x);
	  if ( (angle < 91.0) && (norme == 1) && (n == 1) ) angle = 0.0;
	  _DIST2D_(z, y+theNeighbors[j].y, x+theNeighbors[j].x);
	  dx += theNeighbors[j].x;
	  dy += theNeighbors[j].y;
	  dist = sqrt( dx*dx + dy*dy  );
	  _CHOIX_
	  if ( r > theResult[z][y][x] ) { 
	    theResult[z][y][x] = (r32)( r );
	    if ( theNeighbors[j].x == -1 ) {
	      theNeighbor[z][y][x] = _HHB_;
	    } else if ( theNeighbors[j].x ==  1 ) {
	      theNeighbor[z][y][x] = _HHF_;
	    } else if ( theNeighbors[j].y == -1 ) {
	      theNeighbor[z][y][x] = _HBH_;
	    } else if ( theNeighbors[j].y ==  1 ) {
	      theNeighbor[z][y][x] = _HFH_;	 
	    }
	  }
	  
	}

      }
    }
    break;

  default :
    fprintf( stderr, "%s: can not deal with such type\n", proc );
    return( -1 );
  }

  return( 1 );
}
				    
