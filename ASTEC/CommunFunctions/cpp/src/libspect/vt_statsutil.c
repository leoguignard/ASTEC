/*************************************************************************
 * vt_statsutil.c -
 *
 * $Id: vt_statsutil.c,v 1.2 2000/06/29 16:48:21 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Thu Mar 23 11:59:08 MET 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vt_statsutil.h>



static void _InitListOfTypePoint( typePoint *theList,
				  int first, int last )
{
  int i;
  for ( i=first; i<=last; i++ ) {

    theList[i].sum0 = 0;
    theList[i].sum1 = 0;
    theList[i].sum2 = 0;
    theList[i].sum3 = 0;
    theList[i].sum4 = 0;

    theList[i].moy = 0;
    theList[i].variance = 0;
    theList[i].ecarttype = 0;
    theList[i].skewness = 0;
    theList[i].kurtosis = 0;
 
    theList[i].med = 0;
    theList[i].indMax = 0;
    theList[i].valMax = 0;

    theList[i].gauss[0] = 0;
    theList[i].gauss[1] = 0;
    theList[i].gauss[2] = 0;

    theList[i].nsgauss[0] = 0;
    theList[i].nsgauss[1] = 0;
    theList[i].nsgauss[2] = 0;
    theList[i].nsgauss[3] = 0;

  }
}










void _ComputeStatsValues( vt_image *imHisto,
			  typePoint *theXenon,
			  typePoint *theHmpao,
			  int theMinXenon, int theMaxXenon,
			  int theMinHmpao, int theMaxHmpao )
{
  int minXenon = theMinXenon;
  int maxXenon = theMaxXenon;
  int minHmpao = theMinHmpao;
  int maxHmpao = theMaxHmpao;
  
  float ***theHisto = (float***)NULL;
  int hmpao, xenon;
  int x, y;

  double sum, max;

  double *theX, *theY, *theC, *theS;
  double *theAllocatedBuffer = (double*)NULL;

  int i;
  

  if ( imHisto->type != FLOAT ) return;
  theHisto = (float***)imHisto->array;



  x = imHisto->dim.y;
  if ( x < imHisto->dim.x ) x = imHisto->dim.x;

  theAllocatedBuffer = (double*)malloc( 4 * x * sizeof(double) );
  if ( theAllocatedBuffer == (double*)NULL ) {
    return;
  }
  
  theX = theY = theC = theS = theAllocatedBuffer;
  theY += x;
  theC += 2*x;
  theS += 3*x;




  if ( minXenon < 0 ) minXenon = 0;
  if ( maxXenon >= imHisto->dim.x ) maxXenon = imHisto->dim.x-1;

  if ( minHmpao < 0 ) minHmpao = 0;
  if ( maxHmpao >= imHisto->dim.y ) maxHmpao = imHisto->dim.y-1;



  _InitListOfTypePoint( theXenon, 0, imHisto->dim.x-1 );
  _InitListOfTypePoint( theHmpao, 0, imHisto->dim.y-1 );






  /* Calcul de la moyenne, de la mediane et du maximum
     
     variance = 1/N sum( n_i * x_i * x_i ) - x_moy * x_moy
     ecart-type = sigma = (standard deviation) = sqrt(variance)           
     
     variance = moment centre d'ordre 2 = E[ (x-m)^2 ]

     coefficient d'asymetrie (skewness) = E[ (x-m)^3 ] / sigma^3
     
     coefficient d'aplatissement (kurtosis) = E[ (x-m)^4 ] / sigma^4

     Rayon de courbure pour une courbe y=f(x)
     R = epsilon (1+y'^2)^(3/2)/y'')
  */
  



  for ( hmpao = 0; hmpao < imHisto->dim.y; hmpao++ ) {
    
    /* moments  de la distribution
     */
    
    for ( xenon = 0; xenon < imHisto->dim.x; xenon++ ) {
      theHmpao[hmpao].sum0 += theHisto[0][hmpao][xenon];
      theHmpao[hmpao].sum1 += xenon * theHisto[0][hmpao][xenon];
      theHmpao[hmpao].sum2 += xenon * xenon * theHisto[0][hmpao][xenon];
      theHmpao[hmpao].sum3 += xenon * xenon * xenon * theHisto[0][hmpao][xenon];
      theHmpao[hmpao].sum4 += xenon * xenon * xenon * xenon * theHisto[0][hmpao][xenon];
    }

    /* grandeurs statistiques
     */

    if ( theHmpao[hmpao].sum0 > _MINVAL_ ) {
      theHmpao[hmpao].moy = theHmpao[hmpao].sum1 / theHmpao[hmpao].sum0;
      theHmpao[hmpao].variance = theHmpao[hmpao].sum2 / theHmpao[hmpao].sum0 
	- theHmpao[hmpao].moy * theHmpao[hmpao].moy;
      if ( theHmpao[hmpao].variance < 0.0 ) theHmpao[hmpao].variance = 0.0;
      theHmpao[hmpao].ecarttype = sqrt( theHmpao[hmpao].variance );

      if ( theHmpao[hmpao].ecarttype > 0 ) {
	theHmpao[hmpao].skewness = theHmpao[hmpao].sum3 / theHmpao[hmpao].sum0
	  - 3 * theHmpao[hmpao].moy * theHmpao[hmpao].sum2 / theHmpao[hmpao].sum0
	  + 2 * theHmpao[hmpao].moy * theHmpao[hmpao].moy * theHmpao[hmpao].moy;
	theHmpao[hmpao].skewness /= theHmpao[hmpao].ecarttype * theHmpao[hmpao].ecarttype * theHmpao[hmpao].ecarttype;
	theHmpao[hmpao].kurtosis = theHmpao[hmpao].sum4 / theHmpao[hmpao].sum0
	  - 4 * theHmpao[hmpao].moy * theHmpao[hmpao].sum3 / theHmpao[hmpao].sum0
	  + 6 * theHmpao[hmpao].moy * theHmpao[hmpao].sum2 / theHmpao[hmpao].sum0
	  - 3 * theHmpao[hmpao].moy * theHmpao[hmpao].moy * theHmpao[hmpao].moy * theHmpao[hmpao].moy;
	theHmpao[hmpao].kurtosis /= theHmpao[hmpao].ecarttype * theHmpao[hmpao].ecarttype;
	theHmpao[hmpao].kurtosis /= theHmpao[hmpao].ecarttype * theHmpao[hmpao].ecarttype;
      }
    }

    /* mediane
     */
    if ( theHmpao[hmpao].sum0 > _MINVAL_ ) {
      sum = 0.0;
      for ( xenon = 0; xenon < imHisto->dim.x - 1; xenon++ ) {
	sum += theHisto[0][hmpao][xenon];
	if ( 2.0*sum < theHmpao[hmpao].sum0 
	     && 2.0*(sum+theHisto[0][hmpao][xenon+1]) >= theHmpao[hmpao].sum0  ) {
	  theHmpao[hmpao].med = xenon + ( theHmpao[hmpao].sum0/2.0 - sum) / theHisto[0][hmpao][xenon+1];
	}
      }
    }
    
    /* maximum
     */
    max = theHisto[0][hmpao][0];
    x = 0;
    for ( xenon = 1; xenon < imHisto->dim.x; xenon++ ) {
      if ( max < theHisto[0][hmpao][xenon] ) {
	max = theHisto[0][hmpao][xenon];
	x = xenon;
      }
    }
    theHmpao[hmpao].indMax = x;
    theHmpao[hmpao].valMax = max; 
    
    
    /* fit de distributions
       sur la distribution de xenon pour un  hmpao fixe
     */
    if ( hmpao < minHmpao || hmpao > maxHmpao ) continue;
    
    for ( i=0, xenon = 0; xenon < imHisto->dim.x; xenon++ ) {
      if ( theHisto[0][hmpao][xenon] > _MINVAL_  ) {
	theX[i] = xenon;
	theY[i] = theHisto[0][hmpao][xenon];
	theC[i] = 1.0;
	theS[i] = 1.0;
	i++;
      }
    }
    
    theHmpao[hmpao].gauss[0] = theHisto[0][hmpao][ (int)( theHmpao[hmpao].moy+0.5 ) ];
    theHmpao[hmpao].gauss[1] = theHmpao[hmpao].moy;
    theHmpao[hmpao].gauss[2] = theHmpao[hmpao].ecarttype;
    (void)VT_Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, theC, DOUBLE, theS, DOUBLE, i,
					  theHmpao[hmpao].gauss, 3, _GaussianForLM );
    
    theHmpao[hmpao].nsgauss[0] = theHisto[0][hmpao][ (int)( theHmpao[hmpao].moy+0.5 ) ];
    theHmpao[hmpao].nsgauss[1] = theHmpao[hmpao].moy;
    theHmpao[hmpao].nsgauss[2] = theHmpao[hmpao].ecarttype;
    theHmpao[hmpao].nsgauss[3] = theHmpao[hmpao].ecarttype;
    (void)VT_Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, theC, DOUBLE, theS, DOUBLE, i,
					  theHmpao[hmpao].nsgauss, 4, _NonSymetricGaussianForLM );
  }



  
  /* XENON Calcul de la moyenne, de la mediane et du maximum
   */


  for ( xenon = 0; xenon < imHisto->dim.x; xenon++ ) {
    
    /* moments  de la distribution
     */
    
    for ( hmpao = 0; hmpao < imHisto->dim.y; hmpao++ ) {
      theXenon[xenon].sum0 += theHisto[0][hmpao][xenon];
      theXenon[xenon].sum1 += hmpao * theHisto[0][hmpao][xenon];
      theXenon[xenon].sum2 += hmpao * hmpao * theHisto[0][hmpao][xenon];
      theXenon[xenon].sum3 += hmpao * hmpao * hmpao * theHisto[0][hmpao][xenon];
      theXenon[xenon].sum4 += hmpao * hmpao * hmpao * hmpao * theHisto[0][hmpao][xenon];
    }
   
    /* grandeurs statistiques
     */

    if ( theXenon[xenon].sum0 > _MINVAL_ ) {
      theXenon[xenon].moy = theXenon[xenon].sum1 / theXenon[xenon].sum0;
      theXenon[xenon].variance = theXenon[xenon].sum2 / theXenon[xenon].sum0 
	- theXenon[xenon].moy * theXenon[xenon].moy;
      if ( theXenon[xenon].variance < 0 ) theXenon[xenon].variance = 0.0;
      theXenon[xenon].ecarttype = sqrt( theXenon[xenon].variance );
      
      if ( theXenon[xenon].ecarttype > 0 ) {
	theXenon[xenon].skewness = theXenon[xenon].sum3 / theXenon[xenon].sum0
	  - 3 * theXenon[xenon].moy * theXenon[xenon].sum2 / theXenon[xenon].sum0
	  + 2 * theXenon[xenon].moy * theXenon[xenon].moy * theXenon[xenon].moy;
	theXenon[xenon].skewness /= theXenon[xenon].ecarttype * theXenon[xenon].ecarttype * theXenon[xenon].ecarttype;
	theXenon[xenon].kurtosis = theXenon[xenon].sum4 / theXenon[xenon].sum0
	  - 4 * theXenon[xenon].moy * theXenon[xenon].sum3 / theXenon[xenon].sum0
	  + 6 * theXenon[xenon].moy * theXenon[xenon].sum2 / theXenon[xenon].sum0
	  - 3 * theXenon[xenon].moy * theXenon[xenon].moy * theXenon[xenon].moy * theXenon[xenon].moy;
	theXenon[xenon].kurtosis /= theXenon[xenon].ecarttype * theXenon[xenon].ecarttype;
	theXenon[xenon].kurtosis /= theXenon[xenon].ecarttype * theXenon[xenon].ecarttype;
      }
    }

    /* mediane
     */
    if ( theXenon[xenon].sum0 > _MINVAL_ ) {
      sum = 0.0;
      for ( hmpao = 0; hmpao < imHisto->dim.y-1; hmpao++ ) {
	sum += theHisto[0][hmpao][xenon];
	if ( 2.0*sum < theXenon[xenon].sum0 
	     && 2.0*(sum+theHisto[0][hmpao+1][xenon]) >= theXenon[xenon].sum0 ) {
	  theXenon[xenon].med = hmpao + ( theXenon[xenon].sum0/2.0 - sum) / theHisto[0][hmpao+1][xenon];
	}
      }
    }

    /* maximum
     */
    max = theHisto[0][0][xenon];
    y = 0;
    for ( hmpao = 1; hmpao < imHisto->dim.y; hmpao++ ) {
      if ( max < theHisto[0][hmpao][xenon] ) {
	max = theHisto[0][hmpao][xenon];
	y = hmpao;
      }
    }
    theXenon[xenon].indMax = y;
    theXenon[xenon].valMax = max;    
    


    /* fit de distributions
       sur la distribution de xenon pour un  hmpao fixe
     */
    if ( xenon < minXenon || xenon > maxXenon ) continue;

    for ( i=0, hmpao = 0; hmpao < imHisto->dim.y; hmpao++ ) {
      if ( theHisto[0][hmpao][xenon] > _MINVAL_  ) {
	theX[hmpao] = hmpao;
	theY[hmpao] = theHisto[0][hmpao][xenon];
	theC[hmpao] = 1.0;
	theS[hmpao] = 1.0;
	i ++;
      }
    }

    theXenon[xenon].gauss[0] = theHisto[0][ (int)( theXenon[xenon].moy+0.5 ) ][xenon];
    theXenon[xenon].gauss[1] = theXenon[xenon].moy;
    theXenon[xenon].gauss[2] = theXenon[xenon].ecarttype;
    (void)VT_Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, theC, DOUBLE, theS, DOUBLE, i,
					  theXenon[xenon].gauss, 3, _GaussianForLM );
    
    theXenon[xenon].nsgauss[0] = theHisto[0][ (int)( theXenon[xenon].moy+0.5 ) ][xenon];
    theXenon[xenon].nsgauss[1] = theXenon[xenon].moy;
    theXenon[xenon].nsgauss[2] = theXenon[xenon].ecarttype;
    theXenon[xenon].nsgauss[3] = theXenon[xenon].ecarttype;
    (void)VT_Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, theC, DOUBLE, theS, DOUBLE, i,
					  theXenon[xenon].nsgauss, 4, _NonSymetricGaussianForLM );
    
    if ( 0 ) fprintf( stdout, " XENON #%2d : AMP = %f MOY = %f EC1 =%f EC2 =%f\n",
	     xenon, theXenon[xenon].nsgauss[0], theXenon[xenon].nsgauss[1],
	     theXenon[xenon].nsgauss[2], theXenon[xenon].nsgauss[3] );

  }
  
  free( theAllocatedBuffer );
}
