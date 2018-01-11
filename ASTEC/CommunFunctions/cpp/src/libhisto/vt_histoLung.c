/*************************************************************************
 * vt_histoLung.c -
 *
 * $Id: vt_histoLung.c,v 1.3 2000/07/04 10:38:04 greg Exp $
 *
 * Copyright (c) INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon Jun 26 17:50:05 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <vt_histoLung.h>

void _ProcessHistoLung( vt_image *theHisto,
			double *theGaussParam )
{
  char *proc = "_ProcessHistoLung";
  double *theX = (double*)NULL;
  double *theY = (double*)NULL;
  double *theC = (double*)NULL;
  double *theS = (double*)NULL;
  double *tmpBuf = (double*)NULL;
  int xmax = theHisto->dim.x-1;
  int min, max;
  int i, j;
  int offset;
  double x;
  

  if ( theHisto->dim.y != 1 || theHisto->dim.z != 1 ) {
    VT_Error( "l'entree n'est pas un histogramme", proc );
    return;
  }







  switch( theHisto->type ) {
  default :
    VT_Error( "such input type is not handled yet", proc );
    return;
  case SINT :
    {
      int *theBuf = (int*)theHisto->buf;
      while ( xmax > 0 && theBuf[xmax] == 0 )
	xmax --;
    }
    break;
  case FLOAT :
    {
      float *theBuf = (float*)theHisto->buf;
      while ( xmax > 0 && theBuf[xmax] == 0 )
	xmax --;
    }
    break;
  case DOUBLE :
    {
      double *theBuf = (double*)theHisto->buf;
      while ( xmax > 0 && theBuf[xmax] == 0 )
	xmax --;
    }
    break;
  }




  tmpBuf = (double*)VT_Malloc( 4 * theHisto->dim.x * sizeof(double) );
  if ( tmpBuf == (double*)NULL ) return;
  theX = theY = theC = theS = tmpBuf;
  theY +=   theHisto->dim.x;
  theC += 2*theHisto->dim.x;
  theS += 3*theHisto->dim.x;
  for ( i=0; i<theHisto->dim.x; i++ ) {
    theX[i] = theY[i] = 0.0;
    theC[i] = theS[i] = 1.0;
  }












  /* premier pic vers 250
   */
  min = 100;
  max = 400;
  offset = 0;
  theGaussParam[3*offset+1] = 250;
  theGaussParam[3*offset+2] = 10;

  switch ( theHisto->type ) {
  default :
    VT_Error( "such input type is not handled yet", proc );
    return;
  case SINT :
    {
      int *theBuf = (int*)theHisto->buf;
      for ( i=0, j=min; j<=max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
      theGaussParam[3*offset+0] = theBuf[ (int)(theGaussParam[3*offset+1]) ];
    }
    break;
  case FLOAT :
    {
      float *theBuf = (float*)theHisto->buf;
      for ( i=0, j=min; j<=max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
      theGaussParam[3*offset+0] = theBuf[ (int)(theGaussParam[3*offset+1]) ];
    }
    break;
  }
  
  printf( "\n" );
  printf( " --- INIT #%d : amp=%f   moy=%f   sig=%f\n", offset,
	 theGaussParam[3*offset+0], theGaussParam[3*offset+1], 
	 theGaussParam[3*offset+2] );

  (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				     theC, DOUBLE, theS, DOUBLE, 
				     max-min+1,
				     &(theGaussParam[3*offset]), 3, 
				     _GaussianForLM );
  
  printf( " --- RES. #%d : amp=%f   moy=%f   sig=%f\n", offset,
	 theGaussParam[3*offset+0], theGaussParam[3*offset+1], 
	 theGaussParam[3*offset+2] );
  printf( "\n" );
  





  /* deuxieme pic vers 350
   */
  min = 300;
  max = 400;
  offset = 1;
  theGaussParam[3*offset+1] = 350;
  theGaussParam[3*offset+2] = 10;

  switch ( theHisto->type ) {
  default :
    VT_Error( "such input type is not handled yet", proc );
    return;
  case SINT :
    {
      int *theBuf = (int*)theHisto->buf;
      for ( i=0, j=min; j<=max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
	x = j - theGaussParam[1];
	theY[i] -= theGaussParam[0] * exp( - x * x / ( 2.0 * theGaussParam[2] ) );
      }
      theGaussParam[3*offset+0] = theBuf[ (int)(theGaussParam[3*offset+1]) ];
    }
    break;
  case FLOAT :
    {
      float *theBuf = (float*)theHisto->buf;
      for ( i=0, j=min; j<=max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
	x = j - theGaussParam[1];
	theY[i] -= theGaussParam[0] * exp( - x * x / ( 2.0 * theGaussParam[2] ) );
      }
      theGaussParam[3*offset+0] = theBuf[ (int)(theGaussParam[3*offset+1]) ];
    }
    break;
  }
  theGaussParam[3*offset+0] = 100.0;
  printf( "\n" );
  printf( " --- INIT #%d : amp=%f   moy=%f   sig=%f\n", offset,
	 theGaussParam[3*offset+0], theGaussParam[3*offset+1], 
	 theGaussParam[3*offset+2] );

  (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				     theC, DOUBLE, theS, DOUBLE, 
				     max-min+1,
				     &(theGaussParam[3*offset]), 3, 
				     _GaussianForLM );
  
  printf( " --- RES. #%d : amp=%f   moy=%f   sig=%f\n", offset,
	 theGaussParam[3*offset+0], theGaussParam[3*offset+1], 
	 theGaussParam[3*offset+2] );
  printf( "\n" );
  







  /* les 2 gaussiennes */
  min = 100;
  max = 1400;


  switch ( theHisto->type ) {
  default :
    VT_Error( "such input type is not handled yet", proc );
    return;
  case SINT :
    {
      int *theBuf = (int*)theHisto->buf;
      for ( i=0, j=min; j<=max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
    }
    break;
  case FLOAT :
    {
      float *theBuf = (float*)theHisto->buf;
      for ( i=0, j=min; j<=max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
    }
    break;
  }
  
  /*
  (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				     theC, DOUBLE, theS, DOUBLE, 
				     max-min+1,
				     theGaussParam, 6, 
				     _MixtureOf2GaussiansForLM );

  */

  VT_Free( (void**)&tmpBuf );
  return;









  /* troisieme pic vers 1300
   */
  min = 100;
  max = 300;
  offset = 2;
  theGaussParam[3*offset+1] = 200;
  theGaussParam[3*offset+2] = 10;

  switch ( theHisto->type ) {
  default :
    VT_Error( "such input type is not handled yet", proc );
    return;
  case SINT :
    {
      int *theBuf = (int*)theHisto->buf;
      for ( i=0, j=min; j<=max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
      theGaussParam[3*offset+0] = theBuf[ (int)(theGaussParam[3*offset+1]) ];
    }
    break;
  case FLOAT :
    {
      float *theBuf = (float*)theHisto->buf;
      for ( i=0, j=min; j<=max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
      theGaussParam[3*offset+0] = theBuf[ (int)(theGaussParam[3*offset+1]) ];
    }
    break;
  }
  
  printf( "\n" );
  printf( " --- INIT #%d : amp=%f   moy=%f   sig=%f\n", offset,
	 theGaussParam[3*offset+0], theGaussParam[3*offset+1], 
	 theGaussParam[3*offset+2] );

  (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				     theC, DOUBLE, theS, DOUBLE, 
				     max-min+1,
				     &(theGaussParam[3*offset]), 3, 
				     _GaussianForLM );
  
  printf( " --- RES. #%d : amp=%f   moy=%f   sig=%f\n", offset,
	 theGaussParam[3*offset+0], theGaussParam[3*offset+1], 
	 theGaussParam[3*offset+2] );
  printf( "\n" );
  



  for ( offset=0; offset<3; offset++ ) {
    if ( theGaussParam[3*offset+2] < 0.0 )
      theGaussParam[3*offset+2] *= -1.0;
  }




  







  /* les 3 gaussiennes */
  min = 100;
  max = 1400;


  switch ( theHisto->type ) {
  default :
    VT_Error( "such input type is not handled yet", proc );
    return;
  case SINT :
    {
      int *theBuf = (int*)theHisto->buf;
      for ( i=0, j=min; j<=max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
    }
    break;
  case FLOAT :
    {
      float *theBuf = (float*)theHisto->buf;
      for ( i=0, j=min; j<=max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
    }
    break;
  }
  
  
  (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				     theC, DOUBLE, theS, DOUBLE, 
				     max-min+1,
				     theGaussParam, 9, 
				     _MixtureOf3GaussiansForLM );

  for ( offset=0; offset<3; offset++ ) {
    printf( " --- RES. #%d : amp=%f   moy=%f   sig=%f\n", offset,
	   theGaussParam[3*offset+0], theGaussParam[3*offset+1], 
	   theGaussParam[3*offset+2] );
  }
  printf( "\n" );
















  VT_Free( (void**)&tmpBuf );
}
