#include <vt_recUS2D.h>

typedef struct {
  /* coordonnees cartesiennes voxels entieres
   */
  int xve;
  int yve;
  /* coordonnees cartesiennes voxels reelles
   */
  double xvr;
  double yvr;
  /* coordonnees millimetriques reelles
   */
  double xmr;
  double ymr;

  /* coordonnees polaires (voxels images) entieres
   */
  int rve;
  int tve;

  /* coordonnees polaires (voxels images) reelles
   */
  double rvr;
  double tvr;

  /* coordonnees polaires millimetriques et angulaires (degres)
   */
  double rmr;
  double tdr;

} typePoint;

 
int _ReconstructUS2D( vt_image *theUs,
		       vt_image *theRec,
		       typeGeometryUS2D *theGeom )
{
  char *proc = "_ReconstructUS2D";
  unsigned char ***bufRec = (unsigned char ***)theRec->array;
  unsigned char ***bufUs  = (unsigned char ***)theUs->array;
  /* c'est en voxels
     le point d'origine des rayons est (xMiddle,0)
   */
  double xMiddle = ((double)theRec->dim.x-1)/2.0;
  double xRelative;

  double vx = theRec->siz.x * theRec->siz.x;

  typePoint cur, v[4];
  int i;
  double sum, coef[4];


  if ( theUs->type != UCHAR || theRec->type != UCHAR ) {
    VT_Error( "bad images type", proc );
    return( -1 );
  }

  
  for ( cur.yve = 0; cur.yve < theRec->dim.y; cur.yve ++ )
  for ( cur.xve = 0; cur.xve < theRec->dim.x; cur.xve ++ ) {

    bufRec[0][cur.yve][cur.xve] = 0;

    xRelative = cur.xve - xMiddle;

    /* c'est en mm
     */
    cur.xmr = cur.xve * theRec->siz.x;
    cur.ymr = cur.yve * theRec->siz.y;

    cur.rmr = sqrt( xRelative * xRelative * vx + cur.ymr * cur.ymr );
    if ( cur.rmr < theGeom->radiusMin || cur.rmr > theGeom->radiusMax ) continue;

    /* c'est en degres
     */
    cur.tdr = asin( (cur.xmr - xMiddle * theRec->siz.x) / cur.rmr ) * 180.0/3.1415927;
    if ( cur.tdr < theGeom->thetaMin || cur.tdr > theGeom->thetaMax ) continue;
    
    /* on passe en coordonnees voxels reelles
     */
    cur.rvr = (theUs->dim.x - 1) * (cur.rmr - theGeom->radiusMin) / 
      (theGeom->radiusMax - theGeom->radiusMin);
    cur.tvr = (theUs->dim.y - 1) * (cur.tdr - theGeom->thetaMin) / 
      (theGeom->thetaMax - theGeom->thetaMin);
    
    /* coin "superieur droit"
       v[0] = (tve, rve)                 v[2] = (tve+1, rve) 
                           cur (tvr,rvr)
       v[1] = (tve, rve+1)               v[3] = (tve+1, rve+1) 
     */
    v[0].rve = (int)cur.rvr;
    v[0].tve = (int)cur.tvr;
    if ( v[0].rve < 0 || v[0].rve >= theUs->dim.x-1 ) continue;
    if ( v[0].tve < 0 || v[0].tve >= theUs->dim.y-1 ) continue;
    
    v[1] = v[0];
    v[1].rve ++;
    
    v[2] = v[0];
    v[2].tve ++;
    
    v[3] = v[1];
    v[3].tve ++;


    /* les coefficients sont ceux de l'interpolation bilineaire
       dans l'espace polaire
    */
    coef[0] = (1 -(cur.rvr - v[0].rve)) * (1 -(cur.tvr - v[0].tve));
    coef[1] = (cur.rvr - v[0].rve)      * (1 -(cur.tvr - v[0].tve));
    coef[2] = (1 -(cur.rvr - v[0].rve)) * (cur.tvr - v[0].tve);
    coef[3] = (cur.rvr - v[0].rve)      * (cur.tvr - v[0].tve);
    
    sum = 0.0;
    for ( i=0; i<4; i++ )
      sum += coef[i] *bufUs[0][v[i].tve][v[i].rve];

    bufRec[0][cur.yve][cur.xve] = (int)( sum + 0.5 );

  }

  return( 1 );
}
