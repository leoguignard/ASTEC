
#include <vt_gaussienne.h>

#ifndef NO_PROTO
static void   FitGaussienne( r32 *h, double *K, double *sigma, double *moy, int deb, int fin );
static double Gaussienne( double K, double sigma, double moy, double x );
static double DGaussienneDSigma( double K, double sigma, double moy, double x );
static double DGaussienneDK( double sigma, double moy, double x );
static double DGaussienneDMoy( double K, double sigma, double moy, double x );
static double ChiGaussienne( r32 *h, double K, double sigma, double moy, int deb, int fin );
static void   DeriveChiGaussienne( r32 *h, double K, double sigma, double moy, double Beta[3], int deb, int fin );
static void   Derive2ChiGaussienne( double K, double sigma, double moy, double Alpha[9], int deb, int fin );
static int    GaussienneDsDkDm( double Alpha[9], double Beta[3], double *ds, double *dk, double *dm );
#else 
static void   FitGaussienne();
static double Gaussienne();
static double DGaussienneDSigma();
static double DGaussienneDK();
static double DGaussienneDMoy();
static double ChiGaussienne();
static void   DeriveChiGaussienne();
static void   Derive2ChiGaussienne();
static int    GaussienneDsDkDm();
#endif

static double g_pi = 3.141592653589793238462643383279;

#ifndef NO_PROTO
int VT_FitGaussienne( r32 *tab, int l, double *K, double *sigma, double *moy )
#else
int VT_FitGaussienne( tab, l, K, sigma, moy )
r32 *tab;
int l;
double *K;
double *sigma;
double *moy;
#endif
{
    int i, i_min, i_max;
    float val_max;
    double lK, lS, lM;
    char *proc="VT_FitGaussienne";

    *K = *sigma = *moy = 0.0;

    /*--- calcul des bornes min et max ---*/
    i_min = 0;   i_max = l-1;
    i = 0;     while ( (tab[i] <= 0.1) && (i<l-1) ) i ++;
    i_min = i;
    i = l-1;   while ( (tab[i] <= 0.1) && (i>0) )   i --;
    i_max = i;
    if ( i_min >= i_max ) {
	VT_Error( "array filled with 0", proc );
	return( -1 );
    }

    /*--- calcul de la valeur max ---*/
    val_max = tab[i_min];
    for (i = i_min+1; i <= i_max; i ++ )
	if ( tab[i] > val_max ) val_max = tab[i];
    
    /*--- initialisation de la gaussienne ---*/
    lM = (double)(i_min + i_max) / 2.0;
    /*--- 99 % des points de la gaussienne sont entre deb et fin ---*/
    lS = (double)(i_max - i_min) / 2.0 / 2.58;
    /*--- hauteur de la gaussienne ---*/
    lK = VT_Sqrt( 2.0 * g_pi ) * lS * (double)(val_max);

    FitGaussienne( tab, &lK, &lS, &lM, i_min, i_max );
    
    *K = lK;
    *sigma = lS;
    *moy = lM;

    return( 1 );
}

#ifndef NO_PROTO
static void FitGaussienne( r32 *h, double *K, double *sigma, double *moy, int deb, int fin )
#else
static void FitGaussienne( h, K, sigma, moy, deb, fin )
r32 *h;
double *K, *sigma, *moy;
int deb, fin;
#endif
{
	double kk, ss, mm;
	double lambda, ChiA, ChiAdA;
	double Beta[3], Alpha[9], Alph[9], ds, dk, dm;
	int stop;

	stop = 0;
	lambda = 0.001;
	kk = *K;
	ss = *sigma;
	mm = *moy;
	
	/*--- fit d'une courbe gaussienne    ---*/
	/*--- methode de Levenberg-Marquardt ---*/
	while( stop == 0 ) {
		ChiA = ChiGaussienne( h, kk, ss, mm, deb, fin );
		DeriveChiGaussienne( h, kk, ss, mm, Beta, deb, fin );
		Derive2ChiGaussienne( kk, ss, mm, Alpha, deb, fin );
		Alph[0] = Alpha[0] * (1.0 + lambda);
		Alph[4] = Alpha[4] * (1.0 + lambda);
		Alph[8] = Alpha[8] * (1.0 + lambda);
		Alph[3] = Alph[1] = Alpha[1];
		Alph[6] = Alph[2] = Alpha[2];
		Alph[7] = Alph[5] = Alpha[5];
		if ( GaussienneDsDkDm( Alph, Beta, &ds, &dk, &dm) ) {
			ds += ss;
			dk += kk;
			dm += mm;
			ChiAdA = ChiGaussienne( h, dk, ds, dm, deb, fin );
			if ( ChiAdA >= ChiA ) { lambda *= 10; }
			else {
				lambda /= 10;
				ss = ds;
				kk = dk;
				mm = dm;
				stop = ( (ChiA - ChiAdA)/ChiA < 1e-3 );
			}
			if ( lambda > 1e+10 ) stop = 1;
		}
		else
			stop = 1;
	}
	/*--- retour ---*/
	*K = kk;
	*sigma = ss;
	*moy = mm;
}

#ifndef NO_PROTO
static double Gaussienne( double K, double sigma, double moy, double x )
#else
static double Gaussienne( K, sigma, moy, x )
double K, sigma, moy, x;
#endif
{
  double X;
  double res = 0.0;
  
  if (sigma != 0.0) {
    X = ( (x-moy) * (x-moy) ) / ( 2.0 * sigma * sigma );
    res = (K * exp(-X))/(VT_Sqrt(2.0*g_pi)*sigma) ;
  }
  return(res);
}

#ifndef NO_PROTO
static double DGaussienneDSigma( double K, double sigma, double moy, double x )
#else
static double DGaussienneDSigma( K, sigma, moy, x )
double K, sigma, moy, x;
#endif
{
	double X;
	double res = 0.0;

	if (sigma != 0.0) {	  
	  X = Gaussienne( K, sigma, moy, x );
	  X *= ((x-moy)*(x-moy))/(sigma*sigma) - 1.0;
	  res = X/sigma;
	}
	return(res);
}

#ifndef NO_PROTO
static double DGaussienneDK( double sigma, double moy, double x )
#else
static double DGaussienneDK( sigma, moy, x )
double sigma, moy, x;
#endif
{
	double X;
	double res = 0.0;
	
	if (sigma != 0.0) {
	  X = ( (x-moy) * (x-moy) ) / ( 2.0 * sigma * sigma );
	  res = (exp(-X))/(VT_Sqrt(2.0*g_pi)*sigma) ;
	}
	return (res);
}

#ifndef NO_PROTO
static double DGaussienneDMoy( double K, double sigma, double moy, double x )
#else
static double DGaussienneDMoy( K, sigma, moy, x )
double K, sigma, moy, x;
#endif
{
	double X;
	double res = 0.0;
	
	if (sigma != 0.0) {
	  X = Gaussienne( K, sigma, moy, x );
	  res = (X * (x-moy))/(sigma*sigma) ;
	}
	return (res);
}

#ifndef NO_PROTO
static double ChiGaussienne( r32 *h, double K, double sigma, double moy, int deb, int fin )
#else 
static double ChiGaussienne( h, K, sigma, moy, deb, fin )
r32 *h;
double K, sigma, moy;
int deb, fin;
#endif
{
	register int i;
	double val, sum;

	sum = 0;
	for (i=deb; i <= fin; i++) {
		val = (double)h[i] - Gaussienne( K, sigma, moy, (double)i);
		sum += val*val;
	}
	return( sum );
}

#ifndef NO_PROTO
static void DeriveChiGaussienne( r32 *h, double K, double sigma, double moy, double Beta[3], int deb, int fin )
#else 
static void DeriveChiGaussienne( h, K, sigma, moy, Beta, deb, fin )
r32 *h;
double K, sigma, moy, Beta[3];
int deb, fin;
#endif
{
	register int i;
	double val;
	
	Beta[0] = Beta[1] = Beta[2] = 0.0;
	for (i=deb; i <= fin; i++) {
	    val = (double)h[i] - Gaussienne( K, sigma, moy, (double)i);
	    Beta[0] += val * DGaussienneDSigma( K, sigma, moy, (double)i );
	    Beta[1] += val * DGaussienneDK( sigma, moy, (double)i );
	    Beta[2] += val * DGaussienneDMoy( K, sigma, moy, (double)i );
	}
}

#ifndef NO_PROTO
static void Derive2ChiGaussienne( double K, double sigma, double moy, double Alpha[9], int deb, int fin )
#else 
static void Derive2ChiGaussienne( K, sigma, moy, Alpha, deb, fin )
double K, sigma, moy, Alpha[9];
int deb, fin;
#endif
{
	register int i;
	double ds, dk, dm;
	
	Alpha[0] = Alpha[1] = Alpha[2] = Alpha[3] = Alpha[4] = 0.0;
	Alpha[5] = Alpha[6] = Alpha[7] = Alpha[8] = 0.0;
	for (i=(int)(deb); i <= (int)(fin); i++) {
		ds = DGaussienneDSigma( K, sigma, moy, (double)i );
		dk = DGaussienneDK( sigma, moy, (double)i );
		dm = DGaussienneDMoy( K, sigma, moy, (double)i );
		Alpha[0] += ds*ds;
		Alpha[1] += ds*dk;
		Alpha[2] += ds*dm;
		Alpha[4] += dk*dk;
		Alpha[5] += dk*dm;
		Alpha[8] += dm*dm;
	}
	Alpha[3] = Alpha[1];
	Alpha[6] = Alpha[2];
	Alpha[7] = Alpha[5];
}
#ifndef NO_PROTO
static int GaussienneDsDkDm( double Alpha[9], double Beta[3], double *ds, double *dk, double *dm )
#else
static int GaussienneDsDkDm( Alpha, Beta, ds, dk, dm )
double Alpha[9], Beta[3];
double *ds, *dk, *dm;
#endif
{
	double delta, com[9];

	*ds = *dk = *dm = 0.0;
		
	delta = Alpha[0]*Alpha[4]*Alpha[8] + Alpha[3]*Alpha[7]*Alpha[2] + Alpha[6]*Alpha[1]*Alpha[5];
	delta -= Alpha[2]*Alpha[4]*Alpha[6] + Alpha[5]*Alpha[7]*Alpha[0] + Alpha[8]*Alpha[1]*Alpha[3];
	if (delta == 0.0) return( 0 );
	
	com[0] = Alpha[4]*Alpha[8] - Alpha[5]*Alpha[7];
	com[1] = - (Alpha[3]*Alpha[8] - Alpha[5]*Alpha[6]);
	com[2] = Alpha[3]*Alpha[7] - Alpha[4]*Alpha[6];
	com[3] = - (Alpha[1]*Alpha[8] - Alpha[2]*Alpha[7]);
	com[4] = Alpha[0]*Alpha[8] - Alpha[2]*Alpha[6];
	com[5] = - (Alpha[0]*Alpha[7] - Alpha[1]*Alpha[6]);
	com[6] = Alpha[1]*Alpha[5] - Alpha[2]*Alpha[4];
	com[7] = - (Alpha[0]*Alpha[5] - Alpha[2]*Alpha[3]);
	com[8] = Alpha[0]*Alpha[4] - Alpha[1]*Alpha[3];

	*ds = com[0]*Beta[0] + com[3]*Beta[1] + com[6]*Beta[2];
	*dk = com[1]*Beta[0] + com[4]*Beta[1] + com[7]*Beta[2];
	*dm = com[2]*Beta[0] + com[5]*Beta[1] + com[8]*Beta[2];

	*ds /= delta;
	*dk /= delta;
	*dm /= delta;
	return( 1 );
}
