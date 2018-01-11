/***************************************************************************
 * powell.c
 *
 * $Author: greg $
 * $Revision: 1.1 $
 * $Log: powell.c,v $
 * Revision 1.1  2002/10/18 18:02:07  greg
 * *** empty log message ***
 *
 * Revision 1.4  2001/12/11 09:13:47  greg
 * Simplification des interfaces
 *
 * Revision 1.3  2001/12/07 15:47:52  greg
 * Ajout du rigid 2D / compilation en C++
 *
 *
 *
 * $Id: powell.c,v 1.1 2002/10/18 18:02:07 greg Exp $
 ***************************************************************************/


/* ----------------------------------------------------------------------

   Implementation de la methode de MINIMISATION de Powell. Adapte des 
   Numerical Recipies.
   
   Alexis Roche, 1997-98.

   ----------------------------------------------------------------------- */

 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



#include <powell.h>

/*--- variables statiques
  ---*/

static int _VERBOSE_ = 1;
static FILE* mystdout = NULL;



double Brent (double ax, 
		  double bx, 
		  double cx, 
		  double tol,
		  double *xmin,
		  double *pt, /* le point */
		  double *di, /* la direction */
		  double *aux, 
		  int n,
		  double (*func)(double *, void *),
		  void *par
		  )
{
  int i, iter;
  double a,b;
  double x,v,w;
  double u;
  double xm;
  double d;
  double etemp,fu,fv,fw,fx,p,q,r,tol1,tol2;
  double e=0.0;
  
  /* --- Significations des variables :
     a,b   : intervalle dans lequel est recherche le min
     x,v,w : points realisant les trois plus petites valeurs 
             de la fonction
     w     : deuxieme du classement derriere x
     v     : precedente valeur de w
     u     : point de la derniere evaluation de la fonction
     d     : deplacement de x pour obtenir u
     xm    : milieu de a et b
     
     tol   : On n'evalue jamais la fonction en y si elle a deja
             ete evaluee en x avec |x-y| < tol
	     Raison d'etre : au voisinage du minimum, la fonction f 
	     est assimilable a une parabole (du moins on l'espere!)
	     Il est donc inutile et couteux en temps de calcul 
	     d'exiger pour l'argmin une precision egale a la precision 
	     en double (1e-15)
	     => Un bon choix pour tol = racine carree de la precision
	     en double, i.e. 3e-8
  --- */

  if ( mystdout == NULL )
    PowellVerboseOnStderr( );

  /* --- Precaution en principe inutile ... */
  d=0;
  
  /* --- On ordonne a et b si necessaire --- */
  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  
  /* --- Initialisations --- */
  x = w = v = bx;
  for ( i=0; i<n; i++ ) aux[i] = pt[i] + x * di[i];
  fx = (*func)( aux, par );
  fw = fv = fx;

  /* --- Boucle iterative ... --- */
  if ( _VERBOSE_)
    fprintf(mystdout,"\n");

  for (iter=1;iter<=ITMAXBRENT;iter++) {
    xm = 0.5*(a+b);
    tol2 = 2.0 * ( tol1=tol*fabs(x)+ZEPS );
    
    if ( _VERBOSE_ >= 2 ) {
      fprintf(mystdout,"  %d BRENT - it n.", _VERBOSE_);
      fprintf(mystdout,"%3d  ", iter);
    }

    /* --- Test de sortie --- */
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return (fx) ; 
    }

    /* --- 1er cas : l'approximation parabolique entre les poins
       (x,fx), (w,fw), (v,fv) est PEUT-ETRE acceptable. 
       Reste a le verifier plus amplement; si elle ne l'est 
       pas, on calcule le nouveau point par l'antique mais 
       toujours utile methode de la section d'or (due aux 
       pythagoriciens!)
       --- */
    if (fabs(e) > tol1) {
      
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) 
	  || p <= q*(a-x) 
	  || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=Sign(tol1,xm-x);
      }
    }
    
    /* --- 2e cas : on est deja sur de ne pas pouvoir 
       effectuer une approximation parabolique; on va 
       calculer le nouveau point par la methode de la 
       section d'or
       --- */
    else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }


    u=(fabs(d) >= tol1 ? x+d : x+Sign(tol1,d));
    /* fu=(*f)(u); */
    for ( i=0; i<n; i++ ) aux[i] = pt[i] + u * di[i];
    fu = (*func)( aux, par );

		
    /* --- Si u est "meilleur" que x, on construit un nouveau
       triplet (v,w,x) en abandonnant la valeur actuelle
       de v. Par ailleurs on restreint l'intervalle de
       recherche [a,b] --- */
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      Shft(v,w,x,u);
      Shft(fv,fw,fx,fu);
    }

    /* --- Sinon, on teste si u est "meilleur" que w puis, dans
       la negative, si u est "meilleur" que v. 
       Meme dans le pire cas, ie fu > fv, il se passe
       au moins quelquechose car on restreint l'intervalle
       de recherche [a,b] --- */
    else {
      if (u < x) a=u; else b=u;
      
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      }
      
      else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
      
    }
    
  }
  
  /* --- Le programme arrive ici si il a atteint le nombre d'iterations
     maximal sans satisfaire au test de sortie... 
     --- */
  
  if ( _VERBOSE_)
    fprintf( mystdout, "brent: Too many iterations in Brent - exit\n" );
  *xmin = x;
  return (fx);
}










/* -----------------------------------------------------------------------
     MNBRAK : routine pour encadrer le minimum d'une fonction a variable
     reelle. Indispensable avant toute recherche de minimum.
     *ax et *bx sont les 2 bornes d'un intervalle

     on renvoie
     *ax : 1ere borne
     *bx : minimum estime
     *cx : 2eme borne
     
     
   ----------------------------------------------------------------------- */

void Mnbrak(double *ax, 
	    double *bx, 
	    double *cx, 
	    double *pt, /* le point */
	    double *di, /* la direction */
	    double *aux,
	    int n,
	    double (*func)(double *, void *),
	    void *par )
{
  int i;
  double fa;
  double fb;
  double fc; 
  double fu;
  double tmp;
  
  double ulim, u, r, q; 


  /*
  fa=(*func)(*ax);
  fb=(*func)(*bx);
  */

  /* on ordonne a et b d'apres la valeur de la fonction;
     on aura fb <= fa
   */
  for ( i=0; i<n; i++ ) aux[i] = pt[i] + (*ax) * di[i];
  fa = (*func)( aux, par );
  for ( i=0; i<n; i++ ) aux[i] = pt[i] + (*bx) * di[i];
  fb = (*func)( aux, par );

  if (fb > fa) {
    tmp = *ax; *ax = *bx; *bx = tmp;
    tmp = fa;  fa = fb;   fb = tmp;
  }

  /* on cherche c : estimee initiale
     c est de l'autre cote de b par rapport a a
   */
  *cx = (*bx) + GOLD * (*bx-*ax);
  /* fc=(*func)(*cx); */
  for ( i=0; i<n; i++ ) aux[i] = pt[i] + (*cx) * di[i];
  fc = (*func)( aux, par );


  /* on cherche encore 
     on veut fc plus petite que fb 
   */
  while (fb > fc) {

    /* u -> parabolic extrapolation de la fonction
       en utilisant a, b et c */
    r = (*bx-*ax) * (fb-fc);
    q = (*bx-*cx) * (fb-fa);
    u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*Sign(Max(fabs(q-r),TINY),q-r));

    ulim=(*bx)+GLIMIT*(*cx-*bx);


    if ( (*bx-u)*(u-*cx) > 0.0) {
      /* u est entre b et c */
      /* fu=(*func)(u); */
      for ( i=0; i<n; i++ ) aux[i] = pt[i] + u * di[i];
      fu = (*func)( aux, par );

      if (fu < fc) {
	/* minimum u entre b et c :
	   a <- b
	   b <- u 
	*/
	*ax = (*bx);
	*bx = u; 
	/* fa = fb; 
	   fb = fu;
	*/
	return;

      } else if (fu > fb) {
	/* minimum b entre a et u
	   c <- u
	 */
	*cx = u;
	/* fc = fu; */
	return;
      }
      
      /* ici (fb >= fu >= fc)
	 les points sont tjourous ordonnes
	 on prend u de l'autre cote de c par rapport a b
      */
      u=(*cx)+GOLD*(*cx-*bx);
      /* fu=(*func)(u); */
      for ( i=0; i<n; i++ ) aux[i] = pt[i] + u * di[i];
      fu = (*func)( aux, par );
      
    } /* fin de if ( (*bx-u)*(u-*cx) > 0.0) */


    else if ((*cx-u)*(u-ulim) > 0.0) {
      /* u est entre c et la limite pour u */
      /* fu=(*func)(u); */
      for ( i=0; i<n; i++ ) aux[i] = pt[i] + u * di[i];
      fu = (*func)( aux, par );

      if (fu < fc) {
	/* u est apres c et fu < fc
	   b <- c
	   c <- u
	   u <- c + GOLD* (c-b)
	*/
	*bx = *cx;   *cx = u;   u = *cx+GOLD*(*cx-*bx);
	fb = fc;     fc = fu;
	for ( i=0; i<n; i++ ) aux[i] = pt[i] + u * di[i];
	fu = (*func)( aux, par );
      }
    } 
    
    else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      /* la limite pour u est entre u et c
	 u <- la limite pour u
      */
      u = ulim;
      for ( i=0; i<n; i++ ) aux[i] = pt[i] + u * di[i];
      fu = (*func)( aux, par );
    } 

    else {
      u = (*cx)+GOLD*(*cx-*bx);
      /* fu=(*func)(u); */
      for ( i=0; i<n; i++ ) aux[i] = pt[i] + u * di[i];
      fu = (*func)( aux, par );
    }
    
    /* a <- b
       b <- c
       c <- u
    */
    *ax = *bx;   *bx = *cx;   *cx = u;
    fa = fb;     fb = fc;     fc = fu;

  } /* while ( fb > fc ) */
}








/* ------------------------------------------------------------------------
      LINMIN : minimise une fonction dans une direction donnee.
      Appelle MNBRAK, BRENT et FUNC1D.
   ------------------------------------------------------------------------ */

void Linmin (double *p, /* le point */
	     double *xi, /* la direction */
	     double *aux, 
	     int n, /* nombre de parametres */
	     double *fret, 
	     double tol,
	     double (*func)(double *, void *),
	     void *par )
{
  int j;
  double xx, xmin, bx, ax;


  if ( mystdout == NULL )
    PowellVerboseOnStderr( );


  /* valeurs initiales */
  ax = INITAX;
  xx = INITXX;
  
  if ( _VERBOSE_ >= 3 )
    fprintf(mystdout, "\n  BRAKETTING\n");

  Mnbrak(&ax, &xx, &bx, p, xi, aux, n, func, par );
  /* le point (p + xx * xi)
     est un minimum entre (p + ax * xi) et (p + bx * xi) */

  *fret = Brent( ax, xx, bx, tol, &xmin, 
		 p, xi, aux, n, func, par );
  
  /* mise a jour */
  for (j=0; j<n; j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
}











/* ------------------------------------------------------------------------

      POWELL : fonction calculant le MINIMUM d'une fonction de plusieurs 
      variables en utilisant la methode des directions de Powell
      alpha :     vecteurs de parametres
      n :         dimension de l'espace de recherche (nb de parametres)
      tol, ftol : tolerances pour l'algo de Powell
      iter :      nb d'iterations qui seront faites
      measure :   valeur finale courante du critere
      func : pointeur sur la fonction de calcul de similarite 
      par : les parametres pour la minimisation

   ------------------------------------------------------------------------ */

void Powell (double *p, 
		 int n, 
		 double tol, 
		 double ftol, 
		 int *iter, 
		 double *fret,
		 double (*func)(double *, void *),
		 void *par )
{
  int i, j;
  int ibig;
  
  double del;
  double fp;
  double fptt;
  double t;
  void *allocatedBuffer = (void*)NULL;
  char *aux;
  double *pt=(double*)NULL;
  double *ptt=(double*)NULL;
  double *auxForLinmin = (double*)NULL;
  double *xit=(double*)NULL;    /* direction courante selon laquelle on minimise */
  double **xi = (double**)NULL; /* tableau des directions pour la minisation */

       

  /* --- Signification de ces variables ...
     
     p       le point a trouver a partir du point de depart
     fret    valeur de func en p
     xi      la matrice de l'ensemble de directions (initialement, 
     on prend en general les n vecteurs de la base canonique)
     
     iter    nombre d'iterations realisees
     ftol    precision souhaitee sur les valeurs de func
     ---- */
  
  
  if ( mystdout == NULL )
    PowellVerboseOnStderr( );


  /*--- Allocations
    allocatedBuffer contient 
    - double pt[n]
    - double ptt[n]
    - double xit[n]
    - double xi[n][n]

    Il est important d'allouer les pointeurs APRES les doubles
    car sur une machine 32 bits comme une station sun, les 
    pointeurs sont donc sur 4 octets alors que les doubles 
    sont sur 8 octets.
    Si on met les pointeurs AVANT des doubles et que 'n' est impair
    alors les frontieres entre doubles ne sont plus sur des 
    multiples de 8 et ... Bus Error.
    (Gregoire Malandain Fri Dec  7 16:42:34 MET 2001)
   */
  allocatedBuffer = (void*)malloc( n * sizeof(double) +
				   n * sizeof(double) +
				   n * sizeof(double) +
				   n * sizeof(double) +
				   n * n * sizeof(double) +
				   n * sizeof(double*) );
  if ( allocatedBuffer == (void*)NULL ) {
    if ( _VERBOSE_ ) 
      fprintf( stderr, "Powell: unable to allocate auxiliary buffer (1)\n" );
    return;
  }
  
  /*--- initialisation des pointeurs
    mise a zero de xi[][] par memset
    et initialisation a la base canonique
  */
  pt = (double*)allocatedBuffer;
  ptt = (double*)allocatedBuffer; ptt += n;
  xit = (double*)allocatedBuffer; xit += 2*n;
  auxForLinmin = (double*)allocatedBuffer; auxForLinmin += 3*n;
					     
  aux = (char*)allocatedBuffer;   aux += 4*n * sizeof(double);
  xi = (double**)(aux + n * n * sizeof(double) );

  (void)memset( (void*)aux, 0, n * n* sizeof(double) );
  for ( i=0; i<n; i++, aux+=n * sizeof(double) ) {
    xi[i] = (double*)aux;
    xi[i][i] = 1;
  }
	
  /* --- Calcul du critere au point initial --- */
  /* *fret=(*func)(p); */
  *fret = (*func)( p, par );

	
  /* --- Sauvegarde du point initial --- */
  for (j=0; j<n; j++) pt[j] = p[j];





  /* --- Debut de la boucle --- */
  for (*iter=1 ; ; ++(*iter) ) {
	        
    fp = (*fret);
    ibig = 0;
    del = 0.0;



    /* --- On minimise dans chaque direction --- */
    for (i=0; i<n; i++) {
      
      for (j=0; j<n; j++) xit[j]  = xi[j][i];
      
      if ( _VERBOSE_ ) {
	fprintf(mystdout,"\nPOWELL (Iteration %d - Direction %d)", *iter, i+1);
	if ( _VERBOSE_ > 1 ) {
	  fprintf(mystdout," [ DIR = ");
	  for (j=0; j<n; j++) fprintf(mystdout," %5.3f", xit[j] );
	  fprintf(mystdout," ] ");
	}
      }
      


      fptt = (*fret);
      Linmin( p, xit, auxForLinmin, n, fret, tol, func, par );

      if ( _VERBOSE_ ) {
	fprintf(mystdout,"[VAL : %lf -> %lf]\n", fptt, *fret );
	fprintf(mystdout,"\n" );
      }

      /* est-ce la direction de plus grande variation
	 (positive ou negative) ?
      */
	 
      if (fabs(fptt-(*fret)) > del) {
	del=fabs(fptt-(*fret));
	ibig=i;
      }
    } /* fin de minimisation selon chaque direction */ 

    
    
    /* --- Test d'arret --- */
    if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {

      free( allocatedBuffer );
      return;
    }
    


    /* --- Test d'erreur pour cause de trop grand nombre
       d'iterations --- */
    if (*iter == ITMAXPOWELL) {
      if ( _VERBOSE_)
	fprintf( mystdout, "powell: Too many iterations in Powell - exit\n" );
      free( allocatedBuffer );
      return;
    }
    

    
    /* --- Calcule le point extrapole et la direction moyenne
       dans laquelle on a evolue. 
       Sauvegarde l'ancien point de depart --- */
    for (j=0; j<n; j++) {
      ptt[j]=2.0*p[j]-pt[j];
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
    }
    
    /* --- Calcule le minimum dans la nouvelle direction sous
       certaines conditions et sauvegarde la nouvelle
       direction --- */
    
    /*  fptt=(*func)(ptt); */
    fptt = (*func)( ptt, par );

    if (fptt < fp) {
      t = 2.0*(fp-2.0*(*fret)+fptt)*Sqr(fp-(*fret)-del)
	-del*Sqr(fp-fptt);
      if (t < 0.0) {
	
	if ( _VERBOSE_)
	  fprintf(mystdout,"\n\t New direction ");
	
	/* linmin(p,xit,n,fret,func,tol); */
	Linmin( p, xit, auxForLinmin, n, fret, tol, func, par );

	for (j=0; j<n; j++) {
	  xi[j][ibig]=xit[j];
	  
	  /*   xi[j][ibig]=xi[j][n-1]; 
	       xi[j][n-1]=xit[j];
	  */   
	}
      }
    }

  } /* fin de la boucle des iterations */
  
  
  free ( allocatedBuffer );
  return;
}


/* Descente de gradient du 1er ordre avec regle de Cauchy pour calculer le pas de descente. Cet algorithme 
   revient a une minimisation par lignes dans la direction du gradient.
   Le buffer scl sert a mettre les parametres a la meme echelle. */
void SteepestDescent (double *p, 
	     int n, 
	     double tol, 
	     double ftol, 
	     int *iter, 
	     double *fret,
	     double (*func)(double *, void *),
	     void (*dfunc)(double *, double *, void *),
	     double * scl,
	     void *par )

{

  int i;
  /* double epsilon, aux, auxcut, step; */
  double fp, normg, invsqrt_n, epsilon;
  double *pt = (double*)NULL; 
  double *g = (double*)NULL; 
  double *xit = (double*)NULL;
  double *auxForLinmin = (double*)NULL;

  /* Allocations */
  pt  = (double *) calloc ( n, sizeof (double) );
  g   = (double *) calloc ( n, sizeof (double) );
  xit = (double *) calloc ( n, sizeof (double) );
  auxForLinmin = (double *) calloc ( n, sizeof (double) );

  /* Initialisation */
  invsqrt_n = 1.0/sqrt((double)n);
  *fret = (*func)(p, par);
  (*dfunc)(p, g, par);
   

  /* Debut de la boucle */
  for (*iter=1; *iter <= ITMAXCAUCHY; ++(*iter) ) {

    fp = (*fret);
    
    /* Baratin */
    fprintf(stderr,"\nGRADIENT DESCENT (Iteration %d)\n", *iter);
    fprintf(stderr,"param = [");
    for (i=0; i<n; i++) fprintf(stderr," %5.3f", p[i] );
    fprintf(stderr, "  ]\n" );
    fprintf(stderr, "grad  = [" );
    for ( i = 0; i < n; i ++ ) fprintf(stderr, "  %f", g[i] );
    fprintf(stderr, "  ]\n" );
    
    /* Calcul de la direction de descente. L'idee est de descendre dans une direction:
       - S*g, ou S est une matrice diagonale definie positive. */
    normg = 0.0;
    for (i=0; i<n; i++)
      normg += Sqr(scl[i] * g[i]);
    normg = sqrt(normg);
    
    /* Si ||g||=0, on se deplace dans une direction fixee */
    if ( normg == 0.0 ) {
      break;
      /* for (i=0; i<n; i++)
	 xit[i] = invsqrt_n; */
    } 
    else {
      for (i=0; i<n; i++)
	xit[i] = - scl[i]*g[i]/normg;
    }
    
    /* On stocke la position courante */
    for (i=0; i<n; i++) 
      pt[i]=p[i];

    /* On cherche un minimum local du critere dans la direction
       du gradient. La nouvelle position est:
       p <- (p - x g)
    */
    Linmin( p, xit, auxForLinmin, n, fret, tol, func, par );
    
      
    /* Calcul de la variation en position */
    epsilon = 0;
    for ( i = 0; i < n; i ++ )
      epsilon = Max ( epsilon, Abs(p[i]-pt[i]) );

    /* Test d'arret */
    /*    if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) */
    if (epsilon <= ftol)
      break;
    
    /* Mise a jour du gradient */
    (*dfunc)(p, g, par);
    
  } /* fin de la boucle iterative */
  
  /* Si on arrive ici, c'est p-e qu'il y a eu trop d'iterations */
  if ( *iter > ITMAXCAUCHY ) {
    if ( _VERBOSE_)
      fprintf( stderr, "powell: Too many iterations in Powell - exit\n" );
  }

  /* Liberation memoire */
  /*  free(pt); 
      free(g);
      free(xit);
      free(auxForLinmin);
  */
  return;
}

/* Methode du gradient conjugue. Les directions successives sont construites par la
   methode de Polak-Ribiere */

void ConjugateGradient (double *p, 
			int n, 
			double tol, 
			double ftol, 
			int *iter, 
			double *fret,
			double (*func)(double *, void *),
			void (*dfunc)(double *, double *, void *),
			double * scl,
			void *par )

{

  int i;
  /* double epsilon, aux, auxcut, step; */
  double fp, normg, invsqrt_n, epsilon, gamma;
  double *pt = (double*)NULL; 
  double *g = (double*)NULL; 
  double *dg = (double*)NULL; 
  double *xit = (double*)NULL;
  double *auxForLinmin = (double*)NULL;

  /* Allocations */
  pt  = (double *) calloc ( n, sizeof (double) );
  g   = (double *) calloc ( n, sizeof (double) );
  dg  = (double *) calloc ( n, sizeof (double) );
  xit = (double *) calloc ( n, sizeof (double) );
  auxForLinmin = (double *) calloc ( n, sizeof (double) );

  /* Initialisation */
  invsqrt_n = 1.0/sqrt((double)n);
  *fret = (*func)(p, par);
  (*dfunc)(p, g, par);
  for (i=0; i<n; i++)
    xit[i] = - g[i];

  /* Debut de la boucle */
  for (*iter=1; *iter <= ITMAXCAUCHY; ++(*iter) ) {

    fp = (*fret);
    
    /* Baratin */
    fprintf(stderr,"\nGRADIENT DESCENT (Iteration %d)\n", *iter);
    fprintf(stderr,"param = [");
    for (i=0; i<n; i++) fprintf(stderr," %5.3f", p[i] );
    fprintf(stderr, "  ]\n" );
    fprintf(stderr, "grad  = [" );
    for ( i = 0; i < n; i ++ ) fprintf(stderr, "  %f", g[i] );
    fprintf(stderr, "  ]\n" );
    fprintf(stderr, "dir   = [" );
    for ( i = 0; i < n; i ++ ) fprintf(stderr, "  %f", xit[i] );
    fprintf(stderr, "  ]\n" );
    
    /* On stocke la position courante et le gradient courant */
    for (i=0; i<n; i++) {
      pt[i] = p[i];
      dg[i] = -g[i]; 
    }

    /* On cherche un minimum local du critere dans la direction
       du gradient. La nouvelle position est:
       p <- (p - x g)
    */
    Linmin( p, xit, auxForLinmin, n, fret, tol, func, par );
    
      
    /* Calcul de la variation en position */
    epsilon = 0.0;
    for ( i = 0; i < n; i ++ )
      epsilon = Max ( epsilon, Abs(p[i]-pt[i]) );

    /* Test d'arret */
    /*    if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) */
    if (epsilon <= ftol)
      break;
    
    /* Mise a jour du gradient */
    (*dfunc)(p, g, par);

     /* Variation du gradient et calcul du gamma */
    gamma = 0.0;
    normg = 0.0;
    for (i=0; i<n; i++) {
      normg += Sqr(dg[i]);
      dg[i] += g[i];     
      gamma += dg[i]*g[i];
      /* gamma += Sqr(g[i]); Fletcher-Reeves */
    }
    if (normg == 0.0)
      break;
    gamma = gamma / normg;

    /* Mise a jour de la direction de descente */
    for (i=0; i<n; i++)
      xit[i] = gamma*xit[i] - g[i];
    
  } /* fin de la boucle iterative */
  
  /* Si on arrive ici, c'est p-e qu'il y a eu trop d'iterations */
  if ( *iter > ITMAXCAUCHY ) {
    if ( _VERBOSE_)
      fprintf( stderr, "powell: Too many iterations in Powell - exit\n" );
  }

  /* Liberation memoire */
  /*  free(pt); 
      free(dg);
      free(g);
      free(xit);
      free(auxForLinmin);
  */
  return;
}



/* Methode de QuasiNewton : minimisation par lignes dans la direction H*g, 
   ou g est le gradient du critere et H une approximation du hessien inverse
   (on utilise la formule BFGS).

   Le buffer scl sert a mettre les parametres a la meme echelle, et par 
   la meme occasion, a initialiser H. */

void QuasiNewton (double *p, 
	     int n, 
	     double tol, 
	     double ftol, 
	     int *iter, 
	     double *fret,
	     double (*func)(double *, void *),
	     void (*dfunc)(double *, double *, void *),
	     double * scl,
	     void *par )

{

  int i,j;
  double fp, invsqrt_n, epsilon;
  double *dp  = (double*)NULL; 
  double *g   = (double*)NULL; 
  double *dg  = (double*)NULL; 
  double *xit = (double*)NULL;
  double *auxForLinmin = (double*)NULL;
  double ** H;
  /* double trace_H; */

  /* Variables pour le calcul du pseudo-hessien */
  double dpdg, dp_Hdgdg;
  double * dp_Hdg;
  double ** Haux1, ** Haux2, ** Haux3;
  double c1, c2;

   /* ---
     Signification de ces variables :
     g        -> gradient du critere
     dg       -> variation du gradient
     H        -> estimee courante du hessien inverse
     --- */  

  /* Allocations */
  dp  = (double *) calloc ( n, sizeof (double) );
  g   = (double *) calloc ( n, sizeof (double) );
  dg  = (double *) calloc ( n, sizeof (double) );
  xit = (double *) calloc ( n, sizeof (double) );
  auxForLinmin = (double *) calloc ( n, sizeof (double) );
  H       = _Powell_Matrix ( n, n );
  
  dp_Hdg  = (double *) calloc ( n, sizeof(double) );
  Haux1   = _Powell_Matrix ( n, n );
  Haux2   = _Powell_Matrix ( n, n );
  Haux3   = _Powell_Matrix ( n, n );   



  /* Initialisation */
  invsqrt_n = 1.0/sqrt((double)n);
  *fret = (*func)(p, par);
  (*dfunc)(p, g, par);
   
  for ( i = 0; i < n; i ++ )
    for ( j = 0; j < n; j++ ) {
      if ( i == j ) 
	H [i][j] = scl[i];
      else 
	H [i][j] = 0.0;
    }
  
  /* Debut de la boucle */
  for (*iter=1; *iter <= ITMAXCAUCHY; ++(*iter) ) {
    
    fp = (*fret);
    
    /* Baratin */
    fprintf(stderr,"\nGRADIENT DESCENT (Iteration %d)\n", *iter);
    fprintf(stderr,"param = [");
    for (i=0; i<n; i++) fprintf(stderr," %5.3f", p[i] );
    fprintf(stderr, "  ]\n" );
    fprintf(stderr, "grad  = [" );
    for ( i = 0; i < n; i ++ ) fprintf(stderr, "  %f", g[i] );
    fprintf(stderr, "  ]\n" );
    
    /* Calcul de la direction de descente.
       L'idee est de descendre dans la direction H*g.
    */
    for ( i=0; i<n; i ++ ) {
      xit[i] = 0.0;
      for ( j=0; j<n; j ++ )
	xit[i] -= H[i][j]*g[j];
    }
    
    /* On stocke la position courante et le gradient courant */
    for (i=0; i<n; i++) {
      dp[i] = -p[i];
      dg[i] = -g[i]; 
    }
    /* On cherche un minimum local du critere dans la direction
       choisie. La nouvelle position est:
       p <- (p - x xit)
    */
    Linmin( p, xit, auxForLinmin, n, fret, tol, func, par );
   
      
    /* Calcul de la variation en position */
    epsilon = 0;
    for ( i = 0; i < n; i ++ ) {
      dp[i] += p[i];
      epsilon = Max ( epsilon, Abs(dp[i]) );
    }

    /* Test d'arret */
    if (epsilon <= ftol)
      break;
    
    /* Mise a jour du gradient */
    (*dfunc)(p, g, par);
    
    /* Mise a jour du hessien inverse. On utilise la formule de BFGS, voir
       Numerical Recipes. */

    /* Variation du gradient */
    for (i=0; i<n; i++) 
      dg[i] += g[i]; 
    

    /* Calcul du vecteur (dp - H dg)  et des produits scalaires:
       dp^t dg, et (dp - H dg)^t dg */
    _Powell_MultMatVect ( H, dg, dp_Hdg, n, n ); 

    dpdg = 0;
    dp_Hdgdg = 0;
    for ( i = 0; i < n; i ++ ) {
      dp_Hdg [i] = dp [i] - dp_Hdg [i];
      
      dpdg += dp [i] * dg [i]; 
      dp_Hdgdg += dp_Hdg [i] * dg [i];
    }
    
    /* On met a jour le hessien si dpdg est suffisamment grand */
    if ( fabs (dpdg) > 1e-5 ) {
      
      _Powell_OuterProduct ( dp_Hdg, dp, Haux1, n ); 
      _Powell_OuterProduct ( dp, dp_Hdg, Haux2, n ); 
      _Powell_OuterProduct ( dp, dp, Haux3, n );     
      
      /*  On a alors la formule:
	  H = H 
	  + (1/dpdg)*Haux1 
	  + (1/dpdg)*Haux2 
	  + ( (dp_Hdgdg)/(dpdg^2) )*Haux3 */
      
      c1 = 1.0 / dpdg;
      c2 = dp_Hdgdg * c1 * c1; 
      
      for ( i = 0; i < n; i ++ )
	for ( j = 0; j < n; j ++ ) 
	  H[i][j] = H[i][j] + c1*Haux1[i][j] + c1*Haux2[i][j] - c2*Haux3[i][j]; 
    }
    
  } /* fin de la boucle iterative */
  
  /* Si on arrive ici, c'est p-e qu'il y a eu trop d'iterations */
  if ( *iter > ITMAXCAUCHY ) {
    if ( _VERBOSE_)
      fprintf( stderr, "powell: Too many iterations in Powell - exit\n" );
  }

  
  fprintf( stderr, "\n\nHessien = [ \n" );
  for ( i = 0; i < n; i ++ ) {
    for ( j = 0; j < n; j ++ ) 
      fprintf(stderr, "%f\t", H[i][j] );
     fprintf(stderr, "\n");
    }
  fprintf(stderr, "]\n");

  
  /* Liberation memoire */
  /*  free(pt); 
      free(g);
      free(dg);
      free(xit);
      free(auxForLinmin);
  */
  return;
}



/* Allocation d'une matrice p x q */

double ** _Powell_Matrix ( int p, int q )
{

  int i;
  double ** H;

  H = (double **) calloc ( p, sizeof(double *) );
  for ( i = 0; i < p; i ++ )
    H [i] = (double *) calloc ( q, sizeof(double) );

  return H;
}

/* Liberation d'une matrice allouee par _Powell_Matrix */

void _Powell_FreeMatrix ( double ** H, int p, int q )
{
  int i;

  for ( i = 0; i < p; i ++ )
    free ( H [i] );

  free (H);
  
  return;
}

/* Multiplication d'un vecteur par une matrice carree */

void _Powell_MultMatVect ( double ** A, double * x, double * ax, int p, int q )
{

  int i, j;
  double aux;

  for ( i = 0; i < p; i ++ ) {
    aux = 0;
    for ( j = 0; j < q; j ++ ) {
      aux += A [i][j] * x [j];
    }
    ax[i] = aux;
  }
  
  return;
}



/* Outer product of two vectors */

void _Powell_OuterProduct ( double * x, double * y, double ** A, int p )
{

  int i, j;
  
  for ( i = 0; i < p; i ++ ) 
    for ( j = 0; j < p; j ++ ) 
      A [i][j] = x [i] * y [j];
  
  return;
}




void PowellVerboseON()
{
  if ( _VERBOSE_ <= 0 ) _VERBOSE_ = 1;
  else _VERBOSE_ ++;
}

void PowellSetVerbose( int v )
{
  _VERBOSE_ = v;
}

void PowellVerboseOFF()
{
  _VERBOSE_ = 0;
}

void PowellVerboseOnStdout( )
{
  mystdout = stdout;
}

void PowellVerboseOnStderr( )
{
  mystdout = stderr;
}



