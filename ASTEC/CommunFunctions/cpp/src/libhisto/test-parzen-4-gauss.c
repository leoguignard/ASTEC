#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

double drand48(void);
void srand48(long int seedval);
double erf(double x);

/*----------------------------------------------------------------------*/
/*        Random Init of random generator                               */ 
/*----------------------------------------------------------------------*/
void Rnd_init()
{
    srand48( (long) time( NULL) );
}

/*----------------------------------------------------------------------*/
/*       1D gaussian law of variance 1 and mean 0                       */
/*----------------------------------------------------------------------*/
double Rnd_normal0()
{
    register double u1,u2;

    u1 = drand48();
    u2 = drand48();
    return ( sqrt(-2.0*log(u1))*cos(2* 3.14159265358979323846 *u2));
}

double Rnd_normal( double mu, double sigma )
{
    return ( sigma * Rnd_normal0( ) + mu );
}







#define NB 100

int main( int argc, char *argv[] )
{
  char out[256] = "testmatlab";
  char longname[256];

  int fd = 0;
  FILE *fp = NULL;

  double sigma = 1.0;
  double mu = 10.0;

  double a = 1.2;
  double b = 2.5;

  double xmin = mu - 4.0 * sigma;
  double xmax = a*mu+b + 4.0 * sigma * sqrt(a);

  double x, dx = 0.01;
  double bin = 1.0;

  int  nx0 = (int)(25 / bin)+1;
  int nx = (xmax-xmin)/dx + 1;

  double bin2 = 0.5;
  double *x0, *y1, *y2;
  double *z1, *z2;
  double *x1, *x2;
  double *py1;
  double *py2;

  double parzen_sigma = 0.5;

  double ds = 0.10;
  double t, c, s;
  
  int i, j, k;
  double sample1[NB];
  double sample2[NB];
  
  long int init = time(NULL);

  printf("init = %d\n", init );
  srand48( init );

  sprintf( longname, "%s.raw", out );
  fd = creat( longname, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  sprintf( longname, "%s.m", out );
  fp = fopen( longname, "w" );
  
  fprintf( fp, "\n" );
  k = strlen( out );
  for ( i = k-1; i >= 0 && out[i] != '/' ; i-- )
      ;
  fprintf( fp, "fid = fopen('%s.raw', 'r' );\n", &(out[i+1]) );
  fprintf( fp, "\n" );

  
  for (i=0; i<NB; i++ ) sample1[i] = Rnd_normal( mu, sigma );

  write( fd, sample1, NB * sizeof(double) );
  fprintf( fp, "\n" );
  fprintf( fp, "[SAMPLES1, NS1] = fread( fid, %d, 'float%d' );\n", NB, 8*sizeof(double) );
  
  for (i=0; i<NB; i++ ) sample2[i] = a*sample1[i] + b;
  write( fd, sample2, NB * sizeof(double) );
  fprintf( fp, "[SAMPLES2, NS2] = fread( fid, %d, 'float%d' );\n", NB, 8*sizeof(double) );
  fprintf( fp, "\n" );



 
  x0 = (double*)malloc( nx0*sizeof(double) );
  y1 = (double*)malloc( nx0*sizeof(double) );
  y2 = (double*)malloc( nx0*sizeof(double) );
  for (k=0;k<nx0;k++) {
    x0[k] = bin/2.0 + k * bin;
    y1[k] = y2[k] = 0.0;
  }

  for (i=0; i<NB; i++ ) {
    y1[(int)(sample1[i]/bin)] ++;
    y2[(int)(sample2[i]/bin)] ++;
  }

  for (k=0;k<nx0;k++) {
    y1[k]/= NB;
    y2[k]/= NB;
  }

  write( fd, x0, nx0 * sizeof(double) );
  write( fd, y1, nx0 * sizeof(double) );
  write( fd, y2, nx0 * sizeof(double) );

  fprintf( fp, "\n" );
  fprintf( fp, "[X0, NX0] = fread( fid, %d, 'float%d' );\n", nx0, 8*sizeof(double) );
  fprintf( fp, "[Y1, NY1] = fread( fid, %d, 'float%d' );\n", nx0, 8*sizeof(double) );
  fprintf( fp, "[Y2, NY2] = fread( fid, %d, 'float%d' );\n", nx0, 8*sizeof(double) );
  fprintf( fp, "\n" );
  


  py1 = (double*)malloc( nx*sizeof(double) );
  py2 = (double*)malloc( nx*sizeof(double) );
  for (k=0;k<nx;k++) py1[k] = py2[k] = 0.0;

  c = 1.0 / ( parzen_sigma * sqrt( 2.0 * 3.1415926536 ) );
  s = 1.0 / ( 2.0 * parzen_sigma * parzen_sigma );

  for (k=0;k<nx;k++) {
    x = xmin + k * dx;
    for (i=0; i<nx0; i++ ) {
      t = x0[i] - x;
      py1[k] += y1[i] * c * exp( - t * t  * s );
    }
  }

  for (k=0;k<nx;k++) {
    x = xmin + k * dx;
    for (i=0; i<nx0; i++ ) {
      t = x0[i] - x;
      py2[k] += y2[i] * c * exp( - t * t  * s );
    }
  }

  write( fd, py1, nx * sizeof(double) );
  write( fd, py2, nx * sizeof(double) );
  fprintf( fp, "\n" );
  fprintf( fp, "[PY1, NPY1] = fread( fid, %d, 'float%d' );\n", nx, 8*sizeof(double) );
  fprintf( fp, "[PY2, NPY2] = fread( fid, %d, 'float%d' );\n", nx, 8*sizeof(double) );
  fprintf( fp, "\n" );



  z2 = (double*)malloc( nx0*sizeof(double) );

  for (k=0;k<nx0;k++) z2[k] = 0.0;
  c = 1.0 / ( parzen_sigma * sqrt( 2.0 ) );

  for (k=0;k<nx0;k++) {
    for (i=0; i<nx0; i++ ) {
      z2[k] += y2[i] * 0.5 * ( erf( c * ( a*(x0[k]+bin/2.0)+b - x0[i] ) ) 
			       - erf( c * ( a*(x0[k]-bin/2.0)+b - x0[i]) ) );
      
    }
  }

  z1 = (double*)malloc( nx0*sizeof(double) );

  for (k=0;k<nx0;k++) z1[k] = 0.0;
  c = 1.0 / ( parzen_sigma * sqrt( 2.0 ) );

  for (k=0;k<nx0;k++) {
    for (i=0; i<nx0; i++ ) {
      z1[k] += y1[i] * 0.5 * ( erf( c * ( x0[k]+bin/2.0 - x0[i] ) ) 
			       - erf( c * ( x0[k]-bin/2.0 - x0[i]) ) );
      
    }
  }

  write( fd, z1, nx0 * sizeof(double) );
  write( fd, z2, nx0 * sizeof(double) );
  fprintf( fp, "\n" );
  fprintf( fp, "[Z1, NZ1] = fread( fid, %d, 'float%d' );\n", nx0, 8*sizeof(double) );
  fprintf( fp, "[Z2, NZ2] = fread( fid, %d, 'float%d' );\n", nx0, 8*sizeof(double) );
  fprintf( fp, "\n" );


  fprintf( fp, "X=[%f:%f:%f];\n",xmin, dx, xmax );
  fprintf( fp, "GX1=1/(sqrt(2*pi)*%f) * exp( - (X-%f) .* (X-%f) / (2.0*%f*%f) );\n",
	   sigma, mu, mu, sigma, sigma );
  fprintf( fp, "GX2=1/(sqrt(2*pi)*%f) * exp( - (X-%f) .* (X-%f) / (2.0*%f*%f) );\n",
	   sigma*sqrt(a), a*mu+b, a*mu+b, sigma*sqrt(a), sigma*sqrt(a) );
  
  fprintf( fp, "\n" );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );
  fprintf( fp, "for i = 1:%d\n", NB );
  fprintf( fp, "  plot( [SAMPLES1(i),SAMPLES1(i)], [0,0.05], 'k' );\n" );
  fprintf( fp, "end\n" );
  fprintf( fp, "stem(X0, Y1);\n" );
  fprintf( fp, "plot(X,GX1,'k','Linewidth',2);\n" );
  fprintf( fp, "plot(X,PY1, 'Linewidth',2);\n" );
  fprintf( fp, "axis([5 20 0 0.5]);\n" );
  fprintf( fp, "hold off;\n" );

  fprintf( fp, "\n" );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );
  fprintf( fp, "for i = 1:%d\n", NB );
  fprintf( fp, "  plot( [SAMPLES2(i),SAMPLES2(i)], [0,0.05], 'k' );\n" );
  fprintf( fp, "end\n" );
  fprintf( fp, "stem(X0, Y2);\n" );
  fprintf( fp, "plot(X,GX2,'k','Linewidth',2);\n" );
  fprintf( fp, "plot(X,PY2, 'Linewidth',2);\n" );
  fprintf( fp, "axis([5 20 0 0.5]);\n" );
  fprintf( fp, "hold off;\n" );


  fprintf( fp, "\n" );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );
  fprintf( fp, "stem(X0, Y1);\n" );
  fprintf( fp, "plot(X0,Z2, 'rs-', 'Linewidth',2);\n" );
  fprintf( fp, "plot(X0,Z1, 'ko-', 'Linewidth',1);\n" );
  fprintf( fp, "axis([5 20 0 0.5]);\n" );
  fprintf( fp, "hold off;\n" );



  fclose( fp );
  close( fd );
}
  
