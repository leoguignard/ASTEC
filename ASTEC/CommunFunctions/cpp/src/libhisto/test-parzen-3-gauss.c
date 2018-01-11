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

  double xmin = mu - 4.0 * sigma;
  double xmax = mu + 4.0 * sigma;
  double x, dx = 0.01;
  double bin1 = 1.0;
  double bin2 = 0.5;
  double *y1, *y2;
  double *z1, *z2;
  double *x1, *x2;
  int nx1, nx2;
  int nx = (xmax-xmin)/dx + 1;
  double *py1;
  double *py2;

  double parzen_sigma = 0.5;

  double ds = 0.10;
  double t, c, s;
  
  int i, j, k;
  double sample[NB];
  
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

  
  for (i=0; i<NB; i++ ) sample[i] = Rnd_normal( mu, sigma );
   write( fd, sample, NB * sizeof(double) );
  fprintf( fp, "\n" );
   fprintf( fp, "[SAMPLES, NS] = fread( fid, %d, 'float%d' );\n", NB, 8*sizeof(double) );

  nx1 = (int)(20 / bin1)+1;
  nx2 = (int)(20 / bin2)+1;

  x1 = (double*)malloc( nx1*sizeof(double) );
  y1 = (double*)malloc( nx1*sizeof(double) );
  z1 = (double*)malloc( nx1*sizeof(double) );
  x2 = (double*)malloc( nx2*sizeof(double) );
  y2 = (double*)malloc( nx2*sizeof(double) );
  z2 = (double*)malloc( nx2*sizeof(double) );

  for (k=0;k<nx1;k++) {
    x1[k] = bin1/2.0 + k * bin1;
    y1[k] = 0.0;
  }
  for (k=0;k<nx2;k++) {
    x2[k] = bin2/2.0 + k * bin2;
    y2[k] = 0.0;
  }

  for (i=0; i<NB; i++ ) {
    y1[(int)(sample[i]/bin1)] ++;
    y2[(int)(sample[i]/bin2)] ++;
  }

  for (k=0;k<nx1;k++) {
    y1[k]/= NB;
  }
  for (k=0;k<nx2;k++) {
    y2[k]/= NB;
  }








  write( fd, x1, nx1 * sizeof(double) );
  write( fd, y1, nx1 * sizeof(double) );
  fprintf( fp, "\n" );
  fprintf( fp, "[X1, NX1] = fread( fid, %d, 'float%d' );\n", nx1, 8*sizeof(double) );
  fprintf( fp, "[Y1, NY1] = fread( fid, %d, 'float%d' );\n", nx1, 8*sizeof(double) );
  
  write( fd, x2, nx2 * sizeof(double) );
  write( fd, y2, nx2 * sizeof(double) );
  fprintf( fp, "\n" );
  fprintf( fp, "[X2, NX2] = fread( fid, %d, 'float%d' );\n", nx2, 8*sizeof(double) );
  fprintf( fp, "[Y2, NY2] = fread( fid, %d, 'float%d' );\n", nx2, 8*sizeof(double) );
  fprintf( fp, "\n" );




  py1 = (double*)malloc( nx*sizeof(double) );
  py2 = (double*)malloc( nx*sizeof(double) );

  


  for (k=0;k<nx;k++) py1[k] = py2[k] = 0.0;

  c = 1.0 / ( parzen_sigma * sqrt( 2.0 * 3.1415926536 ) );
  s = 1.0 / ( 2.0 * parzen_sigma * parzen_sigma );

  for (k=0;k<nx;k++) {
    x = xmin + k * dx;
    for (i=0; i<nx1; i++ ) {
      t = x1[i] - x;
      py1[k] += y1[i] * c * exp( - t * t  * s );
    }
  }

  for (k=0;k<nx;k++) {
    x = xmin + k * dx;
    for (i=0; i<nx2; i++ ) {
      t = x2[i] - x;
      py2[k] += y2[i] * c * exp( - t * t  * s );
    }
  }


  write( fd, py1, nx * sizeof(double) );
  fprintf( fp, "\n" );
  fprintf( fp, "[PY1, NPY1] = fread( fid, %d, 'float%d' );\n", nx, 8*sizeof(double) );
  write( fd, py2, nx * sizeof(double) );
  fprintf( fp, "[PY2, NPY2] = fread( fid, %d, 'float%d' );\n", nx, 8*sizeof(double) );
  fprintf( fp, "\n" );


  c = 1.0 / ( parzen_sigma * sqrt( 2.0 ) );

  for (k=0;k<nx1;k++) z1[k] = 0.0;
  for (k=0;k<nx1;k++) {
    for (i=0; i<nx2; i++ ) {
      z1[k] += y2[i] * 0.5 * ( erf( c * (x1[k]+bin1/2.0 - x2[i] ) ) - erf( c * (x1[k]-bin1/2.0 - x2[i]) ) );
    }
  }

  for (k=0;k<nx2;k++) z2[k] = 0.0;
  for (k=0;k<nx2;k++) {
    for (i=0; i<nx1; i++ ) {
      z2[k] += y1[i] * 0.5 * ( erf( c * (x2[k]+bin2/2.0 - x1[i] ) ) - erf( c * (x2[k]-bin2/2.0 - x1[i]) ) );
    }
  }

  write( fd, z1, nx1 * sizeof(double) );
  fprintf( fp, "\n" );
  fprintf( fp, "[Z1, NZ1] = fread( fid, %d, 'float%d' );\n", nx1, 8*sizeof(double) );
  write( fd, z2, nx2 * sizeof(double) );
  fprintf( fp, "[Z2, NZ2] = fread( fid, %d, 'float%d' );\n", nx2, 8*sizeof(double) );
  fprintf( fp, "\n" );


  fprintf( fp, "X=[%f:%f:%f];\n",xmin, dx, xmax );
  fprintf( fp, "GX=1/(sqrt(2*pi)*%f) * exp( - (X-%f) .* (X-%f) / (2.0*%f*%f) );\n",
	   sigma, mu, mu, sigma, sigma );
  
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );

  fprintf( fp, "for i = 1:%d\n", NB );
  fprintf( fp, "  plot( [SAMPLES(i),SAMPLES(i)], [0,0.05], 'k' );\n" );
  fprintf( fp, "end\n" );
  fprintf( fp, "stem(X1, Y1);\n" );
  fprintf( fp, "plot(X,GX,'k','Linewidth',2);\n" );
  fprintf( fp, "axis([5 15 0 0.5]);\n" );
  fprintf( fp, "hold off;\n" );

  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );

  fprintf( fp, "for i = 1:%d\n", NB );
  fprintf( fp, "  plot( [SAMPLES(i),SAMPLES(i)], [0,0.05], 'k' );\n" );
  fprintf( fp, "end\n" );
  fprintf( fp, "stem(X2, Y2);\n" );
  fprintf( fp, "plot(X,GX,'k','Linewidth',2);\n" );
  fprintf( fp, "axis([5 15 0 0.5]);\n" );
  fprintf( fp, "hold off;\n" );

  fprintf( fp, "\n" );
  fprintf( fp, "\n" );

  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );

  fprintf( fp, "stem(X1, Y1);\n" );
  fprintf( fp, "plot(X,GX,'k','Linewidth',1);\n" );
  fprintf( fp, "plot(X,PY1, 'Linewidth',2);\n" );
  fprintf( fp, "plot(X2, Z2, 'sr-');\n" );
  fprintf( fp, "axis([5 15 0 0.5]);\n" );
  fprintf( fp, "hold off;\n" );

  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );

  fprintf( fp, "stem(X2, Y2);\n" );
  fprintf( fp, "plot(X,GX,'k','Linewidth',1);\n" );
  fprintf( fp, "plot(X,PY2, 'Linewidth',2);\n" );
  fprintf( fp, "plot(X1, Z1, 'sr-');\n" );
  fprintf( fp, "axis([5 15 0 0.5]);\n" );
  fprintf( fp, "hold off;\n" );

  fprintf( fp, "\n" );
  fprintf( fp, "\n" );

  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );

  fprintf( fp, "stem(X1, Y1);\n" );
  fprintf( fp, "plot(X1, Z1, 'sr-', 'Linewidth',2);\n" );
  fprintf( fp, "axis([5 15 0 0.5]);\n" );
  fprintf( fp, "hold off;\n" );

  fprintf( fp, "\n" );
  fprintf( fp, "\n" );

  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );

  fprintf( fp, "stem(X2, Y2);\n" );
  fprintf( fp, "plot(X2, Z2, 'sr-', 'Linewidth',2);\n" );
  fprintf( fp, "axis([5 15 0 0.5]);\n" );
  fprintf( fp, "hold off;\n" );

  fprintf( fp, "fclose(fid);\n" );




  fclose( fp );
  close( fd );
}
  
