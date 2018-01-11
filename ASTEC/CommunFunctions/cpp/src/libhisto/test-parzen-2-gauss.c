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







#define NB  100
#define NB1 75
#define NB  500
#define NB1 375

int main( int argc, char *argv[] )
{
  char out[256] = "testmatlab";
  char longname[256];

  int fd = 0;
  FILE *fp = NULL;

  double sigma1 = 1.0;
  double mu1 = 0.0;
  double sigma2 = 0.5;
  double mu2 = 3.0;

  double xmin = mu1 - 4.0 * sigma1;
  double xmax = mu2 + 4.0 * sigma2;
  double x, dx = 0.01;
  double *y;
  int nx = (xmax-xmin)/dx + 1;
  
  double parzen_sigma = 0.20;
  double ds = 0.10;
  double t, c, s;
  
  int i, j, k;
  double sample[NB];
  
  Rnd_init();


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

  
  for (i=0; i<NB1; i++ )  sample[i] = Rnd_normal( mu1, sigma1 );
  for (i=NB1; i<NB; i++ ) sample[i] = Rnd_normal( mu2, sigma2 );
  write( fd, sample, NB * sizeof(double) );
  fprintf( fp, "\n" );
  fprintf( fp, "[SAMPLES, NS] = fread( fid, %d, 'float%d' );\n", NB, 8*sizeof(double) );


  y = (double*)malloc( nx*sizeof(double) );

  for ( j=0; j<4; j++, parzen_sigma += ds ) {

    for (k=0;k<nx;k++) y[k] = 0.0;

    c = 1.0 / ( parzen_sigma * sqrt( 2.0 * 3.1415926536 ) );
    s = 1.0 / ( 2.0 * parzen_sigma * parzen_sigma );
    for (k=0;k<nx;k++) {
      x = xmin + k * dx;
      for (i=0; i<NB; i++ ) {
	t = sample[i] - x;
	y[k] += c/(double)NB * exp( - t * t  * s );
      }
    }
    
    printf( "sigma = %f\n",  parzen_sigma );

    write( fd, y, nx * sizeof(double) );
    fprintf( fp, "\n" );
    fprintf( fp, "[PARZEN%d, NP] = fread( fid, %d, 'float%d' );\n", j, nx, 8*sizeof(double) );
    
  }

    


  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );

  fprintf( fp, "\n" );
  fprintf( fp, "X=[%f:%f:%f];\n",xmin, dx, xmax );

  fprintf( fp, "\n" );
  fprintf( fp, "plot(X,PARZEN0,'r--','Linewidth',2);\n" );
  fprintf( fp, "plot(X,PARZEN1,'g--','Linewidth',2);\n" );
  fprintf( fp, "plot(X,PARZEN2,'c--','Linewidth',2);\n" );
  fprintf( fp, "plot(X,PARZEN3,'m--','Linewidth',2);\n" );

  fprintf( fp, "legend('sigma = 0.2','sigma = 0.3','sigma = 0.4','sigma = 0.5' );\n" );
  fprintf( fp, "for i = 1:%d\n", NB );
  fprintf( fp, "  plot( [SAMPLES(i),SAMPLES(i)], [0,0.1], 'k' );\n" );
  fprintf( fp, "end\n" );

  fprintf( fp, "GX= %f/(sqrt(2*pi)*%f) * exp( - (X-%f) .* (X-%f) / (2.0*%f*%f)) + %f/(sqrt(2*pi)*%f) * exp( - (X-%f) .* (X-%f) / (2.0*%f*%f));\n",
	   (double)NB1/(double)NB, sigma1, mu1, mu1, sigma1, sigma1,
	   (double)(NB-NB1)/(double)NB, sigma2, mu2, mu2, sigma2, sigma2 );
  fprintf( fp, "plot(X,GX,'b','Linewidth',2);\n" );


  fprintf( fp, "axis([-4 5 0 0.4]);\n" );
  fprintf( fp, "hold off;\n" );
  fprintf( fp, "fclose(fid);\n" );
  fprintf( fp, "\n" );
  fclose( fp );
  close( fd );
}
