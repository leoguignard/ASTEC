#include <math.h>

/******************************************
 Procedures de manipulation des rotations
*******************************************/
#define TH_TINY 1e-6

void Rvector_To_Rmatrix ( double R[4][4], double * r )
{
  register double n[3];
  register double theta, aux1, aux2, auxnx2, auxny2, auxnz2, theta2;
  int i,j;
  
  
  /* --- Calcul de l'angle (defini a 2*pi pres) ... --- */
  theta2 =  r[0]*r[0] + r[1]*r[1] + r[2]*r[2] ;
  theta = sqrt (theta2);
  
  /*--- ... et de l'axe unitaire --- */
  if (theta > TH_TINY) { 
    {
      n[0]=r[0]/theta;
      n[1]=r[1]/theta;
      n[2]=r[2]/theta;
    }

  /* --- On va calculer la matrice R par la formule de Rodrigues :
          R = Id + sin(theta)*Sn + (1-cos(theta))*Sn^2  
       
	  avec 
	               0  -nz  ny
	       Sn =    nz   0 -nx
	              -ny  nx   0

   --- */

  /* --- Terme identite --- */
    for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      R[i][j]=0;

    R[0][0]=1;
    R[1][1]=1;
    R[2][2]=1;
    
    /* if (theta==0) return;*/
    
    /* --- Ajout du terme en sin(theta) --- */
    aux1=sin(theta);
    
    aux2=aux1*n[2];
    R[0][1] += -aux2;
    R[1][0] +=  aux2;
    
    aux2=aux1*n[1];
    R[0][2] +=  aux2;
    R[2][0] += -aux2;
    
    aux2=aux1*n[0];
    R[1][2] += -aux2;
    R[2][1] += aux2;
    
    
    /* --- Ajout du terme en (1-cos(theta)) 
       
       sachant que 
       
                -ny^2-nz^2   nxny         nxnz
       Sn^2 =    nxny       -nx^2-nz^2    nynz
                 nxnz        nynz        -nx^2-n^y^2

		 
    --- */
    aux1=1-cos(theta);
    
    auxnx2=n[0]*n[0];
    auxny2=n[1]*n[1];
    auxnz2=n[2]*n[2];
    
    R[0][0] += -aux1*(auxny2+auxnz2);
    R[1][1] += -aux1*(auxnx2+auxnz2);
    R[2][2] += -aux1*(auxnx2+auxny2);
    
    aux2=aux1*n[0]*n[1];
    R[0][1] += aux2;
    R[1][0] += aux2;
    
    aux2=aux1*n[0]*n[2];
    R[0][2] += aux2;
    R[2][0] += aux2;
    
    aux2=aux1*n[1]*n[2];
    R[2][1] += aux2;
    R[1][2] += aux2;
  }
  
  else {
    /* --- Developpement limite autour de theta=0 en utilisant 
           la formule:
	   
	   R = I + sin(theta)/theta Sr + (1-cos(theta))/theta2 Sr^2

	   D'ou:
	   
	   R ~ I + (1-theta2/6)*Sr + (1/2-theta2/24)*Sr^2
    
       --- */

    /* --- Terme identite --- */
    for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      R[i][j]=0;

    R[0][0]=1;
    R[1][1]=1;
    R[2][2]=1;
    
    if (theta==0) return;
    
    /* --- Ajout du terme en 1-theta2/6  --- */
    
    /* theta2=theta*theta;*/
    
    aux1=1-theta2/6;
    
    aux2=aux1*r[2];
    R[0][1] += -aux2;
    R[1][0] +=  aux2;
    
    aux2=aux1*r[1];
    R[0][2] +=  aux2;
    R[2][0] += -aux2;
  
    aux2=aux1*r[0];
    R[1][2] += -aux2;
    R[2][1] += aux2;
  

    /* --- Ajout du terme en (1/2 - theta2/24) --- */
    aux1=0.5 - theta2/24;
    
    auxnx2=r[0]*r[0];
    auxny2=r[1]*r[1];
    auxnz2=r[2]*r[2];
    
    R[0][0] += -aux1*(auxny2+auxnz2);
    R[1][1] += -aux1*(auxnx2+auxnz2);
    R[2][2] += -aux1*(auxnx2+auxny2);
    
    aux2=aux1*r[0]*r[1];
    R[0][1] += aux2;
    R[1][0] += aux2;
    
    aux2=aux1*r[0]*r[2];
    R[0][2] += aux2;
    R[2][0] += aux2;
    
    aux2=aux1*r[1]*r[2];
    R[2][1] += aux2;
    R[1][2] += aux2;
    
  }
}





void Rmatrix_To_Rvector( double R[4][4], double * r )
{
  double xx, yy, zz, xy, yz, zx;
  double f, g;
  double cosinus, theta;
  
  cosinus = 0.5 * ( R[0][0] + R[1][1] +  R[2][2] -1.0 );
  
  if ( cosinus > 1.0 ) cosinus = 1.0;
  else if ( cosinus < -1.0 )  cosinus = -1.0;
  
  /* theta est dans [0, pi] 
   */
  theta = acos( (double)cosinus );
	
  if ( theta == 0.0 ) {

    r[0] = 0.0 ;
    r[1] = 0.0 ;	
    r[2] = 0.0 ;
    return;
  }


  if ( 1 + cosinus > sqrt(2.0)/2.0 ) {

    /* l'angle theta n'est pas voisin de pi
       on calcule f = sin( theta ) / theta
    */
    if ( theta > TH_TINY )
      f = sin( theta ) / theta;
    else 
      f = 1 - theta * theta * ( 1.0 / 6.0 - theta * theta / 120 );
    r[0] = ( R[2][1] - R[1][2] ) / ( 2.0 * f );
    r[1] = ( R[0][2] - R[2][0] ) / ( 2.0 * f );
    r[2] = ( R[1][0] - R[0][1] ) / ( 2.0 * f );

    return;
  }

  /* l'angle theta est voisin de pi,
     f tend vers 0, on ne peut plus l'inverser
  */
  g = (1.0 - cosinus ) / (theta * theta) ;
  
  xx = (R[0][0] - cosinus) / g;
  yy = (R[1][1] - cosinus) / g;
  zz = (R[2][2] - cosinus) / g;
  
  r[0] = sqrt( xx );
  r[1] = sqrt( yy );
  r[2] = sqrt( zz );
  
  if ( cosinus == -1.0 ) {

    /* l'angle est pi, 
       on a le choix sur sens de r
       on choisit ri > 0 avec ri = max (rx ry rz)
     */

    yz = (R[2][1] + R[1][2]) / (2.0 * g);
    zx = (R[0][2] + R[2][0]) / (2.0 * g);
    xy = (R[1][0] + R[0][1]) / (2.0 * g);
    
    if ( xx >= yy && xx >= zz ) {
      if ( xy < 0.0 ) r[1] *= -1.0;
      if ( zx < 0.0 ) r[2] *= -1.0;
    }
    else if ( yy >= zz ) {
      if ( xy < 0.0 ) r[0] *= -1.0;
      if ( yz < 0.0 ) r[2] *= -1.0;
    }
    else {
      if ( zx < 0.0 ) r[0] *= -1.0;
      if ( yz < 0.0 ) r[1] *= -1.0;
    }

    return;
  }

  /* il y a ambiguite sur le module de r
     ( theta ou 2 pi - theta ) et sur le signe des composantes
     on choisit theta dans [0 pi]
  */
  if ( R[2][1] - R[1][2] < 0.0 ) r[0] *= -1.0;
  if ( R[0][2] - R[2][0] < 0.0 ) r[1] *= -1.0;
  if ( R[1][0] - R[0][1] < 0.0 ) r[2] *= -1.0;

}
