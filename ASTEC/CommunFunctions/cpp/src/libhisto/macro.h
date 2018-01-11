/***************************************************************************
 * macro.h
 *
 * $Author: greg $
 * $Revision: 1.1 $
 * $Log: macro.h,v $
 * Revision 1.1  2002/10/18 18:02:07  greg
 * *** empty log message ***
 *
 * Revision 1.2  2001/12/07 15:48:11  greg
 * En-tete
 *
 *
 *
 * $Id: macro.h,v 1.1 2002/10/18 18:02:07 greg Exp $
 ***************************************************************************/

#ifndef MACRO_H
#define MACRO_H
 
#ifdef __cplusplus
extern "C" {
#endif


/* --- Constantes à moi que j'ai --- */
#define INFINITECUTOFF 1e10  



/* --- Quelques constantes universelles --- */

#define Pi             3.14159265358979      /* pi */
#define InvSqrt2Pi     0.39894228040143      /* 1/sqrt(2*pi) */ 
#define Exp1	       2.71828182845905      /* e */ 
#define Log2           0.69314718055995      /* log(2) */ 


/* --- Quelques macros utiles dans la vie --- */

#define Max(A,B) ( (A) > (B) ? (A) : (B) ) 
#define Min(A,B) ( (A) < (B) ? (A) : (B) ) 

/* Warning = le résultat de ces macros dépend du compilo ou de l'architecture */
#define Floor(a) ((a)>0.0 ? (int)(a):( ((int)(a)-a)!= 0.0 ? (int)(a)-1 : (int)(a)))
#define Round(a) (Floor(a+0.5))
#define Ceil(a)  (-(Floor(-(a))))

#define uFloor(a) ( (int)(a) )
#define uRound(a) ( (int)(a+0.5) )
#define uCeil(a) ( ( (int)(a)-a )!=0.0 ? (int)(a+1) : (int)(a) )


#define IsPositive(a) (( (a) < 0 ) ? -1 : 1)
#define Abs(a) ( (a) > 0.0 ? (a) : (-(a)) )
/*#define Sign(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))*/
#define Sign(a,b) ((b) > 0.0 ? Abs(a) : -Abs(a))
#define Shft(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)


#define Sqr(a) ( (a)*(a) )
#define Cube(a) ( (a)*(a)*(a) )

#define Sinc(x)    \
 ( (x)==0.0 ? 1 : sin((x))/(x) )

#define Gaussian(x,scale)     \
 exp ( - (x)*(x) / (scale) )


  /* Fonctions de cout et pondérations associées à 
     différents M-estimateurs. 
     
     Remarque: par convention, toutes les fonctions de cout sont 
     de la forme:
     
     Distance(x,c) = c^2 f(x/c).
  */
 
#define HuberDistance(x,c) \
  ( ( (x) < (c) ) && ( (x) > (-(c)) ) ? (0.5*Sqr(x)) : ( (c)*(Abs(x)-0.5*(c)) ) )
#define HuberWeight(x,c) \
  ( ( (x) < (c) ) && ( (x) > (-(c)) ) ? 1.0 : ( (c)/(Abs(x)) ) )
#define HuberStdCutOff 1.345
#define HuberKvalue 0.467


#define FairDistance(x,c) \
  ( Sqr((c)) * ( Abs((x)/(c)) - log(1.0 + Abs((x)/(c))) ) )
#define FairWeight(x,c) \
  ( 1.0 / (1.0 + Abs((x)/(c)) ) )
#define FairStdCutOff 1.3998
#define FairKvalue 0.3007

    
#define CauchyDistance(x,c) \
  ( 0.5 * Sqr((c)) * log( 1.0 + Sqr((x)/(c)) ) )
#define CauchyWeight(x,c) \
  ( 1.0 / ( 1.0 + Sqr((x)/(c)) ) )
#define CauchyStdCutOff 2.385
#define CauchyKvalue 0.349


#define GemanDistance(x,c) \
  ( 0.5 * Sqr((x)) / ( 1.0 + Sqr((x)/(c)) ) )
#define GemanWeight(x,c) \
  ( 1.0 / ( Sqr( 1.0 + Sqr((x)/(c)) ) ) )
#define GemanStdCutOff 3.648
#define GemanKvalue 0.416

#define WelschDistance(x,c) \
  ( 0.5 * Sqr((c)) * ( 1.0 - exp( - Sqr((x)/(c)) ) ) )
#define WelshWeight(x,c) \
  ( exp( - Sqr((x)/(c)) ) ) 
#define WelshStdCutOff 4.685
#define WelshKvalue 0.468

#define TukeyDistance(x,c) \
  ( ( Abs(x) < (c) ) ? ( (Sqr((c))/6.0) * ( 1.0 - Cube( 1.0 - Sqr((x)/(c)))) ) : (Sqr((c))/6.0) )
#define TukeyWeight(x,c) \
  ( ( Abs(x) < (c) ) ? ( Sqr( 1.0 - Sqr( (x)/(c))) ) : 0.0 )  
#define TukeyStdCutOff 4.685
#define TukeyKvalue 0.437

#define NaiveDistance(x,c) \
  ( ( (x) < (c) ) && ( (x) > (-(c)) ) ? (0.5*Sqr(x)) : (0.5*Sqr(c)) ) 
#define NaiveWeight(x,c) \
  ( ( (x) < (c) ) && ( (x) > (-(c)) ) ? 1.0 : 0.0 )
#define NaiveStdCutOff 2.795
#define NaiveKvalue 0.4952

  /* Efficacité de ces estimateurs.
     Valeur de cut-off donnant 95% d'efficacité et si nécessaire, 
     entre parenthèses, celle qui correspond au point de rupture à 50%
     du S-esimateur associé.
 
       Huber => c = 1.345
       Fair => c = 1.3998
       Cauchy => c = 2.3849
       Geman => c = 3.648 (valeur approchée par Alexis)
       Welsh => c = 4.6851 
       Naive => c = 2.7955 (1.041)
       Tukey => c = 4.685 (1.547)

  */


  /*
    #define MDistance(x,c) \
    ( GemanDistance(x,c) )
    
    #define MWeight(x,c) \
    ( GemanWeight(x,c) )
    
    #define MStdCutOff GemanStdCutOff

    #define MName "Geman-McClure"

  */

 #define MDistance(x,c) \
     ( GemanDistance(x,c) ) 
  
 #define MWeight(x,c) \
     ( GemanWeight(x,c) ) 
  
 #define MStdCutOff GemanStdCutOff 

 #define MKvalue GemanKvalue 
  
 #define MName "Geman-McClure"


  /*#define MDistance(x,c) \
    ( HuberDistance(x,c) )
    
    #define MWeight(x,c) \
    ( HuberWeight(x,c) )
    
    #define MStdCutOff HuberStdCutOff
    
    #define MKvalue HuberKvalue
    
    #define MName "Huber"*/


#ifdef __cplusplus
}
#endif
 
#endif  
