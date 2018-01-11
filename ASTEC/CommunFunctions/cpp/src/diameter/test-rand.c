
#include <stdio.h>
#include <time.h>
#include <rand.h>

#define MAX 10

int main()
{
  int i, j;
  int t[MAX];
  int n = 100000;
  int init = time(0);
  
  for (i=0; i<MAX; i++ ) t[i] =0;

  _SetRandomSeed(init);

  for (i=0; i<n; i++ ) {
    j = _GetRandomIntNb( 0, MAX-1 );
    t[j] ++;
  }

  for (i=0; i<MAX; i++ ) printf( " %6d", i );
  printf( "\n" );
  for (i=0; i<MAX; i++ ) printf( " %6d", t[i] );
  printf( "\n" );


}
