/*
These routines use a full hash-table in hash.c and hash.h.  It also calls the
random number generator in marsaglia.c.  To set the hash-table up, we
must initialise it from main before the routines here are called.  The
number of contacts in the walk is initialised to 0 since the new configuration
is a straight line. 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ising_saw3d_h.h"
#include "hash_chain.h"
#include "marsagliazo.h"
#include "allo.h"
static int **store;
static int *iv,*iy;
static int neigh[3][6];

void setising(int n)
{
 int i,j;
  
  iv = (int *)malloc(3 * sizeof(int));
  iy = (int *)malloc(3 * sizeof(int));
  store   = (int **)alloint(n,3);
  for(i=0;i<6;i++) for(j=0;j<3;j++) neigh[j][i] = 0;
  for(i=0;i<3;i++){neigh[i][i] = 1;
		   neigh[i][i+3] = -1;}

return;
}


int ising_saw3d_h(iw, n, htabfull, mm,cen,sigma, beta,h,iconf, Delta, K)
 int ***iw;
 int n;
 int ***htabfull;
 int mm;
 int *cen;
 int **sigma;
 double beta;
 double h;
 int iconf;
 double Delta;
 double K;
{
 int i,p,i1,j;
 int label;
 int iv[3];
 int eno,enn,dener,dsener;
 double randomNumber;
 int newState = 0;
 dsener = *cen;
 for (p=0;p<n;p++)
{
  eno = enn = 0;
  for(i=0;i<3;i++) iv[i] = iw[p][i][iconf];
  
  do
  {
  		randomNumber = ran1real_();
   		if(randomNumber < 0.333){newState = -1;}
		else if (randomNumber < 0.666 && randomNumber >= 0.333 ) {newState = 0;}
		else {newState = 1;}
	} while (newState == sigma[p][iconf]);
	

    for(i1=0;i1<6;i1++){
     for(j=0;j<3;j++) iy[j] = iv[j] + neigh[j][i1];
     if(inhash_c(iy,htabfull,mm,&label,iconf))
     {
      eno = eno + sigma[p][iconf]*sigma[label][iconf] + K*sigma[p][iconf]*sigma[label][iconf]*sigma[p][iconf]*sigma[label][iconf];
      enn = enn + newState*sigma[label][iconf] + K*newState*sigma[label][iconf]*newState*sigma[label][iconf];
     }
    }
    eno = eno - Delta*sigma[p][iconf]*sigma[p][iconf];
    enn = enn - Delta*newState*newState;
   dener = enn-eno;
 if(ran1real_() <= pow( beta , (double) (dener)))
  {sigma[p][iconf] = newState; dsener = dsener + dener;}
  else {dsener = dsener;}
 
 }

 *cen = dsener; 
 return 1;
 }


