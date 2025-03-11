/*
			     Program pivot_ete.c

This program accepts a walk in iw[][].  It selects a vertex on iw, and
performs a pivot on the shortest segment.  The pivot move is chosen random
from the octahedral group.

IN

   ***iw      -   walk
   n       -   number of edges in walk
   ***htabfull  -    full hash table
   beta      -   the contact potential
   *c       -   the number of contacts
   bet[k]      -    value of the temperature for the k-th chain
   iconf       -    index of the configuration finita in the k-th chain

OUT

    0      -  unsuccessful
    1      -  identity
    2,3... -  other types of moves
     ***iw       -    new polygon

local variables:
     neigh[3][6] -    neighbours of vertex
     piv         -    pivot
     **iu        -    work space
     **ix        -    list of lables
     *iv         -    work-space
     length      -    length of shorter path

These routines use a full hash-table in hash.c and hash.h.  It also calls the
random number generator in marsaglia.c.  To set the hash-table up, we
must initialise it from main before the routines here are called.  The
number of contacts in the walk is initialised to 0 since the new configuration
is a straight line. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pivot_chain_ete.h"
#include "hash_chain.h"
#include "marsagliazo.h"
#include "allo.h"
static int **store;
static int *iv,*iy;
static int neigh[3][6];

void setpivot(int n)
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


int pivot_chain_ete(iw, n, htabfull, mm,c,cen,sigma, beta,iconf, K)
 int ***iw;
 int n;
 int ***htabfull;
 int mm;
 int *c;
 int *cen;
 int **sigma;
 double beta;
 int iconf;
 double K;
{
 int i,p,sp,inc,i1,j;
 int label;
 int cont,cold,cnew,ce;
 int rot[3],sym[3],iv[3],point[3];

 for(i=0;i<3;i++) rot[i] = i;

/* Get an element of the octahedral group */
 octahedral(rot,sym);

/* Find the pivot */
 p = ranint_(&n);

/* Starting near the pivot, rotate the vertices and check for
intersections */
 if( p > (int) ( (float) n / 2. ) ) inc = 1;
 else inc = -1;
 for(i=0;i<3;i++) point[i] = iw[p][i][iconf];
 sp = p;
 p += inc;
 while( (p < n)&&(p > -1) ){
 for(i=0;i<3;i++) iv[i] = iw[p][i][iconf];

/* Rotate iv now */
 rotate(iv,rot,sym,point);

/* Is iv already in the hash-table? */
 if(inhash_c(iv,htabfull,mm,&label,iconf)) return 0;
 for(i=0;i<3;i++) store[p][i] = iv[i];
 p += inc;
 }

/* Update the hash table to find the number of contacts */
 cont = 0;
 ce = 0;
 p = sp;
 p += inc;
 while( (p < n)&&(p > -1) ){
  for(i=0;i<3;i++) iv[i] = iw[p][i][iconf];
  remhash_c(iv,htabfull,mm,p,iconf);
    for(i1=0;i1<6;i1++){
     for(j=0;j<3;j++) iy[j] = iv[j] + neigh[j][i1];
     if(inhash_c(iy,htabfull,mm,&label,iconf))
      {
	cont--;
        ce = ce - sigma[label][iconf]*sigma[p][iconf] - K*sigma[label][iconf]*sigma[p][iconf]*sigma[label][iconf]*sigma[p][iconf];
      }
     }
  p += inc;
  }

 p = sp;
 p += inc;
 while( (p < n)&&(p > -1) ){
  for(i=0;i<3;i++) iv[i] = store[p][i];
  addhash_c(iv,htabfull,mm,p,iconf);
    for(i1=0;i1<6;i1++){
     for(j=0;j<3;j++) iy[j] = iv[j] + neigh[j][i1];
     if(inhash_c(iy,htabfull,mm,&label,iconf))
      {
       cont++;
       ce = ce + sigma[label][iconf]*sigma[p][iconf] + K*sigma[label][iconf]*sigma[p][iconf]*sigma[label][iconf]*sigma[p][iconf];
      }
   }
  p += inc;
  }

 if(ran1real_() >= pow( beta , (double) (ce)))
 {
 p = sp;
 p += inc;
 while( (p < n)&&(p > -1) ){
  for(i=0;i<3;i++) iv[i] = store[p][i];
  remhash_c(iv,htabfull,mm,p,iconf);
  for(i=0;i<3;i++) iv[i] = iw[p][i][iconf];
  addhash_c(iv,htabfull,mm,p,iconf);
  p += inc;
  }
 return 0;
 }

/* The new configuration is accepted. We have to update 
   the old walk */
 p = sp;
 p += inc;
 while( (p < n)&&(p > -1) ){
 for(i=0;i<3;i++) iw[p][i][iconf] = store[p][i];
 p += inc;
 }
 *c += cont;
 *cen += ce;
 return 1;
 }

void octahedral(int rot[], int sym[])
{
 int i,j,k,l;
 k = 2;
 for(i=0;i<3;i++)
 {
  j = ranint_(&k);
  if(j)sym[i] = 1;
  else sym[i] = -1;
 }
 k = 3;
 for(i=0;i<2;i++)
 {
  j = ranint_(&k) + i;
  k -= 1;
  l = rot[j];
  rot[j] = rot[i];
  rot[i] = l;
 }
 return;
}

void rotate(int iv[], int rot[], int sym[], int point[])
{
 int iz[3],i;
 for(i=0;i<3;i++) iz[i] = iv[i];
 for(i=0;i<3;i++) iv[i] = sym[i] * (iz[rot[i]] - point[rot[i]] )
							  + point[i];
 return;
}

