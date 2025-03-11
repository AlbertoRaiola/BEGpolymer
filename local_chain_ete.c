/*                          Routine Local_ete

These routines perform Local moves on a walk.
It is sligtly different from the local.c file because 
in this case we have to handle walks and not polygons, so
we have to pay attention to the extreme vertices of the walk.

IN:  ***iw       -   walk [n][3]
     n          -   current number of vertices in iw
     ***htabfull -   filled hash-table 
     mm         -   length of hash-tables
     *c         -   number of contacts in polygon or walk
     beta      -    contact temperature
     iconf      -   indice della configurazione finita nella k-th chain

OUT: 0  -  unsuccessful attempt
     1  -  successful attempt
     c  - new contact number

The walk and the filled hash-table are also updated.
The moves are of to types:

		 _
 I     __| <-> _|     Corner move
 
        _
 II    | | <-> |_|    Crankshaft move

 The vertices in the local segment are v0,v1,v2,v3, and the edges between
 these are ed01,ed12,ed23.

 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "local_chain_ete.h"
#include "marsagliazo.h"
#include "hash_chain.h"


static int neigh[3][6];

void setlocal(int n)
{
 int i,j;
  for(i=0;i<6;i++) for(j=0;j<3;j++) neigh[j][i] = 0;
  for(i=0;i<3;i++){neigh[i][i] = 1;
		   neigh[i][i+3] = -1;}
  return;
}



int local_chain_ete(iw, n, htabfull, mm,c,cen,sigma,beta,iconf, K) 
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
 int v0[3],v1[3],v2[3],v3[3];
 int ed01[3],ed12[3],ed23[3];
 int contact,label0,label1,label2,label3;
 int ce;
 int label;
 int i,j,k;
 int lg0112,lg0123,lg1223;
 int newv1[3],newv2[3];
 int vertex[3];
 int cold,cnew;
 
/* Find all the labels and their vertices randomly */
 do{label1 = ranint_(&n);}
 while(label1<3 || label1 >n-3);
       /* In the case of polygons we don't need the while */
       label0 = label1 -1;
       label2 = label1 +1;
       label3 = label2 +1;

 for(i=0;i<3;i++)
 {
    v0[i] = iw[label0][i][iconf];
    v1[i] = iw[label1][i][iconf];
    v2[i] = iw[label2][i][iconf];
    v3[i] = iw[label3][i][iconf];
  }

/* Compute the edges (this gives the direction of the edges) */
 for(i=0;i<3;i++)
 {
    ed01[i] = v1[i] - v0[i];
    ed12[i] = v2[i] - v1[i];
    ed23[i] = v3[i] - v2[i];
  }

/* Determine the type of situation.
The logical variables assume 0 if edges are perpendicular, 1 if parallel and
-1 if anti-parallel */
 lg0112 = ed01[0]*ed12[0] + ed01[1]*ed12[1] + ed01[2]*ed12[2];
 lg0123 = ed01[0]*ed23[0] + ed01[1]*ed23[1] + ed01[2]*ed23[2];
 lg1223 = ed12[0]*ed23[0] + ed12[1]*ed23[1] + ed12[2]*ed23[2];

/* Return if we have a straight segment --- */
if((lg0112 == 1) && (lg0123 == 1) && (lg1223 == 1))return(0);

/* Case __|, vertex v2, corner move */
if((lg0112 == 1) && (lg1223 == 0))
{
   for(i=0;i<3;i++) newv2[i] = v1[i] + v3[i] - v2[i];
   if(inhash_c(newv2,htabfull,mm,&label,iconf))return(0);

 /* Count contacts */
   contact = 0;
   ce = 0;
   for(i=0;i<6;i++)
   {
    for(j=0;j<3;j++) vertex[j] = v2[j] + neigh[j][i];
    if(inhash_c(vertex,htabfull,mm,&label,iconf)) 
     {
      contact--;
      ce = ce - sigma[label][iconf]*sigma[label2][iconf] - K*pow(sigma[label][iconf]*sigma[label2][iconf], 2);
     }
   }
   for(i=0;i<6;i++)
   {
    for(j=0;j<3;j++) vertex[j] = newv2[j] + neigh[j][i];
    if(inhash_c(vertex,htabfull,mm,&label,iconf)) 
    {
    contact++;
    ce = ce + sigma[label][iconf]*sigma[label2][iconf] + K*pow(sigma[label][iconf]*sigma[label2][iconf], 2);
    }
   }

    if(pow( beta, (double)(ce)) < ran1double_()) return(0);

 /* Accept */
   remhash_c(v2,htabfull,mm,label2,iconf); 
   addhash_c(newv2,htabfull,mm,label2,iconf);
   for(i=0;i<3;i++) iw[label2][i][iconf] = newv2[i];
   *c += contact;
   *cen += ce;
   return 1;
  }


/* Case |__, vertex v1, corner move */
if((lg0112 == 0) && (lg1223 == 1)){
 for(i=0;i<3;i++) newv1[i] = v0[i] + v2[i] - v1[i];
 if(inhash_c(newv1,htabfull,mm,&label,iconf))return(0);
 
 /* Count contacts */
 contact = 0;
 ce = 0;
 for(i=0;i<6;i++){
  for(j=0;j<3;j++) vertex[j] = v1[j] + neigh[j][i];
  if(inhash_c(vertex,htabfull,mm,&label,iconf))
   {
   contact--;
   ce = ce - sigma[label1][iconf]*sigma[label][iconf] - K*pow(sigma[label1][iconf]*sigma[label][iconf], 2);
   }
  }
 for(i=0;i<6;i++){
  for(j=0;j<3;j++) vertex[j] = newv1[j] + neigh[j][i];
  if(inhash_c(vertex,htabfull,mm,&label,iconf))
    {
    contact++;
    ce = ce + sigma[label1][iconf]*sigma[label][iconf] + K*pow(sigma[label1][iconf]*sigma[label][iconf], 2);
    } 
   }

    if(pow( beta, (double)(ce)) < ran1double_()) return(0);

 /* Accept */
 remhash_c(v1,htabfull,mm,label1,iconf); 
 addhash_c(newv1,htabfull,mm,label1,iconf);
 for(i=0;i<3;i++) iw[label1][i][iconf] = newv1[i];
 *c += contact;
 *cen += ce;
 return 2;}

/*      _       */
/* Case  |_, either vertex v1 or v2, corner moves */

if((lg0112 == 0) && (lg1223 == 0) && (lg0123 != -1))
{
  i=2; 
  if(ranint_(&i))
  {
 for(i=0;i<3;i++) newv1[i] = v0[i] + v2[i] - v1[i];
 if(inhash_c(newv1,htabfull,mm,&label,iconf))return(0);

 /* Count contacts */
 contact = 0;
 ce = 0;
 for(i=0;i<6;i++){
  for(j=0;j<3;j++) vertex[j] = v1[j] + neigh[j][i];
  if(inhash_c(vertex,htabfull,mm,&label,iconf))
    {
    contact--;
    ce = ce - sigma[label][iconf]*sigma[label1][iconf] - K*pow(sigma[label][iconf]*sigma[label1][iconf], 2);
    }
    }
 for(i=0;i<6;i++){
  for(j=0;j<3;j++) vertex[j] = newv1[j] + neigh[j][i];
  if(inhash_c(vertex,htabfull,mm,&label,iconf))
   {
     contact++;
    ce = ce + sigma[label][iconf]*sigma[label1][iconf] + K*pow(sigma[label][iconf]*sigma[label1][iconf], 2);
   }
   }

    if(pow( beta, (double)(ce)) < ran1double_()) return(0);

 /* Accept */
 remhash_c(v1,htabfull,mm,label1,iconf); 
 addhash_c(newv1,htabfull,mm,label1,iconf);
 for(i=0;i<3;i++) iw[label1][i][iconf] = newv1[i];
 *c += contact;
 *cen += ce;
 return 3;
   }
   else
   {
   for(i=0;i<3;i++) newv2[i] = v1[i] + v3[i] - v2[i];
   if(inhash_c(newv2,htabfull,mm,&label,iconf))return(0);

 /* Count contacts */
   contact = 0;
   ce =0;
   for(i=0;i<6;i++){
   for(j=0;j<3;j++) vertex[j] = v2[j] + neigh[j][i];
   if(inhash_c(vertex,htabfull,mm,&label,iconf))
     {  
       contact--;
       ce = ce - sigma[label2][iconf]*sigma[label][iconf] - K*pow(sigma[label2][iconf]*sigma[label][iconf], 2);
     }
   }
   for(i=0;i<6;i++){
   for(j=0;j<3;j++) vertex[j] = newv2[j] + neigh[j][i];
   if(inhash_c(vertex,htabfull,mm,&label,iconf))
     {
       contact++;
       ce = ce + sigma[label2][iconf]*sigma[label][iconf] + K*pow(sigma[label2][iconf]*sigma[label][iconf], 2);
     }
   }

    if(pow( beta, (double)(ce)) < ran1double_()) return(0);

 /* Accept */
 remhash_c(v2,htabfull,mm,label2,iconf); 
 addhash_c(newv2,htabfull,mm,label2,iconf);
 for(i=0;i<3;i++) iw[label2][i][iconf] = newv2[i];
 *c += contact;
 *cen += ce;
 return 4;

   }
}

/* Case |_|, crankshaft */
if((lg0112 == 0) && (lg1223 == 0) && (lg0123 == -1)){

/* Pick one of four directions perpendicular to ed12 */
do{ for(i=0;i<3;i++) vertex[i] = 0; i=3; i=ranint_(&i); j=2;
		     vertex[i] = 2*ranint_(&j)-1;}while(ed12[i] != 0);

for(i=0;i<3;i++) newv1[i] = v0[i] + vertex[i];
for(i=0;i<3;i++) newv2[i] = v3[i] + vertex[i];
 if(inhash_c(newv1,htabfull,mm,&label,iconf)) return 0;
 if(inhash_c(newv2,htabfull,mm,&label,iconf)) return 0;

/*  Count contacts */ 
 contact = 0;
 ce = 0;
 for(i=0;i<6;i++){
  for(j=0;j<3;j++) vertex[j] = v1[j] + neigh[j][i];
  if(inhash_c(vertex,htabfull,mm,&label,iconf))
    {
      contact--;
      ce = ce - sigma[label1][iconf]*sigma[label][iconf] - K*pow(sigma[label1][iconf]*sigma[label][iconf], 2);
    }
  }
 for(i=0;i<6;i++){
  for(j=0;j<3;j++) vertex[j] = newv1[j] + neigh[j][i];
  if(inhash_c(vertex,htabfull,mm,&label,iconf))
    {
      contact++;
      ce = ce + sigma[label1][iconf]*sigma[label][iconf] + K*pow(sigma[label1][iconf]*sigma[label][iconf], 2);
    }
  }
 for(i=0;i<6;i++){
  for(j=0;j<3;j++) vertex[j] = v2[j] + neigh[j][i];
  if(inhash_c(vertex,htabfull,mm,&label,iconf))
    { 
      contact--;
      ce = ce - sigma[label2][iconf]*sigma[label][iconf] - K*pow(sigma[label2][iconf]*sigma[label][iconf], 2);
    }
  }
 for(i=0;i<6;i++){
  for(j=0;j<3;j++) vertex[j] = newv2[j] + neigh[j][i];
  if(inhash_c(vertex,htabfull,mm,&label,iconf))
    {
       contact++;
      ce = ce + sigma[label2][iconf]*sigma[label][iconf] + K*pow(sigma[label2][iconf]*sigma[label][iconf], 2);
    }
  }

  contact = contact + 2;
  ce = ce + 2*sigma[label1][iconf]*sigma[label2][iconf] + 2*K*pow(sigma[label1][iconf]*sigma[label2][iconf], 2);

    if(pow( beta, (double)(ce)) < ran1double_()) return(0);

  /* Accept */  
 remhash_c(v1,htabfull,mm,label1,iconf); 
 remhash_c(v2,htabfull,mm,label2,iconf);
 addhash_c(newv1,htabfull,mm,label1,iconf); 
 addhash_c(newv2,htabfull,mm,label2,iconf);
 for(i=0;i<3;i++) iw[label1][i][iconf] = newv1[i];
 for(i=0;i<3;i++) iw[label2][i][iconf] = newv2[i];
 *c += contact;
 *cen += ce;
 return 5;} 

return 10;

}
