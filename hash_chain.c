/*                   Routines hash_chain.c

The routines in this file operate on a global hash-table htab in
3 dimensions.  The table is externally defined as htab[n][3][nchain], and
contains the vertices of the i-th vertex in htab[m][0][nchain] and
htab[m][1][nchain], where we find m via the hash-function islot.
nchain is the number of chains considered in the catene program.


The following routines are present:

void sethash_c:  initialises htab and the hash parameters in kt.
                 Call this routine once at the beginning of run.

int islot:  returns the slot of the vertex.

int inhash_c:  queries iv in htab[iconf], returns a 0 if false.

void addhash_c:  adds iv to htab[iconf].

void remhash_c:  removes iv from htab[iconf].

int inhash:  queries iv in htab, returns a 0 if false.

void addhash:  adds iv to htab.

void remhash:  removes iv from htab.
*/

#include <math.h>
#include "hash_chain.h"
#include "string.h"
#include "stdlib.h"

#define MAXINT -2147483648

      static int *kt;

/*
    Set the htable and the hash-parameters up.
*/

      void sethash_c( int ***htab, int mm, int nchain )
      {
        int i,j,k;
        float x;
        kt = (int *)malloc(3 * sizeof(int));
        for(i=0;i<3;i++)
        {
          x = (float) (i+1) / 4.;
          kt[i] = (int) (pow( (float) mm , x ));
        }
	for(k=0;k<nchain;k++)
	{
        for(i=0;i<mm;i++) for(j=0;j<4;j++) htab[i][j][k] = MAXINT;
	}
        return ;
      }

/*
    Calculate the slot of vertex iv in a table of length mm.
*/

      int islot( int *iv, int mm )
      {
        int i,j;
        i = 0;
        for(j=0;j<3;j++) i += iv[j] * kt[j];
        return (i + 10000000)%mm;
      }

/*
    Queries iv in htab[iconf]?
*/

      int inhash_c( int *iv, int ***htab , int mm, int *label,int iconf )
      {
        int i,j;
        i = islot(iv,mm);
        while(htab[i][0][iconf] != MAXINT)
        {
          j = 0;
          while(j < 3)
          {
            if(htab[i][j][iconf] == iv[j])
            {
              if(j == 2)
              {
                *label = htab[i][3][iconf];
                return 1;
              }
              else j = j + 1;
            }
            else break;
          }
          i = (i+1)%mm;
        }
        return 0;
      }

/*
    Add iv to the table htab.
*/

      int addhash_c( int *iv , int ***htab , int mm , int label, int iconf)
      {
        int i,j;
        i = islot(iv,mm);
        while(htab[i][0][iconf] != MAXINT) i = (i+1)%mm;
        for(j=0;j<3;j++) htab[i][j][iconf] = iv[j];
        htab[i][3][iconf] = label;
        return i;
      }

/*
    Remove the vertex iv from htab.
*/

      void remhash_c( int *iv , int ***htab , int mm , int label, int iconf)
      {
      int i,j,k,l;
      int im[10];
        i = islot(iv,mm);
        if(htab[i][0][iconf] == MAXINT ) return;
        j = 0;
        while(j < 3)
        {
          if((htab[i][j][iconf] != iv[j]) || (htab[i][3][iconf] != label))
          {
            i = (i+1)%mm;
            j = 0;
          }
          else j = j + 1;
        }
        htab[i][0][iconf] = MAXINT;
        j = (i+1)%mm;
        while(htab[j][0][iconf] != MAXINT)
        {
          for(k=0;k<3;k++) im[k] = htab[j][k][iconf];
          l = islot(im,mm);
          if(((i < l)&&(l <= j))||((i < l)&&(j < i))||
                 ((l <= j)&&(j < i))) j = (j+1)%mm;
          else
          {
            for(k=0;k<4;k++) htab[i][k][iconf] = htab[j][k][iconf];
            htab[j][0][iconf] = MAXINT;
            i = j;
            j = (j+1)%mm;
          }
        }
        return;
      }




/*
    Queries iv in htab?
*/

      int inhash( int *iv, int **htab , int mm, int *label)
      {
        int i,j;
        i = islot(iv,mm);
        while(htab[i][0] != MAXINT)
        {
          j = 0;
          while(j < 3)
          {
            if(htab[i][j] == iv[j])
            {
              if(j == 2)
              {
                *label = htab[i][3];
                return 1;
              }
              else j = j + 1;
            }
            else break;
          }
          i = (i+1)%mm;
        }
        return 0;
      }

/*
    Add iv to the table htab.
*/

      int addhash( int *iv , int **htab , int mm , int label)
      {
        int i,j;
        i = islot(iv,mm);
        while(htab[i][0] != MAXINT) i = (i+1)%mm;
        for(j=0;j<3;j++) htab[i][j] = iv[j];
        htab[i][3] = label;
        return i;
      }

/*
    Remove the vertex iv from htab.
*/

      void remhash( int *iv , int **htab , int mm , int label)
      {
      int i,j,k,l;
      int im[10];
        i = islot(iv,mm);
        if(htab[i][0] == MAXINT ) return;
        j = 0;
        while(j < 3)
        {
          if((htab[i][j] != iv[j]) || (htab[i][3] != label))
          {
            i = (i+1)%mm;
            j = 0;
          }
          else j = j + 1;
        }
        htab[i][0] = MAXINT;
        j = (i+1)%mm;
        while(htab[j][0] != MAXINT)
        {
          for(k=0;k<3;k++) im[k] = htab[j][k];
          l = islot(im,mm);
          if(((i < l)&&(l <= j))||((i < l)&&(j < i))||
                 ((l <= j)&&(j < i))) j = (j+1)%mm;
          else
          {
            for(k=0;k<4;k++) htab[i][k] = htab[j][k];
            htab[j][0] = MAXINT;
            i = j;
            j = (j+1)%mm;
          }
        }
        return;
      }
