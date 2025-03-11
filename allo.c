/*                          Program Allofloat

We use these routines to allocate memory for floating point arrays
with two indices, of the kind float [nx][ny].  The sintax of the
call to allofloat is of the form:

float **imatrix;

imatrix = (float **)allofloat(nx,ny);

  There are three routines:

    allofloat   -   allocates an floating point matrix
    freefloat   -   free memory assigned to a matrix
    frerror   -   error shooting in memory allocation
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "allo.h"

float **allofloat( int nx , int ny )
{
 int i;
 float **m;

/* Allocate pointers to the rows */

 m = (float **)malloc(nx * sizeof(float *));
 if(!m) frerror("allocation failure 1 in allofloat");

/* Allocate memory to the matrix */

 m[0] = (float *)malloc(nx * ny * sizeof(float));
 if(!m[0]) frerror("allocation failure 2 in allofloat");

/* Assign row pointers */

 for(i=0;i<nx;i++) m[i] = m[0] + i * ny;

/* Return m to the main program */

 return m;
}

void freefloat(float **m)
{
 free((void *)m[0]);
 free((void *)m);
}

void frerror(char *error_text)
{
 fprintf(stderr,"Run time error ...\n");
 fprintf(stderr,"%s\n",error_text);
 fprintf(stderr,"...now exiting to system...\n");
 exit(1);
}

/*                          Program Allodouble

We use these routines to allocate memory for floating point arrays
with two indices, of the kind float [nx][ny].  The sintax of the
call to allofloat is of the form:

double **imatrix;

imatrix = (double **)allodouble(nx,ny);

  There are three routines:

    allodouble   -   allocates an floating point matrix
    freedouble   -   free memory assigned to a matrix
    frerror   -   error shooting in memory allocation
*/


double **allodouble( int nx , int ny )
{
 int i;
 double **m;

/* Allocate pointers to the rows */

 m = (double **)malloc(nx * sizeof(double *));
 if(!m) frerror("allocation failure 1 in allodouble");

/* Allocate memory to the matrix */

 m[0] = (double *)malloc(nx * ny * sizeof(double));
 if(!m[0]) frerror("allocation failure 2 in allodouble");

/* Assign row pointers */

 for(i=0;i<nx;i++) m[i] = m[0] + i * ny;

/* Return m to the main program */

 return m;
}

void freedouble(double **m)
{
 free((void *)m[0]);
 free((void *)m);
}

/*                          Program Alloint

We use these routines to allocate memory for integer arrays
with two indices, of the kind int [nx][ny].  The sintax of the
call to alloint is of the form:

int **imatrix;

imatrix = (int **)alloint(nx,ny);

  There are three routines:

    alloint   -   allocates an integer matrix
    freeint   -   free memory assigned to a matrix
    frerror   -   error shooting in memory allocation
*/

int **alloint( int nx , int ny )
{
 int i;
 int **m;

/* Allocate pointers to the rows */

 m = (int **)malloc(nx * sizeof(int *));
 if(!m) frerror("allocation failure 1 in alloint");

/* Allocate memory to the matrix */

 m[0] = (int *)malloc(nx * ny * sizeof(int));
 if(!m[0]) frerror("allocation failure 2 in alloint");

/* Assign row pointers */

 for(i=0;i<nx;i++) m[i] = m[0] + i * ny;

/* Return m to the main program */

 return m;
}

void freeint(int **m)
{
 free((void *)m[0]);
 free((void *)m);
}

