/******************************************************************************
* FILE: omp_bug6.c
* DESCRIPTION:
*   This program compiles and runs fine, but produces the wrong result.
*   Compare to omp_orphan.c.
* AUTHOR: Blaise Barney  6/05
* LAST REVISED: 06/30/05
******************************************************************************/
/******************************************************************************
 * Bugs:
 * 1.Calling dotprod() in parallel construct in main creates multiple stack frames
 *  for dotprod() all of which have a local sum variable. Hence, the compilation
 *  error since sum is now private for each thread. 
 *  We actually want to parallelize the for loop inside dotprod(), so we can either
 *  move the loop into main or move the parallel clause to dotprod() which is what
 *  I have done.
 * 2. Another issue was that dotprod didn't return sum.
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define VECLEN 100

float a[VECLEN], b[VECLEN];

float dotprod ()
{
int i,tid;
float sum;

#pragma omp parallel private(tid)
{
  tid = omp_get_thread_num();
  #pragma omp for
  for (i=0; i < VECLEN; i++)
    {
    sum = sum + (a[i]*b[i]);
    printf("  tid= %d i=%d\n",tid,i);
    }
}
return sum; // bug, original program did not return
}


int main (int argc, char *argv[]) {
int i;
float sum;

for (i=0; i < VECLEN; i++)
  a[i] = b[i] = 1.0 * i;
sum = 0.0;

// #pragma omp parallel shared(sum) // bug
  sum = dotprod();

printf("Sum = %f\n",sum);

}

