/******************************************************************************
* FILE: omp_bug2.c
* DESCRIPTION:
*   Another OpenMP program with a bug. 
* AUTHOR: Blaise Barney 
* LAST REVISED: 04/06/05 
******************************************************************************/
/*****************************************************************************
 * Bugs:
 * 1. tid should be thread-private, so should i and nthreads, since they are 
 *    overwritten (nthreads is written to by only one thread)
 * 2. since total is shared, it should be initialized before entering
 *    parallel construct
 * 3. total should be reduced using +, otherwise threads will overwrite the shared
 *    variable total.
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) 
{
int nthreads, i, tid;
float total = 0.0; // correct initialization

/*** Spawn parallel region ***/

// #pragma omp parallel // bug here

#pragma omp parallel private(tid, i, nthreads) // total is shared by default
  {
  /* Obtain thread number */
  tid = omp_get_thread_num();
  /* Only master thread does this */
  if (tid == 0) {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
  }
  printf("Thread %d is starting...\n",tid);

  #pragma omp barrier

  /* do some work */
  // total = 0.0; // bug here
  // #pragma omp for schedule(dynamic,10) // bug here
  #pragma omp for schedule(dynamic,10) reduction(+:total)
  for (i=0; i<1000000; i++) 
     total = total + i*1.0;

  printf ("Thread %d is done! Total= %e\n",tid,total);

  } /*** End of parallel region ***/
}
