/******************************************************************************
* FILE: omp_bug5.c
* DESCRIPTION:
*   Using SECTIONS, two threads initialize their own array and then add
*   it to the other's array, however a deadlock occurs.
* AUTHOR: Blaise Barney  01/29/04
* LAST REVISED: 04/06/05
******************************************************************************/
/******************************************************************************
 * Deadlock occurs because the thread executing section 1 (thread 1) acquires
 * locka while thread 2 acquires lockb which leads to both threads waiting for
 * a resource the other thread holds.
 * An initial attempt at solving this might be to release locka in thread 1 before
 * holding locka again. The same thing is done in thread 2 with lockb. This way,
 * both threads acquire both locks in the same order (please see the section
 * marked sol1). But this leads to inconsistent output in the case that some thread
 * acquires both locks before the other can initialize its data. We want to make
 * sure both threads initialize their arrays before adding the two arrays. To do
 * this, we'll have to enforce some serialization. We can do this by seperating the
 * initialization and addition into two "sections" constructs as done in the final
 * solution.
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define N 1000000
#define PI 3.1415926535
#define DELTA .01415926535

int main (int argc, char *argv[]) 
{
int nthreads, tid, i;
float a[N], b[N];
omp_lock_t locka, lockb;

/* Initialize the locks */
omp_init_lock(&locka);
omp_init_lock(&lockb);

/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel shared(a, b, nthreads, locka, lockb) private(tid)
  {

  /* Obtain thread number and number of threads */
  tid = omp_get_thread_num();
  #pragma omp master
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d starting...\n", tid);
  #pragma omp barrier

  #pragma omp sections
  {
    // here, we don't need to acquire locks since only one thread executes each section
    #pragma omp section
    {
      for (i=0; i<N; i++)
        a[i] = i * DELTA;
    }

    #pragma omp section
    {
      for (i=0; i<N; i++)
        b[i] = i * PI;
    }

  }

  // both threads wait before moving forward
  #pragma omp sections
  {
    #pragma omp section
    {
      omp_set_lock(&locka);
      omp_set_lock(&lockb);

      printf("Thread %d adding a[] to b[]\n",tid);
      for (i=0; i<N; i++)
        b[i] += a[i];

      omp_unset_lock(&locka);
      omp_unset_lock(&lockb);
    }

    #pragma omp section
    {
      omp_set_lock(&locka);
      omp_set_lock(&lockb);

      printf("Thread %d adding b[] to a[]\n",tid);
      for (i=0; i<N; i++)
        a[i] += b[i];

      
      omp_unset_lock(&locka);
      omp_unset_lock(&lockb);
    }
  }
  
  /* sol1
  #pragma omp sections nowait
    {
    #pragma omp section
      {
      printf("Thread %d initializing a[]\n",tid);
      omp_set_lock(&locka);
      for (i=0; i<N; i++)
        a[i] = i * DELTA;
      omp_unset_lock(&locka);
      omp_set_lock(&locka);
      omp_set_lock(&lockb);
      printf("Thread %d adding a[] to b[]\n",tid);
      for (i=0; i<N; i++)
        b[i] += a[i];
      omp_unset_lock(&lockb);
      omp_unset_lock(&locka);
      }

    #pragma omp section
      {
      printf("Thread %d initializing b[]\n",tid);
      omp_set_lock(&lockb); // bug
      for (i=0; i<N; i++)
        b[i] = i * PI;
      omp_unset_lock(&lockb);
      omp_set_lock(&locka); // bug
      omp_set_lock(&lockb);
      printf("Thread %d adding b[] to a[]\n",tid);
      for (i=0; i<N; i++)
        a[i] += b[i];
      omp_unset_lock(&lockb);
      omp_unset_lock(&locka);
      // omp_unset_lock(&locka);
      // omp_unset_lock(&lockb);
      }
    }*/ /* end of sections */
  }  /* end of parallel region */

  printf("Printing final value of a and b:\n");
  printf("a[N-1]: %f\tb[N-1]: %f\n", a[N-1], b[N-1]);

}

