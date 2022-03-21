#ifndef COMPUTE_RESIDUAL_HPP
#define COMPUTE_RESIDUAL_HPP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#define IS_PARALLEL 1
#else
#define IS_PARALLEL 0
#endif

double residual_serial(int, double*, double*);
double residual_parallel(int, double*, double*);
void print_solutions(int, double*);
double compute_residual(int, double*, double*);

#endif