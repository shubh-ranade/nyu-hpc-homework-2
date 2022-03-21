#include "compute-residual.hpp"

// void print_solution(int N, double *u) {
//     int i,j,index;
//     for(i = 0; i < N; i++) {
//         for(j = 0; j < N; j++) {
//             index = N * i + j;
//             printf("%3f  ", u[index]);
//         }
//         printf("\n");
//     }
// }

double compute_residual(int N, double* u, double* f) {
    if(IS_PARALLEL)
        return residual_parallel(N, u, f);
    else
        return residual_serial(N, u, f); 
}

double residual_parallel(int N, double *u, double *f) {
    double res = 0.0;
    double diff;
    double h = 1.0 / (N + 1);
    double h2 = h * h;
    int i,j, index;
    #pragma omp parallel shared(h2, h) private(i,j,index,diff)\
        reduction(+:res)
    {
        #pragma omp sections nowait
        {
            #pragma omp section
            {
                // first row:
                diff = f[0] + (-4 * u[0] + u[1] + u[N]) / h2;
                res += diff * diff;
                for(j = 1; j < N-1; j++) {
                    diff = f[j] + (-4 * u[j] + u[j-1] + u[j+1] + u[j+N])\
                                    / h2;
                    res += diff * diff;
                }
                diff = f[N-1] + (-4 * u[N-1] + u[N-2] + u[2*N-1]) / h2;
                res += diff * diff;
            }

            #pragma omp section
            {
                // last row:
                index = N*(N-1);
                diff = f[index] + (u[index + 1] + u[index - N] - 4 * u[index])\
                                    / h2;
                res += diff * diff;
                for(j = 1; j < N-1; j++) {
                    index = N*(N-1) + j;
                    diff = f[index] + (u[index - 1] + u[index + 1]\
                                + u[index - N] - 4 * u[index]) / h2;
                    res += diff * diff;
                }
                index = N * N - 1;
                diff = f[index] + (u[index - 1] + u[index - N] - 4 * u[index])\
                                / h2;
                res += diff * diff;
            }

        }
        // inner rows:
        #pragma omp for
        for(i = 1; i < N - 1; i++){
            index = N * i;
            diff = f[index] + (u[index + 1] + u[index + N] + u[index - N] \
                        - 4 * u[index]) / h2;
            res += diff * diff;
            for(j = 1; j < N - 1; j++) {
                index = N * i + j;
                diff = f[index] + (u[index + 1] + u[index - 1] + u[index + N] \
                            + u[index - N] - 4 * u[index]) / h2;
                res += diff * diff;
            }
            index = N * i + (N - 1);
            diff = f[index] + (u[index - 1] + u[index - N] + u[index + N] \
                        - 4 * u[index]) / h2;
            res += diff * diff;
        }
    }
    return sqrt(res);
}

double residual_serial(int N, double* u, double *f) {
    double res = 0.0;
    double diff;
    double h = 1.0 / (N + 1);
    double h2 = h * h;
    int i,j, index;
    // first row:
    diff = f[0] + (-4 * u[0] + u[1] + u[N]) / h2;
    res += diff * diff;
    for(j = 1; j < N-1; j++) {
        diff = f[j] + (-4 * u[j] + u[j-1] + u[j+1] + u[j+N]) / h2;
        res += diff * diff;
    }
    diff = f[N-1] + (-4 * u[N-1] + u[N-2] + u[2*N-1]) / h2;
    res += diff * diff;

    // inner rows:
    for(i = 1; i < N - 1; i++){
        index = N * i;
        diff = f[index] + (u[index + 1] + u[index + N] + u[index - N] \
                    - 4 * u[index]) / h2;
        res += diff * diff;
        for(j = 1; j < N - 1; j++) {
            index = N * i + j;
            diff = f[index] + (u[index + 1] + u[index - 1] + u[index + N] \
                        + u[index - N] - 4 * u[index]) / h2;
            res += diff * diff;
        }
        index = N * i + (N - 1);
        diff = f[index] + (u[index - 1] + u[index - N] + u[index + N] \
                    - 4 * u[index]) / h2;
        res += diff * diff;
    }

    // last row:
    index = N*(N-1);
    diff = f[index] + (u[index + 1] + u[index - N] - 4 * u[index]) / h2;
    res += diff * diff;
    for(j = 1; j < N-1; j++) {
        index = N*(N-1) + j;
        diff = f[index] + (u[index - 1] + u[index + 1] + u[index - N] \
                    - 4 * u[index]) / h2;
        res += diff * diff;
    }
    index = N * N - 1;
    diff = f[index] + (u[index - 1] + u[index - N] - 4 * u[index]) / h2;
    res += diff * diff;
    
    return sqrt(res);
}
