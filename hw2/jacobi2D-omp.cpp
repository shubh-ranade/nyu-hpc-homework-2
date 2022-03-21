#include "compute-residual.hpp"
#include "utils.h"

#define MARGIN 1e6

void run_jacobi_serial(double*, double*, double*, int);
void run_jacobi_parallel(double*, double*, double*, int);
double jacobi(double*, int, int);
void jacobi_helper(double*, double*, double*, int);
int nthreads = 1;

int main(int argc, char **argv){
    if(argc != 3) {
        printf("Usage: %s gridsize max_iter\n", argv[0]);
        return EXIT_SUCCESS;
    }
    const char* run_type;
    if(IS_PARALLEL)
        run_type = "parallel";
    else
        run_type = "serial";
    
    int i,j;
    int N = atoi(argv[1]);
    int max_iter = atoi(argv[2]);
    double *f = (double*) malloc(N * N * sizeof(double));
    
    // initialize f
    for(i = 0; i < N; i++) {
        for(j = 0; j < N; j++) {
            f[N*i+j] = 1;
        }
    }

    Timer t;
    t.tic();
    double res = jacobi(f, N, max_iter);
    double time = t.toc();
    
    printf("Final jacobi residual: %f\n", res);
    printf("Time taken by %s jacobi method with nthreads %d, gridsize %d, and max_iter %d: %f\n", run_type, nthreads, N, max_iter, time);
    
    free(f);
    return EXIT_SUCCESS;
}

void jacobi_helper(double *u, double *u_new, double *f, int N) {
    if(IS_PARALLEL)
        run_jacobi_parallel(u, u_new, f, N);
    else
        run_jacobi_serial(u, u_new, f, N);
}

void run_jacobi_serial(double *u, double *u_new, double *f, int N) {
    int i,j, ind;
    double h = 1.0 / (N + 1);
    double h2 = h * h;
    u_new[0] = (h2 * f[0] + u[1] + u[N]) / 4.0;
    for(j = 1; j < N-1; j++) {
        u_new[j] = (h2 * f[j] + u[j-1] + u[j+1] + u[j+N]) / 4.0;
    }
    u_new[N-1] = (h2 * f[N-1] + u[N-2] + u[2*N-1]) / 4.0;

    for(i = 1; i < N - 1; i++) {
        u_new[N*i] = (h2 * f[N*i] + u[N * (i - 1)] + u[N * (i + 1)]\
                     + u[N * i + 1]) / 4.0;
        for(j = 1; j < N - 1; j++) {
            ind = N * i + j;
            u_new[ind] = (h2 * f[ind] + u[ind + 1] + u[ind - 1]\
                                + u[ind + N] + u[ind - N]) / 4.0 ;
        }
        u_new[N*(i + 1) - 1] = (h2 * f[N*(i + 1) - 1] \
                    + u[N*(i + 1) - 2]  + u[N*i - 1] + u[N*(i + 2) - 1]) / 4.0;
    }

    u_new[N*(N-1)] = (h2 * f[N*(N-1)] + u[N * (N-1) + 1] \
                            + u[N*(N-2)]) / 4.0;
    for(j = 1; j < N-1; j++) {
        ind = N*(N-1) + j;
        u_new[ind] = (h2 * f[ind] + u[ind-1] + u[ind+1] \
                            + u[ind - N]) / 4.0;
    }
    u_new[N*N-1] = (h2 * f[N*N-1] + u[N*N-2] + u[N*(N-1)-1]) / 4.0;
}

void run_jacobi_parallel(double *u, double *u_new, double* f, int N) {
    int i,j,ind;
    double h = 1.0 / (N+1);
    double h_sqr = h * h;
    #pragma omp parallel private(i,j,ind) shared(u, u_new, h, h_sqr, N, f)
    {
        #ifdef _OPENMP
        if(omp_get_thread_num() == 0)
            nthreads = omp_get_num_threads();
        #endif
        #pragma omp sections nowait
        {
            #pragma omp section
            {
                //First row
                u_new[0] = 0.25 * (h_sqr * f[0] + u[1] + u[N]);
                for(j = 1; j < N-1; j++) {
                    u_new[j] = 0.25 * (h_sqr * f[j] + u[j-1] + u[j+1] + u[j+N]);
                }
                u_new[N-1] = 0.25 * (h_sqr * f[N-1] + u[N-2] + u[2*N-1]);
            }

            #pragma omp section
            {
                //Last row
                ind = N * (N-1);
                u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind + 1] + u[ind - N]);
                for(j = 1; j < N-1; j++) {
                    ind = N*(N-1) + j;
                    u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind - 1] \
                                            + u[ind + 1] + u[ind - N]);
                }
                ind = N*N - 1;
                u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind - 1]\
                                + u[ind - N]);
            }
        }

        #pragma omp for
        for(i=1; i < N-1; i++) {
            ind = N*i;
            u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind - N] + u[ind + N] \
                                + u[ind + 1]);
            for(j=1; j<N-1; j++) {
                ind = N*i + j;
                u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind-1] \
                                    + u[ind+1] + u[ind-N] + u[ind+N]);
            }
            ind = N * (i + 1) - 1;
            u_new[ind] = 0.25 * (h_sqr * f[ind] + u[ind - 1] + u[ind - N]\
                                + u[ind + N]);
        }
    }
}

double jacobi(double *f, int N, int max_iter) {
    int iter_num;
    double *u, *u_2;
    double curr_res, init_res;

    // initialize u and u_2
    u = (double*) malloc(N * N * sizeof(double));
    u_2 = (double*) malloc(N * N * sizeof(double));

    //Compute the initial residual
    init_res = compute_residual(N, u, f);
    curr_res = 1.0;
    
    for(iter_num = 0; iter_num < max_iter && (init_res / curr_res) < MARGIN; iter_num++) {
        
        jacobi_helper(u, u_2, f, N);
        double *swp_tmp = u;
        u = u_2;
        u_2 = swp_tmp;

        curr_res = compute_residual(N, u, f);
        // printf("Jacobi iter: %04d norm of residual: %f\n", iter_num+1, curr_res);
    }
    
    free(u);
    free(u_2);
    return curr_res;
}
