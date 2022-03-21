#include "compute-residual.hpp"
#include "utils.h"

#define MARGIN 1e6

void gs_helper(double*, double*, int);
void gs_serial(double*, double*, int);
void gs_parallel(double*, double*, int);
double gs(double*, int, int);
int nthreads = 1;

int main(int argc, char **argv) {
    if(argc != 3) {
        printf("Usage: %s gridsize max_iterations\n", argv[0]);
        return EXIT_FAILURE;
    }
    const char* run_type;
    if(IS_PARALLEL)
        run_type = "parallel";
    else
        run_type = "serial";

    int N = atoi(argv[1]);
    int max_iter = atoi(argv[2]);
    double *f = (double*) malloc(N * N * sizeof(double));
    int ind;
    double elapsed;
    // timestamp_type t_start, t_end;
	Timer t;
	double time;
    for(ind = 0; ind < N * N; ind++) {
        f[ind] = 1;
    }
    // get_timestamp(&t_start);
	t.tic();
    double res = gs(f, max_iter, N);
    // get_timestamp(&t_end);
	time = t.toc();

    printf("Final gauss-seidel residual: %f\n", res);
    printf("Time taken by %s gauss-seidel method with nthreads %d, gridsize %d, and max_iter %d: %f\n", run_type, nthreads, N, max_iter, time);

    free(f);
    return EXIT_SUCCESS;
}

double gs(double* f, int max_iter, int N) {
    int iter;
    double curr_res, init_res;
    double *u = (double*) calloc(N*N, sizeof(double));
    init_res = compute_residual(N, u, f);
    curr_res = 1;
    for(iter = 0; iter < max_iter && (init_res / curr_res) < MARGIN; iter++) {
        gs_helper(u, f, N);
        curr_res = compute_residual(N, u, f);
    }
    free(u);
    return curr_res;
}

void gs_helper(double* u, double* f, int N) {
    if(IS_PARALLEL)
        gs_parallel(u, f, N);
    else
        gs_serial(u, f, N);
}

void gs_serial(double* u, double* f, int N) {
    // Run a red-black GS method
    int i,j,color, ind;
    double h2 = 1.0 / ((N+1) * (N+1));
    for(color = 0; color < 2; color++) {
        //Handle the corners separately 
        if(color == 0) { //Red
            ind = 0;
            u[ind] = (h2 * f[ind] + u[ind+1] + u[ind+N]) / 4;
            ind = N*N - 1;
            u[ind] = (h2 * f[ind] + u[ind-1] + u[ind-N]) / 4;
            if(N%2 == 1) {
                //All corners red
                ind = N - 1;
                u[ind] = (h2 * f[ind] + u[ind - 1] + u[ind + N]) / 4;
                ind = N * (N-1);
                u[ind] = (h2 * f[ind] + u[ind + 1] + u[ind - N]) / 4;
            }
        } else { //Black
            if(N%2 == 0) {
                //Two black corners
                ind = N - 1;
                u[ind] = (h2 * f[ind] + u[ind - 1] + u[ind + N]) / 4;
                ind = N * (N-1);
                u[ind] = (h2 * f[ind] + u[ind + 1] + u[ind - N]) / 4;
            }
        }
        
        //Handle the top row
        for(ind = 2-color; ind < N - 1; ind+=2) {
            u[ind] = (h2 * f[ind] + u[ind - 1] + u[ind + 1]\
                    + u[ind + N]) / 4;
        }

        //Handle the bottom row
        for(ind = N*(N-1) + 2 - (color+N)%2; ind < N*N-1; ind+=2) {
            u[ind] = (h2 * f[ind] + u[ind - 1] + u[ind + 1] \
                        + u[ind - N]) / 4;
        }

        //Handle the interior rows
        for(i = 1; i < N-1; i++) {
            //Handle the edges
            if(i % 2 == color) {
                //Left edge
                ind = N*i;
                u[ind] = (h2 * f[ind] + u[ind+1] + u[ind+N] + u[ind-N]) / 4;
            }
            if((i + N - 1) % 2 == color) {
                //Right edge
                ind = N *  (i+1) - 1;
                u[ind] = (h2 * f[ind] + u[ind-1] + u[ind+N] + u[ind-N]) / 4;
            }
            for(j = 2 - (i+color)%2; j < N-1; j+=2) {
                ind = N*i + j;
                u[ind] = (h2 * f[ind] + u[ind-1] + u[ind+1] + u[ind+N]\
                            + u[ind-N]) / 4;
            }
        }
    }
}

void gs_parallel(double* u, double* f, int N) {
    // Run a red-black GS method
    int i,j,color, ind;
    double h2 = 1.0 / ((N+1) * (N+1));
    for(color = 0; color < 2; color++) {
        #pragma omp parallel shared(u, f, N, color, h2) private(i,j,ind)
        {
            #ifdef _OPENMP
            if(omp_get_thread_num() == 0)
                nthreads = omp_get_num_threads();
            #endif
            #pragma omp sections nowait
            {
                #pragma omp section
                {
                    //Handle the corners separately 
                    if(color == 0) { //Red
                        ind = 0;
                        u[ind] = (h2 * f[ind] + u[ind+1] + u[ind+N]) / 4;
                        ind = N*N - 1;
                        u[ind] = (h2 * f[ind] + u[ind-1] + u[ind-N]) / 4;
                        if(N%2 == 1) {
                            //All corners red
                            ind = N - 1;
                            u[ind] = (h2 * f[ind] + u[ind - 1] + u[ind + N])\
                                    / 4;
                            ind = N * (N-1);
                            u[ind] = (h2 * f[ind] + u[ind + 1] + u[ind - N])\
                                     / 4;
                        }
                    } else { //Black
                        if(N%2 == 0) {
                            //Two black corners
                            ind = N - 1;
                            u[ind] = (h2 * f[ind] + u[ind - 1] + u[ind + N])\
                                     / 4;
                            ind = N * (N-1);
                            u[ind] = (h2 * f[ind] + u[ind + 1] + u[ind - N])\
                                     / 4;
                        }
                    }
                }
                
                #pragma omp section
                {
                    //Handle the top row
                    for(ind = 2-color; ind < N - 1; ind+=2) {
                        u[ind] = (h2 * f[ind] + u[ind - 1] + u[ind + 1]\
                                + u[ind + N]) / 4;
                    }
                }

                #pragma omp section
                {
                    //Handle the bottom row
                    for(ind = N*(N-1) + 2 - (color+N)%2; ind < N*N-1; ind+=2) {
                        u[ind] = (h2 * f[ind] + u[ind - 1] + u[ind + 1] \
                                    + u[ind - N]) / 4;
                    }
                }
            }

            //Handle the interior rows
            #pragma omp for
            for(i = 1; i < N-1; i++) {
                //Handle the edges
                if(i % 2 == color) {
                    //Left edge
                    ind = N*i;
                    u[ind] = (h2 * f[ind] + u[ind+1] + u[ind+N] + u[ind-N])\
                            / 4;
                }
                if((i + N - 1) % 2 == color) {
                    //Right edge
                    ind = N * (i+1) - 1;
                    u[ind] = (h2 * f[ind] + u[ind-1] + u[ind+N] + u[ind-N])\
                            / 4;
                }
                for(j = 2 - (i+color)%2; j < N-1; j+=2) {
                    ind = N*i + j;
                    u[ind] = (h2 * f[ind] + u[ind-1] + u[ind+1] + u[ind+N]\
                                + u[ind-N]) / 4;
                }
            }
        }
    }
}
