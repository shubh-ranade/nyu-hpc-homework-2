all: mmult1 jacobi gs jacobi-no-omp gs-no-omp omp_solved2 omp_solved3 omp_solved4 omp_solved5 omp_solved6 val_test01_solved val_test02_solved

mmult1: MMult1.cpp
	g++ -std=c++11 -O2 -march=native -fopenmp MMult1.cpp -o MMult1

mmult1-no-omp: MMult1.cpp
	g++ -std=c++11 -O3 -march=native MMult1.cpp -o MMult1-no-omp

jacobi: jacobi2D-omp.cpp compute-residual.cpp compute-residual.hpp utils.h
	g++ -std=c++11 -O3 -fopenmp jacobi2D-omp.cpp compute-residual.cpp -o jacobi2D-omp

gs: gs2D-omp.cpp compute-residual.cpp compute-residual.hpp utils.h
	g++ -std=c++11 -O3 -fopenmp gs2D-omp.cpp compute-residual.cpp -o gs2D-omp

jacobi-no-omp: jacobi2D-omp.cpp compute-residual.cpp compute-residual.hpp utils.h
	g++ -std=c++11 -O3 jacobi2D-omp.cpp compute-residual.cpp -o jacobi2D-no-omp

gs-no-omp: gs2D-omp.cpp compute-residual.cpp compute-residual.hpp utils.h
	g++ -std=c++11 -O3 gs2D-omp.cpp compute-residual.cpp -o gs2D-no-omp

omp_solved2: omp_solved2.c
	g++ -std=c++11 -O3 -fopenmp omp_solved2.c -o omp_solved2

omp_solved3: omp_solved3.c
	g++ -std=c++11 -O3 -fopenmp omp_solved3.c -o omp_solved3

omp_solved4: omp_solved4.c
	g++ -std=c++11 -O3 -fopenmp omp_solved4.c -o omp_solved4

omp_solved5: omp_solved5.c
	g++ -std=c++11 -O3 -fopenmp omp_solved5.c -o omp_solved5

omp_solved6: omp_solved6.c
	g++ -std=c++11 -O3 -fopenmp omp_solved6.c -o omp_solved6

val_test01_solved:
	g++ -std=c++11 val_test01_solved.cpp -o val_test01_solved

val_test02_solved:
	g++ -std=c++11 val_test02_solved.cpp -o val_test02_solved

clean:
	rm -rf *.out
	rm jacobi2D-omp
	rm gs2D-omp
	rm jacobi2D-no-omp
	rm gs2D-no-omp
	rm MMult1
	rm omp_solved2
	rm omp_solved3
	rm omp_solved4
	rm omp_solved5
	rm omp_solved6
	rm val_test01_solved
	rm val_test02_solved
