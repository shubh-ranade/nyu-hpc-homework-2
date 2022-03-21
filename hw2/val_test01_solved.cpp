# include <cstdlib>
# include <iostream>

using namespace std;

int main ( );
void f ( int n );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST01.
//
//  Discussion:
//
//    TEST01 calls F, which has a memory "leak".  This memory leak can be
//    detected by VALGRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2011
//
//  Bugs:
//    1. the first issue is that the fibonacci loop tries to access x[n] which is out of
//       bounds for an array of size n. Here we get invalid read/write error because cout reads
//       x[n] and we assign x[n] = x[n-1] + x[n-2]. I solved this by changing i <= n to i < n.
//    2. secondly, there is a mismatched free/delete[]/delete error, which arises because we
//       allocate x using malloc be free using delete []. To solve this, we can either allocate using
//       new int[n] or free using free(x) instead of delete [] x.
{
  int n = 10;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  C++ version.\n";
  cout << "  A sample code for analysis by VALGRIND.\n";

  f ( n );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST01\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void f ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    F computes N+1 entries of the Fibonacci sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2011
//
{
  int i;
  int *x;

  x = ( int * ) malloc ( n * sizeof ( int ) ); // error
  // x = new int[n]; // fix mismatched free/delete[]/delete error by initializing x using new int[]

  x[0] = 1;
  cout << "  " << 0 << "  " << x[0] << "\n";

  x[1] = 1;
  cout << "  " << 1 << "  " << x[1] << "\n";

  // for ( i = 2; i <= n; i++ ) // error 1
  for ( i = 2; i < n; i++ )
  {
    x[i] = x[i-1] + x[i-2];
    cout << "  " << i << "  " << x[i] << "\n";
  }

  // delete [] x; // if x is initialized with new int[], use delete []
  free(x); // if x is initialized with malloc, use free()

  return;
}
