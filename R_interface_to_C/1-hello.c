#include "R.h"

// First, there is no main function. The presence of one in either C or C++ programs will 
// create problems in compiling the shared libraries that are needed for using such code in R. 
// For creating a shared library that can be called from R, the function that serves the purpose
// of main in C/C++ code must
// 1. have void as a return type and
// 2. have arguments that are pointers.
void hello(int* n) {
  int i;
  for (i = 0; i < *n; ++i)
    Rprintf("Hi ya'll %d times!\n", i);
}
