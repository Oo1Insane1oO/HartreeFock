#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <mpi.h>

#ifdef TESTS
    #include "../tests/test_main.cpp"
#endif

#include "hartreefocksolver.h"

int main(int argc, char *argv[]) {
    /* main function */
    // let eigen use threads
    Eigen::initParallel();

    #ifdef TESTS
        test_main();
        return 0;
    #endif
} // end main
