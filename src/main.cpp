#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <chrono>
#include <yaml-cpp/yaml.h>

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
