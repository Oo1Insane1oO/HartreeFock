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

    // write basis to file and exit
    #ifdef GENERATEBASIS
        Cartesian* basis = new Cartesian();
        basis->setup(std::atoi(argv[1]), std::atoi(argv[2]));
        basis->restructureStates();
        basis->writeToFile(argv[3]);
        delete basis;
        exit(0);
    #endif

    #ifdef TESTS
        test_main();
        return 0;
    #endif

    // TODO: fix input file and parser

    // dimensions, cutoff, numParticles
    double w = 0.5;
    HartreeFockSolver* HFS = new HartreeFockSolver(2, 2, 2);
    HFS->getIntegralObj()->initializeParameters(w);
    double E = HFS->iterate(100, 1e-8);

    std::cout << E << std::endl;

    delete HFS;
} // end main
