#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <mpi.h>

#ifdef TESTS
    #include "../tests/test_main.cpp"
#endif

#include "hartreefocksolver.h"
#include "hermite/hexpander.h"

int main(int argc, char *argv[]) {
    /* main function */
    // let eigen use threads
    Eigen::initParallel();
// 
//     unsigned int dim = 3;
//     Hexpander* h = new Hexpander(6, dim, 1, 1, Eigen::VectorXd::Zero(dim),
//             Eigen::VectorXd::Zero(dim));
//     std::cout << h->coeff(0,1,0,0) << std::endl;
//     delete h;
//     exit(1);

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
    #ifdef GAUSSHERMITE
        double w = 1.0;
        HartreeFockSolver* HFS = new HartreeFockSolver(2, 42, 6);
//         HartreeFockSolver* HFS = new HartreeFockSolver(3, 2, 2);
        HFS->getIntegralObj()->initializeParameters(w);
        double E = HFS->iterate(100, 1e-8);
    #endif
    
    #ifdef STYPEGAUSSIAN
        Eigen::VectorXd scalingVec = Eigen::VectorXd::Zero(4);
        Eigen::MatrixXd centralMatrix = Eigen::MatrixXd::Zero(4,2);
        scalingVec << 0.25, 0.5, 1.0, 1.5;
        centralMatrix << 0, 0,
                         0, 0,
                         0, 0,
                         0, 0;

        HartreeFockSolver* HFS = new HartreeFockSolver(2, 6, 2);
        HFS->getIntegralObj()->initializeParameters(scalingVec, centralMatrix);
        double E = HFS->iterate(100, 1e-8);
    #endif

    std::cout << "E0 = " << E << std::endl;

    delete HFS;
} // end main
