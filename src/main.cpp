#define NDEBUG
#define EIGEN_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>

#include <Eigen/Dense>
#include <mpi.h>
#include <yaml-cpp/yaml.h>

#ifdef GAUSSHERMITE
    #include "integrals/gaussianintegrals.h"
#endif

#ifdef DOUBLEWELL
    #include "integrals/doublewell.h"
#endif

#ifdef TESTS
    #include "../tests/test_main.cpp"
#endif

#include "hermite/hexpander.h"

YAML::Node argsParser(const char* inputFile) {
    /* parse arguments and set non-optional parameters. return a YAML node with
     * the parsed arguments */

    // create YAML node map
    YAML::Node inputs = YAML::LoadFile(inputFile);

    // check file
    #ifdef GAUSSHERMITE
        if (!inputs["omega"]) {
            std::cout << "Omega missing" << std::endl;
        } // end if
    #endif
    #ifdef DOUBLEWELL
        if (!inputs["R"]) {
            std::cout << "R missing" << std::endl;
        } // end if
    #endif

    if (!inputs["numParticles"] || !inputs["numBasis"] || !inputs["dim"] ||
        !inputs["maxIter"]) {
        /* check that non-optional parameters are given in input file */
        std::cout << "Input file incorrectly setup" << std::endl;
    } // end if

    if (!inputs["progress"]) {
        inputs["progress"] = true;
    } // end if
    if (!inputs["filename"]) {
        inputs["filename"] = "";
    } // end if

    return inputs;
} // end function argsParser

int main(int argc, char *argv[]) {
    /* main function */

    // initialize MPI
    int myRank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    // let eigen use threads
    Eigen::initParallel();
    
    // set inputs
    YAML::Node inputs = argsParser(argv[1]);
    unsigned int progress = (inputs["progress"].as<bool>() ? exp(fmod(4.5,
                    inputs["maxIter"].as<unsigned int>())) : 0);

    #ifdef TESTS
        test_main();
        return 0;
    #endif

    // TODO: fix input file and parser

    // dimensions, cutoff, numParticles
    #ifdef GAUSSHERMITE
        GaussianIntegrals* HFS = new
            GaussianIntegrals(inputs["dim"].as<unsigned int>(),
                    inputs["numBasis"].as<unsigned int>(),
                    inputs["numParticles"].as<unsigned int>());
        std::string message = HFS->initializeParameters(inputs["omega"] .
                as<double>());
        if (message.compare("")) {
            if (myRank == 0) {
                std::cout << message << std::endl;
            } // end if
            delete HFS;
            MPI_Finalize();
            return 0;
        } // end if

        auto start = std::chrono::high_resolution_clock::now();
        double E = HFS->iterate(inputs["maxIter"].as<unsigned int>(), 1e-6,
                progress);
        auto end = std::chrono::high_resolution_clock::now();
        if (myRank == 0) {
            if (inputs["filename"].as<std::string>().compare("")) {
                HFS->writeCoefficientsToFile(inputs["filename"] .
                        as<std::string>(),
                        std::to_string(inputs["omega"].as<double>()));
            } else {
                std::cout << Methods::stringPos(numProcs, 3) << std::endl;
                std::chrono::duration<double> time = end - start;
                std::cout << "Time: " << time.count() << "s" << std::endl;
                std::cout << "Iter: " << HFS->getIterations() << std::endl;
                std::cout << std::setprecision(15) << "E0 = " << E << std::endl;
            } // end ifelse
        } // end if
    #endif
    
//     #ifdef STYPEGAUSSIAN
//         Eigen::VectorXd scalingVec = Eigen::VectorXd::Zero(4);
//         Eigen::MatrixXd centralMatrix = Eigen::MatrixXd::Zero(4,2);
//         scalingVec << 0.25, 0.5, 1.0, 1.5;
//         centralMatrix << 0, 0,
//                          0, 0,
//                          0, 0,
//                          0, 0;
// 
//         HartreeFockSolver* HFS = new HartreeFockSolver(2, 6, 2);
//         HFS->getIntegralObj()->initializeParameters(scalingVec, centralMatrix);
//         double E = HFS->iterate(100, 1e-10);
//     #endif

    #ifdef DOUBLEWELL
        DoubleWell* HFS = new
            DoubleWell(inputs["dim"].as<unsigned int>(),
                    inputs["numBasis"].as<unsigned int>(),
                    inputs["numParticles"].as<unsigned int>());
        std::string message = HFS->initializeParameters(inputs["R"] .
                as<double>());
        if (message.compare("")) {
            if (myRank == 0) {
                std::cout << message << std::endl;
            } // end if
            delete HFS;
            MPI_Finalize();
            return 0;
        } // end if
        
        auto start = std::chrono::high_resolution_clock::now();
        double E = HFS->iterate(inputs["maxIter"].as<unsigned int>(), 1e-8,
                progress);
        auto end = std::chrono::high_resolution_clock::now();
        if (myRank == 0) {
            if (inputs["filename"].as<std::string>().compare("")) {
                HFS->writeCoefficientsToFile(inputs["filename"] .
                        as<std::string>(),
                        std::to_string(inputs["R"].as<double>()));
            } // end if
            std::cout << Methods::stringPos(numProcs, 3) << std::endl;
            std::chrono::duration<double> time = end - start;
            std::cout << "Time: " << time.count() << "s" << std::endl;
            std::cout << std::setprecision(15) << "E0 = " << E << std::endl;
        } // end if
    #endif

    // free HFS object and end MPI
    delete HFS;
    MPI_Finalize();

    return 0;
} // end main
