#include "hartreefocksolver.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

HartreeFockSolver::HartreeFockSolver(const unsigned int dimension, unsigned int
        cut, const unsigned int numParticles) : Integrals(dimension, cut) {
    /* set dimensions, cutoff and number of particles and initialize basis and
     * integrals */
    m_dim = dimension;
    m_numParticles = numParticles;
    interaction = true;

    // grab info from default communicator 
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
} // end constructor

HartreeFockSolver::~HartreeFockSolver() {
} // end deconstructor

Integrals* HartreeFockSolver::getIntegralObj() {
    /* return a pointer to Integrals */
    return dynamic_cast<Integrals*>(this);
} // end function getIntegralObj

void HartreeFockSolver::setInteraction(bool a) {
    /* set interaction on (if a=true) or false (if a=false) */
    interaction = a;
} // end function setInteraction

inline unsigned int HartreeFockSolver::dIndex(const unsigned int& N, const
        unsigned int& i, const unsigned int& j, const unsigned int& k, const
        unsigned int& l) {
    /* calculate offset for 4d-matrix (square case) for indices (i,j,k,l) */
    return i + N * (j + N * (k + N*l));
} // end function dIndex

inline void HartreeFockSolver::assemble(unsigned int progressDivider) {
    /* assemble integral elements (with symmetries) */

    // class Cartesian (in derived basis) creates full-shell with both
    // spins(each one-body function is created for both spin up and down),
    // divide by two to take the restricted form of the Hartree-Fock equation.
    m_numStates = Integrals::getSize()/2;

    // matrix containing elements <i|h|j> (one-body elements) and elements
    // <i|j> (overlap elements)
    if (myRank == 0) {
        /* only calculate oneBody and overlap elements in root */
        oneBodyElements = Eigen::MatrixXd::Zero(m_numStates, m_numStates);
        overlapElements = Eigen::MatrixXd::Zero(m_numStates, m_numStates);

        // set one-body (uncoupled) elements and overlap elements
        for (unsigned int p = 0; p < m_numStates; ++p) {
            for (unsigned int q = p; q < m_numStates; ++q) {
                overlapElements(p,q) = Integrals::overlapElement(p,q);
                oneBodyElements(p,q) = Integrals::oneBodyElement(p,q);
                if (p != q) {
                    /* only off-diagonal symmetric elements need to be set */
                    overlapElements(q,p) = overlapElements(p,q);
                    oneBodyElements(q,p) = oneBodyElements(p,q);
                } // end if
            } // end forq
        } // end forp
    } // end if

    // set two-body coupled (Coulomb) integral elements
    if (interaction) {
        // create matrix containing pairs (p,q)
        int subSize = m_numStates * (m_numStates+1);
        subSize /= 2;
        Eigen::ArrayXXi pqMap(subSize*subSize,4);
        int j = 0;
        for (unsigned int p = 0; p < m_numStates; ++p) {
            for (unsigned int q = p; q < m_numStates; ++q) {
                for (unsigned int r = 0; r < m_numStates; ++r) {
                    for (unsigned int s = r; s < m_numStates; ++s) {
                        pqMap(j,0) = p;
                        pqMap(j,1) = q;
                        pqMap(j,2) = r;
                        pqMap(j,3) = s;
                        j++;
                    } // end fors
                } // end forr
            } // end forq
        } // end forp
        
        // distribute sizes
        Eigen::ArrayXd pqrsElements, pqsrElements;
        Eigen::ArrayXi displ(numProcs);
        Eigen::ArrayXi sizes(numProcs);
        if (myRank == 0) {
            pqrsElements = Eigen::ArrayXd::Zero(subSize*subSize);
            pqsrElements = Eigen::ArrayXd::Zero(subSize*subSize);
            for (int p = 0; p < numProcs; ++p) {
                sizes(p) = Methods::divider(p, subSize*subSize, numProcs);
                displ(p) = sizes.head(p).sum();
            } // end forp

            // weight the sizes based on the sum of all (p,q,r,s) in process
            Eigen::ArrayXi sums = Eigen::ArrayXi::Zero(numProcs);
            auto findSum = [&]() {
                /* find the sum of all pqrs elements in each process chunk */
                for (unsigned int i = 0; i < sums.size(); ++i) {
                    sums(i) = pqMap.block(displ(i), 0, sizes(i),4).sum();
                } // end fori
            }; // end lambda findSum
            findSum();
            int originalMean = floor(sums.mean());
            // make sure total sum of pqrs within process chunk is within
            // +-1 of originalMean
            for (int i = 0; i < numProcs-1; ++i) {
                for (int j = displ(i); j < displ(i)+sizes(i+1); ++j) {
                    /* iterate over elements in next process */
                    int newSum = sums(i) + pqMap.row(j).sum();
                    if (newSum < originalMean) {
                        /* take pqrs element j from next process */
                        sizes(i) += 1;
                        sizes(i+1) -= 1;
                        displ(i+1) += 1;
                    } else if ((newSum == originalMean-1) || (newSum ==
                                originalMean+1) || (newSum > originalMean)) {
                        /* break if total sum is within +-1 of originalMean
                         * */
                        break;
                    } // end ifeif

                    // find sum of current chunks
                    findSum();
                } // end forj
            } // end fori
        } // end if
        MPI_Bcast(sizes.data(), numProcs, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(displ.data(), numProcs, MPI_INT, 0, MPI_COMM_WORLD);

        // array containing two-body elements <ij|1/r_12|kl> for subset (r,s)
        // of set of range of (p,q) in each process
        Eigen::ArrayXd myTmpTwoBody(sizes(myRank)),
            myTmpTwoBodyAS(sizes(myRank));
        int pqstart = displ(myRank);
        int pqEnd = pqstart+sizes(myRank);
        if (progressDivider) {
            progressDivider = (int)exp(fmod(4.5, pqEnd));
        } // end if
        int rs = 0;
        // save first part of progress bar
        std::string progressPosition, progressBuffer;
        if (progressDivider) {
            progressPosition = Methods::stringPos(myRank, 3) + "Progress: [";
        } // end if
        for (unsigned int pq = pqstart; pq < pqEnd; ++pq) {
            myTmpTwoBody(rs) = Integrals::coulombElement(pqMap(pq,0),
                    pqMap(pq,1), pqMap(pq,2), pqMap(pq,3));
            myTmpTwoBodyAS(rs) = Integrals::coulombElement(pqMap(pq,0),
                    pqMap(pq,1), pqMap(pq,3), pqMap(pq,2));
            
            // print progress
            if (progressDivider) {
                /* show progress if given */
                int progressStart = pq - pqstart;
                if (!(static_cast<int>(fmod(progressStart,
                                    Methods::divider(progressStart, sizes(myRank),
                                        progressDivider))))) {
                    /* print only a few times */
                    progressBuffer = progressPosition;
                    Methods::printProgressBar(progressBuffer,
                            (float)((progressStart==sizes(myRank)-1) ?
                                progressStart : (progressStart+1)) /
                            sizes(myRank), 55, "Two-Body");
                } // end if
            } // end if

            rs++;
        } // end forpq

        // gather subresults from slaves into complete matrix in root
        MPI_Gatherv(myTmpTwoBody.data(), sizes(myRank), MPI_DOUBLE,
                pqrsElements.data(), sizes.data(), displ.data(), MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
        MPI_Gatherv(myTmpTwoBodyAS.data(), sizes(myRank), MPI_DOUBLE,
                pqsrElements.data(), sizes.data(), displ.data(), MPI_DOUBLE, 0,
                MPI_COMM_WORLD);

        // set symmetric values in full two-body matrix used in Hartree-Fock
        // algorithm
        Eigen::ArrayXd tmpTwoBody;
        if (myRank == 0) {
            /* allocate complete matrix for root only */
            tmpTwoBody = Eigen::ArrayXd::Zero(m_numStates * m_numStates *
                    m_numStates * m_numStates);

            int pqrs = 0;
            for (unsigned int p = 0; p < m_numStates; ++p) {
                for (unsigned int q = p; q < m_numStates; ++q) {
                    for (unsigned int r = 0; r < m_numStates; ++r) {
                        for (unsigned int s = r; s < m_numStates; ++s) {
                            double value = pqrsElements(pqrs);
                            tmpTwoBody(dIndex(m_numStates, p,q,r,s)) = value;
                            tmpTwoBody(dIndex(m_numStates, r,q,p,s)) = value;
                            tmpTwoBody(dIndex(m_numStates, r,s,p,q)) = value;
                            tmpTwoBody(dIndex(m_numStates, p,s,r,q)) = value;
                            tmpTwoBody(dIndex(m_numStates, q,p,s,r)) = value;
                            tmpTwoBody(dIndex(m_numStates, s,p,q,r)) = value;
                            tmpTwoBody(dIndex(m_numStates, s,r,q,p)) = value;
                            tmpTwoBody(dIndex(m_numStates, q,r,s,p)) = value;

                            double asvalue = pqsrElements(pqrs);
                            tmpTwoBody(dIndex(m_numStates, p,q,s,r)) = asvalue;
                            tmpTwoBody(dIndex(m_numStates, q,p,r,s)) = asvalue;

                            pqrs++;
                        } // end fors
                    } // end forr
                } // end forq
            } // end forp
    
            // array containing antisymmetric elements <ij|1/r_12|kl>_AS =
            // 2<ij|1/r_12|kl>_- <ij|1/r_12|lk>
            twoBodyElements = Eigen::ArrayXd::Zero(m_numStates * m_numStates *
                    m_numStates * m_numStates);
            for (unsigned int p = 0; p < m_numStates; ++p) {
                for (unsigned int q = 0; q < m_numStates; ++q) {
                    for (unsigned int r = 0; r < m_numStates; ++r) {
                        for (unsigned int s = 0; s < m_numStates; ++s) {
                            twoBodyElements(dIndex(m_numStates, p,q,r,s)) =
                                2*tmpTwoBody(dIndex(m_numStates, p,r,q,s)) -
                                tmpTwoBody(dIndex(m_numStates, p,r,s,q));
                        } // end fors
                    } // end forr
                } // end forq
            } // end forp
        } // end if
    } // end if
} // end function assemble

inline void HartreeFockSolver::setDensityMatrix() {
    /* set density matrix in HartreeFock */
    for (unsigned int c = 0; c < coefficients.rows(); ++c) {
        for (unsigned int d = 0; d < coefficients.cols(); ++d) {
            densityMatrix(c,d) = 0;
            for (unsigned int i = 0; i < m_numParticles/2; ++i) {
                densityMatrix(c,d) += coefficients(c,i) * coefficients(d,i);
            } // end fori
        } // end ford
    } // end forc
} // end function setDensityMatrix

inline void HartreeFockSolver::setFockMatrix() {
    /* set Hartree-Fock matrix */
    FockMatrix.setZero();
    for (unsigned int p = 0; p < m_numStates; ++p) {
        for (unsigned int q = p; q < m_numStates; ++q) {
            FockMatrix(p,q) = oneBodyElements(p,q);
            for (unsigned int r = 0; r < m_numStates; ++r) {
                for (unsigned int s = 0; s < m_numStates; ++s) {
                    FockMatrix(p,q) += densityMatrix(r,s) *
                        twoBodyElements(dIndex(m_numStates, p,q,r,s));
                } // end fors
            } // end forr

            // matrix is symmetric by definition
            FockMatrix(q,p) = FockMatrix(p,q);
        } // end forq
    } // end forp
} // end function sethartreeFockMatrix

double HartreeFockSolver::iterate(const unsigned int& maxIterations, const
        double& eps, const unsigned int progressDivider) {
    /* run Hartree-Fock algorithm for finding coefficients and energy until
     * threshold convergence or until maxIterations is reached */

    // pre-calculate one- and two-body matrix-elements and set initial density
    // matrix with coefficient matrix set to identity
    assemble(progressDivider);
    if (myRank == 0) {
        coefficients = Eigen::MatrixXd::Identity(m_numStates, m_numStates);
        densityMatrix = Eigen::MatrixXd::Zero(m_numStates, m_numStates);
        setDensityMatrix();
        FockMatrix = Eigen::MatrixXd::Zero(m_numStates, m_numStates);
        Eigen::VectorXd previousEnergies = Eigen::VectorXd::Zero(m_numStates);

        // initialize eigenvalue/vector solver for hermitian matrix (Fock
        // matrix is build to be hermitian)
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver;

        // run Hartree-Fock algorithm
        std::string progressPosition, progressBuffer;
        progressPosition = Methods::stringPos(myRank, 3) + "Progress: [";
        for (unsigned int count = 0; count < maxIterations; ++count) {
            /* run for maxIterations or until convergence is reached */

            // set HF-matrix with current coefficients
            setFockMatrix();

            // find eigenvalues and eigenvector (HartreeFock-energies and
            // coefficients respectively)
            eigenSolver.compute(FockMatrix, overlapElements);
            coefficients = eigenSolver.eigenvectors();

            // set density matrix with new coefficients
            setDensityMatrix();

            // check for convergence with RMS of difference between previous
            // and current energies 
            double diff = sqrt((eigenSolver.eigenvalues() -
                        previousEnergies).squaredNorm() / m_numStates);
            if (diff < eps) {
                break;
            } // end if

            // update previous energies
            previousEnergies = eigenSolver.eigenvalues();

            // print progress
            if (progressDivider) {
                /* show progress if given */
                if (!(static_cast<int>(fmod(count, Methods::divider(count,
                                        maxIterations, progressDivider))))) {
                    /* print only a few times */
                    progressBuffer = progressPosition;
                    Methods::printProgressBar(progressBuffer,
                            (float)((count==maxIterations-1) ? count :
                                (count+1)) / maxIterations, 55, "HF");
                } // end if
            } // end if
        } // end forcount

        // find estimate for ground state energy for m_numParticles
        groundStateEnergy = 2*eigenSolver.eigenvalues().segment(0,
                m_numParticles/2).sum();
        for (unsigned int a = 0; a < m_numStates; ++a) {
            for (unsigned int b = 0; b < m_numStates; ++b) {
                for (unsigned int c = 0; c < m_numStates; ++c) {
                    for (unsigned int d = 0; d < m_numStates; ++d) {
                        groundStateEnergy -= densityMatrix(a,c) *
                            densityMatrix(b,d) *
                            twoBodyElements(dIndex(m_numStates, a,c,b,d));
                    } // end ford
                } // end forc
            } // end forb
        } // end fora
        return groundStateEnergy;
    } // end if
    return 0;
} // end function iterate

void HartreeFockSolver::writeCoefficientsToFile(const std::string& filename,
        const std::string& omega) {
    /* write coefficients (in order basis functions were given in integral
     * object) to YAML file */
    std::ofstream outFile(filename + ".yaml");
    YAML::Node info;
    info["E0"] = groundStateEnergy;
    info["omega"] = omega;
    info["dim"] = m_dim;
    info["numbasis"] = m_numStates;
    info["numparticles"] = m_numParticles;
    std::vector<double> tmpCol(m_numStates);
    for (unsigned int p = 0; p < m_numParticles/2; ++p) {
        for (unsigned int i = 0; i < m_numStates; ++i) {
            tmpCol[i] = coefficients(i,p);
        } // end fori
        info["coeffs"].push_back(tmpCol);
    } // end forp
    if (outFile.is_open()) {
        outFile << info;
    } else {
        std::cout << "File : " + filename + "could not be opened." <<
            std::endl;
    } // end ifelse
} // end function writeCoefficientsToFile
