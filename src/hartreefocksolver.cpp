#include "hartreefocksolver.h"

#include <iostream>

HartreeFockSolver::HartreeFockSolver(const unsigned int dimension, unsigned int
        cut, const unsigned int numParticles) : Integrals(dimension, cut) {
    /* set dimensions, cutoff and number of particles and initialize basis and
     * integrals */
    m_dim = dimension;
    m_numParticles = numParticles;
    interaction = true;
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

inline void HartreeFockSolver::assemble() {
    /* assemble integral elements (with symmetries) */
    m_numStates = Integrals::getBasis()->getSize()/2;
//     m_numStates = (m_numStates%2==0 ? m_numStates/2 : (m_numStates+1)/2);

    // array containing elements <ij|1/r_12|ij>_AS 
    twoBodyElements = Eigen::ArrayXd::Zero(m_numStates * m_numStates *
            m_numStates * m_numStates);

    // matrix containing elements <i|h|j>
    overlapElements = Eigen::MatrixXd::Zero(m_numStates, m_numStates);
    oneBodyElements = Eigen::MatrixXd::Zero(m_numStates, m_numStates);

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

    // set two-body coupled (Coulomb) integral elements
    for (unsigned int p = 0; p < m_numStates; ++p) {
        for (unsigned int q = p; q < m_numStates; ++q) {
            for (unsigned int r = 0; r < m_numStates; ++r) {
                for (unsigned int s = r; s < m_numStates; ++s) {
                    twoBodyElements(dIndex(m_numStates, p,q,r,s)) =
                        Integrals::coulombElement(p,q,r,s);
                } // end fors
            } // end forq
        } // end forr
    } // end forp

    // set remaning symmetric elements
    for (unsigned int p = 0; p < m_numStates; ++p) {
        for (unsigned int r = 0; r < m_numStates; ++r) {
            for (unsigned int q = p; q < m_numStates; ++q) {
                for (unsigned int s = r; s < m_numStates; ++s) {
                    double value = twoBodyElements(dIndex(m_numStates,
                                p,r,q,s));
                    twoBodyElements(dIndex(m_numStates, q,s,p,r)) = value;
                    twoBodyElements(dIndex(m_numStates, q,r,p,s)) = value;
                    twoBodyElements(dIndex(m_numStates, p,s,q,r)) = value;
                    twoBodyElements(dIndex(m_numStates, r,p,s,q)) = value;
                    twoBodyElements(dIndex(m_numStates, s,p,r,q)) = value;
                    twoBodyElements(dIndex(m_numStates, r,q,s,p)) = value;
                    twoBodyElements(dIndex(m_numStates, s,q,r,p)) = value;
                } // end fors
            } // end forq
        } // end forr
    } // end forp
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
    for (unsigned int i = 0; i < m_numStates; ++i) {
        for (unsigned int j = i; j < m_numStates; ++j) {
            FockMatrix(i,j) = oneBodyElements(i,j);
            for (unsigned int k = 0; k < m_numStates; ++k) {
                for (unsigned int l = 0; l < m_numStates; ++l) {
                    FockMatrix(i,j) += densityMatrix(k,l) *
                        (2*twoBodyElements(dIndex(m_numStates, i,k,j,l)) -
                         twoBodyElements(dIndex(m_numStates, i,k,l,j)));
                } // end forl
            } // end fork

            // matrix is symmetric by definition
            FockMatrix(j,i) = FockMatrix(i,j);
        } // end forj
    } // end fori
} // end function sethartreeFockMatrix

double HartreeFockSolver::iterate(const unsigned int& maxIterations, const
        double& eps) {
    /* run Hartree-Fock algorithm for finding coefficients and energy until
     * threshold convergence or until maxIterations is reached */

    // pre-calculate one- and two-body matrix-elements and set initial density
    // matrix with coefficient matrix set to identity
    assemble();
    coefficients = Eigen::MatrixXd::Identity(m_numStates, m_numStates);
    densityMatrix = Eigen::MatrixXd::Zero(m_numStates, m_numStates);
    setDensityMatrix();
    FockMatrix = Eigen::MatrixXd::Zero(m_numStates, m_numStates);
    Eigen::VectorXd previousEnergies = Eigen::VectorXd::Zero(m_numStates);

    // initialize eigenvalue/vector solver for hermitian matrix (Fock matrix is
    // build to be hermitian)
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver;

    // run Hartree-Fock algorithm
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

        // check for convergence with RMS of difference between previous and
        // current energies 
        double diff = sqrt((eigenSolver.eigenvalues() -
                    previousEnergies).squaredNorm() / m_numStates);
        if (diff < eps) {
            break;
        } // end if

        // update previous energies
        previousEnergies = eigenSolver.eigenvalues();
    } // end forcount

    std::cout << FockMatrix << std::endl;

    // find estimate for ground state energy for m_numParticles
    double groundStateEnergy = 2*eigenSolver.eigenvalues().segment(0,
            m_numParticles/2).sum();
    for (unsigned int a = 0; a < m_numStates; ++a) {
        for (unsigned int b = 0; b < m_numStates; ++b) {
            for (unsigned int c = 0; c < m_numStates; ++c) {
                for (unsigned int d = 0; d < m_numStates; ++d) {
                    groundStateEnergy -= densityMatrix(a,c) *
                        densityMatrix(b,d) *
                        (2*twoBodyElements(dIndex(m_numStates, a,b,c,d)) -
                         twoBodyElements(dIndex(m_numStates, a,b,d,c)));
                } // end ford
            } // end forc
        } // end forb
    } // end fora

    return groundStateEnergy;
} // end function iterate
