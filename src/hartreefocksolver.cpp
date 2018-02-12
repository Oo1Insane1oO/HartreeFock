#include "hartreefocksolver.h"

HartreeFockSolver::HartreeFockSolver(const unsigned int dimension, unsigned int
        cut, const unsigned int numParticles) : Integrals(dimension, cut) {
        m_dim = dimension;
        m_numStates = Integrals::getBasis()->Cartesian::getStates().rows();
        m_numParticles = numParticles;
} // end constructor

HartreeFockSolver::~HartreeFockSolver() {
} // end deconstructor

inline unsigned int HartreeFockSolver::dIndex(const unsigned int& N, const
        unsigned int& i, const unsigned int& j, const unsigned int& k, const
        unsigned int& l) {
    /* calculate offset for 4d-matrix (square case) for indices (i,j,k,l) */
    return i + N * (j + N * (k + N*l));
} // end function dIndex

inline void HartreeFockSolver::assemble() {
    /* assemble integral elements (with symmetries) */

    // array containing elements <ij|1/r_12|ij>_AS 
    twoBodyElements = Eigen::ArrayXd::Zero(m_numStates * m_numStates *
            m_numStates * m_numStates);

    // matrix containing elements <i|h|j>
    overlapElements = Eigen::MatrixXd::Zero(m_numStates, m_numStates);
    oneBodyElements = Eigen::MatrixXd::Zero(m_numStates, m_numStates);

    // set values based on conserved numbers (skip calculation if total angular
    // momentum and spin orthogonality is not satisfied)
    int orbitalp, orbitalr, spinp, spinq, spinr, spins, orbitalSumpq,
        spinSumpq;
    for (unsigned int p = 0; p < m_numStates; ++p) {
        orbitalp = Integrals::getBasis()->Cartesian::getSumn(p);
        spinp = *(Integrals::getBasis()->Cartesian::getStates(p)(m_dim));
        for (unsigned int q = p; q < m_numStates; ++q) {
            overlapElements(p,q) = Integrals::overlapElement(p,q);
            overlapElements(q,p) = overlapElements(p,q);
            oneBodyElements(p,q) = Integrals::kineticElement(p,q) +
            Integrals::potentialElement(p,q);
            oneBodyElements(q,p) = oneBodyElements(p,q);

            orbitalSumpq = orbitalp +
                Integrals::getBasis()->Cartesian::getSumn(q);
            spinq = *(Integrals::getBasis()->Cartesian::getStates(q)(m_dim));
            spinSumpq = spinp + spinq;
            for (unsigned int r = 0; r < m_numStates; ++r) {
                orbitalr = Integrals::getBasis()->Cartesian::getSumn(r);
                spinr = *(Integrals::getBasis()->
                        Cartesian::getStates(r)(m_dim));
                for (unsigned int s = r; s < m_numStates; ++s) {
                    spins = *(Integrals::getBasis()->
                            Cartesian::getStates(s)(m_dim));
                    if ((orbitalSumpq == (orbitalr + Integrals::getBasis()->
                                    Cartesian::getSumn(s))) && (spinSumpq ==
                                (spinr + spins))) {
                        /* make sure sum quantum numbers are conserved */
                        unsigned int idx = dIndex(m_numStates, p,q,r,s);
                        if ((spinp == spinr) && (spinq == spins)) {
                            /* make sure spin orthogonality is satisfied */
                            twoBodyElements(idx) +=
                                Integrals::coulombElement(p, q, r, s);
                        } // end if
                        if ((spinp == spins) && (spinq == spinr)) {
                            /* make sure spin orthogonality is satisfied for
                             * antisymmetric term */
                            twoBodyElements(idx) -=
                                Integrals::coulombElement(p, q, s, r);
                        } // end if
                        twoBodyElements(dIndex(m_numStates, p,q,s,r)) =
                            -twoBodyElements(idx);
                    } // end if
                } // end forp
            } // end forq
        } // end forr
    } // end fors
} // end function assemble

inline void HartreeFockSolver::setDensityMatrix() {
    /* set density matrix in HartreeFock */
    densityMatrix.setZero();
    for (unsigned int c = 0; c < coefficients.rows(); ++c) {
        for (unsigned int d = 0; d < coefficients.cols(); ++d) {
            for (unsigned int i = 0; i < m_numStates; ++i) {
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
            FockMatrix(i,j) += oneBodyElements(i,j); 
            for (unsigned int k = 0; k < m_numStates; ++k) {
                for (unsigned int l = 0; l < m_numStates; ++l) {
                    FockMatrix(i,j) += densityMatrix(k,l) *
                        twoBodyElements(dIndex(m_numStates, i, k, j, l));
                } // end forl
            } // end fork

            // matrix is symmetric by definition
            FockMatrix(j,i) = FockMatrix(i,j);
        } // end forj
    } // end fori
} // end function sethartreeFockMatrix

void HartreeFockSolver::iterate(const unsigned int& maxIterations, const
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
    unsigned int count = 0;
    do {
        /* run for maxIterations or until convergence is reached */

        // set HF-matrix with current coefficients
        setFockMatrix();

        // find eigenvalues and eigenvector (HartreeFock-energies and
        // coefficients respectively)
        eigenSolver.compute(FockMatrix, overlapElements);
        coefficients = eigenSolver.eigenvectors();

        // set density matrix with new coefficients
        setDensityMatrix();

        // update previous energies and increment count
        previousEnergies = eigenSolver.eigenvalues();
        count++;
    } while ((count < maxIterations) && (fabs((eigenSolver.eigenvalues() -
                        previousEnergies).norm()) > eps));

    // find estimate for gound state energy for m_numParticles
    double groundStateEnergy = eigenSolver.eigenvalues().segment(0,
            m_numParticles).sum();
    for (unsigned int a = 0; a < m_numStates; ++a) {
        for (unsigned int b = 0; b < m_numStates; ++b) {
            for (unsigned int c = 0; c < m_numStates; ++c) {
                for (unsigned int d = 0; d < m_numStates; ++d) {
                    groundStateEnergy -= 0.5 * densityMatrix(a,c) *
                        densityMatrix(b,d) *
                        twoBodyElements(dIndex(m_numStates, a,b,c,d));
                } // end ford
            } // end forc
        } // end forb
    } // end fora

    std::cout << "E: " << groundStateEnergy << " Iter: " << count << std::endl;
} // end function iterate
