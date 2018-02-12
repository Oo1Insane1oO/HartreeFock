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

    // set size of integral array
    integralElements = Eigen::ArrayXd::Zero(m_numStates * m_numStates *
            m_numStates * m_numStates);

    // set values based on conserved numbers (skip calculation if total angular
    // momentum and spin orthogonality is not satisfied)
    int orbitalp, orbitalr, spinp, spinq, spinr, spins, orbitalSumpq,
        spinSumpq;
    double overlapTermAndkineticTerm;
    for (unsigned int p = 0; p < m_numStates; ++p) {
        orbitalp = Integrals::getBasis()->Cartesian::getSumn(p);
        spinp = *(Integrals::getBasis()->Cartesian::getStates(p)(m_dim));
        for (unsigned int q = p; q < m_numStates; ++q) {
            orbitalSumpq = orbitalp +
                Integrals::getBasis()->Cartesian::getSumn(q);
            spinq = *(Integrals::getBasis()->Cartesian::getStates(q)(m_dim));
            spinSumpq = spinp + spinq;
            overlapTermAndkineticTerm = Integrals::overlapElement(p,q) +
                Integrals::kineticElement(p,q);
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
                        integralElements(idx) = overlapTermAndkineticTerm;
                        if ((spinp == spinr) && (spinq == spins)) {
                            /* make sure spin orthogonality is satisfied */
                            integralElements(idx) +=
                            Integrals::coulombElement(p, q, r, s);
                        } // end if
                        if ((spinp == spins) && (spinq == spinr)) {
                            /* make sure spin orthogonality is satisfied for
                             * antisymmetric term */
                            integralElements(idx) -=
                            Integrals::coulombElement(p, q, s, r);
                        } // end if
                        integralElements(dIndex(m_numStates, p,q,s,r)) =
                            -integralElements(idx);
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

inline void HartreeFockSolver::setHartreeFockMatrix() {
    /* set Hartree-Fock matrix */
    HartreeFockMatrix.setZero();
    for (unsigned int i = 0; i < m_numStates; ++i) {
        for (unsigned int j = i; j < m_numStates; ++j) {
            for (unsigned int k = 0; k < m_numStates; ++k) {
                for (unsigned int l = 0; l < m_numStates; ++l) {
                    HartreeFockMatrix(i,j) += densityMatrix(k,l) *
                        integralElements(dIndex(m_numStates, i, k, j, l));
                } // end forl
            } // end fork

            // matrix is symmetric by definition
            HartreeFockMatrix(j,i) = HartreeFockMatrix(i,j);
        } // end forj
    } // end fori
} // end function sethartreeFockMatrix

void HartreeFockSolver::iterate(const unsigned int& maxIterations, const
        double& eps) {
    /* run Hartree-Fock algorithm for finding coefficients and energy until
     * threshold convergence or until maxIterations is reached */

    // pre-calculate matrix-elements, allocate coefficient matrix and
    // Hartree-Fock matrix and set density matrix with unity
    assemble();
    coefficients = Eigen::MatrixXd::Identity(m_numStates, m_numStates);
    densityMatrix = Eigen::MatrixXd::Zero(m_numStates, m_numStates);
    setDensityMatrix();
    HartreeFockMatrix = Eigen::MatrixXd::Zero(m_numStates, m_numStates);
    Eigen::VectorXd previousEnergies = Eigen::VectorXd::Zero(m_numStates);

    // set eigen solver for matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver;

//     for (unsigned int i = 0; i < integralElements.size(); ++i) {
//         if (integralElements(i) != 0) {
//             std::cout << integralElements(i) << std::endl;
//         }
//     }

    // run Hartree-Fock algorithm
    unsigned int count = 0;
    double diff = 1. + eps;
    while (count < maxIterations && diff > eps) {
        /* run for maxIterations or until convergence is reached */

        // set HF-matrix with current coefficients
        setHartreeFockMatrix();

        // find eigenvalues and eigenvector (HF-energies and coefficients
        // respectively)
        eigenSolver.compute(HartreeFockMatrix);
        coefficients = eigenSolver.eigenvectors();

        // set density matrix with new coefficients
        setDensityMatrix();

        // set difference between previes and current values for convergence
        // test
        diff = fabs((eigenSolver.eigenvalues() - previousEnergies).mean());

        // update previous energies and increment count
        previousEnergies = eigenSolver.eigenvalues();
        count++;
    } // end while count

    // find estimate for gound state energy for m_numParticles
    double groundStateEnergy = eigenSolver.eigenvalues().segment(0,
            m_numParticles).sum();
    for (unsigned int a = 0; a < m_numStates; ++a) {
        for (unsigned int b = 0; b < m_numStates; ++b) {
            for (unsigned int c = 0; c < m_numStates; ++c) {
                for (unsigned int d = 0; d < m_numStates; ++d) {
                    groundStateEnergy -= 0.5 * densityMatrix(a,c) *
                        densityMatrix(b,d) *
                        integralElements(dIndex(m_numStates, a,b,c,d));
                } // end ford
            } // end forc
        } // end forb
    } // end fora

    std::cout << "E: " << groundStateEnergy << " Iter: " << count << std::endl;
} // end function iterate
