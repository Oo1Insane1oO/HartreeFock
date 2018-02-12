#ifndef HARTREEFOCKSOLVER_H
#define HARTREEFOCKSOLVER_H

#include <Eigen/Dense>

#include "integrals/gaussianintegrals.h"

class GaussianIntegrals;
using Integrals = GaussianIntegrals;

class HartreeFockSolver : public Integrals {
    private:
        unsigned int m_dim, m_numStates, m_numParticles;

        Eigen::ArrayXd twoBodyElements;
        Eigen::MatrixXd oneBodyElements;

        Eigen::MatrixXd HartreeFockMatrix, densityMatrix, coefficients;

        inline void assemble();
        
        inline void setDensityMatrix();
        inline void setHartreeFockMatrix();

        inline unsigned int dIndex(const unsigned int&, const unsigned int&,
                const unsigned int&, const unsigned int&, const unsigned int&);

    public:
        HartreeFockSolver (const unsigned int, unsigned int, const unsigned
                int);
        virtual ~HartreeFockSolver ();

        void iterate(const unsigned int&, const double&);
};

#endif /* HARTREEFOCKSOLVER_H */
