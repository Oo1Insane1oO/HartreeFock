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
        Eigen::MatrixXd oneBodyElements, overlapElements;

        Eigen::MatrixXd FockMatrix, densityMatrix, coefficients;

        inline void assemble();
        
        inline void setDensityMatrix();
        inline void setFockMatrix();

        inline unsigned int dIndex(const unsigned int&, const unsigned int&,
                const unsigned int&, const unsigned int&, const unsigned int&);

    public:
        HartreeFockSolver (const unsigned int, unsigned int, const unsigned
                int);
        virtual ~HartreeFockSolver ();

        Integrals* getIntegralObj();

        double iterate(const unsigned int&, const double&);
};

#endif /* HARTREEFOCKSOLVER_H */
