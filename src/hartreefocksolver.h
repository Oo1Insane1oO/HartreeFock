#ifndef HARTREEFOCKSOLVER_H
#define HARTREEFOCKSOLVER_H

#include <Eigen/Dense>

#ifdef GAUSSHERMITE
    #include "integrals/gaussianintegrals.h"

    class GaussianIntegrals;
    using Integrals = GaussianIntegrals;
#endif

#ifdef STYPEGAUSSIAN
    #include "integrals/gaussianstypeintegrals.h"

    class GaussianStypeIntegrals;
    using Integrals = GaussianStypeIntegrals;
#endif

class HartreeFockSolver : public Integrals {
    private:
        unsigned int m_dim, m_numStates, m_numParticles;

        bool interaction;

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

        void setInteraction(bool);
};

#endif /* HARTREEFOCKSOLVER_H */
