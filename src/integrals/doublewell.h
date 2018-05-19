#ifndef DOUBLEWELL_H
#define DOUBLEWELL_H

#include <Eigen/Dense>

#include "gaussianintegrals.h"
#include "../basisfunctions/dwc.h"
#include "../hartreefocksolver.h"

class GaussianIntegrals;

class DoubleWell : 
    public HartreeFockSolver<DoubleWell>,
    private DWC,
    private GaussianIntegrals {
    friend class HartreeFockSolver<DoubleWell>;
    using HartreeFockSolver<DoubleWell>::m_dim;
    private:
        unsigned int m_axis, m_numBasis;
        double R, RsqrdFactor;

        double potentialDWElement(const unsigned int&, const unsigned int&);
        double potDWSum(const int&, const int&);

    public:
        using HartreeFockSolver<DoubleWell>::iterate;
        using HartreeFockSolver<DoubleWell>::writeCoefficientsToFile;
        DoubleWell (const unsigned int, unsigned int, const unsigned int);
        virtual ~DoubleWell ();
        
        std::string initializeParameters(double, unsigned int=0);

        double overlapElement(const unsigned int&, const unsigned int&);
        double oneBodyElement(const unsigned int&, const unsigned int&);
        double coulombElement(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
};

#endif /* DOUBLEWELL_H */
