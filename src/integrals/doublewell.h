#ifndef DOUBLEWELL_H
#define DOUBLEWELL_H

#include <Eigen/Dense>

#include "gaussianintegrals.h"
#include "../basisfunctions/dwc.h"

class GaussianIntegrals;
class HartreeFockSolver;

class DoubleWell : public GaussianIntegrals, private DWC {
    private:
        unsigned int m_axis, m_numBasis;
        double R, RsqrdFactor;

        double potentialDWElement(const unsigned int&, const unsigned int&);
        double potDWSum(const int&, const int&);

        void assembleTwoBodyElementsHarmonicOscillator();

    public:
        DoubleWell (HartreeFockSolver*, const unsigned int, unsigned int);
        virtual ~DoubleWell ();

        unsigned int getSize();
        
        std::string initializeParameters(double, unsigned int=0);

        double overlapElement(const unsigned int&, const unsigned int&);
        double oneBodyElement(const unsigned int&, const unsigned int&);
        double coulombElement(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
};

#endif /* DOUBLEWELL_H */
