#ifndef DOUBLEWELL_H
#define DOUBLEWELL_H

#include <Eigen/Dense>

#include "gaussianintegrals.h"
#include "../basisfunctions/dwc.h"

class GaussianIntegrals;

class DoubleWell : public GaussianIntegrals, private DWC {
    private:
        unsigned int m_axis, m_numBasis;
        double R, RsqrdFactor;

        Eigen::ArrayXd normalizationFactors;

        double potentialDWElement(const unsigned int&, const unsigned int&);
        double potDWSum(const int&, const int&);

        const double& normalizationFactor(const unsigned int&) const;

        void setNormalizations();

//         double expansionResult(const unsigned int&, const unsigned int&,
//                 double(GaussianIntegrals::*)(const unsigned int&, const
//                     unsigned int&));

    public:
        DoubleWell (const unsigned int, unsigned int);
        virtual ~DoubleWell ();

        unsigned int getSize();
        
        std::string initializeParameters(double, unsigned int=0);

        double overlapElement(const unsigned int&, const unsigned int&);
        double oneBodyElement(const unsigned int&, const unsigned int&);
        double coulombElement(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
};

#endif /* DOUBLEWELL_H */
