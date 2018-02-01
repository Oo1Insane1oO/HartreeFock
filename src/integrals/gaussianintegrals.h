#ifndef GAUSSIANINTEGRALS_H
#define GAUSSIANINTEGRALS_H

#include "../basis/gaussianbasis.h"
#include "../basisfunctions/gaussiancontractedbasis.h"

class GaussianBasis;

class GaussianIntegrals : public GaussianBasis  {
    private:
        unsigned int m_dim;
        double expScaleFactor, sqrtFactor, F0, F1, xScale, sqrtScale, powScale;

        Eigen::ArrayXd normalizationFactors;

        inline double incompleteOverlapIntegral(const unsigned int&);
        inline double incompleteByPartsFactorG(const unsigned int&);
        inline double incompleteByPartsFactorF(const unsigned int&);

        inline double overlapd(const unsigned int&, const unsigned int&); 

        const double& normalizationFactor(const unsigned int&) const;
        
        void setNormalizations();
        void setF0();
        void setF1();

    public:
        GaussianIntegrals(const unsigned int, unsigned int, double=1);
        virtual ~GaussianIntegrals();

        GaussianBasis* getBasis();

        double overlapElement(const unsigned int&, const unsigned int&);
        double kineticElement(const unsigned int&, const unsigned int&);
        double coulombElement(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
};

#endif /* GAUSSIANINTEGRALS_H */
