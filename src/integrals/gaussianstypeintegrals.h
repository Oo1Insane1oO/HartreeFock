#ifndef GAUSSIANSTYPEINTEGRALS_H
#define GAUSSIANSTYPEINTEGRALS_H

#include "../basis/stypebasis.h"

class StypeBasis;

class GaussianStypeIntegrals : public StypeBasis {
    private:
        unsigned int m_dim;
        double sqrt2pi, sqrtpi;
        bool isotropic; 
        Eigen::ArrayXd normalizationFactors;

        void setNormalizations();

        double overlapd(const unsigned int&, const unsigned int&);

    public:
        GaussianStypeIntegrals (const unsigned int, const unsigned int,
                double=1);
        virtual ~GaussianStypeIntegrals ();

        StypeBasis* getBasis();

        void initializeParameters(const Eigen::VectorXd&);
        void initializeParameters(const Eigen::MatrixXd&);
        
        double overlapElement(const unsigned int&, const unsigned int&);
        double kineticElement(const unsigned int&, const unsigned int&);
        double coulombElement(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
        double potentialElement(const unsigned int&, const unsigned int&);
};

#endif /* GAUSSIANSTYPEINTEGRALS_H */
