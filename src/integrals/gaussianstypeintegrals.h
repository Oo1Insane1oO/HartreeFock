#ifndef GAUSSIANSTYPEINTEGRALS_H
#define GAUSSIANSTYPEINTEGRALS_H

#include "../basis/stypebasis.h"

class StypeBasis;

class GaussianStypeIntegrals : public StypeBasis {
    private:
        unsigned int m_dim;
        double sqrt2pi, sqrtpi, sqrtpiDim;
        bool isotropic; 
        Eigen::ArrayXd normalizationFactors;

        void setNormalizations();

        double overlapd(const unsigned int&, const unsigned int&);
        
        inline double laplacianElement(const unsigned int&, const unsigned
                int&);
        inline double potentialElement(const unsigned int&, const unsigned
                int&);
        inline double coulombElement2D(const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&);
        inline double coulombElement3D(const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&);
        double (GaussianStypeIntegrals::*coulombElementFunc)(const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&);

    public:
        GaussianStypeIntegrals (const unsigned int, const unsigned int,
                double=1);
        virtual ~GaussianStypeIntegrals ();

        StypeBasis* getBasis();

        void initializeParameters(const Eigen::VectorXd&, const
                Eigen::MatrixXd&);
        void initializeParameters(const Eigen::MatrixXd&, const
                Eigen::MatrixXd&);
        
        double overlapElement(const unsigned int&, const unsigned int&);
        double oneBodyElement(const unsigned int&, const unsigned int&);
        double coulombElement(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
};

#endif /* GAUSSIANSTYPEINTEGRALS_H */
