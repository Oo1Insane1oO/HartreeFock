#ifndef GAUSSIANINTEGRALS_H
#define GAUSSIANINTEGRALS_H

#include "../basis/gaussianbasis.h"
#include "../hermite/hexpander.h"

class GaussianBasis;

class GaussianIntegrals : public GaussianBasis  {
    private:
        unsigned int m_dim;
        double expScaleFactor, sqrtFactor, F0, F1, xScale, sqrtScale, powScale,
               sqrtScale1, xScaleHalf, coulomb2DFactor, coulomb3DFactor;

        Eigen::ArrayXd normalizationFactors;

        std::unique_ptr<Hexpander> coeffs;
        
        inline double laplacianElement(const unsigned int&, const unsigned
                int&);
        inline double potentialElement(const unsigned int&, const unsigned
                int&);

        inline double ddexpr(const int&, const int&,
                double(GaussianIntegrals::*)(const int&, const int&));
        inline double ddexprOverlap(const int&, const int&);
        inline double ddexprLaplacian(const int&, const int&);
        inline double ddexprPotential(const int&, const int&);

        inline double overlapd(const unsigned int&, const unsigned int&); 

        const double& normalizationFactor(const unsigned int&) const;

        double (GaussianIntegrals::*coulombFunc)(const unsigned int&,
                const unsigned int&, const unsigned int&, const unsigned int&);
        inline double coulombElement2D(const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&);
        inline double coulombElement3D(const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&, const unsigned int&);
        inline double coulomb2D(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
        inline double coulomb3D(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
        
        void setNormalizations();

    public:
        GaussianIntegrals(const unsigned int, unsigned int, double=1);
        virtual ~GaussianIntegrals();

        GaussianBasis* getBasis();

        void initializeParameters(double);

        double overlapElement(const unsigned int&, const unsigned int&);
        double oneBodyElement(const unsigned int&, const unsigned int&);
        double coulombElement(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
};

#endif /* GAUSSIANINTEGRALS_H */
