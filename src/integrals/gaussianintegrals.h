#ifndef GAUSSIANINTEGRALS_H
#define GAUSSIANINTEGRALS_H

#include "../basis/gaussianbasis.h"
#include "../hermite/hexpander.h"
#include "../hartreefocksolver.h"

#include <string>

class GaussianBasis;

class GaussianIntegrals : 
    public HartreeFockSolver<GaussianIntegrals>, 
    protected GaussianBasis {
    using HartreeFockSolver<GaussianIntegrals>::m_dim;
    private:
        double expScaleFactor, sqrtFactor, F0, F1, xScale, sqrtScale, powScale,
               sqrtScale1, xScaleHalf, coulomb2DFactor, coulomb3DFactor;

        Eigen::ArrayXd normalizationFactors;
        const double& normalizationFactor(const unsigned int&) const;

        std::unique_ptr<Hexpander> coeffs2D, coeffs3D;
        
        inline double laplacianElement(const unsigned int&, const unsigned
                int&);
        inline double potentialElement(const unsigned int&, const unsigned
                int&);

        inline double ddexprLaplacian(const int&, const int&);
        inline double ddexprPotential(const int&, const int&);

        inline double overlapd(const unsigned int&, const unsigned int&); 

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

    protected:
        unsigned int m_cutOff;
        
        double ddexpr(const int&, const int&,
                double(GaussianIntegrals::*)(const int&, const int&));
        double ddexprOverlap(const int&, const int&);
        void setNormalizations();

    public:
        using HartreeFockSolver<GaussianIntegrals>::iterate;
        using HartreeFockSolver<GaussianIntegrals>::writeCoefficientsToFile;
        GaussianIntegrals(const unsigned int, unsigned int, const unsigned int,
                double=1);
        virtual ~GaussianIntegrals();

        std::string initializeParameters(double);

        double overlapElement(const unsigned int&, const unsigned int&);
        double oneBodyElement(const unsigned int&, const unsigned int&);
        double coulombElement(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
};

#endif /* GAUSSIANINTEGRALS_H */
