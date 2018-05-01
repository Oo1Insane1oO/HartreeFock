#ifndef DOUBLEWELL_H
#define DOUBLEWELL_H

#include <Eigen/Dense>

#include "gaussianintegrals.h"
class GaussianIntegrals;

class DoubleWell : public GaussianIntegrals {
    private:
        unsigned int m_dim, m_axis;
        double R, RsqrdFactor;

        double potentialDWElement(const unsigned int&, const unsigned int&);
        double potDWSum(const int&, const int&);

    public:
        DoubleWell (const unsigned int, unsigned int);
        virtual ~DoubleWell ();
        
        std::string initializeParameters(double, Eigen::MatrixXd, unsigned int=0);

        virtual double overlapElement(const unsigned int&, const unsigned int&)
            = 0;
        double oneBodyElement(const unsigned int&, const unsigned int&);
        virtual double coulombElement(const unsigned int&, const unsigned int&,
                const unsigned int&, const unsigned int&) = 0;
};

#endif /* DOUBLEWELL_H */
