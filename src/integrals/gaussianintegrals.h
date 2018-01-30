#ifndef GAUSSIANINTEGRALS_H
#define GAUSSIANINTEGRALS_H

#include "../basis/gaussianbasis.h"
#include "../basisfunctions/gaussiancontractedbasis.h"

class GaussianBasis;

class GaussianIntegrals : public GaussianBasis  {
    private:
        unsigned int m_dim;
        double expScaleFactor; 

    public:
        GaussianIntegrals(const unsigned int, unsigned int, double=1);
        virtual ~GaussianIntegrals();

        GaussianBasis* getBasis();

        double element(const unsigned int&, const unsigned int&, const unsigned
                int&, const unsigned int&);
};

#endif /* GAUSSIANINTEGRALS_H */
