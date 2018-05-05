#ifndef GAUSSIANBASIS_H
#define GAUSSIANBASIS_H

#include "../basisfunctions/cartesian.h"
#include "../methods.h"

class GaussianBasis : public Cartesian {
    private:
        unsigned int m_dim;
        double scaling;

    public:
        GaussianBasis ();
        GaussianBasis (unsigned int, unsigned int);
        virtual ~GaussianBasis ();

        unsigned int getSize();

        void setup(unsigned int, unsigned int);
};

#endif /* GAUSSIANBASIS_H */
