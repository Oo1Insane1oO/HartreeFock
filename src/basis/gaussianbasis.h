#ifndef GAUSSIANBASIS_H
#define GAUSSIANBASIS_H

#include "../basisfunctions/cartesian.h"
#include "../methods.h"
#include "../basisfunctions/gaussiancontractedbasis.h"

class GaussianBasis : public Cartesian, public GaussianContractedBasis {
    private:
        unsigned int m_dim;

    public:
        GaussianBasis ();
        GaussianBasis (unsigned int, unsigned int, double);
        virtual ~GaussianBasis ();

        void setup(unsigned int, unsigned int, double=1);
        void setPrimitives(const double scaling);

        const Eigen::Ref<const Eigen::VectorXi> getExpVec(const unsigned int&)
            const;
};

#endif /* GAUSSIANBASIS_H */
