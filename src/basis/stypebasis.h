#ifndef STYPEBASIS_H
#define STYPEBASIS_H

#include <Eigen/Dense>

#include "../basisfunctions/gaussiancontractedbasis.h"

class StypeBasis : public GaussianContractedBasis {
    private:
        unsigned int m_dim;

    public:
        StypeBasis (unsigned int);
        virtual ~StypeBasis ();

        unsigned int getSize();
        
        void setPrimitives(const double);
};

#endif /* STYPEBASIS_H */
