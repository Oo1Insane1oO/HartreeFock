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

        GaussianContractedBasis* getBasis();
        
        void setPrimitives(const double);
        void setPrimitives(const Eigen::VectorXd&, const Eigen::MatrixXd&);
        void setPrimitives(const Eigen::MatrixXd&, const Eigen::MatrixXd&);
};

#endif /* STYPEBASIS_H */
