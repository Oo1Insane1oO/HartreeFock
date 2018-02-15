#ifndef GAUSSIANCONTRACTEDBASIS_H
#define GAUSSIANCONTRACTEDBASIS_H

#include <vector>
#include <Eigen/Dense>

#include "gaussianprimitivebasis.h"
#include <memory>

class GaussianContractedBasis {
    private:
        std::vector<std::unique_ptr<GaussianPrimitiveBasis>> m_primitives;

    public:
        GaussianContractedBasis();
        virtual ~GaussianContractedBasis();

        void addPrimitive(const Eigen::VectorXd&, const Eigen::VectorXi&);
        const GaussianPrimitiveBasis* getPrimitive(const unsigned int&) const;
        unsigned int getNumPrimitives();
};

#endif /* GAUSSIANCONTRACTEDBASIS_H */
