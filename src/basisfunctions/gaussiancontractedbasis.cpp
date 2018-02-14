#include "gaussiancontractedbasis.h"

GaussianContractedBasis::GaussianContractedBasis() {
} // end constructor

GaussianContractedBasis::~GaussianContractedBasis() {
} // end deconstructor

void GaussianContractedBasis::addPrimitive(const double weight, const
        Eigen::VectorXi& expVec) {
    /* add primitive */
    m_primitives.push_back(std::make_unique<GaussianPrimitiveBasis>(weight,
                expVec));
} // end function addPrimitives

const GaussianPrimitiveBasis* GaussianContractedBasis::getPrimitive(const
        unsigned int& j) const {
    return m_primitives.at(j).get();
} // end function getPrimitive

unsigned int GaussianContractedBasis::getNumPrimitives() {
    /* get number of primitives (size of basis) */
    return m_primitives.size();
} // end function getNumPrimitives
