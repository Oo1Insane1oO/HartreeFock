#include "gaussianstypeintegrals.h"

GaussianStypeIntegrals::GaussianStypeIntegrals(const unsigned int dim, unsigned
        int numParameters, double scaling) : StypeBasis(dim) {
    m_dim = dim;
}// end constructor

GaussianStypeIntegrals::~GaussianStypeIntegrals() {
} // end deconstructor

StypeBasis* GaussianStypeIntegrals::getBasis() {
    /* return a pointer to basis class */
    return dynamic_cast<StypeBasis*>(this);
} // end function getBasis

void GaussianStypeIntegrals::initializeParameters(const Eigen::VectorXd&
        scalingVector) {
    /* initialize scaling factors and primitives for isotropic case */
    StypeBasis::setPrimitives(scalingVector);
} // end function initializeParameters

void GaussianStypeIntegrals::initializeParameters(const Eigen::MatrixXd&
        scalingMatrix) {
    /* initialize scaling factors and primitives for non-isotropic case */
    StypeBasis::setPrimitives(scalingMatrix);
} // end function initializeParameters
        
double GaussianStypeIntegrals::overlapElement(const unsigned int& i, const
        unsigned int& j) {
    return 0.0;
} // end function overlapElement

double GaussianStypeIntegrals::kineticElement(const unsigned int& i, const
        unsigned int& j) {
    return 0.0;
} // end function kineticElement

double GaussianStypeIntegrals::potentialElement(const unsigned int&, const
        unsigned int&) {
    return 0.0;
} // end function potentialElement

double GaussianStypeIntegrals::coulombElement(const unsigned int& i, const
        unsigned int& j, const unsigned int& k, const unsigned int& l) {
    return 0.0;
}  // end function coulombElement
