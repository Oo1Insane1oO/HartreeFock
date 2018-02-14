#include "gaussianstypeintegrals.h"

GaussianStypeIntegrals::GaussianStypeIntegrals(const unsigned int dim, unsigned
        int numParameters, double scaling) : StypeBasis(dim) {
    /* TODO: fix for general scaling */
    m_dim = dim;
}// end constructor

GaussianStypeIntegrals::~GaussianStypeIntegrals() {
} // end deconstructor

StypeBasis* GaussianStypeIntegrals::getBasis() {
    /* return a pointer to basis class */
    return dynamic_cast<StypeBasis*>(this);
} // end function getBasis

void GaussianStypeIntegrals::initializeParameters(double w) {
    /* TODO: fix for general scaling */
} // end function initializzeParameters
        
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
