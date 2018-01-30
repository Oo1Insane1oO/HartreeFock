#include "gaussianintegrals.h"

GaussianIntegrals::GaussianIntegrals(const unsigned int dim, unsigned int
        cutOff, double scaling) : GaussianBasis() {
    m_dim = dim;
    expScaleFactor = scaling;
    GaussianBasis::setup(cutOff, dim, expScaleFactor);
} // end constructor

GaussianIntegrals::~GaussianIntegrals() {
} // end deconstructor

GaussianBasis* GaussianIntegrals::getBasis() {
    /* return a pointer to GaussianBasis */
    // TODO: check that caller actually gets GaussianBasis
    return dynamic_cast<GaussianBasis*>(this);
} // end function getBasis

double GaussianIntegrals::element(const unsigned int& i, const unsigned int& j,
        const unsigned int& k, const unsigned int& l) {
    /* calculate and return the integral element <ij|H|kl> */
    return 0;
} // end function element
