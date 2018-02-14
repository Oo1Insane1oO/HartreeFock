#include "stypebasis.h"

StypeBasis::StypeBasis(unsigned int dim) : GaussianContractedBasis() {
    m_dim = dim;
} // end constructor

StypeBasis::~StypeBasis() {
} // end deconstructor

unsigned int StypeBasis::getSize() {
    return GaussianContractedBasis::getNumPrimitives();
} // end function getSize

void StypeBasis::setPrimitives(const double scaling) {
    /* create and set primitive functions in contracted basis */
} // end function setPrimitives
