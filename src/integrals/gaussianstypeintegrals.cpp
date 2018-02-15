#include "gaussianstypeintegrals.h"

#include <boost/math/special_functions/erf.hpp>
#include <iostream>

GaussianStypeIntegrals::GaussianStypeIntegrals(const unsigned int dim, unsigned
        int numParameters, double scaling) : StypeBasis(dim) {
    m_dim = dim;
    sqrt2pi = sqrt(2*M_PI);
    sqrtpi = sqrt(M_PI);
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
    setNormalizations();
    isotropic = true;
} // end function initializeParameters

void GaussianStypeIntegrals::initializeParameters(const Eigen::MatrixXd&
        scalingMatrix) {
    /* initialize scaling factors and primitives for non-isotropic case */
    StypeBasis::setPrimitives(scalingMatrix);
    setNormalizations();
    isotropic = false;
} // end function initializeParameters

void GaussianStypeIntegrals::setNormalizations() {
    /* calculate and set normalization factors for all basis functions */
    normalizationFactors =
        Eigen::ArrayXd::Zero(StypeBasis::getBasis()->getNumPrimitives());
    for (unsigned int i = 0; i < StypeBasis::getBasis()->getNumPrimitives();
            ++i) {
        double prod = 1.0;
        for (unsigned int d = 0; d < m_dim; ++d) {
            prod *= overlapd(i,d) * overlapd(i,d);
        } // end ford
        normalizationFactors(i) = 1.0 / prod;
    } // end fori
    normalizationFactors = normalizationFactors.cwiseSqrt();
} // end function setNormalizations

double GaussianStypeIntegrals::overlapd(const unsigned int& i, const unsigned
        int& d) {
    /* calculate and return the 1D overlap element for state i */
    return sqrt2pi / sqrt(StypeBasis::getBasis()->getPrimitive(i)->scaling(d));
} // end function overlapd
        
double GaussianStypeIntegrals::overlapElement(const unsigned int& i, const
        unsigned int& j) {
    double prod = 1.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        prod *= sqrt2pi /
            sqrt(StypeBasis::getBasis()->getPrimitive(i)->scaling(d) +
                    StypeBasis::getBasis()->getPrimitive(j)->scaling(d));
    } // end ford
    return prod * normalizationFactors(i) * normalizationFactors(j);
} // end function overlapElement

double GaussianStypeIntegrals::kineticElement(const unsigned int& i, const
        unsigned int& j) {
    /* evaluate and return preset kineticFunc */
    double sum = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        double a = StypeBasis::getBasis()->getPrimitive(i)->scaling(d);
        double b = StypeBasis::getBasis()->getPrimitive(j)->scaling(d);
        double prod = a*b * (sqrtpi / (2*pow((a+b), 1.5)));
        for (unsigned int dd = 0; dd < m_dim; ++dd) {
            if (dd != d) {
                prod *= sqrtpi / sqrt(StypeBasis::getBasis()->getPrimitive(i)
                        -> scaling(dd) +
                        StypeBasis::getBasis()->getPrimitive(j) ->
                        scaling(dd));
            } // end if
        } // end fordd
        sum += prod;
    } // end ford
    return 0.5 * sum * normalizationFactors(i) * normalizationFactors(j);
} // end function kineticElement

double GaussianStypeIntegrals::potentialElement(const unsigned int& i, const
        unsigned int& j) {
    double prod = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        prod *= sqrt2pi /
            pow(StypeBasis::getBasis()->getPrimitive(i)->scaling(d) +
                    StypeBasis::getBasis()->getPrimitive(j)->scaling(d), 1.5);
    } // end ford
    return 0.5 * prod * normalizationFactors(i) * normalizationFactors(j);
} // end function potentialElement

double GaussianStypeIntegrals::coulombElement(const unsigned int& i, const
        unsigned int& j, const unsigned int& k, const unsigned int& l) {
    return 0.0;
}  // end function coulombElement
