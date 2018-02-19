#include "gaussianstypeintegrals.h"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <iostream>

GaussianStypeIntegrals::GaussianStypeIntegrals(const unsigned int dim, unsigned
        int numParameters, double scaling) : StypeBasis(dim) {
    m_dim = dim;
    sqrt2pi = sqrt(2*M_PI);
    sqrtpi = sqrt(M_PI);
    sqrtpiDim = pow(sqrtpi, dim);
    if (m_dim==2) {
        coulombElementFunc = &GaussianStypeIntegrals::coulombElement2D;
    } else {
        coulombElementFunc = &GaussianStypeIntegrals::coulombElement3D;
    } // end if
}// end constructor

GaussianStypeIntegrals::~GaussianStypeIntegrals() {
} // end deconstructor

StypeBasis* GaussianStypeIntegrals::getBasis() {
    /* return a pointer to basis class */
    return dynamic_cast<StypeBasis*>(this);
} // end function getBasis

void GaussianStypeIntegrals::initializeParameters(const Eigen::VectorXd&
        scalingVector, const Eigen::MatrixXd& centralMatrix) {
    /* initialize scaling factors and primitives for isotropic case */
    StypeBasis::setPrimitives(scalingVector, centralMatrix);
    setNormalizations();
    isotropic = true;
} // end function initializeParameters

void GaussianStypeIntegrals::initializeParameters(const Eigen::MatrixXd&
        scalingMatrix, const Eigen::MatrixXd& centralMatrix) {
    /* initialize scaling factors and primitives for non-isotropic case */
    StypeBasis::setPrimitives(scalingMatrix, centralMatrix);
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
    double expSum = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        double sid = StypeBasis::getBasis()->getPrimitive(i)->scaling(d);
        double sjd = StypeBasis::getBasis()->getPrimitive(j)->scaling(d);
        double sijd = sid + sjd;
        prod *= sijd;
        expSum += sid*sjd/sijd *
            pow((StypeBasis::getBasis()->getPrimitive(i)->central(d) -
                        StypeBasis::getBasis()->getPrimitive(j)->central(d)),
                    2);
    } // end ford
    return sqrtpiDim / sqrt(prod) * exp(-expSum) * normalizationFactors(i) *
        normalizationFactors(j);
} // end function overlapElement

inline double GaussianStypeIntegrals::laplacianElement(const unsigned int& i,
        const unsigned int& j) {
    double sum = 0.0;
    double prod2 = 1.0;
    double expSum = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        double sid = StypeBasis::getBasis()->getPrimitive(i)->scaling(d);
        double sjd = StypeBasis::getBasis()->getPrimitive(j)->scaling(d);
        double proddd = 0.5*sid*sjd*sqrtpi/pow(sqrt(sid+sjd), 3/2.);
        double cid = StypeBasis::getBasis()->getPrimitive(i)->central(d);
        double cjd = StypeBasis::getBasis()->getPrimitive(j)->central(d);
        double sdProd = sid*sjd/(sid+sjd);
        expSum +=  sdProd * pow(cid - cjd, 2);
        for (unsigned int dd = 0; dd < m_dim; ++dd) {
            if (dd != d) {
                proddd *= sqrtpi /
                    sqrt(StypeBasis::getBasis()->getPrimitive(i)->scaling(dd) +
                            StypeBasis::getBasis()->getPrimitive(j)->scaling(dd));
            } // end if
        } // end fordd
        prod2 *= sqrtpi * sdProd / (sid+sjd) * (cid*cid - cjd*cjd);
        sum += proddd;
    } // end ford
    return 4 * exp(-expSum) * sum * prod2 * normalizationFactors(i) *
        normalizationFactors(j);
} // end function kineticElement

inline double GaussianStypeIntegrals::potentialElement(const unsigned int& i,
        const unsigned int& j) {
    /* calculate HO potential term */
    double prod = 1.0;
    double expSum = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        double sid = StypeBasis::getBasis()->getPrimitive(i)->scaling(d);
        double sjd = StypeBasis::getBasis()->getPrimitive(j)->scaling(d);
        double sijd = sid + sjd;
        double cid = StypeBasis::getBasis()->getPrimitive(i)->central(d);
        double cjd = StypeBasis::getBasis()->getPrimitive(j)->central(d);
        expSum += sid*sjd/sijd * pow((cid - cjd), 2);
        double prodd = pow((sid*cid+sjd*cjd)*sijd, 2);
        for (unsigned int dd = 0; dd < m_dim; ++dd) {
            if (dd != d) {
                prodd /=
                    pow((StypeBasis::getBasis()->getPrimitive(i)->scaling(dd) +
                                StypeBasis::getBasis()->getPrimitive(j)->scaling(dd)),
                            2);
            } // end if
        } // end fordd
        prod *= sqrtpi * (prodd / sqrt(sijd) + 0.5 / pow(sqrt(sijd), 3));
    } // end ford
    return 0.5 * prod * exp(-expSum) * normalizationFactors(i) * normalizationFactors(j);
} // end function potentialElement

double GaussianStypeIntegrals::oneBodyElement(const unsigned int& i, const
        unsigned int& j) {
    /* calculate and return oneBodyElement <i|h|k> = <i|K|j> + <i|P|j>, where K
     * is the kinetic part and P is the potential part */
    return -0.5*laplacianElement(i,j) + potentialElement(i,j);
} // end function oneBodyElement

inline double GaussianStypeIntegrals::coulombElement2D(const unsigned int& i,
        const unsigned int& j, const unsigned int& k, const unsigned int& l) {
    /* coulomb 2D isotropic */
    double expSum1 = 0.0;
    double expSum2 = 0.0;
    double mu = StypeBasis::getBasis()->getPrimitive(i)->scaling(0) +
        StypeBasis::getBasis()->getPrimitive(k)->scaling(0);
    double nu = StypeBasis::getBasis()->getPrimitive(j)->scaling(0) +
        StypeBasis::getBasis()->getPrimitive(l)->scaling(0);
    double sigma = (mu+nu) / (4*mu*nu);
    double DeltaSq = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        double sid1 = StypeBasis::getBasis()->getPrimitive(i)->scaling(d);
        double sjd1 = StypeBasis::getBasis()->getPrimitive(k)->scaling(d);
        double sijd1 = sid1 + sjd1;
        double cid1 = StypeBasis::getBasis()->getPrimitive(i)->central(d);
        double cjd1 = StypeBasis::getBasis()->getPrimitive(k)->central(d);
        expSum1 += sid1*sjd1/sijd1 * pow((cid1 - cjd1), 2);
        double sid2 = StypeBasis::getBasis()->getPrimitive(j)->scaling(d);
        double sjd2 = StypeBasis::getBasis()->getPrimitive(l)->scaling(d);
        double sijd2 = sid2 + sjd2;
        double cid2 = StypeBasis::getBasis()->getPrimitive(j)->central(d);
        double cjd2 = StypeBasis::getBasis()->getPrimitive(l)->central(d);
        expSum2 += sid2*sjd2/sijd2 * pow((cid2 - cjd2), 2);
        double prodd1 = pow((sid2*cid2+sjd2*cjd2)*sijd2, 2);
        double prodd2 = pow((sid2*cid2+sjd2*cjd2)*sijd2, 2);
        for (unsigned int dd = 0; dd < m_dim; ++dd) {
            if (dd != d) {
                prodd1 /=
                    pow((StypeBasis::getBasis()->getPrimitive(i)->scaling(dd) +
                                StypeBasis::getBasis()->getPrimitive(k)->scaling(dd)),
                            2);
                prodd2 /=
                    pow((StypeBasis::getBasis()->getPrimitive(j)->scaling(dd) +
                                StypeBasis::getBasis()->getPrimitive(l)->scaling(dd)),
                            2);
            } // end if
        } //  end fordd
        DeltaSq += prodd1 - prodd2;
    } // end ford
    double DeltaFactor = -DeltaSq/(8*sigma);
    return normalizationFactors(i) * normalizationFactors(j) *
        normalizationFactors(k) * normalizationFactors(l) * exp(-expSum1) *
        exp(-expSum2) * pow(M_PI, 5/2.) / (2*sqrt(sigma)*mu*nu) *
        exp(DeltaFactor)*boost::math::cyl_bessel_i(0, DeltaFactor);
}  // end function coulombElement

inline double GaussianStypeIntegrals::coulombElement3D(const unsigned int& i,
        const unsigned int& j, const unsigned int& k, const unsigned int& l) {
    /* coulomb 3D isotropic */
    double expSum1 = 0.0;
    double expSum2 = 0.0;
    double mu = StypeBasis::getBasis()->getPrimitive(i)->scaling(0) +
        StypeBasis::getBasis()->getPrimitive(k)->scaling(0);
    double nu = StypeBasis::getBasis()->getPrimitive(j)->scaling(0) +
        StypeBasis::getBasis()->getPrimitive(l)->scaling(0);
    double sigma = (mu+nu) / (4*mu*nu);
    double Delta = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        double sid1 = StypeBasis::getBasis()->getPrimitive(i)->scaling(d);
        double sjd1 = StypeBasis::getBasis()->getPrimitive(k)->scaling(d);
        double sijd1 = sid1 + sjd1;
        double cid1 = StypeBasis::getBasis()->getPrimitive(i)->central(d);
        double cjd1 = StypeBasis::getBasis()->getPrimitive(k)->central(d);
        expSum1 += sid1*sjd1/sijd1 * pow((cid1 - cjd1), 2);
        double sid2 = StypeBasis::getBasis()->getPrimitive(j)->scaling(d);
        double sjd2 = StypeBasis::getBasis()->getPrimitive(l)->scaling(d);
        double sijd2 = sid2 + sjd2;
        double cid2 = StypeBasis::getBasis()->getPrimitive(j)->central(d);
        double cjd2 = StypeBasis::getBasis()->getPrimitive(l)->central(d);
        expSum2 += sid2*sjd2/sijd2 * pow((cid2 - cjd2), 2);
        double prodd1 = pow((sid2*cid2+sjd2*cjd2)*sijd2, 2);
        double prodd2 = pow((sid2*cid2+sjd2*cjd2)*sijd2, 2);
        for (unsigned int dd = 0; dd < m_dim; ++dd) {
            if (dd != d) {
                prodd1 /=
                    pow((StypeBasis::getBasis()->getPrimitive(i)->scaling(dd) +
                                StypeBasis::getBasis()->getPrimitive(k)->scaling(dd)),
                            2);
                prodd2 /=
                    pow((StypeBasis::getBasis()->getPrimitive(j)->scaling(dd) +
                                StypeBasis::getBasis()->getPrimitive(l)->scaling(dd)),
                            2);
            } // end if
        } //  end fordd
        Delta += sqrt(prodd1 - prodd2);
    } // end ford
    return normalizationFactors(i) * normalizationFactors(j) *
        normalizationFactors(k) * normalizationFactors(l) * exp(-expSum1) *
        exp(-expSum2) * pow(M_PI, 3) / (4*sigma*pow(mu*nu,3/2.)) *
        boost::math::erf(Delta*sqrt(4*sigma))/Delta;
} // end function coulombElement3D

double GaussianStypeIntegrals::coulombElement(const unsigned int& i, const
        unsigned int& j, const unsigned int& k, const unsigned int& l) {
    return (this->*coulombElementFunc)(i,j,k,l);
} // end function coulombElement
