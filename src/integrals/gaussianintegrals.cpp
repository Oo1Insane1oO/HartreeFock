#include "gaussianintegrals.h"
#include "../hermite/hermite.h"
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>

GaussianIntegrals::GaussianIntegrals(const unsigned int dim, unsigned int
        cutOff, double scaling) : GaussianBasis(cutOff, dim, scaling) {
    m_dim = dim;
    expScaleFactor = scaling;
    sqrtFactor = sqrt(scaling);

    xScale = 1.0; // omega in HO case FIXME: generalize this 
    sqrtScale = 2./sqrt(pow(xScale, m_dim));
    powScale = pow(xScale, 2*m_dim);

    setF0();
    setF1();

    setNormalizations();
} // end constructor

GaussianIntegrals::~GaussianIntegrals() {
} // end deconstructor

GaussianBasis* GaussianIntegrals::getBasis() {
    /* return a pointer to GaussianBasis */
    // TODO: check that caller actually gets GaussianBasis
    return dynamic_cast<GaussianBasis*>(this);
} // end function getBasis

void GaussianIntegrals::setNormalizations() {
    /* calculate and set normalization factors for all basis functions */
    normalizationFactors =
        Eigen::ArrayXd::Constant(GaussianBasis::Cartesian::getNumberOfStates(),
                1.0);
    for (unsigned int i = 0; i < GaussianBasis::Cartesian::getNumberOfStates();
            ++i) {
        for (unsigned int d = 0; d < m_dim; ++d) {
            int n = *(GaussianBasis::Cartesian::getStates(i)(d));
            normalizationFactors(i) *= 1.0 / sqrt(pow(2,n) *
                    boost::math::factorial<double>(n)*sqrtFactor);
        } // end ford
    } // end fori
} // end function setNormalizations

void GaussianIntegrals::setF0() {
    /* set first incomplete overlap integral factor */
    F0 = sqrt(M_PI/expScaleFactor);
} // end function setF0

void GaussianIntegrals::setF1() {
    /* set second incomplete overlap integral factor */
    F1 = 0;
} // end function setF1

const double& GaussianIntegrals::normalizationFactor(const unsigned int& n)
    const {
    /* normalization for Gauss-Hermite of order n */
    return normalizationFactors(n); 
} // end function normalizationFactor

inline double GaussianIntegrals::overlapd(const unsigned int& n, const unsigned
        int& m) {
    /* calculate and return <g_n|g_m> */
    int s = n + m;
    if (s%2 || s<=-1) {
        return 0.0;
    } // end if
    return boost::math::tgamma<double>((s+1)/2.) / sqrt(xScale);
//     return sqrt(M_PI/xScale) *
//         boost::math::factorial<double>(n)/boost::math::factorial<double>(n/2) *
//         pow(0.5, s);
} // end function overlapd

double GaussianIntegrals::overlapElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return the overlap integral element <i|j> */
    double prod = sqrtScale;
    for (unsigned int d = 0; d < m_dim; ++d) {
        int nid = *(GaussianBasis::Cartesian::getStates(i)(d));
        int njd = *(GaussianBasis::Cartesian::getStates(j)(d));
        prod *= normalizationFactor(nid) * normalizationFactor(njd) *
            overlapd(nid,njd);
    } // end ford

    return prod * normalizationFactor(i) * normalizationFactor(j);
} // end function overlap

double GaussianIntegrals::kineticElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return the kinetic integral element -<i|nabla|j> */
    double sums = 0.0;
    for (unsigned int p = 0; p < i; ++p) {
        for (unsigned int q = 0; q < j; ++q) {
            Eigen::Array3d sumsd = Eigen::Array3d::Zero(3);
            for (unsigned int d = 0; d < m_dim; ++d) {
                Eigen::Array3d tmpProdsdd = Eigen::Array3d::Constant(3,1.0);
                for (unsigned int dd = 0; dd < m_dim; ++dd) {
                    int ndd = *(GaussianBasis::Cartesian::getStates(p)(dd));
                    int mdd = *(GaussianBasis::Cartesian::getStates(q)(dd));
                    double CpCq = HC(p)[ndd]*HC(q)[mdd];
                    if (dd != d) {
                        tmpProdsdd *= CpCq*xScale*overlapd(ndd, mdd);
                    } else {
                        tmpProdsdd(0) *= CpCq*xScale *
                            mdd*(mdd-1)*overlapd(ndd,mdd-2);
                        tmpProdsdd(1) *= CpCq*xScale *
                            (2*mdd+1)*overlapd(ndd,mdd);
                        tmpProdsdd(2) *= CpCq*xScale *
                            overlapd(ndd,mdd+2);
                    } // end if
                } // end fordd
                sumsd += tmpProdsdd;
            } // end ford
            sums += sumsd[0] - sumsd[1] + sumsd[2];
        } // end forq
    } // end forp
    return -0.5*sums * normalizationFactor(i) * normalizationFactor(j); 
} // end function kinetic

double GaussianIntegrals::potentialElement(const unsigned int& i, const
        unsigned int& j) {
    /* calculate and return <i|0.5w^2r^2|j> */
    double sums = 0.0;
    for (unsigned int p = 0; p < i; ++p) {
        for (unsigned int q = 0; q < j; ++q) {
            double sumsd = 0.0;
            for (unsigned int d = 0; d < m_dim; ++d) {
                double tmpProdsdd = 1.0; 
                for (unsigned int dd = 0; dd < m_dim; ++dd) {
                    int ndd = *(GaussianBasis::Cartesian::getStates(p)(dd));
                    int mdd = *(GaussianBasis::Cartesian::getStates(q)(dd));
                    double CpCq = HC(p)[ndd]*HC(q)[mdd];
                    if (dd != d) {
                        tmpProdsdd *= CpCq*xScale*overlapd(ndd,mdd);
                    } else {
                        tmpProdsdd *= CpCq*overlapd(ndd,mdd+2);
                    } // end if
                } // end fordd
                sumsd += tmpProdsdd;
            } // end ford
            sums += sumsd;
        } // end forq
    } // end forp
    return 0.5*xScale*sums * normalizationFactor(i) * normalizationFactor(j);
} // end function potentialElement

double GaussianIntegrals::coulombElement(const unsigned int& i, const unsigned
        int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> */
    return 0;
} // end function coulomb

inline double GaussianIntegrals::incompleteOverlapIntegral(const unsigned int&
        l) {
    /* calculate incomplete integral used in overlap integral element with
     * precalculated F0 and F1 */
    return 0;
} // end function incompleteOverlapIntegral

inline double GaussianIntegrals::incompleteByPartsFactorG(const unsigned int&
        m) {
    /* calculate first part of incomplete overlap integral */
    return 0;
} // end function incompleteByPartsFactor

inline double GaussianIntegrals::incompleteByPartsFactorF(const unsigned int&
        l) {
    /* calculate second part of incomplete overlap integral */
    if (l%2==0) {
        /* even case */
        return 0;
    } // end if

    return 2*sqrt(M_PI) * 2*boost::math::factorial<double>(l) /
        boost::math::factorial<double>(l/2) * pow(1/(2*sqrt(expScaleFactor)),
                l+1);
} // end function incompleteByPartsFactorF
