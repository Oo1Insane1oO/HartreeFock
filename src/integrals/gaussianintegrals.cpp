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
    sqrtScale1 = sqrt(pow(xScale, m_dim));
    sqrtScale = 1./sqrtScale1;
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
    /* calculate and return <g_n|g_m> (overlap) in 1 dimension */
    int s = n + m;
    if (s%2 || s<=-1) {
        return 0.0;
    } // end if
    return boost::math::tgamma<double>((s+1)/2.) / sqrt(xScale);
} // end function overlapd

inline double GaussianIntegrals::ddexpr(const int& ndd, const int& mdd,
        double(GaussianIntegrals::* f)(const int&, const int&)) {
    /* expression for sum over contracted functions */
    double sums = 0.0;
    int p = 0;
    int q = 0;
    do {
        /* run atleast once to take care of ndd=0 */
        do {
            /* run atleast once to take care of mdd=0 */
            sums += HC(ndd)[p]*HC(mdd)[q] * (this->*f)(p,q);
            q++;
        } while(q<mdd);
        p++;
    } while(p<ndd);
    return sums;
} // end function ddexpr

inline double GaussianIntegrals::ddexprOverlap(const int& p, const int& q) {
    /* expression for 1D overlap element */
    return xScale * overlapd(p,q);
} // end function ddexpr1

inline double GaussianIntegrals::ddexprLaplacian(const int& p, const int& q) {
    /* expression for 1D laplacian element */
    return xScale * (q*(q-1)*overlapd(p,q-2) - (2*q+1)*overlapd(p,q) +
            overlapd(p,q+2));
} // end function ddexpr2

inline double GaussianIntegrals::ddexprPotential(const int& p, const int& q) {
    /* expression for 1D potential element */
    return xScale * overlapd(p,q+2);
} // end function ddexpr1

double GaussianIntegrals::overlapElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return the overlap integral element <i|j> */
    double prod = 1.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        prod *= ddexpr(*(GaussianBasis::Cartesian::getStates(i)(d)),
                *(GaussianBasis::Cartesian::getStates(j)(d)),
                &GaussianIntegrals::ddexprOverlap);
    } // end ford

//     return (i==j ? sqrtScale1 : 0.0);
    return prod * sqrtScale1 * normalizationFactor(i) * normalizationFactor(j);
} // end function overlapElement

double GaussianIntegrals::kineticElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return the kinetic integral element -0.5<i|nabla|j> */
    double sums = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        double tmpProdsd = 1.0;
        for (unsigned int dd = 0; dd < m_dim; ++dd) {
            int ndd = *(GaussianBasis::Cartesian::getStates(i)(dd));
            int mdd = *(GaussianBasis::Cartesian::getStates(j)(dd));
            if (dd != d) {
                tmpProdsd *= ddexpr(ndd, mdd,
                        &GaussianIntegrals::ddexprOverlap);
            } else {
                tmpProdsd *= ddexpr(ndd, mdd,
                        &GaussianIntegrals::ddexprLaplacian);
            } // end ifelse
        } // end fordd
        sums += tmpProdsd;
    } // end ford

    return -0.5*sums * normalizationFactor(i) * normalizationFactor(j);
} // end function kineticElement

double GaussianIntegrals::potentialElement(const unsigned int& i, const
        unsigned int& j) {
    /* calculate and return the kinetic integral element -<i|nabla|j> */
    double sums = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        double tmpProdsd = 1.0;
        for (unsigned int dd = 0; dd < m_dim; ++dd) {
            int ndd = *(GaussianBasis::Cartesian::getStates(i)(dd));
            int mdd = *(GaussianBasis::Cartesian::getStates(j)(dd));
            if (dd != d) {
                tmpProdsd *= ddexpr(ndd, mdd,
                        &GaussianIntegrals::ddexprOverlap);
            } else {
                tmpProdsd *= ddexpr(ndd, mdd,
                        &GaussianIntegrals::ddexprPotential);
            } // end ifelse
        } // end fordd
        sums += tmpProdsd;
    } // end ford

    return 0.5*xScale*sums * normalizationFactor(i) * normalizationFactor(j);
} // end function kinetic

double GaussianIntegrals::coulombElement(const unsigned int& i, const unsigned
        int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> */
    return 0.0;
} // end function coulombElement

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
