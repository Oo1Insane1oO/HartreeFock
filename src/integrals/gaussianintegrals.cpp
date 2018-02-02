#include "gaussianintegrals.h"
#include "../hermite/hermite.h"
#include <boost/math/special_functions/factorials.hpp>

GaussianIntegrals::GaussianIntegrals(const unsigned int dim, unsigned int
        cutOff, double scaling) : GaussianBasis() {
    m_dim = dim;
    expScaleFactor = scaling;
    sqrtFactor = sqrt(scaling);
    GaussianBasis::setup(cutOff, dim, expScaleFactor);

    xScale = 1.0; // omega in HO case TODO: generalize this 
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
        Eigen::ArrayXd::Zero(GaussianBasis::Cartesian::getn().size());
    for (unsigned int i = 0; i < GaussianBasis::Cartesian::getn().size(); ++i)
    {
        int n = GaussianBasis::Cartesian::getn(i);
        normalizationFactors(i) =
            1./sqrt(pow(2,n)*boost::math::factorial<double>(n)*sqrtFactor);
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
    double sum = 0;
    for (unsigned int p = 0; p < n; ++p) {
        for (unsigned int q = 0; q < m; ++q) {
            int pq = p + q;
            if (pq%2==0) {
                /* integral is zero for odd powers */
                int s = pq + 1;
                sum += HC(n)[p]*HC(m)[q]/s *
                    boost::math::factorial<double>(s/2);
            } // end if
        } // end forq
    } // end forp

    return sum;
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

    return prod;
} // end function overlap

double GaussianIntegrals::kineticElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return the kinetic integral element -<i|nabla|j> */
    return 0;
} // end function kinetic

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
