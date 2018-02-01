#include "gaussianintegrals.h"
#include <boost/math/special_functions/factorials.hpp>

GaussianIntegrals::GaussianIntegrals(const unsigned int dim, unsigned int
        cutOff, double scaling) : GaussianBasis() {
    m_dim = dim;
    expScaleFactor = scaling;
    sqrtFactor = sqrt(M_PI/expScaleFactor);
    GaussianBasis::setup(cutOff, dim, expScaleFactor);

    xScale = 1.0; // omega in HO case TODO: generalize this 

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
        normalizationFactors(i) =
            1./sqrt(pow(2,i)*boost::math::factorial<double>(i)*sqrtFactor);
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

double GaussianIntegrals::overlapElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return the overlap integral element <i|j> */
    double sum = 0;
    for (unsigned int p = 0; p < i; ++p) {
        for (unsigned int q = 0; q < j; ++q) {
            double prod = 1;
            for (unsigned int d = 0; d < m_dim; ++d) {
                int pq = *(GaussianBasis::Cartesian::getStates(p)(d)) +
                    *(GaussianBasis::Cartesian::getStates(q)(d));
                if (pq%2==0) {
                    /* integral is zero for odd */
                    int s = pq + 1;
                    sum += 1./sqrt(xScale) * 2./s *
                        boost::math::factorial<double>(s/2.);
                } // end if
            } // end ford
        } // end forq
    } // end forp

    return sum;
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
