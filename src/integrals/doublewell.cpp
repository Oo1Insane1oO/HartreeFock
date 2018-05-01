#include "doublewell.h"
#include "../hermite/hermite.h"
#include "../methods.h"

#include <boost/math/special_functions/gamma.hpp>

DoubleWell::DoubleWell(const unsigned int dim, unsigned int cutOff) :
    GaussianIntegrals(dim, cutOff), DWC() {
} // end constructor

DoubleWell::~DoubleWell() {
} // end deconstructor

std::string DoubleWell::initializeParameters(double _R, unsigned int axis) {
    /* initialize GaussianIntegrals and grab well-separation parameter R */

    // default to shifting well in x-direction
    m_axis = axis;
    R = _R;
    RsqrdFactor = 1.0/8.0 * R*R;
    std::string message = GaussianIntegrals::initializeParameters(1.0);

    m_numBasis = GaussianIntegrals::GaussianBasis::getSize()/2;

    setNormalizations();

    return message;
} // end function initializeParameters

unsigned int DoubleWell::getSize() {
    return 2*m_numBasis;
} // end function getSize

void DoubleWell::setNormalizations() {
    /* precalculate normalizations in arraty */
    normalizationFactors = Eigen::ArrayXd::Constant(m_numBasis, 1.0);
    for (unsigned int k = 0; k  < m_numBasis; ++k) {
        normalizationFactors(k) = 1.0 / sqrt(overlapElementNonNormal(k,k));
    } // end fork
} // end function setNormalizations

const double& DoubleWell::normalizationFactor(const unsigned int& n) const {
    /* normalization for Gauss-Hermite of order n */
    return normalizationFactors(n);
} // end function normalizationFactor

double DoubleWell::overlapElementNonNormal(const unsigned int& i, const
        unsigned int& j) {
    /* calculate and return the overlap integral element <i|j> (not normalized,
     * for use in setNormalizations function) */
    double res = 0.0;
    for (unsigned int p = 0; p < m_numBasis; ++p) {
        for (unsigned int q = 0; q < m_numBasis; ++q) {
            res += DWC::C(i,p) * DWC::C(j,q) *
                GaussianIntegrals::overlapElement(p,q);
        } // end forq
    } // end forp

    return res;
} // end function overlapElement

double DoubleWell::overlapElement(const unsigned int& i, const unsigned int& j)
{
    /* calculate and return the overlap integral element <i|j> */
    double res = 0.0;
    for (unsigned int p = 0; p < m_numBasis; ++p) {
        for (unsigned int q = 0; q < m_numBasis; ++q) {
            res += DWC::C(i,p) * DWC::C(j,q) *
                GaussianIntegrals::overlapElement(p,q);
        } // end forq
    } // end forp

    return res * normalizationFactor(i) * normalizationFactor(j);
} // end function overlapElement

double DoubleWell::oneBodyElement(const unsigned int& i, const unsigned int& j)
{
    /* calculate and return oneBodyElement <i|h|k> = <i|K|j> + <i|P|j>, where K
     * is the kinetic part and P is the potential part. For this case the K+P
     * part is the HO-potential part taken from GaussianIntegrals with the
     * added DW part */
    double res = 0.0;
    for (unsigned int p = 0; p < m_numBasis; ++p) {
        res += RsqrdFactor +
            *(GaussianIntegrals::GaussianBasis::Cartesian::getStates(p)(m_dim +
                        2));
    } // end forp

    return res + potentialDWElement(i,j);
} // end function oneBodyElement

double DoubleWell::potentialDWElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate -0.5*w*R*abs(axis) */
    double sum = 0.0;
    for (unsigned int p = 0; p < m_numBasis; ++p) {
        for (unsigned int q = 0; q < m_numBasis; ++q) {
            double res = 1.0;
            for (unsigned int d = 0; d < m_dim; ++d) {
                const int& nd =
                    GaussianIntegrals::GaussianBasis::Cartesian::getn(p,d);
                const int& md =
                    GaussianIntegrals::GaussianBasis::Cartesian::getn(q,d);
                if (d != m_axis) {
                    res *= GaussianIntegrals::ddexpr(nd, md,
                            &DoubleWell::ddexprOverlap);
                } else {
                    res *= potDWSum(nd, md);
                } // end ifselse
            } // end ford
            sum += DWC::C(i,p) * DWC::C(j,q) * res;
        } // end forq 
    } // end forp

    return -0.5 * R * sum * normalizationFactor(i) * normalizationFactor(j);
} // end function potentialDWElement

double DoubleWell::potDWSum(const int& ndd, const int& mdd) {
    /* expression for sum over contracted functions */
    auto potDW = [](const int& p, const int& q) {
        /* analytic expression for integrals over abs(axis) times a gaussian */
        int s = p + q;
        if ((s<=-1) or (s%2==1)) {
            return static_cast<double>(0.0);
        } // end if
        
        return boost::math::tgamma<double>((s+2)/2.);
    }; // end lambda potDw

    double sums = 0.0;
    for (int p = 0; p <= ndd; ++p) {
        for (int q = 0; q <= mdd; ++q) {
            sums += HC(ndd)[p]*HC(mdd)[q] * potDW(p,q);
        } // end forq
    } // end forp

    return sums;
} // end function ddexpr


double DoubleWell::coulombElement(const unsigned int& i, const unsigned int& j,
        const unsigned int& k, const unsigned int& l) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> */
    double res = 0.0;
    for (unsigned int p = 0; p < m_numBasis; ++p) {
        for (unsigned int q = 0; q < m_numBasis; ++q) {
            for (unsigned int r = 0; r < m_numBasis; ++r) {
                for (unsigned int s = 0; s < m_numBasis; ++s) {
                    res += DWC::C(i,p)*DWC::C(j,q)*DWC::C(k,r)*DWC::C(l,s) *
                        GaussianIntegrals::coulombElement(p,q,r,s);
                } // end fors
            } // end forr
        } // end forq
    } // end forp

    return res * normalizationFactor(i) * normalizationFactor(j) *
        normalizationFactor(k) * normalizationFactor(l) ;
} // end function coulombElement 
