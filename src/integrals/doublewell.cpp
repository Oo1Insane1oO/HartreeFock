#include "doublewell.h"
#include "../hermite/hermite.h"
#include "../methods.h"
#include "../hartreefocksolver.h"

#include <boost/math/special_functions/gamma.hpp>

DoubleWell::DoubleWell(const unsigned int dim, unsigned int cutOff) :
    GaussianIntegrals(hfsIn, dim, cutOff), DWC() {
    m_numBasis = cutOff;
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

    // increase size of basis
    GaussianIntegrals::GaussianBasis::setup(2*DWC::C.rows(), m_dim);

    return message;
} // end function initializeParameters

unsigned int DoubleWell::getSize() {
    /* return number of states (by two because of spin) */
    return 2*m_numBasis;
} // end function getSize

double DoubleWell::overlapElement(const unsigned int& i, const unsigned int& j)
{
    /* calculate and return the overlap integral element <i|j> */
    double res = 0.0;
    for (unsigned int p = 0; p < DWC::C.rows(); ++p) {
        for (unsigned int q = 0; q < DWC::C.rows(); ++q) {
            res += DWC::C(p,i) * DWC::C(q,j) *
                GaussianIntegrals::overlapElement(p,q);
        } // end forq
    } // end forp

    return res;
} // end function overlapElement

double DoubleWell::oneBodyElement(const unsigned int& i, const unsigned int& j)
{
    /* calculate and return oneBodyElement <i|h|k> = <i|K|j> + <i|P|j>, where K
     * is the kinetic part and P is the potential part. For this case the K+P
     * part is the HO-potential part taken from GaussianIntegrals with the
     * added DW part */
    double res = 0.0;
    if (i == j) {
        /* Add diagonal part */
        res += DWC::epsDW(i);
    } // end if

    return res;
} // end function oneBodyElement

double DoubleWell::potentialDWElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate -0.5*w*R*abs(axis) */
    double sum = 0.0;
    for (unsigned int p = 0; p < DWC::C.rows(); ++p) {
        for (unsigned int q = 0; q < DWC::C.rows(); ++q) {
            double res = 1.0;
            for (unsigned int d = 0; d < m_dim; ++d) {
                const int& nd =
                    GaussianIntegrals::GaussianBasis::Cartesian::getn(p,d);
                const int& md =
                    GaussianIntegrals::GaussianBasis::Cartesian::getn(q,d);
                if (d == m_axis) {
                    res *= potDWSum(nd, md);
                } else {
                    res *= GaussianIntegrals::ddexpr(nd, md,
                            &DoubleWell::ddexprOverlap);
                } // end ifselse
            } // end ford
            sum += DWC::C(p,i) * DWC::C(q,j) * res;
        } // end forq 
    } // end forp

    return -0.5 * R * sum;
} // end function potentialDWElement

double DoubleWell::potDWSum(const int& ndd, const int& mdd) {
    /* expression for sum over contracted functions */
    double sums = 0.0;
    for (int p = 0; p <= ndd; ++p) {
        for (int q = 0; q <= mdd; ++q) {
            int s = p + q;
            if (s%2!=1) {
                sums += HC(ndd)[p]*HC(mdd)[q] *
                    boost::math::tgamma<double>((s+2.)/2.);
            } // end if
        } // end forq
    } // end forp

    return sums;
} // end function potDWSum

double DoubleWell::coulombElement(const unsigned int& i, const unsigned int& j,
        const unsigned int& k, const unsigned int& l) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> */
    double res = 0.0;
    for (unsigned int p = 0; p < DWC::C.rows(); ++p) {
        for (unsigned int q = 0; q < DWC::C.rows(); ++q) {
            for (unsigned int r = 0; r < DWC::C.rows(); ++r) {
                for (unsigned int s = 0; s < DWC::C.rows(); ++s) {
                    res += DWC::C(p,i)*DWC::C(q,j)*DWC::C(r,k)*DWC::C(s,l) *
                        GaussianIntegrals::coulombElement(p,q,r,s);
                } // end fors
            } // end forr
        } // end forq
    } // end forp

    return res;
} // end function coulombElement
