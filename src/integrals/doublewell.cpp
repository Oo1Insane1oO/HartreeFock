#include "doublewell.h"
#include "../hermite/hermite.h"

#include <boost/math/special_functions/gamma.hpp>

DoubleWell::DoubleWell(const unsigned int dim, unsigned int cutOff) :
    GaussianIntegrals(dim, cutOff), DWC() {
} // end constructor

DoubleWell::~DoubleWell() {
} // end deconstructor

std::string DoubleWell::initializeParameters(double _R, unsigned int axis) {
    /* initialize GaussianIntegrals  and grab well-separation parameter R */

    // default to shifting well in x-direction
    m_axis = axis;
    R = _R;
    RsqrdFactor = 1.0/8.0 * R*R;
    std::string message = GaussianIntegrals::initializeParameters(1.0);

    return message;
} // end function initializeParameters

double DoubleWell::overlapElement(const unsigned int& i, const unsigned int& j)
{
    /* calculate and return the overlap integral element <i|j> */
    double res = 0.0;
    for (unsigned int p = 0; p < DWC::C.col(i).size(); ++p) {
        for (unsigned int q = 0; q < DWC::C.col(j).size(); ++q) {
            res += DWC::C(p,i) * C(q,j) *
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
    for (unsigned int p = 0; p < DWC::C.col(i).size(); ++p) {
        for (unsigned int q = 0; q < DWC::C.col(j).size(); ++q) {
            res += DWC::C(p,i) * C(q,j) *
                GaussianIntegrals::oneBodyElement(p,q);
        } // end forq
    } // end forp

    // add overlap part
    if (i == j) {
        res += RsqrdFactor;
    } // end if

    return res + potentialDWElement(i,j);
} // end function oneBodyElement

double DoubleWell::potentialDWElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate -0.5*w*R*abs(axis) */
    double sum = 0.0;
    for (unsigned int p = 0; p < DWC::C.col(i).size(); ++p) {
        for (unsigned int q = 0; q < DWC::C.col(j).size(); ++q) {
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
            sum += res * normalizationFactor(p) * normalizationFactor(q);
        } // end forq 
    } // end forp

    return -0.5 * R * sum;
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
    for (unsigned int p = 0; p < DWC::C.col(i).size(); ++p) {
        for (unsigned int q = 0; q < DWC::C.col(j).size(); ++q) {
            for (unsigned int r = 0; r < DWC::C.col(k).size(); ++r) {
                for (unsigned int s = 0; s < DWC::C.col(l).size(); ++s) {
                    res += DWC::C(p,i)*DWC::C(q,j)*DWC::C(r,k)*DWC::C(s,l) *
                        GaussianIntegrals::coulombElement(p,q,r,s);
                } // end fors
            } // end forr
        } // end forq
    } // end forp

    return res;
} // end function coulombElement 
