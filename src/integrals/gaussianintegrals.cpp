#include "gaussianintegrals.h"
#include "../hermite/hermite.h"
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>

GaussianIntegrals::GaussianIntegrals(const unsigned int dim, unsigned int
        cutOff, double scaling) : GaussianBasis(cutOff, dim) {
    m_dim = dim;
    expScaleFactor = scaling;
    sqrtFactor = sqrt(scaling);
    coeffs = std::make_unique<Hexpander>();
} // end constructor

GaussianIntegrals::~GaussianIntegrals() {
} // end deconstructor

GaussianBasis* GaussianIntegrals::getBasis() {
    /* return a pointer to GaussianBasis */
    return dynamic_cast<GaussianBasis*>(this);
} // end function getBasis

void GaussianIntegrals::initializeParameters(double omega) {
    /* set value of oscillator frequency */
    xScale = omega;
    xScaleHalf = 0.5*xScale;
    sqrtScale1 = sqrt(pow(xScale, m_dim));
    sqrtScale = 1./sqrtScale1;
    powScale = pow(xScale, 2*m_dim);

    // choose coulombElement function for 2D or 3D
    if (m_dim == 2) {
        coulombFunc = &GaussianIntegrals::coulomb2D;
    } else {
        coulombFunc = &GaussianIntegrals::coulomb3D;
    } // end if

    // calculate and set normalization factors to array
    setNormalizations();
} // end function setPositionScaling

void GaussianIntegrals::setNormalizations() {
    /* calculate and set normalization factors for all basis functions */
    normalizationFactors =
        Eigen::ArrayXd::Constant(GaussianBasis::Cartesian::getNumberOfStates(),
                1.0);
    for (unsigned int i = 0; i < GaussianBasis::Cartesian::getNumberOfStates();
            ++i) {
        for (unsigned int d = 0; d < m_dim; ++d) {
            int n = *(GaussianBasis::Cartesian::getStates(i)(d));
            normalizationFactors(i) *= 1.0 / ddexpr(n,n,
                    &GaussianIntegrals::ddexprOverlap);
        } // end ford
    } // end fori
    normalizationFactors = normalizationFactors.cwiseSqrt();
} // end function setNormalizations

const double& GaussianIntegrals::normalizationFactor(const unsigned int& n)
    const {
    /* normalization for Gauss-Hermite of order n */
    return normalizationFactors(n); 
} // end function normalizationFactor

inline double GaussianIntegrals::overlapd(const unsigned int& n, const unsigned
        int& m) {
    /* calculate and return <g_n|g_m> (overlap) in 1 dimension */
    int s = n + m;
    if ((s<=-1) || (s%2)) {
        return 0.0;
    } // end if

    return boost::math::tgamma<double>((s+1)/2.) / sqrt(xScale);
} // end function overlapd

inline double GaussianIntegrals::ddexpr(const int& ndd, const int& mdd,
        double(GaussianIntegrals::* f)(const int&, const int&)) {
    /* expression for sum over contracted functions */
    double sums = 0.0;
    for (int p = 0; p <= ndd; ++p) {
        for (int q = 0; q <= mdd; ++q) {
            sums += HC(ndd)[p]*HC(mdd)[q] * (this->*f)(p,q);
        } // end forq
    } // end forp

    return sums;
} // end function ddexpr

inline double GaussianIntegrals::ddexprOverlap(const int& p, const int& q) {
    /* expression for 1D overlap element */
    return overlapd(p,q);
} // end function ddexpr1

inline double GaussianIntegrals::ddexprLaplacian(const int& p, const int& q) {
    /* expression for 1D laplacian element */
    return xScale * (q*(q-1)*overlapd(p,q-2) - (2*q+1)*overlapd(p,q) +
            overlapd(p,q+2));
} // end function ddexpr2

inline double GaussianIntegrals::ddexprPotential(const int& p, const int& q) {
    /* expression for 1D potential element */
    return overlapd(p,q+2);
} // end function ddexpr1

inline double GaussianIntegrals::laplacianElement(const unsigned int& i, const
        unsigned int& j) {
    /* calculate and return the laplacian integral element <i|nabla|j> */
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

    return sums * normalizationFactor(i) * normalizationFactor(j);
} // end function laplacianElement 

inline double GaussianIntegrals::potentialElement(const unsigned int& i, const
        unsigned int& j) {
    /* calculate and return the HO potential integral element <i|0.5wr^2|j> */
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

    return xScaleHalf*sums * normalizationFactor(i) * normalizationFactor(j);
} // end function potentialElement

inline double GaussianIntegrals::coulombElement2D(const unsigned int& ix, const
        unsigned int& iy, const unsigned int& jx, const unsigned int& jy, const
        unsigned int& kx, const unsigned int& ky, const unsigned int& lx, const
        unsigned int& ly) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> in 2D for level (i,k,j,l) (calculate 1D integral
     * numerically) */
    double sum = 0.0;
//     static Eigen::VectorXd centerVec = Eigen::VectorXd::Constant(m_dim, 0.0);
    for (unsigned int px = 0; px <= ix+kx; ++px) {
        for (unsigned int py = 0; py <= iy+ky; ++py) {
            for (unsigned int qx = 0; qx <= jx+lx; ++qx) {
                for (unsigned int qy = 0; qy <= jy+ly; ++qy) {
                    double tmp = coeffs->coeff(ix,kx,px) *
                        coeffs->coeff(iy,ky,py) * coeffs->coeff(jx,lx,qx) *
                        coeffs->coeff(jy,ly,qy) *
                        coeffs->auxiliary2D(0,px+qx,py+qy);
                    // fix sign in (-1)^(qx + qy) part
                    sum += (((qx+qy)%2==0) ? tmp : -tmp);
//                     Methods::sepPrint(px, py, qx, qy);
                } // end forqy
            } // end forqx
        } // end forpy
    } // end forpx
    return sum;
} // end function coulombElement2D

inline double GaussianIntegrals::coulombElement3D(const unsigned int& ix, const
        unsigned int& iy, const unsigned int& iz, const unsigned int& jx, const
        unsigned int& jy, const unsigned int& jz, const unsigned int& kx, const
        unsigned int& ky, const unsigned int& kz, const unsigned int& lx, const
        unsigned int& ly, const unsigned int& lz) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> in 3D for level (i,k,j,l) (calculate 1D integral
     * numerically) */
    double sum = 0.0;
    static Eigen::VectorXd centerVec = Eigen::VectorXd::Constant(m_dim, 0.0);
    for (unsigned int px = 0; px <= ix+kx; ++px) {
        for (unsigned int py = 0; py <= iy+ky; ++py) {
            for (unsigned int pz = 0; pz <= iz+kz; ++pz) {
                for (unsigned int qx = 0; qx <= jx+lx; ++qx) {
                    for (unsigned int qy = 0; qy <= jy+ly; ++qy) {
                        for (unsigned int qz = 0; qz <= jz+lz; ++qz) {
                            double tmp = coeffs->coeff(ix, kx, px) *
                                coeffs->coeff(iy, ky, py) * coeffs->coeff(iz,
                                        kz, pz) * coeffs->coeff(jx, lx, qx) *
                                coeffs->coeff(jy, ly, qy) * coeffs->coeff(jz,
                                        lz, qz) *
                                GaussianBasis::Hexpander::auxiliary3D(px+qx,
                                        py+qy, pz+qz, 0, xScaleHalf, centerVec,
                                        0.0);
                            // fix sign in (-1)^(qx + qy + qz) part
                            sum += (((qx+qy+qz)%2==0) ? tmp : -tmp);
                        } // end forqz
                    } // end forqy
                } // end forqx
            } // end forpz
        } // end forpy
    } // end forpx
    return sum * 2*pow(M_PI, 2.5) * pow(xScaleHalf, 1.5);
} // end function coulombElement3D

inline double GaussianIntegrals::coulomb2D(const unsigned int& i, const unsigned
        int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate full coulomb integral in 2D case */

    // grab hermite coefficients
    const std::vector<int>& HCIx =
        HC(*(GaussianBasis::Cartesian::getStates(i)(0)));
    const std::vector<int>& HCIy =
        HC(*(GaussianBasis::Cartesian::getStates(i)(1)));
    const std::vector<int>& HCJx =
        HC(*(GaussianBasis::Cartesian::getStates(j)(0)));
    const std::vector<int>& HCJy =
        HC(*(GaussianBasis::Cartesian::getStates(j)(1)));
    const std::vector<int>& HCKx =
        HC(*(GaussianBasis::Cartesian::getStates(k)(0)));
    const std::vector<int>& HCKy =
        HC(*(GaussianBasis::Cartesian::getStates(k)(1)));
    const std::vector<int>& HCLx =
        HC(*(GaussianBasis::Cartesian::getStates(l)(0)));
    const std::vector<int>& HCLy =
        HC(*(GaussianBasis::Cartesian::getStates(l)(1)));

    // find the maximum index (upper limit for coeffs and integrals needed)
    unsigned int pmax = Methods::max(HCIx.size(), HCKx.size(), HCJx.size(),
            HCLx.size(), HCIy.size(), HCKy.size(), HCJy.size(), HCLy.size());
    unsigned int auxMax = 2*Methods::max(HCIx.size()+HCKx.size(),
            HCJx.size()+HCLx.size(), HCIy.size()+HCKy.size(),
            HCJy.size()+HCLy.size());
    static Eigen::VectorXd centerVec = Eigen::VectorXd::Constant(m_dim, 0.0);

    // set all coefficients and integrals needed
    coeffs->setCoefficients(pmax, pmax, xScaleHalf, xScaleHalf, 0.0);
    coeffs->setAuxiliary2D(auxMax, auxMax, xScaleHalf, xScaleHalf, xScaleHalf,
            xScaleHalf, centerVec);

    double sum = 0.0;
    for (unsigned int ix = 0; ix < HCIx.size(); ++ix) {
        if ((HCIx.size()%2==0 && ix%2==0) || (HCIx.size()%2!=0 && ix%2!=0)) {
            /* ignore calculation when Hermite Coefficient is zero */
            continue;
        } // end if
        for (unsigned int iy = 0; iy < HCIy.size(); ++iy) {
            if ((HCIy.size()%2==0 && iy%2==0) || (HCIy.size()%2!=0 && iy%2!=0))
            {
                /* ignore calculation when Hermite Coefficient is zero */
                continue;
            } // end if
            for (unsigned int jx = 0; jx < HCJx.size(); ++jx) {
                if ((HCJx.size()%2==0 && jx%2==0) || (HCJx.size()%2!=0 &&
                            jx%2!=0)) {
                    /* ignore calculation when Hermite Coefficient is zero */
                    continue;
                } // end if
                for (unsigned int jy = 0; jy < HCJy.size(); ++jy) {
                    if ((HCJy.size()%2==0 && jy%2==0) || (HCJy.size()%2!=0 &&
                                jy%2!=0)) {
                        /* ignore calculation when Hermite Coefficient is zero
                         * */
                        continue;
                    } // end if
                    for (unsigned int kx = 0; kx < HCKx.size(); ++kx) {
                        if ((HCKx.size()%2==0 && kx%2==0) || (HCKx.size()%2!=0 &&
                                    kx%2!=0)) {
                            /* ignore calculation when Hermite Coefficient is
                             * zero */
                            continue;
                        } // end if
                        for (unsigned int ky = 0; ky < HCKy.size(); ++ky) {
                            if ((HCKy.size()%2==0 && ky%2==0) ||
                                    (HCKy.size()%2!=0 && ky%2!=0)) {
                                /* ignore calculation when Hermite Coefficient
                                 * is zero */
                                continue;
                            } // end if
                            for (unsigned int lx = 0; lx < HCLx.size(); ++lx) {
                                if ((HCLx.size()%2==0 && lx%2==0) ||
                                        (HCLx.size()%2!=0 && lx%2!=0)) {
                                    /* ignore calculation when Hermite
                                     * Coefficient is zero */
                                    continue;
                                } // end if
                                for (unsigned int ly = 0; ly < HCLy.size();
                                        ++ly) {
                                    if ((HCLy.size()%2==0 && ly%2==0) ||
                                            (HCLy.size()%2!=0 && ly%2!=0)) {
                                        /* ignore calculation when Hermite
                                         * Coefficient is zero */
                                        continue;
                                    } // end if
                                    sum += HCIx[ix] * HCIy[iy] * HCKx[kx] *
                                        HCKy[ky] * HCJx[jx] * HCJy[jy] *
                                        HCLx[lx] * HCLy[ly] *
                                        coulombElement2D(ix,iy, jx,jy, kx,ky,
                                                lx,ly);
                                } // end forly
                            } // end forlx
                        } // end forky
                    } // end forkx
                } // end forjy
            } // end forjx
        } // end foriy
    } // end forix

    return sum * 2 * pow(M_PI, 1.5) / (4 * sqrt(xScaleHalf));
} // end function coulomb2D

inline double GaussianIntegrals::coulomb3D(const unsigned int& i, const unsigned
        int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate full coulomb integral in 2D case */
    const std::vector<int>& HCIx =
        HC(*(GaussianBasis::Cartesian::getStates(i)(0)));
    const std::vector<int>& HCIy =
        HC(*(GaussianBasis::Cartesian::getStates(i)(1)));
    const std::vector<int>& HCIz =
        HC(*(GaussianBasis::Cartesian::getStates(i)(2)));
    const std::vector<int>& HCKx =
        HC(*(GaussianBasis::Cartesian::getStates(k)(0)));
    const std::vector<int>& HCKy =
        HC(*(GaussianBasis::Cartesian::getStates(k)(1)));
    const std::vector<int>& HCKz =
        HC(*(GaussianBasis::Cartesian::getStates(k)(2)));
    const std::vector<int>& HCJx =
        HC(*(GaussianBasis::Cartesian::getStates(j)(0)));
    const std::vector<int>& HCJy =
        HC(*(GaussianBasis::Cartesian::getStates(j)(1)));
    const std::vector<int>& HCJz =
        HC(*(GaussianBasis::Cartesian::getStates(j)(2)));
    const std::vector<int>& HCLx =
        HC(*(GaussianBasis::Cartesian::getStates(l)(0)));
    const std::vector<int>& HCLy =
        HC(*(GaussianBasis::Cartesian::getStates(l)(1)));
    const std::vector<int>& HCLz =
        HC(*(GaussianBasis::Cartesian::getStates(l)(2)));
    
    unsigned int pmax = Methods::max(HCIx.size(), HCKx.size(), HCJx.size(),
            HCLx.size(), HCIy.size(), HCKy.size(), HCJy.size(), HCLy.size(),
            HCIz.size(), HCKz.size(), HCJz.size(), HCLz.size());
    coeffs->setCoefficients(pmax, pmax, xScaleHalf, xScaleHalf, 0.0);

    double sum = 0.0;
    for (unsigned int ix = 0; ix < HCIx.size(); ++ix) {
        for (unsigned int iy = 0; iy < HCIy.size(); ++iy) {
            for (unsigned int iz = 0; iy < HCIz.size(); ++iz) {
                for (unsigned int kx = 0; kx < HCKx.size(); ++kx) {
                    for (unsigned int ky = 0; ky < HCKy.size(); ++ky) {
                        for (unsigned int kz = 0; iz < HCKz.size(); ++kz) {
                            for (unsigned int jx = 0; jx < HCJx.size(); ++jx) {
                                for (unsigned int jy = 0; jy < HCJy.size();
                                        ++jy) {
                                    for (unsigned int jz = 0; jz < HCJz.size();
                                            ++jz) {
                                        for (unsigned int lx = 0; lx <
                                                HCLx.size(); ++lx) {
                                            for (unsigned int ly = 0; ly <
                                                    HCLy.size(); ++ly) {
                                                for (unsigned int lz = 0; lz <
                                                        HCLz.size(); ++lz) {
                                                    sum += 
                                                        HCIx[ix] * HCIy[iy] *
                                                        HCIz[iz] * HCKx[kx] *
                                                        HCKy[ky] * HCKz[kz] *
                                                        HCJx[jx] * HCJy[jy] *
                                                        HCJz[jz] * HCLx[lx] *
                                                        HCLy[ly] * HCLz[lz] *
                                                        coulombElement3D(
                                                                ix,iy,iz,
                                                                jx,jy,jz,
                                                                kx,ky,kz,
                                                                lx,ly,lz);
                                                } // end forlz
                                            } // end forly
                                        } // end forlx
                                    } // end forjz
                                } // end forjy
                            } // end forjx
                        } // end forkz
                    } // end forky
                } // end forkx
            } // end for iz
        } // end foriy
    } // end forix

    return sum;
} // end function coulomb3D

double GaussianIntegrals::overlapElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return the overlap integral element <i|j> */
    double prod = 1.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        prod *= ddexpr(*(GaussianBasis::Cartesian::getStates(i)(d)),
                *(GaussianBasis::Cartesian::getStates(j)(d)),
                &GaussianIntegrals::ddexprOverlap);
    } // end ford

    return prod * normalizationFactor(i) * normalizationFactor(j);
} // end function overlapElement

double GaussianIntegrals::oneBodyElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return oneBodyElement <i|h|k> = <i|K|j> + <i|P|j>, where K
     * is the kinetic part and P is the potential part */

    return -0.5*laplacianElement(i,j) + potentialElement(i,j);
} // end function oneBodyElements

double GaussianIntegrals::coulombElement(const unsigned int& i, const unsigned
    int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> */
//     int nixMod = (*(GaussianBasis::Cartesian::getStates(i)(0)))%2;
//     int niyMod = (*(GaussianBasis::Cartesian::getStates(i)(1)))%2;
//     int nkxMod = (*(GaussianBasis::Cartesian::getStates(k)(0)))%2;
//     int nkyMod = (*(GaussianBasis::Cartesian::getStates(k)(1)))%2;
//     int njxMod = (*(GaussianBasis::Cartesian::getStates(j)(0)))%2;
//     int njyMod = (*(GaussianBasis::Cartesian::getStates(j)(1)))%2;
//     int nlxMod = (*(GaussianBasis::Cartesian::getStates(l)(0)))%2;
//     int nlyMod = (*(GaussianBasis::Cartesian::getStates(l)(1)))%2;
//     bool integrandIsEven = ((nixMod == nkxMod) && (niyMod == nkyMod)) ==
//         ((njxMod == nlxMod) && (njyMod == nlyMod));

    const Eigen::Array<int*, Eigen::Dynamic, 1>& ni =
        GaussianBasis::Cartesian::getStates(i).segment(0, m_dim);
    const Eigen::Array<int*, Eigen::Dynamic, 1>& nj =
        GaussianBasis::Cartesian::getStates(j).segment(0, m_dim);
    const Eigen::Array<int*, Eigen::Dynamic, 1>& nk =
        GaussianBasis::Cartesian::getStates(k).segment(0, m_dim);
    const Eigen::Array<int*, Eigen::Dynamic, 1>& nl =
        GaussianBasis::Cartesian::getStates(l).segment(0, m_dim);

    Eigen::ArrayXi nSum = Eigen::ArrayXi::Zero(m_dim);
    Methods::refSum(nSum, ni, nj, nk, nl);
    nSum = nSum.unaryExpr([](const int m) {return m%2;});

    bool integrandIsEven = false;
    for (unsigned int i = 0; i < nSum.size(); ++i) {
        for (unsigned int j = 0; j < nSum.size(); ++j) {
            if (nSum(i) == nSum(j)) {
                integrandIsEven = true;
            } else { 
                integrandIsEven = false;
                break;
            } // end if
        } // end forj
    } // end fori

    if (integrandIsEven) {
        /* make sure integrand is even (odd integrand yields zero) */
        return normalizationFactors(i) * normalizationFactors(j) *
            normalizationFactors(k) * normalizationFactors(l) *
            (this->*coulombFunc)(i,j,k,l);
    } else {
        return 0.0;
    } // end ifselse
} // end function coulombElement
