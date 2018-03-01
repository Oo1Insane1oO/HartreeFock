#include "gaussianintegrals.h"
#include "../hermite/hermite.h"
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>

GaussianIntegrals::GaussianIntegrals(const unsigned int dim, unsigned int
        cutOff, double scaling) : GaussianBasis(cutOff, dim) {
    m_dim = dim;
    expScaleFactor = scaling;
    sqrtFactor = sqrt(scaling);
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
    sqrtScale1 = sqrt(pow(xScale, m_dim));
    sqrtScale = 1./sqrtScale1;
    powScale = pow(xScale, 2*m_dim);

    // choose coulombElement function for 2D or 3D
    if (m_dim == 2) {
        coulombElementFunc = &GaussianIntegrals::coulombElement2D;
    } else {
        coulombElementFunc = &GaussianIntegrals::coulombElement3D;
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

    return 0.5*xScale*sums * normalizationFactor(i) * normalizationFactor(j);
} // end function potentialElement

inline double GaussianIntegrals::coulombElement2D(const unsigned int& i, const
        unsigned int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> in 2D for level (i,k,j,l) (calculate 1D integral
     * numerically) */
    double sum = 0.0;
    int Ix = *(GaussianBasis::Cartesian::getStates(i)(0));
    int Iy = *(GaussianBasis::Cartesian::getStates(i)(1));
    int Kx = *(GaussianBasis::Cartesian::getStates(k)(0));
    int Ky = *(GaussianBasis::Cartesian::getStates(k)(1));
    int Jx = *(GaussianBasis::Cartesian::getStates(j)(0));
    int Jy = *(GaussianBasis::Cartesian::getStates(j)(1));
    int Lx = *(GaussianBasis::Cartesian::getStates(l)(0));
    int Ly = *(GaussianBasis::Cartesian::getStates(l)(1));
    for (int ix = 0; ix < Ix; ++ix) {
        for (int iy = 0; iy < Iy; ++iy) {
            for (int jx = 0; jx < Jx; ++jx) {
                for (int jy = 0; jy < Jy; ++jy) {
                    for (int kx = 0; kx < Kx; ++kx) {
                        for (int ky = 0; ky < Ky; ++ky) {
                            for (int lx = 0; lx < Lx; ++lx) {
                                for (int ly = 0; ly < Ly; ++ly) {
                                    sum += HC(Ix)[ix]*HC(Iy)[iy] *
                                        HC(Jx)[jx]*HC(Jy)[jy] *
                                        HC(Kx)[kx]*HC(Ky)[ky] *
                                        HC(Lx)[lx]*HC(Ly)[ly];
                                    // TODO: fix numerical part
                                } // end forly
                            } // end forlx
                        } // end forky
                    } // end forkx
                } // end forjy
            } // end forjx
        } // end foriy
    } // end foriy
    return sum;
} // end function coulombElement2D

inline double GaussianIntegrals::coulombElement3D(const unsigned int& i, const
        unsigned int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> in 3D for level (i,k,j,l) (calculate 1D integral
     * numerically) */
    double sum = 0.0;
    int Ix = *(GaussianBasis::Cartesian::getStates(i)(0));
    int Iy = *(GaussianBasis::Cartesian::getStates(i)(1));
    int Iz = *(GaussianBasis::Cartesian::getStates(i)(2));
    int Kx = *(GaussianBasis::Cartesian::getStates(k)(0));
    int Ky = *(GaussianBasis::Cartesian::getStates(k)(1));
    int Kz = *(GaussianBasis::Cartesian::getStates(k)(2));
    int Jx = *(GaussianBasis::Cartesian::getStates(j)(0));
    int Jy = *(GaussianBasis::Cartesian::getStates(j)(1));
    int Jz = *(GaussianBasis::Cartesian::getStates(z)(2));
    int Lx = *(GaussianBasis::Cartesian::getStates(l)(0));
    int Ly = *(GaussianBasis::Cartesian::getStates(l)(1));
    int Lz = *(GaussianBasis::Cartesian::getStates(l)(2));
    for (int ix = 0; ix < Ix; ++ix) {
        for (int iy = 0; iy < Iy; ++iy) {
            for (int iz = 0; iz < Iz; ++iz) {
                for (int jx = 0; jx < Jx; ++jx) {
                    for (int jy = 0; jy < Jy; ++jy) {
                        for (int jz = 0; jz < Jz; ++jz) {
                            for (int kx = 0; kx < Kx; ++kx) {
                                for (int ky = 0; ky < Ky; ++ky) {
                                    for (int kz = 0; kz < Kz; ++kz) {
                                        for (int lx = 0; lx < Lx; ++lx) {
                                            for (int ly = 0; ly < Ly; ++ly) {
                                                for (int lz = 0; lz < Lz; ++lz)
                                                {
                                                    sum +=
                                                        HC(Ix)[ix]*HC(Iy)[iy] *
                                                        HC(Jx)[jx]*HC(Jy)[jy] *
                                                        HC(Kx)[kx]*HC(Ky)[ky] *
                                                        HC(Lx)[lx]*HC(Ly)[ly];
                                                    // TODO: fix numerical part
                                                } // end forlz
                                            } // end forly
                                        } // end forlx
                                    } // end forkz
                                } // end forky
                            } // end forkx
                        } // end forjz
                    } // end forjy
                } // end forjx
            } // end foriz
        } // end foriy
    } // end forix
    return sum;
} // end function coulombElement3D

double GaussianIntegrals::overlapElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return the overlap integral element <i|j> */
    if (*(GaussianBasis::Cartesian::getStates(i)(m_dim+1)) !=
            *(GaussianBasis::Cartesian::getStates(j)(m_dim+1))) {
        /* respect spin orthogonality */
        return 0.0;
    } // end if
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
    if (*(GaussianBasis::Cartesian::getStates(i)(m_dim+1)) !=
            *(GaussianBasis::Cartesian::getStates(j)(m_dim+1))) {
        /* respect spin orthogonality */
        return 0.0;
    } // end if

    return -0.5*laplacianElement(i,j) + potentialElement(i,j);
} // end function oneBodyElements

double GaussianIntegrals::coulombElement(const unsigned int& i, const unsigned
        int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> */
    if (((GaussianBasis::Cartesian::getSumn(i) +
                    GaussianBasis::Cartesian::getSumn(j)) !=
                (GaussianBasis::Cartesian::getSumn(k) +
                 GaussianBasis::Cartesian::getSumn(l))) ||
            ((*(GaussianBasis::Cartesian::getStates(i)(m_dim+1)) +
              *(GaussianBasis::Cartesian::getStates(j)(m_dim+1))) !=
             *((GaussianBasis::Cartesian::getStates(k)(m_dim+1)) +
                 *(GaussianBasis::Cartesian::getStates(l)(m_dim)))) ||
            (*(GaussianBasis::Cartesian::getStates(i)(m_dim+1)) !=
             *(GaussianBasis::Cartesian::getStates(k)(m_dim+1))) ||
            (*(GaussianBasis::Cartesian::getStates(j)(m_dim+1)) !=
             *(GaussianBasis::Cartesian::getStates(l)(m_dim+1)))) {
        /* make sure total angular momentum and total spin is conserved */
        return 0.0;
    } // end if

    return (this->*coulombElementFunc)(i,j,k,l);
} // end function coulombElement
