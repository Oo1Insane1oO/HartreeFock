#include "hexpander.h"

#include "../gaussianquadrature.h"

#include <cmath>
#include <iostream>

#include <boost/math/special_functions/gamma.hpp>

Hexpander::Hexpander() {
    /* default constructor */
} // end constructor

Hexpander::~Hexpander() {
} // end deconstructor

double Hexpander::boysIntegrand(double u, const unsigned int& n, const double&
        pRR) {
    /* Integrand assuming Gauss-Hermite method is used (exponential part is
     * absorbed) */
    return pow(u, 2*n) * exp(-u*u*pRR); // * exp(-u*u
} // end function boysIntegrand

double Hexpander::modifiedIntegrand(double u, const unsigned int& n, const
        double& pRR) {
    /* integrand assuming Gauss-Chebyshev method is used (the 1/sqrt(1-u^2)
     * part is absorbed) */
    return pow(u, 2*n) * exp(-u*u*pRR); // * 1/sqrt(1-u*u);
} // end function modifiedIntegrand

double Hexpander::boys(const unsigned int& n, const double& pRR) {
    /* calculate boys function using boost library */
//     return GaussianQuadrature::gaussHermiteQuad(n+1, this,
//             &Hexpander::boysIntegrand, n, pRR*pRR);
    if (fabs(pRR) <= 1e-15) {
        return 1./n;
    } else {
        return 0.5/pow(pRR*pRR, n) * boost::math::tgamma_lower<double>(n+0.5,
                pRR*pRR);
    } // end ifelse
} // end function boys

double Hexpander::modified(const unsigned int& n, const double& pRR) {
    /* calculate modified function using Simpson's rule */
    return GaussianQuadrature::gaussChebyshevQuad(100, this,
            &Hexpander::modifiedIntegrand, n, pRR*pRR);
} // end function boys

double Hexpander::auxiliary3D(const unsigned int& ix, const unsigned int& iy,
        const unsigned int& iz, const unsigned int& n, const unsigned int& p,
        const Eigen::VectorXd& P, const double& R) {
    /* calculate auxiliary integral (Boys function) */
    double val = 0.0;
    if ((ix==0) && (iy==iz) && (iz==ix)) {
        val += pow(-2*p, n) * boys(n, p*R*R);
    } else if ((ix==iy) && (iy==0)) {
        if (iz > 1) {
            val += (iz-1) * auxiliary3D(ix,iy,iz-2,n+1,p,P,R);
        } // end if
        val += P(2) * auxiliary3D(ix,iy,iz-1,n+1,p,P,R);
    } else if (ix==0) {
        if (iy > 1) {
            val += (iy-1) * auxiliary3D(ix,iy-2,iz,n+1,p,P,R);
        } // end if
        val += P(1) * auxiliary3D(ix,iy-1,iz,n+1,p,P,R);
    } else {
        if (ix > 1) {
            val += (ix-1) * auxiliary3D(ix-2,iy,iz,n+1,p,P,R);
        } // end if
        val += P(0) * auxiliary3D(ix-1,iy,iz,n+1,p,P,R);
    } // end ifeifeifelse

    return val;
} // end function auxiliary3D

double Hexpander::auxiliary2D(const unsigned int& ix, const unsigned int& iy,
        const unsigned int& n, const unsigned int& p, const Eigen::VectorXd& P,
        const double& R) {
    /* calculate auxiliary integral (Modified Bessel function) */
    double val = 0.0;
    if ((ix==0) && (iy==ix)) {
        val += pow(-2*p, n) * modified(n, p*R*R);
    } else if (ix==0) {
        if (iy > 1) {
            val += (iy-1) * auxiliary2D(ix,iy-2,n+1,p,P,R);
        } // end if
        val += P(1) * auxiliary2D(ix,iy-1,n+1,p,P,R);
    } else {
        if (ix > 1) {
            val += (ix-1) * auxiliary2D(ix-2,iy,n+1,p,P,R);
        } // end if
        val += P(0) * auxiliary2D(ix-1,iy,n+1,p,P,R);
    } // end ifeifeifelse
    return val;
} // end function auxiliary3D

double Hexpander::coeff(const int& i, const int& j, const int& t, const double&
        a, const double& b, const double& Qx) {
    /* calculate Hermite coefficient */
    double p = a+b;
    double q = a*b/p;
    if ((t<0) || (t>(i+j))) {
        /* bounds for t */
        return 0.0;
    } else if (((i==0) || (j==0)) && (i==j) && (i==t)) {
        /* initial coefficient */
        return exp(-q*Qx*Qx);
    } else if (j==0) {
        /* decrement level i */
        return 0.5/p * coeff(i-1,j,t-1,a,b,Qx) - q*Qx/a * coeff(i-1,j,t,a,b,Qx)
            + (t+1)
            * coeff(i-1,j,t+1,a,b,Qx);
    } else {
        /* decrement level j */
        return 0.5/p * coeff(i,j-1,t-1,a,b,Qx) - q*Qx/b * coeff(i,j-1,t,a,b,Qx)
            + (t+1)
            * coeff(i,j-1,t+1,a,b,Qx);
    } // end ifeifeifelse
} // end function coeff
