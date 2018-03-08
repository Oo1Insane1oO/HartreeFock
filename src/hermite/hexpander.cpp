#include "hexpander.h"

#include "../gaussianquadrature.h"
#include "../methods.h"

#include <cmath>
#include <iostream>

#include <boost/math/special_functions/gamma.hpp>

Hexpander::Hexpander() {
    /* default constructor */
} // end constructor

Hexpander::~Hexpander() {
} // end deconstructor

void Hexpander::setCoefficients(unsigned int iMax, unsigned int jMax, double a,
        double b, double Qx) {
    /* find coefficients */
    double p = a + b;
    double q = a*b/p;

    unsigned int tMax = iMax+jMax+1;
    coefficients = EigenMatArrXd::Constant(iMax+1, jMax+1,
            Eigen::ArrayXd::Zero(tMax));
    coefficients(0,0)(0) = exp(-q*Qx*Qx);
    for (unsigned int i = 0; i < coefficients.cols(); ++i) {
        for (unsigned int t = 0; t < tMax; t++) {
            if ((i==0) && (t==0)) {
                continue;
            } // end if
            int im = i - 1;
            int tm = t - 1;
            int tp = t + 1;

            double Em = 0;
            if(checkIndices(0, im, tm)) {
                Em = coefficients(0,im)(tm);
            } // end if
            double Ed = 0;
            if(checkIndices(0, im, t)) {
                Ed = coefficients(0,im)(t);
            } // end if
            double Ep = 0;
            if(checkIndices(0, im, tp)) {
                Ep = coefficients(0,im)(tp);
            } // end if
            coefficients(0,i)(t) = 0.5/p*Em + q*Qx/b*Ed + tp*Ep;
        } // end fort
    } // end fori
    for (unsigned int i = 1; i < coefficients.rows(); ++i) {
        for (unsigned int j = 0; j < coefficients.cols(); ++j) {
            for (unsigned int t = 0; t < tMax; ++t) {
                int im = i - 1;
                int tm = t - 1;
                int tp = t + 1;

                double Em = 0;
                if(checkIndices(im, j, tm)) {
                    Em = coefficients(im,j)(tm);
                } // end if
                double Ed = 0;
                if(checkIndices(im, j, t)) {
                    Ed = coefficients(im,j)(t);
                } // end if
                double Ep = 0;
                if(checkIndices(im, j, tp)) {
                    Ep = coefficients(im,j)(tp);
                } // end if
                coefficients(i,j)(t) = 0.5/p*Em + q*Qx/a*Ed + tp*Ep;
            } // end fort
        } // end forj
    } // end fori
} // end function setCoefficients

void Hexpander::setAuxiliary2D(unsigned int xMax, unsigned int yMax, double a,
        double b, double c, double d, const Eigen::VectorXd& PQ) {
    /* setup integral elements in 2D */
    unsigned int nMax = xMax + yMax;
    double p = (a+c)*(b+d) / (a+c+b+d);
    integrals = EigenArrMatXd::Constant(nMax+1, Eigen::MatrixXd::Zero(xMax+1,
                yMax+1));
    double powVal = 1;
    for (unsigned int n = 0; n <= nMax; ++n) {
        /* calculate initial integrals */
        integrals(n)(0,0) = powVal * GaussianQuadrature::gaussChebyshevQuad(50,
                this, &Hexpander::modifiedIntegrand, n, p*PQ.squaredNorm());
        powVal *= -2*p;
    } // end forn

    for (unsigned int ts = 1; ts <= nMax; ++ts) {
        for (unsigned int n = 0; n <= nMax-ts; ++n) {
            for (unsigned int i = 0; i <= xMax; ++i) {
                for (unsigned int j = 0; j <= yMax; ++j) {
                    if (i+j != ts || i+j == 0) {
                        /* out of bounds */
                        continue;
                    } // end if
                    unsigned int ijMax = Methods::max(i,j);
                    int i2 = i;
                    int j2 = j;
                    int i3 = i;
                    int j3 = j;
                    int factor = i;
                    double XPQ = PQ(0);
                    if (ijMax == i) {
                        i2 = i-2;
                        i3 = i-1;
                        factor = i3;
                        XPQ = PQ(0);
                    } else {
                        j2 = j-2;
                        j3 = j-1;
                        factor = j3;
                        XPQ = PQ(1);
                    } // end ifelse
                    double currValue = 0;
                    if (i2 >= 0 && j2 >= 0) {
                        currValue += factor * integrals(n+1)(i2,j2);
                    } // end if
                    if (i3 >= 0 && j3 >= 0) {
                        currValue += XPQ * integrals(n+1)(i3,j3);
                    } // end if
                    integrals(n)(i,j) = currValue;
                } // end forj
            } // end fori
        } // end forn
    } // end forts
} // end function setIntegrals

bool Hexpander::checkIndices(const int& ia, const int& ib, const int& t) {
    if (t < 0 || t > (ia+ib) || ia < 0 || ib < 0) {
        return false;
    } else {
        return true;
    } // end if
} // end function checkIndices

double Hexpander::boysIntegrand(double u, const unsigned int& n, const double&
        pRR) {
    /* Integrand of Boys function */
    return pow(u, 2*n) * exp(-u*u*pRR);
} // end function boysIntegrand

double Hexpander::modifiedIntegrand(double u, const unsigned int& n, const
        double& pRR) {
    /* integrand assuming Gauss-Chebyshev method is used (the 1/sqrt(1-u^2)
     * part is absorbed) */
    return pow(u, 2*n) * exp(-u*u*pRR); // * 1/sqrt(1-u*u);
} // end function modifiedIntegrand

double Hexpander::boys(const unsigned int& n, const double& pRR) {
    /* calculate boys function using boost library */
    if (fabs(pRR) <= 1e-15) {
        /* case centering is in zero */
        return 1./n;
    } else {
        /* calculate full integral */
        return 0.5/pow(pRR*pRR, n+0.5) *
            boost::math::tgamma_lower<double>(n+0.5, pRR*pRR);
    } // end ifelse
} // end function boys

double Hexpander::auxiliary3D(const unsigned int& ix, const unsigned int& iy,
        const unsigned int& iz, const unsigned int& n, const double& p,
        const Eigen::VectorXd& P, const double& R) {
    /* calculate auxiliary integral (Boys function) */
    double val = 0.0;
    if ((ix==0) && (iy==0) && (iz==0)) {
        val += pow(-2*p,n) * boys(n, p*R*R);
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

const double& Hexpander::coeff(const unsigned int& i, const unsigned int& j,
        const unsigned int& t) const {
    /* return coefficient E^{ij}_t */
    return coefficients(i,j)(t);
} // end coeff

const double& Hexpander::auxiliary2D(const unsigned int& n, const unsigned int&
        i, const unsigned int& j) {
    /* return integral R^{ij}_n */
    return integrals(n)(i,j);
} // end function integral2D
