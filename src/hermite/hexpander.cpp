#include "hexpander.h"
#include <cmath>
#include <iostream>

Hexpander::Hexpander() {
    /* default constructor */
} // end constructor

Hexpander::Hexpander(const size_t& N, const size_t& dim, const double& a, const
        double& b, const Eigen::VectorXd& A, const Eigen::VectorXd& B) {
    /* grab maximum angular momentum N (aka number of states) and fill
     * coefficients matrix */
    setup(N, dim, a, b, A, B);
} // end constructor

Hexpander::~Hexpander() {
} // end deconstructor

void Hexpander::setup(const size_t& N, const size_t& dim, const double& a,
        const double& b, const Eigen::VectorXd& A, const Eigen::VectorXd& B) {
    /* calculate and set coefficients in coefficient matrix */

    Eigen::VectorXd diffAB = A - B;
    coeffs = Eigen::Matrix<Eigen::MatrixXd, Eigen::Dynamic,
           Eigen::Dynamic>(N,N);
    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
            coeffs(i,j) = Eigen::MatrixXd(i+j+1,dim);
            for (unsigned int t = 0; t < coeffs(i,j).rows(); ++t) {
                for (unsigned int d = 0; d < dim; ++d) {
                    coeffs(i,j)(t,d) = calculateCoeff(i, j, t, a, b,
                            diffAB(d));
                } // end ford
            } // end fort
        } // end forj
    } // end fori
} // end function setup

double Hexpander::calculateCoeff(const int& i, const int& j, const int& t,
        const double& a, const double& b, const double& Qx) {
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
        return 0.5/p * calculateCoeff(i-1,j,t-1,a,b,Qx) - q*Qx/a *
            calculateCoeff(i-1,j,t,a,b,Qx) + (t+1)
            * calculateCoeff(i-1,j,t+1,a,b,Qx);
    } else {
        /* decrement level j */
        return 0.5/p * calculateCoeff(i,j-1,t-1,a,b,Qx) - q*Qx/b *
            calculateCoeff(i,j-1,t,a,b,Qx) + (t+1)
            * calculateCoeff(i,j-1,t+1,a,b,Qx);
    } // end ifeifeifelse
} // end function calculateCoeff

const double& Hexpander::coeff(const size_t& i, const size_t& j, const size_t&
        t, const size_t& d) const {
    /* get coefficient for t'th Hermite function in dimension d at level i+j */
    return coeffs(i,j)(t,d);
} // end function getCoefficient
