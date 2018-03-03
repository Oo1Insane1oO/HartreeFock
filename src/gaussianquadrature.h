#ifndef GAUSSIANQUADRATURE_H
#define GAUSSIANQUADRATURE_H

#include "hermite/hermite.h"

#include <boost/math/special_functions/factorials.hpp>
#include <Eigen/Dense>
#include <iostream>

class GaussianQuadrature {
    private:
        GaussianQuadrature() {};
        virtual ~GaussianQuadrature() {};

        static void setWeightsAndPoints(const size_t&, Eigen::VectorXd&,
                Eigen::VectorXd&);
        static void setPoints(const size_t&, Eigen::VectorXd&);
        static void generateRandom(double*);
        static void findRootsHermite(const size_t&, Eigen::VectorXd&);

    public:
        template<typename U, typename F, typename... Args> static inline double
            gaussHermiteQuad(const size_t n, const U& obj, F f, Args... args) {
            /* calculate integral of f using Gauss-Hermite Quadrature */
            Eigen::VectorXd points, weights;
            
            /* set points, weights and normalization constant */
            setWeightsAndPoints(n, weights, points);
            double A = 1./(pow(2,n) * boost::math::factorial<double>(n) *
                    sqrt(M_PI));

            double sum = 0;
            for (unsigned int i = 0; i < points.size(); ++i) {
                double hi = H(points(i), i);
                sum += weights(i) * A * hi*hi * (obj->*f)(points(i), args...);
            } // end fori

            return sum;
        } // end function gaussHermiteQuad
        
        template<typename U, typename F, typename... Args> static inline double
            gaussChebyshevQuad(const size_t n, const U& obj, F f, Args... args) {
            /* calculate integral of f using Gauss-Chebyshev Quadrature */
            Eigen::ArrayXd points(n), weights(n);
            for (unsigned int i = 0; i < n; ++i) {
                weights(i) = M_PI/n;
                points(i) = cos(weights(i) * (i-0.5));
            } // end fori

            double sum = 0.0;
            for (unsigned int i = 0; i < points.size(); ++i) {
                sum += weights(i) * (obj->*f)(points(i), args...);
            } // end fori

            return sum;
        } // end function gaussChebyshevQuad
};

#endif /* GAUSSIANQUADRATURE_H */
