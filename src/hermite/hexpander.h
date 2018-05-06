#ifndef HEXPANDER_H
#define HEXPANDER_H

#include <cmath>

#include <Eigen/Dense>

#include "../gaussianquadrature.h"
#include "../methods.h"

using EigenMatArrXd = Eigen::Matrix<Eigen::ArrayXd, Eigen::Dynamic,
      Eigen::Dynamic>;
using EigenArrMatXd = Eigen::Array<Eigen::MatrixXd, Eigen::Dynamic, 1>;
using EigenArrCubeXd = Eigen::Array<Eigen::Array<Eigen::MatrixXd,
      Eigen::Dynamic, 1>, Eigen::Dynamic, 1>;

class Hexpander {
    private:
        EigenMatArrXd coefficients;
        EigenArrMatXd integrals2D;
        EigenArrCubeXd integrals3D;

        double boys(const unsigned int&, const double&);

        double boysIntegrand(double, const unsigned int&, const double&);
        double modifiedIntegrand(double, const unsigned int, const double);

        static constexpr auto modifiedIntegrandStatic = [](double u, const
                unsigned int n, const double pRR) {
            double powu = 1.0;
            for (unsigned int i = 0; i < 2*n; ++i) {
                powu *= u;
            } // end fori

            return powu * exp(-u*u*pRR);
        };

        template<const unsigned int xMax, const unsigned int yMax, const
            unsigned int nMax=xMax+yMax>
        static constexpr std::array<std::array<std::array<double, xMax+1>,
                         yMax+1>, nMax+1> I2Dvalues(const double a, const
                                 double b, const double c, const double d,
                                 const std::array<double,2> PQ) {
//             unsigned int nMax = xMax + yMax;
            double p = (a+c)*(b+d) / (a+c+b+d);
            std::array<std::array<std::array<double, xMax+1>, yMax+1>, nMax+1>
                values{};
            double powVal = 1;
            double s = 50;
            double pis = M_PI/s; // all weights are pi/n
            for (unsigned int n = 0; n <= nMax; ++n) {
                /* calculate initial integrals */
                values[n][0][0] = powVal; 
                double tmpVal = 0.0;
                for (unsigned int i = 0; i < s; ++i) {
                    tmpVal += modifiedIntegrandStatic(cos(pis * (i-0.5)),
                            n, p*(PQ[0]*PQ[0] + PQ[1]*PQ[1]));
                } // end fori
                values[n][0][0] *= tmpVal * pis;
                powVal *= -2*p;
            } // end forn

            for (unsigned int ts = 1; ts <= nMax; ++ts) {
                for (unsigned int n = 0; n <= nMax-ts; ++n) {
                    for (unsigned int i = 0; i <= xMax; ++i) {
                        for (unsigned int j = 0; j <= yMax; ++j) {
                            if ((i+j != ts) || (i+j == 0)) {
                                /* out of bounds */
                                continue;
                            } // end if
                            unsigned int ijMax = Methods::max(i,j);
                            int i2 = i;
                            int j2 = j;
                            int i3 = i;
                            int j3 = j;
                            int factor = i;
                            double XPQ = PQ[0];
                            if (ijMax == i) {
                                i2 = i-2;
                                i3 = i-1;
                                factor = i3;
                                XPQ = PQ[0];
                            } else {
                                j2 = j-2;
                                j3 = j-1;
                                factor = j3;
                                XPQ = PQ[1];
                            } // end ifelse
                            double currValue = 0;
                            if (i2 >= 0 && j2 >= 0) {
                                currValue += factor * values[n+1][i2][j2];
                            } // end if
                            if (i3 >= 0 && j3 >= 0) {
                                currValue += XPQ * values[n+1][i3][j3];
                            } // end if
                            values[n][i][j] = currValue;
                        } // end forj
                    } // end fori
                } // end forn
            } // end forts

            return values;
        }

        bool checkIndices(const int&, const int&, const int&);

    public:
        Hexpander();
        virtual ~Hexpander ();
        
        void setCoefficients(unsigned int, unsigned int, double, double, double);
        void setAuxiliary2D(unsigned int, unsigned int, double, double, double,
                double, const Eigen::VectorXd&);
        void setAuxiliary3D(unsigned int, unsigned int, unsigned int, double,
                double, double, double, const Eigen::VectorXd&);

        const double& coeff(const unsigned int&, const unsigned int&, const
                unsigned int&) const;
        const double& auxiliary2D(const unsigned int&, const unsigned int&,
                const unsigned int&) const;
        const double& auxiliary3D(const unsigned int&, const unsigned int&,
                const unsigned int&, const unsigned int&) const;
};

#endif /* HEXPANDER_H */
