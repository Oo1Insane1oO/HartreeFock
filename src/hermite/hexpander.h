#ifndef HEXPANDER_H#define HEXPANDER_H

#include <cmath>

#include <Eigen/Dense>

using EigenMatArrXd = Eigen::Matrix<Eigen::ArrayXd, Eigen::Dynamic,
      Eigen::Dynamic>;
// using EigenArrMatXd = Eigen::Array<Eigen::MatrixXd, Eigen::Dynamic, 1>;
// using EigenArrCubeXd = Eigen::Array<Eigen::Array<Eigen::MatrixXd,
//       Eigen::Dynamic, 1>, Eigen::Dynamic, 1>;

class HexpanderCoeffs {
    private:
        bool checkIndices(const int&, const int&, const int&);

    protected:
        EigenMatArrXd coefficients;

        HexpanderCoeffs() {
        } // end constructor
        virtual ~HexpanderCoeffs() {
        } // end deconstructor

    public:
        const double& coeff(const unsigned int&, const unsigned int&, const
                unsigned int&) const;

        void setCoefficients(unsigned int, unsigned int, double, double,
                double);
}; // end class HexpanderCoeffs

template<const unsigned int XM, const unsigned int YM, 
    const double A, const double B, const double C, const double D, 
    const std::array<double, 2> PQ, 
    const unsigned int NM=XM+YM>
class Hexpander2D : public HexpanderCoeffs {
    private:
        static constexpr auto I2Dvalues = []() {
            /* 1D integral values */

            constexpr auto modifiedIntegrand = [](double u, const
                    unsigned int n, const double pRR) {
                /* modified integrand (using Gauss-Chebyshev method) */
                double powu = 1.0;
                for (unsigned int i = 0; i < 2*n; ++i) {
                    powu *= u;
                } // end fori

                return powu * exp(-u*u*pRR);
            };

            std::array<std::array<std::array<double, XM+1>, YM+1>, NM+1>
                values{};
            double p = (A+C)*(B+D) / (A+C+B+D);
            double powVal = 1;
            double s = 50;
            double pis = M_PI/s; // all weights are pi/n
            for (unsigned int n = 0; n <= NM; ++n) {
                /* calculate initial integrals */
                values[n][0][0] = powVal; 
                double tmpVal = 0.0;
                for (unsigned int i = 0; i < s; ++i) {
                    tmpVal += modifiedIntegrand(cos(pis * (i-0.5)), n,
                            p*(PQ[0]*PQ[0] + PQ[1]*PQ[1]));
                } // end fori
                values[n][0][0] *= tmpVal * pis;
                powVal *= -2*p;
            } // end forn

            for (unsigned int ts = 1; ts <= NM; ++ts) {
                for (unsigned int n = 0; n <= NM-ts; ++n) {
                    for (unsigned int i = 0; i <= XM; ++i) {
                        for (unsigned int j = 0; j <= YM; ++j) {
                            if ((i+j != ts) || (i+j == 0)) {
                                /* out of bounds */
                                continue;
                            } // end if
                            unsigned int ijMax = (i > j ? i : j);
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
        }(); // end lambda I2Dvalues

    public:
        Hexpander2D() {
        } // end constructor
        virtual ~Hexpander2D (){
        } // end deconstructor

//         void setAuxiliary2D(unsigned int, unsigned int, double, double, double,
//                 double, const Eigen::VectorXd&);

        const double& auxiliary2D(const unsigned int& n, const unsigned int& i,
                const unsigned int& j) const {
            /* return integral R^{ij}_n */
            return I2Dvalues[n][i][j];
        } // end function integral2D
};

template<const unsigned int XM, const unsigned int YM, const unsigned int ZM,
    const double A, const double B, const double C, const double D, 
    const std::array<double, 3> PQ,
    const unsigned int NM=XM+YCM+ZM>
class Hexpander3D : public HexpanderCoeffs {
    private:
        static constexpr auto I3Dvalues = []() {
            /* 1D integral values */

            constexpr auto boysIntegrand = [](double u, const unsigned int n,
                    const double pRR) {
                /* Integrand of Boys function (using Gauss-Chebyshev method) */
                double powu = 1;
                for (unsigned int i = 0; i < 2*n; ++i) {
                    powu *= u;
                } // end fori

                return powu * sqrt(1 - u*u) * exp(-u*u*pRR);
            };

            std::array<std::array<std<<array<std::array<double, XM+1>, YM+1>,
            ZM+1>, NM+1> values{}; 
            double p = (A+C)*(B+D) / (A+C+B+D);
            double powVal = 1;
            double s = 50;
            double pis = M_PI/s; // all weights are pi/n
            for (unsigned int n = 0; n < nMax; ++n) {
                values[n][0][0][0] = powVal;
                double tmpVal = 0.0;
                for (unsigned int i = 0; i < s; ++i) {
                    tmpVal += boysIntegrand(cos(pis * (i-0.5)), n,
                            p*(PQ[0]*PQ[0] + PQ[1]*PQ[1] + PQ[2]*PQ[2]));
                } // end fori
                values[n][0][0][0] *= tmpVal * pis;
                powVal *= -2*p;
            } // end forn

            for (unsigned int ts = 1; ts <= NM; ++ts) {
                for (unsigned int n = 0; n <= NM - ts; ++n) {
                    for (unsigned int i = 0; i <= XM; ++i) {
                        for (unsigned int j = 0; j <= YM; ++j) {
                            for (unsigned int k = 0; k <= ZM; ++k) {
                                if (i+j+k != ts || i+j+k == 0) {
                                    /* out of bounds */
                                    continue;
                                } // end if
                                unsigned int ijkMax = ((i>j) ? ((i>k) ? i : k)
                                        : ((j>k) ? j : k));
                                int i2 = i;
                                int j2 = j;
                                int k2 = k;
                                int i3 = i;
                                int j3 = j;
                                int k3 = k;
                                int factor = i;
                                double XPQ = PQ[0];
                                if (ijkMax == i) {
                                    i2 = i - 2;
                                    i3 = i - 1;
                                    factor = i3;
                                    XPQ = PQ[0];
                                } else if (ijkMax == j) {
                                    j2 = j - 2;
                                    j3 = j - 1;
                                    factor = j3;
                                    XPQ = PQ[1];
                                } else {
                                    k2 = k - 2;
                                    k3 = k - 1;
                                    factor = k3;
                                    XPQ = PQ[2];
                                } // end ifeifelse
                                double  currValue = 0;
                                if (i2 >= 0 && j2 >= 0 && k2 >= 0) {
                                    currValue += factor * values[n+1][i2][j2][k2];
                                } // end if
                                if (i3 >= 0 && j3 >= 0 && k3 >= 0) {
                                    currValue += XPQ * values[n+1][i3][j3][k3];
                                } // end if
                                values[n][i][j][k] = currValue;
                            } // end fork
                        } // end forj
                    } // end fori
                } // end forn
            } // end forts
        }(); // end lambda I3Dvalues

    public:
        Hexpander3D() {
        } // end constructor
        virtual ~Hexpander3D (){
        } // end deconstructor
        
//         void setAuxiliary3D(unsigned int, unsigned int, unsigned int, double,
//                 double, double, double, const Eigen::VectorXd&);

        const double& auxiliary3D(const unsigned int& n, const unsigned int& i,
                const unsigned int& j, const unsigned int& k) const {
            /* return integral R^{ij}_n */
            return I3Dvalues[n][i][j][k];
        } // end function integral2D
};

#endif /* HEXPANDER_H */
