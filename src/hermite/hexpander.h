#ifndef HEXPANDER_H
#define HEXPANDER_H

#include <Eigen/Dense>

using EigenMatArrXd = Eigen::Matrix<Eigen::ArrayXd, Eigen::Dynamic,
      Eigen::Dynamic>;
using EigenArrMatXd = Eigen::Array<Eigen::MatrixXd, Eigen::Dynamic, 1>;
using EigenArrCubeXd = Eigen::Array<Eigen::Array<Eigen::MatrixXd,
      Eigen::Dynamic, 1>, Eigen::Dynamic, 1>;

class Hexpander {
    private:
        EigenArrMatXd integrals2D;
        EigenArrCubeXd integrals3D;

        Eigen::ArrayXd coefficients;
    
        unsigned int iMp1, jMp1, nMp1;

        double boys(const unsigned int&, const double&);

        double boysIntegrand(double, const unsigned int&, const double&);
        double modifiedIntegrand(double, const unsigned int&, const double&);

        bool checkIndices(const int&, const int&, const int&);

        inline unsigned int cidx(const unsigned int& i, const unsigned int& j,
                const unsigned int& n) const {
            /* calculate index in coefficients matrix */
            return n + nMp1 * (j + jMp1 * i);
        } // end function cidx

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
