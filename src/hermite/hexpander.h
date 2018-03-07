#ifndef HEXPANDER_H
#define HEXPANDER_H

#include <Eigen/Dense>

using EigenMatArrXd = Eigen::Matrix<Eigen::ArrayXd, Eigen::Dynamic,
      Eigen::Dynamic>;
using EigenArrMatXd = Eigen::Array<Eigen::MatrixXd, Eigen::Dynamic, 1>;

class Hexpander {
    private:
        EigenMatArrXd coefficients;
        EigenArrMatXd integrals;

        double boys(const unsigned int&, const double&);

        double boysIntegrand(double, const unsigned int&, const double&);
        double modifiedIntegrand(double, const unsigned int&, const double&);

        bool checkIndices(const int&, const int&, const int&);

    public:
        Hexpander();
        virtual ~Hexpander ();
        
        void setCoefficients(unsigned int, unsigned int, double, double, double);
        void setAuxiliary2D(unsigned int, unsigned int, double, double, double,
                double, const Eigen::VectorXd&);

        const double& coeff(const unsigned int&, const unsigned int&, const
                unsigned int&) const;
        const double& auxiliary2D(const unsigned int&, const unsigned int&,
                const unsigned int&);
        
        double auxiliary2D(const unsigned int&, const unsigned int&, const
                unsigned int&, const double&, const Eigen::VectorXd&, const
                double&);
        double auxiliary3D(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&, const double&, const
                Eigen::VectorXd&, const double&);
};

#endif /* HEXPANDER_H */
