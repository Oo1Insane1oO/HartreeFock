#ifndef HEXPANDER_H
#define HEXPANDER_H

#include <Eigen/Dense>

using EigMatMatXd = Eigen::Matrix<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic>;

class Hexpander {
    private:
        EigMatMatXd coeffs;

        double calculateCoeff(const int&, const int&, const int&, const
                double&, const double&, const double&);

    public:
        Hexpander();
        Hexpander(const size_t&, const size_t&, const double&, const double&,
                const Eigen::VectorXd&, const Eigen::VectorXd&);
        virtual ~Hexpander ();
        
        void setup(const size_t&, const size_t&, const double&, const double&,
                const Eigen::VectorXd&, const Eigen::VectorXd&);

        const double& coeff(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&) const;
};

#endif /* HEXPANDER_H */
