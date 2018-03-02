#ifndef HEXPANDER_H
#define HEXPANDER_H

#include <Eigen/Dense>

using EigMatMatXd = Eigen::Matrix<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic>;

class Hexpander {
    private:
        EigMatMatXd coeffs;
        double boys(const unsigned int&, const double&);
        double modified(const unsigned int&, const double&);

        double boysIntegrand(double, const unsigned int&, const double&);
        double modifiedIntegrand(double, const unsigned int&, const double&);

        template<typename U, typename F, typename... Args> inline double
            simpsons(const double& limMin, const double &limMax, const U& u, F
                    f, Args...  args) {
            /* evaluate the definite integral of function f with Simpson's rule
             * */
            return (limMax - limMin) / 6. * ((u->*f)(limMin, args...) +
                    4*(u->*f)(0.5*(limMin+limMax), args...) + (u->*f)(limMax,
                        args...));
        } // end function simpsons

    public:
        Hexpander();
        Hexpander(const size_t&, const size_t&, const double&, const double&,
                const Eigen::VectorXd&, const Eigen::VectorXd&);
        virtual ~Hexpander ();
        
        double coeff(const int&, const int&, const int&, const double&, const
                double&, const double&);
        double auxiliary2D(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&, const Eigen::VectorXd&,
                const double&);
        double auxiliary3D(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&, const unsigned int&, const
                Eigen::VectorXd&, const double&);
};

#endif /* HEXPANDER_H */
