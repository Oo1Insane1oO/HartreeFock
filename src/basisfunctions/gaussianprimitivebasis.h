#ifndef GAUSSIANPRIMITIVEBASIS_H
#define GAUSSIANPRIMITIVEBASIS_H

#include <Eigen/Dense>

class GaussianPrimitiveBasis {
    private:
        double m_weight;
        Eigen::VectorXi m_exponentVector;
    
    public:
        GaussianPrimitiveBasis ();
        GaussianPrimitiveBasis (double, Eigen::VectorXi);
        virtual ~GaussianPrimitiveBasis ();

        int dExponent(const unsigned int&) const;
        const Eigen::VectorXi &expVec() const;
        void setExponent(int, const unsigned int);
        void setExponent(const Eigen::VectorXi&);

        void setWeight(double);

        int exponentMax() const;
        int angularMomentum() const;
};

#endif /* GAUSSIANPRIMITIVEBASIS_H */
