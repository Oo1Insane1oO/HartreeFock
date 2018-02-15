#ifndef GAUSSIANPRIMITIVEBASIS_H
#define GAUSSIANPRIMITIVEBASIS_H

#include <Eigen/Dense>

class GaussianPrimitiveBasis {
    private:
        double m_weight;
        Eigen::VectorXi m_exponentVector;
        Eigen::VectorXd m_scalingVector;
    
    public:
        GaussianPrimitiveBasis ();
        GaussianPrimitiveBasis (Eigen::VectorXd, Eigen::VectorXi);
        virtual ~GaussianPrimitiveBasis ();

        int dExponent(const unsigned int&);
        const Eigen::VectorXi& expVec() const;
        const Eigen::VectorXd& scalingVec() const;
        void setExponent(int, const unsigned int);
        void setExponent(const Eigen::VectorXi&);

        void setScaling(const Eigen::VectorXd&);

        int exponentMax() const;
        int angularMomentum() const;
};

#endif /* GAUSSIANPRIMITIVEBASIS_H */
