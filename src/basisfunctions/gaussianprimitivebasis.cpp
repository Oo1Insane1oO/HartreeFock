#include "gaussianprimitivebasis.h"

GaussianPrimitiveBasis::GaussianPrimitiveBasis() {
    /* default constructor */
} // end constructor

GaussianPrimitiveBasis::GaussianPrimitiveBasis(Eigen::VectorXd scalingVec,
        Eigen::VectorXi expVec) {
    /* initialize */
    m_scalingVector = scalingVec;
    m_exponentVector = expVec;
} // end constructor

GaussianPrimitiveBasis::~GaussianPrimitiveBasis() {
} // end deconstructor

int GaussianPrimitiveBasis::dExponent(const unsigned int& dIdx) {
    /* return exponent for dimension dIdx */
    return m_exponentVector(dIdx);
} // end function dExponent

const Eigen::VectorXi& GaussianPrimitiveBasis::expVec() const {
    /* return exponent for all dimension */
    return m_exponentVector;
} // end function dExponent

const Eigen::VectorXd& GaussianPrimitiveBasis::scalingVec() const {
    /* return scaling for all dimension */
    return m_scalingVector;
} // end function dExponent

void GaussianPrimitiveBasis::setScaling(const Eigen::VectorXd& scalingVector) {
    /* set m_scalingVector to scalingVector */
    m_scalingVector = scalingVector;
} // end functionsetScaling 

void GaussianPrimitiveBasis::setExponent(int value, const unsigned int dIdx) {
    /* set exponent in dimension dIdx */
    m_exponentVector(dIdx) = value;
} // end function setExponent

void GaussianPrimitiveBasis::setExponent(const Eigen::VectorXi& expVec) {
    /* set exponent value */
    m_exponentVector = expVec;
} // end function setExponent

int GaussianPrimitiveBasis::exponentMax() const {
    /* return maximum exponent value */
    return m_exponentVector.maxCoeff();
} // end function exponentMax

int GaussianPrimitiveBasis::angularMomentum() const {
    /* return angular momentum (sum of exponents) */
    return m_exponentVector.sum();
} // end function angularMomentum
