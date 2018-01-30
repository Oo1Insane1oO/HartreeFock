#include "gaussianprimitivebasis.h"

GaussianPrimitiveBasis::GaussianPrimitiveBasis() {
    /* default constructor */
} // end constructor

GaussianPrimitiveBasis::GaussianPrimitiveBasis(double weight, Eigen::VectorXi
        expVec) {
    m_weight = weight;
    m_exponentVector = expVec;
} // end constructor

GaussianPrimitiveBasis::~GaussianPrimitiveBasis() {
} // end deconstructor

int GaussianPrimitiveBasis::dExponent(const unsigned int& dIdx) const {
    /* return exponent for dimension dIdx */
    return m_exponentVector(dIdx);
} // end function dExponent

const Eigen::VectorXi &GaussianPrimitiveBasis::expVec() const {
    /* return exponent for dimension dIdx */
    return m_exponentVector;
} // end function dExponent

void GaussianPrimitiveBasis::setWeight(double weight) {
    /* set m_weight to weight */
    m_weight = weight;
} // end function setWeight

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
