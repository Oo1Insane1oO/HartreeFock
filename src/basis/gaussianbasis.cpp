#include "gaussianbasis.h"

GaussianBasis::GaussianBasis() : Cartesian::Cartesian(),
    GaussianContractedBasis::GaussianContractedBasis() {
        /* default constructor */
} // end constructor

GaussianBasis::GaussianBasis(unsigned int cut, unsigned int dimension) :
    GaussianContractedBasis::GaussianContractedBasis() {
    GaussianBasis::setup(cut, dimension);
} // end constructor

GaussianBasis::~GaussianBasis() {
} // end deconstructor

void GaussianBasis::setup(unsigned int cut, unsigned int dimension, double
        scaling) {
    /* set number of dimensions */
    m_dim = dimension;
    Cartesian::setup(cut, m_dim);
    Cartesian::restructureStates();
    setPrimitives(scaling);
} // end function setDim

void GaussianBasis::setPrimitives(const double scaling) {
    /* create and set primitive functions in contracted basis */
    for (unsigned int i = 0; i < Cartesian::getStates().rows(); ++i) {
        Eigen::VectorXi exponents = Eigen::VectorXi::Zero(m_dim);
        for (unsigned int d = 0; d < m_dim; ++d) {
            exponents(d) = *(Cartesian::getStates()(i,d));
        } // end ford
        GaussianContractedBasis::addPrimitive(scaling, exponents);
    } // end fori
} // end function setPrimitives

const Eigen::Ref<const Eigen::VectorXi> GaussianBasis::getExpVec(const unsigned
        int& j) const { 
    return GaussianContractedBasis::getPrimitive(j)->expVec();
} // end function getExpVec
