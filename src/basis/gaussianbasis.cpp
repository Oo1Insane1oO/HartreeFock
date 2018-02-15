#include "gaussianbasis.h"

GaussianBasis::GaussianBasis() : Cartesian::Cartesian() {
    /* default constructor */
} // end constructor

GaussianBasis::GaussianBasis(unsigned int cut, unsigned int dimension) {
    setup(cut, dimension);
} // end constructor

GaussianBasis::~GaussianBasis() {
} // end deconstructor

unsigned int GaussianBasis::getSize() {
    /* return number of states (basis functions) */
    return Cartesian::getStates().rows();
} // end function getSize

void GaussianBasis::setup(unsigned int cut, unsigned int dimension) {
    /* set number of dimensions */
    m_dim = dimension;
    Cartesian::setup(cut, m_dim);
    Cartesian::restructureStates();
} // end function setDim
