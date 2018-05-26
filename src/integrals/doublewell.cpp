#include "doublewell.h"
#include "../hermite/hermite.h"
#include "../methods.h"

#include <boost/math/special_functions/gamma.hpp>

DoubleWell::DoubleWell(const unsigned int dim, unsigned int cutOff, const
        unsigned int numParticles) :
    HartreeFockSolver<DoubleWell>(this, dim, 2*cutOff, numParticles),
    DWC(dim),
    GaussianIntegrals(dim, 2*DWC::C.rows(), numParticles) {
    /* construct */
    m_numBasis = cutOff;
} // end constructor

DoubleWell::~DoubleWell() {
} // end deconstructor

std::string DoubleWell::initializeParameters(double _R, unsigned int axis) {
    /* initialize GaussianIntegrals and grab well-separation parameter R
     * (axis=0 by default) */

    // default to shifting well in x-direction
    m_axis = axis;
    R = _R;
    RsqrdFactor = 1.0/8.0 * R*R;
    std::string message = GaussianIntegrals::initializeParameters(1.0);

    // precalculate two-body elements over HO-functions elements
    GaussianIntegrals::HartreeFockSolver<GaussianIntegrals>::assemble(1);
    
    DoubleWell::HartreeFockSolver<DoubleWell>::setReadFile(false);
    DoubleWell::HartreeFockSolver<DoubleWell>::setIndexWeigthing(false);

    return message;
} // end function initializeParameters

double DoubleWell::overlapElement(const unsigned int& i, const unsigned int& j)
{
    /* calculate and return the overlap integral element <i|j> */
    return ((i==j) ? 1.0 : 0.0);
} // end function overlapElement

double DoubleWell::oneBodyElement(const unsigned int& i, const unsigned int& j)
{
    /* calculate and return oneBodyElement <i|h|k> = <i|K|j> + <i|P|j>, where K
     * is the kinetic part and P is the potential part. For this case the K+P
     * part is the HO-potential part taken from GaussianIntegrals with the
     * added DW part */
    return ((i==j) ? DWC::epsDW(i) : 0.0);
} // end function oneBodyElement

double DoubleWell::coulombElement(const unsigned int& i, const unsigned int& j,
        const unsigned int& k, const unsigned int& l) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> */
    double res = 0.0;
    for (Eigen::SparseMatrix<double>::InnerIterator pIt(DWC::C, i); pIt; ++pIt)
    for (Eigen::SparseMatrix<double>::InnerIterator qIt(DWC::C, j); qIt; ++qIt)
    for (Eigen::SparseMatrix<double>::InnerIterator rIt(DWC::C, k); rIt; ++rIt)
    for (Eigen::SparseMatrix<double>::InnerIterator sIt(DWC::C, l); sIt; ++sIt)
    {
        res += pIt.value()*qIt.value()*rIt.value()*sIt.value() *
            GaussianIntegrals::HartreeFockSolver::getTwoBodyElement(pIt.row(),
                    qIt.row(), rIt.row(), sIt.row());
    } // end for pIt,qIt,rIt,sIt

    return res;
} // end function coulombElement
