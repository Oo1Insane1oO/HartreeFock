#include "stypebasis.h"

StypeBasis::StypeBasis(unsigned int dim) : GaussianContractedBasis() {
    m_dim = dim;
} // end constructor

StypeBasis::~StypeBasis() {
} // end deconstructor

unsigned int StypeBasis::getSize() {
    return GaussianContractedBasis::getNumPrimitives();
} // end function getSize

GaussianContractedBasis* StypeBasis::getBasis() {
    /* return contracted basis set */
    return dynamic_cast<GaussianContractedBasis*>(this);
} // end function getPrimitive

void StypeBasis::setPrimitives(const Eigen::VectorXd& scalingVector) {
    /* create and set primitive functions in contracted basis (isotropic) */
    for (unsigned int i = 0; i < scalingVector.size(); ++i) {
        GaussianContractedBasis::addPrimitive(Eigen::VectorXd::Constant(m_dim,
                    scalingVector(i)), Eigen::VectorXi::Zero(m_dim));
    } // end fori
} // end function setPrimitives

void StypeBasis::setPrimitives(const Eigen::MatrixXd& scalingMatrix) {
    /* create and set primitive functions in contracted basis (non-isotropic)
     * */
    for (unsigned int i = 0; i < scalingMatrix.rows(); ++i) {
        GaussianContractedBasis::addPrimitive(scalingMatrix.row(i),
                Eigen::VectorXi::Zero(m_dim));
    } // end fori
} // end function setPrimitives
