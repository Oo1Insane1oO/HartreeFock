#ifdef STYPEGAUSSIAN
#include <UnitTest++/UnitTest++.h>

#include <memory>

#include <Eigen/Dense>

#include "../src/hartreefocksolver.h"

#include <iostream>

SUITE(HFStype) {
    /* group for testing setup and energy of HFStype */

    class HFStypeFix {
        public:
            std::unique_ptr<StypeBasis> sb;
            double eps;

            HFStypeFix() {
                sb = std::make_unique<StypeBasis>(2);
                eps = 1e-14;
            } // end constructor
            ~HFStypeFix() {
            } // end deconstructor
    }; // end class HFStypeFix

    TEST_FIXTURE(HFStypeFix, checkScalingVec2D) {
        /* check for isotropic case */
        sb->~StypeBasis();
        new (&*sb) StypeBasis(2);
        Eigen::VectorXd scalingVec = Eigen::VectorXd::Zero(5);
        scalingVec << 0.25, 0.5, 1.0, 1.5, 2.0;
        sb->setPrimitives(scalingVec);

        for (unsigned int i = 0; i < sb->getSize(); ++i) {
            for (unsigned int d = 0; d < 2; ++d) {
                CHECK_EQUAL(scalingVec(i),
                        sb->getBasis()->getPrimitive(i)->scalingVec()(d));
            } // end ford
        } // end fori
    } // end TEST_FIXTURE checkScalingVec2D
    
    TEST_FIXTURE(HFStypeFix, checkScalingVec3D) {
        /* check for isotropic case */
        sb->~StypeBasis();
        new (&*sb) StypeBasis(3);
        Eigen::VectorXd scalingVec = Eigen::VectorXd::Zero(5);
        scalingVec << 0.25, 0.5, 1.0, 1.5, 2.0;
        sb->setPrimitives(scalingVec);

        for (unsigned int i = 0; i < sb->getSize(); ++i) {
            for (unsigned int d = 0; d < 3; ++d) {
                CHECK_EQUAL(scalingVec(i),
                        sb->getBasis()->getPrimitive(i)->scalingVec()(d));
            } // end ford
        } // end fori
    } // end TEST_FIXTURE checkScalingVec3D

    TEST_FIXTURE(HFStypeFix, checkScalingMat2D) {
        /* check for non-isotropic case */
        sb->~StypeBasis();
        new (&*sb) StypeBasis(2);
        Eigen::MatrixXd scalingMat = Eigen::MatrixXd::Random(5,2);
        sb->setPrimitives(scalingMat);

        for (unsigned int i = 0; i < sb->getSize(); ++i) {
            for (unsigned int d = 0; d < 2; ++d) {
                CHECK_EQUAL(scalingMat(i,d),
                        sb->getBasis()->getPrimitive(i)->scalingVec()(d));
            } // end ford
        } // end fori
    } // end TEST_FIXTURE checkScalingMat2D
    
    TEST_FIXTURE(HFStypeFix, checkScalingMat3D) {
        /* check for non-isotropic case */
        sb->~StypeBasis();
        new (&*sb) StypeBasis(3);
        Eigen::MatrixXd scalingMat = Eigen::MatrixXd::Random(5,3);
        sb->setPrimitives(scalingMat);

        for (unsigned int i = 0; i < sb->getSize(); ++i) {
            for (unsigned int d = 0; d < 3; ++d) {
                CHECK_EQUAL(scalingMat(i,d),
                        sb->getBasis()->getPrimitive(i)->scalingVec()(d));
            } // end ford
        } // end fori
    } // end TEST_FIXTURE checkScalingMat3D
} // end SUITE HFStype
#endif
