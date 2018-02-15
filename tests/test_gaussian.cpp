#ifdef GAUSSHERMITE
#include <UnitTest++/UnitTest++.h>

#include <memory>

#include <Eigen/Dense>

#include "../src/hartreefocksolver.h"

SUITE(HFGaussUnperturbed) {
    /* group for testing that energies are correct for non-interacting case in
     * harmonic osciallator potential */

    class HFGaussUnperturbedFix {
        public:
            std::unique_ptr<HartreeFockSolver> HFS;
            double eps;

            HFGaussUnperturbedFix() {
                HFS = std::make_unique<HartreeFockSolver>(2,2,2);
                eps = 1e-14;
            } // end constructor
            ~HFGaussUnperturbedFix() {
            } // end deconstructor
    }; // end class HFGaussUnperturbedFix

    TEST_FIXTURE(HFGaussUnperturbedFix, check2DN2) {
        /* check energies in 2D for N=2 */
        HFS->~HartreeFockSolver();
        new (&*HFS) HartreeFockSolver(2,6,2);

        HFS->getIntegralObj()->initializeParameters(1.0);
        CHECK_CLOSE(2, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.5);
        CHECK_CLOSE(1, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.25);
        CHECK_CLOSE(0.5, HFS->iterate(100, 1e-8), 1e-14);
    } // end TEST_FIXTURE check2Dw1N2
    
    TEST_FIXTURE(HFGaussUnperturbedFix, check2DN6) {
        /* check energies in 2D for N=6 */
        HFS->~HartreeFockSolver();
        new (&*HFS) HartreeFockSolver(2,12,6);

        HFS->getIntegralObj()->initializeParameters(1.0);
        CHECK_CLOSE(10, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.5);
        CHECK_CLOSE(5, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.25);
        CHECK_CLOSE(2.5, HFS->iterate(100, 1e-8), 1e-14);
    } // end TEST_FIXTURE check2Dw1N6
    
    TEST_FIXTURE(HFGaussUnperturbedFix, check2DN12) {
        /* check energies in 2D for N=12 */
        HFS->~HartreeFockSolver();
        new (&*HFS) HartreeFockSolver(2,30,12);

        HFS->getIntegralObj()->initializeParameters(1.0);
        CHECK_CLOSE(28, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.5);
        CHECK_CLOSE(14, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.25);
        CHECK_CLOSE(7, HFS->iterate(100, 1e-8), 1e-14);
    } // end TEST_FIXTURE check2Dw1N12
    
    TEST_FIXTURE(HFGaussUnperturbedFix, check2DN20) {
        /* check energies in 2D for N=12 */
        HFS->~HartreeFockSolver();
        new (&*HFS) HartreeFockSolver(2,42,20);

        HFS->getIntegralObj()->initializeParameters(1.0);
        CHECK_CLOSE(60, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.5);
        CHECK_CLOSE(30, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.25);
        CHECK_CLOSE(15, HFS->iterate(100, 1e-8), 1e-14);
    } // end TEST_FIXTURE check2Dw1N20
    
    TEST_FIXTURE(HFGaussUnperturbedFix, check3DN2) {
        /* check energies in 3D for N=2 */
        HFS->~HartreeFockSolver();
        new (&*HFS) HartreeFockSolver(3,8,2);

        HFS->getIntegralObj()->initializeParameters(1.0);
        CHECK_CLOSE(3, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.5);
        CHECK_CLOSE(1.5, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.25);
        CHECK_CLOSE(0.75, HFS->iterate(100, 1e-8), 1e-14);
    } // end TEST_FIXTURE check3Dw1N2
    
    TEST_FIXTURE(HFGaussUnperturbedFix, check3DN8) {
        /* check energies in 3D for N=8 */
        HFS->~HartreeFockSolver();
        new (&*HFS) HartreeFockSolver(3,20,8);

        HFS->getIntegralObj()->initializeParameters(1.0);
        CHECK_CLOSE(18, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.5);
        CHECK_CLOSE(9, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.25);
        CHECK_CLOSE(4.5, HFS->iterate(100, 1e-8), 1e-14);
    } // end TEST_FIXTURE check3Dw1N8
    
    TEST_FIXTURE(HFGaussUnperturbedFix, check3DN20) {
        /* check energies in 3D for N=20 */
        HFS->~HartreeFockSolver();
        new (&*HFS) HartreeFockSolver(3,30,20);

        HFS->getIntegralObj()->initializeParameters(1.0);
        CHECK_CLOSE(60, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.5);
        CHECK_CLOSE(30, HFS->iterate(100, 1e-8), 1e-14);
        
        HFS->getIntegralObj()->initializeParameters(0.25);
        CHECK_CLOSE(15, HFS->iterate(100, 1e-8), 1e-14);
    } // end TEST_FIXTURE check3Dw1N20
} // end SUITE HFGauss
#endif
