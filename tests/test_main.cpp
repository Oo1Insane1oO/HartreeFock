
#ifdef GAUSSHERMITE
    #include "test_gaussian.cpp"
#endif

#ifdef STYPEGAUSSIAN
    #include "test_stype.cpp"
#endif

int test_main() {
    return UnitTest::RunAllTests();
} // end function test_main
