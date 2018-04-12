#ifndef HERMITE_H
#define HERMITE_H

#include <vector>
#include <cmath>

static std::vector<long int> HC0() {
    std::vector<long int> coeffs = std::vector<long int>{1};

    return coeffs;
}
static std::vector<long int> HC1() {
    std::vector<long int> coeffs = std::vector<long int>{0,2};

    return coeffs;
}
static std::vector<long int> HC2() {
    std::vector<long int> coeffs = std::vector<long int>{-2,0,4};

    return coeffs;
}
static std::vector<long int> HC3() {
    std::vector<long int> coeffs = std::vector<long int>{0,-12,0,8};

    return coeffs;
}
static std::vector<long int> HC4() {
    std::vector<long int> coeffs = std::vector<long int>{12,0,-48,0,16};

    return coeffs;
}
static std::vector<long int> HC5() {
    std::vector<long int> coeffs = std::vector<long int>{0,120,0,-160,0,32};

    return coeffs;
}
static std::vector<long int> HC6() {
    std::vector<long int> coeffs = std::vector<long int>{-120,0,720,0,-480,0,64};

    return coeffs;
}
static std::vector<long int> HC7() {
    std::vector<long int> coeffs = std::vector<long int>{0,-1680,0,3360,0,-1344,0,128};

    return coeffs;
}
static std::vector<long int> HC8() {
    std::vector<long int> coeffs = std::vector<long int>{1680,0,-13440,0,13440,0,-3584,0,256};

    return coeffs;
}
static std::vector<long int> HC9() {
    std::vector<long int> coeffs = std::vector<long int>{0,30240,0,-80640,0,48384,0,-9216,0,512};

    return coeffs;
}
static std::vector<long int> HC10() {
    std::vector<long int> coeffs = std::vector<long int>{-30240,0,302400,0,-403200,0,161280,0,-23040,0,1024};

    return coeffs;
}
static std::vector<long int> HC11() {
    std::vector<long int> coeffs = std::vector<long int>{0,-665280,0,2217600,0,-1774080,0,506880,0,-56320,0,2048};

    return coeffs;
}
static std::vector<long int> HC12() {
    std::vector<long int> coeffs = std::vector<long int>{665280,0,-7983360,0,13305600,0,-7096320,0,1520640,0,-135168,0,4096};

    return coeffs;
}
static std::vector<long int> HC13() {
    std::vector<long int> coeffs = std::vector<long int>{0,17297280,0,-69189120,0,69189120,0,-26357760,0,4392960,0,-319488,0,8192};

    return coeffs;
}
static std::vector<long int> HC14() {
    std::vector<long int> coeffs = std::vector<long int>{-17297280,0,242161920,0,-484323840,0,322882560,0,-92252160,0,12300288,0,-745472,0,16384};

    return coeffs;
}
static std::vector<long int> HC15() {
    std::vector<long int> coeffs = std::vector<long int>{0,-518918400,0,2421619200,0,-2905943040,0,1383782400,0,-307507200,0,33546240,0,-1720320,0,32768};

    return coeffs;
}
static std::vector<long int> HC16() {
    std::vector<long int> coeffs = std::vector<long int>{518918400,0,-8302694400,0,19372953600,0,-15498362880,0,5535129600,0,-984023040,0,89456640,0,-3932160,0,65536};

    return coeffs;
}
static std::vector<long int> HC17() {
    std::vector<long int> coeffs = std::vector<long int>{0,17643225600,0,-94097203200,0,131736084480,0,-75277762560,0,20910489600,0,-3041525760,0,233963520,0,-8912896,0,131072};

    return coeffs;
}
static std::vector<long int> HC(int n) {
   if (n > 17) {
       return std::vector<long int>{-1};
   }
   switch(n) {
       case 0: return HC0();
       case 1: return HC1();
       case 2: return HC2();
       case 3: return HC3();
       case 4: return HC4();
       case 5: return HC5();
       case 6: return HC6();
       case 7: return HC7();
       case 8: return HC8();
       case 9: return HC9();
       case 10: return HC10();
       case 11: return HC11();
       case 12: return HC12();
       case 13: return HC13();
       case 14: return HC14();
       case 15: return HC15();
       case 16: return HC16();
       case 17: return HC17();
       default: return std::vector<long int>{0};
   }
}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
 template<typename T> static T H0(T x) {return 1;}
 template<typename T> static T H1(T x) {return 2*x;}
 template<typename T> static T H2(T x) {return 4*pow(x, 2) - 2;}
 template<typename T> static T H3(T x) {return 8*pow(x, 3) - 12*x;}
 template<typename T> static T H4(T x) {return 16*pow(x, 4) - 48*pow(x, 2) + 12;}
 template<typename T> static T H5(T x) {return 32*pow(x, 5) - 160*pow(x, 3) + 120*x;}
 template<typename T> static T H6(T x) {return 64*pow(x, 6) - 480*pow(x, 4) + 720*pow(x, 2) - 120;}
 template<typename T> static T H7(T x) {return 128*pow(x, 7) - 1344*pow(x, 5) + 3360*pow(x, 3) - 1680*x;}
 template<typename T> static T H8(T x) {return 256*pow(x, 8) - 3584*pow(x, 6) + 13440*pow(x, 4) - 13440*pow(x, 2) + 1680;}
 template<typename T> static T H9(T x) {return 512*pow(x, 9) - 9216*pow(x, 7) + 48384*pow(x, 5) - 80640*pow(x, 3) + 30240*x;}
 template<typename T> static T H10(T x) {return 1024*pow(x, 10) - 23040*pow(x, 8) + 161280*pow(x, 6) - 403200*pow(x, 4) + 302400*pow(x, 2) - 30240;}
 template<typename T> static T H11(T x) {return 2048*pow(x, 11) - 56320*pow(x, 9) + 506880*pow(x, 7) - 1774080*pow(x, 5) + 2217600*pow(x, 3) - 665280*x;}
 template<typename T> static T H12(T x) {return 4096*pow(x, 12) - 135168*pow(x, 10) + 1520640*pow(x, 8) - 7096320*pow(x, 6) + 13305600*pow(x, 4) - 7983360*pow(x, 2) + 665280;}
 template<typename T> static T H13(T x) {return 8192*pow(x, 13) - 319488*pow(x, 11) + 4392960*pow(x, 9) - 26357760*pow(x, 7) + 69189120*pow(x, 5) - 69189120*pow(x, 3) + 17297280*x;}
 template<typename T> static T H14(T x) {return 16384*pow(x, 14) - 745472*pow(x, 12) + 12300288*pow(x, 10) - 92252160*pow(x, 8) + 322882560*pow(x, 6) - 484323840*pow(x, 4) + 242161920*pow(x, 2) - 17297280;}
 template<typename T> static T H15(T x) {return 32768*pow(x, 15) - 1720320*pow(x, 13) + 33546240*pow(x, 11) - 307507200*pow(x, 9) + 1383782400*pow(x, 7) - 2905943040*pow(x, 5) + 2421619200*pow(x, 3) - 518918400*x;}
 template<typename T> static T H16(T x) {return 65536*pow(x, 16) - 3932160*pow(x, 14) + 89456640*pow(x, 12) - 984023040*pow(x, 10) + 5535129600*pow(x, 8) - 15498362880*pow(x, 6) + 19372953600*pow(x, 4) - 8302694400*pow(x, 2) + 518918400;}
#pragma GCC diagnostic pop
template<typename T> static T H(T x, int n) {
   if (n > 16) {
       return -1;
   }
   switch(n) {
       case 0: return H0(x);
       case 1: return H1(x);
       case 2: return H2(x);
       case 3: return H3(x);
       case 4: return H4(x);
       case 5: return H5(x);
       case 6: return H6(x);
       case 7: return H7(x);
       case 8: return H8(x);
       case 9: return H9(x);
       case 10: return H10(x);
       case 11: return H11(x);
       case 12: return H12(x);
       case 13: return H13(x);
       case 14: return H14(x);
       case 15: return H15(x);
       case 16: return H16(x);
       default: return 0;
   }
}
#endif /* HERMITE_H */
