#ifndef HERMITE_H
#define HERMITE_H

#include <vector>
#include <cmath>

static std::vector<int> HC0() {
    std::vector<int> coeffs = std::vector<int>{1};

    return coeffs;
}
static std::vector<int> HC1() {
    std::vector<int> coeffs = std::vector<int>{0,2};

    return coeffs;
}
static std::vector<int> HC2() {
    std::vector<int> coeffs = std::vector<int>{-2,0,4};

    return coeffs;
}
static std::vector<int> HC3() {
    std::vector<int> coeffs = std::vector<int>{0,-12,0,8};

    return coeffs;
}
static std::vector<int> HC4() {
    std::vector<int> coeffs = std::vector<int>{12,0,-48,0,16};

    return coeffs;
}
static std::vector<int> HC5() {
    std::vector<int> coeffs = std::vector<int>{0,120,0,-160,0,32};

    return coeffs;
}
static std::vector<int> HC6() {
    std::vector<int> coeffs = std::vector<int>{-120,0,720,0,-480,0,64};

    return coeffs;
}
static std::vector<int> HC7() {
    std::vector<int> coeffs = std::vector<int>{0,-1680,0,3360,0,-1344,0,128};

    return coeffs;
}
static std::vector<int> HC8() {
    std::vector<int> coeffs = std::vector<int>{1680,0,-13440,0,13440,0,-3584,0,256};

    return coeffs;
}
static std::vector<int> HC9() {
    std::vector<int> coeffs = std::vector<int>{0,30240,0,-80640,0,48384,0,-9216,0,512};

    return coeffs;
}
static std::vector<int> HC10() {
    std::vector<int> coeffs = std::vector<int>{-30240,0,302400,0,-403200,0,161280,0,-23040,0,1024};

    return coeffs;
}
static std::vector<int> HC11() {
    std::vector<int> coeffs = std::vector<int>{0,-665280,0,2217600,0,-1774080,0,506880,0,-56320,0,2048};

    return coeffs;
}
static std::vector<int> HC12() {
    std::vector<int> coeffs = std::vector<int>{665280,0,-7983360,0,13305600,0,-7096320,0,1520640,0,-135168,0,4096};

    return coeffs;
}
static std::vector<int> HC13() {
    std::vector<int> coeffs = std::vector<int>{0,17297280,0,-69189120,0,69189120,0,-26357760,0,4392960,0,-319488,0,8192};

    return coeffs;
}
static std::vector<int> HC(int n) {
   if (n > 13) {
       return std::vector<int>{-1};
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
       default: return std::vector<int>{0};
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
#pragma GCC diagnostic pop
template<typename T> static T H(T x, int n) {
   if (n > 12) {
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
       default: return 0;
   }
}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
 template<typename T> static T dH0(T x) {return 0;}
 template<typename T> static T dH1(T x) {return 2;}
 template<typename T> static T dH2(T x) {return 8*x;}
 template<typename T> static T dH3(T x) {return 24*pow(x, 2) - 12;}
 template<typename T> static T dH4(T x) {return 64*pow(x, 3) - 96*x;}
 template<typename T> static T dH5(T x) {return 160*pow(x, 4) - 480*pow(x, 2) + 120;}
 template<typename T> static T dH6(T x) {return 384*pow(x, 5) - 1920*pow(x, 3) + 1440*x;}
 template<typename T> static T dH7(T x) {return 896*pow(x, 6) - 6720*pow(x, 4) + 10080*pow(x, 2) - 1680;}
 template<typename T> static T dH8(T x) {return 2048*pow(x, 7) - 21504*pow(x, 5) + 53760*pow(x, 3) - 26880*x;}
 template<typename T> static T dH9(T x) {return 4608*pow(x, 8) - 64512*pow(x, 6) + 241920*pow(x, 4) - 241920*pow(x, 2) + 30240;}
 template<typename T> static T dH10(T x) {return 10240*pow(x, 9) - 184320*pow(x, 7) + 967680*pow(x, 5) - 1612800*pow(x, 3) + 604800*x;}
 template<typename T> static T dH11(T x) {return 22528*pow(x, 10) - 506880*pow(x, 8) + 3548160*pow(x, 6) - 8870400*pow(x, 4) + 6652800*pow(x, 2) - 665280;}
 template<typename T> static T dH12(T x) {return 49152*pow(x, 11) - 1351680*pow(x, 9) + 12165120*pow(x, 7) - 42577920*pow(x, 5) + 53222400*pow(x, 3) - 15966720*x;}
#pragma GCC diagnostic pop
template<typename T> static T dH(T x, int n) {
   if (n > 12) {
       return -1;
   }
   switch(n) {
       case 0: return dH0(x);
       case 1: return dH1(x);
       case 2: return dH2(x);
       case 3: return dH3(x);
       case 4: return dH4(x);
       case 5: return dH5(x);
       case 6: return dH6(x);
       case 7: return dH7(x);
       case 8: return dH8(x);
       case 9: return dH9(x);
       case 10: return dH10(x);
       case 11: return dH11(x);
       case 12: return dH12(x);
       default: return 0;
   }
}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
 template<typename T> static T ddH0(T x) {return 0;}
 template<typename T> static T ddH1(T x) {return 0;}
 template<typename T> static T ddH2(T x) {return 8;}
 template<typename T> static T ddH3(T x) {return 48*x;}
 template<typename T> static T ddH4(T x) {return 192*pow(x, 2) - 96;}
 template<typename T> static T ddH5(T x) {return 640*pow(x, 3) - 960*x;}
 template<typename T> static T ddH6(T x) {return 1920*pow(x, 4) - 5760*pow(x, 2) + 1440;}
 template<typename T> static T ddH7(T x) {return 5376*pow(x, 5) - 26880*pow(x, 3) + 20160*x;}
 template<typename T> static T ddH8(T x) {return 14336*pow(x, 6) - 107520*pow(x, 4) + 161280*pow(x, 2) - 26880;}
 template<typename T> static T ddH9(T x) {return 36864*pow(x, 7) - 387072*pow(x, 5) + 967680*pow(x, 3) - 483840*x;}
 template<typename T> static T ddH10(T x) {return 92160*pow(x, 8) - 1290240*pow(x, 6) + 4838400*pow(x, 4) - 4838400*pow(x, 2) + 604800;}
 template<typename T> static T ddH11(T x) {return 225280*pow(x, 9) - 4055040*pow(x, 7) + 21288960*pow(x, 5) - 35481600*pow(x, 3) + 13305600*x;}
 template<typename T> static T ddH12(T x) {return 540672*pow(x, 10) - 12165120*pow(x, 8) + 85155840*pow(x, 6) - 212889600*pow(x, 4) + 159667200*pow(x, 2) - 15966720;}
#pragma GCC diagnostic pop
template<typename T> static T ddH(T x, int n) {
   if (n > 12) {
       return -1;
   }
   switch(n) {
       case 0: return ddH0(x);
       case 1: return ddH1(x);
       case 2: return ddH2(x);
       case 3: return ddH3(x);
       case 4: return ddH4(x);
       case 5: return ddH5(x);
       case 6: return ddH6(x);
       case 7: return ddH7(x);
       case 8: return ddH8(x);
       case 9: return ddH9(x);
       case 10: return ddH10(x);
       case 11: return ddH11(x);
       case 12: return ddH12(x);
       default: return 0;
   }
}
#endif /* HERMITE_H */
