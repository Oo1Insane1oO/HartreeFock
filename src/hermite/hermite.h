#ifndef HERMITE_H
#define HERMITE_H
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
