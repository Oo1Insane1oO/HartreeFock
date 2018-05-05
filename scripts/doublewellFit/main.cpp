#include "cartesian.h"
#include <iostream>
#include <iomanip>
#include <string>

int main(int argc, char *argv[]) {
    int cut = atoi(argv[1]);
    int dim = atoi(argv[2]);
    std::string fname = argv[3];

    Cartesian* cart = new Cartesian();
    cart->setup(cut, dim);
    cart->restructureStates();
    cart->writeToFile(fname);

    delete cart;

    return 0;
} // end main
