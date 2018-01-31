#ifndef METHODS_H
#define METHODS_H

#include <Eigen/Dense>
#include <string>
#include <iostream>
#include <iomanip>
#include <boost/math/special_functions/factorials.hpp>

class Methods {
    private:
        Methods ();
        virtual ~Methods ();

    public:
        static void printProgressBar(std::string&, const float&, const int,
                const std::string& extra="");
        static std::string stringPos(const unsigned int&, const int&);
        static int divider(const unsigned int&, const unsigned int&, const
                int&);
        static void setDisplacement(int[], int[], const int[], int);
        
        template <typename T> static inline T refSum(const Eigen::Matrix<T*,
                Eigen::Dynamic, 1>& vec) {
            /* sum values of a vector of pointers(de-reference and sum) */
            T sum = 0;
            for (unsigned int i = 0; i < vec.size(); ++i) {
                sum += *(vec(i));
            } // end fori
            return sum;
        } // end function refSum
};

#endif /* METHODS_H */
