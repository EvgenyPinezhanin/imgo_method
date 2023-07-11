#ifndef MY_MATH_H_
#define MY_MATH_H_

#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <vector>

using std::vector;

double euclideanDistance(const vector<double> &val1, const vector<double> &val2);

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif // MY_MATH_H_
