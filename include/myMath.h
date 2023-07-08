#ifndef MY_MATH_H_
#define MY_MATH_H_

#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

inline double euclideanDistance(vector<double> val1, vector<double> val2) {
    double res = 0.0;
    size_t size = val1.size();
    for (int i = 0; i < size; i++) {
        res += (val1[i] - val2[i]) * (val1[i] - val2[i]);
    }
    return sqrt(res);
}

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif // MY_MATH_H_
