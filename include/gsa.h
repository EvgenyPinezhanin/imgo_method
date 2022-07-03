#ifndef GSA_H
#define GSA_H

#include<vector>
#include<algorithm>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include<math.h>
#else
    #include<cmath>
#endif

using namespace std;

struct trial {
    double x, z;

    trial(double _x, double _z) : x(_x), z(_z) {};
};

class gsa_method {
private:
    double a, b;
    double eps;
    double r;
    double (*f)(double);

    vector<trial> trial_points;

    void addInSort(double x);
    double searchMinX();
    double selectNewPoint(size_t &t);

public:
    gsa_method(double (*_f)(double), double _a, double _b, double _eps, double _r = 2.0);
    
    void setFunc(double (*_f)(double));
    void setA(double _a);
    void setB(double _b);
    void setEps(double _eps);
    
    double solve(int &n);
};

#endif // GSA_H
