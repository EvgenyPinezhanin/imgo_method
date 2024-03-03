#ifndef MY_MATH_H_
#define MY_MATH_H_

#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <vector>

struct point {
    std::vector<double> x;

    point(double _x) : x(_x) {};
    point(double _x, double _y) : x({ _x, _y }) {};
    point(std::vector<double> _x) : x(_x) {};
    point(std::vector<double> _x, double _y) : x(_x) { x.push_back(_y); };

    size_t getDimension() const { return x.size(); };
};

class mnk {
private:
    std::vector<std::vector<double>> A;
    std::vector<double> b;

public:
    mnk() : A(), b() {};
    mnk(const std::vector<std::vector<double>>& _A, const std::vector<double>& _b)
        : A(_A), b(_b) {};

    void set_A(const std::vector<std::vector<double>>& _A) { A = _A; };
    void set_b(const std::vector<double>& _b) { b = _b; };

    void solve(std::vector<double> &x);
};

double euclideanDistance(const std::vector<double> &val1, const std::vector<double> &val2);

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif // MY_MATH_H_
