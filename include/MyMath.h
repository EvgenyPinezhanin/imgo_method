#ifndef _MY_MATH_H_
#define _MY_MATH_H_

#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <vector>

struct point {
    std::vector<double> x;

    point() : x() {};
    point(double _x) : x(_x) {};
    point(double _x, double _y) : x({ _x, _y }) {};
    point(std::vector<double> _x) : x(_x) {};
    point(std::vector<double> _x, double _y) : x(_x) { x.push_back(_y); };

    size_t getDimension() const { return x.size(); };
};

class mnk {
private:
    std::vector<std::vector<double>> A;
    std::vector<double> B;

public:
    mnk() : A(), B() {};
    mnk(const std::vector<std::vector<double>>& _A, const std::vector<double>& _B)
        : A(_A), B(_B) {};

    void setA(const std::vector<std::vector<double>>& _A) { A = _A; };
    void setB(const std::vector<double>& _B) { B = _B; };

    void solve(std::vector<double> &x) const;
};

double euclideanDistance(const std::vector<double> &firstValue, const std::vector<double> &secondValue);
double euclideanDistance(double firstValue, double secondValue);

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif // _MY_MATH_H_
