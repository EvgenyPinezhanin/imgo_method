#ifndef IMGO_H
#define IMGO_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <peano.h>

using namespace std;

void addInSort(vector<trial> &vec, trial tr);
double searchMinXTrial(vector<trial> &trials, int m);

class imgo_method : public optimization_method {
protected:
    double (*f_single)(double, int);
    double r, d;
    int den, key;

    int M;
    vector<vector<trial>> I;
    vector<bool> calc_I;
    vector<double> mu;
    vector<double> z_star;

    trial newTrial_single(double x);
    trial newTrial(double x);
    double deltaX(double x1, double x2);
    double newPoint(size_t t);
    double selectNewPoint(size_t &t, trial last_trial);

public:
    // IMGO method
    imgo_method(double (*_f)(double, int), int _m = 0, double _a = 0.0, double _b = 10.0, double _r = 2.0, double _d = 0.01, double _eps = 0.0001, int _Nmax = 1000)
        : optimization_method(nullptr, 1, _m, vector<double>{_a}, vector<double>{_b}, _eps, _Nmax), 
        f_single(_f), r(_r), d(_d), I(m + 1), calc_I(m + 1), mu(m + 1), z_star(m + 1) {}

    // MGGSA method
    imgo_method(double (*_f)(vector<double>, int), int _n, size_t _m = 0, vector<double> _A = vector<double>(), 
    vector<double> _B = vector<double>(), double _r = 2.0, double _d = 0.01, double _eps = 0.0001, int _Nmax = 1000, int _den = 10, int _key = 1)
        : optimization_method(_f, _n, _m, _A, _B, _eps, _Nmax), r(_r), d(_d), den(_den), key(_key), M(0), I(m + 1), calc_I(m + 1), mu(m + 1), z_star(m + 1) {}
    
    // IMGO method
    void setFunc(double (*_f)(double, int)) { f_single = _f; };
    // MGGSA method
    void setFunc(double (*_f)(vector<double>, int)) { optimization_method::setF(_f); };
    // IMGO method
    void setA(double _a) { optimization_method::setA(vector<double>{_a}); };
    void setB(double _b) { optimization_method::setB(vector<double>{_b}); };
    // MGGSA method
    void setA(const vector<double> &_A) { optimization_method::setA(_A); };
    void setB(const vector<double> &_B) { optimization_method::setB(_B); };
    // IMGO method
    void setAB(double _a, double _b) { optimization_method::setAB(vector<double>{_a}, vector<double>{_b}); };
    // MGGSA method
    void setAB(const vector<double> &_A, const vector<double> &_B) { optimization_method::setAB(_A, _B); };
    void setM(int _m);
    void setR(double _r) { r = _r; };
    void setD(double _d) { d = _d; };
    void setDen(int _den) { den = _den; };
    void setKey(int _key) { key = _key; };

    void getPoints(vector<vector<double>> &points_vec);

    // MGGSA method
    void y(double x, vector<double> &X);
    
    // IMGO method
    double solve(int &count, Stop stop = ACCURACY);
    // MGGSA method
    void solve(int &count, vector<double> &X, Stop stop = ACCURACY);
    // IMGO method
    bool solve_test(double x_opt, int &count, Stop stop = ACCURACY);
    // MGGSA method
    bool solve_test(vector<double> x_opt, int &count, Stop stop = ACCURACY);
};

#endif // IMGO_H
