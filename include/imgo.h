#ifndef IMGO_H
#define IMGO_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <peano.h>

using namespace std;

enum Stop {ACCURACY, NUMBER, ACCURNUMBER};

struct trial {
    double x, z;
    size_t nu;

    trial(double _x = 0.0, double _z = 0.0, int _nu = 0) : x(_x), z(_z), nu(_nu) {}
};

void addInSort(vector<trial> &vec, trial tr);
double searchMinXTrial(vector<trial> &trials, int m);

class optimization_method {
protected:
    double (*f)(vector<double>, int);
    int n, m;
    vector<double> A, B;
    double eps;
    int Nmax;

    vector<trial> trial_points;

    virtual trial newTrial(double x) = 0;
    virtual double newPoint(size_t t) = 0;
    virtual double selectNewPoint(size_t &t, trial last_trial) = 0;

public:
    optimization_method(double (*_f)(vector<double>, int), int _n, int _m, const vector<double> &_A, const vector<double> &_B, double _eps, int _Nmax) 
        : f(_f), n(_n), m(_m), A(_A), B(_B), eps(_eps), Nmax(_Nmax) {}

    void setF(double (*_f)(vector<double>, int)) { f = _f; };
    void setN(int _n) { n = _n; };
    void setM(int _m) { m = _m; };
    void setNmax(int _Nmax) { Nmax = _Nmax; };
    void setEps(double _eps) { eps = _eps; };
    void setA(const vector<double> &_A) { A = _A; };
    void setB(const vector<double> &_B) { B = _B; };
    void setAB(const vector<double> &_A, const vector<double> &_B) { A = _A; B = _B; };

    void getTrialPoints(vector<trial> &trial_vec) const { trial_vec = trial_points; };

    virtual void solve(int &n, vector<double> &X, Stop stop = ACCURACY) = 0;
};

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
    imgo_method(double (*_f)(double, int), int _m = 0, double _a = 0.0, double _b = 10.0, double _r = 2.0, double _d = 0.01, double _eps = 0.0001, int _Nmax = 1000)
        : optimization_method(nullptr, 1, _m, vector<double>{_a}, vector<double>{_b}, _eps, _Nmax), 
        f_single(_f), r(_r), d(_d), I(m + 1), calc_I(m + 1), mu(m + 1), z_star(m + 1) {}

    imgo_method(double (*_f)(vector<double>, int), int _n, size_t _m = 0, vector<double> _A = vector<double>(), 
    vector<double> _B = vector<double>(), double _r = 2.0, double _d = 0.01, double _eps = 0.0001, int _Nmax = 1000, int _den = 10, int _key = 1)
        : optimization_method(_f, _n, _m, _A, _B, _eps, _Nmax), r(_r), d(_d), den(_den), key(_key), M(0), I(m + 1), calc_I(m + 1), mu(m + 1), z_star(m + 1) {}
    
    void setFunc(double (*_f)(double, int)) { f_single = _f; };
    void setFunc(double (*_f)(vector<double>, int)) { optimization_method::setF(_f); };
    void setA(double _a) { optimization_method::setA(vector<double>{_a}); };
    void setB(double _b) { optimization_method::setB(vector<double>{_b}); };
    void setA(const vector<double> &_A) { optimization_method::setA(_A); };
    void setB(const vector<double> &_B) { optimization_method::setB(_B); };
    void setAB(double _a, double _b) { optimization_method::setAB(vector<double>{_a}, vector<double>{_b}); };
    void setAB(const vector<double> &_A, const vector<double> &_B) { optimization_method::setAB(_A, _B); };
    void setM(int _m);
    void setR(double _r) { r = _r; };
    void setD(double _d) { d = _d; };
    void setDen(int _den) { den = _den; };
    void setKey(int _key) { key = _key; };

    void getPoints(vector<vector<double>> &points_vec);

    void y(double x, vector<double> &X);
    
    double solve(int &count, Stop stop = ACCURACY);
    void solve(int &count, vector<double> &X, Stop stop = ACCURACY);
    bool solve_test(double x_opt, int &count, Stop stop = ACCURACY);
    bool solve_test(vector<double> x_opt, int &count, Stop stop = ACCURACY);
};

#endif // IMGO_H
