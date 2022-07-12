#ifndef MGGSA_H
#define MGGSA_H

#include <vector>

#include <opt_method.h>

using namespace std;

class mggsa_method : public optimization_method_constrained {
private:
    double r, d;
    int den, key;

    vector<vector<trial_constr>> I;
    vector<bool> calc_I;
    vector<double> mu;
    vector<double> z_star;

    int M;

    trial_constr newTrial(double x);
    double newPoint(int t);
    double selectNewPoint(int &t, trial_constr last_trial);

public:
    mggsa_method(double (*_f)(vector<double>, int), int _n = 2, int _m = 0, vector<double> _A = vector<double>(), 
                 vector<double> _B = vector<double>(), double _r = 2.0, double _d = 0.01, int _den = 10, 
                 int _key = 1, double _eps = 0.0001, int _Nmax = 1000)
                : optimization_method_constrained(_f, _n, _m, _A, _B, _eps, _Nmax), r(_r), d(_d), den(_den), key(_key), 
                  M(0), I(m + 1), calc_I(m + 1), mu(m + 1), z_star(m + 1) {}
    
    void setM(int _m);
    void setR(double _r) { r = _r; };
    void setD(double _d) { d = _d; };
    void setDen(int _den) { den = _den; };
    void setKey(int _key) { key = _key; };

    void getPoints(vector<vector<double>> &points_vec);

    void y(double x, vector<double> &X);
    
    void solve(int &count, vector<double> &X, Stop stop = ACCURACY);
    bool solve_test(vector<double> x_opt, int &count, Stop stop = ACCURACY);
};

#endif // MGGSA_H
