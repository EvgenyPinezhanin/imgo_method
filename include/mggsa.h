#ifndef MGGSA_H
#define MGGSA_H

#include <vector>
#include <functional>

#include <opt_method.h>

using namespace std;

class mggsa_method : public optimization_method_constrained {
private:
    double r, d;
    int den, key;
    int incr;

    vector<trial_constr> last_trials;
    vector<int> last_trials_pos;

    vector<vector<trial_constr>> I;
    vector<double> h_nu;
    vector<bool> calc_I;
    vector<double> mu;
    vector<double> z_star;

    int M;

    trial_constr newTrial(double x) override;
    double newPoint(int t) override;
    double selectNewPoint(int &t) override;

    double calc_h(double x_k_1);
    bool check_density(double h_j);
    void y(double x, vector<double> &X);
    void x(const vector<double> &P, vector<double> &X);

public:
    mggsa_method(function<double(vector<double>, int)> _f, int _n = 2, int _m = 0, vector<double> _A = vector<double>(), 
                 vector<double> _B = vector<double>(), double _r = 2.0, double _d = 0.01, int _den = 10, 
                 int _key = 1, double _eps = 0.0001, int _Nmax = 1000, int _incr = 0)
                : optimization_method_constrained(_f, _n, _m, _A, _B, _eps, _Nmax), r(_r), d(_d), den(_den), key(_key), incr(_incr), 
                  last_trials(1), last_trials_pos(1), M(0), I((size_t)m + 1), h_nu(0), calc_I((size_t)m + 1), 
                  mu((size_t)m + 1), z_star((size_t)m + 1) {}
    
    void setM(int _m);
    void setR(double _r) { r = _r; };
    void setD(double _d) { d = _d; };
    void setDen(int _den) { den = _den; };
    void setKey(int _key) { key = _key; };
    void setIncr(int _incr) { incr = _incr; };

    int getM() const { return m; };
    double getR() const { return r; };
    double getD() const { return d; };

    void getPoints(vector<vector<double>> &points_vec);
    int getCountPoints() const { return (int)trial_points.size() - 2; };
    void getMu(vector<double> &mu_vec) const { mu_vec = mu; };
    
    void solve(int &count, vector<double> &X, Stop stop = Stop::ACCURACY);
    bool solve_test(vector<double> X_opt, int &count, Stop stop = Stop::ACCURACY);
};

#endif // MGGSA_H
