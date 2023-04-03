#ifndef MGGSA_H
#define MGGSA_H

#include <vector>
#include <functional>

#include <opt_method.h>

using namespace std;

enum class TypeSolve { SOLVE, RESOLVE };

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
    mggsa_method(function<double(vector<double>, int)> _f, int _n = 2, int _numberConstraints = 0,
                 vector<double> _A = vector<double>(), vector<double> _B = vector<double>(), double _r = 2.0,
                 double _d = 0.01, int _den = 10, int _key = 1, double _eps = 0.0001, int _maxIters = 1000,
                 int _maxEvals = 1000, int _incr = 0) : optimization_method_constrained(_f, _n, _numberConstraints,
                 _A, _B, _eps, _maxIters, _maxEvals), r(_r), d(_d), den(_den), key(_key), incr(_incr), last_trials(1),
                 last_trials_pos(1), M(0), I((size_t)numberConstraints + 1), h_nu(0), calc_I((size_t)numberConstraints + 1),
                 mu((size_t)numberConstraints + 1), z_star((size_t)numberConstraints + 1) {}
    
    void setNumberConstraints(int _numberConstraints);
    void setR(double _r) { r = _r; };
    void setD(double _d) { d = _d; };
    void setDen(int _den) { den = _den; };
    void setKey(int _key) { key = _key; };
    void setIncr(int _incr) { incr = _incr; };

    double getR() const { return r; };
    double getD() const { return d; };
    int getDen() const { return den; };
    int getKey() const { return key; };
    int getIncr() const { return incr; };

    void getPoints(vector<vector<double>> &points);
    void getLambda(vector<double> &lambdas) const { lambdas = mu; };
    
    void solve(int &countIters, int &countTrials, int &countEvals, vector<double> &X, TypeSolve type);
    void solve(int &countIters, int &countTrials, int &countEvals, vector<double> &X) override;

    bool solve_test(vector<double> X_opt, int &countIters, int &countTrials, int &countEvals);
};

#endif // MGGSA_H
