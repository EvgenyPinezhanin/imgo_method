#ifndef IMGO_H
#define IMGO_H

#include <vector>
#include <functional>

#include <opt_method.h>
#include <task.h>

using namespace std;

class imgo_method : public optimization_method_constrained {
private:
    function<double(double, int)> f;
    double r, d;

    trial_constr last_trial;
    int last_trial_pos;

    vector<vector<trial_constr>> I;
    vector<bool> calc_I;
    vector<double> mu;
    vector<double> z_star;

    trial_constr newTrial(double x) override;
    double newPoint(int t) override;
    double selectNewPoint(int &t) override;

public:
    imgo_method(function<double(double, int)> _f, int _m = 0, double _a = 0.0, double _b = 10.0, double _r = 2.0, double _d = 0.0, 
                double _eps = 0.0001, int _Nmax = 1000) : optimization_method_constrained(nullptr, 1, _m, vector<double>{_a}, 
                vector<double>{_b}, _eps, _Nmax), f(_f), r(_r), d(_d), last_trial(0.0, 0.0, 0), last_trial_pos(0), I((size_t)m + 1),
                calc_I((size_t)m + 1), mu((size_t)m + 1), z_star((size_t)m + 1) {}
    
    void setF(const function<double(double, int)> &_f) { f = _f; };
    void setA(double _a) { optimization_method::setA(vector<double>{_a}); };
    void setB(double _b) { optimization_method::setB(vector<double>{_b}); };
    void setAB(double _a, double _b) { optimization_method::setAB(vector<double>{_a}, vector<double>{_b}); };
    void setM(int _m);
    void setR(double _r) { r = _r; };
    void setD(double _d) { d = _d; };

    function<double(double, int)> getF() const { return f; };
    double getA() const { return A[0]; };
    double getB() const { return B[0]; };
    int getM() const { return m; };
    double getR() const { return r; };
    double getD() const { return d; };

    void getMu(vector<double> &mu_vec) const { mu_vec = mu; };

    void solve(int &count, double &x, Stop stop = Stop::ACCURACY);
    void solve(int &count, vector<double> &X, Stop stop = Stop::ACCURACY);
    
    bool solve_test(double x_opt, int &count, Stop stop = Stop::ACCURACY);
};

#endif // IMGO_H
