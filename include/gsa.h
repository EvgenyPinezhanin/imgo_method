#ifndef GSA_H
#define GSA_H

#include <vector>

#include <opt_method.h>

using namespace std;

class gsa_method : public optimization_method_non_constrained {
private:
    double r;
    double (*f)(double);

    double m; // parameter m

    trial last_trial;

    trial newTrial(double x) override;
    double newPoint(int t) override;
    double selectNewPoint(int &t) override;

public:
    gsa_method(double (*_f)(double), double _a = 0.0, double _b = 10.0, double _r = 2.0, double _eps = 0.001, int _Nmax = 1000)
        : optimization_method_non_constrained(nullptr, 1, vector<double>{_a}, vector<double>{_b}, _eps, _Nmax), 
          f(_f), r(_r), m(0), last_trial(0.0, 0.0) {}
    
    void setF(double (*_f)(double)) { f = _f; };
    void setA(double _a) { optimization_method::setA(vector<double>{_a}); };
    void setB(double _b) { optimization_method::setB(vector<double>{_b}); };
    void setAB(double _a, double _b) { optimization_method::setAB(vector<double>{_a}, vector<double>{_b}); };
    void setR(double _r) { r = _r; };
    
    void solve(int &count, double &x, Stop stop = Stop::ACCURACY);
    void solve(int &count, vector<double> &X, Stop stop = Stop::ACCURACY);
};

#endif // GSA_H
