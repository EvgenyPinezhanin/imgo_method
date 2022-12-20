#ifndef OPT_METHOD_H
#define OPT_METHOD_H

#include <vector>

using namespace std;

enum class Stop { ACCURACY, NUMBER, ACCURNUMBER };

class optimization_method {
protected:
    int n; // dimension
    vector<double> A, B; // area of search
    double eps; // accuracy
    int Nmax; // maximum number of iterations

public:
    optimization_method(int _n, const vector<double> &_A, const vector<double> &_B, double _eps, int _Nmax) 
        : n(_n), A(_A), B(_B), eps(_eps), Nmax(_Nmax) {}

    void setN(int _n) { n = _n; };
    void setNmax(int _Nmax) { Nmax = _Nmax; };
    void setEps(double _eps) { eps = _eps; };
    void setA(const vector<double> &_A) { A = _A; };
    void setB(const vector<double> &_B) { B = _B; };
    void setAB(const vector<double> &_A, const vector<double> &_B) { A = _A; B = _B; };

    int getN() { return n; };
    int getNmax() { return Nmax; };
    double getEps() { return eps; };
    vector<double> getA() { return A; };
    vector<double> getB() { return B; };

    virtual void solve(int &count, vector<double> &X, Stop stop = Stop::ACCURACY) = 0;
};

struct trial {
    double x, z;

    trial(double _x = 0.0, double _z = 0.0) : x(_x), z(_z) {};
};

class optimization_method_non_constrained : public optimization_method {
protected:
    double (*f)(vector<double>); // target function

    vector<trial> trial_points;

    virtual trial newTrial(double x) = 0;
    virtual double newPoint(int t) = 0;
    virtual double selectNewPoint(int &t) = 0;

public:
    optimization_method_non_constrained(double (*_f)(vector<double>), int _n, const vector<double> &_A, const vector<double> &_B, double _eps, int _Nmax) 
        : optimization_method(_n, _A, _B, _eps, _Nmax), f(_f) {}

    void setF(double (*_f)(vector<double>)) { f = _f; };

    void getTrialPoints(vector<trial> &trial_vec) const { trial_vec = trial_points; };
};

struct trial_constr {
    double x, z;
    int nu;

    trial_constr(double _x = 0.0, double _z = 0.0, int _nu = 0) : x(_x), z(_z), nu(_nu) {}
};

class optimization_method_constrained : public optimization_method {
protected:
    double (*f)(vector<double>, int); // target function
    int m; // number of constraints

    vector<trial_constr> trial_points;

    virtual trial_constr newTrial(double x) = 0;
    virtual double newPoint(int t) = 0;
    virtual double selectNewPoint(int &t) = 0;

public:
    optimization_method_constrained(double (*_f)(vector<double>, int), int _n, int _m, const vector<double> &_A, const vector<double> &_B, double _eps, int _Nmax) 
        : optimization_method(_n, _A, _B, _eps, _Nmax), f(_f), m(_m) {}

    void setF(double (*_f)(vector<double>, int)) { f = _f; };
    void setM(int _m) { m = _m; };

    void getTrialPoints(vector<trial_constr> &trial_vec) const { trial_vec = trial_points; };
};

#endif // OPT_METHOD_H
