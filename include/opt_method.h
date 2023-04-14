#ifndef OPT_METHOD_H
#define OPT_METHOD_H

#include <vector>
#include <functional>

using namespace std;

class optimization_method {
protected:
    int n; // dimension
    vector<double> A, B; // area of search
    double eps; // accuracy
    int maxIters; // maximum number of iterations
    int maxEvals; // maximum number of calls to the target function
    int countEvals; // number of calls to the target function

public:
    optimization_method(int _n, const vector<double> &_A, const vector<double> &_B, double _eps, int _maxIters, int _maxEvals) 
                        : n(_n), A(_A), B(_B), eps(_eps), maxIters(_maxIters), maxEvals(_maxEvals), countEvals(0) {}

    void setN(int _n) { n = _n; };
    void setMaxIters(int _maxIters) { maxIters = _maxIters; };
    void setMaxEvals(int _maxEvals) { maxEvals = _maxEvals; };
    void setEps(double _eps) { eps = _eps; };
    void setA(const vector<double> &_A) { A = _A; };
    void setB(const vector<double> &_B) { B = _B; };
    void setAB(const vector<double> &_A, const vector<double> &_B) { A = _A; B = _B; };

    int getN() { return n; };
    int getMaxIters() { return maxIters; };
    int getMaxEvals() { return maxEvals; };
    double getEps() { return eps; };
    vector<double> getA() { return A; };
    vector<double> getB() { return B; };

    virtual void solve(int &countIters, int &countEvals, vector<double> &X) = 0;
};

struct trial {
    double x, z;

    trial(double _x = 0.0, double _z = 0.0) : x(_x), z(_z) {};
};

class optimization_method_non_constrained : public optimization_method {
protected:
    function<double(vector<double>)> f; // target function

    vector<trial> trial_points;

    virtual trial newTrial(double x) = 0;
    virtual double newPoint(int t) = 0;
    virtual double selectNewPoint(int &t) = 0;

public:
    optimization_method_non_constrained(function<double(vector<double>)> _f, int _n, const vector<double> &_A,
                                        const vector<double> &_B, double _eps, int _max_iters, int _max_evals)
                                        : optimization_method(_n, _A, _B, _eps, _max_iters, _max_evals), f(_f) {}

    void setF(function<double(vector<double>)> _f) { f = _f; };

    function<double(vector<double>)> getF() const { return f; };

    void getTrialPoints(vector<trial> &trial_points) const { trial_points = this->trial_points; };
    int getNumberTrialPoints() const { return trial_points.size(); };
};

struct trial_constr {
    double x, z;
    int nu;

    trial_constr(double _x = 0.0, double _z = 0.0, int _nu = 0) : x(_x), z(_z), nu(_nu) {}
};

class optimization_method_constrained : public optimization_method {
protected:
    function<double(vector<double>, int)> f; // target function
    int numberConstraints; // number of constraints

    vector<trial_constr> trial_points;

    virtual trial_constr newTrial(double x) = 0;
    virtual double newPoint(int t) = 0;
    virtual double selectNewPoint(int &t) = 0;

public:
    optimization_method_constrained(function<double(vector<double>, int)> _f, int _n, int _numberConstarints,
                                    const vector<double> &_A, const vector<double> &_B, double _eps, int _max_iters,
                                    int _max_evals) : optimization_method(_n, _A, _B, _eps, _max_iters, _max_evals),
                                    f(_f), numberConstraints(_numberConstarints) {}

    void setF(function<double(vector<double>, int)> _f) { f = _f; };
    void setNumberConstraints(int _numberConstraints) { numberConstraints = _numberConstraints; };

    function<double(vector<double>, int)> getF() const { return f; };
    int getNumberConstraints() const { return numberConstraints; };

    void getTrialPoints(vector<trial_constr> &trial_points) const { trial_points = this->trial_points; };
    int getNumberTrialPoints() const { return trial_points.size(); };
};

#endif // OPT_METHOD_H
