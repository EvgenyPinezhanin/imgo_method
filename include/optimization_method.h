#ifndef OPTIMIZATION_METHOD_H_
#define OPTIMIZATION_METHOD_H_

#include <vector>
#include <functional>

using namespace std;

class OptimizationMethod {
protected:
    int n; // dimension
    vector<double> A, B; // area of search
    double eps; // accuracy
    int maxTrials; // maximum number of trials
    int maxFevals; // maximum number of calls to the target function
    int numberFevals; // number of calls to the target function

public:
    OptimizationMethod(int _n, const vector<double> &_A, const vector<double> &_B, double _eps, int _maxTrials,
                       int _maxFevals) : n(_n), A(_A), B(_B), eps(_eps), maxTrials(_maxTrials), maxFevals(_maxFevals),
                       numberFevals(0) {}

    void setN(int _n) { n = _n; };
    void setMaxTrials(int _maxTrials) { maxTrials = _maxTrials; };
    void setMaxFevals(int _maxFevals) { maxFevals = _maxFevals; };
    void setEps(double _eps) { eps = _eps; };
    void setA(const vector<double> &_A) { A = _A; };
    void setB(const vector<double> &_B) { B = _B; };
    void setAB(const vector<double> &_A, const vector<double> &_B) { A = _A; B = _B; };

    int getN() const { return n; };
    int getMaxTrials() const { return maxTrials; };
    int getMaxFevals() const { return maxFevals; };
    double getEps() const { return eps; };
    void getA(vector<double>& _A) const { _A = A; };
    void getB(vector<double>& _B) const { _B = B; };
    void getAB(vector<double>& _A, vector<double>& _B) const { _A = A; _B = B; };

    virtual void solve(int &numberTrials, int &numberFevals, vector<double> &X) = 0;
};

struct Trial {
    double x, z;

    Trial(double _x = 0.0, double _z = 0.0) : x(_x), z(_z) {};
};

class OptimizationMethodNoConstrained : public OptimizationMethod {
protected:
    function<double(vector<double>)> f; // target function

    vector<Trial> trialPoints;

    virtual Trial newTrial(double x) = 0;
    virtual double newPoint(int t) = 0;
    virtual double selectNewPoint(int &t) = 0;

public:
    OptimizationMethodNoConstrained(function<double(vector<double>)> _f, int _n, const vector<double> &_A,
                                     const vector<double> &_B, double _eps, int _maxTrials, int _maxFevals)
                                     : OptimizationMethod(_n, _A, _B, _eps, _maxTrials, _maxFevals), f(_f) {}

    void setF(function<double(vector<double>)> _f) { f = _f; };

    function<double(vector<double>)> getF() const { return f; };

    void getTrialPoints(vector<Trial> &trialPoints) const { trialPoints = this->trialPoints; };
    int getNumberTrialPoints() const { return trialPoints.size(); };
};

struct TrialConstrained : public Trial {
    int nu;

    TrialConstrained(double _x = 0.0, double _z = 0.0, int _nu = 0) : Trial(_x, _z), nu(_nu) {}
};

class OptimizationMethodConstrained : public OptimizationMethod {
protected:
    function<double(vector<double>, int)> f; // target function
    int numberConstraints; // number of constraints

    vector<TrialConstrained> trialPoints;

    virtual TrialConstrained newTrial(double x) = 0;
    virtual double newPoint(int t) = 0;
    virtual double selectNewPoint(int &t) = 0;

public:
    OptimizationMethodConstrained(function<double(vector<double>, int)> _f, int _n, int _numberConstarints,
                                  const vector<double> &_A, const vector<double> &_B, double _eps, int _maxTrials,
                                  int _maxFevals) : OptimizationMethod(_n, _A, _B, _eps, _maxTrials, _maxFevals),
                                  f(_f), numberConstraints(_numberConstarints) {}

    void setF(function<double(vector<double>, int)> _f) { f = _f; };
    void setNumberConstraints(int _numberConstraints) { numberConstraints = _numberConstraints; };

    function<double(vector<double>, int)> getF() const { return f; };
    int getNumberConstraints() const { return numberConstraints; };

    void getTrialPoints(vector<TrialConstrained> &trialPoints) const { trialPoints = this->trialPoints; };
    int getNumberTrialPoints() const { return trialPoints.size(); };
};

#endif // OPTIMIZATION_METHOD_H_
