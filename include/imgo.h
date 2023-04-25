#ifndef IMGO_H_
#define IMGO_H_

#include <vector>
#include <functional>

#include <optimization_method.h>
#include <task.h>

using namespace std;

class ImgoMethod : public OptimizationMethodConstrained {
private:
    function<double(double, int)> f;
    double r, d;

    TrialConstrained lastTrial;
    int lastTrialPos;

    vector<vector<TrialConstrained>> I;
    vector<bool> calcI;
    vector<double> mu;
    vector<double> zStar;

    TrialConstrained newTrial(double x) override;
    double newPoint(int t) override;
    double selectNewPoint(int &t) override;

public:
    ImgoMethod(function<double(double, int)> _f = nullptr, int _numberConstraints = 0, double _a = 0.0, double _b = 1.0, 
               double _r = 2.0, double _d = 0.0, double _eps = 0.0001, int _maxTrials = 1000, int _maxFevals = 1000)
               : OptimizationMethodConstrained(nullptr, 1, _numberConstraints, vector<double>{_a}, vector<double>{_b},
               _eps, _maxTrials, _maxFevals), f(_f), r(_r), d(_d), lastTrial(0.0, 0.0, 0), lastTrialPos(0),
               I((size_t)numberConstraints + 1), calcI((size_t)numberConstraints + 1), mu((size_t)numberConstraints + 1),
               zStar((size_t)numberConstraints + 1) {}
    
    void setF(const function<double(double, int)> &_f) { f = _f; };
    void setA(double _a) { OptimizationMethod::setA(vector<double>{ _a }); };
    void setB(double _b) { OptimizationMethod::setB(vector<double>{ _b }); };
    void setAB(double _a, double _b) { OptimizationMethod::setAB(vector<double>{ _a }, vector<double>{ _b }); };
    void setNumberConstraints(int _numberConstraints);
    void setR(double _r) { r = _r; };
    void setD(double _d) { d = _d; };

    function<double(double, int)> getF() const { return f; };
    double getA() const { return A[0]; };
    double getB() const { return B[0]; };
    double getR() const { return r; };
    double getD() const { return d; };

    void getLambda(vector<double> &lambdas) const { lambdas = mu; };

    void solve(int &numberTrials, int &numberFevals, double &x);
    void solve(int &numberTrials, int &numberFevals, vector<double> &X) override;

    bool solveTest(double xOpt, int &numberTrials, int &numberFevals);
};

#endif // IMGO_H_
