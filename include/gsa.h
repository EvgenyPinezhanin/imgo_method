#ifndef GSA_H_
#define GSA_H_

#include <vector>
#include <functional>

#include <optimization_method.h>
#include <task.h>

using namespace std;

class GsaMethod : public OptimizationMethodNoConstrained {
private:
    function<double(double)> f;
    double r;

    double m; // parameter m

    Trial lastTrial;
    int lastTrialPos;

    Trial newTrial(double x) override;
    double newPoint(int t) override;
    double selectNewPoint(int &t) override;

public:
    GsaMethod(function<double(double)> _f = nullptr, double _a = 0.0, double _b = 10.0, double _r = 2.0, double _eps = 0.001, 
              int _maxTrials = 1000, int _maxFevals = 1000) : OptimizationMethodNoConstrained(nullptr, 1,
              vector<double>{ _a }, vector<double>{ _b }, _eps, _maxTrials, _maxFevals), f(_f), r(_r), m(0),
              lastTrial(0.0, 0.0), lastTrialPos(0) {};
    
    void setF(function<double(double)> _f) { f = _f; };
    void setA(double _a) { OptimizationMethod::setA(vector<double>{ _a }); };
    void setB(double _b) { OptimizationMethod::setB(vector<double>{ _b }); };
    void setAB(double _a, double _b) { OptimizationMethod::setAB(vector<double>{ _a }, vector<double>{ _b }); };
    void setR(double _r) { r = _r; };

    function<double(double)> getF() const { return f; };
    double getA() const { return A[0]; };
    double getB() const { return B[0]; };
    double getR() const { return r; };

    double getLambda() const { return m; };
    
    void solve(int &numberTrials, int &numberFevals, double &x);
    void solve(int &numberTrials, int &numberFevals, vector<double> &X) override;
};

#endif // GSA_H_
