#ifndef MGGSA_H
#define MGGSA_H

#include <vector>
#include <functional>

// #include <AbstractOptimizationMethods/IConstrainedStronginOptimizationMethod.h>

/* using namespace std;

enum class TypeSolve { SOLVE, RESOLVE };

class MggsaMethod : public OptimizationMethodConstrained {
private:
    double r, d;
    int den, key;
    int incr;

    vector<vector<double>> points;

    vector<TrialConstrained> lastTrials;
    vector<int> lastTrialsPos;

    vector<vector<TrialConstrained>> I;
    vector<double> hNu;
    vector<bool> calcI;
    vector<double> mu;
    vector<double> zStar;

    int M;

    bool coincideX;

    TrialConstrained newTrial(double x) override;
    double newPoint(int t) override;
    double selectNewPoint(int &t) override;

    double calcH(double xNew);
    bool checkDensity(double h);
    void y(double x, vector<double> &X);
    void x(const vector<double> &P, vector<double> &X);

public:
    MggsaMethod(function<double(vector<double>, int)> _f = nullptr, int _n = 2, int _numberConstraints = 0,
                const vector<double>& _A = vector<double>(), const vector<double>& _B = vector<double>(),
                double _r = 2.0, double _d = 0.01, int _den = 10, int _key = 1, double _eps = 0.0001,
                int _maxTrials = 1000, int _maxFevals = 1000, int _incr = 0) : OptimizationMethodConstrained(_f,
                _n, _numberConstraints, _A, _B, _eps, _maxTrials, _maxFevals), r(_r), d(_d), den(_den), key(_key),
                incr(_incr), points(0), lastTrials(1), lastTrialsPos(1), M(0), coincideX(false),
                I((size_t)numberConstraints + 1), hNu(0), calcI((size_t)numberConstraints + 1),
                mu((size_t)numberConstraints + 1), zStar((size_t)numberConstraints + 1) {}

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

    bool getCoincideX() const { return coincideX; };
    void getPoints(vector<vector<double>> &points) const { points = this->points; };
    void getL(vector<double> &L) const { L = mu; };

    void solve(int &numberTrials, int &numberFevals, vector<double> &X, TypeSolve type);
    void solve(int &numberTrials, int &numberFevals, vector<double> &X) override;

    bool solveTest(vector<double> XOpt, int &numberTrials, int &numberFevals, TypeSolve type);
    bool solveTest(vector<double> XOpt, int &numberTrials, int &numberFevals);
}; */

#endif // MGGSA_H
