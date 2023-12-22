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

/* void printResultMggsa(string taskName, int dimension, int numberConstraints, const vector<double> &A, const vector<double> &B,
                      const vector<double> &lipschitzConst, const vector<double> &xOpt, double optimalF, int maxTrials, int maxFevals,
                      double eps, double r, double d, int den, int key, int incr, int numberTrials, int numberFevals,
                      const vector<double> &estLipschitzConst, const vector<double> &X, double f) {
    const auto defaultPrecision = cout.precision();
    cout << setprecision(8);

    cout << "Function: " << taskName << "\n";
    cout << "Dimension = " << dimension << "\n";
    cout << "Number of constraints = " << numberConstraints << "\n";
    cout << "[A; B] = [(";
    for (int i = 0; i < dimension - 1; i++) {
        cout << A[i] << ", ";
    }
    cout << A[dimension - 1] << "); (";
    for (int i = 0; i < dimension - 1; i++) {
        cout << B[i] << ", ";
    }
    cout << B[dimension - 1] << ")]" << "\n";
    if (!lipschitzConst.empty()) {
        cout << "Lipschitz constant:" << "\n";
        cout << "L*(f) = " << lipschitzConst[numberConstraints] << "\n";
        for (int j = 0; j < numberConstraints; j++) {
            cout << "L*(g" << j + 1 << ") = " << lipschitzConst[j] << "\n";
        }
    }
    cout << "X* = (";
    for (int i = 0; i < dimension - 1; i++) {
        cout << xOpt[i] << ", ";
    }
    cout << xOpt[dimension - 1] << ")" << "\n";
    cout << "f(X*) = " << optimalF << "\n";
    cout << "Parameters for method:" << "\n";
    cout << "Maximum of trials = " << maxTrials << "\n";
    cout << "Maximum of fevals = " << maxFevals << "\n";
    cout << "eps = " << eps << " r = " << r << " d = " << d << "\n";
    cout << "Parameters for constructing the Peano curve:" << "\n";
    cout << "m = " << den << " key = " << key;
    if (incr >= 0) cout << " incr = " << incr << "\n";
        else cout << "\n";
    cout << "Trials result:" << "\n";
    cout << "Number of trials = " << numberTrials << "\n";
    cout << "Number of fevals = " << numberFevals << "\n";
    cout << "Estimation of the Lipschitz constant:" << "\n";
    cout << "L(f) = " << estLipschitzConst[numberConstraints] << "\n";
    for (int j = 0; j < numberConstraints; j++) {
        cout << "L(g" << j + 1 << ") = " << estLipschitzConst[j] << "\n";
    }
    cout << "X = (";
    for (int i = 0; i < dimension - 1; i++) {
        cout << X[i] << ", ";
    }
    cout << X[dimension - 1] << ")" << "\n";
    cout << "f(X) = " << f << "\n";
    double sum = 0.0;
    for (int i = 0; i < dimension; i++) {
        sum += (xOpt[i] - X[i]) * (xOpt[i] - X[i]);
    }
    cout << "||X* - X|| = " << sqrt(sum) << "\n";
    cout << "|f(X*) - f(X)| = " << abs(optimalF - f) << "\n";
    cout << endl;

    cout << setprecision(defaultPrecision);
}

void printResultDirect(string taskName, int dimension, const vector<double> &A, const vector<double> &B, const vector<double> &xOpt,
                       double optimalF, int maxIters, int maxFevals, double magicEps, double volumeReltol, double sigmaReltol,
                       direct_algorithm algorithm, int numberFevals, const vector<double> &X, double f) {
    cout << "Function: " << taskName << "\n";
    cout << "Dimension = " << dimension << "\n";
    cout << "[A; B] = [(";
    for (int i = 0; i < dimension - 1; i++) {
        cout << A[i] << ", ";
    }
    cout << A[dimension - 1] << "); (";
    for (int i = 0; i < dimension - 1; i++) {
        cout << B[i] << ", ";
    }
    cout << B[dimension - 1] << ")]" << "\n";
    cout << "X* = (";
    for (int i = 0; i < dimension - 1; i++) {
        cout << xOpt[i] << ", ";
    }
    cout << xOpt[dimension - 1] << ")" << "\n";
    cout << "f(X*) = " << optimalF << "\n";
    cout << "Parameters for method:" << "\n";
    cout << "Maximum of iters = " << maxIters << "\n";
    cout << "Maximum of fevals = " << maxFevals << "\n";
    cout << "Magic eps = " << magicEps << "\n";
    cout << "Volume reltol = " << volumeReltol << " Sigma reltol = " << sigmaReltol << "\n";
    cout << "Type of algorithm: " << ((algorithm == DIRECT_ORIGINAL) ? "DIRECT_ORIGINAL" : "DIRECT_GABLONSKY") << "\n";
    cout << "Trials result:" << "\n";
    cout << "Number of fevals = " << numberFevals << "\n";
    cout << "X = (";
    for (int i = 0; i < dimension - 1; i++) {
        cout << X[i] << ", ";
    }
    cout << X[dimension - 1] << ")" << "\n";
    cout << "f(X) = " << f << "\n";
    double sum = 0.0;
    for (int i = 0; i < dimension; i++) {
        sum += (xOpt[i] - X[i]) * (xOpt[i] - X[i]);
    }
    cout << "||X* - X|| = " << sqrt(sum) << "\n";
    cout << "|f(X*) - f(X)| = " << abs(optimalF - f) << "\n";
    cout << endl;
} */

#ifndef MGGSA_RESULT_METHOD_H_
#define MGGSA_RESULT_METHOD_H_

#include <vector>

/* #include <base_structures/result_methods/GeneralNumericalResultMethod.h>

using std::vector;

enum class MggsaStopCriteria { accuracy, error, maxTrials, maxFevals, coincidePoints };

using MggsaResultMethod = GeneralNumericalResultMethod<MggsaStopCriteria, vector<double>>; */

#endif // MGGSA_RESULT_METHOD_H_
