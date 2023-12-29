#ifndef MGGSA_H
#define MGGSA_H

#include <vector>
#include <string>
#include <functional>

#include <base_structures/trials/IndexTrial.h>

enum class StopCriteria { accuracy, error, maxTrials, maxFevals, coincidePoints };

struct Task {
    std::string name;
    int n;
    std::vector<double> A, B, X_opt, L;
    double eps;
    int maxTrials, maxFevals;

    bool used;

    Task(std::string _name, int _n, std::vector<double> _A, std::vector<double> _B,
         std::vector<double> _X_opt, std::vector<double> _L, double _eps, int _maxTrials,
         int _maxFevals, bool _used = true)
        : name(_name), n(_n), A(_A), B(_B), X_opt(_X_opt), L(_L), eps(_eps), maxTrials(_maxTrials),
          maxFevals(_maxFevals), used(_used) {};
};

struct TaskMggsa : public Task {
    double (*f)(std::vector<double>, int);
    int numberConstraints;
    double r, d;
    int den, key, incr;

    TaskMggsa(double (*_f)(std::vector<double>, int), std::string _name, int _n, int _numberConstraints, std::vector<double> _A,
              std::vector<double> _B, std::vector<double> _XOpt, std::vector<double> _L, double _eps, int _maxTrials, int _maxFevals,
              double _r, double _d, int _den, int _key, int _incr, bool _used = true) : Task(_name, _n, _A, _B, _XOpt, _L, _eps,
              _maxTrials, _maxFevals, _used), f(_f), numberConstraints(_numberConstraints), r(_r), d(_d), den(_den),
              key(_key), incr(_incr) {};
};

enum class TypeSolve { SOLVE, RESOLVE };

class MggsaMethod {
private:
    std::function<double(std::vector<double>, int)> f;
    int numberConstraints;

    int dimension;
    std::vector<double> A, B; // area of search
    double eps; // accuracy
    int maxTrials, maxFevals;
    int numberFevals;

    double r, d;
    int den, key;
    int incr;

    std::vector<opt::IndexTrial> trialPoints;
    std::vector<std::vector<double>> points;

    std::vector<opt::IndexTrial> lastTrials;
    std::vector<int> lastTrialsPos;

    std::vector<std::vector<opt::IndexTrial>> I;
    std::vector<double> hNu;
    std::vector<bool> calcI;
    std::vector<double> mu;
    std::vector<double> zStar;

    int M;

    bool coincideX;

    opt::IndexTrial newTrial(double x);
    double newPoint(int t);
    double selectNewPoint(int &t);

    double calcH(double xNew);
    bool checkDensity(double h);
    void y(double x, std::vector<double> &X);
    void x(const std::vector<double> &P, std::vector<double> &X);

public:
    MggsaMethod(std::function<double(std::vector<double>, int)> _f = nullptr, int _n = 2, int _numberConstraints = 0,
                const std::vector<double>& _A = std::vector<double>(), const std::vector<double>& _B = std::vector<double>(),
                double _r = 2.0, double _d = 0.01, int _den = 10, int _key = 1, double _eps = 0.0001,
                int _maxTrials = 1000, int _maxFevals = 1000, int _incr = 0) : f(_f), dimension(_n),
                numberConstraints(_numberConstraints), A(_A), B(_B), eps(_eps), maxTrials(_maxTrials),
                maxFevals(_maxFevals), r(_r), d(_d), den(_den), key(_key),
                incr(_incr), points(0), lastTrials(1), lastTrialsPos(1), M(0), coincideX(false),
                I((size_t)numberConstraints + 1), hNu(0), calcI((size_t)numberConstraints + 1),
                mu((size_t)numberConstraints + 1), zStar((size_t)numberConstraints + 1) {}

    void setF(std::function<double(std::vector<double>, int)> _f) { f = _f; };
    void setM(int _m) { numberConstraints = _m; };
    void setN(int _n) { dimension = _n; };
    void setMaxTrials(int _maxTrials) { maxTrials = _maxTrials; };
    void setMaxFevals(int _maxFevals) { maxFevals = _maxFevals; };
    void setEps(double _eps) { eps = _eps; };
    void setA(const std::vector<double> &_A) { A = _A; };
    void setB(const std::vector<double> &_B) { B = _B; };
    void setAB(const std::vector<double> &_A, const std::vector<double> &_B) { A = _A; B = _B; };
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
    void getPoints(std::vector<std::vector<double>> &points) const { points = this->points; };
    void getL(std::vector<double> &L) const { L = mu; };

    void solve(int &numberTrials, int &numberFevals, std::vector<double> &X, TypeSolve type);
    void solve(int &numberTrials, int &numberFevals, std::vector<double> &X);

    bool solveTest(std::vector<double> XOpt, int &numberTrials, int &numberFevals, TypeSolve type);
    bool solveTest(std::vector<double> XOpt, int &numberTrials, int &numberFevals);
};

void printResultMggsa(std::string taskName, int dimension, int numberConstraints, const std::vector<double> &A, const std::vector<double> &B,
                      const std::vector<double> &lipschitzConst, const std::vector<double> &xOpt, double optimalF, int maxTrials, int maxFevals,
                      double eps, double r, double d, int den, int key, int incr, int numberTrials, int numberFevals,
                      const std::vector<double> &estLipschitzConst, const std::vector<double> &X, double f);

#endif // MGGSA_H

/*
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
