#ifndef SCANNING_REPORT_H_
#define SCANNING_REPORT_H_

#include <ostream>

#include <base_classes/IGeneralReport.h>
#include <base_structures/GeneralParametersNumericalOptMethod.h>
#include <opt_methods/ScanningMethod.h>
#include <opt_problems/OneDimensionalProblem.h>
#include <result_methods/ResultMethod.h>
#include <tasks/ScanningTask.h>

class ScanningReport: public opt::IGeneralReport<ScanningTask, OneDimensionalProblem, ResultMethod> {
protected:
    void printOptProblem(std::ostream &stream, const OneDimensionalProblem &optProblem) const override;
    void printMethodParameters(std::ostream &stream,
                               const opt::GeneralParametersNumericalOptMethod &parameters) const override;
    void printResultMethod(std::ostream &stream, const ResultMethod &result) const override;
    void printErrorEstimate(std::ostream &stream, const OneDimensionalProblem &optProblem,
                            const ResultMethod &result) const override;

public:
    ScanningReport():
        opt::IGeneralReport<ScanningTask, OneDimensionalProblem, ResultMethod>()
    {};
};

#endif // SCANNING_REPORT_H_

/* void printResultGsa(string taskName, double a, double b, double lipschitzConst, double xOpt, double optimalF,
                    int maxTrials, int maxFevals, double eps, double r, int numberTrials, int numberFevals,
                    double estLipschitzConst, double x, double f) {
    const auto defaultPrecision = cout.precision();
    cout << setprecision(8);

    cout << "Function: " << taskName << "\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< "\n";
    cout << "Lipschitz constant = " << lipschitzConst << "\n";
    cout << "X* = " << xOpt << "\n";
    cout << "f(X*) = " << optimalF << "\n";
    cout << "Parameters for method:" << "\n";
    cout << "Maximum of trials = " << maxTrials << "\n";
    cout << "Maximum of fevals = " << maxFevals << "\n";
    cout << "eps = " << eps << " r = " << r << "\n";
    cout << "Trials result:" << "\n";
    cout << "Number of trials = " << numberTrials << "\n";
    cout << "Number of fevals = " << numberFevals << "\n";
    cout << "Estimation of the Lipschitz constant = " << estLipschitzConst << "\n";
    cout << "X = " << x << "\n";
    cout << "f(X) = " << f << "\n";
    cout << "|X* - X| = " << abs(xOpt - x) << "\n";
    cout << "|f(X*) - f(X)| = " << abs(optimalF - f) << "\n";
    cout << endl;

    cout << setprecision(defaultPrecision);
}

void printResultImgo(string taskName, int numberConstraints, double a, double b, const vector<double> &lipschitzConst, double xOpt,
                     double optimalF, int maxTrials, int maxFevals, double eps, double r, double d, int numberTrials, int numberFevals,
                     const vector<double> &estLipschitzConst, double x, double f) {
    const auto defaultPrecision = cout.precision();
    cout << setprecision(8);

    cout << "Function: " << taskName << "\n";
    cout << "Number of constraints = " << numberConstraints << "\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< "\n";
    cout << "Lipschitz constant:" << "\n";
    cout << "L*(f) = " << lipschitzConst[numberConstraints] << "\n";
    for (int j = 0; j < numberConstraints; j++) {
        cout << "L*(g" << j + 1 << ") = " << lipschitzConst[j] << "\n";
    }
    cout << "X* = " << setprecision(8) << xOpt << "\n";
    cout << "f(X*) = " << setprecision(8) << optimalF << "\n";
    cout << "Parameters for method:" << endl;
    cout << "Maximum of trials = " << maxTrials << "\n";
    cout << "Maximum of fevals = " << maxFevals << "\n";
    cout << "eps = " << eps << " r = " << r << " d = " << d << "\n";
    cout << "Trials result:" << "\n";
    cout << "Number of trials = " << numberTrials << "\n";
    cout << "Number of fevals = " << numberFevals << "\n";
    cout << "Estimation of the Lipschitz constant:" << "\n";
    cout << "L(f) = " << estLipschitzConst[numberConstraints] << "\n";
    for (int j = 0; j < numberConstraints; j++) {
        cout << "L(g" << j + 1 << ") = " << estLipschitzConst[j] << "\n";
    }
    cout << "X = " << x << "\n";
    cout << "f(X) = " << f << "\n";
    cout << "|X* - X| = " << abs(xOpt - x) << "\n";
    cout << "|f(X*) - f(X)| = " << abs(optimalF - f) << "\n";
    cout << endl;

    cout << setprecision(defaultPrecision);
}

void printResultMggsa(string taskName, int dimension, int numberConstraints, const vector<double> &A, const vector<double> &B,
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