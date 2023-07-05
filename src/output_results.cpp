/* #include <output_results.h>

#include <iostream>
#include <iomanip>
#include <cmath>

void printResultGsa(string taskName, double a, double b, double lipschitzConst, double xOpt, double fXOpt,
                    int maxTrials, int maxFevals, double eps, double r, int numberTrials, int numberFevals,
                    double estLipschitzConst, double x, double fX) {
    const auto defaultPrecision = cout.precision();
    cout << setprecision(8);

    cout << "Function: " << taskName << "\n";
    cout << "[a; b] = [" << a << "; " << b << "]"<< "\n";
    cout << "Lipschitz constant = " << lipschitzConst << "\n";
    cout << "X* = " << xOpt << "\n";
    cout << "f(X*) = " << fXOpt << "\n";
    cout << "Parameters for method:" << "\n";
    cout << "Maximum of trials = " << maxTrials << "\n";
    cout << "Maximum of fevals = " << maxFevals << "\n";
    cout << "eps = " << eps << " r = " << r << "\n";
    cout << "Trials result:" << "\n";
    cout << "Number of trials = " << numberTrials << "\n";
    cout << "Number of fevals = " << numberFevals << "\n";
    cout << "Estimation of the Lipschitz constant = " << estLipschitzConst << "\n";
    cout << "X = " << x << "\n";
    cout << "f(X) = " << fX << "\n";
    cout << "|X* - X| = " << abs(xOpt - x) << "\n";
    cout << "|f(X*) - f(X)| = " << abs(fXOpt - fX) << "\n";
    cout << endl;

    cout << setprecision(defaultPrecision);
}

void printResultImgo(string taskName, int numberConstraints, double a, double b, const vector<double> &lipschitzConst, double xOpt,
                     double fXOpt, int maxTrials, int maxFevals, double eps, double r, double d, int numberTrials, int numberFevals,
                     const vector<double> &estLipschitzConst, double x, double fX) {
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
    cout << "f(X*) = " << setprecision(8) << fXOpt << "\n";
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
    cout << "f(X) = " << fX << "\n";
    cout << "|X* - X| = " << abs(xOpt - x) << "\n";
    cout << "|f(X*) - f(X)| = " << abs(fXOpt - fX) << "\n";
    cout << endl;

    cout << setprecision(defaultPrecision);
}

void printResultMggsa(string taskName, int dimension, int numberConstraints, const vector<double> &A, const vector<double> &B,
                      const vector<double> &lipschitzConst, const vector<double> &xOpt, double fXOpt, int maxTrials, int maxFevals,
                      double eps, double r, double d, int den, int key, int incr, int numberTrials, int numberFevals,
                      const vector<double> &estLipschitzConst, const vector<double> &X, double fX) {
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
    cout << "f(X*) = " << fXOpt << "\n";
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
    cout << "f(X) = " << fX << "\n";
    double sum = 0.0;
    for (int i = 0; i < dimension; i++) {
        sum += (xOpt[i] - X[i]) * (xOpt[i] - X[i]);
    }
    cout << "||X* - X|| = " << sqrt(sum) << "\n";
    cout << "|f(X*) - f(X)| = " << abs(fXOpt - fX) << "\n";
    cout << endl;

    cout << setprecision(defaultPrecision);
}

void printResultDirect(string taskName, int dimension, const vector<double> &A, const vector<double> &B, const vector<double> &xOpt,
                       double fXOpt, int maxIters, int maxFevals, double magicEps, double volumeReltol, double sigmaReltol,
                       direct_algorithm algorithm, int numberFevals, const vector<double> &X, double fX) {
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
    cout << "f(X*) = " << fXOpt << "\n";
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
    cout << "f(X) = " << fX << "\n";
    double sum = 0.0;
    for (int i = 0; i < dimension; i++) {
        sum += (xOpt[i] - X[i]) * (xOpt[i] - X[i]);
    }
    cout << "||X* - X|| = " << sqrt(sum) << "\n";
    cout << "|f(X*) - f(X)| = " << abs(fXOpt - fX) << "\n";
    cout << endl;
}

void addPointGnuplot(ofstream &ofstr, double x, double f) {
    ofstr << x << " " << f << "\n";
    ofstr << "\n" << endl;
}

void addPointGnuplot(ofstream &ofstr, const vector<double> &X, double f) {
    size_t dimensionX = X.size();
    for (int i = 0; i < dimensionX; i++) {
        ofstr << X[i] << " ";
    }
    ofstr << f << "\n";
    ofstr << "\n" << endl;
}

void addPointGnuplot(ofstream &ofstr, const vector<double> &X) {
    size_t dimensionX = X.size();
    for (int i = 0; i < dimensionX - 1; i++) {
        ofstr << X[i] << " ";
    }
    ofstr << X[dimensionX - 1] << "\n";
    ofstr << "\n" << endl;
}

void addPointGnuplot(ofstream &ofstr, Trial trial) {
    ofstr << trial.x << " " << trial.z << "\n";
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<double> &X, const vector<double> &f) {
    size_t numberPoints = X.size();
    for (int i = 0; i < numberPoints; i++) {
        ofstr << X[i] << " " << f[i] << "\n";
    }
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X, const vector<double> &f) {
    size_t sizeX = X.size(), dimensionX = X[0].size();
    for (int i = 0; i < sizeX; i++) {
        for (int j = 0; j < dimensionX; j++) {
            ofstr << X[i][j] << " ";
        }
        ofstr << f[i] << "\n";
    }
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X) {
    size_t sizeX = X.size(), dimensionX = X[0].size();
    for (int i = 0; i < sizeX; i++) {
        for (int j = 0; j < dimensionX - 1; j++) {
            ofstr << X[i][j] << " ";
        }
        ofstr << X[i][dimensionX - 1] << "\n";
    }
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<Trial> &trials) {
    size_t numberPoints = trials.size();
    for (int i = 0; i < numberPoints; i++) {
        ofstr << trials[i].x << " " << trials[i].z << "\n";
    }
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<TrialConstrained> &trials) {
    size_t numberPoints = trials.size();
    for (int i = 0; i < numberPoints; i++) {
        ofstr << trials[i].x << " " << trials[i].z << "\n";
    }
    ofstr << "\n" << endl;
}

void addPointsGnuplot(ofstream &ofstr, const vector<vector<double>> &X, const vector<TrialConstrained> &trials) {
    size_t sizeX = X.size(), dimensionX = X[0].size();
    for (int i = 0; i < sizeX; i++) {
        for (int j = 0; j < dimensionX; j++) {
            ofstr << X[i][j] << " ";
        }
        ofstr << trials[i].z << "\n";
    }
    ofstr << "\n" << endl;
}

void initArrayGnuplot(ofstream &ofstr, string nameArray, int numberElements) {
    ofstr << "array " << nameArray << "[" << numberElements << "]" << endl;
}

void setValueInArrayGnuplot(ofstream &ofstr, string nameArray, int index, string element) {
    ofstr << nameArray << "[" << index << "] = \"" << element << "\"" << endl;
}

void drawGraphGnuplot(string nameScript) {
    drawGraphGnuplot(nameScript, vector<int>{});
}

void drawGraphGnuplot(string nameScript, int arg) {
    drawGraphGnuplot(nameScript, vector<int>{ arg });
}

void drawGraphGnuplot(string nameScript, const vector<int> &args) {
    string inputString;

#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    inputString = "chmod +x " + nameScript;
    if (system(inputString.c_str()))
        cerr << "Error chmod" << endl;
#endif

    inputString = "gnuplot -c " + nameScript;
    for (int i = 0; i < args.size(); i++) {
        inputString += " " + to_string(args[i]);
    }
    if (system(inputString.c_str()))
        cerr << "Error gnuplot" << endl;
}
 */