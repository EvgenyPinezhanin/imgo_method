#ifndef OUTPUT_RESULTS_H_
#define OUTPUT_RESULTS_H_

#include <fstream>
#include <string>
#include <vector>

#include <optimization_method.h>
#include <direct_method.h>

using namespace std;

void printResultGsa(string taskName, double a, double b, double lipschitzConst, double xOpt, double fXOpt, int maxTrials, int maxFevals,
                    double eps, double r, int numberTrials, int numberFevals, double estLipschitzConst, double x, double fX);

void printResultImgo(string taskName, int numberConstraints, double a, double b, const vector<double> &lipschitzConst, double xOpt,
                     double fXOpt, int maxTrials, int maxFevals, double eps, double r, double d, int numberTrials, int numberFevals,
                     const vector<double> &estLipschitzConst, double x, double fX);

void printResultMggsa(string taskName, int dimension, int numberConstraints, const vector<double> &A, const vector<double> &B,
                      const vector<double> &lipschitzConst, const vector<double> &xOpt, double fXOpt, int maxTrials, int maxFevals,
                      double eps, double r, double d, int den, int key, int incr, int numberTrials, int numberFevals,
                      const vector<double> &estLipschitzConst, const vector<double> &X, double fX);

void printResultDirect(string taskName, int dimension, const vector<double> &A, const vector<double> &B, const vector<double> &xOpt,
                       double fXOpt, int maxIters, int maxFevals, double magicEps, double volumeReltol, double sigmaReltol,
                       direct_algorithm algorithm, int numberFevals, const vector<double> &X, double fX);

void addPointGnuplot(ofstream &ofstr, double x, double f);
void addPointGnuplot(ofstream &ofstr, Trial trial);

void addPointGnuplot(ofstream &ofstr, const vector<double> &X, double f);

void addPointsGnuplot(ofstream &ofstr, const vector<double> &X, vector<double> f);
void addPointsGnuplot(ofstream &ofstr, vector<Trial> trials);
void addPointsGnuplot(ofstream &ofstr, vector<TrialConstrained> trials);

void addPointsGnuplot(ofstream &ofstr, vector<vector<double>> X, vector<double> f);
void addPointsGnuplot(ofstream &ofstr, vector<vector<double>> X, vector<TrialConstrained> trials);

template<typename T>
void setVariableGnuplot(ofstream &ofstr, string nameVariable, T value, bool stringExpression = true) {
    ofstr << nameVariable <<  (stringExpression ? " = \"" : " = ") << value << (stringExpression ? "\"" : "") << endl;
}

void initArrayGnuplot(ofstream &ofstr, string nameArray, int numberElements);

template<typename T>
void setValueInArrayGnuplot(ofstream &ofstr, string nameArray, int index, T value, bool stringExpression = true) {
    ofstr << nameArray << "[" << index << (stringExpression ? "] = \"" : "] = ") <<
                                 value << (stringExpression ? "\"" : "") << endl;
}

void drawGraphGnuplot(string nameScript);
void drawGraphGnuplot(string nameScript, int arg);
void drawGraphGnuplot(string nameScript, vector<int> args);

#endif // OUTPUT_RESULTS_H_
