#ifndef I_GENERAL_OPTIMIZATION_METHOD_H_
#define I_GENERAL_OPTIMIZATION_METHOD_H_

#include <vector>

#include <ResultMethod.h>

using std::vector;

template <typename ObjectiveFunction>
class IGeneralOptimizationMethod {
protected:
    ObjectiveFunction objFunction;
    int dimension;
    vector<double> lowerBound, upBound; // area of search

    int maxFevals, numberFevals;

    virtual double compute(vector<double> X) const = 0;

public:
    IGeneralOptimizationMethod(ObjectiveFunction _objFunction, int _dimension, const vector<double> &_lowerBound,
                               const vector<double> &_upBound, int _maxFevals)
                               : objFunction(_objFunction), dimension(_dimension), lowerBound(_lowerBound),
                               upBound(_upBound), maxFevals(_maxFevals), numberFevals(0) {};

    void setF(ObjectiveFunction _objFunction) { objFunction = _objFunction; };
    ObjectiveFunction getF() const { return objFunction; };

    void setDimension(int _dimension) { dimension = _dimension; };
    int getDimension() const { return dimension; };

    void setLowerBound(const vector<double> &_lowerBound) { lowerBound = _lowerBound; };
    void getLowerBound(vector<double>& _lowerBound) const { _lowerBound = lowerBound; };

    void setUpBound(const vector<double> &_upBound) { upBound = _upBound; };
    void getUpBound(vector<double>& _upBound) const { _upBound = upBound; };

    void setBounds(const vector<double> &_lowerBound, const vector<double> &_upBound) {
        lowerBound = _lowerBound; upBound = _upBound;
    };
    void getBounds(vector<double>& _lowerBound, vector<double>& _upBound) const {
        _lowerBound = lowerBound; _upBound = upBound;
    };

    void setMaxFevals(int _maxFevals) { maxFevals = _maxFevals; };
    int getMaxFevals() const { return maxFevals; };

    virtual void solve(ResultMethod *result) = 0;
    virtual bool solveTest(vector<double> XOpt, ResultMethod *result) = 0;
};

#endif // I_GENERAL_OPTIMIZATION_METHOD_H_
