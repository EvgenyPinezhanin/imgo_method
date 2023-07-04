#ifndef TASK_OPTIMIZATION_METHOD_H_
#define TASK_OPTIMIZATION_METHOD_H_

#include <vector>

using std::vector;

template <typename ObjectiveFunction>
class TaskOptimizationMethod {
private:
    ObjectiveFunction objFunction;
    int dimension;
    vector<double> lowerBound, upBound; // area of search

    virtual double compute(vector<double> X) const = 0;

public:
    TaskOptimizationMethod(ObjectiveFunction _objFunction, int _dimension, const vector<double> &_lowerBound,
                           const vector<double> &_upBound)
                           : objFunction(_objFunction), dimension(_dimension), lowerBound(_lowerBound),
                           upBound(_upBound) {};

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
};

#endif // TASK_OPTIMIZATION_METHOD_H_
