#ifndef MULTI_DIMENSIONAL_TASK_OPTIMIZATION_METHOD_H_
#define MULTI_DIMENSIONAL_TASK_OPTIMIZATION_METHOD_H_

#include <vector>

using std::vector;

template <typename ObjectiveFunctionType>
class MultiDimensionalTaskOptimizationMethod {
protected:
    ObjectiveFunctionType objFunction;
    int dimension;
    vector<double> lowerBound, upBound; // area of search

public:
    MultiDimensionalTaskOptimizationMethod(ObjectiveFunctionType _objFunction, int _dimension,
                                           const vector<double> &_lowerBound, const vector<double> &_upBound)
                                           : objFunction(_objFunction), dimension(_dimension),
                                           lowerBound(_lowerBound), upBound(_upBound) {};

    void setF(ObjectiveFunctionType _objFunction) { objFunction = _objFunction; };
    ObjectiveFunctionType getF() const { return objFunction; };

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

    virtual double computeObjFunction(const vector<double> X) const = 0;
};

#endif // MULTI_DIMENSIONAL_TASK_OPTIMIZATION_METHOD_H_
