#ifndef I_MULTI_DIMENSIONAL_OPTIMIZATION_TASK_H_
#define I_MULTI_DIMENSIONAL_OPTIMIZATION_TASK_H_

#include <vector>

#include <base_classes/opt_tasks/IMultiDimensionalOptimizationTask.h>

using std::vector;

template <typename ObjectiveFunctionType>
class IMultiDimensionalOptimizationTask {
protected:
    int numberConstraints;

public:
    IMultiDimensionalOptimizationTask(ObjectiveFunctionType _objFunction, int _dimension,
                                      const vector<double> &_lowerBound, const vector<double> &_upBound)
                                      : objFunction(_objFunction), dimension(_dimension), lowerBound(_lowerBound),
                                      upBound(_upBound) {};

    void setNumberConstraints(double _numberConstraints) { numberConstraints = _numberConstraints; };
    int getNumberConstraints() const { return numberConstraints; };

    virtual double computeConstraints(const vector<double> &X, int index) const = 0;
};

#endif // I_MULTI_DIMENSIONAL_OPTIMIZATION_TASK_H_
