#ifndef I_ONE_DIMENSIONAL_TASK_H_
#define I_ONE_DIMENSIONAL_TASK_H_

#include <functional>

#include <base_classes/opt_tasks/IOneDimensionalOptimizationTask.h>

using std::function;

class OneDimensionalTask : public IOneDimensionalOptimizationTask<function<double(double)>> {
public:
    OneDimensionalTask(function<double(double)> _objFunction = nullptr, double _lowerBound = 0.0, double _upBound = 1.0)
                       : IOneDimensionalOptimizationTask<function<double(double)>>(_objFunction, _lowerBound, _upBound) {};

    double computeObjFunction(double x) const override { return objFunction(x); };
};

#endif // I_ONE_DIMENSIONAL_OPTIMIZATION_TASK_H_
