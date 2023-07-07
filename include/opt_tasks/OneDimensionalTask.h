#ifndef I_ONE_DIMENSIONAL_TASK_H_
#define I_ONE_DIMENSIONAL_TASK_H_

#include <functional>

#include <base_classes/opt_tasks/IOptimizationTask.h>
#include <base_classes/search_areas/OneDimensionalSearchArea.h>

using std::function;

class OneDimensionalTask : public IOptimizationTask<function<double(double)>, OneDimensionalSearchArea, double> {
public:
    OneDimensionalTask(const function<double(double)> &_objFunction = nullptr,
                       const OneDimensionalSearchArea &_area = OneDimensionalSearchArea(0.0, 1.0),
                       double _optPoint = 0.0)
                       : IOptimizationTask<function<double(double)>, OneDimensionalSearchArea, double>(_objFunction,
                       _area, _optPoint) {};

    double computeObjFunction(double x) const override { return objFunction(x); };
};

#endif // I_ONE_DIMENSIONAL_OPTIMIZATION_TASK_H_
