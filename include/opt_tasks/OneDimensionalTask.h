#ifndef ONE_DIMENSIONAL_TASK_H_
#define ONE_DIMENSIONAL_TASK_H_

#include <functional>

#include <base_classes/opt_tasks/IGeneralOptTask.h>
#include <base_classes/search_areas/OneDimensionalSearchArea.h>

using std::function;

class OneDimensionalTask : public IGeneralOptTask<function<double(double)>, OneDimensionalSearchArea, double> {
private:
    double lipschitzConstant;

public:
    OneDimensionalTask(const function<double(double)> &_objFunction = nullptr,
                       const OneDimensionalSearchArea &_area = OneDimensionalSearchArea(0.0, 1.0),
                       double _optPoint = 0.0, double _lipschitzConstant = -1.0)
                      : IGeneralOptTask<function<double(double)>, OneDimensionalSearchArea, double>(_objFunction,
                      _area, _optPoint), lipschitzConstant(_lipschitzConstant) {};

    void setLipschitzConstant(double _lipschitzConstant) { lipschitzConstant = _lipschitzConstant; };
    double getLipschitzConstant() const { return lipschitzConstant; };

    double computeObjFunction(double x) const override { return objFunction(x); };
};

#endif // ONE_DIMENSIONAL_OPTIMIZATION_TASK_H_
