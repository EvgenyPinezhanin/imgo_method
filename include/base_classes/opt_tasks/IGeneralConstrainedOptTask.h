#ifndef I_GENERAL_CONSTRAINED_OPT_TASK_H_
#define I_GENERAL_CONSTRAINED_OPT_TASK_H_

#include <base_classes/opt_tasks/IGeneralOptTask.h>

template <typename ObjectiveFunctionType, typename SearchAreaType, typename OptimalPointType>
class IConstrainedOptimizationTask : public IOptimizationTask<ObjectiveFunctionType, SearchAreaType, OptimalPointType> {
protected:
    int numberConstraints;

public:
    IGeneralConstrainedOptTask(const ObjectiveFunctionType &_objFunction, int _numberConstraints,
                               const SearchAreaType &_area, const OptimalPointType &_optPoint)
                               : IOptimizationTask<ObjectiveFunctionType>(_objFunction, _area, _optPoint),
                               numberConstraints(_numberConstraints) {};

    void setNumberConstraints(double _numberConstraints) { numberConstraints = _numberConstraints; };
    int getNumberConstraints() const { return numberConstraints; };

    virtual double computeConstraints(double x, int index) const = 0;
};

#endif // I_GENERAL_CONSTRAINED_OPT_TASK_H_
