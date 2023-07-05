#ifndef I_ONE_DIMENSIONAL_CONSTRAINED_OPTIMIZATION_TASK_H_
#define I_ONE_DIMENSIONAL_CONSTRAINED_OPTIMIZATION_TASK_H_

#include <base_classes/opt_tasks/IOneDimensionalOptimizationTask.h>

template <typename ObjectiveFunctionType>
class IOneDimensionalConstrainedOptimizationTask : public IOneDimensionalOptimizationTask<ObjectiveFunctionType> {
protected:
    int numberConstraints;

public:
    IOneDimensionalConstrainedOptimizationTask(ObjectiveFunctionType _objFunction, int _numberConstraints,
                                               double _lowerBound, double _upBound)
                                               : IOneDimensionalOptimizationTask<ObjectiveFunctionType>(_objFunction,
                                               _lowerBound, _upBound), numberConstraints(_numberConstraints) {};

    void setNumberConstraints(double _numberConstraints) { numberConstraints = _numberConstraints; };
    int getNumberConstraints() const { return numberConstraints; };

    virtual double computeConstraints(double x, int index) const = 0;
};

#endif // I_ONE_DIMENSIONAL_CONSTRAINED_OPTIMIZATION_TASK_H_
