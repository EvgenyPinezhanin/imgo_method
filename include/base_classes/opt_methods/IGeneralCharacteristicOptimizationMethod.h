#ifndef I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_
#define I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_

#include <vector>

#include <base_classes/opt_methods/IGeneralNumericalOptimizationMethod.h>

using std::vector;

template <typename SolutionType, typename TrialType, typename TaskOptimizationMethodType, typename ResultMethodType, typename PointType>
class IGeneralCharacteristicOptimizationMethod
    : public IGeneralNumericalOptimizationMethod<SolutionType, TrialType, TaskOptimizationMethodType, ResultMethodType, PointType> {
protected:
    double t;

    virtual void calcCharacteristic() = 0;

public:
    IGeneralCharacteristicOptimizationMethod(const TaskOptimizationMethodType &_task, double _accuracy, int _maxTrials, int _maxFevals)
            : IGeneralNumericalOptimizationMethod<SolutionType, TrialType, TaskOptimizationMethodType, ResultMethodType, PointType>(
                                             _task, _accuracy, _maxTrials, _maxFevals), t(0) {};
};

#endif // I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_
