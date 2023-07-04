#ifndef I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_
#define I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralNumericalOptimizationMethod.h>

using std::vector;

template <typename TaskOptimizationMethodType, typename TrialType>
class IGeneralCharacteristicOptimizationMethod : public IGeneralNumericalOptimizationMethod<TaskOptimizationMethodType, TrialType> {
protected:
    double t;

    virtual void calcCharacteristic() = 0;

public:
    IGeneralCharacteristicOptimizationMethod(const TaskOptimizationMethodType &_task, double _accuracy, int _maxTrials, int _maxFevals)
                                             : IGeneralNumericalOptimizationMethod<TaskOptimizationMethodType, TrialType>(_task,
                                             _accuracy, _maxTrials, _maxFevals), t(0) {};
};

#endif // I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_
