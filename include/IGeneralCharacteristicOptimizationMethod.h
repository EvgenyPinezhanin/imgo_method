#ifndef I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_
#define I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralNumericalOptimizationMethod.h>

using std::vector;

template <typename ObjectiveFunction, typename Trial> 
class IGeneralCharacteristicOptimizationMethod : public IGeneralNumericalOptimizationMethod<ObjectiveFunction, Trial> {
protected:
    double t;

    virtual void calcCharacteristic() = 0;

public:
    IGeneralCharacteristicOptimizationMethod(ObjectiveFunction _objFunction, int _dimension, const vector<double> &_lowerBound,
                                             const vector<double> &_upBound, double _accuracy, int _maxTrials, int _maxFevals)
                                             : IGeneralNumericalOptimizationMethod<ObjectiveFunction, Trial>(_objFunction,
                                             _dimension, _lowerBound, _upBound, _accuracy, _maxTrials, _maxFevals), t(0) {};
};

#endif // I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_
