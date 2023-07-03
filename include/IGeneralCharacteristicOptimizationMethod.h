#ifndef I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_
#define I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralNumericalOptimizationMethod.h>

using std::vector;

template <typename Trial> 
class IGeneralCharacteristicOptimizationMethod : public IGeneralNumericalOptimizationMethod<Trial> {
protected:
    double t;

    virtual void calcCharacteristic() = 0;

public:
    IGeneralCharacteristicOptimizationMethod(int _dimension, const vector<double> &_lowerBound, const vector<double> &_upBound,
                                             double _accuracy, int _maxTrials, int _maxFevals) : IGeneralNumericalOptimizationMethod<Trial>(
                                             _dimension, _lowerBound, _upBound, _accuracy, _maxTrials, _maxFevals), t(0) {};
};

#endif // I_GENERAL_CHARACTERISTIC_OPTIMIZATION_METHOD_H_
