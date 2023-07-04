#ifndef I_STRONGIN_OPTIMIZATION_METHOD_H_
#define I_STRONGIN_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralStronginOptimizationMethod.h>
#include <Trial.h>

using std::vector;

template <typename ObjectiveFunction>
class IStronginOptimizationMethod : public IGeneralStronginOptimizationMethod<ObjectiveFunction, Trial> {
protected:
    double r;
    double constantEstimation;

public:
    IStronginOptimizationMethod(ObjectiveFunction _objFunction, int _dimension, const vector<double> &_lowerBound,
                                const vector<double> &_upBound, double _r, double _accuracy, int _maxTrials,
                                int _maxFevals)
                                : IGeneralStronginOptimizationMethod<ObjectiveFunction, Trial>(_objFunction, _dimension,
                                _lowerBound, _upBound, _accuracy, _maxTrials, _maxFevals), r(_r), constantEstimation(0) {};

    void setR(double _r) { r = _r; };
    double getR() const { return r; };

    double getConstantEstimation() const { return constantEstimation; };
};

#endif // I_STRONGIN_OPTIMIZATION_METHOD_H_
