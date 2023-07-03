#ifndef I_STRONGIN_OPTIMIZATION_METHOD_H_
#define I_STRONGIN_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralStronginOptimizationMethod.h>
#include <IObjectiveFunction.h>
#include <Trial.h>

using std::vector;

template <typename ObjectiveFunction>
class IStronginOptimizationMethod : public IGeneralStronginOptimizationMethod<Trial>, public IObjectiveFunction<ObjectiveFunction> {
protected:
    double r;

public:
    IStronginOptimizationMethod(int _dimension, const vector<double> &_lowerBound, const vector<double> &_upBound,
                                double _r, double _accuracy, int _maxTrials, int _maxFevals) : IGeneralStronginOptimizationMethod(
                                _dimension, _lowerBound, _upBound, _accuracy, _maxTrials, _maxFevals), r(_r) {};

    void setR(double _r) { r = _r; };

    double getR() const { return r; };
};

#endif // I_STRONGIN_OPTIMIZATION_METHOD_H_
