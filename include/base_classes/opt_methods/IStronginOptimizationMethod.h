#ifndef I_STRONGIN_OPTIMIZATION_METHOD_H_
#define I_STRONGIN_OPTIMIZATION_METHOD_H_

#include <vector>

#include <base_classes/opt_methods//IGeneralStronginOptimizationMethod.h>
#include <Trials/Trial.h>

using std::vector;

template <typename SolutionType, typename TaskOptimizationMethodType, typename ResultMethodType, typename PointType>
class IStronginOptimizationMethod
    : public IGeneralStronginOptimizationMethod<SolutionType, Trial, TaskOptimizationMethodType, ResultMethodType, PointType> {
protected:
    double r;
    double constantEstimation;

public:
    IStronginOptimizationMethod(const TaskOptimizationMethodType &_task, double _r, double _accuracy, int _maxTrials, int _maxFevals)
                 : IGeneralStronginOptimizationMethod<SolutionType, Trial, TaskOptimizationMethodType, ResultMethodType, PointType>(
                                _task, _accuracy, _maxTrials, _maxFevals), r(_r), constantEstimation(0) {};

    void setR(double _r) { r = _r; };
    double getR() const { return r; };

    double getConstantEstimation() const { return constantEstimation; };
};

#endif // I_STRONGIN_OPTIMIZATION_METHOD_H_
