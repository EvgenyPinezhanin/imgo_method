#ifndef I_CONSTRAINED_STRONGIN_OPTIMIZATION_METHOD_H_
#define I_CONSTRAINED_STRONGIN_OPTIMIZATION_METHOD_H_

#include <vector>

#include <base_classes/opt_methods/IGeneralStronginOptimizationMethod.h>
#include <Trials/ConstrainedTrial.h>

using std::vector;

template <typename SolutionType, typename OptimizationTaskType, typename PointType>
class IConstrainedStronginOptimizationMethod
    : public IGeneralStronginOptimizationMethod<SolutionType, ConstrainedTrial, OptimizationTaskType, PointType> {
protected:
    vector<double> r;
    vector<double> constantsEstimation;

    vector<vector<ConstrainedTrial>> I;
    vector<bool> calcI;
    vector<double> mu;
    vector<double> zStar;

public:
    IConstrainedStronginOptimizationMethod(const OptimizationTaskType &_task, vector<double> _r, double _accuracy,
                                           int _maxTrials, int _maxFevals)
            : IGeneralStronginOptimizationMethod<SolutionType, ConstrainedTrial, OptimizationTaskType, PointType>(
                                           _task, _accuracy, _maxTrials, _maxFevals, 0), r(_r), constantsEstimation(0) {};

    void setR(const vector<double> &_r) { r = _r; };
    void getR(vector<double>& _r) const { _r = r; };

    void getConstantEstimation(const vector<double> &_constantsEstimation) const { _constantsEstimation = constantsEstimation; };
};

#endif // I_CONSTRAINED_STRONGIN_OPTIMIZATION_METHOD_H_
