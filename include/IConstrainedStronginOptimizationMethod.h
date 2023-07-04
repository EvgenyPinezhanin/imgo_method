#ifndef I_CONSTRAINED_STRONGIN_OPTIMIZATION_METHOD_H_
#define I_CONSTRAINED_STRONGIN_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralStronginOptimizationMethod.h>
#include <ConstrainedTrial.h>

using std::vector;

template <typename TaskOptimizationMethodType>
class IConstrainedStronginOptimizationMethod
    : public IGeneralStronginOptimizationMethod<TaskOptimizationMethodType, ConstrainedTrial> {
protected:
    vector<double> r;
    vector<double> constantsEstimation;

    vector<vector<ConstrainedTrial>> I;
    vector<bool> calcI;
    vector<double> mu;
    vector<double> zStar;

public:
    IConstarinedStronginOptimizationMethod(const TaskOptimizationMethodType &_task, vector<double> _r, double _accuracy,
                                           int _maxTrials, int _maxFevals)
                                           : IGeneralStronginOptimizationMethod<TaskOptimizationMethodType, ConstarinedTrial>(_task,
                                           _accuracy, _maxTrials, _maxFevals), r(_r), constantEstimation(0) {};

    void setR(const vector<double> &_r) { r = _r; };
    void getR(vector<double>& _r) const { _r = r; };

    void getConstantEstimation(vector<double> _constantsEstimation) const { _constantEstimation = constantsEstimation; };
};

#endif // I_CONSTRAINED_STRONGIN_OPTIMIZATION_METHOD_H_
