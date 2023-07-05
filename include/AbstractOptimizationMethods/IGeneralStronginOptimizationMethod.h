#ifndef I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_
#define I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_

#include <vector>

#include <AbstractOptimizationMethods/IGeneralCharacteristicOptimizationMethod.h>

using std::vector;

template <typename SolutionType, typename TrialType, typename TaskOptimizationMethodType, typename ResultMethodType, typename PointType>
class IGeneralStronginOptimizationMethod
    : public IGeneralCharacteristicOptimizationMethod<SolutionType, TrialType, TaskOptimizationMethodType, ResultMethodType, PointType> {
protected:
    vector<TrialType> lastTrials;
    vector<int> lastTrialsPos;

    virtual PointType newPoint() = 0;

    virtual int insertInSorted(TrialType trial) = 0;
    virtual double searchMin() const = 0;

public:
    IGeneralStronginOptimizationMethod(const TaskOptimizationMethodType &_task, double _accuracy, int _maxTrials, int _maxFevals)
    : IGeneralCharacteristicOptimizationMethod<SolutionType, TrialType, TaskOptimizationMethodType, ResultMethodType, PointType>(
                                       _task, _accuracy, _maxTrials, _maxFevals), lastTrials(0), lastTrialsPos(0) {};
};

#endif // I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_
