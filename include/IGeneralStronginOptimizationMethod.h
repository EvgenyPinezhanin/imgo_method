#ifndef I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_
#define I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_

#include <vector>

#include <IGeneralCharacteristicOptimizationMethod.h>

using std::vector;

template <typename TaskOptimizationMethodType, typename TrialType>
class IGeneralStronginOptimizationMethod : public IGeneralCharacteristicOptimizationMethod<TaskOptimizationMethodType, TrialType> {
protected:
    vector<TrialType> lastTrials;
    vector<int> lastTrialsPos;

    virtual int insertInSorted(vector<TrialType> &trials, TrialType trial) = 0;
    virtual double searchMin() const = 0;

public:
    IGeneralStronginOptimizationMethod(const TaskOptimizationMethodType &_task, double _accuracy, int _maxTrials, int _maxFevals)
                                       : IGeneralCharacteristicOptimizationMethod<TaskOptimizationMethodType, TrialType>(_task,
                                       _accuracy, _maxTrials, _maxFevals), lastTrials(0), lastTrialsPos(0) {};
};

#endif // I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_
