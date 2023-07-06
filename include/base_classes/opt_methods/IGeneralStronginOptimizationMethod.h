#ifndef I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_
#define I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_

#include <vector>

#include <base_classes/opt_methods/IGeneralCharacteristicOptimizationMethod.h>
#include <result_methods/StronginResultMethod.h>

using std::vector;

template <typename SolutionType, typename TrialType, typename OptimizationTaskType, typename PointType>
class IGeneralStronginOptimizationMethod
    : public IGeneralCharacteristicOptimizationMethod<SolutionType, TrialType, OptimizationTaskType, StronginResultMethod<SolutionType>, PointType> {
protected:
    vector<TrialType> lastTrials;
    vector<int> lastTrialsPos;

    StopCriteria stopCriteria;

    virtual PointType newPoint() = 0;

    virtual int insertInSorted(TrialType trial) = 0;
    virtual double searchMin() const = 0;

    void setDataInResultData(StronginResultMethod<SolutionType> &result) const {
        result.numberFevals = IGeneralOptimizationMethod<OptimizationTaskType, StronginResultMethod<SolutionType>, PointType>::numberFevals;
        result.numberTrials = IGeneralNumericalOptimizationMethod<SolutionType, TrialType, OptimizationTaskType, StronginResultMethod<SolutionType>, PointType>::numberTrials;
        result.stopCriteria = stopCriteria;
        result.solution = IGeneralNumericalOptimizationMethod<SolutionType, TrialType, OptimizationTaskType, StronginResultMethod<SolutionType>, PointType>::estimateSolution();
    }

public:
    IGeneralStronginOptimizationMethod(const OptimizationTaskType &_task, double _accuracy, int _maxTrials, int _maxFevals,
                                       int numberLastTrials)
    : IGeneralCharacteristicOptimizationMethod<SolutionType, TrialType, OptimizationTaskType, StronginResultMethod<SolutionType>, PointType>(
                                       _task, _accuracy, _maxTrials, _maxFevals), lastTrials(numberLastTrials),
                                       lastTrialsPos(numberLastTrials), stopCriteria(StopCriteria::accuracy) {};
};

#endif // I_GENERAL_STRONGIN_OPTIMIZATION_METHOD_H_
