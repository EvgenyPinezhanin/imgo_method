#ifndef STRONGIN_RESULT_METHOD_H_
#define STRONGIN_RESULT_METHOD_H_

#include <base_structures/result_methods/ResultMethod.h>

enum class StopCriteria { accuracy, maxTrials, maxFevals, coincidePoints };

template <typename SolutionType>
struct StronginResultMethod : public ResultMethod<SolutionType> {
    int numberTrials;
    StopCriteria stopCriteria;

    StronginResultMethod(int _numberTrials = 0, int _numberFevals = 0, StopCriteria _stopCriteria = StopCriteria::accuracy,
                         SolutionType _solution = SolutionType{})
                         : ResultMethod<SolutionType>(_numberFevals, _solution), numberTrials(_numberTrials),
                         stopCriteria(_stopCriteria) {};
};

#endif // STRONGIN_RESULT_METHOD_H_
