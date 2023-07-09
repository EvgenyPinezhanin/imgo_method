#ifndef GENERAL_NUMERICAL_RESULT_METHOD_H_
#define GENERAL_NUMERICAL_RESULT_METHOD_H_

#include <base_structures/result_methods/GeneralResultMethod.h>

enum class GeneralNumericalStoppingCondition { accuracy, error, maxTrials, maxFevals };

template <typename StoppingConditionType, typename SolutionType>
struct GeneralNumericalResultMethod : GeneralResultMethod<SolutionType> {
    int numberTrials, numberFevals;

    StoppingConditionType stoppingCondition;

    GeneralNumericalResultMethod(SolutionType _solution = SolutionType(), int _numberTrials = 0, int _numberFevals = 0,
                                 StoppingConditionType _stoppingCondition = StoppingConditionType())
                                : GeneralResultMethod<SolutionType>(_solution), numberTrials(_numberTrials),
                                numberFevals(_numberFevals), stoppingCondition(_stoppingCondition) {};
};

#endif // GENERAL_NUMERICAL_RESULT_METHOD_H_
