#ifndef GENERAL_NUMERICAL_RESULT_METHOD_H_
#define GENERAL_NUMERICAL_RESULT_METHOD_H_

#include <base_structures/result_methods/GeneralResultMethod.h>

namespace opt {
    enum class GeneralNumericalStoppingCondition { accuracy, error, maxTrials, maxFevals };

    template <typename StoppingConditionType, typename PointType>
    struct GeneralNumericalResultMethod : GeneralResultMethod<PointType> {
        int numberTrials, numberFevals;

        StoppingConditionType stoppingCondition;

        GeneralNumericalResultMethod(const PointType &_point = PointType(), double _value = 0.0, int _numberTrials = 0,
                                     int _numberFevals = 0, StoppingConditionType _stoppingCondition = StoppingConditionType())
                                    : GeneralResultMethod<PointType>(_point, _value), numberTrials(_numberTrials),
                                    numberFevals(_numberFevals), stoppingCondition(_stoppingCondition) {};
    };
}

#endif // GENERAL_NUMERICAL_RESULT_METHOD_H_
