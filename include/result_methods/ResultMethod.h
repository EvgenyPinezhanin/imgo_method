#ifndef RESULT_METHOD_H_
#define RESULT_METHOD_H_

#include <base_structures/result_methods/GeneralNumericalResultMethod.h>

using StoppingCondition = opt::GeneralNumericalStoppingCondition;
using ResultMethod = opt::GeneralNumericalResultMethod<StoppingCondition, double>;

#endif // RESULT_METHOD_H_
